#![allow(unused_parens)]
use std::ops::Range;
use std::ops::Mul;
use std::time::Instant;

use std::sync::Arc;
use std::sync::Mutex;
use std::sync::mpsc;

extern crate crossbeam;

use self::crossbeam::thread;
use self::crossbeam::atomic::AtomicCell;


use mat;
use mat::Material;
use crate::space::*;
use crate::scn::*;

pub fn compose(scene: Scene,width: usize,height: usize,samples: usize,threads: usize) -> Vec<u8> {
    let accel_struct = {
        let t0 = Instant::now();
        println!("accelerating...");
        let mut tris = Vec::new();
        let mut spheres = Vec::new();
        scene.mesh.tris.iter().for_each(|tri|tris.push(tri));
        scene.mesh.spheres.iter().for_each(|sphere|spheres.push(sphere));
        let accel = accelerate(Geometry {
            aabb: AABB { min: Vec3::new(-10.0,-10.0,-10.0), max: Vec3::new(10.0,10.0,10.0) },
            tris: tris,
            spheres: spheres,
        });
        println!("spent {}s accelerating",t0.elapsed().as_secs_f32());
        accel
    }; 
    let pixels = render_com(&scene,&accel_struct,width,height,samples,threads);
    render_norm(&scene,&accel_struct,width,height);
    render_depth(&scene,&accel_struct,width,height);
    
    let mut data = Vec::with_capacity(width * height);
    for pixel in pixels {
        let c: Color = pixel.into();
        data.push(c.r);
        data.push(c.g);
        data.push(c.b);
    }
    return data;
}

pub fn render_com<'a>(scene: &Scene, accel_struct: &AccelStruct, width: usize, height: usize, samples: usize, threads: usize) -> Vec<Vec3> {
    //higher is much better for large scenes
    const CHUNK_SIZE: usize = 256;

    let t0 = Instant::now();
    println!("rendering...");

    let num_pixels = width * height;
    let chunks = num_pixels/CHUNK_SIZE;
    
    let dirs = scene.camera.dirs(width, height);
    let accel_struct = accel_struct;
    let camera_origin = scene.camera.origin;
    let mut pixels: Vec<AtomicCell<[Option<Vec3>;CHUNK_SIZE]>> = Vec::with_capacity(chunks);
    
    let mut chunk_status: Vec<u8> = Vec::new(); //0=unrendered, 1=in progress, 2=done
    for _ in 0..chunks {
        pixels.push(AtomicCell::new([None;CHUNK_SIZE]));
        chunk_status.push(0); 
    }
    if num_pixels % CHUNK_SIZE != 0 { //chunk at the end for leftover pixels. this won't be totally filled
        pixels.push(AtomicCell::new([None;CHUNK_SIZE]));
    }

    //channel threads use to communicate that they finished their task.
    //main thread will start a new thread occupied with an unrendered chunk
    //upon recieving the message.
    let (tx, rx) = mpsc::channel();

    //threads are started initially by sending the message that all threads
    //have finished doing nothing, and need to be given work.
    for _ in 0..threads { tx.send(None).unwrap(); }

    thread::scope(|s| {
    loop {
        let done = rx.recv().unwrap(); //loop waits to recieve message that a thread is done
        let new_chunk = { //get next chunk to render
            if done.is_some() {
                chunk_status[done.unwrap()] = 2;
            }
            let mut new_chunk: Option<usize> = None;
            for i in 0..chunks {
                if chunk_status[i] == 0 {
                    chunk_status[i] = 1;
                    new_chunk = Some(i);
                    break;
                }
            }
            new_chunk
        };
        { //progress indicator
            let aspect = (height as f32) / (width as f32);
            let line_length: usize = ((chunks as f32) / aspect).sqrt() as usize;
            if new_chunk.is_none() || new_chunk.unwrap() > 0 {
                print!("\x1b[{}A\n",chunks/line_length+2);
            }
            let mut done_chunks = 0;
            for i in 0..chunks {
                if chunk_status[i] == 2 { done_chunks +=1; }
            }
            print!("rendering on {} threads... {}/{}\n",threads,done_chunks,chunks);
            for i in 0..chunks {
                match chunk_status[i] {
                    0 => { print!("░░"); }
                    1 => { print!("▒▒"); }
                    _ => { print!("▓▓"); }
                }
                if (i+1) % line_length == 0 && i+1 < line_length * (chunks/line_length) {
                    print!("\n");
                }
            }
            print!("\n");
        }


        if new_chunk.is_some() { //render the chunk in a new thread, meanwhile restart the loop
            let chunk_index = new_chunk.unwrap();
            let dirs = &dirs;
            let accel_struct = &accel_struct;
            let pixels = &pixels[chunk_index];
            let tx = tx.clone();
            s.spawn(move |_| { //actual rendering code here:
                let mut new_pixels = [None;CHUNK_SIZE];
                for i in 0..CHUNK_SIZE {
                    let dir = dirs[(chunk_index * CHUNK_SIZE) + i];
                    let ray = Ray::new(scene.camera.origin, dir);
                    let col = accel_struct.trace(&ray);
                    if col.is_none() {
                        let c = mat::background(ray.dir);
                        new_pixels[i] = Some(c);
                    }
                    else {
                        let col = col.unwrap();
                        let c = scene.mats[col.mat()].shade(&ray,col,&accel_struct,&scene,1,samples);
                        new_pixels[i] = Some(c);
                    }
                }
                pixels.store(new_pixels);
                tx.send(Some(chunk_index)).unwrap();

            });
        }
        //loop keeps running even after there are no chunks to assign,
        //but it has to stop when they are all finished rendering
        else if !(chunk_status.contains(&1)) {
            break;
        }
    }
    });

    //pixels is currently a vector of arrays of pixels,
    //merge it into a single vector of pixels:
    let mut output = Vec::new(); 
    for thread in pixels {
    for pixel in thread.load() {
        if output.len() >= num_pixels { break; }
        if pixel.is_some() {
            output.push(pixel.unwrap());
        }
        else {
            output.push(Vec3::ZERO);
        }
    }}
    println!("spent {}s rendering",t0.elapsed().as_secs_f32());
    return output;
}

pub fn render_norm<'a>(scene: &Scene, accel_struct: &AccelStruct, width: usize, height: usize) -> Vec<Option<Vec3>> {
    let num_pixels = width * height;
    
    let dirs = scene.camera.dirs(width, height);
    let camera_origin = scene.camera.origin;
    let mut pixels: Vec<Option<Vec3>> = Vec::with_capacity(num_pixels);
    
    for i in 0..num_pixels {
        let ray = Ray::new(scene.camera.origin, dirs[i]);
        let col = accel_struct.trace(&ray);
        if col.is_some() {
            pixels.push(Some(col.unwrap().normal(&ray)));
        } else {
            pixels.push(None);
        }
    }
    return pixels;
}
pub fn render_depth<'a>(scene: &Scene, accel_struct: &AccelStruct, width: usize, height: usize) -> Vec<f64> {
    let num_pixels = width * height;
    
    let dirs = scene.camera.dirs(width, height);
    let camera_origin = scene.camera.origin;
    let mut pixels: Vec<f64> = Vec::with_capacity(num_pixels);
    
    for i in 0..num_pixels {
        let ray = Ray::new(scene.camera.origin, dirs[i]);
        let col = accel_struct.trace(&ray);
        if col.is_some() {
            pixels.push(col.unwrap().depth(&ray));
        } else {
            pixels.push(f64::MAX);
        }
    }
    return pixels;
}


struct AABB { min: Vec3, max: Vec3, }
impl AABB {
    fn subdiv(&self,subdivs: usize) -> Vec<AABB> {
        let s = subdivs;
        let subdivs = subdivs as f64;
        let size = self.max - self.min;
        let sub_size = size/subdivs;

        let mut grid = vec![];
        for x in 0..s { for y in 0..s { for z in 0..s {
            let origin = Vec3::new(x as f64,y as f64,z as f64) * sub_size;
            grid.push(AABB {
                min: self.min + origin,
                max: self.min + origin + Vec3::ONE * sub_size * 1.01,
            });
        }}}
        return grid;
    }
}

fn accelerate<'a>(geom: Geometry<'a>) -> AccelStruct<'a> {
    const SUBDIVS: usize = 2;
    const OBJS_PER: usize = 70;
    let mut children: Vec<Box<dyn AccelNode+Send+Sync+'a>> = Vec::new();
    for subdiv in geom.aabb.subdiv(SUBDIVS) {
        let mut tris = Vec::new();
        let mut spheres = Vec::new();
        for tri in &geom.tris {
            if tri_aabb(&tri.verts[0..3],&subdiv) {
                tris.push(*tri)
            }
        }
        for sphere in &geom.spheres {
            if sphere_aabb(sphere.origin,sphere.radius,&subdiv) {
                spheres.push(*sphere);
            }
        }
        if tris.is_empty() && spheres.is_empty() { continue; }
        let num_objects = tris.len() + spheres.len();
        let new_geom = Geometry {
            aabb: subdiv,
            tris: tris,
            spheres: spheres,
        };
        if num_objects <= OBJS_PER {
            children.push(Box::new(new_geom));
        } else {
            children.push(Box::new(accelerate(new_geom)));
        }
    }
    AccelStruct {
        aabb: geom.aabb,
        children: children,
    }
}

pub trait AccelNode<'a> {
    fn trace(&self, ray: &Ray) -> Option<Box<dyn Collision+Send+Sync+'a>>;
}
pub struct AccelStruct<'a> {
    aabb: AABB,
    children: Vec<Box<dyn AccelNode<'a>+Send+Sync+'a>>,
}
pub struct Geometry<'a> {
    aabb: AABB,
    tris: Vec<&'a Tri>,
    spheres: Vec<&'a Sphere>,
}
impl<'a> AccelNode<'a> for AccelStruct<'a> {
    fn trace(&self, ray: &Ray) -> Option<Box<dyn Collision+Send+Sync+'a>> {
        let mut cols = Vec::new();
        if ray_aabb(ray,&self.aabb) {
            for child in &self.children {
                let col = child.trace(ray);
                if col.is_some() { cols.push(col.unwrap()) };
            }
            cols.into_iter().min_by_key(|col| {
                let depth = col.depth(ray);
                f64_ord(depth)
            })
        } else {
            None
        }
    }
}
impl<'a> AccelNode<'a> for Geometry<'a> {
    fn trace(&self, ray: &Ray) -> Option<Box<dyn Collision+Send+Sync+'a>> {
        if ray_aabb(ray,&self.aabb) {
            let mut cols: Vec<Box<dyn Collision+Send+Sync+'a>> = Vec::new();
            for tri in &self.tris {
                let hit = ray_tri(ray, &tri.verts);
                if hit.is_some() {
                    let hit = hit.unwrap();
                    cols.push(Box::new(TriCol {
                        t: hit.0,
                        u: hit.1,
                        v: hit.2,
                        tri: tri,
                    }));
                }
            }
            for sphere in &self.spheres {
                let hit = ray_sphere(ray, sphere.origin, sphere.radius);
                if hit {
                    cols.push(Box::new(SphereCol{ sphere: sphere}));
                }
            }
            cols.into_iter().min_by_key(|col| {
                let depth = col.depth(ray);
                f64_ord(depth)
            })
        } else {
            None
        }
    }
}

#[derive(Clone,Copy,Debug)]
pub struct Ray {
    pub start: Vec3,
    pub dir: Vec3,
}
impl Ray {
    pub fn new(start: Vec3, dir: Vec3) -> Ray {
        Ray {
            start: start,
            dir: dir,
        }
    }
}
impl Mul<f64> for Ray {
    type Output = Vec3;
    fn mul(self, other: f64) -> Vec3 {
        self.start + (self.dir * other)
    }
}

pub trait Collision {
    fn depth(&self, ray: &Ray) -> f64;
    fn normal(&self, ray: &Ray) -> Vec3 { Vec3::UP }
    fn mat(&self) -> usize { 0 }
}
pub struct SphereCol<'a> { 
    sphere: &'a Sphere
} impl<'a> Collision for SphereCol<'a> {
    fn depth(&self, ray: &Ray) -> f64 { 
        (self.sphere.origin - ray.start).magn()
    }
    fn mat(&self) -> usize {
        self.sphere.mat
    }
}
pub struct TriCol<'a> {
    t: f64, u: f64, v: f64, tri: &'a Tri
} impl<'a> Collision for TriCol<'a> {
    fn depth(&self, ray: &Ray) -> f64 { self.t }
    fn normal(&self, ray: &Ray) -> Vec3 {
        (self.tri.norms[1] * self.u) +
        (self.tri.norms[2] * self.v) +
        (self.tri.norms[0] * (1.0 - self.u - self.v))
    }
    fn mat(&self) -> usize {
        self.tri.mat
    }
}



fn ray_tri(ray: &Ray, tri: &[Vec3;3]) -> Option<(f64,f64,f64)> {
    //Moller-Trumbore algorithm:
    const EPSILON: f64 = 0.000001;
    
    let edge0 = tri[1] - tri[0];
    let edge1 = tri[2] - tri[0];
    
    let h = ray.dir.cross(edge1);
    let a = edge0.dot(h);
    
    if a < EPSILON { return None; }
    
    let f = 1.0/a;
    let s = ray.start - tri[0];
    let u = f * (s.dot(h));
    
    if u < 0.0 || u > 1.0 { return None; }

    
    let q = s.cross(edge0);
    let v = f * (ray.dir.dot(q));
    
    if v < 0.0 || (u + v) > 1.0 { return None; }
    
    let t = f * (edge1.dot(q));
   
    if t > EPSILON {
        return Some(( t, u, v ));
    }
    else { return None; }
}




fn ray_aabb(ray: &Ray, aabb: &AABB) -> bool {
    let mut tmin: f64 = f64::MIN;
    let mut tmax: f64 = f64::MAX;

    for i in 0..3 {
        let dir_rcp = 1.0/ray.dir[i];
        let mut t0 = (aabb.min[i] - ray.start[i]) * dir_rcp;
        let mut t1 = (aabb.max[i] - ray.start[i]) * dir_rcp;
        if dir_rcp < 0.0 { let tmp = t1; t1 = t0; t0 = tmp; }
        tmin = tmin.max(t0.min(t1));
        tmax = tmax.min(t1.max(t0));
        if tmax < tmin { return false; }
    }
    return tmax > 0.0;

}
fn ray_sphere(ray: &Ray, origin: Vec3, radius:f64) -> bool {

    let a = ray.dir.dot(ray.dir);
    let b = (ray.dir * 2.0).dot(ray.start - origin);
    let c = origin.dot(origin) + 
            ray.start.dot(ray.start) - 
            2.0 * origin.dot(ray.start) - 
            radius * radius;

    let dsc = b * b - (4.0 * a * c);

    dsc >= 0.0
}
fn tri_aabb(tri: &[Vec3], aabb: &AABB) -> bool {
    let aabb_norms = [Vec3::RIGHT, Vec3::UP, Vec3::FORWARD];
    for i in 0..3 {
        let (tri_min,tri_max) = project(tri,aabb_norms[i]);
        if tri_max < aabb.min[i] || tri_min > aabb.max[i] { return false; }
    }

    let aabb_verts = [
        aabb.min,
        Vec3::new(aabb.min.x,aabb.max.y,aabb.min.z),
        Vec3::new(aabb.min.x,aabb.min.y,aabb.max.z),
        Vec3::new(aabb.min.x,aabb.max.y,aabb.max.z),
        Vec3::new(aabb.max.x,aabb.min.y,aabb.min.z),
        Vec3::new(aabb.max.x,aabb.max.y,aabb.min.z),
        Vec3::new(aabb.max.x,aabb.min.y,aabb.max.z),
        aabb.max
    ];
    let tri_edges = [tri[0] - tri[1], tri[1] - tri[2], tri[2] - tri[0]];
    let tri_norm = (tri_edges[0]).cross(tri_edges[1]);
    {
        let tri_offset = tri_norm.dot(tri[0]);
        let (aabb_min, aabb_max) = project(&aabb_verts,tri_norm);
        if aabb_max < tri_offset || aabb_min > tri_offset { return false; }
    }
    

    for tri_edge in tri_edges {
    for aabb_norm in aabb_norms {
        let axis = tri_edge.cross(aabb_norm);
        let (aabb_min, aabb_max) = project(&aabb_verts,axis);
        let (tri_min, tri_max) = project(tri,axis);
        if aabb_max < tri_min || aabb_min > tri_max { return false; }
    }}

    true
}
fn sphere_aabb(origin: Vec3, radius: f64, aabb: &AABB) -> bool {
    let r2 = radius * radius;
    let mut dmin = 0.0;
    for i in 0..3 {
        if origin[i] < aabb.min[i] {
            let n = origin[i] - aabb.min[i];
            dmin += n * n;
        } else if origin[i] > aabb.max[i] {
            let n = origin[i] - aabb.max[i];
            dmin += n * n;
        }
    }
    dmin <= r2
}

fn project(points: &[Vec3], axis: Vec3) -> (f64,f64) {
    let mut min = f64::MAX;
    let mut max = f64::MIN;
    for point in points {
        let v = axis.dot(*point);
        min = min.min(v);
        max = max.max(v);
    }
    (min,max)
}

