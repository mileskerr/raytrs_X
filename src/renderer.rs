#![allow(unused_parens)]
use std::ops::Range;
use crate::space::*;
use crate::scn::*;
use std::time::Instant;
use mat::Material;
use std::ops::Mul;




pub fn render(scene: Scene, width: usize, height: usize) -> Vec<u8> {
    let aabb = AABB::from_verts(&scene.mesh.verts);
    let t0 = Instant::now();
    println!("accelerating...");
    let accel_struct = accelerate(&scene,aabb,2);
    let mut data = vec![];
    println!("spent {}s accelerating",t0.elapsed().as_secs_f32());

    let t0 = Instant::now();
    println!("rendering...");
    for dir in scene.camera.dirs(width, height) {
        let ray = Ray::new(scene.camera.origin, dir);
        let cols = accel_struct.check_ray(&ray);
        if cols.is_empty() {
            data.push(0);
            data.push(0);
            data.push(0);
        }
        else {
            let col = get_nearest(cols);
            if let Collision::Tri(hit) = col {
                let c: Color = get_material(
                    &scene.mats,&scene.mesh.mats,hit.i
                ).shade(
                    &ray,&hit,&scene,&accel_struct
                );
                data.push(c.r);
                data.push(c.g);
                data.push(c.b);
            } else {
                data.push(255);
                data.push(255);
                data.push(255);
            }

        }
    }
    println!("spent {}s rendering",t0.elapsed().as_secs_f32());
    return data;
}

pub fn get_nearest(cols: Vec<Collision>) -> Collision {
    cols.into_iter().min_by_key(|col| {
        match col {
            Collision::Tri(hit) => { 
                (hit.t * 1000.0) as i32 
            }
            _ => { 
                i32::MAX 
            }
        }
    }).unwrap()
}

//the way materials are stored is kinda wacky, so this function exists
fn get_material<'a>(mats: &'a Vec<Box<dyn Material>>, mat_inds: &Vec<(Range<usize>,usize)>, index: usize) -> &'a Box<dyn Material> {
    for ind in mat_inds.iter().filter(|j| j.0.contains(&index)) {
        return &mats[ind.1];
    }
    return &mats[0];
}


fn accelerate<'a>(scene: &'a Scene, aabb: AABB, iters: usize) -> AccelStruct<'a> { //super unoptimized
    const SUBDIVS: usize = 2;

    let mut children = vec![];
    for subdiv in aabb.subdiv(SUBDIVS) {
        if iters == 0 {
            let tris = subdiv.get_tris(&scene.mesh);
            let lights = subdiv.get_lights(&scene.lights);
            if tris.is_empty() && lights.is_empty() { continue; }
            if !tris.is_empty() {
                let slice = MeshSlice {
                    mesh: &scene.mesh,
                    inds: subdiv.get_tris(&scene.mesh),
                };
                children.push(
                    AccelStruct {
                        aabb: subdiv.clone(),
                        child: Box::new(slice),
                    }
                );
            } 
            if !lights.is_empty() {
                children.push(
                    AccelStruct {
                        aabb: subdiv.clone(),
                        child: Box::new(lights),
                    }
                );
            }
        } else {
            children.push(
                accelerate(scene, subdiv, iters-1)
            );
        }
    }
    AccelStruct {
        aabb: aabb,
        child: Box::new(children),
    }
}

#[derive(Clone)]
struct AABB { //axis aligned bounding box
    min: Vec3,
    max: Vec3,
}

impl<'a> AABB {
    fn from_verts(verts: &'a Vec<Vec3>) -> AABB {
        let mut aabb = AABB { min: Vec3::MAX, max: Vec3::MIN };
        for vert in verts {
            aabb.min = aabb.min.min(vert);
            aabb.max = aabb.max.max(vert);
        }
        return aabb;
    }
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
    fn get_tris(&self, mesh: &Mesh) -> Vec<usize> {
        let mut ok_tris = vec![];
        for i in 0..mesh.tris.len() {
            let tri = mesh.get_verts(i);
            if tri_aabb(&tri,&self) { ok_tris.push(i); }
        }
        return ok_tris;
    }
    fn get_lights(&self, lights: &Vec<Light>) -> Vec<Box<dyn AccelNode + 'a>> {
        let mut ok_lights: Vec<Box<dyn AccelNode + 'a>> = vec![];
        for light in lights {
            if sphere_aabb(light.origin,light.radius,&self) { ok_lights.push(Box::new(*light)); }
        }
        return ok_lights;
    }
}



pub trait AccelNode {
    fn check_ray(&self, ray: &Ray) -> Vec<Collision>;
}
pub enum Collision<'a> {
    Tri(TriHit),
    Light(&'a Light),
}
#[derive(Clone,Copy,Debug)]
pub struct TriHit {
    pub t: f64, pub u: f64, pub v: f64, pub i: usize
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

pub struct MeshSlice<'a> {
    pub mesh: &'a Mesh,
    pub inds: Vec<usize>,
} impl<'a> AccelNode for MeshSlice<'a> {
    fn check_ray(&self, ray: &Ray) -> Vec<Collision> {
        let mut hit = vec![];
        for i in &self.inds {
            let tri = &self.mesh.get_verts(*i);
            ray_tri(ray, tri).map(|h| {
                hit.push(Collision::Tri(TriHit {
                    t: h.0, u: h.1, v: h.2, i: *i,
                })); 0
            });
        }
        return hit;
    }
}


pub struct AccelStruct<'a> {
    aabb: AABB,
    child: Box<dyn AccelNode + 'a>,
}
impl<'a> AccelNode for Light {
    fn check_ray(&self, ray: &Ray) -> Vec<Collision> {
        if ray_sphere(ray, self.origin, self.radius) {
            vec![Collision::Light(&self)]
        } else {
            vec![]
        }
    }
}
impl<'a> AccelNode for Vec<Box<dyn AccelNode>> {
    fn check_ray(&self, ray: &Ray) -> Vec<Collision> {
        let mut hit = vec![];
        for node in self {
            hit.append(&mut node.check_ray(ray));
        }
        return hit;
    }
}
impl<'a, T> AccelNode for Vec<T>
where T: AccelNode {
    fn check_ray(&self, ray: &Ray) -> Vec<Collision> {
        let mut hit = vec![];
        for node in self {
            hit.append(&mut node.check_ray(ray));
        }
        return hit;
    }
}
impl<'a> AccelNode for AccelStruct<'a> {
    fn check_ray(&self, ray: &Ray) -> Vec<Collision> {
        let mut hit = vec![];
        if ray_aabb(ray, &self.aabb) {
            hit.append(&mut self.child.check_ray(ray));
        }
        return hit;
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
fn tri_aabb(tri: &[Vec3;3], aabb: &AABB) -> bool {
    let aabb_norms = [Vec3::RIGHT, Vec3::UP, Vec3::FORWARD];
    let mut tri_min = 0.0; let mut tri_max = 0.0;
    let mut aabb_min = 0.0; let mut aabb_max = 0.0;
    for i in 0..3 {
        (tri_min,tri_max) = project(tri,aabb_norms[i]);
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
        (aabb_min, aabb_max) = project(&aabb_verts,tri_norm);
        if aabb_max < tri_offset || aabb_min > tri_offset { return false; }
    }
    

    for tri_edge in tri_edges {
    for aabb_norm in aabb_norms {
        let axis = tri_edge.cross(aabb_norm);
        (aabb_min, aabb_max) = project(&aabb_verts,axis);
        (tri_min, tri_max) = project(tri,axis);
        if aabb_max < tri_min || aabb_min > tri_max { return false; }
    }}

    true
}
fn sphere_aabb(origin: Vec3, radius: f64, aabb: &AABB) -> bool {
    let r2 = radius * radius;
    let mut dmin = 0.0;
    for i in 0..3 {
        if origin[i] < aabb.min[i] {
            dmin += (origin[i] - aabb.min[i]).sqrt();
        } else if origin[i] > aabb.max[i] {
            dmin += (origin[i] - aabb.max[i]).sqrt();
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

