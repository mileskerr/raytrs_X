extern crate png;

use std::fs;
use crate::space::*;
use crate::extras::*;
use std::path::Path;

pub fn render(mesh: Mesh, camera: Camera, width: usize, height: usize) -> Vec<u8> {
    let aabb = AABB::from_verts(&mesh.verts);
    let accel_struct = accelerate(&mesh,aabb);
    let mut data = vec![];

    for dir in camera.dirs(width, height) {
        let ray = (camera.origin, dir);
        let hit = accel_struct.check_ray(&ray);
        if hit.is_empty() {
            data.push(0);
            data.push(0);
            data.push(0);
        }
        else {
            let mut nearest = hit[0].clone();
            for h in hit { if h.0.t < nearest.0.t { nearest = h.clone(); } }
            let hit = nearest;


            let c: Color = get_normal(&mesh.norms, &mesh.tris[hit.1], hit.0.u, hit.0.v).into();
            data.push(c.r);
            data.push(c.g);
            data.push(c.b);
        }
    }
    return data;
}

fn get_normal<'a> (
    normals: &'a Vec<Vec3>, tri: &'a [usize;3], u: f64, v: f64
) -> Vec3 {
    (normals[tri[1]] * u) +
    (normals[tri[2]] * v) +
    (normals[tri[0]] * (1.0 - u - v))
}

fn accelerate<'a>(mesh: &'a Mesh, aabb: AABB) ->
AccelStruct<Vec<AccelStruct<MeshSlice>>> {
    const subdivs: usize = 4;

    let mut children = vec![];
    for subdiv in aabb.subdiv(subdivs) {
        let slice = MeshSlice {
            mesh: &mesh,
            inds: subdiv.get_tris(&mesh),
        };
        if !slice.inds.is_empty() {
            children.push(
                AccelStruct {
                    aabb: subdiv,
                    child: slice,
                }
            );
        }
    }
    AccelStruct {
        aabb: aabb,
        child: children,
    }
}


pub struct AABB {
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
                max: self.min + origin + (Vec3::ONE * sub_size),
            });
        }}}
        return grid;
    }
    fn get_tris(&self, mesh: &Mesh) -> Vec<usize> {
        let mut ok_tris = vec![];
        for i in 0..mesh.tris.len() {
            let tri = mesh.tris[i];
            if (
                mesh.verts[tri[0]] < self.max && mesh.verts[tri[0]] > self.min ||
                mesh.verts[tri[1]] < self.max && mesh.verts[tri[1]] > self.min ||
                mesh.verts[tri[2]] < self.max && mesh.verts[tri[2]] > self.min
            ) {
                ok_tris.push(i);
            }
        }
        return ok_tris;
    }
}






pub trait AccelNode {
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<(TriHit,usize)>;
}

pub struct Mesh {
    pub verts: Vec<Vec3>,
    pub norms: Vec<Vec3>,
    pub tris: Vec<[usize;3]>,
} 

pub struct MeshSlice<'a> {
    pub mesh: &'a Mesh,
    pub inds: Vec<usize>,
} impl<'a> AccelNode for MeshSlice<'a> {
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<(TriHit,usize)> {
        let mut hit = vec![];
        for i in &self.inds {
            let inds = &self.mesh.tris[*i];
            let tri = &[self.mesh.verts[inds[0]],self.mesh.verts[inds[1]],self.mesh.verts[inds[2]]];
            ray_tri(ray, tri).map(|h| { hit.push((h,*i)); 0});
        }
        return hit;
    }
}


pub struct AccelStruct<T>
where T: AccelNode {
    aabb: AABB,
    child: T,
}

impl<T> AccelNode for Vec<AccelStruct<T>>
where T: AccelNode {
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<(TriHit,usize)> {
        let mut hit = vec![];
        for accel_struct in self {
            hit.append(&mut accel_struct.check_ray(ray));
        }
        return hit;
    }
}
impl<T> AccelNode for AccelStruct<T>
where T: AccelNode {
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<(TriHit,usize)> {
        let mut hit = vec![];
        if ray_aabb(ray, &self.aabb) {
            hit.append(&mut self.child.check_ray(ray));
        }
        return hit;
    }
}

#[derive(Clone)]
pub struct TriHit {
    pub t: f64, pub u: f64, pub v: f64,
}
fn ray_tri(ray: &(Vec3,Vec3), tri: &[Vec3;3]) -> Option<TriHit> {
    //Moller-Trumbore algorithm:
    const EPSILON: f64 = 0.000001;
    
    let edge0 = tri[1] - tri[0];
    let edge1 = tri[2] - tri[0];
    
    let h = ray.1.cross(edge1);
    let a = edge0.dot(h);
    
    if a < EPSILON { return None; }
    
    let f = 1.0/a;
    let s = ray.0 - tri[0];
    let u = f * (s.dot(h));
    
    if u < 0.0 || u > 1.0 { return None; }

    
    let q = s.cross(edge0);
    let v = f * (ray.1.dot(q));
    
    if v < 0.0 || (u + v) > 1.0 { return None; }
    
    let t = f * (edge1.dot(q));
    
    if t > EPSILON {
        return Some(TriHit { t: t, u: u, v: v });
    }
    else { return None; }
}




fn ray_aabb(ray: &(Vec3,Vec3), aabb: &AABB) -> bool {
    let mut tmin: f64 = f64::MIN;
    let mut tmax: f64 = f64::MAX;

    let ray_start = ray.0.x;
    let ray_dir = ray.1.x;
    let aabb_min = aabb.min.x;
    let aabb_max = aabb.max.x;
    {
        let dir_rcp = 1.0/ray_dir;
        let mut t0 = (aabb_min - ray_start) * dir_rcp;
        let mut t1 = (aabb_max - ray_start) * dir_rcp;

        if dir_rcp < 0.0 { let tmp = t1; t1 = t0; t0 = tmp; }

        tmin = tmin.max(t0.min(t1));
        tmax = tmax.min(t1.max(t0));
        if tmax < tmin { return false; }
    }
    let ray_start = ray.0.y;
    let ray_dir = ray.1.y;
    let aabb_min = aabb.min.y;
    let aabb_max = aabb.max.y;
    {
        let dir_rcp = 1.0/ray_dir;
        let mut t0 = (aabb_min - ray_start) * dir_rcp;
        let mut t1 = (aabb_max - ray_start) * dir_rcp;

        if dir_rcp < 0.0 { let tmp = t1; t1 = t0; t0 = tmp; }

        tmin = tmin.max(t0.min(t1));
        tmax = tmax.min(t1.max(t0));
        if tmax < tmin { return false; }
    }
    let ray_start = ray.0.z;
    let ray_dir = ray.1.z;
    let aabb_min = aabb.min.z;
    let aabb_max = aabb.max.z;
    {
        let dir_rcp = 1.0/ray_dir;
        let mut t0 = (aabb_min - ray_start) * dir_rcp;
        let mut t1 = (aabb_max - ray_start) * dir_rcp;

        if dir_rcp < 0.0 { let tmp = t1; t1 = t0; t0 = tmp; }

        tmin = tmin.max(t0.min(t1));
        tmax = tmax.min(t1.max(t0));
        if tmax < tmin { return false; }
    }
    return true;

}

