#![allow(unused_parens)]
extern crate png;

use crate::space::*;
use crate::scn::*;

pub trait Material {
    fn shade( &self, hit: &TriHit, scene: &Scene ) -> Color;
}
pub struct NormalMaterial;
impl Material for NormalMaterial {
    fn shade( &self, hit: &TriHit, scene: &Scene ) -> Color {
        let normals = &scene.mesh.norms;
        let tri = &scene.mesh.tris[hit.i];
        (
            (normals[tri[4]] * hit.u) +
            (normals[tri[5]] * hit.v) +
            (normals[tri[3]] * (1.0 - hit.u - hit.v))
        ).into()
    }
}



pub fn render(scene: Scene, width: usize, height: usize) -> Vec<u8> {
    let aabb = AABB::from_verts(&scene.mesh.verts);
    let accel_struct = accelerate(&scene.mesh,aabb);
    let mut data = vec![];

    for dir in scene.camera.dirs(width, height) {
        let ray = (scene.camera.origin, dir);
        let hit = accel_struct.check_ray(&ray);
        if hit.is_empty() {
            data.push(0);
            data.push(0);
            data.push(0);
        }
        else {
            let mut nearest = hit[0].clone();
            for h in hit { if h.t < nearest.t { nearest = h.clone(); } }
            let hit = nearest;


            let c: Color = scene.get_material(hit.i).shade(&hit,&scene);
            data.push(c.r);
            data.push(c.g);
            data.push(c.b);
        }
    }
    return data;
}


fn accelerate<'a>(mesh: &'a Mesh, aabb: AABB) ->
AccelStruct<Vec<AccelStruct<MeshSlice>>> {
    const SUBDIVS: usize = 6;

    let mut children = vec![];
    for subdiv in aabb.subdiv(SUBDIVS) {
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
        let fudge = Vec3::uniform(0.02); //TODO: figure out why this is neccessary

        let mut ok_tris = vec![];
        for i in 0..mesh.tris.len() {
            let tri = mesh.tris[i];
            let max = self.max+fudge;
            let min = self.min-fudge;
            if (
                mesh.verts[tri[0]] <= max && mesh.verts[tri[0]] >= min ||
                mesh.verts[tri[1]] <= max && mesh.verts[tri[1]] >= min ||
                mesh.verts[tri[2]] <= max && mesh.verts[tri[2]] >= min
            ) {
                ok_tris.push(i);
            }
        }
        return ok_tris;
    }
}



pub trait AccelNode {
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<TriHit>;
}
#[derive(Clone)]
pub struct TriHit { //depth, UV, index of hit triangle
    pub t: f64, pub u: f64, pub v: f64, pub i: usize
}

pub struct MeshSlice<'a> {
    pub mesh: &'a Mesh,
    pub inds: Vec<usize>,
} impl<'a> AccelNode for MeshSlice<'a> {
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<TriHit> {
        let mut hit = vec![];
        for i in &self.inds {
            let inds = &self.mesh.tris[*i];
            let tri = &[self.mesh.verts[inds[0]],self.mesh.verts[inds[1]],self.mesh.verts[inds[2]]];
            ray_tri(ray, tri).map(|h| {
                hit.push(TriHit {
                    t: h.0, u: h.1, v: h.2, i: *i,
                }); 0
            });
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
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<TriHit> {
        let mut hit = vec![];
        for accel_struct in self {
            hit.append(&mut accel_struct.check_ray(ray));
        }
        return hit;
    }
}
impl<T> AccelNode for AccelStruct<T>
where T: AccelNode {
    fn check_ray(&self, ray: &(Vec3,Vec3)) -> Vec<TriHit> {
        let mut hit = vec![];
        if ray_aabb(ray, &self.aabb) {
            hit.append(&mut self.child.check_ray(ray));
        }
        return hit;
    }
}

fn ray_tri(ray: &(Vec3,Vec3), tri: &[Vec3;3]) -> Option<(f64,f64,f64)> {
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
        return Some(( t, u, v ));
    }
    else { return None; }
}




fn ray_aabb(ray: &(Vec3,Vec3), aabb: &AABB) -> bool {
    let mut tmin: f64 = f64::MIN;
    let mut tmax: f64 = f64::MAX;

    if ray.1.x != 0.0 {
        let dir_rcp = 1.0/ray.1.x;
        let mut t0 = (aabb.min.x - ray.0.x) * dir_rcp;
        let mut t1 = (aabb.max.x - ray.0.x) * dir_rcp;

        if dir_rcp < 0.0 { let tmp = t1; t1 = t0; t0 = tmp; }

        tmin = tmin.max(t0.min(t1));
        tmax = tmax.min(t1.max(t0));
        if tmax < tmin { return false; }
    }
    if ray.1.y != 0.0 {
        let dir_rcp = 1.0/ray.1.y;
        let mut t0 = (aabb.min.y - ray.0.y) * dir_rcp;
        let mut t1 = (aabb.max.y - ray.0.y) * dir_rcp;

        if dir_rcp < 0.0 { let tmp = t1; t1 = t0; t0 = tmp; }

        tmin = tmin.max(t0.min(t1));
        tmax = tmax.min(t1.max(t0));
        if tmax < tmin { return false; }
    }
    if ray.1.z != 0.0 {
        let dir_rcp = 1.0/ray.1.z;
        let mut t0 = (aabb.min.z - ray.0.z) * dir_rcp;
        let mut t1 = (aabb.max.z - ray.0.z) * dir_rcp;

        if dir_rcp < 0.0 { let tmp = t1; t1 = t0; t0 = tmp; }

        tmin = tmin.max(t0.min(t1));
        tmax = tmax.min(t1.max(t0));
        if tmax < tmin { return false; }
    }
    return tmax > 0.0;

}

