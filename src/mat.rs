#![allow(unused_parens)]

use renderer::{Ray,TriHit,AccelStruct,AccelNode};
use space::*;
use scn::*;

pub trait Material {
    fn shade( &self, ray: &Ray, hit: &TriHit, scene: &Scene, accel_struct: &AccelStruct ) -> Color;
}
pub struct NormalMaterial;
impl Material for NormalMaterial {
    fn shade( &self, _: &Ray, hit: &TriHit, scene: &Scene, _: &AccelStruct ) -> Color {
        get_normal(hit, scene).into()
    }
}
pub struct Unlit(pub Color);
impl Material for Unlit {
    fn shade( &self, _: &Ray, _: &TriHit, _: &Scene, _: &AccelStruct ) -> Color {
        self.0.clone()
    }
}
pub struct Simple(pub Color, pub Color);
impl Material for Simple {
    fn shade( &self, ray: &Ray, hit: &TriHit, scene: &Scene, accel_struct: &AccelStruct ) -> Color {
        const S0: f64 = 0.3;
        const S1: f64 = 0.32;
        
        const AOSIZE: f64 = 0.5;
        const AOSTR: f64 = 0.5;

pub struct Lambertian {
    color: Color,
} impl Lambertian {
    pub fn new(color: Color) -> Lambertian {
        Lambertian { color: color }
    }
} impl Material for Lambertian {
    fn shade<'a> (
        &self, in_ray: &Ray, col: Box<dyn Collision+'a>,
        accel_struct: &AccelStruct, scene: &Scene, iter: u8, samples: usize
    ) -> Vec3 {

        const LIGHT_FACTOR: f64 = 2.5; //setting to anything other than one VIOLATES CONSERVATION OF ENERGY!
        if iter == 0 { return Vec3::ZERO; } 

        let pos = *in_ray * col.depth(in_ray);

        let mut light = Vec3::ZERO;

        let normal = col.normal(in_ray);
        let y_unit = normal;
        let x_unit = Vec3::RIGHT.cross(y_unit).unit();
        let z_unit = y_unit.cross(x_unit).unit();

        let view_matrix = Matrix3::new(x_unit,y_unit,z_unit);
        for _ in 0..samples {
            let dir = rand_in_cap(PI);
            let real_dir = view_matrix * dir;
            let out_ray = Ray::new(pos, real_dir);
            let col = accel_struct.trace(&out_ray);
            let mut new_light = Vec3::ZERO;
            let dot = normal.dot(real_dir).clamp(0.0,1.0);
            if col.is_some() {
                let col = col.unwrap();
                new_light = scene.mats[col.mat()].shade(&out_ray,col,accel_struct,scene,iter-1,samples/2);
            } else {
                new_light = background(dir);
            }
            light = light + new_light * dot * LIGHT_FACTOR;
        }
        let hdr: Vec3 = self.color.clone().into();
        return hdr * light/(samples as f64);
    }
}


fn get_normal( hit: &TriHit, scene: &Scene ) -> Vec3 {
    let normals = &scene.mesh.norms;
    let tri = &scene.mesh.tris[hit.i];
    (
        (normals[tri[4]] * hit.u) +
        (normals[tri[5]] * hit.v) +
        (normals[tri[3]] * (1.0 - hit.u - hit.v))
    )
}

pub fn background(dir: Vec3) -> Vec3 {
    Vec3::ONE * dir.y * 0.4
}

