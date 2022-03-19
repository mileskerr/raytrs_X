#![allow(unused_parens)]

use renderer::*;
use space::*;
use scn::*;
use std::f64::consts::{TAU,PI};

pub trait Material {
    fn shade<'a> (
        &self, in_ray: &Ray, col: Box<dyn Collision + 'a>,
        accel_struct: &AccelStruct, scene: &Scene, iter: u8, samples: usize
    ) -> Vec3;
}

pub struct NormalMaterial;
impl Material for NormalMaterial {
    fn shade<'a> (
        &self, in_ray: &Ray, col: Box<dyn Collision+'a>, _: &AccelStruct, _: &Scene, _: u8 ,_: usize
    ) -> Vec3 {
        col.normal(in_ray)
    }
}

pub struct Emissive {
    color: Color,
    intensity: f64
} impl Emissive {
    pub fn new(color: Color, intensity: f64) -> Emissive {
        Emissive { color: color, intensity: intensity, }
    }
} impl Material for Emissive {
    fn shade<'a> (
        &self, ray: &Ray, col: Box<dyn Collision+'a>, _: &AccelStruct, _: &Scene, _: u8, _: usize
    ) -> Vec3 {
        let hdr: Vec3 = self.color.clone().into();
        hdr * self.intensity
    }
}

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

        const LIGHT_FACTOR: f64 = 1.5; //setting to anything other than one VIOLATES CONSERVATION OF ENERGY!
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


fn rand_in_cap(theta: f64) -> Vec3 {
    let cos_theta = theta.cos();

    let z = rand::random::<f64>() * (1.0-cos_theta) + cos_theta;
    let phi = rand::random::<f64>() * TAU;
    let a = (1.0 - (z * z)).sqrt();

    Vec3::new(a * phi.cos(), a * phi.sin(), z)
}

pub fn background(dir: Vec3) -> Vec3 {
    let color0 = Vec3::new(0.5,0.7,1.0);
    let color1 = Vec3::new(0.4,0.3,0.4);
    let t = (dir.y * 0.5)+0.5;
    color0 * t + color1 * (1.0-t)
}

