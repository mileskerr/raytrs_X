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
    hemi: Vec<Vec3>,
} impl Lambertian {
    pub fn new(color: Color) -> Lambertian {
        const SAMPLES: usize = 64;

        let mut hemi = Vec::new();
        for i in 0..SAMPLES {
            hemi.push(rand_in_cap(PI));
        }
        Lambertian { color: color, hemi: hemi }
    }
} impl Material for Lambertian {
    fn shade<'a> (
        &self, in_ray: &Ray, col: Box<dyn Collision+'a>,
        accel_struct: &AccelStruct, scene: &Scene, iter: u8, samples: usize
    ) -> Vec3 {
        if iter == 0 { return Vec3::ZERO; } 

        let hdr: Vec3 = self.color.clone().into();
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
            if col.is_some() {
                let col = col.unwrap();
                let col_pos = out_ray * col.depth(&out_ray);
                let dot = normal.dot((col_pos - pos).unit()).clamp(0.0,1.0);
                let new_light = scene.mats[col.mat()].shade(&out_ray,col,accel_struct,scene,iter-1,samples/2);
                light = light + new_light * dot * 2.0;
            }
        }
        return light/(self.hemi.len() as f64);
    }
}


fn rand_in_cap(theta: f64) -> Vec3 {
    let cos_theta = theta.cos();

    let z = rand::random::<f64>() * (1.0-cos_theta) + cos_theta;
    let phi = rand::random::<f64>() * TAU;
    let a = (1.0 - (z * z)).sqrt();

    Vec3::new(a * phi.cos(), a * phi.sin(), z)
}

