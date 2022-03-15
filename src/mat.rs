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

        let mut lightness = 0.0;
        let mut highlight = 0.0;

        let pos = ray.clone() * hit.t;
        let normal = get_normal(hit, scene);
        for l in &scene.lights {
            match l {
                Light::Point{origin,intensity} => {
                    let light_dir = (*origin - pos).unit();
                    let shadow_hit = accel_struct.check_ray(&Ray::new(ray.clone() * hit.t,light_dir));
                    if shadow_hit.is_empty() {
                        let refl_dir = ray.dir.reflect(normal);
                        let lin = normal.dot(light_dir).clamp(0.0,1.0);
                        lightness += lin * lin * intensity;
                        highlight += ((refl_dir.dot(light_dir) -0.5-S0)/S1).clamp(0.0,1.0)
                    } else {
                        let h = &shadow_hit[0];
                        lightness = ((h.t/AOSIZE).clamp(0.0,1.0)-1.0)*AOSTR
                    }
                }
                _ => {}
            }
        }
        if lightness >= 0.0 {
            (self.0.clone() * lightness) +
            (self.1.clone() * (1.0-lightness)) +
            Color::WHITE * highlight
        } else {
            self.1.clone() * (lightness+1.0)
        }
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
