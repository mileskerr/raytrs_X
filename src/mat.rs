#![allow(unused_parens)]

use renderer::TriHit;
use renderer::Ray;
use space::*;
use scn::*;

pub trait Material {
    fn shade( &self, ray: &Ray, hit: &TriHit, scene: &Scene ) -> Color;
}
pub struct NormalMaterial;
impl Material for NormalMaterial {
    fn shade( &self, _: &Ray, hit: &TriHit, scene: &Scene ) -> Color {
        get_normal(hit, scene).into()
    }
}
pub struct Unlit(pub Color);
impl Material for Unlit {
    fn shade( &self, _: &Ray, _: &TriHit, _: &Scene ) -> Color {
        self.0.clone()
    }
}
pub struct Simple(pub Color);
impl Material for Simple {
    fn shade( &self, ray: &Ray, hit: &TriHit, scene: &Scene ) -> Color {
        let mut lightness = 0.0;
        let pos = ray.clone() * hit.t;
        let normal = get_normal(hit, scene);
        for l in &scene.lights {
            match l {
                Light::Point{origin,intensity} => {
                    let dir = (*origin - pos).unit();
                    let lin = normal.dot(dir).clamp(0.0,1.0);
                    lightness += lin * lin * intensity;
                }
                _ => {}
            }
        }
        self.0.clone() * lightness
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
