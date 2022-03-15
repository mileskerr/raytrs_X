#![allow(unused_parens)]

use renderer::*;
use space::*;
use scn::*;

pub trait Material {
    fn shade( &self, ray: &Ray, hit: &TriHit, scene: &Scene, accel_struct: &AccelStruct ) -> Color;
}
pub struct NormalMaterial;
impl Material for NormalMaterial {
    fn shade( &self, _: &Ray, hit: &TriHit, scene: &Scene, _: &AccelStruct ) -> Color {
        return get_normal(hit, scene).into();
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
