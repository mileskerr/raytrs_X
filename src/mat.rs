use renderer::TriHit;
use space::*;
use scn::*;

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
pub struct Unlit(pub Color);
impl Material for Unlit {
    fn shade( &self, _: &TriHit, _: &Scene ) -> Color {
        self.0.clone()
    }
}
