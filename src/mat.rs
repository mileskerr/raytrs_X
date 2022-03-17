#![allow(unused_parens)]

use renderer::*;
use space::*;
use scn::*;

pub trait Material {
    fn shade<'a>( &self, ray: &Ray, col: Box<dyn Collision + 'a>, accel_struct: &AccelStruct ) -> Color;
}
pub struct NormalMaterial;
impl Material for NormalMaterial {
    fn shade<'a>( &self, ray: &Ray, col: Box<dyn Collision+'a>, _: &AccelStruct ) -> Color {
        return col.normal(ray).into();
    }
}


