#![allow(unused_parens)]

use renderer::*;
use space::*;
use scn::*;

pub trait Material {
    fn shade<'a>( &self, ray: &Ray, col: Box<dyn Collision + 'a>, accel_struct: &AccelStruct ) -> Vec3;
}
pub struct NormalMaterial;
impl Material for NormalMaterial {
    fn shade<'a>( &self, ray: &Ray, col: Box<dyn Collision+'a>, _: &AccelStruct ) -> Vec3 {
        return col.normal(ray).into();
    }
}

pub struct Emissive {
    color: Color,
    intensity: f64
}
impl Emissive {
    pub fn new(color: Color, intensity: f64) -> Emissive {
        Emissive {
            color: color,
            intensity: intensity,
        }
    }
}
impl Material for Emissive {
    fn shade<'a>( &self, ray: &Ray, col: Box<dyn Collision+'a>, _: &AccelStruct ) -> Vec3 {
        let hdr: Vec3 = self.color.clone().into();
        hdr * self.intensity
    }
}


