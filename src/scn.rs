use crate::space::*;
use mat::Material;


pub struct Scene {
    pub mesh: Mesh,
    pub mats: Vec<(usize,Box<dyn Material>)>,
    pub lights: Vec<Light>,
    pub camera: Camera,
}

pub enum Light {
    Point{ origin: Vec3, intensity: f64 },
    Sphere{ origin: Vec3, radius: f64, intensity: f64 },
}

pub struct Mesh {
    pub verts: Vec<Vec3>,
    pub norms: Vec<Vec3>,
    pub tris: Vec<[usize;9]>,
}
pub struct Camera {
    pub origin: Vec3,
    pub direction: Vec3,
    pub length: f64,
}
impl Camera {
    pub fn new( origin: Vec3, direction: Vec3, length: f64) -> Camera {
        Camera { origin: origin, direction: direction, length: length }
    }
    pub fn dirs(&self, width: usize, height: usize) -> Vec<Vec3> {

        //first do matrix math to transform the easy-to-understand
        //camera properties into something that's actually useful:
        let z_unit = self.direction.unit();
        let x_unit = Vec3::UP.cross(z_unit).unit();
        let y_unit = z_unit.cross(x_unit).unit();

        let view_matrix = Matrix3::new(x_unit,y_unit,z_unit);

        let aspect = (width as f64) / (height as f64);
        let half = aspect/2.0;

        let upper_left  = Vec3::new(-half, 0.5,self.length);
        let upper_right = Vec3::new(half,  0.5,self.length);
        let lower_right = Vec3::new(half, -0.5,self.length);



        //iterate through pixels and calculate their direction:
        let mut dir = upper_left;

        let dx = (upper_right.x - upper_left.x)/(width as f64);
        let dy = (upper_right.y - lower_right.y)/(height as f64);

        let num_pixels = width * height;
        let mut dirs = Vec::with_capacity(num_pixels);

        for _y in 0..height {
            dir.x = upper_left.x;
            for _x in 0..width {
                dirs.push(view_matrix * dir);
                dir.x += dx;
            }
            dir.y -= dy;
        }
        return dirs;

    }
}
