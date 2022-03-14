use crate::space::*;
use std::ops::Range;
use mat::Material;


pub struct Scene {
    pub mesh: Mesh,
    pub mats: Vec<Box<dyn Material>>,
    pub lights: Vec<Light>,
    pub camera: Camera,
}
pub enum Light {
    Point { origin: Vec3, intensity: f64 },
    Sphere { origin: Vec3, radius: f64, intensity: f64 },
}
pub struct Mesh {
    pub mats: Vec<(Range<usize>,usize)>,
    pub verts: Vec<Vec3>,
    pub norms: Vec<Vec3>,
    pub txs: Vec<Vec3>,
    pub tris: Vec<[usize;9]>,
}
impl Mesh {
    pub fn join(&mut self, mut other: Mesh) {
        let offs = [self.verts.len(),self.norms.len(),self.txs.len()];
        let mut offs_tris = vec![];
        for tri in other.tris {
            let mut new_tri = [0;9];
            for i in 0..3 {
                new_tri[i] = tri[i] + offs[0];
            } for i in 3..6 {
                new_tri[i] = tri[i] + offs[1];
            } for i in 6..9 {
                new_tri[i] = tri[i] + offs[2];
            }
            offs_tris.push(new_tri)
        }
        let mut offs_mats = vec![];
        for i in 0..other.mats.len() {
            offs_mats.push((
                (other.mats[i].0.start+self.tris.len()..other.mats[i].0.end+self.tris.len()),
                other.mats[i].1
            ));
        }
        self.mats.append(&mut offs_mats);
        self.verts.append(&mut other.verts);
        self.norms.append(&mut other.norms);
        self.txs.append(&mut other.txs);
        self.tris.append(&mut offs_tris);
    }
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
