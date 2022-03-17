use std::ops::Mul;
use std::ops::Add;
use std::ops::Sub;
use std::ops::Neg;
use std::ops::Div;
use std::ops::Index;
use std::cmp::Ordering;

#[derive(Clone,Copy,Debug,PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Vec3 {
    pub const ZERO: Vec3 = Vec3 {x:0.0,y:0.0,z:0.0};
    pub const ONE: Vec3 = Vec3 {x:1.0,y:1.0,z:1.0};
    pub const RIGHT: Vec3 = Vec3 {x:1.0,y:0.0,z:0.0};
    pub const UP: Vec3 = Vec3 {x:0.0,y:1.0,z:0.0};
    pub const FORWARD: Vec3 = Vec3 {x:0.0,y:0.0,z:1.0};
    pub const MAX: Vec3 = Vec3 {x:f64::MAX,y:f64::MAX,z:f64::MAX};
    pub const MIN: Vec3 = Vec3 {x:f64::MIN,y:f64::MIN,z:f64::MIN};

    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3{ x: x as f64, y: y as f64, z: z as f64 }
    }
    pub fn uniform(value: f64) -> Vec3 {
        Vec3{ x: value, y: value, z: value }
    }
    pub fn dot(self, other: Vec3) -> f64 {
        &self.x * other.x +
        &self.y * other.y +
        &self.z * other.z
    }
    pub fn cross(self, other: Vec3) -> Vec3 {
        let x = self.y * other.z - self.z * other.y;
        let y = self.z * other.x - self.x * other.z;
        let z = self.x * other.y - self.y * other.x;
        Vec3::new(x,y,z)
    }
    pub fn reflect(self, axis: Vec3) -> Vec3 {
        (axis * (axis.dot(self))) * -2.0 + self
    }
    pub fn magn(&self) -> f64 {
        self.dot(*self).sqrt()
    }
    pub fn unit(self) -> Vec3 {
        self/self.magn()
    }
    pub fn max(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.x.max(other.x),
            y: self.y.max(other.y),
            z: self.z.max(other.z),
        }
    }
    pub fn min(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.x.min(other.x),
            y: self.y.min(other.y),
            z: self.z.min(other.z),
        }
    }
}
impl Index<usize> for Vec3 {
    type Output = f64;
    fn index(&self, i: usize) -> &f64 {
        match i {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vec3 index out of range")
        }
    }
}
impl Mul<Vec3> for Vec3 { //this implimentation makes mathmaticians cry
    type Output = Vec3;
    fn mul(self, other: Vec3) -> Vec3 {
        let x = self.x * other.x;
        let y = self.y * other.y;
        let z = self.z * other.z;
        Vec3::new(x,y,z)
    }
}
impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, other: f64) -> Vec3 {
        let x = self.x * other;
        let y = self.y * other;
        let z = self.z * other;
        Vec3::new(x,y,z)
    }
}
impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, other: f64) -> Vec3 {
        let factor = 1.0 / other;
        let x = self.x * factor;
        let y = self.y * factor;
        let z = self.z * factor;
        Vec3::new(x,y,z)
    }
}
impl Add for Vec3 {
    type Output = Vec3;
    fn add(self, other: Vec3) -> Vec3 {
        let x = self.x + other.x;
        let y = self.y + other.y;
        let z = self.z + other.z;
        Vec3::new(x,y,z)
    }
}
impl Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, other: Vec3) -> Vec3 {
        let x = self.x - other.x;
        let y = self.y - other.y;
        let z = self.z - other.z;
        Vec3::new(x,y,z)
    }
}
impl Neg for Vec3 {
    type Output = Vec3;
    fn neg(self) -> Vec3 {
        let x = -self.x;
        let y = -self.y;
        let z = -self.z;
        Vec3::new(x,y,z)
    }
}
impl From<Color> for Vec3 {
    fn from(other: Color) -> Self {
        let r = ((other.r-128)as f64)/ 128.0;
        let g = ((other.g-128)as f64)/ 128.0;
        let b = ((other.b-128)as f64)/ 128.0;
        Vec3::new(r,g,b)
    }
}
impl PartialOrd<Vec3> for Vec3 {
    fn partial_cmp(&self, _: &Vec3) -> Option<Ordering> {
        None
    }
    fn lt(&self, other: &Vec3) -> bool {
        self.x < other.x && self.y < other.y && self.z < other.z
    }
    fn gt(&self, other: &Vec3) -> bool {
        self.x > other.x && self.y > other.y && self.z > other.z
    }
    fn ge(&self, other: &Vec3) -> bool {
        self.x >= other.x && self.y >= other.y && self.z >= other.z
    }
    fn le(&self, other: &Vec3) -> bool {
        self.x <= other.x && self.y <= other.y && self.z <= other.z
    }
}

impl Div<Vec3> for f64 {
    type Output = Vec3;
    fn div(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self/other.x,
            y: self/other.y,
            z: self/other.z,
        }
    }
}

#[derive(Clone,Copy,Debug)]
pub struct Matrix3 {
    pub a: Vec3,
    pub b: Vec3,
    pub c: Vec3,
}
impl Matrix3 {
    pub fn new( a: Vec3, b: Vec3, c: Vec3) -> Matrix3 {
        Matrix3 { a: a, b: b, c: c }
    }
}
impl Mul<Vec3> for Matrix3 {
    type Output = Vec3;
    fn mul(self, other: Vec3) -> Vec3 {
        let x = self.a * other.x;
        let y = self.b * other.y;
        let z = self.c * other.z;
        x + y + z
    }
}


#[derive(Clone,Debug)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}
impl Color {
    pub const WHITE: Color = Color { r: 255, g: 255, b: 255 };
    pub const BLACK: Color = Color { r: 0, g: 0, b: 0 };
    pub fn new(r: u8, g: u8, b: u8) -> Color
    { Color{ r: r, g: g, b: b } }
}
impl Add for Color {
    type Output = Color;
    fn add(self, other: Color) -> Color {
        let mut r = self.r as u16 + other.r as u16;
        let mut g = self.g as u16 + other.g as u16;
        let mut b = self.b as u16 + other.b as u16;
        if r > 255
        { r = 255; }
        if g > 255
        { g = 255; }
        if b > 255
        { b = 255; }
        Color::new(r as u8, g as u8, b as u8)
    }
}
impl Mul<f64> for Color {
    type Output = Color;
    fn mul(self, other: f64) -> Color {

        let mut r = self.r as f64 * other;
        let mut g = self.g as f64 * other;
        let mut b = self.b as f64 * other;
        if r > 255.0
        { r = 255.0; }
        if g > 255.0
        { g = 255.0; }
        if b > 255.0
        { b = 255.0; }
        Color::new(r as u8, g as u8, b as u8)
    }
}
impl From<Vec3> for Color {
    fn from(other: Vec3) -> Self {
        let r: u8 = (other.x * 128.0 + 128.0) as u8;
        let g: u8 = (other.y * 128.0 + 128.0) as u8;
        let b: u8 = (other.z * 128.0 + 128.0) as u8;
        Color::new(r,g,b)
    }
}

#[allow(non_camel_case_types)]
#[derive(PartialEq,PartialOrd)]
pub struct f64_ord(pub f64); //wrapper type for f64 since the default only impliments PartialOrd
impl Eq for f64_ord {}
impl Ord for f64_ord {
    fn cmp(&self, other: &Self) -> Ordering {
        let lhs = &self.0;
        let rhs = &other.0;
        match lhs.partial_cmp(rhs) {
            Some(ordering) => ordering,
            None => {
                if lhs.is_nan() {
                    if rhs.is_nan() {
                        Ordering::Equal
                    } else {
                        Ordering::Greater
                    }
                } else {
                    Ordering::Less
                }
            }
        }
    }
}
