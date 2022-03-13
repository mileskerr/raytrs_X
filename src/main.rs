pub mod space;
mod extras;
mod renderer;
use crate::renderer::*;
use crate::space::*;
use crate::extras::*;
use std::fs;
use std::path::Path;
use std::io::BufWriter;
use std::time::{Duration,Instant};

fn main() {
    render();
}

const OBJ_PATH: &str = "utah_teapot.obj";
const WIDTH: usize = 256;
const HEIGHT: usize = 256;

fn render() {
    let t0 = Instant::now();
    let camera = Camera::new(Vec3::new(0.0,3.0,-10.0),Vec3::new(0.0,0.0,1.0),1.0);
    let mesh = read_obj(&fs::read_to_string(OBJ_PATH).unwrap(), Vec3::ZERO, Vec3::ONE);

    let data = renderer::render(mesh, camera, WIDTH, HEIGHT);

    let file = fs::File::create("render.png").unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, WIDTH as u32, HEIGHT as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().unwrap();

    writer.write_image_data(&data).unwrap();
    println!("done in {}s",t0.elapsed().as_secs_f32());
}







fn read_obj(contents: &str, offset: Vec3, scale: Vec3) -> Mesh {
    let mut verts: Vec<Vec3> = vec![];
    let mut norms: Vec<Vec3> = vec![];
    let mut tris: Vec<[usize;3]> = vec![];
    
    for line in contents.lines() {
        let is_vert = line.find("v ");
        if is_vert.is_some() { 
            let values: Vec<&str> = line.split(' ').collect();
            verts.push(Vec3 {
                x: values[1].parse::<f64>().unwrap(),
                y: values[2].parse::<f64>().unwrap(),
                z: values[3].parse::<f64>().unwrap(),
            } * scale + offset);
        }
        let is_norm = line.find("vn ");
        if is_norm.is_some() { 
            let values: Vec<&str> = line.split(' ').collect();
            norms.push(Vec3 {
                x: values[1].parse::<f64>().unwrap(),
                y: values[2].parse::<f64>().unwrap(),
                z: values[3].parse::<f64>().unwrap(),
            });
        }
        
        let is_face = line.find("f ");
        if is_face.is_some() { 
            let values: Vec<&str> = line.split(' ').collect();
            let mut i = Vec::new();
            for value in &values[1..] {
                if value.is_empty() == false {
                    let ind: Vec<&str> = value.split('/').collect();
                    i.push( ind[0].parse::<usize>().unwrap()-1 );
                }
            }
            tris.push([i[0], i[1], i[2]]);
            if i.len() > 3 { //quad
                tris.push([i[0], i[2], i[3]]);
            }
        }
    }
    Mesh {
        verts: verts,
        norms: norms,
        tris: tris,
    }
}
