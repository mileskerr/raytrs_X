/*
 * space impliments the really basic types I use everywhere
 * scn contains more advanced types for building scenes
 * i'll let you guess what renderer does
 *
*/
pub mod space;
mod mat;
mod scn;
mod renderer;

use crate::space::*;
use crate::scn::*;
use std::fs;
use std::io::BufWriter;
use std::time::Instant;

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

    let scene = Scene {
        mesh: mesh,
        camera: camera,
        mats: vec![
            Box::new(mat::Simple(Color::new(255,0,0))),
            Box::new(mat::Simple(Color::new(0,255,0)))
        ],
        lights: vec![
            Light::Point{ origin: Vec3::new(1.0,5.0,-2.0), intensity: 1.0 }
        ],
    };
    let data = renderer::render(scene, WIDTH, HEIGHT);

    let file = fs::File::create("render.png").unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, WIDTH as u32, HEIGHT as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().unwrap();

    writer.write_image_data(&data).unwrap();
    println!("done in {}s",t0.elapsed().as_secs_f32());
}





fn read_obj(contents: &str, offset: Vec3, scale: Vec3) -> Mesh { //TODO: handle errors and return result
    let mut verts: Vec<Vec3> = vec![];
    let mut norms: Vec<Vec3> = vec![];
    let mut tris: Vec<[usize;9]> = vec![];
    
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
            let mut vi = Vec::new();
            let mut ni = Vec::new();
            let mut ti = Vec::new();
            for value in &values[1..] {
                if value.is_empty() == false {
                    let ind: Vec<&str> = value.split('/').collect();
                    vi.push( ind[0].parse::<usize>().unwrap()-1 );
                    ni.push( ind[2].parse::<usize>().unwrap()-1 );
                    ti.push( ind[1].parse::<usize>().unwrap()-1 );
                }
            }
            tris.push([vi[0],vi[1],vi[2], ni[0],ni[1],ni[2], ti[0],ti[1],ti[2]]);
            if ni.len() > 3 { //quad
                tris.push([vi[0],vi[2],vi[3], ni[0],ni[2],ni[3], ti[0],ti[2],ti[3]]);
            }
        }
    }
    Mesh {
        mats: vec![((0..tris.len()),0)],
        verts: verts,
        norms: norms,
        tris: tris,
    }
}
