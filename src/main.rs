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
use mat::Material;

const OBJ_PATH: &str = "utah_teapot.obj";
const WIDTH: usize = 256;
const HEIGHT: usize = 256;

fn main() {
    let num_pixels = WIDTH * HEIGHT;
    let t0 = Instant::now();
    let camera = Camera::new(Vec3::new(0.0,3.0,-10.0),Vec3::new(0.0,0.0,1.0),1.0);
    let mut mesh1 = read_obj(&fs::read_to_string(OBJ_PATH).unwrap(), Vec3::ZERO, Vec3::ONE);
    let mesh2 = read_obj(&fs::read_to_string(OBJ_PATH).unwrap(), Vec3::new(-1.0,2.0,2.0), Vec3::ONE);
    let mesh3 = Mesh {
        tris: vec![
            Tri::new(
                Vec3::new(3.0,0.0,3.0),
                Vec3::new(3.0,0.0,-3.0),
                Vec3::new(-3.0,0.0,3.0),
                0,
            ),
            Tri::new(
                Vec3::new(3.0,0.0,-3.0),
                Vec3::new(-3.0,0.0,-3.0),
                Vec3::new(-3.0,0.0,3.0),
                0,
            )
        ],
        spheres: vec![
			Sphere::new(
				Vec3::new(1.0,4.0,-2.0),
				0.75,
                1,
			)
		],
    };
    mesh1.join(mesh2);
    mesh1.join(mesh3);

    let scene: Scene = Scene {
        mesh: mesh1,
        camera: camera,
        mats: vec![
            Box::new(mat::Lambertian::new(Color::WHITE)),
            Box::new(mat::Emissive::new(Color::WHITE,10.0)),
        ],
    };
    let pixels = renderer::render(scene, WIDTH, HEIGHT, 2000, 16);
    let mut data = Vec::new();
    for i in 0..num_pixels {
        let mut sur = Vec::new();
        sur.push((i, 25.0));
        if i > WIDTH+1 {
            sur.push((i-WIDTH, 4.0));
            sur.push((i-WIDTH+1, 1.0));
            sur.push((i-WIDTH-1, 1.0));
        } if i > (WIDTH * 2) {
            sur.push((i-WIDTH * 2, 0.5));
        } if i+WIDTH+1 < num_pixels {
            sur.push((i+WIDTH, 4.0));
            sur.push((i+WIDTH+1, 1.0));
            sur.push((i+WIDTH-1, 1.0));
        } if i+(2 * WIDTH) < num_pixels {
            sur.push((i+(2 * WIDTH), 0.5));
        } if i > 0 {
            sur.push((i-1, 4.0));
        } if i+1<num_pixels {
            sur.push((i+1, 4.0));
        }
        let mut pixel = Vec3::ZERO;
        let mut total = 0.0;
        for s in sur {
            let diff = (0.2/((pixels.1[i] - pixels.1[s.0]).abs()+0.01)-2.0).clamp(0.0,1.0);
            let hdr: Vec3 = pixels.0[s.0].into();
            let weight = diff * s.1;
            pixel = pixel + (hdr * weight);
            total += weight;
        }
        let c: Color = (pixel / total).into();
        data.push(c.r);
        data.push(c.g);
        data.push(c.b);
    }
    let file = fs::File::create("render.png").unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, WIDTH as u32, HEIGHT as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().unwrap();

    writer.write_image_data(&data).unwrap();

    let mut data = Vec::new();
    for pixel in pixels.0 {
        data.push(pixel.r);
        data.push(pixel.g);
        data.push(pixel.b);
    }
    let file = fs::File::create("noise.png").unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, WIDTH as u32, HEIGHT as u32);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().unwrap();

    writer.write_image_data(&data).unwrap();
    println!("done in {}s",t0.elapsed().as_secs_f32());


}





fn read_obj(contents: &str, offset: Vec3, scale: Vec3) -> Mesh { //TODO: handle errors and return result
    let mut tris: Vec<Tri> = vec![];
    let mut norms: Vec<Vec3> = vec![];
    let mut verts: Vec<Vec3> = vec![];
    
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
            tris.push(Tri {
                    verts: [verts[vi[0]],verts[vi[1]],verts[vi[2]]], 
                    norms: [norms[ni[0]],norms[ni[1]],norms[ni[2]]],
                    mat: 0,
            });
            if ni.len() > 3 { //quad
                tris.push(Tri {
                        verts: [verts[vi[0]],verts[vi[2]],verts[vi[3]]], 
                        norms: [norms[ni[0]],norms[ni[2]],norms[ni[3]]],
                        mat: 0,
                });
            }
        }
    }
    Mesh {
        tris: tris,
        spheres: Vec::new(),
    }
}
