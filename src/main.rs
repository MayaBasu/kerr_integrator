mod trajectory_graph;
mod functions;
mod numeric_solvers;
mod constants;
use crate::numeric_solvers::root_hunt_and_peck;
use crate::functions::*;
use crate::constants::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use csv::Writer;

mod derivatives;

use roots::{Roots, find_roots_quartic};
use egui;

use std::io;
use serde::Serialize;

fn main() -> Result<(), Box<dyn Error>>{

    println!("Please select a functionality: show radial derivative (1), theta derivative (2).");
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read line");
    let input: u32 = input.trim().parse().expect("Please type 1 or 2");

    match input {
        1 => {
            println!("Setting up the radial derivative ...");
            let radialgraph = derivatives::RadialGraph::new();

            let j = serde_json::to_string(&radialgraph)?;

            // Print, write to a file, or send to an HTTP server.
            println!("{}", j);

           // let mut wtr = Writer::from_path("foo.csv")?;
         //   wtr.serialize(radialgraph)?;
            let mut file = File::create("output.json")?;
            file.write_all(j.as_bytes())?;
            Ok(())

        }

        2 => {
            println!("Setting up the theta derivative ...");
            let thetagraph = derivatives::ThetaGraph::new();

            let j = serde_json::to_string(&thetagraph)?;

            // Print, write to a file, or send to an HTTP server.
           // println!("{}", j);

            // let mut wtr = Writer::from_path("foo.csv")?;
            //   wtr.serialize(radialgraph)?;
            let mut file = File::create("output.json")?;
            file.write_all(j.as_bytes())?;
            Ok(())

        }

            _=> {println!("Please select 1 or 2"); panic!();}

    }






}



fn setup(r_coefficients:[f64;5], theta_coefficients:[f64;5], simple: bool) -> ((f64, f64, f64), (f64, f64, f64)) { //((r start, rmin, rmax),(theta start,theta min,theta max))


    let L_sqr = C + LZ.powi(2);

    let beta = (A.powi(2))*(1.0 - E.powi(2));

    let (z_plus, z_minus) = if A != 0.0{
        (((L_sqr + beta) + ((L_sqr + beta).powi(2) - 4.0*C*beta).sqrt())/(2.0*beta),
        ((L_sqr + beta) - ((L_sqr + beta).powi(2) - 4.0*C*beta).sqrt())/(2.0*beta))

    } else {
        (0.0,0.0)
    };


   println!("z plus and minus {}   {}",z_plus,z_minus);



    let multiple_roots = find_roots_quartic(r_coefficients[4], r_coefficients[3], r_coefficients[2], r_coefficients[1], r_coefficients[0]);
    let mut root_list: Vec<f64> = Vec::new();

    match multiple_roots {
        Roots::Four(roots) =>{
            for root in roots{
                if root >=0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::Three(roots) =>{
            for root in roots{
                if root >=0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::Two(roots) =>{
            for root in roots{
                if root >=0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::One(roots) =>{
            for root in roots{
                if root >=0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::No(roots)=>{

            for root in roots{
                if root >0.0{
                    root_list.push(root);

                }
            }
        }
    }
    println!("The positive roots of the radial derivative function are at {:?}",root_list);
    if root_list.len() < 2{
        panic!("There were less than 2 roots :( ")
    }

    for _i in 0..root_list.len()-1 {
       // println!("sample intermediate value: {}",coefficients_to_poly((root_list[i+1]-root_list[i])/2.0+root_list[i],r_coefficients));
    }
    let (r_start,r_min,r_max)  = if simple {
        println!("Implementing simple setup...");
        let r_start = root_list[root_list.len() -2];
        let r_min = root_list[root_list.len()-2];
        let r_max = root_list[root_list.len()-1];
        println!("Chose r initial {}, r min {} r max {}",r_start,r_min,r_max);
        (r_start,r_min,r_max)

    }else {
        let mut r_start = String::new();
        println!("Select an initial starting value for r:");

        io::stdin()
            .read_line(&mut r_start)
            .expect("Failed to read line");
        let r_start = r_start.trim().parse().expect("Please type a number!");
        let mut r_min = 0.0;
        let mut r_max = 0.0;
        println!("roots {:?}", root_list);
        for root in root_list {
            if root <= r_start {
                r_min = root;
            }
            if root > r_start {
                r_max = root;
                break
            }
        };
        (r_start,r_min,r_max)
    };
    let (theta_start,theta_min,theta_max) = if simple {

        let theta_start =  std::f64::consts::PI/2.0;//1.0471975511965979;
        let graph = |x: f64| -> f64 {
            coefficients_to_poly(x.cos(), theta_coefficients)
        };

        let (theta_min,theta_max) = if A != 0.0{root_hunt_and_peck(theta_start,graph)}else {(0.0,0.0)};
        println!("selecting theta initial{} theta  min {} and theta max {}",theta_start,theta_min,theta_max);

        (theta_start,theta_min,theta_max)

    } else{
        let mut theta_start =COS_I.acos();
        println!("Starting with theta = {}",theta_start);


        let graph = |x: f64| -> f64 {
            coefficients_to_poly(x.cos(), theta_coefficients)
        };
        let (theta_min,theta_max) = if A != 0.0{root_hunt_and_peck(theta_start,graph)}else {(0.0,0.0)};

        println!("Then the bounds are {} to {} for theta",theta_min,theta_max);
        (theta_start,theta_min,theta_max)

    };


    ((r_start,r_min,r_max),(theta_start,theta_min,theta_max))


}



