mod trajectory_graph;
mod functions;
mod numeric_solvers;
mod constants;
use crate::numeric_solvers::root_hunt_and_peck;
use crate::functions::*;
use crate::constants::*;

mod radial_derivative;

use roots::{Roots, find_roots_quartic};
use egui;
use plotters::prelude::*;
use std::io;


fn main() -> eframe::Result<()> {

    println!("Please select a functionality: show radial derivative (1), plot (2).");
    let mut func = String::new();
    io::stdin()
        .read_line(&mut func)
        .expect("Failed to read line");
    let func: u32 = func.trim().parse().expect("Please type 1 or 2");

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size((1600.0, 800.0)),
        ..eframe::NativeOptions::default()
    };

    match func {
        1 => {
            println!("Setting up the radial derivative ...");
            eframe::run_native(
                "Visualizer",
                native_options,
                Box::new(|cc| Ok(Box::new(radial_derivative::RadialGraph::new(cc)))),
            )


        }
        2=> {
            println!("Plotting the trajectory");
            eframe::run_native(
                "Visualizer",
                native_options,
                Box::new(|cc| Ok(Box::new(trajectory_graph::Graph::new(cc, setup(radial_coefficients(LZ, E, C), theta_coefficients(LZ, E, C), true))))),
            )
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

    for i in 0..root_list.len()-1 {
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
        let mut theta_start =String::new();
        println!("Select an initial starting value for theta:");

        io::stdin()
            .read_line(&mut theta_start)
            .expect("Failed to read line");
        let theta_start = theta_start.trim().parse().expect("Please type a number!");

        let graph = |x: f64| -> f64 {
            coefficients_to_poly(x.cos(), theta_coefficients)
        };
        let (theta_min,theta_max) = if A != 0.0{root_hunt_and_peck(theta_start,graph)}else {(0.0,0.0)};

        println!("Then the bounds are {} to {} for theta",theta_min,theta_max);
        (theta_start,theta_min,theta_max)

    };


    ((r_start,r_min,r_max),(theta_start,theta_min,theta_max))


}



