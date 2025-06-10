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

            println!("{}", j);

            let mut file = File::create("output.json")?;
            file.write_all(j.as_bytes())?;
            Ok(())

        }

        2 => {
            println!("Setting up the theta derivative ...");
            let thetagraph = derivatives::ThetaGraph::new();
            let j = serde_json::to_string(&thetagraph)?;
            let mut file = File::create("output.json")?;
            file.write_all(j.as_bytes())?;
            Ok(())

        }

        3 => {

            parameter_converter(root_finder(radial_coefficients(LZ,E,C)));
            Ok(())

        }

            _=> {println!("Please select 1 or 2"); panic!();}

    }






}



