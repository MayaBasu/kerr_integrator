
mod functions;

mod numeric_integrators;
use std::fs::File;
use std::io::Write;
mod constants;

use crate::functions::*;
use crate::constants::*;
use crate::numeric_integrators::*;

use std::error::Error;

mod derivatives;
mod structs;

use serde::Serialize;
use crate::structs::Graph;
use serde_arrays;

fn main() -> Result<(), Box<dyn Error>>{


    let radialgraph = Graph::new(LZ,E,C);

    let j = serde_json::to_string(&radialgraph)?;

    let mut file = File::create("outpu2t.json")?;
    file.write_all(j.as_bytes())?;
    println!("{:?}",find_radial_parameters(LZ,E,C));
    println!("{}  {}   {}", E,LZ,C);
    Ok(())


        /*

    println!("Please select a functionality: show radial derivative (1), theta derivative (2).");
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read line");
    let input: u32 = input.trim().parse().expect("Please type 1 or 2");

    match input {

        3 => {



        }

            _=> {println!("Please select 1 or 2"); panic!();}

    }

         */






}



