
mod functions;
mod numeric_integrators;
mod constants;

use crate::functions::*;
use crate::constants::*;
use crate::numeric_integrators::*;
use std::error::Error;

mod derivatives;



use serde::Serialize;

fn main() -> Result<(), Box<dyn Error>>{

    let roots = root_finder(radial_coefficients(LZ,E,C));
    println!("roots are{:?}",roots);
    let (p,e,p3,p4) = parameter_converter(roots);
    println!("parameters:  {:?}",(p,e,p3,p4));
    let array = integrate_r(100.0,p,e,p3,p4);

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



