use egui_plot::{Plot,PlotPoints,Line};
use crate::constants::{E, LZ, C};
use crate::functions::{r_derivative_magnitude, radial_coefficients};
use std::io;

use serde::{Deserialize, Serialize};
use std::error::Error;
use csv::Writer;


const RADIAL_DERIVATIVE_GRANULATION: u32 = 100;

use serde_json::Result;

#[derive(Serialize, Deserialize)]
#[derive(Default)]
pub(crate) struct RadialGraph {

    upper_bound:u32,

    pub(crate) data:Vec<[f64;2]>,

}




impl RadialGraph {



    pub(crate) fn new() -> Self {

        let mut upper_bound = String::new();
        println!("Upper Bound of Graph: (u32)");
        io::stdin()
            .read_line(&mut upper_bound)
            .expect("Failed to read line");
        let upper_bound: u32 = upper_bound.trim().parse().expect("Please type a u32");
        let coefficients = radial_coefficients(LZ,E,C);
        println!("The radial coefficients are");
        println!("{:?}",coefficients);
        let data:Vec<[f64;2]> = (0..upper_bound*RADIAL_DERIVATIVE_GRANULATION).map(|i|
        [(i as f64)/(RADIAL_DERIVATIVE_GRANULATION as f64), r_derivative_magnitude((i as f64)/(RADIAL_DERIVATIVE_GRANULATION as f64), coefficients)]
        ).collect();

        Self {
            upper_bound,
            data,
        }
    }
}

