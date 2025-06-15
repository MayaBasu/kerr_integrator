
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


    let graph = Graph::new(LZ,E,C);
    let graph_json = serde_json::to_string(&graph)?;
    let mut file = File::create("outpu2t.json")?;
    file.write_all(graph_json.as_bytes())?;

    println!("{:?}",find_radial_parameters(LZ,E,C));
    println!("{}  {}   {}", E,LZ,C);
    Ok(())
}



