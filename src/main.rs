
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
use crate::structs::{Graph, Stream};

fn main() -> Result<(), Box<dyn Error>>{
    println!("{:?}",find_radial_parameters(LZ,E,C));
    println!("E: {}  \n L_z: {}  \n C: {} \n L: {} \n cos I: {}", E,LZ,C,(C+LZ*LZ).sqrt(), C.sqrt()/((C+LZ*LZ).sqrt()));


    let graph = Graph::new(LZ,E,C);
  //  println!("{:?}",graph.phi);



    let stream = Stream{
        h:integrate_H(&graph,0.0,0.92).to_vec() //2.35
    };

    let graph_json = serde_json::to_string(&graph)?;
    let mut file = File::create("ballistic_graph.json")?;
    file.write_all(graph_json.as_bytes())?;


    let graph_stream = serde_json::to_string(&stream)?;
    let mut file = File::create("stream_width.json")?;
    file.write_all(graph_stream.as_bytes())?;





    Ok(())
}



