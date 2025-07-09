
mod functions;

mod numeric_integrators;
use std::fs::File;
use std::io::Write;
mod constants;

use crate::functions::*;
use crate::tests::test_tetrads;
use crate::constants::*;
use crate::numeric_integrators::*;

use std::error::Error;
use crate::derivatives::{r_derivative_propertime, theta_derivative};
use crate::tetrads::*;
mod derivatives;
mod structs;
mod tetrads;
mod tests;
mod stream_intersection;

use crate::structs::{Graph, Stream};

fn main() -> Result<(), Box<dyn Error>>{

    println!("{:?}",find_radial_parameters(LZ,E,C));
    println!("E: {}  \nL_z: {}  \nC: {} \nL: {} \ncos I: {}", E,LZ,C,(C+LZ*LZ).sqrt(), C.sqrt()/((C+LZ*LZ).sqrt()));


    let graph = Graph::new(LZ,E,C);
    let stream = Stream{
        h:integrate_H(&graph,0.0,0.90).to_vec() //2.35
    };







    let graph_json = serde_json::to_string(&graph)?;
    let mut file = File::create("ballistic_graph.json")?;
    file.write_all(graph_json.as_bytes())?;


    let graph_stream = serde_json::to_string(&stream)?;
    let mut file = File::create("stream_width.json")?;
    file.write_all(graph_stream.as_bytes())?;



    test_tetrads(graph,900);





    Ok(())
}



