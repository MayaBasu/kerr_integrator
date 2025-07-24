use crate::derivatives::{r_derivative_propertime, theta_derivative};
use crate::functions::find_theta_parameters;
use crate::structs::{ StellarParams};
mod functions;
mod numeric_integrators;
mod derivatives;
mod structs;
mod tetrads;
mod tests;
mod star_initialization;


pub const M: f64 = 1.0;
pub const A: f64 = 0.9;


const STEP_SIZE: f64 = 0.0001;
const NUM_STEPS: usize = 400;

fn main() {

    let e = 0.9999;
    let lz = 0.5*6.5;
    let c = 6.5*6.5*(1.0-0.5*0.5);






    let stellar_params:StellarParams = StellarParams {lz,e,c};

    let r_deriv = theta_derivative(1000.0, 0.0,false,stellar_params);


    println!("{:?}",r_deriv);




/*

    let mut graph = GeodesicGraph::new(stellar_params, 19.265297560836153, 0.5235996057125137,NUM_STEPS,STEP_SIZE); //theta_params.z_minus.sqrt().acos()+0.001
    graph.calculate_stream_width(0.7,0.0);
    graph.serialize("maingraph.txt").unwrap();

    graph.find_self_intersections(1000);

 */






  //  let star = Star::new(stellar_params, 19.265297560836153, 0.5235996057125137, 2.35,NUM_STEPS,STEP_SIZE);

 //   star.serialize("star.json").unwrap()


}



