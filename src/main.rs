use crate::derivatives::{r_derivative_propertime, theta_derivative};
use crate::functions::find_theta_parameters;
use crate::star_initialization::tidal_radius;
use crate::structs::{GeodesicGraph, Star, StellarParams};
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
const NUM_STEPS: usize = 40000;
const num_phi_bins: usize=100;
fn main() {

    let e = 0.9999;
    let lz = 0.5*6.5;
    let c = 6.5*6.5*(1.0-0.5*0.5);
    let stellar_radius = 0.47;
    let stellar_mass = 10.0_f64.powi(-6);
    let tidal_radius = tidal_radius(stellar_radius,stellar_mass);
    let stellar_params:StellarParams = StellarParams {lz,e,c};


/*
    let mut graph = GeodesicGraph::new(stellar_params, 9950.0, 1.5235996057125137,NUM_STEPS,STEP_SIZE); //theta_params.z_minus.sqrt().acos()+0.001
    graph.calculate_stream_width(0.7,0.0);
    graph.clone().find_intersections();
    graph.serialize("ballistic_graph.json").unwrap()

 */

    println!("The stellar mass is {:?}tidal radius is at {:?}",stellar_mass,tidal_radius);









    let star = Star::new(stellar_params, tidal_radius, 800.0,stellar_radius,NUM_STEPS,STEP_SIZE,1);

    star.serialize("star.json").unwrap()


}



