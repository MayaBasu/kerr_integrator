
use rayon::prelude::*;
use crate::functions::radial_roots;
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

pub const NUM_PHI_BINS:usize =  1000;



fn main() {

    let step_size =  0.0001;
    let num_steps = 80000;
    let num_slices = 100;

    let cos_inclination = 0.5;
    let l_mag:f64 = 6.5_f64;
    let e = 0.9999;
    let stellar_radius = 0.47;
    let stellar_mass = 10.0_f64.powi(-6);


    let lz = cos_inclination*l_mag;
    let c = l_mag.powi(2)*(1.0-cos_inclination.powi(2));
    let stellar_params:StellarParams = StellarParams {lz,e,c};
    let tidal_radius = tidal_radius(stellar_radius,stellar_mass);
    let (periapsis,apoapsis) = (radial_roots(stellar_params)[2],radial_roots(stellar_params)[3]);
    //check that the star is disrupted within it's orbit

    //calculate_total_angular_momentum(periapsis,apoapsis,stellar_params,num_steps,step_size,tidal_radius,stellar_radius,num_slices)
    let mut graph = GeodesicGraph::new(stellar_params, 9000.0, 1.5, num_steps, step_size);
    graph.calculate_stream_width(0.7,0.0);
    graph.serialize("src/maingraph.txt").unwrap()
}


fn collect_data_for_fig_oneandtwo(stellar_params: StellarParams,num_steps:usize, step_size:f64,tidal_radius:f64){
    println!("the tidal radius is {tidal_radius}");

    let graph = GeodesicGraph::new(stellar_params, tidal_radius, std::f64::consts::PI/2.0, num_steps, step_size);

    graph.serialize("maingraph.txt").unwrap()
}
fn calculate_total_angular_momentum(periapsis:f64,apoapsis:f64,stellar_params: StellarParams,num_steps:usize,step_size:f64,tidal_radius:f64,stellar_radius:f64,num_slices:usize){
    println!("The tidal radius is {tidal_radius}, the periapsis is {periapsis} apoapsis is {apoapsis}");
    if periapsis > tidal_radius{
        println!("No Disruption");
        panic!()
    }

    let mut graph = GeodesicGraph::new(stellar_params, 9000.0, 1.5,num_steps,step_size);
    graph.calculate_stream_width(0.7,0.0);
    let lag = graph.lag(tidal_radius,9000.0,1.5).2;
    let stop_time = graph.find_stop_time()+lag;

    println!("Stop time is {stop_time} The lag is {:?}",lag);
    println!("The first self-intersection time is {:?}, the stream was calculated out to {:?}",stop_time,graph.t_graph.last());
    let stop_index = graph.return_phi_at_t(stop_time).0;
    graph = graph.trim(stop_index);
    graph.serialize("maingraph.txt").unwrap();


    let step_size =  0.01;
    let num_steps = 1000000;

    let mut star = Star::new_trimmed(stellar_params, 47.0, 1.5,stellar_radius,num_steps,step_size,num_slices,stop_time);

    println!("Calculating the total angular momentum as a function of time...");
    let angular_momentums = star.find_angular_momentum_over_time(stop_time,10);
    println!("Angular momentums: {:?}",angular_momentums);

    star.serialize("src/stars/star_1.0_6.5_6.txt").unwrap();

}



