
use rayon::prelude::*;
use serde::Serialize;
use crate::functions::radial_roots;
use crate::star_initialization::tidal_radius;
use crate::structs::{GeodesicGraph, Star, StellarParams};
mod functions;
mod numeric_integrators;
mod derivatives;
mod structs;
mod star_initialization;


pub const M: f64 = 1.0;
pub const A: f64 = 0.9;

pub const NUM_PHI_BINS:usize =  1000;



fn main() {
    let cos_inclinations = [0.50];

    for cos_inclination in cos_inclinations.into_iter(){
        let step_size =  0.00001;
        let num_steps = 500000;
        let num_slices = 100;

      //  let cos_inclination = 0.49;
        let l_mag:f64 = 6.50_f64;
        let e = 0.9999;
        let stellar_radius = 0.47;
        let stellar_mass = 10.0_f64.powi(-6);


        let lz = cos_inclination*l_mag;
        let c = l_mag.powi(2)*(1.0-cos_inclination.powi(2));
        let stellar_params:StellarParams = StellarParams {lz,e,c};
        let tidal_radius = tidal_radius(stellar_radius,stellar_mass);
        let (periapsis,apoapsis) = (radial_roots(stellar_params)[2],radial_roots(stellar_params)[3]);
        //check that the star is disrupted within it's orbit
        let label= make_label(cos_inclination,l_mag,stellar_mass.log(10.0));
        //calculate_total_angular_momentum_alt(periapsis,apoapsis,stellar_params,num_steps,step_size,tidal_radius,stellar_radius,num_slices,label)
        // let mut star = Star::new_trimmed(stellar_params, 47.0, 1.5, stellar_radius,  num_steps, step_size, 100,5000000.0);

       //  calculate_single_geodesic(StellarParams{lz:stellar_params.lz,c:stellar_params.c,e:0.9999},num_steps,step_size,true);
        calculate_single_geodesic(StellarParams{lz:stellar_params.lz,c:stellar_params.c,e:0.9996894},num_steps,step_size,true);

        //star.serialize("star.json").unwrap
    }






}

fn make_label(cos_I:f64, l:f64, exp:f64) -> String {

    let mut label = String::new();
    label.push_str(&format!("src/stars/star_demo_{}_{}_{}.txt", cos_I, l, exp.round().abs()));

    println!("Made the label {label} from the parameters cos_I = {cos_I}, l = {l}, exp = {exp}");

    label
}


fn collect_data_for_fig_oneandtwo(stellar_params: StellarParams,num_steps:usize, step_size:f64,tidal_radius:f64){
    println!("the tidal radius is {tidal_radius}");

    let graph = GeodesicGraph::new(stellar_params, tidal_radius, std::f64::consts::PI/2.0, num_steps, step_size);

    graph.serialize("maingraph.txt").unwrap()
}
fn calculate_total_angular_momentum(periapsis:f64,apoapsis:f64,stellar_params: StellarParams,num_steps:usize,step_size:f64,tidal_radius:f64,stellar_radius:f64,num_slices:usize,path_label:String){
    println!("The tidal radius is {tidal_radius}, the periapsis is {periapsis} apoapsis is {apoapsis}");
    if periapsis > tidal_radius{
        println!("No Disruption");
        panic!()
    }

    let mut graph = GeodesicGraph::new(stellar_params, 9900.0, 1.57,num_steps,step_size);
    graph.calculate_stream_width(0.7,0.0);
    let lag = graph.lag(tidal_radius,9000.0,1.57).2;
    let stop_time = graph.find_stop_time()+lag;

    println!("Stop time is {stop_time} The lag is {:?}",lag);
    println!("The first self-intersection time is {:?}, the stream was calculated out to {:?}",stop_time,graph.t_graph.last());
    let stop_index = graph.return_phi_at_t(stop_time).0;
    graph = graph.trim(stop_index);
    graph.serialize("maingraph.txt").unwrap();


    let step_size =  0.0001;
    let num_steps = 200000;

    let mut star = Star::new_trimmed(stellar_params, 47.0, 1.5,stellar_radius,num_steps,step_size,num_slices,stop_time);

    println!("Calculating the total angular momentum as a function of time...");
    let angular_momentums = star.find_angular_momentum_over_time(stop_time,10);
    println!("Angular momentums: {:?}",angular_momentums);

    star.serialize(&*path_label).unwrap();

}

fn calculate_total_angular_momentum_alt(periapsis:f64,apoapsis:f64,stellar_params: StellarParams,num_steps:usize,step_size:f64,tidal_radius:f64,stellar_radius:f64,num_slices:usize,path_label:String){
    println!("The tidal radius is {tidal_radius}, the periapsis is {periapsis} apoapsis is {apoapsis}");
    if periapsis > tidal_radius{
        println!("No Disruption");
        panic!()
    }


    let mut graph = GeodesicGraph::new(StellarParams{lz:stellar_params.lz,c:stellar_params.c,e:0.9996894}, 9900.0, 1.57,num_steps,step_size);

    println!("calculated tiny stream");
    graph.calculate_stream_width(0.07,0.0);
    println!("calculated tiny stream width");
    let lag = graph.lag(tidal_radius,9000.0,1.57).2;
    let stop_time = graph.find_stop_time()+lag;

    println!("Stop time is {stop_time} The lag is {:?}",lag);
    println!("The first self-intersection time is {:?}, the stream was calculated out to {:?}",stop_time,graph.t_graph.last());
    let stop_index = graph.return_phi_at_t(stop_time).0;
    graph = graph.trim(stop_index);
    graph.serialize("maingraph.txt").unwrap();


    let step_size =  0.0001;
    let num_steps = 200000;

    let mut star = Star::new_trimmed(stellar_params, 47.0, 1.5,stellar_radius,num_steps,step_size,num_slices,stop_time);

    println!("Calculating the total angular momentum as a function of time...");
    let angular_momentums = star.find_angular_momentum_over_time(stop_time,10);
    println!("Angular momentums: {:?}",angular_momentums);

    star.serialize(&*path_label).unwrap();

}


fn calculate_single_geodesic(stellar_params: StellarParams,num_steps:usize,step_size:f64,find_point:bool){
    let mut geodesic = GeodesicGraph::new(stellar_params,47.0,1.57,num_steps,step_size);
    if find_point{
        geodesic.calculate_stream_width(0.07,0.0);
    }


  //  geodesic.calculate_stream_width(0.7,0.0);
    //let lag = geodesic.lag(tidal_radius,9000.0,1.5).2;
  //  let stop_time = geodesic.find_stop_time();//+lag
    //let stop_index = geodesic.return_phi_at_t(stop_time).0;
   // geodesic = geodesic.trim(stop_index);

   // geodesic.wrapped_data = geodesic.wrap_data().to_vec();
    //  println!("geodesidc wrapped data is {:?}",geodesic.wrapped_data);
   // println!("Stop time  is {stop_time}");
    geodesic.serialize("tiny.txt").unwrap();
}