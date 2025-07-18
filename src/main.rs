use crate::functions::find_theta_parameters;
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

fn main() {
    let e = 0.9999;
    let lz = 0.5*6.5;
    let c = 6.5*6.5*(1.0-0.5*0.5);

    let stellar_params:StellarParams = StellarParams::new(lz,e,c);
    let theta_params = find_theta_parameters(stellar_params);




    let mut graph = GeodesicGraph::new(stellar_params, 19.265297560836153, 0.5235996057125137); //theta_params.z_minus.sqrt().acos()+0.001
    graph.calculate_stream_width(0.7,0.0);
    graph.serialize("maingraph.txt").unwrap();




    let star = Star::new(stellar_params, 19.265297560836153, 0.5235996057125137, 2.35);

    star.serialize("star.json").unwrap()


}



