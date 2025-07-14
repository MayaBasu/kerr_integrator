use crate::structs::{GeodesicGraph, StellarParams,StarChunk};


pub fn initialize_star_chunks(stellar_params: StellarParams, r_cm:f64, r_stellar:f64, num_slices:usize, theta_initial:f64) -> Vec<StarChunk>{
    //find the binding energy of slices of a star which is disrupted at a certain radius with the freezing approximation
    //in units of G = M = c = 1, the stellar radius is 0.47
    // so taylor expanding the newtonian potential we have delta e = GMr/r_cm^2 where r is how far out and r_cm is the radiout out fromhe center of mass
    //first most basic model, is that of an even density sphere.
    let mut chunk_data = Vec::with_capacity(num_slices);
    let delta_r = r_stellar/num_slices as f64;
    for i in 0..num_slices {
        let r_out = i as f64 * delta_r;
        let delta_binding_energy= r_out / (r_cm.powi(2));
        let binding_energy = stellar_params.e - delta_binding_energy;
        let fraction_of_star= fractional_stellar_slice_volume(r_out, delta_r, r_stellar);
        let chunk_stellar_parameters = StellarParams::new(stellar_params.lz, binding_energy, stellar_params.c);
        let chunk_geodesic_graph = GeodesicGraph::new(chunk_stellar_parameters,r_cm,theta_initial);

        chunk_data.push(
            StarChunk {
                geodesic_graph: chunk_geodesic_graph,
                fraction_of_star,
                binding_energy
            })
    }
    chunk_data
}

pub fn fractional_stellar_slice_volume(r_out:f64,delta_r:f64,r_star:f64) ->f64{
    let stellar_chunk_radius = (r_star.powi(2)-r_out.powi(2)).sqrt();
    let stellar_volume = (4.0/3.0)* std::f64::consts::PI *r_star.powi(3);
    let stellar_chunk_volume = std::f64::consts::PI*stellar_chunk_radius.powi(2)*delta_r;
    stellar_chunk_volume/stellar_volume
}


