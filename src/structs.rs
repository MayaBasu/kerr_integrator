use serde::{Serialize};
use std::fs::File;
use std::io::Write;
use std::error::Error;
use serde_json::Value::Array;
use crate::functions::{delta, distance, lower_distance_bound, radial_roots, sigma};
use crate::numeric_integrators::{integrate_H, integrate_phi, integrate_r, integrate_t, integrate_theta};
use crate::star_initialization::initialize_star_chunks;

#[derive(Clone, Copy, Debug)]
pub struct RadialParams{
    pub p:f64,
    pub e:f64,
    pub p3:f64,
    pub p4:f64
}
#[derive(Clone, Copy)]
pub struct ThetaParams{
    pub beta:f64,
    pub z_plus:f64,
    pub z_minus:f64,
}
#[derive(Debug, Clone, Copy,Serialize)]


pub struct StellarParams{
    pub(crate) lz:f64,
    pub(crate) c:f64,
    pub(crate) e:f64
}


pub struct WrappedDataPoint {
    pub index:usize,
    pub wrap:usize,
    pub phi:f64,
    pub r:f64,
    pub theta:f64,
    pub h:f64

}

#[derive(Serialize,Clone)]
pub struct GeodesicGraph {
    pub(crate) stellar_params: StellarParams,

    pub(crate) num_steps: usize,
    pub(crate) step_size: f64,

    pub(crate) t_graph: Vec<f64>,
    pub(crate) radial_graph: Vec<f64>,
    pub(crate) theta_graph: Vec<f64>,
    pub(crate) phi_graph: Vec<f64>,

    pub(crate) stream_height: Vec<f64>,

    pub(crate) intersections: Vec<(usize,usize,bool)>,
    pub(crate) bl_t_of_first_intersection: (f64,f64),
}

impl GeodesicGraph {
    pub fn new(stellar_params: StellarParams, r_initial: f64, theta_initial: f64, num_steps:usize,step_size:f64) -> Self {
        println!("Initializing a new Graph with: \n\
         (e,lz,c) = {:?} \n\
         (r_0,theta_0) = ({},{}) \n\
         It has radial roots at {:?}\
         \n Calculating trajectory for {num_steps} in Mino time spaced at {step_size}",
                 stellar_params, r_initial, theta_initial, radial_roots(stellar_params));

        let radial_graph = integrate_r(r_initial, stellar_params,num_steps,step_size).0;
        let theta_graph = integrate_theta(theta_initial,stellar_params,num_steps,step_size).0;
        let t_graph = integrate_t(radial_graph.clone(), theta_graph.clone(), stellar_params,num_steps,step_size);
        let phi_graph = integrate_phi(radial_graph.clone(),theta_graph.clone(), stellar_params,num_steps,step_size);
        println!(" phi values are {:?}",phi_graph);

        Self {
            stellar_params,

            num_steps,
            step_size,

            radial_graph,
            theta_graph,
            phi_graph,
            t_graph,

            stream_height: Vec::new(),
            intersections: Vec::new(),

            bl_t_of_first_intersection: (0.0,0.0),
        }
    }
    pub fn serialize(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let graph_json = serde_json::to_string(&self)?;
        let mut file = File::create(file_path)?;
        file.write_all(graph_json.as_bytes())?;
        Ok(())
    }
    pub fn calculate_stream_width(&mut self, h_initial: f64, h_dot_initial: f64) {
        self.stream_height = integrate_H(&self, h_dot_initial, h_initial, self.stellar_params,self.num_steps,self.step_size);
    }

    pub fn return_phi_at_t(self, global_time: f64) -> (usize,f64){
        assert!(self.t_graph.last().unwrap() > &global_time);

        let mut past_difference = (self.t_graph[0] - global_time).abs();

        for i in 1..self.t_graph.len() {
            let current_difference = (self.t_graph[i] - global_time).abs();

            if past_difference < current_difference {
                let best_index = i - 1;
                return (best_index,self.phi_graph[best_index])
            } else {
                past_difference = current_difference;
            }
        }
        panic!("Didn't find a closest match")
    }

    pub fn wrap_data(&self,num_phi_bins:usize)->Array<Vec<WrappedDataPoint>>{
        let mut wrapped_data = Array::new();
        for (index, phi) in self.phi_graph.iter().enumerate(){
            let remainder = phi % (2.0*std::f64::consts::PI);
            let wrap = ((phi-remainder)/(2.0*std::f64::consts::PI)).round() as i32;
            let bin_number = (remainder/num_phi_bins).floor() as i32;
            wrapped_data[bin_number].push(WrappedDataPoint{
                index,
                wrap,
                phi,
                r:self.radial_graph(index),
                theta:self.theta_graph(index),
                h:self.stream_height(index)})
        }
        wrapped_data
    }

    fn find_intersections(mut self,num_phi_bins:usize) -> Vec<((usize,usize), bool, f64,f64,bool)> {
        let wrapped_data = self.wrap_data(num_phi_bins);
        let mut intersections = Vec::new();
        for phi_bin in wrapped_data{
            if wrapped_data[phi_bin].len() < 2{continue} //check that there are at least 2 points at this phi value
            else{
                let sorted_by_wrap = wrapped_data[phi_bin].sort_by(|x| x.wrap.abs());
                if sorted_by_wrap.popfront().wrap == sorted_by_wrap.pop().wrap{continue} //check that there are points from at least two different wraps

                //find the smallest sigma and largest delta
                let sigma_min = wrapped_data[phi_bin].into_iter().map(|point| sigma(point.r,point.theta)).sort().popfont();
                let delta_max = wrapped_data[phi_bin].into_iter().map(|point| delta(point.r)).sort().pop();


                //now we want to find the maximum distance delta r and delta theta between these points
                let min_r_distances_and_corresponding_h_ave =  wrapped_data[phi_bin].into_iter().map(|point|
                                 wrapped_data[phi_bin]
                                     .into_iter()
                                     .map(|second_point|
                                     if second_point.wrap == point.wrap{}
                                     else{
                                         ((point.r-second_point.r).abs(),(point.h+second_point.h)/2.0,(point.index,second_point.index))
                                     }
                                     ).sort_by(|(r_distance,h,indices)| r_distance).pop()
                    );

                let min_theta_distances_and_corresponding_h_ave =  wrapped_data[phi_bin].into_iter().map(|point|
                                wrapped_data[phi_bin]
                                    .into_iter()
                                    .map(|second_point|
                                    if second_point.wrap == point.wrap{}
                                    else{
                                        ((point.theta-second_point.theta).abs(),(point.h+second_point.h)/2.0,(point.index,second_point.index))
                                    }
                                    ).sort_by(|(theta_distance,h,indices)| theta_distance).pop()
                                );

                for (index, (smallest_r_distance,h_r,indices)) in min_r_distances_and_corresponding_h_ave.enumerate() {
                    let (smallest_theta_distance,h_theta,indices) = min_theta_distances_and_corresponding_h_ave[index];

                    let h_max = f64::max(h_r,h_theta);
                    let lower_distance_bound = lower_distance_bound(delta_max,sigma_min,smallest_r_distance,smallest_theta_distance);

                    if h_max>lower_distance_bound{ //if the stream width is wider than this rough lower bound on the distance, there is a chance of collision

                    }


                    let distance = distance(coords_bin[point],coords_bin[second_point],50);


                }






            }



                }
            }
        }
        intersections

    }


}



#[derive(Serialize,Clone)]
pub struct StarChunk {
    pub(crate) geodesic_graph: GeodesicGraph,
    pub(crate) fraction_of_star: f64,
    pub(crate) binding_energy: f64,
}
#[derive(Serialize,Clone)]
pub struct Star {
    star_chunks: Vec<StarChunk>,
    stellar_params: StellarParams,
    r_initial:f64,
    theta_initial:f64,
    weighted_phi_values:Vec<((usize, [f64; 2]), f64)>
}
impl Star{
    pub fn new(stellar_params: StellarParams, r_initial: f64, theta_initial: f64,stellar_radius:f64,num_steps:usize,step_size:f64)-> Self{
        let star_chunks = initialize_star_chunks(stellar_params,r_initial, stellar_radius, 10,theta_initial,num_steps,step_size);
        let weighted_phi_values = Vec::new();
        Self{star_chunks,stellar_params,r_initial,theta_initial,weighted_phi_values}
    }
    /*
    pub fn weighted_phi_values(&mut self, t:f64)->Vec<((usize, [f64; 2]), f64)>{
        let mut phi_values = Vec::new();
        for star_chunk in self.star_chunks.clone(){
            phi_values.push((star_chunk.geodesic_graph.return_phi_at_t(t),star_chunk.fraction_of_star));
        }
        self.weighted_phi_values = phi_values.clone();
        phi_values
    }
    pub fn total_angular_momentum(&mut self, t:f64) ->(f64,f64){
        let phi_values = Star::weighted_phi_values(self, t);
        let mut total_angular_momentum = (0.0,0.0);
        for phi_value in phi_values{
            let (angle,weight) = phi_value;
            let angle = angle.1[1];
            total_angular_momentum = (total_angular_momentum.0 + (angle.cos())*weight,total_angular_momentum.1 + (angle.sin())*weight)
        }
        total_angular_momentum
    }

     */
    pub fn serialize(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let graph_json = serde_json::to_string(&self)?;
        let mut file = File::create(file_path)?;
        file.write_all(graph_json.as_bytes())?;
        Ok(())
    }

}











