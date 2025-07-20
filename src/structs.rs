use serde::{Serialize};
use std::fs::File;
use std::io::Write;
use std::error::Error;
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

    pub fn find_self_intersections(mut self,num_phi_bins:usize){
        let mut intersections = Vec::with_capacity(self.num_steps);
        let binsize = 2.0*std::f64::consts::PI/num_phi_bins as f64;
        let mut indices = Vec::new();
        let mut wraps = Vec::new();
        for bin_number in 0..num_phi_bins{
            let mut points_bin = Vec::new();
            let mut coords_bin = Vec::new();
            for phi_index in 0..self.phi_graph.len(){
                let remainder = self.phi_graph[phi_index] % (2.0*std::f64::consts::PI);
                //println!("NDER {:?}   {:?}",self.phi_graph[phi_index],remainder);

                let phi_value = remainder;
                if (phi_value >= bin_number as f64*binsize) && (phi_value < (bin_number as f64+1.0)*binsize){
                    let wrap = ((self.phi_graph[phi_index] -remainder)/(2.0*std::f64::consts::PI)).round() as i32;
                    wraps.push(wrap);
                    points_bin.push((phi_index,wrap));
                    coords_bin.push((self.radial_graph[phi_index],self.theta_graph[phi_index]));
                }
            }
            if coords_bin.len() > 1{
                intersections.push(Self::test_for_intersection(coords_bin,points_bin.clone(),self.stream_height.clone()));
            }

            indices.push(points_bin)
        }
        print!(" PHI IS {:?}, wraphs are {:?}",self.phi_graph,wraps);
        println!("indices are at {:?}",indices);
        println!("intersections are at {:?}",intersections);


    }

    fn test_for_intersection(coords_bin:Vec<(f64,f64)>, points_bin:Vec<(usize,i32)>, stream_data:Vec<f64>) -> Vec<((usize,usize), bool, f64,f64,bool)> {
        let mut intersections = Vec::new();
        //(point one tested,point two tested), was there a collision?, how far away, stream width, was a more precise integration done,
        let (mut sigma_min, mut delta_max)  = if points_bin.len() == 0{
            panic!("Tried to find stream intersections between an empty set of points.")
        }else{
            (sigma(coords_bin[0].0,coords_bin[0].1),delta(coords_bin[0].0))
        };

        for point in 1..points_bin.len(){
            let sigma = sigma(coords_bin[point].0,coords_bin[point].1);
            let delta = delta(coords_bin[point].0);

            if sigma <sigma_min{
                sigma_min = sigma
            }
            if delta>delta_max{
                delta_max = delta
            }

        }
        println!("max out pf the points {:?} are delta is {delta_max} and sigma is {sigma_min}", coords_bin);
        for point in 0..points_bin.len(){
            for second_point in 0..point{
                if points_bin[point].1 != points_bin[second_point].1{
                    println!("{:?} and {:?}",points_bin[point].1,points_bin[second_point].1);
                    let average_stream_width = (stream_data[point] + stream_data[second_point])/2.0;
                    let delta_r = (coords_bin[point].0-coords_bin[second_point].0).abs();
                    let delta_theta = (coords_bin[point].1-coords_bin[second_point].1).abs();
                    let lower_distance_bound = lower_distance_bound(delta_max,sigma_min,delta_r,delta_theta);

                    if lower_distance_bound > average_stream_width{
                        intersections.push(((point,second_point),false,lower_distance_bound,average_stream_width,false));
                    }else{
                        let distance = distance(coords_bin[point],coords_bin[second_point],50);
                        if distance > average_stream_width{
                            intersections.push(((point,second_point),false,lower_distance_bound,average_stream_width,true));
                        } else{
                            intersections.push(((point,second_point),true,lower_distance_bound,average_stream_width,true));
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









