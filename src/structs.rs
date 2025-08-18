use crate::functions::radial_roots;
use serde::{Serialize};
use std::fs::File;
use std::io::Write;
use std::error::Error;
use std::f32::consts::PI;
use crate::functions::{delta, distance, lower_distance_bound, sigma};
use crate::NUM_PHI_BINS;
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


#[derive(Clone, Copy, Debug,Serialize)]
pub struct WrappedDataPoint {
    pub index:usize,
    pub wrap:usize,
    pub phi:f64,
    pub r:f64,
    pub theta:f64,
    pub h:f64,


}

#[derive(Serialize,Clone)]
pub struct GeodesicGraph {
    pub stellar_params: StellarParams,

    pub num_steps: usize,
    pub step_size: f64,

    pub t_graph: Vec<f64>,
    pub radial_graph: Vec<f64>,
    pub theta_graph: Vec<f64>,
    pub phi_graph: Vec<f64>,

    pub stream_height: Vec<f64>,


    pub wrapped_data: Vec<Vec<WrappedDataPoint>>,

    pub intersections: Vec<(usize,usize,bool)>,
    pub bl_t_of_first_intersection: (f64,f64),
}

impl GeodesicGraph {
    pub fn new(stellar_params: StellarParams, r_initial: f64, theta_initial: f64, num_steps:usize,step_size:f64) -> Self {
        /*
        println!("Initializing a new Graph with: \n\
         (e,lz,c) = {:?} \n\
         (r_0,theta_0) = ({},{}) \n\
         It has radial roots at {:?}\
         \n Calculating trajectory for {num_steps} in Mino time spaced at {step_size}",
                 stellar_params, r_initial, theta_initial, radial_roots(stellar_params));

         */



        let radial_graph = integrate_r(r_initial, stellar_params,num_steps,step_size).0;
        let theta_graph = integrate_theta(theta_initial,stellar_params,num_steps,step_size).0;
        let t_graph = integrate_t(radial_graph.clone(), theta_graph.clone(), stellar_params,num_steps,step_size);
        let phi_graph = integrate_phi(radial_graph.clone(),theta_graph.clone(), stellar_params,num_steps,step_size);


        Self {
            stellar_params,

            num_steps,
            step_size,

            radial_graph,
            theta_graph,
            phi_graph,
            t_graph,

            wrapped_data:[(); NUM_PHI_BINS].map(|_| Vec::new()).to_vec(),

            stream_height: Vec::new(),
            intersections: Vec::new(),

            bl_t_of_first_intersection: (0.0,0.0),
        }
    }

    pub fn trim(mut self,last_index:usize)->Self{
        self.radial_graph.truncate(last_index);
        self.theta_graph.truncate(last_index);
        self.phi_graph.truncate(last_index);
        self.t_graph.truncate(last_index);
        self.stream_height.truncate(last_index);
        self


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

    pub fn return_phi_at_t(&self, global_time: f64) -> (usize,f64){
      //  println!("the last time is {:?}",self.t_graph.last().unwrap());
        //println!("The time calculated up to was {:?}", self.t_graph.last().unwrap());
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

    pub fn lag(&self, tidal_radius:f64,outer_radius:f64,theta_initial:f64)->(usize,f64,f64){
        let graph = GeodesicGraph::new(self.stellar_params, tidal_radius, theta_initial,1000,0.01);

        let mut past_difference = (graph.radial_graph[0] - outer_radius).abs();

        for i in 1..graph.radial_graph.len() {
            let current_difference = (graph.radial_graph[i] - outer_radius).abs();

            if past_difference < current_difference {
                let best_index = i - 1;
                return (best_index,graph.radial_graph[best_index],graph.t_graph[best_index])
            } else {
                past_difference = current_difference;
            }
        }
        panic!("Didn't find a closest match for the radial")


    }



    pub fn wrap_data(&mut self)->[Vec<WrappedDataPoint>; NUM_PHI_BINS]{
        let mut wrapped_data =  [(); NUM_PHI_BINS].map(|_| Vec::new());
        let bin_size = (2.0*std::f64::consts::PI)/NUM_PHI_BINS as f64;
        let mut last_bin = 0;
        for (index, &phi) in self.phi_graph.iter().enumerate(){
            let remainder = phi % (2.0*std::f64::consts::PI);
            let wrap = ((phi-remainder)/(2.0*std::f64::consts::PI)).round() as usize;
            let bin_number = (remainder/ bin_size ).floor() as usize;
            if ((bin_number as f64 -last_bin as f64)).abs()>1.0 && (last_bin !=NUM_PHI_BINS-1){
                println!("Skipped a bin!!!, before was {}, now is {}",last_bin,bin_number)
            }
            last_bin = bin_number;
       //     println!("pushing with bin unumber {:?} since remainder was {:?}",bin_number,remainder);
            wrapped_data[bin_number].push(WrappedDataPoint{
                index,
                wrap,
                phi,
                r:self.radial_graph[index],
                theta:self.theta_graph[index],
                h:self.stream_height[index]})
        }
        self.wrapped_data = wrapped_data.to_vec();

        wrapped_data
    }

    pub(crate) fn find_earliest_intersection(&mut self) -> Vec<(usize, usize)> {
        let mut lowest_upper_index = self.num_steps;
        let mut wrapped_data = self.wrap_data();
        let mut intersections = Vec::new();

      //  println!("{:?}", wrapped_data.len());
        for phi_bin in 0..wrapped_data.len(){
         //   println!("{:?}",wrapped_data[phi_bin]);

            if wrapped_data[phi_bin].len() < 2 {
              //  println!("continuing on {:?}",phi_bin);
                continue
            } //check that there are at least 2 points at this phi value
            else {
                wrapped_data[phi_bin].sort_by_key(|x| x.wrap);
                if wrapped_data[phi_bin][wrapped_data[phi_bin].len()-1].wrap == wrapped_data[phi_bin][0].wrap {
                   // println!("skipped because the wraps are {:?}",wrapped_data[phi_bin]);
                    continue }else { //check that there are points from at least two different wraps
                    //find the smallest sigma and largest delta
                    let mut sigmas = wrapped_data[phi_bin].clone().into_iter().map(|point| sigma(point.r, point.theta)).collect::<Vec<f64>>();
                    sigmas.sort_by(|a,b|a.partial_cmp(b).unwrap());
                    let sigma_min = sigmas[0];


                    let mut deltas = wrapped_data[phi_bin].clone().into_iter().map(|point| delta(point.r)).collect::<Vec<f64>>();
                    deltas.sort_by(|a,b|a.partial_cmp(b).unwrap());
                    let delta_max = deltas[deltas.len()-1];
                    if phi_bin ==40{
                       // println!("sigmas are {:?} and sigma min is {sigma_min}",sigmas);
                       // println!("deltas are {:?} and delta max is {delta_max}",deltas);
                        //println!("the wrapped data is {:?}",wrapped_data[phi_bin]);
                        for wrapped_data_point in &wrapped_data[phi_bin]{
                            //println!("The phi value is at {:?}",wrapped_data_point.phi % (2.0*std::f64::consts::PI))
                        }

                    }

                    //now we want to find the minimum distance delta r and delta theta between these points
                    for point in wrapped_data[phi_bin].clone().into_iter() {
                        if point.index > lowest_upper_index{
                            continue
                        }
                        for second_point in wrapped_data[phi_bin].clone().into_iter() {
                            if second_point.index > lowest_upper_index{
                                continue
                            }
                            if second_point.wrap == point.wrap {
                                continue
                            } else {
                                let delta_r = (point.r - second_point.r).abs();
                                let delta_theta = (point.theta - second_point.theta).abs();
                                let lower_distance_bound = lower_distance_bound(delta_max, sigma_min, delta_r, delta_theta);
                                let h_ave = (point.h + second_point.h) / 2.0;
                                if lower_distance_bound > h_ave {
                                    intersections.push((false, false, (point.index, second_point.index)));
                                } else {
                                    let first_coords = (point.r, point.theta);
                                    let second_coords = (second_point.r, second_point.theta);
                                    let distance = distance(first_coords, second_coords, 50);
                                    if distance > h_ave {
                                        intersections.push((false, true, (point.index, second_point.index)));
                                    } else {
                                        intersections.push((true, true, (point.index, second_point.index)));
                                        lowest_upper_index = point.index.max(second_point.index);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
       // println!("{:?}",intersections);


        let mut positive_intersections = Vec::new();

        for intersection in intersections{
            if intersection.0 == true{

                positive_intersections.push(intersection.2);

            }
        }
        positive_intersections
    }

    pub fn find_stop_time(&mut self)-> f64{
        let intersections = self.clone().find_earliest_intersection();
        let mut times = Vec::new();

        for intersection in intersections{
            let time = self.t_graph[intersection.0].max(self.t_graph[intersection.1]);
            //println!("The phi value is at {:?} and {:?}",self.phi_graph[intersection.0],self.phi_graph[intersection.1]);
            let phi = self.phi_graph[intersection.0];
            let remainder = phi % (2.0*std::f64::consts::PI);
            let bin_size = 2.0*std::f64::consts::PI/NUM_PHI_BINS as f64;
            let bin_number = (remainder/ bin_size ).floor() as usize;
            //println!("The bin is {:?} at time {:?}",bin_number,time);
            times.push(time)
        }

        if times.len()<1{
            println!("No intersection found, so returning a 0");
            return 0.0
        }
        times.sort_by(|a,b|a.partial_cmp(b).unwrap());

        let time_min = times[0];

       // println!("The times of intersection are {:?} and I am taking {:?}",times,time_min);
        time_min

    }

    pub fn angular_momentum(self, bl_time:f64)-> (f64,f64,f64){
        let C = self.stellar_params.c;
        let LZ = self.stellar_params.lz;
        let L = (C+LZ*LZ).sqrt();
        let cos_i = LZ/L;
        let lxy = cos_i.acos().sin()*L;

       // println!("{} {} {} {} {}",C,LZ,L,cos_i,lxy);

        let phi = self.return_phi_at_t(bl_time).1;

        (phi.cos()*lxy,phi.sin()*lxy,LZ)

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
    pub(crate) star_chunks: Vec<StarChunk>,
    stellar_params: StellarParams,
    r_initial:f64,
    theta_initial:f64,
    angular_momentum_over_time:Vec<(f64,f64,f64,f64)>,

}

impl StarChunk{
    pub fn trim(mut self,trim_index:usize)->Self{
        self.geodesic_graph = self.geodesic_graph.trim(trim_index);
        self
    }
}
impl Star{
    pub fn new(stellar_params: StellarParams, r_initial: f64, theta_initial: f64,stellar_radius:f64,num_steps:usize,step_size:f64,num_slices:usize)-> Self{
        let star_chunks = initialize_star_chunks(stellar_params,r_initial, stellar_radius, num_slices,theta_initial,num_steps,step_size);
        let angular_momentum_over_time =Vec::new();
        Self{star_chunks,stellar_params,r_initial,theta_initial,angular_momentum_over_time}
    }


    pub fn new_trimmed(stellar_params: StellarParams, r_initial: f64, theta_initial: f64,stellar_radius:f64,num_steps:usize,step_size:f64,num_slices:usize,t_trim:f64)-> Self{
        let mut star_chunks = initialize_star_chunks(stellar_params,r_initial, stellar_radius, num_slices,theta_initial,num_steps,step_size);
        let mut trimmed_star_chunks = Vec::new();
        let angular_momentum_over_time =Vec::new();

        for (index,star_chunk) in star_chunks.into_iter().enumerate(){


            let trim_index = star_chunk.geodesic_graph.return_phi_at_t(t_trim).0 + 2;
            if star_chunk.geodesic_graph.t_graph[trim_index] < t_trim{
                println!("Trimming at index {:?} which is at time {:?}",trim_index,star_chunk.geodesic_graph.t_graph[trim_index]);
            }

            trimmed_star_chunks.push(star_chunk.trim(trim_index))
        }

        Self{star_chunks:trimmed_star_chunks,stellar_params,r_initial,theta_initial,angular_momentum_over_time}
    }

    pub fn serialize(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let graph_json = serde_json::to_string(&self)?;
        let mut file = File::create(file_path)?;
        file.write_all(graph_json.as_bytes())?;
        Ok(())
    }
    pub fn calculate_total_angular_momentum(&self, time:f64) ->(f64,f64,f64) {
        let mut total_angular_momentum = (0.0,0.0,0.0);
        for star_chunk in self.star_chunks.clone(){
            let fraction_debris = star_chunk.fraction_of_star;
          //  println!("The fraction is {:?}",fraction_debris);
            let chunk_angular_momentum = star_chunk.geodesic_graph.angular_momentum(time);
            //println!("the chunk angular is {:?}",chunk_angular_momentum);
          //  println!("the angular momentum is {:?}",chunk_angular_momentum);

            total_angular_momentum.0 = total_angular_momentum.0 + chunk_angular_momentum.0*fraction_debris;
            total_angular_momentum.1 = total_angular_momentum.1 + chunk_angular_momentum.1*fraction_debris;
            total_angular_momentum.2 = total_angular_momentum.2 + chunk_angular_momentum.2*fraction_debris;
        }
        total_angular_momentum
    }

    pub fn find_angular_momentum_over_time(&mut self, stop_time:f64, divisions:usize) -> Vec<(f64,f64,f64,f64)> {
        let mut angular_momentums = Vec::new();
        for division in 1..divisions+1{

            let time = (stop_time/divisions as f64)*division as f64;

            let total_angular = self.calculate_total_angular_momentum(time);
            println!("Total Angular Momentum at  {:?} is {:?}, division {:?} of {:?}",time,total_angular,division,divisions);
            angular_momentums.push((time,total_angular.0,total_angular.1,total_angular.2))
        }
        self.angular_momentum_over_time = angular_momentums.clone();
        angular_momentums
    }
}











