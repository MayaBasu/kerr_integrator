use serde::{Serialize};
use std::fs::File;
use std::io::Write;
use std::error::Error;
use std::usize;
use itertools::Itertools;
use crate::derivatives::{r_derivative_propertime, theta_derivative};
use crate::functions::{mino_to_bl_time, radial_roots};
use crate::numeric_integrators::{ integrate_H, integrate_geodesic};
use crate::star_initialization::initialize_star_chunks;
use crate::tetrads::lambda_2;

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
impl StellarParams{
    pub fn new(lz:f64,e:f64, c:f64) -> Self{
        Self{e,lz,c}
    }
}

#[derive(Serialize,Clone)]
pub struct GeodesicGraph {
    pub(crate) stellar_params: StellarParams,

    pub(crate) radial_graph: Vec<[f64;2]>,
    pub(crate) theta_graph: Vec<[f64;2]>,
    pub(crate) phi_graph: Vec<[f64;2]>,
    pub(crate) t_graph: Vec<[f64;2]>,

    pub(crate) stream_height: Vec<[f64;2]>,
    pub(crate) possible_self_intersections: Vec<(usize,usize)>,
    pub(crate) distance_at_possible_self_intersections: Vec<(usize,usize,f64)>,
    pub(crate) intersections: Vec<(usize,usize,bool)>,

    pub(crate) bl_t_of_first_intersection: (f64,f64),
}

impl GeodesicGraph {
    pub fn new(stellar_params: StellarParams, r_initial: f64, theta_initial: f64) -> Self {
        println!("Initializing a new Graph with: \n\
         (e,lz,c) = {:?} \n\
         (r_0,theta_0) = ({},{}) \n\
         It has radial roots at {:?}",
                 stellar_params, r_initial, theta_initial, radial_roots(stellar_params));

        let (radial_graph, theta_graph, t_graph, phi_graph) = integrate_geodesic(r_initial, theta_initial, stellar_params);

        let mut possible_self_intersections = Vec::new();
        let mut stream_height = Vec::new();
        let mut distance_at_possible_self_intersections = Vec::new();
        let mut intersections = Vec::new();

        Self {
            radial_graph,
            theta_graph,
            phi_graph,
            t_graph,
            possible_self_intersections,
            stream_height,
            stellar_params,
            distance_at_possible_self_intersections,
            intersections,
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
        self.stream_height = integrate_H(&self, h_dot_initial, h_initial, self.stellar_params)
    }
    pub fn find_possible_self_intersection_points(mut self) -> GeodesicGraph {
        for primary_interval_start_index in 0..self.radial_graph.len() - 1 {
            for secondary_interval_start_index in 0..primary_interval_start_index {
                let (a_phi, b_phi, a_r, b_r) = if self.radial_graph[primary_interval_start_index][1] < self.radial_graph[primary_interval_start_index + 1][1] {
                    // println!("going up");
                    (self.phi_graph[primary_interval_start_index][1],
                     self.phi_graph[primary_interval_start_index + 1][1],
                     self.radial_graph[primary_interval_start_index][1] * self.theta_graph[primary_interval_start_index][1].sin(),
                     self.radial_graph[primary_interval_start_index + 1][1] * self.theta_graph[primary_interval_start_index + 1][1].sin())
                } else {
                    //  println!("goingdown");
                    (self.phi_graph[primary_interval_start_index + 1][1],
                     self.phi_graph[primary_interval_start_index][1],
                     self.radial_graph[primary_interval_start_index + 1][1] * self.theta_graph[primary_interval_start_index + 1][1].sin(),
                     self.radial_graph[primary_interval_start_index][1] * self.theta_graph[primary_interval_start_index][1].sin())
                };

                let (c_phi, d_phi, c_r, d_r) = if self.radial_graph[secondary_interval_start_index][1] < self.radial_graph[secondary_interval_start_index + 1][1] {
                    (self.phi_graph[secondary_interval_start_index][1],
                     self.phi_graph[secondary_interval_start_index + 1][1],
                     self.radial_graph[secondary_interval_start_index][1] * self.theta_graph[secondary_interval_start_index][1].sin(),
                     self.radial_graph[secondary_interval_start_index + 1][1] * self.theta_graph[secondary_interval_start_index + 1][1].sin())
                } else {
                    (self.phi_graph[secondary_interval_start_index + 1][1],
                     self.phi_graph[secondary_interval_start_index][1],
                     self.radial_graph[secondary_interval_start_index + 1][1] * self.theta_graph[secondary_interval_start_index + 1][1].sin(),
                     self.radial_graph[secondary_interval_start_index][1] * self.theta_graph[secondary_interval_start_index][1].sin())
                };

                let a_phi = a_phi % (2.0 * std::f64::consts::PI);
                let b_phi = b_phi % (2.0 * std::f64::consts::PI);
                let c_phi = c_phi % (2.0 * std::f64::consts::PI);
                let d_phi = d_phi % (2.0 * std::f64::consts::PI);

                let delta_phi_1 = a_phi - c_phi;
                let delta_phi_2 = b_phi - d_phi;

                let delta_r_1 = a_r - d_r;
                let delta_r_2 = b_r - c_r;
                // println!("{}",delta_phi_1 * delta_phi_2);

                if (delta_phi_1 * delta_phi_2 < 0.0) && (delta_r_2 * delta_r_1 < 0.0) && (delta_r_2 * delta_r_1 > -1.0) && (delta_phi_1 * delta_phi_2 > -1.0) { // && (delta_r_2*delta_r_1>-1.0) && (delta_phi_1 * delta_phi_2 > -1.0)
                    println!("{}   {}", delta_phi_1, delta_phi_2);
                    self.possible_self_intersections.push((primary_interval_start_index, secondary_interval_start_index));
                }
            }
        }


        //code to search through each arm position end point, and, for each index less than this, I want to compute the deltas of adjacent points

        /*
        for arm_end_position_index in 0..armpositions.len()-1{
            let arm_start = armpositions[arm_end_position_index];
            let arm_end = armpositions[arm_end_position_index+1];

            for primary_interval_start_index in arm_start..arm_end{
            for secondary_interval_start_index in 0..arm_start{
                let delta_phi_1 = phi_graph[secondary_interval_start_index][1]-phi_graph[primary_interval_start_index][1];
                let delta_phi_2 = phi_graph[secondary_interval_start_index+1][1]-phi_graph[primary_interval_start_index+1][1];
                if delta_phi_1*delta_phi_2 < 0.0{
                    intersection_points.push(primary_interval_start_index)
                }


        }
            }
        }

  */

        self
    }
    pub fn find_shortest_distance(mut self) {
        for intersection_point in self.possible_self_intersections {
            let index_1 = intersection_point.0;
            let index_2 = intersection_point.1;

            let coordinate_difference = [
                0.0,
                self.radial_graph[index_1][1] - self.radial_graph[index_2][1],
                self.theta_graph[index_1][1] - self.theta_graph[index_2][1],
                self.phi_graph[index_1][1] - self.phi_graph[index_2][1]];

            let r = self.radial_graph[index_1][1];
            let theta = self.theta_graph[index_1][1];
            let r_dot = r_derivative_propertime(r, theta, true, self.stellar_params);
            let theta_dot = theta_derivative(r, theta, true, self.stellar_params);

            let lambda_2 = lambda_2(r, theta, r_dot, theta_dot, self.stellar_params);

            //take the dot product:
            let mut dot_product = 0.0;

            for i in 0..4 {
                dot_product = dot_product + lambda_2[i] * coordinate_difference[i]
            }
            self.distance_at_possible_self_intersections.push((index_1,index_2,dot_product))
        }
    }

    pub fn find_intersections(mut self) {
        let mut intersections = Vec::new();

        for i in 0..self.distance_at_possible_self_intersections.len(){
            let (index_1,index_2,dot_product) = self.distance_at_possible_self_intersections[i];

            if (self.stream_height[index_1][1]+self.stream_height[index_2][1])/(2.0) > dot_product{
                intersections.push((index_1,index_2,true));
            }
            else{
                intersections.push((index_1,index_2,false));
            }

        }
        self.intersections = intersections;
    }

    pub fn return_phi_at_t(self, global_time: f64) -> (usize,[f64;2]){
        assert!(self.t_graph.last().unwrap()[1] > global_time);

        let mut past_difference = (self.t_graph[0][1] - global_time).abs();

        for i in 1..self.t_graph.len() {
            let current_difference = (self.t_graph[i][1] - global_time).abs();

            if past_difference < current_difference {
                let best_index = i - 1;
                return (best_index,self.phi_graph[best_index])
            } else {
                past_difference = current_difference;
            }
        }
        panic!("Didn't find a closest match")
    }
}
pub struct StarChunk {
    pub(crate) geodesic_graph: GeodesicGraph,
    pub(crate) fraction_of_star: f64,
    pub(crate) binding_energy: f64,
}

pub struct Star {
    star_chunks: Vec<StarChunk>,
    stellar_params: StellarParams,
    r_initial:f64,
    theta_initial:f64,
}

impl Star {
    pub fn new(self, stellar_params: StellarParams, r_initial: f64, theta_initial: f64,stellar_radius:f64)-> Self{
        let star_chunks = initialize_star_chunks(stellar_params,r_initial, stellar_radius, 10,theta_initial);
        Self{star_chunks,stellar_params,r_initial,theta_initial}
    }
}





