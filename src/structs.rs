use serde::{Serialize};
use std::fs::File;
use std::io::Write;

use crate::functions::{find_radial_parameters, find_theta_parameters,mino_to_bl_time};
use crate::numeric_integrators::{integrate_theta, integrate_r, integrate_phi, integrate_H, integrate_t};
use crate::tests::test_tetrads;

#[derive(Serialize)]
pub struct Graph {
    pub(crate) radial: Vec<[f64;2]>,
    pub(crate) theta: Vec<[f64;2]>,
    pub(crate) phi: Vec<[f64;2]>,
    pub(crate) t: Vec<[f64;2]>,
    pub(crate) self_intersections: Vec<(usize,usize)>,
    pub(crate) h: Vec<[f64;2]>
}


#[derive(Serialize)]
pub struct StarChunk {
    pub(crate) binding_energy: f64,
    pub(crate) z_angular_momentum: f64,
    pub(crate) carter_constant: f64,
    pub(crate) phi: Vec<[f64;2]>,

}
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

#[derive(Debug, Clone, Copy)]
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
impl StarChunk {
    pub fn new(lz:f64,e:f64, c:f64) -> Self {
        let radial_params = find_radial_parameters(lz, e, c);
        let theta_params = find_theta_parameters(lz, e, c); //9979
        let (radial_graph,armpositions) = integrate_r( 9979.0, 2.37896861,4.315260870513532, radial_params);

        println!("z minus is {}, and initial theta is then {:?}", theta_params.z_minus, theta_params.z_minus.sqrt().acos());
        let theta_graph = integrate_theta(theta_params.z_minus.sqrt().acos()+0.01, theta_params); //theta_params.z_minus.sqrt().acos()+0.01
        println!("theta initial is {:?}",theta_graph[0]); //theta_params.z_minus.sqrt().acos()+0.01
        let t_graph  = integrate_t(radial_graph.clone(),theta_graph.clone(),lz,e);
        let phi_graph = integrate_phi(radial_graph.clone(),theta_graph.clone(),lz,e);
        let phi_by_t = mino_to_bl_time(t_graph,phi_graph);

        Self{
            binding_energy:e,
            z_angular_momentum:lz,
            carter_constant:c,
            phi:phi_by_t
        }

    }

}
impl Graph {
    pub fn new(stellar_params: StellarParams,r_initial:f64,theta_initial:f64) -> Self {

        println!("Initializing a new Graph with: \n\
         (e,lz,c) = {:?} \n\
         (r_0,theta_0) = ({},{})",
                 stellar_params,r_initial,theta_initial);

        let radial_graph = integrate_r( r_initial, stellar_params);
        let theta_graph = integrate_theta(theta_initial, stellar_params); //theta_params.z_minus.sqrt().acos()+0.01
        let phi_graph = integrate_phi(radial_graph.clone(), theta_graph.clone(), stellar_params);
        let t_graph  = integrate_t(radial_graph.clone(),theta_graph.clone(),stellar_params);
        let mut intersection_points = Vec::new();
        let mut stream_width = Vec::new();

        Self{
            radial: radial_graph,
            theta: theta_graph,
            phi: phi_graph,
            t: t_graph,
            self_intersections:intersection_points,
            stream
        }
    }
    pub fn serialize_trajectory(&self, file_path: String) {
        let graph_json = serde_json::to_string(&self)?;
        let mut file = File::create(file_path)?;
        file.write_all(graph_json.as_bytes())?;
    }



    pub fn find_possible_self_intersection_points(mut self) ->Graph{

        for primary_interval_start_index in 0..self.radial.len()-1 {
            for secondary_interval_start_index in 0..primary_interval_start_index {
                let (A_phi,B_phi,A_r,B_r) = if self.radial[primary_interval_start_index][1] < self.radial[primary_interval_start_index+1][1]{
                   // println!("going up");
                    (self.phi[primary_interval_start_index][1],
                    self.phi[primary_interval_start_index+1][1],
                    self.radial[primary_interval_start_index][1]*self.theta[primary_interval_start_index][1].sin(),
                    self.radial[primary_interval_start_index+1][1]*self.theta[primary_interval_start_index+1][1].sin())


                }else{
                  //  println!("goingdown");
                    (self.phi[primary_interval_start_index+1][1],
                     self.phi[primary_interval_start_index][1],
                     self.radial[primary_interval_start_index+1][1]*self.theta[primary_interval_start_index+1][1].sin(),
                     self.radial[primary_interval_start_index][1]*self.theta[primary_interval_start_index][1].sin())
                };

                let (C_phi,D_phi,C_r,D_r) = if self.radial[secondary_interval_start_index][1] < self.radial[secondary_interval_start_index+1][1]{
                    (self.phi[secondary_interval_start_index][1],
                    self.phi[secondary_interval_start_index+1][1],
                    self.radial[secondary_interval_start_index][1]*self.theta[secondary_interval_start_index][1].sin(),
                    self.radial[secondary_interval_start_index+1][1]*self.theta[secondary_interval_start_index+1][1].sin())
                }else{
                    (self.phi[secondary_interval_start_index+1][1],
                     self.phi[secondary_interval_start_index][1],
                     self.radial[secondary_interval_start_index+1][1]*self.theta[secondary_interval_start_index+1][1].sin(),
                     self.radial[secondary_interval_start_index][1]*self.theta[secondary_interval_start_index][1].sin())
                };

                let A_phi = A_phi%(2.0* std::f64::consts::PI);
                let B_phi = B_phi%(2.0* std::f64::consts::PI);
                let C_phi = C_phi%(2.0* std::f64::consts::PI);
                let D_phi = D_phi%(2.0* std::f64::consts::PI);

                let delta_phi_1 =A_phi-C_phi;
                let delta_phi_2 = B_phi-D_phi;

                let delta_r_1 = A_r-D_r;
                let delta_r_2 = B_r-C_r;
               // println!("{}",delta_phi_1 * delta_phi_2);

                if (delta_phi_1 * delta_phi_2 < 0.0)  && (delta_r_2*delta_r_1<0.0) && (delta_r_2*delta_r_1>-1.0) && (delta_phi_1 * delta_phi_2 > -1.0) { // && (delta_r_2*delta_r_1>-1.0) && (delta_phi_1 * delta_phi_2 > -1.0)
                    println!("{}   {}",delta_phi_1 , delta_phi_2);
                    self.self_intersections.push((primary_interval_start_index,secondary_interval_start_index));



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

            println!("the intersection points are at {:?}", self.self_intersections);
            self

         }
}



