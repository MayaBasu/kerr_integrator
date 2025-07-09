use serde::{Serialize};

use crate::functions::{find_radial_parameters, find_theta_parameters};
use crate::numeric_integrators::{integrate_theta, integrate_r, integrate_phi, integrate_H};

#[derive(Serialize)]
pub struct Graph {
    pub(crate) radial: Vec<[f64;2]>,
    pub(crate) theta: Vec<[f64;2]>,
    pub(crate) phi: Vec<[f64;2]>,
    pub(crate) self_intersections: Vec<(usize,usize)>
}

#[derive(Serialize)]
pub struct Stream {
    pub(crate) h: Vec<[f64;2]>
}

impl Graph {
    pub fn new(lz:f64,e:f64, c:f64) -> Self {
        let radial_params = find_radial_parameters(lz, e, c);
        let theta_params = find_theta_parameters(lz, e, c); //9979

        let (radial_graph,armpositions) = integrate_r( 9979.0, 2.37896861,4.315260870513532, radial_params);

        println!("z minus is {}, and initial theta is then {:?}", theta_params.z_minus, theta_params.z_minus.sqrt().acos());
        let theta_graph = integrate_theta(1.57, theta_params); //theta_params.z_minus.sqrt().acos()+0.01
        println!("theta initial is {:?}",theta_graph[0]); //theta_params.z_minus.sqrt().acos()+0.01
        let phi_graph = integrate_phi(radial_graph.clone(), theta_graph.clone(), lz,e); //theta_params.z_minus.sqrt().acos()+0.01

        let mut intersection_points = Vec::new();
        /*
        for primary_interval_start_index in 0..radial_graph.len()-1 {
            for secondary_interval_start_index in 0..primary_interval_start_index {
                let (A_phi,B_phi,A_r,B_r) = if radial_graph[primary_interval_start_index][1] < radial_graph[primary_interval_start_index+1][1]{
                   // println!("going up");
                    (phi_graph[primary_interval_start_index][1],
                    phi_graph[primary_interval_start_index+1][1],
                    radial_graph[primary_interval_start_index][1]*theta_graph[primary_interval_start_index][1].sin(),
                    radial_graph[primary_interval_start_index+1][1]*theta_graph[primary_interval_start_index+1][1].sin())


                }else{
                  //  println!("goingdown");
                    (phi_graph[primary_interval_start_index+1][1],
                     phi_graph[primary_interval_start_index][1],
                     radial_graph[primary_interval_start_index+1][1]*theta_graph[primary_interval_start_index+1][1].sin(),
                     radial_graph[primary_interval_start_index][1]*theta_graph[primary_interval_start_index][1].sin())
                };

                let (C_phi,D_phi,C_r,D_r) = if radial_graph[secondary_interval_start_index][1] < radial_graph[secondary_interval_start_index+1][1]{
                    (phi_graph[secondary_interval_start_index][1],
                    phi_graph[secondary_interval_start_index+1][1],
                    radial_graph[secondary_interval_start_index][1]*theta_graph[secondary_interval_start_index][1].sin(),
                    radial_graph[secondary_interval_start_index+1][1]*theta_graph[secondary_interval_start_index+1][1].sin())
                }else{
                    (phi_graph[secondary_interval_start_index+1][1],
                     phi_graph[secondary_interval_start_index][1],
                     radial_graph[secondary_interval_start_index+1][1]*theta_graph[secondary_interval_start_index+1][1].sin(),
                     radial_graph[secondary_interval_start_index][1]*theta_graph[secondary_interval_start_index][1].sin())
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
                    intersection_points.push((primary_interval_start_index,secondary_interval_start_index));



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

            println!("the intersection points are at {:?}", intersection_points);

         */






       // println!("phi daata {:?}", phi_graph);
        Self{
            radial: radial_graph.to_vec(),
            theta: theta_graph.to_vec(),
            phi: phi_graph.to_vec(),
            self_intersections:intersection_points
        }
    }
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


