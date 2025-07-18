#![allow(non_snake_case)]
use crate::{M,A};
use crate::structs::{GeodesicGraph, StellarParams};
use crate::functions::{find_radial_parameters,find_theta_parameters};
use crate::derivatives::{psi_derivative, chi_derivative, phi_derivative, H_acceleration, t_derivative};

const STEP_SIZE: f64 = 0.0001;
pub const NUM_STEPS: usize = 40000;

pub fn integrate_geodesic(r_initial:f64, theta_initial:f64, stellar_params: StellarParams) ->(Vec<[f64;2]>,Vec<[f64;2]>,Vec<[f64;2]>,Vec<[f64;2]>){
    let radial_trajectory = integrate_r(r_initial,stellar_params);
    let theta_trajectory = integrate_theta(theta_initial,stellar_params);
    let t_trajectory = integrate_t(radial_trajectory.clone(), theta_trajectory.clone(), stellar_params);
    let phi_trajectory = integrate_phi(radial_trajectory.clone(),theta_trajectory.clone(), stellar_params);

    (radial_trajectory,theta_trajectory,t_trajectory,phi_trajectory)
}
pub fn integrate_H(trajectory_graph:&GeodesicGraph, h_der_initial:f64, h_initial:f64, stellar_params: StellarParams) ->Vec<[f64;2]>{
    let mut h = h_initial;
    let mut h_velocity = h_der_initial;
    let mut phi_graph= Vec::with_capacity(NUM_STEPS);
    let mut steps_since_flip = 0;
    let mut flipped = false;
    for step in 0..155{
        println!("{:?}",trajectory_graph.radial_graph[step][1]);
        phi_graph.push([trajectory_graph.radial_graph[step][1],h]);

    }

    for step in 155..NUM_STEPS{
        steps_since_flip = steps_since_flip + 1;
        let d_tau = STEP_SIZE*(trajectory_graph.radial_graph[step][1].powi(2)+A*A*trajectory_graph.theta_graph[step][1].cos().powi(2));
        let velocity_increment = H_acceleration(trajectory_graph.radial_graph[step][1],
                                                trajectory_graph.theta_graph[step][1],
                                                h,
                                                stellar_params
        )*d_tau;
        h_velocity = h_velocity + velocity_increment;

        (flipped, steps_since_flip,h_velocity)  = if h < 0.0 && steps_since_flip>10 {
            if flipped{
                (false, 0,-h_velocity)}
            else{
                (true,0,-h_velocity)
            }
        }else{
            (flipped, steps_since_flip,h_velocity)
        };

        h = h + h_velocity*d_tau;



        phi_graph.push([trajectory_graph.radial_graph[step][1],h]);
    }
    phi_graph
}


fn integrate_r(r_initial:f64, stellar_params: StellarParams) -> Vec<[f64;2]>{
    //wrapper on the integrate psi function which allows the inputs and outputs of it to be in terms of r, e
    // ven if we integrate wrt psi
    let params = find_radial_parameters(stellar_params);
    let mut psi = r_to_psi(r_initial,params.e, params.p);
    let mut r_graph= Vec::with_capacity(NUM_STEPS);
  //  let mut r_arm_positions = Vec::new();
    //let mut step_since_last_arm_recording = 1000;
    for i in 0..NUM_STEPS{
       // step_since_last_arm_recording = step_since_last_arm_recording+1;
        let x = (i as f64)*STEP_SIZE;
        let increment = psi_derivative(psi,params, stellar_params.e)*STEP_SIZE;
        psi = psi + increment;
        let r = psi_to_r(psi,params.e,params.p);
        r_graph.push([x,r]);
        /*
        if ((r - r_min < 0.000001)||(r_max - r < 0.0001))  && (step_since_last_arm_recording > 1000){
            println!("{step_since_last_arm_recording}   {:?}   {:?} ",r - r_min, r_max-r);
            step_since_last_arm_recording = 0;
            r_arm_positions.push(i);
        }

         */


    }
  //  println!("graph {:?}", r_graph);
    r_graph

} //con
fn integrate_theta(theta_initial:f64, stellar_params: StellarParams) -> Vec<[f64;2]>{

    //wrapper on the integrate chi function which allows the inputs and outputs of it to be in terms of chi, e
    // ven if we integrate wrt chi
    let params = find_theta_parameters(stellar_params);
    let mut chi = theta_to_chi(theta_initial,params.z_minus);
    let mut theta_graph= Vec::with_capacity(NUM_STEPS);

    if A == 0.0{
        for i in 0..NUM_STEPS{
            theta_graph.push([(i as f64)*STEP_SIZE,theta_initial]);
        }
        return theta_graph
    }

    for i in 0..NUM_STEPS{
        let x = (i as f64)*STEP_SIZE;
        let increment = chi_derivative(chi,params)*STEP_SIZE;
        chi = chi + increment;
     //   println!("{:?}",[x,chi_to_theta(chi,params.z_minus)]);
        theta_graph.push([x,chi_to_theta(chi,params.z_minus)]);

    }
    //println!("THE CHANGED BCK GRAPH for chi IS {:?}", graph);
    theta_graph
} //con
fn integrate_phi(radial_graph: Vec<[f64;2]>, theta_graph:Vec<[f64;2]>, stellar_params: StellarParams) -> Vec<[f64;2]>{

    for i in 0..NUM_STEPS{ //check that the x axis are the same for these graphs
        assert_eq!(radial_graph[i][0],theta_graph[i][0]);
    }
    assert_eq!(STEP_SIZE,radial_graph[2][0]-radial_graph[1][0] ); // check that the stepsize is consistant between the graphs

    let mut phi = 0.0;
    let mut phi_graph= Vec::with_capacity(NUM_STEPS);
    for i in 0..NUM_STEPS{
        let increment = phi_derivative(
            radial_graph[i][1],
            theta_graph[i][1], stellar_params.lz,stellar_params.e)*STEP_SIZE;
        phi = phi + increment;
        phi_graph.push([radial_graph[i][0], phi]);

    }

    phi_graph

} //con
fn integrate_t(radial_graph: Vec<[f64;2]>, theta_graph:Vec<[f64;2]>, stellar_params: StellarParams) -> Vec<[f64;2]>{

    for i in 0..NUM_STEPS{ //check that the x axis are the same for these graphs
        assert_eq!(radial_graph[i][0],theta_graph[i][0]);
    }
    assert_eq!(STEP_SIZE,radial_graph[2][0]-radial_graph[1][0] ); // check that the stepsize is consistant between the graphs

    let mut t = 0.0;
    let mut t_graph= Vec::with_capacity(NUM_STEPS);
    for i in 0..NUM_STEPS{
        let increment = t_derivative(
            radial_graph[i][1],
            theta_graph[i][1], stellar_params)*STEP_SIZE;
        t = t + increment;
        t_graph.push([radial_graph[i][0], t]);

    }

    t_graph

} //con

fn psi_to_r(psi:f64, e:f64, p:f64) -> f64{
    p*M/(1.0 + e*psi.cos())
} //con
fn r_to_psi(r:f64,e:f64,p:f64) -> f64{
    ((p*M/r-1.0)/e).acos()

} //con
fn chi_to_theta(chi:f64, zminus:f64) -> f64{
    (zminus.sqrt()*chi.cos()).acos()
} //con
fn theta_to_chi(theta:f64, zminus:f64) -> f64{
   // let z = (theta.cos()).powi(2);
    (theta.cos()/(zminus).sqrt()).acos()
} //con

