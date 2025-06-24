#![allow(non_snake_case)]

use crate::constants::{A, E, M};
use crate::derivatives::{psi_derivative, chi_derivative, phi_derivative, H_acceleration};
use crate::structs::{Graph, RadialParams, ThetaParams};

const STEP_SIZE: f64 = 0.00001;
pub const NUM_STEPS: usize = 400000;


pub fn integrate_r(r_initial:f64, params:RadialParams) -> Vec<[f64;2]>{
    //wrapper on the integrate psi function which allows the inputs and outputs of it to be in terms of r, e
    // ven if we integrate wrt psi
    let mut psi = r_to_psi(r_initial,params.e, params.p);
    let mut r_graph= Vec::with_capacity(NUM_STEPS);

    for i in 0..NUM_STEPS{
        let x = (i as f64)*STEP_SIZE;
        let increment = -psi_derivative(psi,params, E)*STEP_SIZE;
        psi = psi + increment;
        r_graph.push([x,psi_to_r(psi,params.e,params.p)]);

    }
  //  println!("graph {:?}", r_graph);
    r_graph

}
pub fn integrate_theta(theta_initial:f64, params:ThetaParams) -> Vec<[f64;2]>{
    //wrapper on the integrate chi function which allows the inputs and outputs of it to be in terms of chi, e
    // ven if we integrate wrt chi
    println!("theta initial is {theta_initial}");
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
        theta_graph.push([x,chi_to_theta(chi,params.z_minus)]);

    }
    //println!("THE CHANGED BCK GRAPH for chi IS {:?}", graph);
    theta_graph
}

pub fn integrate_phi(radial_graph: Vec<[f64;2]>, theta_graph:Vec<[f64;2]>, lz:f64,e:f64) -> Vec<[f64;2]>{

    for i in 0..NUM_STEPS{ //check that the x axis are the same for these graphs
        assert_eq!(radial_graph[i][0],theta_graph[i][0]);
    }
    assert_eq!(STEP_SIZE,radial_graph[2][0]-radial_graph[1][0] ); // check that the stepsize is consistant between the graphs

    let mut phi = 0.0;
    let mut phi_graph= Vec::with_capacity(NUM_STEPS);
    for i in 0..NUM_STEPS{
        let increment = phi_derivative(
            radial_graph[i][1],
            theta_graph[i][1], lz,e)*STEP_SIZE;
        phi = phi + increment;
        phi_graph.push([radial_graph[i][0], phi]);

    }

    phi_graph

}


fn psi_to_r(psi:f64,e:f64,p:f64) -> f64{
    p*M/(1.0 + e*psi.cos())
}

fn r_to_psi(r:f64,e:f64,p:f64) -> f64{
    ((p*M/r-1.0)/e).acos()

}

fn chi_to_theta(chi:f64, zminus:f64) -> f64{
    (zminus.sqrt()*chi.cos()).acos()
}
fn theta_to_chi(theta:f64, zminus:f64) -> f64{
    let z = (theta.cos()).powi(2);
    ((z/zminus).sqrt()).acos()
}



pub fn integrate_H(trajectory_graph:&Graph,h_der_initial:f64,h_initial:f64)->Vec<[f64;2]>{
    let mut h = h_initial;
    let mut h_velocity = h_der_initial;
    let mut phi_graph= Vec::with_capacity(NUM_STEPS);
    let mut steps_since_flip = 0;
    let mut flipped = false;

    for step in 0..NUM_STEPS{
        steps_since_flip = steps_since_flip + 1;
        let d_tau = STEP_SIZE*(trajectory_graph.radial[step][1].powi(2)+A*A*trajectory_graph.theta[step][1].cos().powi(2));
        let velocity_increment = H_acceleration(trajectory_graph.radial[step][1],
                                                trajectory_graph.theta[step][1],
                                                h
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
        phi_graph.push([trajectory_graph.radial[step][1],h]);
    }
    phi_graph
}

