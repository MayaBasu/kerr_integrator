#![allow(non_snake_case)]
use crate::{M,A};
use crate::structs::{GeodesicGraph, StellarParams};
use crate::functions::{find_radial_parameters,find_theta_parameters};
use crate::derivatives::{psi_derivative, chi_derivative, phi_derivative, H_acceleration, t_derivative};



pub fn integrate_H(trajectory_graph:&GeodesicGraph, h_der_initial:f64, h_initial:f64, stellar_params: StellarParams,num_steps:usize,step_size:f64) ->Vec<f64>{
    let mut h = h_initial;
    let mut h_velocity = h_der_initial;
    let mut phi_graph= Vec::with_capacity(num_steps);
    let mut steps_since_flip = 0;
    let mut flipped = false;
    for step in 0..15{
        println!("{:?}",trajectory_graph.radial_graph[step]);
        phi_graph.push(h);

    }

    for step in 15..num_steps{
        steps_since_flip = steps_since_flip + 1;
        let d_tau = step_size*(trajectory_graph.radial_graph[step].powi(2)+A*A*trajectory_graph.theta_graph[step].cos().powi(2));
        let velocity_increment = H_acceleration(trajectory_graph.radial_graph[step],
                                                trajectory_graph.theta_graph[step],
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



        phi_graph.push(h);
    }
    phi_graph
}
pub fn integrate_t(radial_graph: Vec<f64>, theta_graph:Vec<f64>, stellar_params: StellarParams,num_steps:usize,step_size:f64) -> Vec<f64>{
    let mut t = 0.0;
    let mut t_graph= Vec::with_capacity(num_steps);
    for i in 0..num_steps{
        let increment = t_derivative(
            radial_graph[i],
            theta_graph[i], stellar_params)*step_size;
        t = t + increment;
        t_graph.push(radial_graph[i]);
    }
    t_graph
} //con

pub(crate) fn integrate_r(r_initial:f64, stellar_params: StellarParams, num_steps:usize,step_size:f64) -> (Vec<f64>, f64, usize){
    let params = find_radial_parameters(stellar_params);
    let mut psi = r_to_psi(r_initial,params.e, params.p);
    let mut r_graph= Vec::with_capacity(num_steps);
    for i in 0..num_steps{
        psi = psi + psi_derivative(psi,params, stellar_params.e)*step_size;
        r_graph.push(psi_to_r(psi,params.e,params.p));
    }
    (r_graph,step_size,num_steps)

} //con
pub fn integrate_theta(theta_initial:f64, stellar_params: StellarParams, num_steps:usize,step_size:f64) -> (Vec<f64>,f64,usize){
    let params = find_theta_parameters(stellar_params);
    let mut chi = theta_to_chi(theta_initial,params.z_minus);
    let mut theta_graph= Vec::with_capacity(num_steps);
    if A == 0.0{
        for i in 0..num_steps{
            theta_graph.push(theta_initial);
        }
        return (theta_graph, step_size,num_steps)
    }

    for i in 0..num_steps{
        chi = chi + chi_derivative(chi,params)*step_size;
        theta_graph.push(chi_to_theta(chi,params.z_minus));

    }
    (theta_graph,step_size,num_steps)
} //con
pub fn integrate_phi(radial_graph: Vec<f64>, theta_graph:Vec<f64>, stellar_params: StellarParams, num_steps:usize,step_size:f64) -> Vec<f64>{
    let mut phi = 0.0;
    let mut phi_graph= Vec::with_capacity(num_steps);
    for i in 0..num_steps{
        let increment = phi_derivative(
            radial_graph[i],
            theta_graph[i], stellar_params.lz,stellar_params.e)*step_size;
        phi = phi + increment;
        phi_graph.push(phi);
    }
    phi_graph
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

