#![allow(non_snake_case)]

use std::ops::Deref;
use std::sync::Arc;
use crate::constants::{A, E, M};
use crate::derivatives::{psi_derivative, chi_derivative, phi_derivative, H_acceleration};
use crate::structs::{Graph, RadialParams, ThetaParams};

const STEP_SIZE: f64 = 0.0001;
pub const NUM_STEPS: usize = 40000;


pub fn integrate_r(r_initial:f64, params:RadialParams) -> Arc<[[f64;2]; NUM_STEPS]>{
    //wrapper on the integrate psi function which allows the inputs and outputs of it to be in terms of r, e
    // ven if we integrate wrt psi
    let mut psi = r_to_psi(r_initial,params.e, params.p);
    let mut r_graph= [[0.0;2]; NUM_STEPS];

    for i in 0..NUM_STEPS{
        let x = (i as f64)*STEP_SIZE;
        let increment = -psi_derivative(psi,params, E)*STEP_SIZE;
        psi = psi + increment;

        r_graph[i][0] = x;
        r_graph[i][1] = psi_to_r(psi,params.e,params.p);

    }
    Arc::new(r_graph)
}








pub fn integrate_theta(theta_initial:f64, params:ThetaParams) -> Arc<[[f64;2]; NUM_STEPS]>{
    //wrapper on the integrate chi function which allows the inputs and outputs of it to be in terms of chi, e
    // ven if we integrate wrt chi
    println!("theta initial is {theta_initial}");
    let chi = theta_to_chi(theta_initial,params.z_minus);
 //   println!("chi is {}", chi);
    if A == 0.0{
        let  graph =  Arc::new([[0.0;2]; NUM_STEPS]);
        let mut theta_graph = [[0.0;2]; NUM_STEPS];
        for i in 0..NUM_STEPS{
            theta_graph[i][1] = theta_initial;
            theta_graph[i][0] = (i as f64)*STEP_SIZE;
        }

        return Arc::new(theta_graph)
    }
    let mut graph = integrate_chi(chi, params);
    let mut theta_graph = [[0.0;2]; NUM_STEPS];
    for i in 0..NUM_STEPS{
        theta_graph[i][1] = chi_to_theta(graph[i][1],params.z_minus);
    }
    //println!("THE CHANGED BCK GRAPH for chi IS {:?}", graph);
    Arc::new(theta_graph)
}
pub fn integrate_phi(radial_graph: Arc<[[f64;3]; NUM_STEPS]>, theta_graph:&[[f64;2]; NUM_STEPS], lz:f64,e:f64) -> Arc<[[f64;2]; NUM_STEPS]>{

    for i in 0..NUM_STEPS{ //check that the x axis are the same for these graphs
        assert_eq!(radial_graph[i][0],theta_graph[i][0]);
    }
    assert_eq!(STEP_SIZE,radial_graph[2][0]-radial_graph[1][0] ); // check that the stepsize is consistant between the graphs

    let mut phi = 0.0;
    let mut phi_graph= [[0.0;2]; NUM_STEPS];
    for i in 0..NUM_STEPS{
        let increment = phi_derivative(
            radial_graph[i][1],
            theta_graph[i][1], lz,e)*STEP_SIZE;
        phi = phi + increment;
        phi_graph[i][1] = phi;

        phi_graph[i][0] = radial_graph[i][0];
    }

    Arc::new(phi_graph)

}

fn integrate_chi(chi_initial:f64, params:ThetaParams) ->Arc<[[f64;2];NUM_STEPS]> {

    let mut chi = chi_initial;
   // println!("chi is {chi}");
    let mut graph= [[0.0;2]; NUM_STEPS];


    for step in 0..NUM_STEPS{
        let x = (step as f64)*STEP_SIZE;
        let increment = chi_derivative(chi,params)*STEP_SIZE;
        chi = chi + increment;
      //  println!("increment is {}",increment);
        graph[step][1] = chi;
        graph[step][0] = x;
    }
    Arc::new(graph)
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



pub fn integrate_H(trajectory_graph:&Graph,h_der_initial:f64,h_initial:f64)->Arc<[[f64;2]; NUM_STEPS]>{
    let mut h = h_initial;
    let mut h_velocity = h_der_initial;
 //   println!("h is {h}");
    let mut graph= [[0.0;2]; NUM_STEPS];
    let mut steps_since_flip = 0;
    let mut flipped = false;

    for step in 0..NUM_STEPS{
   //     println!("h is {h}");
        steps_since_flip = steps_since_flip + 1;
        let d_tau = STEP_SIZE*(trajectory_graph.radial[step][1].powi(2)+A*A*trajectory_graph.theta[step][1].cos().powi(2));
        //println!("radial {:?} and theta {:?}", trajectory_graph.radial[step][1], trajectory_graph.theta[step][1]);
        let velocity_increment = H_acceleration(trajectory_graph.radial[step][1],
                                                trajectory_graph.theta[step][1],
                                                h
        )*d_tau;
        h_velocity = h_velocity + velocity_increment;
     //   println!("{} and stepps since {} ",h < 0.002,steps_since_flip);
       (flipped, steps_since_flip,h_velocity)  = if h < 0.0 && steps_since_flip>10 {
          //  println!("flipped");
           if flipped{
               (false, 0,-h_velocity)}
           else{
               (true,0,-h_velocity)
           }

       }else{
           (flipped, steps_since_flip,h_velocity)
       };


        h = h + h_velocity*d_tau;
        //  println!("increment is {}",increment);
        graph[step][1] = h;
        graph[step][0] = trajectory_graph.radial[step][1];

    }
 //   println!("{:?}", graph);

    Arc::new(graph)
}

