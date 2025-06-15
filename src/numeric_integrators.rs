
use crate::constants::{A, E, M};
use crate::derivatives::{psi_derivative, chi_derivative, phi_derivative};
use crate::structs::{RadialParams,ThetaParams};
use crate::functions::delta;
const STEP_SIZE: f64 = 0.0001;
pub const NUM_STEPS: usize = 500000;


pub fn integrate_r(r_initial:f64, params:RadialParams) -> Box<[[f64;3]; NUM_STEPS]>{
    //wrapper on the integrate psi function which allows the inputs and outputs of it to be in terms of r, e
    // ven if we integrate wrt psi
    let psi = r_to_psi(r_initial,params.e, params.p);
    println!("psi is {}", psi);
    let mut graph = integrate_psi(psi, params);
    for i in 0..NUM_STEPS{
        graph[i][1] = psi_to_r(graph[i][1],params.e,params.p);
        graph[i][2] = psi_to_r(graph[i][2],params.e,params.p);
    }
    println!("THE CHANGEkjuihD BCK GRAPH IS {:?}", graph);
    graph
}
pub fn integrate_theta(theta_initial:f64, params:ThetaParams) -> Box<[[f64;2]; NUM_STEPS]>{
    //wrapper on the integrate chi function which allows the inputs and outputs of it to be in terms of chi, e
    // ven if we integrate wrt chi
    let chi = theta_to_chi(theta_initial,params.zminus);
 //   println!("chi is {}", chi);
    if A == 0.0{
        let mut graph =  Box::new([[0.0;2]; NUM_STEPS]);
        for i in 0..NUM_STEPS{
            graph[i][1] = theta_initial;
            graph[i][0] = (i as f64)*STEP_SIZE;
        }

        return graph
    }
    let mut graph = integrate_chi(chi, params);
    for i in 0..NUM_STEPS{
        graph[i][1] = chi_to_theta(graph[i][1],params.zminus);
    }
   // println!("THE CHANGED BCK GRAPH for chi IS {:?}", graph);
    graph
}
pub fn integrate_phi(radial_graph: Box<[[f64;3]; NUM_STEPS]>, theta_graph:&[[f64;2]; NUM_STEPS], lz:f64,e:f64) -> Box<[[f64;2]; NUM_STEPS]>{

    for i in 0..NUM_STEPS{ //check that the x axis are the same for these graphs
        assert_eq!(radial_graph[i][0],theta_graph[i][0]);
    }
    assert_eq!(STEP_SIZE,radial_graph[2][0]-radial_graph[1][0] ); // check that the stepsize is consistant between the graphs

    let mut phi = 0.0;
    let mut phi_graph= Box::new([[0.0;2]; NUM_STEPS]);
    for i in 0..NUM_STEPS{
        let increment = phi_derivative(
            radial_graph[i][1],
            theta_graph[i][1], lz,e)*STEP_SIZE;
        phi = phi + increment;
        phi_graph[i][1] = phi;

        phi_graph[i][0] = radial_graph[i][0];
    }

    phi_graph

}
fn integrate_psi(psi_initial:f64, params:RadialParams) -> Box<[[f64;3];NUM_STEPS]> {

    let mut psi = psi_initial;
    let mut graph= Box::new([[0.0;3]; NUM_STEPS]);
    let mut x = 0.0;

    for step in 0..NUM_STEPS{
        x = (step as f64)*STEP_SIZE;
        let increment = psi_derivative(psi,params, E)*STEP_SIZE;
        psi = psi + increment;
        graph[step][2] = psi_derivative((step as f64)*2.0*(std::f64::consts::PI)/(NUM_STEPS as f64), params, E);
        graph[step][1] = psi;
        graph[step][0] = x;

    }
    println!("{:?}",graph);
    graph
}
fn integrate_chi(chi_initial:f64, params:ThetaParams) ->Box<[[f64;2];NUM_STEPS]> {

    let mut chi = chi_initial;
    let mut graph= Box::new([[0.0;2]; NUM_STEPS]);
    let mut x = 0.0;

    for step in 0..NUM_STEPS{
        x = (step as f64)*STEP_SIZE;
        let increment = chi_derivative(chi,params)*STEP_SIZE;
        chi = chi + increment;
        graph[step][1] = chi;
        graph[step][0] = x;
    }
    graph
}

fn psi_to_r(psi:f64,e:f64,p:f64) -> f64{
    p*M/(1.0 + e*psi.cos())
}

fn r_to_psi(r:f64,e:f64,p:f64) -> f64{
    ((p*M/r-1.0)/e).acos()

}

fn chi_to_theta(chi:f64, zminus:f64) -> f64{
    let z = zminus*(chi.cos()).powi(2);
    (z.sqrt()).acos()
}
fn theta_to_chi(theta:f64, zminus:f64) -> f64{
    let z = (theta.cos()).powi(2);
    ((z/zminus).sqrt()).acos()
}

