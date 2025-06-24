use serde::{Serialize};

use crate::functions::{find_radial_parameters, find_theta_parameters};
use crate::numeric_integrators::{integrate_theta, integrate_r, integrate_phi};

#[derive(Serialize)]
pub struct Graph {
    pub(crate) radial: Vec<[f64;3]>,
    pub(crate) theta: Vec<[f64;2]>,
    pub(crate) phi: Vec<[f64;2]>
}

#[derive(Serialize)]
pub struct Stream {
    pub(crate) h: Vec<[f64;2]>
}

impl Graph {
    pub fn new(lz:f64,e:f64, c:f64) -> Self {
        let radial_params = find_radial_parameters(lz, e, c);
        let theta_params = find_theta_parameters(lz, e, c); //9979

        let radial_graph = integrate_r(9979.0/2.0, radial_params);
        println!("z minus is {}, and initial theta is then {:?}", theta_params.z_minus, theta_params.z_minus.sqrt().acos());
        let theta_graph = integrate_theta(theta_params.z_minus.sqrt().acos()+0.01, theta_params);
        let phi_graph = integrate_phi(radial_graph.clone(), &theta_graph, lz,e);

       // println!("phi daata {:?}", phi_graph);
        Self{
            radial: radial_graph.to_vec(),
            theta: theta_graph.to_vec(),
            phi: phi_graph.to_vec()
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