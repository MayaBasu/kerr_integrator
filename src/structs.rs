use std::ops::Deref;
use serde::{Serialize, Serializer};
use crate::constants::{ COS_I};
use crate::derivatives::phi_derivative;
use crate::functions::{find_radial_parameters, find_theta_parameters};
use crate::NUM_STEPS;
use crate::numeric_integrators::{integrate_theta, integrate_r, integrate_phi};

#[derive(Serialize)]

pub struct Graph {
    pub(crate) radial: Vec<[f64;3]>,
    pub(crate) theta: Vec<[f64;2]>,
    pub(crate) phi: Vec<[f64;2]>

}

impl Graph {
    pub fn new(lz:f64,e:f64, c:f64) -> Self {
        let radialparams = find_radial_parameters(lz, e, c);
        let thetaparams = find_theta_parameters(lz, e, c);

        let radial_graph = integrate_r( 1.6441, radialparams);
        let theta_graph = integrate_theta(COS_I.acos(),thetaparams);
        let phi_graph = integrate_phi(radial_graph.clone(), &theta_graph, lz,e);

        println!("phi daata {:?}", phi_graph);
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
    pub zplus:f64,
    pub zminus:f64,
}