use egui_plot::{PlotPoints,Line,Plot};
use crate::constants::*;
use crate::functions::*;
use crate::numeric_solvers::*;
use egui;
use plotters::prelude::*;


#[derive(Default)]
pub(crate) struct Graph {
    r_initial: f64,
    r_min: f64,
    r_max: f64,

    theta_initial: f64,
    theta_min: f64,
    theta_max: f64,

    data: (PlotPoints,PlotPoints,PlotPoints,Vec<(f64,f64,f64)>),

}
impl Graph { // radial coefficients, theta coefficients
    pub(crate) fn new( vals:((f64, f64, f64), (f64, f64, f64))) -> Self {

        let radial = integrate(vals.0.0,vals.0.1,vals.0.2, Coordinates::Radial, LZ, E, C);
        let angular = integrate(vals.1.0,vals.1.1,vals.1.2, Coordinates::Theta,LZ,E,C);
        println!("Radial is {:?}",radial.points());
        println!("theta is {:?}",angular.points());
        let data =find_phi(radial,angular,LZ,E);

        println!("{:?}",data.3);
        println!("r initial {}, r min {}, r max {}",vals.0.0,vals.0.1,vals.0.2);
        println!("theta initial {}, theta min {}, theta max {}",vals.1.0,vals.1.1,vals.1.2);

        Self {
            r_initial:vals.0.0,
            r_min:vals.0.1,
            r_max:vals.0.2,
            theta_initial:vals.1.0,
            theta_min:vals.1.1,
            theta_max:vals.1.2,
            data,
        }

    }
}