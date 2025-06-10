use egui_plot::PlotPoints;
use crate::constants::*;
use crate::functions::*;


pub fn root_hunt_and_peck<F: Fn(f64)->f64>(y_start: f64, graph: F) -> (f64, f64) { //find a root above and below a starting value given a graph with the desired roots



    let mut lowerbound = y_start;
    let mut upperbound = y_start+0.1;

    let lower_root = {
        if graph(y_start).abs() < 0.001{
            println!("taking lower root");
            y_start
        } else {
            loop {
                // println!("searching between {} and {}",lowerbound,y_start + 0.01);
                match roots::find_root_brent(lowerbound, y_start + 0.01, &graph, &mut 0.001) {
                    Ok(root) => { break root }
                    Err(_message) => {} //println!("{} not found at {}",message,lowerbound)
                };
                lowerbound += -0.01;
                if lowerbound < 0.0 {
                    panic!("No lower root found")
                }
            }
        }
    };


    let upper_root = loop {
        if E > 0.9{
            upperbound += 1.0;
        }
        else {
            upperbound += 0.01;
        }

        match roots::find_root_brent(y_start+0.01,upperbound,& graph,&mut 0.001){
            Ok(root) => { break root}
            Err(_message) => {}
        };
        if upperbound > 2000.0{
            panic!("No upper root found below search criteria");

        }
    };
    // println!("lower {} upper {}",lower_root,upper_root);
    (lower_root,upper_root)

}
pub fn integrate_psi(r_initial:f64, ) {



}





