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
pub fn integrate(y_start:f64, y_min:f64,y_max:f64, coordinate:Coordinates,lz:f64,e:f64,c:f64)->PlotPoints{ //function should

    let coefficients:  [f64;5] = match &coordinate  {
        Coordinates::Radial => radial_coefficients(lz, e, c),
        Coordinates::Theta => theta_coefficients(lz, e, c),
    };

    let derivative= match &coordinate  {  // pick out derivative function
        Coordinates::Radial => r_derivative,
        Coordinates::Theta =>  theta_derivative,
    };

    let ( length_tolerance,time_tolerance)  = match &coordinate{
        Coordinates::Radial => {(1.0,1.0)}
        Coordinates::Theta => {(0.001,0.1)}
    };

    //   let (y_min,y_max)= match &coordinate {
    //        Coordinates::Radial => {( 1.6666666666666665,15.000000000000004)}
    //        Coordinates::Theta => {(1.0471975511965979,2.0943951023841776)}
    //   };



    let mut going_up:bool = true;
    let mut last_switch_location:f64 = 0.0;

    // println!("alsjdflkejlfjasekfjsd {}",y_start);
    let mut y = y_start;

    // println!("{}   {}",y_min,y_max);


    (0..NUM_ITERATIONS).map(|i| {
        if i % 100 == 0 {
            //  println!("{}, {}",i,y);
        }

        let x = i as f64 * DT;

        let increment = {if going_up {derivative(y,coefficients)} else {-derivative(y,coefficients)}}*DT;
        y += increment;

        if ((y-y_min).abs() < length_tolerance || (y-y_max).abs() < length_tolerance) && (x - last_switch_location).abs() > time_tolerance{

            going_up = !going_up;
            last_switch_location = x;
            //    println!("Switched at {}, going up is {}",last_switch_location,going_up)
        }

        [x, y]
    }).collect()

}

pub fn find_phi(r_plot_points: PlotPoints, t_plot_points: PlotPoints,   lz:f64, e:f64) -> (PlotPoints,PlotPoints,PlotPoints,Vec<(f64,f64,f64)>) { //->PlotPoints
    let r_points = r_plot_points.points().to_vec();
    let t_points = t_plot_points.points().to_vec();
    let mut running_phi = 0.0;
    let coords: Vec<_> = (0..r_points.len()).map(|i|{
        let r = r_points[i].y;
        let l = r_points[i].x;
        let theta = t_points[i].y;
        let phi = phi_total(theta,r,lz,e);
        (l, r,theta,phi)
    }).map(|(l,r,theta,phi_der)| {
        running_phi += phi_der*DT;
        (l,r,theta,running_phi, r*theta.sin()*running_phi.cos(),r*theta.sin()*running_phi.sin(),r*theta.cos())
    }
    ).collect();

    let flat_r = (0..r_points.len()).map(|i|{
        [coords[i].0,coords[i].1]
    }).collect();
    let flat_t = (0..r_points.len()).map(|i|{
        [coords[i].0,coords[i].2]
    }).collect();
    let flat_p = (0..r_points.len()).map(|i|{
        [coords[i].0,coords[i].3 % 2.0*std::f64::consts::PI]
    }).collect();

    let cartesian =
        (0..r_points.len()).map(|i|{
            (coords[i].4,coords[i].5,coords[i].6)
        }).collect();

    (flat_r,flat_t,flat_p,cartesian)
}


