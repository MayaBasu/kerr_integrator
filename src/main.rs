
mod functions;

mod numeric_integrators;
use std::fs::File;
use std::io::Write;
mod constants;

use crate::functions::*;
use crate::constants::*;
use crate::numeric_integrators::*;

use std::error::Error;
use crate::derivatives::{r_derivative_propertime, theta_derivative};
use crate::tetrads::*;
mod derivatives;
mod structs;
mod tetrads;

use crate::structs::{Graph, Stream};

fn main() -> Result<(), Box<dyn Error>>{

    println!("{:?}",find_radial_parameters(LZ,E,C));
    println!("E: {}  \nL_z: {}  \nC: {} \nL: {} \ncos I: {}", E,LZ,C,(C+LZ*LZ).sqrt(), C.sqrt()/((C+LZ*LZ).sqrt()));


    let graph = Graph::new(LZ,E,C);
    let stream = Stream{
        h:integrate_H(&graph,0.0,0.90).to_vec() //2.35
    };

    let radial_graph = graph.radial.clone();
    let theta_graph = graph.theta.clone();

    let mut derivatives = Vec::new();

    let phi_graph = graph.phi.clone();

/*
    for intersection_point in graph.self_intersections.clone() {
        let index_1 = intersection_point.0;
        let index_2 = intersection_point.1;

        let coordinate_difference = [
            0.0,
            radial_graph[index_1][1] - radial_graph[index_2][1],
            theta_graph[index_1][1] - theta_graph[index_2][1],
            phi_graph[index_1][1] - phi_graph[index_2][1]];

        let r = radial_graph[index_1][1];
        let theta = theta_graph[index_1][1];
        let r_dot = r_derivative_propertime(r, theta);
        let theta_dot = theta_derivative(r,theta);

        let lambda_2 = lambda_2(r,theta,r_dot,theta_dot);

        //take the dot product:

        let mut dot_product = 0.0;

        for i in 0..4{
            dot_product = dot_product + lambda_2[i]*coordinate_difference[i]
        }
        println!("the dot product is {}",dot_product);
        println!("the stream width is {} and {}",stream.h[index_1][1] ,stream.h[index_2][1])


    }

 */


    let graph_json = serde_json::to_string(&graph)?;
    let mut file = File::create("ballistic_graph.json")?;
    file.write_all(graph_json.as_bytes())?;


    let graph_stream = serde_json::to_string(&stream)?;
    let mut file = File::create("stream_width.json")?;
    file.write_all(graph_stream.as_bytes())?;



    for point in 0..radial_graph.len(){
        if point % 10 == 1{
            derivatives.push(
                [theta_derivative(radial_graph[point][1], theta_graph[point][1],false),
                    ((theta_graph[point][1]-theta_graph[point-1][1])/(0.0001))/(radial_graph[point][1].powi(2)+A*A*theta_graph[point][1].cos().powi(2))]
            )
        }

    }
    println!("the derivatives and approximate derivatives are {:?}",derivatives);

    let testpointindex = 700;
    let r = radial_graph[testpointindex][1];
    let theta = theta_graph[testpointindex][1];
    println!("at point number {testpointindex}  the values are \n r: {:?} \n theta: {:?}  \n r dot: {:?} \n theta dot: {:?}",
             r,
             theta,
        r_derivative_propertime(r,theta,false),
             theta_derivative(r,theta,false)
    );
    println!("delta is {:?} and sigma is {:?}",delta(r),sigma(r, theta));
    println!("the tetrad is \n w0 {:?} ", w_0(r,theta));
    println!("the tetrad is \n w1 {:?} ", w_1(r,theta));
    println!("the tetrad is \n w2 {:?} ", w_2(r,theta));
    println!("the tetrad is \n w3 {:?} ", w_3(r,theta));
    let lambda_2 = lambda_2(r,theta,r_derivative_propertime(r,theta,false),theta_derivative(r,theta,false));
    println!("the lambda 2 in BL coordinates is {:?}",lambda_2);



    Ok(())
}



