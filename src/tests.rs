use crate::derivatives::{r_derivative_propertime, theta_derivative};
use crate::structs::Graph;
use crate::tetrads::{lambda_2, w_0,w_1,w_2,w_3};
use crate::functions::{delta,sigma};
use crate::constants::{A};

pub fn test_tetrads(graph:Graph, testpointindex:u32){
    let radial_graph = graph.radial.clone();
    let theta_graph = graph.theta.clone();

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


}

pub fn test_derivatives(graph:Graph){
    let radial_graph = graph.radial.clone();
    let theta_graph = graph.theta.clone();
    let mut derivatives = Vec::new();

    for point in 0..radial_graph.len(){
        if point % 10 == 1{
            derivatives.push(
                [theta_derivative(radial_graph[point][1], theta_graph[point][1],false),
                    ((theta_graph[point][1]-theta_graph[point-1][1])/(0.0001))/(radial_graph[point][1].powi(2)+A*A*theta_graph[point][1].cos().powi(2))]
            )
        }

    }
    println!("the derivatives and approximate derivatives are {:?}",derivatives)


}
