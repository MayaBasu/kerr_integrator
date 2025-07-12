use crate::derivatives::{r_derivative_propertime, theta_derivative};
use crate::structs::{Graph, Stream};
use crate::tetrads::lambda_2;

pub fn find_shortest_distance(graph:Graph,stream:Stream){
    let radial_graph = graph.radial.clone();
    let theta_graph = graph.theta.clone();
    let phi_graph = graph.phi.clone();

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
        let r_dot = r_derivative_propertime(r, theta,true);
        let theta_dot = theta_derivative(r,theta,true);

        let lambda_2 = lambda_2(r,theta,r_dot,theta_dot);

        //take the dot product:

        let mut dot_product = 0.0;

        for i in 0..4{
            dot_product = dot_product + lambda_2[i]*coordinate_difference[i]
        }
        println!("the dot product is {}",dot_product);
        println!("the stream width is {} and {}",stream.h[index_1][1] ,stream.h[index_2][1])


    }


}

