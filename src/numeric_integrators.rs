use crate::derivatives::psi_derivative;
use crate::functions::*;
const STEP_SIZE: f64 = 0.0001;
const NUM_STEPS: usize = 100;


pub fn integrate_r(r_initial:f64, e:f64, p:f64, p3:f64, p4:f64) -> [[f64;2]; NUM_STEPS]{
    let psi = r_to_psi(r_initial,e, p);
    println!("psi is {}", psi);
    let mut graph = integrate_psi(psi, e, p, p3, p4);
    for i in 0..NUM_STEPS{
        graph[i][1] = psi_to_r(graph[i][1],e,p);
    }
    println!("THE CHANGED BCK GRAPH IS {:?}", graph);
    graph


}



pub fn integrate_psi(psi_initial:f64, e:f64, p:f64, p3:f64, p4:f64) ->[[f64;2];NUM_STEPS] {

    let mut psi = psi_initial;

    println!("psi and psi initial {}, {}", psi, psi_initial);
    let mut x = 0.0;

    let mut graph= [[0.0;2]; NUM_STEPS];

    for step in 0..NUM_STEPS{
        x = (step as f64)*STEP_SIZE;
        println!("inputs are {:?}", (psi,e,p,p3,p4));
        let increment = psi_derivative(psi,e,p,p3,p4)*STEP_SIZE;
        println!("incremnet is {}", increment);
        psi = psi + increment;
        graph[step][1] = psi;
        graph[step][0] = x;
    }
    println!("alskejf {:?}", graph);
    graph
}





