#![allow(non_snake_case)]
use crate::{M,A};
use roots::{find_roots_quartic, Roots, find_roots_quadratic};
use crate::structs::{RadialParams, StellarParams, ThetaParams};

pub(crate) fn find_radial_parameters(stellar_params: StellarParams) ->  RadialParams{
    //This function finds the parameters p,e,p3,p4 from the appendex of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    //using the roots of the radial polynomial - this is to shift from E,Lx,C which determine the radial polynomial, to this other parameterization.

    let [r4,r3,r2,r1] = quartic_root_finder(radial_coefficients(stellar_params));

    let p = 2.0*(r1*r2)/(M*(r1+r2));
    let e = (r1-r2)/(r1+r2);
    let p4 = r4*(1.0+e)/M;
    let p3 = r3*(1.0-e)/M;

    RadialParams{
        p,
        e,
        p3,
        p4,
    }


} //clean
pub(crate) fn find_theta_parameters(stellar_params: StellarParams) ->  ThetaParams{
    //This function finds the parameters p,e,p3,p4 from the appendex of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    //using the roots of the radial polynomial - this is to shift from E,Lx,C which determine the radial polynomial, to this other parameterization.
    if A==0.0{
        return ThetaParams{
            beta:0.0,
            z_plus:0.0,
            z_minus:0.0
        }
    }
    let [zminus,zplus] = quadratic_root_finder(theta_coefficients(stellar_params.clone()));
    println!("z_minus is {zminus} and z_plus is {zplus}");
    let beta = A*A*(1.0-stellar_params.e*stellar_params.e);

    ThetaParams{
        beta,
        z_plus: zplus,
        z_minus: zminus
    }
}  //con
pub fn theta_coefficients(stellar_params: StellarParams) ->[f64;3]{ // coefficients for 0th, cos()^2 and cos()^4
    let c = stellar_params.c;
    let e = stellar_params.e;
    let lz = stellar_params.lz;

    let a2 = 1.0;
    let a1 = -1.0*(
        c+lz*lz+A*A*(1.0-e*e)
    )/(A*A*(1.0-e*e));
    let a0= c/(A*A*(1.0-e*e));
    [a0,a1,a2]
}  //con
pub fn radial_coefficients(stellar_params: StellarParams) ->[f64;5]{
    let lz = stellar_params.lz;
    let e = stellar_params.e;
    let c = stellar_params.c;

    let a4 = e.powi(2)-1.0;
    let a3 = 2.0*M;
    let a2 = e.powi(2)*A.powi(2)-lz.powi(2)-c-A.powi(2);
    let a1 = 2.0*M*(lz.powi(2)+c+A.powi(2)*e.powi(2)-2.0*A*e*lz);
    let a0 = -A.powi(2)*c;
   // println!("{:?}",[a0,a1,a2,a3,a4]);

    [a0,a1,a2,a3,a4]
} //clean
fn quadratic_root_finder(coefficients:[f64;3]) -> [f64; 2] {
    //this function takes in coefficients of a quadratic and outputs a vector of roots
    //it is possible the roots overlap for a certain parameter set, but we want to panic and invesigate this case.

    match find_roots_quadratic(coefficients[2], coefficients[1], coefficients[0]) {
        Roots::Two(roots) =>{
            println!("succesfully found two roots: {:?}", roots);
            roots
        }
        _ => {
            panic!("Two distinct roots were not found!")
        }
    }
}
pub fn radial_roots(stellar_params: StellarParams) -> [f64; 4] {
    quartic_root_finder(radial_coefficients(stellar_params))
}
fn quartic_root_finder(coefficients:[f64;5]) -> [f64; 4] {
    //this function takes in coefficients of a quartic and outputs a vector of roots
    //it is possible the roots overlap for a certain parameter set, but we want to panic and invesigate this case.
    println!("THE COEFFICIENTS ARE {:?}", coefficients);
    match find_roots_quartic(coefficients[4], coefficients[3], coefficients[2], coefficients[1], coefficients[0]) {
        Roots::Four(roots) =>{
            println!("succesfully found four roots: {:?}", roots);
            roots
        }
        Roots::Three(roots) =>{
            println!("succesfully found three roots: {:?}", roots);
            [0.0,roots[0], roots[1],roots[2]]
        }
        Roots::Two(roots) =>{
            println!("succesfully found two roots: {:?}", roots);
            [0.0,0.0, roots[0], roots[1]]
        }

        _ => {
            panic!("Four distinct roots were not found!")
        }
    }
} //clean
pub fn delta(r:f64) -> f64{

    r*r-2.0*M*r+A*A
}
pub fn sigma(r:f64, theta:f64) -> f64{
    r.powi(2) + A.powi(2)*(theta.cos()).powi(2)

}
pub fn K(stellar_params: StellarParams)->f64{  //Batra paper
    let lz = stellar_params.lz;
    let c = stellar_params.c;
    let e = stellar_params.e;
    (lz-A*e)*(lz-A*e)+c

}
pub fn R(theta:f64,stellar_params: StellarParams)-> f64{

    K(stellar_params)-A.powi(2)*(theta.cos()).powi(2)
}
pub fn S(r:f64,stellar_params: StellarParams)-> f64{
    r.powi(2)+K(stellar_params)
}
pub fn P(r:f64,stellar_params: StellarParams)->f64{
    let e = stellar_params.e;
    let lz = stellar_params.lz;

    e*(r*r+A*A)-A*lz
}
pub fn mino_to_bl_time(t_graph:Vec<[f64;2]>, mut other_graph: Vec<[f64;2]>)->Vec<[f64;2]>{
    //function to convert between one graph, such as phi, or r, or theta, in mino time, to a graph in boyer lindquist time
    for point in 0..other_graph.len(){
        other_graph[point][0] = t_graph[point][1]
    }
    other_graph
}
pub fn lower_distance_bound(delta_max:f64,sigma_min:f64,delta_r:f64,delta_theta:f64)->f64{
    ((sigma_min/delta_max)*delta_r.powi(2)+ sigma_min*delta_theta.powi(2)).sqrt()
}

pub fn distance(first_coords:(f64,f64),second_coords:(f64,f64),num_divisions:i32) -> f64{
    let delta_r = first_coords.0-second_coords.0;
    let delta_theta = first_coords.1-second_coords.1;

    let mut distance = 0.0;

    for division in 0..num_divisions{
        let r = first_coords.0+division as f64*(delta_r/num_divisions as f64);
        let theta = first_coords.1+division as f64*(delta_theta/num_divisions as f64);
        distance = distance + lower_distance_bound(delta(r),sigma(r,theta),delta_r,delta_theta);
    }
    distance
}



/*
pub fn rest_mass_squared(theta_min:f64)->f64{
    (
        (C)/((theta_min.cos()).powi(2))
        -(LZ/theta_min.sin()).powi(2)
    )/(A*A) +E*E
}
 */