use crate::A;
use crate::functions::{delta, sigma};
use crate::structs::StellarParams;

pub fn w_0(r:f64, theta:f64) -> [f64;4]{
    let s = (delta(r)/sigma(r,theta)).sqrt();
    [s,0.0,0.0,-A*s*theta.sin().powi(2)]
} //con
pub fn w_1(r:f64,theta:f64)-> [f64;4]{
    let s = (sigma(r,theta)/delta(r)).sqrt();
    [0.0,s,0.0,0.0]
} //con
pub fn w_2(r:f64, theta:f64)-> [f64;4] {
    let s = sigma(r,theta).sqrt();
    [0.0,0.0,s,0.0]
}  //con
pub fn w_3(r:f64, theta:f64)-> [f64;4] {
    let s = theta.sin()/(sigma(r,theta).sqrt());
    [A*s,0.0,0.0,-(r*r+A*A)*s]

}  //con
pub fn lambda_2(r:f64, theta:f64,r_dot:f64,theta_dot:f64,stellar_params: StellarParams) ->[f64;4] {
    let e = stellar_params.e;
    let lz = stellar_params.lz;
    let c = stellar_params.c;
    let sigma = sigma(r,theta);
    let delta = delta(r);

    let lambda_2_0 = (sigma/(c*delta)).sqrt()*A*theta.cos()*r_dot;
    let lambda_2_1 = (1.0/(c*sigma*delta)).sqrt()*A*theta.cos()*(e*(r*r+A*A)-A*lz);
    let lambda_2_2 = -(1.0/(c*sigma)).sqrt()*r*(A*e*theta.sin()-lz/theta.sin());
    let lambda_2_3 = (sigma/c).sqrt()*r*theta_dot;

    println!("lambda_2_0: {lambda_2_0}
    \n lambda_2_1: {lambda_2_1}
    \n lambda_2_2: {lambda_2_2}
    \n lambda_2_3: {lambda_2_3}");

    let symmetric_tetrad = [w_0(r,theta),w_1(r,theta),w_2(r,theta),w_3(r,theta)];
    let tetrad_coefficients = [lambda_2_0,lambda_2_1,lambda_2_2,lambda_2_3];

    let mut lambda_2 = [0.0,0.0,0.0,0.0];
    for boyer_lindquist_coord in 0..4{
        for basis_vector in 0..4{
            lambda_2[boyer_lindquist_coord] =
                lambda_2[boyer_lindquist_coord] +
                    symmetric_tetrad[basis_vector][boyer_lindquist_coord]*
                        tetrad_coefficients[basis_vector]
        }
    }
    lambda_2
} //con

