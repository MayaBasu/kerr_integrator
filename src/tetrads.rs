use crate::constants::{A, C, E, LZ};
use crate::functions::{delta, sigma};

pub fn w_0(r:f64, theta:f64) -> [f64;4]{
    let s = (delta(r)/sigma(r,theta)).sqrt();
    [s,0.0,0.0,-A*s*theta.sin().sqrt()]
}
pub fn w_1(r:f64,theta:f64)-> [f64;4]{
    let s = (sigma(r,theta)/delta(r)).sqrt();
    [0.0,s,0.0,0.0]
}
pub fn w_2(r:f64, theta:f64)-> [f64;4] {
    let s = sigma(r,theta).sqrt();
    [0.0,0.0,s,0.0]
}
pub fn w_3(r:f64, theta:f64)-> [f64;4] {
    let s = theta.sin()/(sigma(r,theta).sqrt());
    [A*s,0.0,0.0,-(r*r+A*A)*s]

}

pub fn lambda_2(r:f64, theta:f64,r_dot:f64,theta_dot:f64) ->[f64;4] {
    let sigma = sigma(r,theta);
    let delta = delta(r);

    let lambda_2_0 = (sigma/(C*delta)).sqrt()*A*theta.cos()*r_dot;
    let lambda_2_1 = (1.0/(C*sigma*delta)).sqrt()*A*theta.cos()*(E*(r*r+A*A)-A*LZ);
    let lambda_2_2 = -(1.0/(C*sigma)).sqrt()*r*(A*E*theta.sin()-LZ/theta.sin());
    let lambda_2_3 = (sigma/C).sqrt()*r*theta_dot;

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
}

