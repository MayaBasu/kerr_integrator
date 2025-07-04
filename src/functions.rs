#![allow(non_snake_case)]
use roots::{find_roots_quartic, Roots, find_roots_quadratic};
use crate::constants::{A, M,K,C,E,LZ};
use crate::structs::{RadialParams, ThetaParams};

//Hellsing, Integra et al (2006)
pub(crate) fn find_radial_parameters(lz:f64, e:f64, c:f64) ->  RadialParams{
    //This function finds the parameters p,e,p3,p4 from the appendex of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    //using the roots of the radial polynomial - this is to shift from E,Lx,C which determine the radial polynomial, to this other parameterization.

    let [r4,r3,r2,r1] = quartic_root_finder(radial_coefficients(lz, e, c));


    let p = 2.0*(r1*r2)/(M*(r1+r2));
    let e = (r1-r2)/(r1+r2);
    let p4 = r4*(1.0+e)/M;
    let p3 = r3*(1.0-e)/M;
    println!("p is {} and e is {} for {} {} {}",p,e,lz,e,c);

    RadialParams{
        p,
        e,
        p3,
        p4,
    }


} //clean
pub(crate) fn find_theta_parameters(lz:f64, e:f64, c:f64) ->  ThetaParams{
    //This function finds the parameters p,e,p3,p4 from the appendex of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    //using the roots of the radial polynomial - this is to shift from E,Lx,C which determine the radial polynomial, to this other parameterization.
    if A==0.0{
        return ThetaParams{
            beta:0.0,
            z_plus:0.0,
            z_minus:0.0
        }
    }
    let [zminus,zplus] = quadratic_root_finder(theta_coefficients(lz, e, c));
    println!("zminus is {zminus} and zplus is {zplus}");
    let beta = A*A*(1.0-e*e);
    println!("theta parameter {}   {}   {}",zminus,zplus,beta);

    ThetaParams{
        beta,
        z_plus: zplus,
        z_minus: zminus
    }
}
pub fn theta_coefficients(lz:f64, e:f64, c:f64) ->[f64;3]{ // coefficients for 0th, cos()^2 and cos()^4
    let a2 = 1.0;
    let a1 = -1.0*(
        c+lz*lz+A*A*(1.0-e*e)
    )/(A*A*(1.0-e*e));
    let a0= c/(A*A*(1.0-e*e));
    [a0,a1,a2]
}
pub fn radial_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{
    let a4 = e.powi(2)-1.0;
    let a3 = 2.0*M;
    let a2 = e.powi(2)*A.powi(2)-lz.powi(2)-c-A.powi(2);
    let a1 = 2.0*M*(lz.powi(2)+c+A.powi(2)*e.powi(2)-2.0*A*e*lz);
    let a0 = -A.powi(2)*c;
   // println!("{:?}",[a0,a1,a2,a3,a4]);

    [a0,a1,a2,a3,a4]
} //clean
fn quadratic_root_finder(coeffs:[f64;3]) -> [f64; 2] {
    //this function takes in coefficients of a quadratic and outputs a vector of roots
    //it is possible the roots overlap for a certain parameter set, but we want to panic and invesigate this case.

    match find_roots_quadratic(coeffs[2], coeffs[1], coeffs[0]) {
        Roots::Two(roots) =>{
            println!("succesfully found two roots: {:?}", roots);
            roots
        }
        _ => {
            panic!("Two distinct roots were not found!")
        }
    }
}
fn quartic_root_finder(coeffs:[f64;5]) -> [f64; 4] {
    //this function takes in coefficients of a quartic and outputs a vector of roots
    //it is possible the roots overlap for a certain parameter set, but we want to panic and invesigate this case.
    println!("THE COEFFICIENTS ARE {:?}", coeffs);
    match find_roots_quartic(coeffs[4], coeffs[3], coeffs[2], coeffs[1], coeffs[0]) {
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
pub fn R(theta:f64)-> f64{

    K-A.powi(2)*(theta.cos()).powi(2)
}
pub fn S(r:f64)-> f64{
    r.powi(2)+K
}

pub fn P(r:f64)->f64{
    E*(r*r+A*A)-A*LZ
}

pub fn rest_mass_squared(theta_min:f64)->f64{
    (
        (C)/((theta_min.cos()).powi(2))
        -(LZ/theta_min.sin()).powi(2)
    )/(A*A) +E*E
}
pub fn w_0(r:f64,theta:f64)-> [f64;4]{
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

