#![allow(non_snake_case)]
use crate::constants::{A,M};
use crate::functions::{delta, R, S, sigma, radial_coefficients, theta_coefficients};
use crate::structs::{RadialParams, ThetaParams};

pub fn psi_derivative(psi:f64, params:RadialParams ,ee:f64) -> f64 {
    // calculate d\psi/d \lambda using A16 from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    let p = params.p;
    let e = params.e;
    let p3 = params.p3;
    let p4 = params.p4;

    M*(
        (1.0-ee*ee)
        *((p-p3)-e*(p+p3*psi.cos()))
        *((p-p4)+e*(p-p4*psi.cos()))
    ).sqrt()/(1.0-e*e)

} //con
pub fn chi_derivative(chi:f64, params:ThetaParams) -> f64 {
    let beta = params.beta;
    let zplus = params.z_plus;
    let zminus = params.z_minus;
    //A3 from  https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015

    (
        beta*
            (zplus-
                zminus* (chi.cos()).powi(2))
    ).sqrt()

} //con
pub fn phi_derivative(r:f64,theta:f64, lz:f64, e:f64) -> f64{
    ((theta.sin()).recip()).powi(2)*lz
        +A*e*(
        (r.powi(2)+A*A)/(delta(r)) -1.0
    )
    -A*A*lz/(delta(r))
}  //con
pub fn r_derivative_propertime(r:f64, theta:f64,negative:bool) -> f64{
    let [a0,a1,a2,a3,a4] = radial_coefficients(LZ, E, C);
 //from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    let sign = if negative{-1.0} else {1.0};
    (
        a4*r.powi(4)+
            a3*r.powi(3)+
            a2*r.powi(2)+
            a1*r.powi(1) +
            a0
    ).sqrt()/(r*r+A*A*theta.cos().powi(2))*sign

} //con
pub fn theta_derivative(r:f64,theta:f64, negative:bool)->f64{ //con
    //from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    let [a0,a1,a2] = theta_coefficients(LZ,E,C);
    let z = theta.cos().powi(2);
    let sign = if negative{-1.0} else {1.0};
    ((a2*z*z + a1*z + a0)*(A*A*(1.0-E*E))/(1.0-z)).sqrt()/(r*r+A*A*theta.cos().powi(2))*sign

} //con
pub fn t_derivative(r:f64,theta:f64,lz:f64,e:f64)->f64{
    //derivative wrt mino time, https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015
    let first_term = e*((r.powi(2)+A.powi(2)).powi(2)/delta(r)
        -A.powi(2)*theta.sin().powi(2)
    );
    let second_term = A*lz*(1.0-
        (r.powi(2)+A.powi(2))/delta(r));
    first_term+second_term
} //con
pub fn H_acceleration(r:f64,theta:f64, H:f64)->f64{
    let sigma = sigma(r,theta);
    let R = R(theta);
    let S = S(r);
    //Equation 8 from ``General relativistic stream crossing in tidal disruption events"
    -H*(
        (r/(sigma.powi(3)))*
            (r.powi(2)-3.0*A*A*theta.cos().powi(2))*
            (1.0+
                3.0*(r.powi(2)*R.powi(2)-A.powi(2)*(theta.cos()).powi(2)*S.powi(2))/
                    (K*sigma.powi(2)))
    )

} // ?

