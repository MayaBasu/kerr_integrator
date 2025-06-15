
use crate::constants::*;
use crate::functions::delta;
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

}



pub fn chi_derivative(chi:f64, params:ThetaParams) -> f64 {
    let beta = params.beta;
    let zplus = params.zplus;
    let zminus = params.zminus;
    //A3 from  https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015

    (
        beta*
            (zplus-
                zminus* (chi.cos()).powi(2))
    ).sqrt()


}


pub fn phi_derivative(r:f64,theta:f64, lz:f64, e:f64) -> f64{
    ((theta.sin()).recip()).powi(2)*lz
        +A*e*(
        (r.powi(2)+A*A)/(delta(r)) -1.0
    )
    -A*A*lz/(delta(r))
}


