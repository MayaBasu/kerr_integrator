
use crate::constants::*;
pub fn psi_derivative(psi:f64, p:f64, e:f64, p3:f64, p4:f64 ) -> f64 {
    // calculate d\psi/d \lambda using A16 from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.044015

    M*(
        (1.0-E*E)
        *((p-p3)-e*(p+p3*psi.cos()))
        *((p-p4)+e*(p-p4*psi.cos()))
    ).sqrt()/(1.0-E*E)

}

