use roots::{find_roots_quartic, Roots};
use crate::constants::{A, M};

pub fn theta_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{ // coefficients for 0th, cos()^2 and cos()^4
    let a0 = c;
    let a2 = -(c+A.powi(2)*(1.0-e.powi(2))+lz.powi(2));
    let a4 = A.powi(2)*(1.0-e.powi(2));
    [a0,0.0,a2,0.0,a4]
}
//Hellsing et al (2006)
pub fn radial_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{
    let a4 = e.powi(2)-1.0;
    let a3 = 2.0*M;
    let a2 = e.powi(2)*2.0*A.powi(2)-2.0*A*lz*e-(A*e-lz).powi(2)-c-A.powi(2);
    let a1 = 2.0*M*(A*e-lz).powi(2)+2.0*M*c;
    let a0 = e.powi(2)*A.powi(4)-2.0*A*lz*e*A.powi(2)+A.powi(2)*lz.powi(2)-A.powi(2)*(A*e-lz).powi(2)-A.powi(2)*c;
   // println!("{:?}",[a0,a1,a2,a3,a4]);
    [a0,a1,a2,a3,a4]
}


pub fn root_finder(coeffs:[f64;5]) -> [f64; 4] {
    //this function takes in coefficients of a quartic and outputs a vector of roots
    //it is possible the roots overlap for a certain parameter set, but we want to panic and invesigate this case.

    match find_roots_quartic(coeffs[4], coeffs[3], coeffs[2], coeffs[1], coeffs[0]) {
        Roots::Four(roots) =>{
            println!("succesfully found four roots: {:?}", roots);
            roots
        }
        _ => {
            panic!("Four distinct roots were not found!")
        }
    }
}


pub(crate) fn parameter_converter(roots: [f64; 4] ) ->  (f64,f64,f64,f64) {

    let [r4,r3,r2,r1] = roots;
    let p = 2.0/(M*(r1.recip()+r2.recip()));
    let e = (p*M)/(r2)-1.0;
    let e2 = 1.0-(p*M)/r1;
    println!("First with r2 is is {} and second with r1 is {}", e,e2);
    let p4 = r4*(1.0+e)/M;
    let p3 = r3*(1.0-e)/M;
    (p,e,p3,p4)
}




pub fn psi_to_r(psi:f64,e:f64,p:f64) -> f64{
    p*M/(1.0 + e*psi.cos())
}

pub fn r_to_psi(r:f64,e:f64,p:f64) -> f64{
    ((p*M/r-1.0)/e).acos()
}



