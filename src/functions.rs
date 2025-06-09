use crate::constants::{A,M};

pub fn theta_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{ // coefficients for 0th, cos()^2 and cos()^4
    let a0 = c;
    let a2 = -(c+A.powi(2)*(1.0-e.powi(2))+lz.powi(2));
    let a4 = A.powi(2)*(1.0-e.powi(2));
    [a0,0.0,a2,0.0,a4]
}
pub fn radial_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{
    let a4 = e.powi(2)-1.0;
    let a3 = 2.0*M;
    let a2 = e.powi(2)*2.0*A.powi(2)-2.0*A*lz*e-(A*e-lz).powi(2)-c-A.powi(2);
    let a1 = 2.0*M*(A*e-lz).powi(2)+2.0*M*c;
    let a0 = e.powi(2)*A.powi(4)-2.0*A*lz*e*A.powi(2)+A.powi(2)*lz.powi(2)-A.powi(2)*(A*e-lz).powi(2)-A.powi(2)*c;
   // println!("{:?}",[a0,a1,a2,a3,a4]);
    [a0,a1,a2,a3,a4]
}

pub fn coefficients_to_poly(x:f64, coefficients:[f64;5])-> f64{
    let mut sum = 0.0;
    for i in 0..coefficients.len() {
        sum += coefficients[i]*x.powi(i as i32)
    }
    sum
}

pub fn r_derivative_magnitude(x:f64, coefficients:[f64;5]) ->f64{ //dr/dlambda
    coefficients_to_poly(x,coefficients).abs().sqrt()
}
pub fn theta_derivative(x:f64, coefficients:[f64;5]) ->f64{ //dtheta/dlambda
    if A == 0.0{
        0.0
    }else{
        let c = x.cos();
        coefficients_to_poly(c,coefficients).abs().sqrt()/x.sin()
    }

}


fn delta(r:f64) -> f64 {
    r.powi(2) - 2.0 * M * r + A.powi(2)

}
fn p(r:f64, lz:f64,e:f64) -> f64 {
    e*(r.powi(2) + A.powi(2))-A*lz
}
fn phi_r(r:f64, lz:f64,e:f64) -> f64{
    let d = delta(r);
    let p = p(r,lz,e);
    if d ==0.0{
        panic!("divided by 0")
    }
    A*p/d
}
fn phi_theta(theta:f64,lz:f64) -> f64 {
   // println!("phi theta is {}, {},{}",lz/(1.0-theta.cos().powi(2)),lz,theta.cos().powi(2));
    lz/(1.0-theta.cos().powi(2))
}
pub fn phi_total(theta:f64,r:f64, lz:f64,e:f64) -> f64 {
   // println!("phi total is {}",phi_r(r,lz,e) + phi_theta(theta,lz) -A*e);
    phi_r(r,lz,e) + phi_theta(theta,lz) -A*e

}


