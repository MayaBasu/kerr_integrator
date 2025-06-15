use roots::{find_roots_quartic, Roots, find_roots_quadratic};
use crate::constants::{A, M};
use crate::structs::{RadialParams, ThetaParams};


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
            zplus:0.0,
            zminus:0.0
        }
    }
    let [zminus,zplus] = quadratic_root_finder(theta_coefficients(lz, e, c));

    let beta = A*A*(1.0-e*e);
    println!("theta parameter {}   {}   {}",zminus,zplus,beta);

    ThetaParams{
        beta,
        zplus,
        zminus
    }
}



fn theta_coefficients(lz:f64, e:f64, c:f64) ->[f64;3]{ // coefficients for 0th, cos()^2 and cos()^4
    let a2 = 1.0;
    let a1 = -1.0*(
        c+lz*lz+A*A*(1.0-e*e)
    )/(A*A*(1.0-e*e));
    let a0= c/(A*A*(1.0-e*e));
    [a0,a1,a2]
}
//Hellsing et al (2006)
fn radial_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{
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
