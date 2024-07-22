
//no hair
pub const M: f64 = 1.0;
pub const A: f64 = 0.99;

//manually adjusted constants
pub const E:f64 =0.95;
pub const COS_I:f64 = 0.5;
pub const L:f64 = 40.0;


// code uses LZ,C,E
pub const LZ:f64 = COS_I*L;
pub const C: f64 = L*L-LZ*LZ;

//step size and number of steps to take during integration
pub const DT :f64 = 0.001;
pub const NUM_ITERATIONS: i32 = 1000;

//sensitivity to mouse moving
pub const MOVE_SCALE: f32 = 0.01;
pub const SCROLL_SCALE: f32 = 0.001;



pub enum Coordinates {  //the two types of coordinates which can be integrated by themselves (wrt Mino times)
    Radial,
    Theta,
}

pub enum Options {
    ViewRadialDerivative,
    Plot,
}