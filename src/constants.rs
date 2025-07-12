
//no hair
pub const M: f64 = 1.0;
pub const A: f64 = 0.9;

//manually adjusted constants
pub const E:f64 =0.9999;


//input parameters:

/*

pub const  LZ:f64 = 1.9;
pub const C: f64 = 1.26;

//derived constants

// code uses LZ,C,E

//Lz  = 6, Q = 3, so L = 6.24
pub const L:f64 = 2.2068076490713913;
pub const COS_I: f64 = 0.5086520415608135
;

//pub const SIN_I:f64 = (1.0-COS_I*COS_I).sqrt();

 */

pub const LZ: f64 = 0.5*6.5;
pub const C:f64 = 6.5*6.5*(1.0-0.5*0.5);






//pub const L:f64 = 6.862943974709396;

//pub const L:f64 = 6.0;


// code uses LZ,C,E


//totoal is 7, so LZ = sqrt(36-9)