
//no hair
pub const M: f64 = 1.0;
pub const A: f64 = 0.99;

//manually adjusted constants
pub const E:f64 =0.9999;

//derived constants

// code uses LZ,C,E


pub const L:f64 = 6.5;
pub const COS_I: f64 = 0.5;

pub const SIN_I:f64 = 1.0-COS_I*COS_I;


pub const  LZ:f64 = SIN_I*L;
pub const C: f64 = L*L-LZ*LZ;



pub const K:f64 = (LZ-A*E).powi(2)+C; //Batra paper


//pub const L:f64 = 6.862943974709396;

//pub const L:f64 = 6.0;


// code uses LZ,C,E


//totoal is 7, so LZ = sqrt(36-9)