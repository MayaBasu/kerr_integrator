
//no hair
pub const M: f64 = 1.0;
pub const A: f64 = 0.9999;

//manually adjusted constants
pub const E:f64 =0.999;
pub const COS_I:f64 = 0.5;
pub const L:f64 = 4.0;


//derived constants
pub const SIN_I:f64 = 1.0-COS_I*COS_I;
// code uses LZ,C,E
pub const LZ:f64 = SIN_I*L;
pub const C: f64 = L*L-LZ*LZ;

