
use crate::functions::{r_derivative_magnitude, radial_coefficients, theta_coefficients, theta_derivative};



const RADIAL_DERIVATIVE_GRANULATION: u32 = 100;
const THETA_DERIVATIVE_GRANULATION: u32 = 1000;

use serde_json::Result;
