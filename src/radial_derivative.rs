use egui_plot::{Plot,PlotPoints,Line};
use crate::constants::{E, LZ, C};
use crate::functions::{r_derivative,radial_coefficients};
use std::io;


const RADIAL_DERIVATIVE_GRANULATION: u32 = 100;


#[derive(Default)]
pub(crate) struct RadialGraph {

    upper_bound:u32,

    data:Vec<[f64;2]>,

}
impl RadialGraph {
    pub(crate) fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let context = &cc.egui_ctx;

        context.tessellation_options_mut(|tess_options| {
            tess_options.feathering = false;
        });

        context.set_visuals(egui::Visuals::light());
        let mut upper_bound = String::new();
        println!("Upper Bound of Graph: (u32)");
        io::stdin()
            .read_line(&mut upper_bound)
            .expect("Failed to read line");
        let upper_bound: u32 = upper_bound.trim().parse().expect("Please type a u32");
        let coefficients = radial_coefficients(LZ,E,C);
        println!("The radial coefficients are");
        println!("{:?}",coefficients);
        let data:Vec<[f64;2]> = (0..upper_bound*RADIAL_DERIVATIVE_GRANULATION).map(|i|
        [(i as f64)/(RADIAL_DERIVATIVE_GRANULATION as f64),r_derivative((i as f64)/(RADIAL_DERIVATIVE_GRANULATION as f64),coefficients)]
        ).collect();
        Self {
            upper_bound,
            data,
        }



    }
}

impl eframe::App for RadialGraph {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_pixels_per_point(1.5);
        egui::CentralPanel::default().show(ctx, |ui| {
            let radial_derivative = Line::new(PlotPoints::new(self.data.clone()));

            Plot::new("Radial Derivative").view_aspect(1.0).show(ui, |plot_ui| { plot_ui.line(radial_derivative);});



        });

        std::thread::sleep(std::time::Duration::from_millis(10));
        ctx.request_repaint();
    }
}