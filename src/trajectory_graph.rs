use egui_plot::{PlotPoints,Line,Plot};
use crate::constants::*;
use crate::functions::*;
use crate::numeric_solvers::*;
use egui;
use plotters::prelude::*;


#[derive(Default)]
pub(crate) struct Graph {
    chart_pitch: f32,
    chart_yaw: f32,
    chart_scale: f32,
    chart_pitch_vel: f32,
    chart_yaw_vel: f32,

    r_initial: f64,
    r_min: f64,
    r_max: f64,

    theta_initial: f64,
    theta_min: f64,
    theta_max: f64,

    data: (PlotPoints,PlotPoints,PlotPoints,Vec<(f64,f64,f64)>),

}
impl Graph {
    pub(crate) fn new(cc: &eframe::CreationContext<'_>, vals:((f64, f64, f64), (f64, f64, f64))) -> Self {
        let context = &cc.egui_ctx;

        context.tessellation_options_mut(|tess_options| {
            tess_options.feathering = false;
        });

        context.set_visuals(egui::Visuals::light());
        let radial = integrate(vals.0.0,vals.0.1,vals.0.2, Coordinates::Radial, LZ, E, C);
        let angular = integrate(vals.1.0,vals.1.1,vals.1.2, Coordinates::Theta,LZ,E,C);
        let data =find_phi(radial,angular,LZ,E);

        println!("{:?}",data.0.points());
        println!("r initial {}, r min {}, r max {}",vals.0.0,vals.0.1,vals.0.2);
        println!("theta initial {}, theta min {}, theta max {}",vals.1.0,vals.1.1,vals.1.2);

        Self {
            chart_pitch: 0.3,
            chart_yaw: 0.9,
            chart_scale: 0.9,
            chart_pitch_vel: 0.0,
            chart_yaw_vel: 0.0,
            r_initial:vals.0.0,
            r_min:vals.0.1,
            r_max:vals.0.2,
            theta_initial:vals.1.0,
            theta_min:vals.1.1,
            theta_max:vals.1.2,
            data,
        }

    }
}
impl eframe::App for Graph {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_pixels_per_point(1.5);
        egui::CentralPanel::default().show(ctx, |ui| {


            let (pitch_delta, yaw_delta, scale_delta) = ui.input(|input| {
                let pointer = &input.pointer;
                let delta = pointer.delta();

                let (pitch_delta, yaw_delta) = match pointer.primary_down() {
                    true => (delta.y * MOVE_SCALE, -delta.x * MOVE_SCALE),
                    false => (self.chart_pitch_vel, self.chart_yaw_vel),
                };

                let scale_delta = input.smooth_scroll_delta.y * SCROLL_SCALE;

                (pitch_delta, yaw_delta, scale_delta)
            });

            self.chart_pitch_vel = pitch_delta;
            self.chart_yaw_vel = yaw_delta;
            self.chart_pitch += self.chart_pitch_vel;
            self.chart_yaw += self.chart_yaw_vel;
            self.chart_scale += scale_delta;

            let root = egui_plotter::EguiBackend::new(ui).into_drawing_area();

            root.fill(&WHITE).unwrap();

            let width = 20.0;

            let x_axis = (-width..width).step(0.1);
            let z_axis = (-width..width).step(0.1);

            let mut chart = ChartBuilder::on(&root)
                .caption(format!("Kerr Geodesic in Cartesian"), (FontFamily::SansSerif, 20))
                .build_cartesian_3d(x_axis, -width..width, z_axis)
                .unwrap();

            chart.with_projection(|mut pb| {
                pb.yaw = self.chart_yaw as f64;
                pb.pitch = self.chart_pitch as f64;
                pb.scale = self.chart_scale as f64;
                pb.into_matrix()
            });

            chart
                .configure_axes()
                .light_grid_style(BLACK.mix(0.15))
                .max_light_lines(3)
                .draw()
                .unwrap();



            chart
                .draw_series(LineSeries::new(
                    self.data.3.clone(),
                    & RED,
                ))
                .unwrap()
                .label("Line");


            root.present().unwrap();
        });
        egui::SidePanel::left("left pannel").show(ctx, |ui| {

            ui.heading("In Mino Time");


            let radial = integrate(self.r_initial, self.r_min,self.r_max,Coordinates::Radial, LZ, E, C);
            let angular = integrate(self.theta_initial, self.theta_min, self.theta_max, Coordinates::Theta, LZ, E, C);
            let data =find_phi(radial, angular, LZ, E);


            let radial_line = Line::new(data.0);
            let theta_line = Line::new(data.1);
            let aziuthal_line  = Line::new(data.2);


            Plot::new("my_plot").view_aspect(2.0).show(ui, |plot_ui|{
                plot_ui.line(radial_line);
                plot_ui.line(theta_line);
                plot_ui.line(aziuthal_line)
            } );


            if ui.button("Unicorns").clicked() {
                std::process::exit(0);
            };
        });
        egui::SidePanel::right("side pannel").show(ctx,
                                                   |ui| {
                                                       ui.heading("Distance between streams");

                                                       ui.label("Graph of closest distance of approach");



                                                       let coeffs = radial_coefficients(LZ,E,C);

                                                       //let line3 = Line::new(self.data.2);
                                                       let linedata: PlotPoints = (0..1000).map(|i| [i as f64/0.1,r_derivative(i as f64/0.1,coeffs)]).collect();
                                                       let line3 = Line::new(linedata);


                                                       Plot::new("my_plot").view_aspect(1.0).show(ui, |plot_ui| { plot_ui.line(line3);});

                                                       if ui.button("Quit").clicked() {
                                                           std::process::exit(0);
                                                       };
                                                   }
        );
        println!("requesting repaint");

        std::thread::sleep(std::time::Duration::from_millis(10));
        ctx.request_repaint();
    }
}