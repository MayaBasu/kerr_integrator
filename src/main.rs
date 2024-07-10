use roots::{Roots, find_roots_quartic, SearchError};
use egui_plot::{Line, Plot, PlotPoints};
use egui;
use plotters::prelude::*;


//no hair
const M: f64 = 1.0;
const A: f64 = 0.999;
const E:f64 = 0.8889917687523382;
const LZ:f64 =  1.92511617871112;
const C: f64 = 1.2876760482236116;
const DT :f64 = 0.01;
const NUM_ITERATIONS: i32 = 1000;
const MOVE_SCALE: f32 = 0.01;
const SCROLL_SCALE: f32 = 0.001;


enum Coordinates {  //the two types of coordinates which can be integrated by themselves (wrt Mino times)
    Radial,
    Theta,
}
fn get_theta_poly_coefficients( lz:f64, e:f64, c:f64)->[f64;5]{ // coefficients for 0th, cos()^2 and cos()^4
    let a0 = c;
    let a2 = -(c+A.powi(2)*(1.0-e.powi(2))+lz.powi(2));
    let a4 = A.powi(2)*(1.0-e.powi(2));
    [a0,0.0,a2,0.0,a4]
}
fn get_radial_poly_coefficients( lz:f64, e:f64, c:f64)->[f64;5]{
    let a4 = e.powi(2)-1.0;
    let a3 = 2.0*M;
    let a2 = e.powi(2)*2.0*A.powi(2)-2.0*A*lz*e-(A*e-lz).powi(2)-c-A.powi(2);
    let a1 = 2.0*M*(A*e-lz).powi(2)+2.0*M*c;
    let a0 = e.powi(2)*A.powi(4)-2.0*A*lz*e*A.powi(2)+A.powi(2)*lz.powi(2)-A.powi(2)*(A*e-lz).powi(2)-A.powi(2)*c;
    //  println!("{:?}",[a0,a1,a2,a3,a4]);
    [a0,a1,a2,a3,a4]
}
fn r_derivative(x:f64, coefficients:[f64;5]) ->f64{ //dr/dlambda
    coefficients_to_poly(x,coefficients).abs().sqrt()
}
fn theta_derivative(x:f64, coefficients:[f64;5]) ->f64{ //dtheta/dlambda
    let c = x.cos();
    coefficients_to_poly(c,coefficients).abs().sqrt()/x.sin()
}
fn coefficients_to_poly(x:f64, coefficients:[f64;5])-> f64{
    let mut sum = 0.0;
    for i in 0..coefficients.len() {
        sum += coefficients[i]*x.powi(i as i32)
    }
    sum
}

fn root_hunt_and_peck<F: Fn(f64)->f64>(y_start: f64, graph: F) -> (f64, f64) { //find a root above and below a starting value given a graph with the desired roots
    let mut lowerbound = y_start;
    let mut upperbound = y_start;

    let lower_root = loop {

        match roots::find_root_brent(lowerbound,y_start+0.01,& graph,&mut 0.001){
            Ok(root) => { break root}
            Err(message) => {} //println!("{} not found at {}",message,lowerbound)
        };
        lowerbound += -0.01;
        if lowerbound < 0.0{
            panic!("No lower root found")
        }
    };

    let upper_root = loop {
        if E > 0.9{
            upperbound += 1.0;
        }
        else {
            upperbound += 0.01;
        }

        match roots::find_root_brent(y_start,upperbound,& graph,&mut 0.001){
            Ok(root) => { break root}
            Err(message) => {}
        };
        if upperbound > 2000.0{
            panic!("No upper root found below search criteria");
            break 1000.0
        }
    };
   // println!("lower {} upper {}",lower_root,upper_root);
    (lower_root,upper_root)

}


fn integrate(y_start:f64, coordinate:Coordinates,lz:f64,e:f64,c:f64)->PlotPoints{ //function should

    let coefficients:  [f64;5] = match &coordinate  {
        Coordinates::Radial => get_radial_poly_coefficients(lz,e,c),
        Coordinates::Theta => get_theta_poly_coefficients(lz,e,c),
    };

    let derivative= match &coordinate  {  // pick out derivative function
        Coordinates::Radial => r_derivative,
        Coordinates::Theta =>  theta_derivative,
    };

    let (y_min,y_max)= match &coordinate {
        Coordinates::Radial => {(2.0,6.0)}
        Coordinates::Theta => {(1.0471975511965979,2.094395102365872)}
    };


    let roots = match &coordinate {
        Coordinates::Radial => {
          //  println!("FOR Radial");
            let graph = |x: f64| -> f64 {
                coefficients_to_poly(x, coefficients)
            };
            let roots = root_hunt_and_peck(y_start,& graph);
            roots
        }
        Coordinates::Theta => {
           // println!("FOR THETA");
            let graph = |x: f64| -> f64 {
                coefficients_to_poly(x.cos(), coefficients)
            };
            let roots = root_hunt_and_peck(y_start,&graph);
            roots
        }
    };


    let mut going_up:bool = true;
    let mut last_switch_location:f64 = 0.0;
    let mut y = y_start;
   // println!("{}   {}",y_min,y_max);


    (0..NUM_ITERATIONS).map(|i| {

        let x = i as f64 * DT;

        let increment = {if going_up {derivative(y,coefficients)} else {-derivative(y,coefficients)}}*DT;
        y += increment;

        if ((y-y_min).abs() < 0.001 || (y-y_max).abs() < 0.001) && (x - last_switch_location).abs() > 1.0{

            going_up = !going_up;
            last_switch_location = x;
        //    println!("Switched at {}, going up is {}",last_switch_location,going_up)
        }
        [x, y]
    }).collect()

}

fn find_phi(r_plot_points: PlotPoints, t_plot_points: PlotPoints,   lz:f64, e:f64) -> (PlotPoints,PlotPoints,PlotPoints,Vec<(f64,f64,f64)>) { //->PlotPoints
    let r_points = r_plot_points.points().to_vec();
    let t_points = t_plot_points.points().to_vec();
    let mut running_phi = 0.0;
    let coords: Vec<_> = (0..r_points.len()).map(|i|{
        let r = r_points[i].y;
        let l = r_points[i].x;
        let theta = t_points[i].y;
        let phi = phi_total(theta,r,lz,e);
        (l, r,theta,phi)
    }).map(|(l,r,theta,phi_der)| {
        running_phi += phi_der*DT;
        (l,r,theta,running_phi, r*theta.sin()*running_phi.cos(),r*theta.sin()*running_phi.sin(),r*theta.cos())
    }
    ).collect();

    let flat_r = (0..r_points.len()).map(|i|{
        [coords[i].0,coords[i].1]
    }).collect();
    let flat_t = (0..r_points.len()).map(|i|{
        [coords[i].0,coords[i].2]
    }).collect();
    let flat_p = (0..r_points.len()).map(|i|{
        [coords[i].0,coords[i].3 % 2.0*std::f64::consts::PI]
    }).collect();

    let cartesian =
        (0..r_points.len()).map(|i|{
            (coords[i].4,coords[i].5,coords[i].6)
        }).collect();

    (flat_r,flat_t,flat_p,cartesian)
}

#[derive(Default)]

struct Graph {
    chart_pitch: f32,
    chart_yaw: f32,
    chart_scale: f32,
    chart_pitch_vel: f32,
    chart_yaw_vel: f32,
    data: (PlotPoints,PlotPoints,PlotPoints,Vec<(f64,f64,f64)>),

}

impl Graph {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let context = &cc.egui_ctx;

        context.tessellation_options_mut(|tess_options| {
            tess_options.feathering = false;
        });

        context.set_visuals(egui::Visuals::light());
        let radial = integrate(2.142857142857143, Coordinates::Radial, LZ, E, C);
        let angular = integrate(1.0471975511965979, Coordinates::Theta,LZ,E,C);
        let data =find_phi(radial,angular,LZ,E);

        Self {
            chart_pitch: 0.3,
            chart_yaw: 0.9,
            chart_scale: 0.9,
            chart_pitch_vel: 0.0,
            chart_yaw_vel: 0.0,
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

            let width = 10.0;

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


            let radial = integrate(3.0, Coordinates::Radial, LZ, E, C);
            let angular = integrate(1.5, Coordinates::Theta,LZ,E,C);
              let data =find_phi(radial,angular,LZ,E);


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
                let coeffs = get_radial_poly_coefficients(LZ,E,C);





                //let line3 = Line::new(self.data.2);
                let linedata = (0..100).map(|i| r_derivative(i as f64/0.1,coeffs));


               // Plot::new("my_plot").view_aspect(1.0).show(ui, |plot_ui| { plot_ui.line(line3);});

                if ui.button("Quit").clicked() {
                    std::process::exit(0);
                };
            }
        );

        std::thread::sleep(std::time::Duration::from_millis(100));
        ctx.request_repaint();
    }
}

fn main() -> eframe::Result<()> {
    let coeffs = get_theta_poly_coefficients(LZ,E,C);
    println!("{}",theta_derivative(2.9,coeffs));
    let pltpnts = integrate(1.0471975511965979,Coordinates::Theta,LZ,E,C);

    println!("{:?}",pltpnts.points());




    let native_options = eframe::NativeOptions {

        viewport: egui::ViewportBuilder::default().with_inner_size((800.0, 800.0)),
        ..eframe::NativeOptions::default()
    };
    eframe::run_native(
        "Visualizer",
        native_options,
        Box::new(|cc| Ok(Box::new(Graph::new(cc)))),
    )

}


fn delta(r:f64) -> f64 {
    r.powi(2) - 2.0 * M * r + A.powi(2)
}
fn sigma(r:f64, c:f64) -> f64 {
    r.powi(2) + A.powi(2)*c.powi(2)
}
fn p(r:f64, lz:f64,e:f64) -> f64 {
    e*(r.powi(2) + A.powi(2))-A*lz
}
fn radial(r:f64, lz:f64, e:f64, c:f64) -> f64{
    p(r,lz,e).powi(2)-delta(r)*(r.powi(2) + (A*e-lz).powi(2)+c)

    //   (p(r,lz,e).powi(2)-delta(r)*(r.powi(2) + (A*e-lz).powi(2)+c))
}

fn phi_r(r:f64, lz:f64,e:f64) -> f64{
    let d = delta(r);
    let p = p(r,lz,e);
    A*p/d
}

fn phi_theta(theta:f64,lz:f64) -> f64 {
    lz/(1.0-theta.cos().powi(2))
}

fn phi_total(theta:f64,r:f64, lz:f64,e:f64) -> f64 {
    phi_r(r,lz,e) + phi_theta(theta,lz) -A*e
}


