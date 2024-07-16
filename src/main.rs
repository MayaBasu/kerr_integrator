use roots::{Roots, find_roots_quartic};
use egui_plot::{Line, Plot, PlotPoints};
use egui;
use plotters::prelude::*;
use std::io;


//no hair
const M: f64 = 1.0;
const A: f64 = 0.999;

//1.0471975511965979 0.0 2.1428571428571463
//1.0471975511965979 0.0 1.5789473684210527
// orbital constants

const E:f64 =0.9693269909705301;
const LZ:f64 = 2.074219118310137;
const C: f64 = 1.4491994255741059;

//const E:f64 = 0.8772863839709705;
//const LZ:f64 = 1.9034666944312097;
//const C: f64 = 1.2652055090196432;


//step size and number of steps to take during integration
const DT :f64 = 0.001;
const NUM_ITERATIONS: i32 = 10000;

//sensitivity to mouse moving
const MOVE_SCALE: f32 = 0.01;
const SCROLL_SCALE: f32 = 0.001;


enum Coordinates {  //the two types of coordinates which can be integrated by themselves (wrt Mino times)
    Radial,
    Theta,
}
#[derive(Default)]
struct Graph {
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
    fn new(cc: &eframe::CreationContext<'_>, vals:((f64,f64,f64),(f64,f64,f64))) -> Self {
        let context = &cc.egui_ctx;

        context.tessellation_options_mut(|tess_options| {
            tess_options.feathering = false;
        });

        context.set_visuals(egui::Visuals::light());
        let radial = integrate(vals.0.0,vals.0.1,vals.0.2, Coordinates::Radial, LZ, E, C);
        let angular = integrate(vals.1.0,vals.1.1,vals.1.2, Coordinates::Theta,LZ,E,C);
        let data =find_phi(radial,angular,LZ,E);

          println!("{:?}",data.1.points());
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
            let angular = integrate(self.theta_initial, self.theta_min,self.theta_max,Coordinates::Theta,LZ,E,C);
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





                                                       //let line3 = Line::new(self.data.2);
                                                       //let linedata = (0..100).map(|i| r_derivative(i as f64/0.1,coeffs));


                                                       // Plot::new("my_plot").view_aspect(1.0).show(ui, |plot_ui| { plot_ui.line(line3);});

                                                       if ui.button("Quit").clicked() {
                                                           std::process::exit(0);
                                                       };
                                                   }
        );

        std::thread::sleep(std::time::Duration::from_millis(10));
        ctx.request_repaint();
    }
}

fn main() -> eframe::Result<()> {


    let native_options = eframe::NativeOptions {

        viewport: egui::ViewportBuilder::default().with_inner_size((1600.0, 800.0)),
        ..eframe::NativeOptions::default()
    };
    eframe::run_native(
        "Visualizer",
        native_options,
        Box::new(|cc| Ok(Box::new(Graph::new(cc, setup(radial_coefficients(LZ, E, C), theta_coefficients(LZ, E, C),true))))),
    )

}
fn setup(r_coefficients:[f64;5], theta_coefficients:[f64;5], simple: bool) -> ((f64, f64, f64), (f64, f64, f64)) { //((r start, rmin, rmax),(theta start,theta min,theta max))


    let L_sqr = C + LZ.powi(2);

    let beta = (A.powi(2))*(1.0 - E.powi(2));


    let z_plus = ((L_sqr + beta) + ((L_sqr + beta).powi(2) - 4.0*C*beta).sqrt())/(2.0*beta);
    let z_minus = ((L_sqr + beta) - ((L_sqr + beta).powi(2) - 4.0*C*beta).sqrt())/(2.0*beta);

   println!("z plus and minus {}   {}",z_plus,z_minus);






    let multiple_roots = find_roots_quartic(r_coefficients[4], r_coefficients[3], r_coefficients[2], r_coefficients[1], r_coefficients[0]);
    let mut root_list: Vec<f64> = Vec::new();

    match multiple_roots {
        Roots::Four(roots) =>{
            for root in roots{
                if root >0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::Three(roots) =>{
            for root in roots{
                if root >0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::Two(roots) =>{
            for root in roots{
                if root >0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::One(roots) =>{
            for root in roots{
                if root >0.0{
                    root_list.push(root);
                }
            }
        }
        Roots::No(roots)=>{

            for root in roots{
                if root >0.0{
                    root_list.push(root);

                }
            }
        }
    }
    println!("The positive roots of the radial derivative function are at {:?}",root_list);
    if root_list.len() < 2{
        panic!("There were less than 2 roots :( ")
    }

    for i in 0..root_list.len()-1 {
       // println!("sample intermediate value: {}",coefficients_to_poly((root_list[i+1]-root_list[i])/2.0+root_list[i],r_coefficients));
    }
    let (r_start,r_min,r_max)  = if simple {
        println!("Implementing simple setup...");
        let r_start = root_list[root_list.len() -2];
        let r_min = root_list[root_list.len()-2];
        let r_max = root_list[root_list.len()-1];
        println!("Chose r initial {}, r min {} r max {}",r_start,r_min,r_max);
        (r_start,r_min,r_max)

    }else {
        let mut r_start = String::new();
        println!("Select an initial starting value for r:");

        io::stdin()
            .read_line(&mut r_start)
            .expect("Failed to read line");
        let r_start = r_start.trim().parse().expect("Please type a number!");
        let mut r_min = 0.0;
        let mut r_max = 0.0;
        println!("roots {:?}", root_list);
        for root in root_list {
            if root <= r_start {
                r_min = root;
            }
            if root > r_start {
                r_max = root;
                break
            }
        };
        (r_start,r_min,r_max)
    };
    let (theta_start,theta_min,theta_max) = if simple {

        let theta_start =  std::f64::consts::PI/2.0;//1.0471975511965979;
        let graph = |x: f64| -> f64 {
            coefficients_to_poly(x.cos(), theta_coefficients)
        };
        let (theta_min,theta_max) = root_hunt_and_peck(theta_start,graph);
        println!("selecting theta initial{} theta  min {} and theta max {}",theta_start,theta_min,theta_max);

        (theta_start,theta_min,theta_max)

    } else{
        let mut theta_start =String::new();
        println!("Select an initial starting value for theta:");

        io::stdin()
            .read_line(&mut theta_start)
            .expect("Failed to read line");
        let theta_start = theta_start.trim().parse().expect("Please type a number!");

        let graph = |x: f64| -> f64 {
            coefficients_to_poly(x.cos(), theta_coefficients)
        };
        let (theta_min,theta_max) = root_hunt_and_peck(theta_start,graph);

        println!("Then the bounds are {} to {} for theta",theta_min,theta_max);
        (theta_start,theta_min,theta_max)

    };


    ((r_start,r_min,r_max),(theta_start,theta_min,theta_max))


}


fn theta_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{ // coefficients for 0th, cos()^2 and cos()^4
    let a0 = c;
    let a2 = -(c+A.powi(2)*(1.0-e.powi(2))+lz.powi(2));
    let a4 = A.powi(2)*(1.0-e.powi(2));
    [a0,0.0,a2,0.0,a4]
}
fn radial_coefficients(lz:f64, e:f64, c:f64) ->[f64;5]{
    let a4 = e.powi(2)-1.0;
    let a3 = 2.0*M;
    let a2 = e.powi(2)*2.0*A.powi(2)-2.0*A*lz*e-(A*e-lz).powi(2)-c-A.powi(2);
    let a1 = 2.0*M*(A*e-lz).powi(2)+2.0*M*c;
    let a0 = e.powi(2)*A.powi(4)-2.0*A*lz*e*A.powi(2)+A.powi(2)*lz.powi(2)-A.powi(2)*(A*e-lz).powi(2)-A.powi(2)*c;
    //  println!("{:?}",[a0,a1,a2,a3,a4]);
    [a0,a1,a2,a3,a4]
}

fn coefficients_to_poly(x:f64, coefficients:[f64;5])-> f64{
    let mut sum = 0.0;
    for i in 0..coefficients.len() {
        sum += coefficients[i]*x.powi(i as i32)
    }
    sum
}

fn r_derivative(x:f64, coefficients:[f64;5]) ->f64{ //dr/dlambda
    coefficients_to_poly(x,coefficients).abs().sqrt()
}
fn theta_derivative(x:f64, coefficients:[f64;5]) ->f64{ //dtheta/dlambda
    let c = x.cos();
    coefficients_to_poly(c,coefficients).abs().sqrt()/x.sin()
}


fn root_hunt_and_peck<F: Fn(f64)->f64>(y_start: f64, graph: F) -> (f64, f64) { //find a root above and below a starting value given a graph with the desired roots



    let mut lowerbound = y_start;
    let mut upperbound = y_start+0.1;

    let lower_root = {
        if graph(y_start).abs() < 0.001{
            println!("taking lower root");
            y_start
        } else {
            loop {
               // println!("searching between {} and {}",lowerbound,y_start + 0.01);
                match roots::find_root_brent(lowerbound, y_start + 0.01, &graph, &mut 0.001) {
                    Ok(root) => { break root }
                    Err(_message) => {} //println!("{} not found at {}",message,lowerbound)
                };
                lowerbound += -0.01;
                if lowerbound < 0.0 {
                    panic!("No lower root found")
                }
            }
        }
    };


    let upper_root = loop {
        if E > 0.9{
            upperbound += 1.0;
        }
        else {
            upperbound += 0.01;
        }

        match roots::find_root_brent(y_start+0.01,upperbound,& graph,&mut 0.001){
            Ok(root) => { break root}
            Err(_message) => {}
        };
        if upperbound > 2000.0{
            panic!("No upper root found below search criteria");

        }
    };
   // println!("lower {} upper {}",lower_root,upper_root);
    (lower_root,upper_root)

}


fn integrate(y_start:f64, y_min:f64,y_max:f64, coordinate:Coordinates,lz:f64,e:f64,c:f64)->PlotPoints{ //function should

    let coefficients:  [f64;5] = match &coordinate  {
        Coordinates::Radial => radial_coefficients(lz, e, c),
        Coordinates::Theta => theta_coefficients(lz, e, c),
    };

    let derivative= match &coordinate  {  // pick out derivative function
        Coordinates::Radial => r_derivative,
        Coordinates::Theta =>  theta_derivative,
    };

    let ( length_tolerance,time_tolerance)  = match &coordinate{
        Coordinates::Radial => {(0.1,1.0)}
        Coordinates::Theta => {(0.001,0.1)}
    };

 //   let (y_min,y_max)= match &coordinate {
//        Coordinates::Radial => {( 1.6666666666666665,15.000000000000004)}
//        Coordinates::Theta => {(1.0471975511965979,2.0943951023841776)}
 //   };



    let mut going_up:bool = true;
    let mut last_switch_location:f64 = 0.0;

   // println!("alsjdflkejlfjasekfjsd {}",y_start);
    let mut y = y_start;

   // println!("{}   {}",y_min,y_max);


    (0..NUM_ITERATIONS).map(|i| {
        if i % 100 == 0 {
          //  println!("{}, {}",i,y);
        }

        let x = i as f64 * DT;

        let increment = {if going_up {derivative(y,coefficients)} else {-derivative(y,coefficients)}}*DT;
        y += increment;

        if ((y-y_min).abs() < length_tolerance || (y-y_max).abs() < length_tolerance) && (x - last_switch_location).abs() > time_tolerance{

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







fn delta(r:f64) -> f64 {
    r.powi(2) - 2.0 * M * r + A.powi(2)
}
fn p(r:f64, lz:f64,e:f64) -> f64 {
    e*(r.powi(2) + A.powi(2))-A*lz
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


