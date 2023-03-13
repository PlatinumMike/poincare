use std::io::Write;
use std::fs::File;

/// Program follows tracer particles, so dr/dt = v, and v = c B, where c is just a constant to convert from tesla to m/s, e.g. c = 1 (m/s)/T for convenience.
/// So the tracer always follows a field line, this can then be used to follow said field lines, and to create a point-care plot.
/// Updating the position is done using simple Runge-Kutta 4.
fn main() {
    //inputs for the magnetic field
    let input = MagInputs {
        minor_rad0: 1.0,
        major_rad0: 3.0,
        z_pos_axis: 0.0,
        mag_field0: 4.5,
    };
    let mut position = [2.7, 0.0, 0.0];
    let max_num_steps: u32 = 1000_000_000; //to avoid infinite loop
    let time_step = 1.0e-4;
    let num_revolutions: u32 = 100; //so you will have num_revolutions+1 points
    let num_elements: usize = 2 * (num_revolutions + 1) as usize;//times 2 because it stores R and Z: (R0,Z0,R1,Z1,...RN-1,ZN-1).

    let mut poincare_points = vec![0.0; num_elements];

    let mut revolution :u32 = 0;
    let mut step :u32 = 0;
    let angle0 = get_angle(&position); //todo: modify check for wrap around if you place the plane at a different position.
    let mut angle_old = angle0;
    let mut angle_new = angle0;
    store_point(&position,& mut poincare_points, revolution); //start already at phi=0, so store a point.

    println!("Starting the simulation");
    while revolution < num_revolutions {
        //interrupt loop if it gets stuck for some reason.
        if step >= max_num_steps {
            println!("Infinite loop, terminating...");
            break;
        }

        //advance by one time step
        update_position(&mut position, &input, time_step);
        step += 1;
        angle_new = get_angle(&position);

        if wrap_around(angle_old, angle_new) {
            revolution += 1;
            store_point(&position, &mut poincare_points, revolution);
        }
        angle_old = angle_new;
    }

    println!("Completed, number of steps used: {}", step);
    println!("Final position: x={}, y={}, z={}", position[0], position[1], position[2]);

    let mut file_handler = File::create("output.csv").expect("Unable to create file");
    for i in 0..num_revolutions + 1 {
        let index = 2 * i as usize;
        write!(file_handler, "{},{}\n", poincare_points[index], poincare_points[index + 1]);
    }
}

struct MagInputs {
    minor_rad0: f64,    //minor radius
    major_rad0: f64,    //major radius
    z_pos_axis: f64, //z position of magnetic axis
    mag_field0: f64,    //on axis field
}

fn update_position(position: &mut [f64; 3], input: &MagInputs, dt: f64) {
    let mut temporary_position = [0.0; 3];
    let mut k1 = [0.0; 3];
    let mut k2 = [0.0; 3];
    let mut k3 = [0.0; 3];
    let mut k4 = [0.0; 3];

    let mut velocity = get_magnetic_field(position, input); //c=1
    for i in 0..3 {
        k1[i] = dt * velocity[i];
        temporary_position[i] = position[i] + 0.5 * k1[i];
    }
    velocity = get_magnetic_field(&temporary_position, input);
    for i in 0..3 {
        k2[i] = dt * velocity[i];
        temporary_position[i] = position[i] + 0.5 * k2[i];
    }
    velocity = get_magnetic_field(&temporary_position, input);
    for i in 0..3 {
        k3[i] = dt * velocity[i];
        temporary_position[i] = position[i] + k3[i];
    }
    velocity = get_magnetic_field(&temporary_position, input);
    for i in 0..3 {
        k4[i] = dt * velocity[i];
        position[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
}

///Simple torus for now
fn get_magnetic_field(position: &[f64; 3], input: &MagInputs) -> [f64; 3] {
    let mut mag_field = [0.0; 3];

    let major_radius = get_major_radius(position);
    let reduced_major = major_radius - input.major_rad0;
    let reduced_z = position[2] - input.z_pos_axis;
    let minor_radius2 = reduced_major * reduced_major + reduced_z * reduced_z;

    let b_toroidal = input.mag_field0 * input.major_rad0 / major_radius;
    // safety factor q = r*Bphi/(R*Btheta)
    // q = 4 (r/a)^2 + 1 = r*bToroidal/(R*bTheta), solve for bTheta, then project to x,y,z.
    let safety_factor = 1.0 + 4.0 * minor_radius2 / (input.minor_rad0 * input.minor_rad0);
    let b_poloidal_over_r = b_toroidal / (major_radius * safety_factor); //bTheta/r

    mag_field[0] = (-position[1] * b_toroidal - position[0] * reduced_z * b_poloidal_over_r) / major_radius;
    mag_field[1] = (position[0] * b_toroidal - position[1] * reduced_z * b_poloidal_over_r) / major_radius;
    mag_field[2] = reduced_major * b_poloidal_over_r;
    mag_field
}

fn get_angle(position: &[f64; 3]) -> f64 {
    position[1].atan2(position[0])
}

fn get_major_radius(position: &[f64; 3]) -> f64 {
    let r2 = position[0] * position[0] + position[1] * position[1];
    r2.sqrt()
}

fn get_vertical_position(position: &[f64; 3]) -> f64 {
    position[2]
}

fn store_point(position: &[f64; 3], destination: &mut [f64], revolution : u32) {
    //todo: instead of storing it directly after a wrap_around has occured, it would be better to interpolate between this and the last position to find the most accurate point at the desired poloidal plane. Using very small time steps would also work.
    let index = (revolution*2) as usize;
    destination[index] = get_major_radius(position);
    destination[index + 1] = get_vertical_position(position);
}

fn wrap_around(angle0: f64, angle1: f64) -> bool {
    //angle is in (-pi,pi], so if it goes beyond pi, it jumps to -pi, so angle1<angle0.
    //but if the particles is moving in direction of decreasing angle, then you should check for angle1>angle0.
    //these can be combined into one check, if the sign flips we have wrapped around, or equivalently, angle0*angle1<0
    angle1 * angle0 < 0.0
}