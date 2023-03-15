use std::io::Write;
use std::fs::File;

/*
Program follows tracer particles, so dr/dt = v, and v = u B, where u is just a scalar to convert from tesla to m/s.
So the tracer always follows a field line, this can then be used to follow said field lines, and to create a point-care plot.
Updating the position is done using simple Runge-Kutta 4.

In addition, the following assumption is made: the field lines of interest wrap around in toroidal direction without reversing.
So B_phi is always positive, or always negative, but never passing through zero. This allows us to parameterize the variables in
terms of toroidal angle phi, instead of time, which is more convenient, as we can do fixed stepsize in phi.
So we get multiple poloidal planes for the price of one, and no interpolation is needed to get a precise phi angle.

So dphi/dt = omega = const. And v_phi = u B_phi = omega R, so u = omega R/B_phi.
This means that dR/dphi = R B_R/B_phi, dZ/dphi = R B_Z/B_phi.
This approach fails very close to the coils because there the field loops back in a small circle around the wires (think of a standard solenoid).
However, we are not interested in the ripple anyway, so limiting ourselves to the main plasma.
*/
fn main() {
    //inputs for the magnetic field
    let input = MagInputs {
        minor_rad0: 0.5,
        major_rad0: 1.0,
        z_pos_axis: 0.0,
        mag_field0: 4.5,
    };
    let mut position = [0.7, 0.0, 0.0]; //cylindrical coordinates are in (R,phi,Z) format
    let num_revolutions: u32 = 100; //except the last point, so each poloidal plane will have num_revolutions points.
    let steps_per_rev: u32 = 1000; //number of steps per revolution todo: add parameter to choose how many of these are saved, right now it saves every step...

    let step_size = 2.0 * std::f64::consts::PI / (steps_per_rev as f64);
    let num_steps = num_revolutions * steps_per_rev;
    let num_elements = (2 * num_revolutions * steps_per_rev) as usize;//times 2 because it stores R and Z: (R0,Z0,R1,Z1,...RN-1,ZN-1) of a field line. To create a poincare plot, use stride=steps_per_rev to cut out just the points of one poloidal plane.
    let mut poincare_points = vec![0.0; num_elements]; //todo: use multiple tracers (10 or so), not just one!

    store_point(&position, &mut poincare_points, 0);

    println!("Starting the simulation");
    for step in 1..num_steps {
        //advance by one time step
        update_position(&mut position, &input, step_size);

        //save data
        store_point(&position, &mut poincare_points, step);
    }

    println!("Completed!");
    println!("Final position: R={}, ϕ={}, Z={}", position[0], position[1], position[2]);

    let mut file_handler = File::create("output.csv").expect("Unable to create file");
    for elem in poincare_points.iter() {
        write!(file_handler, "{}\n", elem).expect("Error writing line to output file");
    }
}

struct MagInputs {
    minor_rad0: f64,    //minor radius
    major_rad0: f64,    //major radius
    z_pos_axis: f64, //z position of magnetic axis
    mag_field0: f64,    //on axis field
}

fn update_position(position: &mut [f64; 3], input: &MagInputs, angle_step: f64) {
    let mut temporary_position = [0.0; 3];
    let mut k1 = [0.0; 3];
    let mut k2 = [0.0; 3];
    let mut k3 = [0.0; 3];
    let mut k4 = [0.0; 3];

    let mut derivative = get_derivative(position, input);
    for i in 0..3 {
        k1[i] = angle_step * derivative[i];
        temporary_position[i] = position[i] + 0.5 * k1[i];
    }
    derivative = get_derivative(&temporary_position, input);
    for i in 0..3 {
        k2[i] = angle_step * derivative[i];
        temporary_position[i] = position[i] + 0.5 * k2[i];
    }
    derivative = get_derivative(&temporary_position, input);
    for i in 0..3 {
        k3[i] = angle_step * derivative[i];
        temporary_position[i] = position[i] + k3[i];
    }
    derivative = get_derivative(&temporary_position, input);
    for i in 0..3 {
        k4[i] = angle_step * derivative[i];
        position[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
}

///Simple torus for now, returns (B_R,B_phi,B_Z)
fn get_magnetic_field(position: &[f64; 3], input: &MagInputs) -> [f64; 3] {
    let mut mag_field = [0.0; 3];

    let major_radius = position[0];
    let reduced_major = major_radius - input.major_rad0;
    let reduced_z = position[2] - input.z_pos_axis;
    let minor_radius2 = reduced_major * reduced_major + reduced_z * reduced_z;

    let b_toroidal = input.mag_field0 * input.major_rad0 / major_radius;
    // safety factor q = r*Bphi/(R*Btheta)
    // q = 4 (r/a)^2 + 1 = r*bToroidal/(R*bTheta), solve for bTheta, then project to x,y,z.
    let safety_factor = 1.0 + 4.0 * minor_radius2 / (input.minor_rad0 * input.minor_rad0);
    let b_poloidal_over_r = b_toroidal / (major_radius * safety_factor); //bTheta/r

    mag_field[0] = -reduced_z * b_poloidal_over_r;
    mag_field[1] = b_toroidal;
    mag_field[2] = reduced_major * b_poloidal_over_r;
    mag_field
}

fn get_derivative(position: &[f64; 3], input: &MagInputs) -> [f64; 3] {
    let mag_field = get_magnetic_field(position, input);
    let mut derivative = [0.0; 3];
    derivative[0] = position[0] * mag_field[0] / mag_field[1]; // dR/dphi = R*B_R/B_phi
    derivative[1] = 1.0; //dphi/dphi = 1
    derivative[2] = position[0] * mag_field[2] / mag_field[1]; // dZ/dphi = R*B_Z/B_phi
    derivative
}

fn store_point(position: &[f64; 3], destination: &mut [f64], step: u32) {
    let index = (step * 2) as usize;
    destination[index] = position[0];
    destination[index + 1] = position[2];
}