# Poincare plotter in rust

Simple Rust code to map out the magnetic flux surfaces in tokamak/stellarator geometry.
Tracer particles follow the field lines, and a  point is placed every time one toroidal revolution is completed.

# Prerequisites

* Rust (for the main program)
* Python 3 (for plotting the data)

# Installation

Simply run `cargo build --release` in the same directory as this README file.

# Running the program

`./target/release/magnetic_field`
Then analyse the output data with python.

# Misc

Disclaimer: my first Rust project.

Mike 2023