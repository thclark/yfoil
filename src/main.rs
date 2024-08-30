use yfoil::geometry::read_geometry_from_file;

use clap::Parser;
use plotly::{Plot, Scatter};

/// Load aerofoil geometry and determine polar and boundary layer characteristics
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path of the file containing aerofoil JSON
    #[arg(long, default_value = "aerofoil.json")]
    file: String,
}

fn main() {
    let args = Args::parse();
    let version = env!("CARGO_PKG_VERSION");

    // Header lines
    println!("----------------------------------------");
    println!("yfoil version {}", version);
    println!("    -- why? because y comes after x...");
    println!("----------------------------------------");

    // Load and parse the aerofoil geometry
    println!("Reading aerofoil input file at {}", args.file);
    let geometry = read_geometry_from_file(args.file).unwrap();

    // Show where the reference is
    println!(
        "Reference x/c, y/c : {}, {}",
        geometry.reference[0], geometry.reference[1]
    );

    // Make a plotly figure

    let mut plot = Plot::new();
    let trace = Scatter::new(geometry.x_c, geometry.y_c);
    plot.add_trace(trace);
    plot.write_html("aerofoil_geometry.html");
}
