use std::sync::Arc;
use std::fs::File;
use std::io::{self, Read};
use array_programming::complex;
use array_programming::accumulate;
use array_programming::view;

fn initialize_buffer(file_path: &str) -> io::Result<Arc<[complex::Complex<f32>]>> {
    let mut file = File::open(file_path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    
    // Convert the byte buffer into Complex<f32> buffer
    // First, ensure the buffer length is a multiple of Complex<f32> size
    if buffer.len() % std::mem::size_of::<complex::Complex<f32>>() != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Buffer size is not a multiple of Complex<f32> size"
        ));
    }
    
    // Create a new vector of Complex<f32>
    let complex_buffer = unsafe {
        // Cast the bytes to Complex<f32> using std::slice::from_raw_parts
        let complex_slice = std::slice::from_raw_parts(
            buffer.as_ptr() as *const complex::Complex<f32>,
            buffer.len() / std::mem::size_of::<complex::Complex<f32>>()
        );
        
        // Clone the slice into a new vector
        Vec::from(complex_slice)
    };
    
    Ok(Arc::from(complex_buffer))
}

use plotters::prelude::*;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let buffer = initialize_buffer("test_dump.txt").expect("Failed to initialize buffer");
    
    println!("Buffer contains {} complex numbers", buffer.len());
    // Access the first complex number
    if !buffer.is_empty() {
            println!("First complex number: {} + {}i", buffer[0].re, buffer[0].im);
        }

    let root = BitMapBackend::new("plotters-doc-data/0.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-0f32..1f32, -1f32..1f32)?;

    chart.configure_mesh().draw()?;

    // TODO : map complex component as well!
    chart
        .draw_series(LineSeries::new(
            buffer.iter().enumerate().map(|(i, x)| (i as f32 / buffer.len() as f32, x.re)),
            &RED,
        ))?
        .label("Real Component of Complex Data")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}
