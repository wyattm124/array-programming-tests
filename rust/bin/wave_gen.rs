use array_programming::complex;
use std::path::Path;
use std::fs::File;
use std::io::{self, Read, Write};
use bincode;

// NOTE : phase shift and time are realtive to the frequency
//  all the calculations are thus relative to the wave period.
//  In this wave, phase shift should be in the range [0, 1) and is a proportion
//  of the wave period, and time is in # of periods since time 0.
fn wave_pt(freq: f32, phase_shift: f32, relative_time: usize) -> complex::Complex<f32> {
    let phase = 2.0 * std::f32::consts::PI * freq * ((relative_time as f32) + phase_shift);
    complex::Complex::new(phase.cos(), phase.sin())
}

// Based on the wave_pt that generates a pt of a wave, populate a buffer with a wave
fn wave_gen(freq: f32, phase_shift: f32, data: &mut [complex::Complex<f32>]) -> () {
    for i in 0..data.len() {
        data[i] += wave_pt(freq, phase_shift, i);
    }
}

// Function to write binary data to file
fn write_wave_binary(
    freq: f32,
    phase_shift: f32,
    num_samples: usize,
    output_path: &Path,
) -> Result<usize, bincode::error::EncodeError> {
    let config = bincode::config::standard();
    let mut buffer = vec![complex::Complex::new(0.0, 0.0); num_samples];
    wave_gen(freq, phase_shift, &mut buffer);
    
    let encoded: Vec<u8> = bincode::encode_to_vec(&buffer, config).unwrap();
    let mut file = File::create(output_path).unwrap();
    bincode::encode_into_std_write(&buffer, &mut file, config)
}

fn read_wave_binary(
    input_path: &Path,
) -> Result<Vec<complex::Complex<f32>>, bincode::error::DecodeError> {
    let mut file = File::open(input_path).unwrap();
    let config = bincode::config::standard();
    bincode::decode_from_std_read::<Vec<complex::Complex<f32>>, _, _>(&mut file, config)
}

fn main() -> std::io::Result<()> {
    let freq = 128.0; // Frequency
    let phase_shift = 0.0; // Phase shift
    let num_samples = 128; // Number of samples
    
    // Write binary file
    let binary_path = Path::new("wave_data.bin");
    write_wave_binary(freq, phase_shift, num_samples, binary_path).unwrap();
    
    // Read binary file
    let read_data = read_wave_binary(binary_path).unwrap();

    // TODO : plot the data
    
    Ok(())
}
