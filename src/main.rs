mod complex;
mod accumulate;
mod view;

fn main() {
    let c1 = complex::Complex::<f32>::new(1.0, 2.0);
    let c2 = complex::Complex::<f32>::new(3.0, 4.0);
    let c3 = c1 + c2;
    println!("{} + {}i", c3.re, c3.im);
}
