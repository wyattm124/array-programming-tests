use crate::complex;
use complex::Complex;
use complex::Field;


fn convolve <T: Field> (a: &[Complex<T>], b: &[Complex<T>], result: &mut [Complex<T>]) -> () {
    for i in 0..a.len() {
        rev_madd(&a[i..i + b.len()], &b, &mut result[i]);
    }
}

fn correlate <T: Field> (a: &[Complex<T>], b: &[Complex<T>], result: &mut [Complex<T>]) -> () {
    for i in 0..a.len() {
        madd(&a[i..i + b.len()], &b, &mut result[i]);
    }
}

// Many Key components of efficient numerical algorithms come down to
//  a good compilation of these madd (multiply and add) loops
// NOTE : this is also a dot product
fn madd <T: Field> (a: &[Complex<T>], b: &[Complex<T>], c: &mut Complex<T>) -> () {
    for i in 0..a.len() {
        *c += a[i] * b[i];
    }
}

fn rev_madd <T: Field> (a: &[Complex<T>], b: &[Complex<T>], c: &mut Complex<T>) -> () {
    let len = a.len();
    for i in 0..len {
        *c += a[i] * b[(len - 1) - i];
    }
}