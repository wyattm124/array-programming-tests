use std::ops;

pub trait Sqrt {
    fn sqrt(self) -> Self;
}

impl Sqrt for f32 {
    fn sqrt(self) -> f32 {
        self.sqrt()
    }
}

pub trait Num : Copy +
                ops::Neg<Output = Self> +
                ops::Add<Output = Self> +
                ops::Sub<Output = Self> +
                ops::Mul<Output = Self> +
                ops::Div<Output = Self> +
                ops::AddAssign +
                ops::SubAssign +
                ops::MulAssign +
                ops::DivAssign {
    // Just need to make sure a Complex is made with an actual number
}

pub trait Field : Num + Sqrt {}

impl Field for f32 {}
impl Num for f32 {}

impl Num for Complex<f32> {}

#[derive(Debug, Clone, Copy)]
#[derive(PartialEq, Eq)]
pub struct Complex<T: Num> {
    pub re: T,
    pub im: T,
}

// Unary Negative
impl<T: Field> ops::Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}

// Binary, element wise addition and subtraction
impl<T: Field> ops::Add<Complex<T>> for Complex<T> {
    type Output = Complex<T>;
    fn add(self, other: Complex<T>) -> Self::Output {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl<T: Field> ops::Sub<Complex<T>> for Complex<T> {
    type Output = Complex<T>;
    fn sub(self, other: Complex<T>) -> Self::Output {
        Complex {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

// Scalar multiplication and division
impl<T: Field> ops::Mul<T> for Complex<T> {
    type Output = Complex<T>;
    fn mul(self, other: T) -> Self::Output {
        Complex {
            re: self.re * other,
            im: self.im * other,
        }
    }
}

impl<T: Field> ops::Div<T> for Complex<T> {
    type Output = Complex<T>;
    fn div(self, other: T) -> Self::Output {
        Complex {
            re: self.re / other,
            im: self.im / other,
        }
    }
}

// Binary Complex multiplaction and Division
impl<T: Field> ops::Mul<Complex<T>> for Complex<T> {
    type Output = Complex<T>;
    fn mul(self, other: Complex<T>) -> Self::Output {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

impl<T: Field> ops::Div<Complex<T>> for Complex<T> {
    type Output = Complex<T>;
    fn div(self, other: Complex<T>) -> Self::Output {
        self * other.inverse()
    }
}

// Assignment Operations
impl<T: Field> ops::AddAssign<Complex<T>> for Complex<T> {
    fn add_assign(&mut self, other: Complex<T>) {
        *self = *self + other;
    }
}

impl<T: Field> ops::SubAssign<Complex<T>> for Complex<T> {
    fn sub_assign(&mut self, other: Complex<T>) {
        *self = *self - other;
    }
}

impl<T: Field> ops::MulAssign<Complex<T>> for Complex<T> {
    fn mul_assign(&mut self, other: Complex<T>) {
        *self = *self * other;
    }
}

impl<T: Field> ops::DivAssign<Complex<T>> for Complex<T> {
    fn div_assign(&mut self, other: Complex<T>) {
        *self = *self / other;
    }
}

// Method specifcs for Complex<T>
impl<T: Field> Complex<T> {
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }

    pub fn conjugate(&self) -> Self {
        Self { re: self.re, im: -self.im }
    }

    pub fn magnitude_squared(&self) -> T {
        self.re * self.re + self.im * self.im
    }

    pub fn magnitude(&self) -> T {
        self.magnitude_squared().sqrt()
    }

    pub fn inverse(&self) -> Complex<T> {
        self.conjugate() / self.magnitude_squared()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_scalar_mult() {
        let a = Complex::<f32>::new(1.0, 2.0);
        assert_eq!(a * 2.0, Complex::<f32>::new(2.0, 4.0));
    }

    #[test]
    fn complex_scalar_div() {
        let a = Complex::<f32>::new(1.0, 2.0);
        assert_eq!(a / 2.0, Complex::<f32>::new(0.5, 1.0));
    }

    #[test]
    fn complex_negate() {
        let a = Complex::<f32>::new(1.0, 2.0);
        assert_eq!(-a, Complex::<f32>::new(-1.0, -2.0));
    }

    #[test]
    fn complex_conjugate() {
        let a = Complex::<f32>::new(1.0, 2.0);
        assert_eq!(a.conjugate(), Complex::<f32>::new(1.0, -2.0));
    }

    #[test]
    fn complex_inverse() {
        let a = Complex::<f32>::new(1.0, 2.0);
        assert_eq!(a.inverse(), Complex::<f32>::new(0.2, -0.4));
    }

    #[test]
    fn complex_magnitude() {
        let a = Complex::<f32>::new(3.0, 4.0);
        let b = Complex::<f32>::new(-4.0, 3.0);
        assert_eq!(a.magnitude_squared(), 25.0);
        assert_eq!(b.magnitude_squared(), 25.0);
        assert_eq!(a.magnitude(), 5.0);
        assert_eq!(b.magnitude(), 5.0);
    }

    #[test]
    fn complex_add() {
        let a = Complex::<f32>::new(1.0, 2.0);
        let b = Complex::<f32>::new(3.0, 4.0);
        assert_eq!(a + b, Complex::<f32>::new(4.0, 6.0));
    }

    #[test]
    fn complex_sub() {
        let a = Complex::<f32>::new(1.0, 2.0);
        let b = Complex::<f32>::new(3.0, 4.0);
        assert_eq!(a - b, Complex::<f32>::new(-2.0, -2.0));
    }

    #[test]
    fn complex_mult() {
        let a = Complex::<f32>::new(1.0, 2.0);
        let b = Complex::<f32>::new(3.0, 4.0);
        assert_eq!(a * b, Complex::<f32>::new(-5.0, 10.0));
    }

    #[test]
    fn complex_div() {
        let a = Complex::<f32>::new(20.0, -4.0);
        let b = Complex::<f32>::new(3.0, 2.0);
        assert_eq!(a / b, Complex::<f32>::new(4.0, -4.0));
    }
}

