use std::ops;

pub trait Sqrt {
    fn sqrt(self) -> Self;
}

impl Sqrt for f32 {
    fn sqrt(self) -> f32 {
        self.sqrt()
    }
}

impl Sqrt for f64 {
    fn sqrt(self) -> f64 {
        self.sqrt()
    }
}

pub trait Num : Copy + 
                Sqrt +
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

impl Num for f32 {}
impl Num for f64 {}

impl Num for Complex<f32> {}
impl Num for Complex<f64> {}

#[derive(Debug, Clone, Copy)]
pub struct Complex<T: Num> {
    pub re: T,
    pub im: T,
}

impl<T: Num> ops::Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl<'a, 'b, T: Num> ops::Add<&'b Complex<T>> for &'a Complex<T> {
    type Output = Complex<T>;
    fn add(self, other: &'b Complex<T>) -> Self::Output {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl<'a, 'b, T: Num> ops::Add<&'b Complex<T>> for &'a Complex<T> {
    type Output = Complex<T>;
    fn add(self, other: &'b Complex<T>) -> Self::Output {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl<'a, 'b, T: Num> ops::Sub<&'b Complex<T>> for &'a Complex<T> {
    type Output = Complex<T>;
    fn sub(self, other: &'b Complex<T>) -> Self::Output {
        Complex {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

impl<'a, 'b, T: Num> ops::Mul<&'b T> for &'a Complex<T> {
    type Output = Complex<T>;
    fn mul(self, other: &'b T) -> Self::Output {
        Complex {
            re: self.re * *other,
            im: self.im * *other,
        }
    }
}

impl<'a, 'b, T: Num> ops::Div<&'b T> for &'a Complex<T> {
    type Output = Complex<T>;
    fn div(self, other: &'b T) -> Self::Output {
        Complex {
            re: self.re / *other,
            im: self.im / *other,
        }
    }
}

impl<'a, 'b, T: Num> ops::Mul<&'b Complex<T>> for &'a Complex<T> {
    type Output = Complex<T>;
    fn mul(self, other: &'b Complex<T>) -> Self::Output {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

impl<'a, 'b, T: Num> ops::Div<&'b Complex<T>> for &'a Complex<T> {
    type Output = Complex<T>;
    fn div(self, other: &'b Complex<T>) -> Self::Output {
        let temp = other.inverse();
        self * &temp
    }
}

impl<'a, 'b, T: Num> ops::AddAssign<&'b Complex<T>> for &'a mut Complex<T> {
    fn add_assign(&mut self, other: &'b Complex<T>) {
        *self = *self + other;
    }
}

impl<'a, 'b, T: Num> ops::SubAssign<&'b Complex<T>> for &'a mut Complex<T> {
    fn sub_assign(&mut self, other: &'b Complex<T>) {
        *self = *self - other;
    }
}

impl<'a, 'b, T: Num> ops::MulAssign<&'b T> for &'a mut Complex<T> {
    fn mul_assign(&mut self, other: &'b T) {
        *self = *self * other;
    }
}

impl<'a, 'b, T: Num> ops::DivAssign<&'b T> for &'a mut Complex<T> {
    fn div_assign(&mut self, other: &'b T) {
        *self = *self / other;
    }
}

impl<T: Num> Complex<T> {
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

    pub fn inverse(&self) -> Self {
        let denominator = self.magnitude_squared();
        Self {
            re: self.re / denominator,
            im: -self.im / denominator,
        }
    }
}

