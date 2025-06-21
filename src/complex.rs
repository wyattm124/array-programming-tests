use std::ops;

#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub re: f32,
    pub im: f32,
}

impl ops::Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl ops::Add for Complex {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

impl ops::Sub for Complex {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

impl ops::Mul<f32> for Complex {
    type Output = Self;
    fn mul(self, other: f32) -> Self {
        Self {
            re: self.re * other,
            im: self.im * other,
        }
    }
}

impl ops::Div<f32> for Complex {
    type Output = Self;
    fn div(self, other: f32) -> Self {
        Self {
            re: self.re / other,
            im: self.im / other,
        }
    }
}

impl ops::Mul for Complex {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

impl ops::Div for Complex {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        self * other.inverse()
    }
}

impl ops::AddAssign for Complex {
    fn add_assign(&mut self, other: Self) {
        self.re += other.re;
        self.im += other.im;
    }
}

impl Complex {
    pub fn new(re: f32, im: f32) -> Self {
        Self { re, im }
    }

    pub fn conjugate(&self) -> Self {
        Self { re: self.re, im: -self.im }
    }

    pub fn magnitude_squared(&self) -> f32 {
        self.re * self.re + self.im * self.im
    }

    pub fn magnitude(&self) -> f32 {
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

