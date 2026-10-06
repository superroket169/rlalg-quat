use crate::{Sqrt, v3f};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quat {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

impl Quat {
    pub const IDENTITY: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
        w: 1.0,
    };

    pub const fn new(x: f32, y: f32, z: f32, w: f32) -> Self {
        Self { x, y, z, w }
    }

    pub fn from_axis_angle(axis: v3f, angle: f32) -> Self {
        let half = angle * 0.5;
        let (s, c) = (libm::sinf(half), libm::cosf(half));

        Self {
            x: axis.x * s,
            y: axis.y * s,
            z: axis.z * s,
            w: c,
        }
    }

    pub fn mag_sq(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w
    }

    pub fn mag(&self) -> f32 {
        self.mag_sq().sqrt()
    }

    pub fn normalize(self) -> Self {
        let m = self.mag();

        Self {
            x: self.x / m,
            y: self.y / m,
            z: self.z / m,
            w: self.w / m,
        }
    }

    pub fn conjugate(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }

    pub fn rotate(self, v: v3f) -> v3f {
        let qv = Self::new(v.x, v.y, v.z, 0.0);
        let r = self * qv * self.conjugate();
        v3f::new(r.x, r.y, r.z)
    }

    pub fn integrate(self, omega: v3f, dt: f32) -> Self {
        let omega_q = Self::new(omega.x, omega.y, omega.z, 0.0);
        let q_dot = self * omega_q;

        Self {
            x: self.x + q_dot.x * 0.5 * dt,
            y: self.y + q_dot.y * 0.5 * dt,
            z: self.z + q_dot.z * 0.5 * dt,
            w: self.w + q_dot.w * 0.5 * dt,
        }
        .normalize()
    }
}

impl Default for Quat {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl core::ops::Mul for Quat {
    type Output = Self;

    fn mul(self, o: Self) -> Self {
        Self {
            w: self.w * o.w - self.x * o.x - self.y * o.y - self.z * o.z,
            x: self.w * o.x + self.x * o.w + self.y * o.z - self.z * o.y,
            y: self.w * o.y - self.x * o.z + self.y * o.w + self.z * o.x,
            z: self.w * o.z + self.x * o.y - self.y * o.x + self.z * o.w,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn identity_rotate_is_noop() {
        let v = v3f::new(1.0, 2.0, 3.0);
        assert_eq!(Quat::IDENTITY.rotate(v), v);
    }

    #[test]
    fn quarter_turn_about_z() {
        let q = Quat::from_axis_angle(v3f::new(0.0, 0.0, 1.0), core::f32::consts::FRAC_PI_2);
        let r = q.rotate(v3f::new(1.0, 0.0, 0.0));

        assert!((r.x).abs() < 1e-5);
        assert!((r.y - 1.0).abs() < 1e-5);
        assert!((r.z).abs() < 1e-5);
    }

    #[test]
    fn normalize_gives_unit_mag() {
        let q = Quat::new(1.0, 2.0, 3.0, 4.0).normalize();
        assert!((q.mag() - 1.0).abs() < 1e-5);
    }
}
