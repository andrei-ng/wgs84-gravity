use nalgebra as na;
use uom::si::f64::{Acceleration, Angle, AngularVelocity, Length, Velocity};
use uom::si::{
    acceleration::meter_per_second_squared, angle::radian, angular_velocity::radian_per_second,
    length::meter, velocity::meter_per_second,
};

pub trait NedComponents<T> {
    fn north(&self) -> T;
    fn east(&self) -> T;
    fn down(&self) -> T;
}

#[derive(Debug, Clone, Copy)]
pub struct Lla {
    lat: Angle,
    lon: Angle,
    alt: Length,
}

impl Lla {
    /// Create a new Lla from Angle and Length values
    pub fn new(lat: Angle, lon: Angle, alt: Length) -> Self {
        Self { lat, lon, alt }
    }

    #[cfg(test)]
    /// Create a new Lla with the given raw values as standard SI units for geodetic position
    pub(crate) fn new_from_raw_si(lat: f64, lon: f64, alt: f64) -> Self {
        Self {
            lat: Angle::new::<radian>(lat),
            lon: Angle::new::<radian>(lon),
            alt: Length::new::<meter>(alt),
        }
    }

    pub fn zeros() -> Self {
        Self {
            lat: Angle::new::<radian>(0.0),
            lon: Angle::new::<radian>(0.0),
            alt: Length::new::<meter>(0.0),
        }
    }

    pub fn lat(&self) -> Angle {
        self.lat
    }

    pub fn lon(&self) -> Angle {
        self.lon
    }

    pub fn alt(&self) -> Length {
        self.alt
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LlaRate {
    lat: AngularVelocity,
    lon: AngularVelocity,
    alt: Velocity,
}

impl LlaRate {
    /// Create a new LlaRate from AngularVelocity and Velocity values
    pub fn new(lat: AngularVelocity, lon: AngularVelocity, alt: Velocity) -> Self {
        Self { lat, lon, alt }
    }

    pub fn zeros() -> Self {
        Self {
            lat: AngularVelocity::new::<radian_per_second>(0.0),
            lon: AngularVelocity::new::<radian_per_second>(0.0),
            alt: Velocity::new::<meter_per_second>(0.0),
        }
    }

    /// Create a new LlaRate with the given raw values as standard SI units for rates
    pub(crate) fn new_from_raw_si(lat: f64, lon: f64, alt: f64) -> Self {
        Self {
            lat: AngularVelocity::new::<radian_per_second>(lat),
            lon: AngularVelocity::new::<radian_per_second>(lon),
            alt: Velocity::new::<meter_per_second>(alt),
        }
    }

    pub fn lat_rate(&self) -> AngularVelocity {
        self.lat
    }

    pub fn lon_rate(&self) -> AngularVelocity {
        self.lon
    }

    pub fn alt_rate(&self) -> Velocity {
        self.alt
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Velocity3d {
    inner: na::Vector3<Velocity>,
}

impl Velocity3d {
    /// Create a new Velocity3d from Velocity values
    pub fn new(x: Velocity, y: Velocity, z: Velocity) -> Self {
        Self {
            inner: na::Vector3::new(x, y, z),
        }
    }

    pub fn zeros() -> Self {
        Self {
            inner: na::Vector3::zeros(),
        }
    }

    #[cfg(test)]
    /// Create a new Velocity3d with the given raw values as standard SI units for velocities
    pub(crate) fn new_from_raw_si(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: na::Vector3::new(
                Velocity::new::<meter_per_second>(x),
                Velocity::new::<meter_per_second>(y),
                Velocity::new::<meter_per_second>(z),
            ),
        }
    }

    /// Create a new Velocity3d from a nalgebra vector of f64 values of velocity in SI units of m/s
    pub(crate) fn from_vector_raw_si(vector: na::Vector3<f64>) -> Self {
        Self {
            inner: na::Vector3::new(
                Velocity::new::<meter_per_second>(vector.x),
                Velocity::new::<meter_per_second>(vector.y),
                Velocity::new::<meter_per_second>(vector.z),
            ),
        }
    }

    /// Convert the Velocity3d to a nalgebra vector of f64 values with standard SI units of m/s
    pub(crate) fn to_vector_raw_si(self) -> na::Vector3<f64> {
        na::Vector3::new(
            self.inner.x.get::<meter_per_second>(),
            self.inner.y.get::<meter_per_second>(),
            self.inner.z.get::<meter_per_second>(),
        )
    }

    pub fn x(&self) -> Velocity {
        self.inner.x
    }

    pub fn y(&self) -> Velocity {
        self.inner.y
    }

    pub fn z(&self) -> Velocity {
        self.inner.z
    }

    /// Get the underlying nalgebra vector
    pub fn inner(&self) -> &na::Vector3<Velocity> {
        &self.inner
    }

    /// Get a mutable reference to the underlying nalgebra vector
    pub fn inner_mut(&mut self) -> &mut na::Vector3<Velocity> {
        &mut self.inner
    }

    /// Consumes the wrapper and returns the underlying nalgebra vector
    pub fn into_inner(self) -> na::Vector3<Velocity> {
        self.inner
    }
}

impl NedComponents<Velocity> for Velocity3d {
    fn north(&self) -> Velocity {
        self.inner.x
    }

    fn east(&self) -> Velocity {
        self.inner.y
    }

    fn down(&self) -> Velocity {
        self.inner.z
    }
}

impl std::ops::Add for Velocity3d {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Velocity3d::from_vector_raw_si(self.to_vector_raw_si() + other.to_vector_raw_si())
    }
}

impl std::ops::Sub for Velocity3d {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Velocity3d::from_vector_raw_si(self.to_vector_raw_si() - other.to_vector_raw_si())
    }
}

impl std::ops::Mul<f64> for Velocity3d {
    type Output = Self;

    fn mul(self, factor: f64) -> Self {
        Velocity3d::from_vector_raw_si(self.to_vector_raw_si() * factor)
    }
}

impl std::ops::Mul<Velocity3d> for f64 {
    type Output = Velocity3d;

    fn mul(self, vec: Velocity3d) -> Velocity3d {
        vec * self
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Acceleration3d {
    inner: na::Vector3<Acceleration>,
}

impl Acceleration3d {
    pub fn new(x: Acceleration, y: Acceleration, z: Acceleration) -> Self {
        Self {
            inner: na::Vector3::new(x, y, z),
        }
    }

    pub fn zeros() -> Self {
        Self {
            inner: na::Vector3::zeros(),
        }
    }

    #[cfg(test)]
    pub(crate) fn new_from_raw_si(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: na::Vector3::new(
                Acceleration::new::<meter_per_second_squared>(x),
                Acceleration::new::<meter_per_second_squared>(y),
                Acceleration::new::<meter_per_second_squared>(z),
            ),
        }
    }

    pub(crate) fn from_vector_raw_si(vector: na::Vector3<f64>) -> Self {
        Self {
            inner: na::Vector3::new(
                Acceleration::new::<meter_per_second_squared>(vector.x),
                Acceleration::new::<meter_per_second_squared>(vector.y),
                Acceleration::new::<meter_per_second_squared>(vector.z),
            ),
        }
    }

    pub(crate) fn to_vector_raw_si(self) -> na::Vector3<f64> {
        na::Vector3::new(
            self.inner.x.get::<meter_per_second_squared>(),
            self.inner.y.get::<meter_per_second_squared>(),
            self.inner.z.get::<meter_per_second_squared>(),
        )
    }

    pub fn x(&self) -> Acceleration {
        self.inner.x
    }

    pub fn y(&self) -> Acceleration {
        self.inner.y
    }

    pub fn z(&self) -> Acceleration {
        self.inner.z
    }

    /// Get the underlying nalgebra vector
    pub fn inner(&self) -> &na::Vector3<Acceleration> {
        &self.inner
    }

    /// Get a mutable reference to the underlying nalgebra vector
    pub fn inner_mut(&mut self) -> &mut na::Vector3<Acceleration> {
        &mut self.inner
    }

    /// Consumes the wrapper and returns the underlying nalgebra vector
    pub fn into_inner(self) -> na::Vector3<Acceleration> {
        self.inner
    }
}

impl NedComponents<Acceleration> for Acceleration3d {
    fn north(&self) -> Acceleration {
        self.inner.x
    }

    fn east(&self) -> Acceleration {
        self.inner.y
    }

    fn down(&self) -> Acceleration {
        self.inner.z
    }
}

impl std::ops::Add for Acceleration3d {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Acceleration3d::from_vector_raw_si(self.to_vector_raw_si() + other.to_vector_raw_si())
    }
}

impl std::ops::Sub for Acceleration3d {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Acceleration3d::from_vector_raw_si(self.to_vector_raw_si() - other.to_vector_raw_si())
    }
}

impl std::ops::Mul<f64> for Acceleration3d {
    type Output = Self;

    fn mul(self, factor: f64) -> Self {
        Acceleration3d::from_vector_raw_si(self.to_vector_raw_si() * factor)
    }
}

impl std::ops::Mul<Acceleration3d> for f64 {
    type Output = Acceleration3d;

    fn mul(self, vec: Acceleration3d) -> Acceleration3d {
        vec * self
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AngularVelocity3d {
    inner: na::Vector3<AngularVelocity>,
}

impl AngularVelocity3d {
    pub fn new(x: AngularVelocity, y: AngularVelocity, z: AngularVelocity) -> Self {
        Self {
            inner: na::Vector3::new(x, y, z),
        }
    }

    pub fn zeros() -> Self {
        Self {
            inner: na::Vector3::zeros(),
        }
    }

    pub(crate) fn new_from_raw_si(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: na::Vector3::new(
                AngularVelocity::new::<radian_per_second>(x),
                AngularVelocity::new::<radian_per_second>(y),
                AngularVelocity::new::<radian_per_second>(z),
            ),
        }
    }

    pub(crate) fn to_vector_raw_si(self) -> na::Vector3<f64> {
        na::Vector3::new(
            self.inner.x.get::<radian_per_second>(),
            self.inner.y.get::<radian_per_second>(),
            self.inner.z.get::<radian_per_second>(),
        )
    }

    pub(crate) fn from_vector_raw_si(vector: na::Vector3<f64>) -> Self {
        Self {
            inner: na::Vector3::new(
                AngularVelocity::new::<radian_per_second>(vector.x),
                AngularVelocity::new::<radian_per_second>(vector.y),
                AngularVelocity::new::<radian_per_second>(vector.z),
            ),
        }
    }

    pub fn x(&self) -> AngularVelocity {
        self.inner.x
    }

    pub fn y(&self) -> AngularVelocity {
        self.inner.y
    }

    pub fn z(&self) -> AngularVelocity {
        self.inner.z
    }

    /// Get the underlying nalgebra vector
    pub fn inner(&self) -> &na::Vector3<AngularVelocity> {
        &self.inner
    }

    /// Get a mutable reference to the underlying nalgebra vector
    pub fn inner_mut(&mut self) -> &mut na::Vector3<AngularVelocity> {
        &mut self.inner
    }

    /// Consumes the wrapper and returns the underlying nalgebra vector
    pub fn into_inner(self) -> na::Vector3<AngularVelocity> {
        self.inner
    }
}

impl NedComponents<AngularVelocity> for AngularVelocity3d {
    fn north(&self) -> AngularVelocity {
        self.inner.x
    }

    fn east(&self) -> AngularVelocity {
        self.inner.y
    }

    fn down(&self) -> AngularVelocity {
        self.inner.z
    }
}

// Mathematical operations for AngularVelocity3d
impl std::ops::Add for AngularVelocity3d {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        AngularVelocity3d::from_vector_raw_si(self.to_vector_raw_si() + other.to_vector_raw_si())
    }
}

impl std::ops::Sub for AngularVelocity3d {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        AngularVelocity3d::from_vector_raw_si(self.to_vector_raw_si() - other.to_vector_raw_si())
    }
}

impl std::ops::Mul<f64> for AngularVelocity3d {
    type Output = Self;

    fn mul(self, factor: f64) -> Self {
        AngularVelocity3d::from_vector_raw_si(self.to_vector_raw_si() * factor)
    }
}

impl std::ops::Mul<AngularVelocity3d> for f64 {
    type Output = AngularVelocity3d;

    fn mul(self, vec: AngularVelocity3d) -> AngularVelocity3d {
        vec * self
    }
}
