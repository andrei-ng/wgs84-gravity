//! This module implements the mechanization equations for Strapdown Inertial Navigation Systems (INS)
//! in the North-East-Down (NED) frame.
//!
//! # References
//!
//! Implementation based on equations in Section 11.2.4 of [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onage&q&f=false)
//!
use uom::si::{angle::radian, length::meter, velocity::meter_per_second};

use crate::gravity;
use crate::inertial;
use crate::physical_types::{Acceleration3d, Lla, LlaRate, NedComponents, Velocity3d};
use tiny_wgs84::WGS84;

/// Navigation frame mechanization (in NED frame)
/// Implementation is based on the equations as described in Section 11.2.4, page 390 - 392 of
/// [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false) page 390-391
#[derive(Debug, Clone, Copy)]
pub struct Mechanization {
    /// Enable/disable rotating earth assumption.
    /// When the flag is enabled, earth's rotation relative to the inertial frame
    /// is ignored in the navigation mechanization equations
    static_earth: bool,
    /// Select gravity model to use in calculation of gravitational acceleration,
    /// see `gravity::Mode` for possible options.
    /// By default the NED model for gravity is used where g = [North, 0.0, Down]
    gravity_mode: gravity::Mode,
}

/// Create a default Mechanization object with
///  - rotating earth modelling
///  - NED decomposition of the normal gravitational acceleration
impl Default for Mechanization {
    fn default() -> Self {
        Self {
            static_earth: false,
            gravity_mode: gravity::Mode::Ned,
        }
    }
}

impl Mechanization {
    /// Enable/disable the rotating earth assumption.
    /// Take into account earth's rotation w.r.t. the inertial frame in navigation mechanization calculations
    pub fn static_earth(&mut self, yes: bool) {
        self.static_earth = yes;
    }

    /// Select gravity model to use in mechanization calculations.
    /// See `gravity::Mode` for possible options.
    pub fn gravity_model(&mut self, mode: gravity::Mode) {
        self.gravity_mode = mode;
    }

    /// Get local value of gravitational acceleration as used in mechanization calculations
    pub fn gravitational_acceleration(&self, pos: &Lla) -> Acceleration3d {
        gravity::gravity(pos.lat(), pos.alt(), self.gravity_mode)
    }

    /// Calculates the rate of change of geodetic position from current position and current velocity in NED frame.
    ///
    /// This implements the position update equation for n-frame (NED)  mechanization,
    /// computing the derivatives of latitude, longitude, and altitude in closed form.
    ///
    /// Takes the current geodetic position [`Lla`] and velocity [`Velocity3d`] in NED frame
    /// and computes the rate of change of position in NED frame [`LlaRate`].
    ///
    /// # References
    ///
    /// Formula (11.50) page 391 of Farrell's "AIDED Navigation - GPS With High Rate Sensors"
    pub fn position_rate(&self, pos: &Lla, vel_ned: &Velocity3d) -> LlaRate {
        // Implements formula (11.50) page 391
        // [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false)

        let rm = WGS84::radius_meridian(pos.lat()).get::<meter>();
        let rn = WGS84::radius_normal(pos.lat()).get::<meter>();

        let v_north = vel_ned.north().get::<meter_per_second>();
        let v_east = vel_ned.east().get::<meter_per_second>();
        let v_down = vel_ned.down().get::<meter_per_second>();
        let lat = pos.lat().get::<radian>();
        let alt = pos.alt().get::<meter>();

        let lat_dot = v_north / (rm + alt);
        let lon_dot = v_east / lat.cos() / (rn + alt);
        let alt_dot = -v_down;

        LlaRate::new_from_raw_si(lat_dot, lon_dot, alt_dot)
    }

    /// Calculates the rate of change of NED velocity from current position, current velocity and specific force in NED frame.
    ///
    /// This implements the velocity update equation for n-frame (NED)  mechanization,
    /// computing the derivatives of velocity in the North, East, and Down directions in closed form.
    ///
    /// Takes the current specific force in the navigation frame (NED), current position in [`Lla`],
    /// current velocity in the navigation frame (NED)[`Velocity3d`], and optionally the local gravitational
    /// acceleration in vector form [`Acceleration3d`] and returns the rate of change of velocity in the navigation
    /// frame (NED)[`Acceleration3d`].
    ///
    /// When the gravitational acceleration is not provided explicitly, it is calculated on the fly using the input  position [`Lla`]
    /// and the build-in gravity model selected at initialization.
    ///
    /// # References
    ///
    /// Formula (11.51) page 391 of Farrell's "AIDED Navigation - GPS With High Rate Sensors"
    ///
    pub fn velocity_rate(
        &mut self,
        specific_force: &Acceleration3d,
        pos: &Lla,
        vel_ned: &Velocity3d,
        g: Option<Acceleration3d>,
    ) -> Acceleration3d {
        let g = match g {
            Some(g) => g,
            None => gravity::gravity(pos.lat(), pos.alt(), self.gravity_mode),
        };

        if self.static_earth {
            let coriolis_accel = inertial::coriolis_acceleration(pos, vel_ned);
            *specific_force + g - coriolis_accel
        } else {
            *specific_force + g
        }
    }
}
