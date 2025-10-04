use crate::physical_types::{Acceleration3d, AngularVelocity3d, Lla, NedComponents, Velocity3d};
use nav_math::linear_algebra;
use uom::si::f64::Angle;
use uom::si::{angle::radian, length::meter, velocity::meter_per_second};

use tiny_wgs84::WGS84;

/// Calculate Earth rotation rate expressed in n-frame (NED) at provided latitude
/// See [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false) page 390-391
pub fn earth_rotation_rate_ned(lat: Angle) -> AngularVelocity3d {
    let lat_rad = lat.get::<radian>();

    // Implements w^n_ie term in formula 11.43 on page 390
    AngularVelocity3d::new_from_raw_si(
        WGS84::W_IE * lat_rad.cos(),
        0.0,
        -WGS84::W_IE * lat_rad.sin(),
    )
}

/// Calculate angular rate of n-frame (NED) relative to e-frame (ECEF) expressed in n-frame
/// This is the transport rate (navigation frame rotation relative to Earth-fixed frame e-frame (ECEF)
/// due to motion over the surface of the Earth)
pub fn angular_velocity_relative_to_e_frame(pos: &Lla, vel: &Velocity3d) -> AngularVelocity3d {
    // Reuse the equations of formula (11.50) page 391 to obtain the same result
    // without relying on pos_dot calculation
    // Also more efficient, since we don't multiply and divide by the same cos(lat) quantity

    let rm = WGS84::radius_meridian(pos.lat()).get::<meter>();
    let rn = WGS84::radius_normal(pos.lat()).get::<meter>();
    let alt = pos.alt().get::<meter>();
    let lat = pos.lat().get::<radian>();
    let v_east = vel.east().get::<meter_per_second>();
    let v_north = vel.north().get::<meter_per_second>();
    let common = v_east / (rn + alt);
    AngularVelocity3d::new_from_raw_si(common, -v_north / (rm + alt), -common * lat.tan())
}

/// Calculate the Coriolis effect in n-frame (NED frame)
/// See [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false) page 391
pub fn coriolis_acceleration(pos: &Lla, vel: &Velocity3d) -> Acceleration3d {
    // Angular velocity vector update formula (11.52) page 391
    let w_n_ie = earth_rotation_rate_ned(pos.lat());
    let w_n_en = angular_velocity_relative_to_e_frame(pos, vel);

    // Total angular velocity: ω_in = ω_en + 2 * ω_ie
    let omega = w_n_en + 2.0 * w_n_ie;
    let omega = omega.to_vector_raw_si();

    let w_skew = linear_algebra::skew_symmetric(&omega);

    Acceleration3d::from_vector_raw_si(w_skew * vel.to_vector_raw_si())
}
