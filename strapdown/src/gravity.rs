use tiny_wgs84::WGS84;

use uom::si::acceleration::meter_per_second_squared;
use uom::si::f64::{Acceleration, Angle, Length};

use crate::physical_types::Acceleration3d;

/// Gravity model selection
#[derive(Default, Debug, Clone, Copy)]
pub enum Mode {
    /// Zero gravity
    /// This is useful when testing or simulating gravity-free environments
    Null,
    /// Use constant gravitational acceleration for Down component, independent of altitude and latitude.
    /// The value is the average of the Equatorial and Polar gravitational acceleration values from
    /// the WGS84 model. For more details see: `tiny-wgs84::WGS84`
    Constant,
    /// Use only the WGS84 normal gravity as the down component and zero for North and East components.
    /// For more details see: `tiny-wgs84::WGS84`
    Normal,
    /// Use NED model for gravity components
    /// `g = [North, 0.0, Down]`
    /// where the North component approximation is based on
    /// [Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems (GNSS Technology and Applications) ](https://books.google.nl/books?id=t94fAgAAQBAJ&printsec=copyrigh) - formulas (2.139) and (2.140).
    /// The Down component is the Normal gravity as given by the WGS84 model
    #[default]
    Ned,
}

pub fn gravity(lat: Angle, alt: Length, mode: Mode) -> Acceleration3d {
    let zero_accel = Acceleration::new::<meter_per_second_squared>(0.0);

    match mode {
        Mode::Null => Acceleration3d::new(zero_accel, zero_accel, zero_accel),
        Mode::Constant => Acceleration3d::new(
            zero_accel,
            zero_accel,
            Acceleration::new::<meter_per_second_squared>((WGS84::G_EQ + WGS84::G_POLE) / 2.0),
        ),
        Mode::Normal => {
            Acceleration3d::new(zero_accel, zero_accel, WGS84::gravity_normal(lat, alt))
        }
        Mode::Ned => Acceleration3d::new(
            WGS84::gravity_north(lat, alt),
            zero_accel,
            WGS84::gravity_down(lat, alt),
        ),
    }
}

#[cfg(test)]
mod tests {
    use crate::physical_types::Acceleration3d;
    use uom::si::f64::{Angle, Length};
    use uom::si::{angle::degree, length::meter};

    use super::*;

    #[test]
    fn gravitational_acceleration() {
        let zero_accel = Acceleration::new::<meter_per_second_squared>(0.0);

        let lat = Angle::new::<degree>(35.0);
        let alt = Length::new::<meter>(10.0);

        assert_eq!(
            gravity(lat, alt, crate::gravity::Mode::Constant).inner(),
            Acceleration3d::new(
                zero_accel,
                zero_accel,
                Acceleration::new::<meter_per_second_squared>(9.806255136883562)
            )
            .inner()
        );

        assert_eq!(
            gravity(lat, alt, crate::gravity::Mode::Normal).inner(),
            Acceleration3d::new(
                zero_accel,
                zero_accel,
                Acceleration::new::<meter_per_second_squared>(9.797305150189478)
            )
            .inner()
        );

        assert_eq!(
            gravity(lat, alt, crate::gravity::Mode::Ned).inner(),
            Acceleration3d::new(
                Acceleration::new::<meter_per_second_squared>(-7.59271637595014e-8),
                zero_accel,
                Acceleration::new::<meter_per_second_squared>(9.797305150189478)
            )
            .inner()
        );
    }
}
