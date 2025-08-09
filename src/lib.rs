/// Gravity Model using the WGS84 reference ellipsoid
pub struct WGS84;

impl WGS84 {
    /// Reciprocal Earth Ellipsoid Flattening
    pub const F_RECIPROCAL: f64 = 298.257_223_563;
    /// Earth Ellipsoid Flattening: 1/F_RECIPROCAL = 0.003_352_810_664_747_480_719_845_528_618_52;
    pub const F: f64 = 1.0 / Self::F_RECIPROCAL;
    /// Earth Equatorial Radius (semi-major axis) \[m\]
    pub const A: f64 = 6_378_137.0;
    /// Earth Ellipsoid semi-minor axis \[m\]
    ///
    /// Computed from `A` and `F` as `B` = `A` * (1 - `F`), resulting in `B` = 6_356_752.314_245_179 \[m\]
    pub const B: f64 = Self::A * (1.0 - Self::F);

    /// First eccentricity
    pub const E1: f64 = 8.181_919_084_262_2e-2;

    /// First eccentricity squared, which can be computed in two ways (implementation uses the first)
    /// ```
    /// use wgs84_gravity::WGS84;
    /// let e1_squared = (2.0 - WGS84::F)*WGS84::F;
    /// // or alternative calculation
    /// let e1_squared = 1.0 - f64::powf(WGS84::B/WGS84::A, 2.0);
    /// ```
    pub const E1_SQ: f64 = 2.0 * Self::F - Self::F * Self::F;

    /// Second eccentricity
    pub const E2: f64 = 8.209_443_794_969_6e-2;

    /// Second eccentricity squared, which can be computed in two ways (implementation uses the first)
    /// ```
    /// use wgs84_gravity::WGS84;
    /// let e2_squared = WGS84::E1_SQ / (1.0 - WGS84::E1_SQ);
    /// // or alternative calculation
    /// let e2_squared = WGS84::F  * (2.0 - WGS84::F) / f64::powf(1.0 - WGS84::F, 2.0);
    /// ```
    pub const E2_SQ: f64 = WGS84::E1_SQ / (1.0 - WGS84::E1_SQ);

    /// Universal Gravitation \[N * m<sup>2</sup> / kg<sup>2</sup>\]
    pub const G: f64 = 6.674 * 1e-11;
    /// Newtonian Geocentric Gravitational Constant \[m<sup>3</sup>/s<sup>2</sup>\]
    pub const GM: f64 = 3.986_004_418 * 1e14;
    /// Earth Rotational Rate relative to the inertial frame \[rad/s\]
    pub const W_IE: f64 = 7.292_115 * 1e-5;
    /// Earth Dynamic Flattening Form Factor
    pub const J2: f64 = 1.082_63 * 1e-3;
    /// Equatorial effective gravity \[m/s<sup>2</sup>\]
    pub const G_EQ: f64 = 9.780_325_335_904_06;
    /// Normal Gravity at Poles \[m/s<sup>2</sup>\]
    pub const G_POLE: f64 = 9.832_184_937_863_065;
    /// Mass of Earth: approx 5.972 * 1e24 \[kg\]
    pub const M_EARTH: f64 = WGS84::W_IE * WGS84::W_IE * WGS84::A * WGS84::A * WGS84::B / WGS84::GM;
    /// North Gravity component scaling factor
    /// Reference:[Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems (GNSS Technology and Applications) ](https://books.google.nl/books?id=t94fAgAAQBAJ&printsec=copyrigh) - formula (2.140)
    const G_NORTH_SCALING: f64 = -8.08 * 1e-9;

    /// Meridian radius of the Earth R<sub>M</sub> \[m\]
    ///
    /// Reference: [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false) - formula (2.6)  
    pub fn radius_meridian(lat: f64) -> f64 {
        let sin_lat_sq = Self::sin_lat_sq(lat);
        WGS84::A * (1.0_f64 - WGS84::E1_SQ) / f64::powf(1.0_f64 - WGS84::E1_SQ * sin_lat_sq, 1.5)
    }

    /// Normal radius of the Earth R<sub>N</sub> \[m\]
    ///
    /// Reference: [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false) - formula (2.7)
    pub fn radius_normal(lat: f64) -> f64 {
        let sin_lat_sq = Self::sin_lat_sq(lat);
        WGS84::A / f64::powf(1.0 - WGS84::E1_SQ * sin_lat_sq, 0.5)
    }

    /// Estimation of Normal Gravity component (down component) [m/s<sup>2</sup>] based on International Gravity Formula including altitude/height dependence
    ///
    /// Based on the International Gravity Formula [Walter D. Lambert](https://earth.geology.yale.edu/~ajs/1945A/360.pdf)
    ///
    /// Truncated Taylor series of the original formulation is used as described [here](https://ahrs.readthedocs.io/en/latest/wgs84.html#ahrs.utils.wgs84.international_gravity)
    /// since for the case of the International Ellipsoid, the third-order terms are negligible.
    pub fn gravity_down(lat: f64, alt: f64) -> f64 {
        let gravity_normal_surface = Self::gravity_surface(lat);
        let k1 = 2.0 * (1.0 + WGS84::F + WGS84::M_EARTH) / WGS84::A;
        let k2 = 4.0 * WGS84::F / WGS84::A;
        let k3 = 3.0 / (WGS84::A * WGS84::A);

        let sin_lat_sq = Self::sin_lat_sq(lat);

        gravity_normal_surface * (1.0 - (k1 - k2 * sin_lat_sq) * alt + k3 * alt * alt)
    }

    /// Same as `gravity_down`
    ///
    /// This is a convenience function to match the naming convention of the AHRS Python package
    pub fn gravity_normal(lat: f64, alt: f64) -> f64 {
        Self::gravity_down(lat, alt)
    }

    /// Estimation of Normal Gravity component (down component) [m/s<sup>2</sup>] at the surface of the reference ellipsoid based on the Somigliana equation
    ///
    /// Reference: [Theoretical_gravity#Somigliana_equation](https://en.wikipedia.org/wiki/Theoretical_gravity#Somigliana_equation)
    pub fn gravity_surface(lat: f64) -> f64 {
        const K: f64 = (WGS84::B * WGS84::G_POLE) / (WGS84::A * WGS84::G_EQ) - 1.0;

        let sin_lat_sq = Self::sin_lat_sq(lat);
        WGS84::G_EQ * (1.0 + K * sin_lat_sq) / f64::sqrt(1.0 - WGS84::E1_SQ * sin_lat_sq)
    }

    /// Estimation of North Gravity component [m/s<sup>2</sup>]
    ///
    /// Reference: [Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems (GNSS Technology and Applications) ](https://books.google.nl/books?id=t94fAgAAQBAJ&printsec=copyrigh) - formula (2.139) and (2.140)
    pub fn gravity_north(lat: f64, alt: f64) -> f64 {
        WGS84::G_NORTH_SCALING * alt * f64::sin(2.0 * lat)
    }

    /// Gravity vector in NED representation [m/s<sup>2</sup>] at provided latitude
    pub fn gravity_ned(lat: f64, alt: f64) -> [f64; 3] {
        [
            WGS84::gravity_north(lat, alt),
            0.0,
            WGS84::gravity_down(lat, alt),
        ]
    }

    /// Earth rotation [rad/s] evaluated in navigation frame (NED frame)
    /// Reference: [AIDED Navigation - GPS With High Rate Sensors, Jay Farrell](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false)
    /// Implements W <sup>n</sup> <sub>ie</sub> term in formula 12.43 on page 390 - 391
    #[rustfmt::skip]
    pub fn earth_rotation_ned(lat: f64) -> [f64; 3] {
        [
            Self::W_IE * lat.cos(),
            0.0,
            -Self::W_IE * lat.sin()
        ]
    }

    /// Useful quantity sin(lat)^2 used throughout WGS84 formulas
    fn sin_lat_sq(lat: f64) -> f64 {
        let sin_lat = f64::sin(lat);
        sin_lat * sin_lat
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, relative_eq};

    const EPS_SMALL: f64 = 1e-21;
    const EPS_NORMAL: f64 = 1e-16;

    #[test]
    fn first_eccentricity_squared() {
        // There are two methods to calculate eccentricity of the reference ellipsoid
        // Either:  e^2 = (2 -f)*f or e^2 = 1 - b^2/a^2
        // We use this test to check that the constants are well defined and give same result
        const BA_RATIO: f64 = WGS84::B / WGS84::A;
        assert_relative_eq!(
            WGS84::E1_SQ,
            (1.0 - BA_RATIO * BA_RATIO),
            epsilon = EPS_NORMAL
        );
    }

    #[test]
    fn second_eccentricity_squared() {
        // There are two methods to calculate eccentricity of the reference ellipsoid
        // Either:  e^2 = (2 -f)*f / (1 - e^2) or e^2 = 1 - b^2/a^2
        // We use this test to check that the constants are well defined and give same result
        assert_relative_eq!(
            WGS84::E2_SQ,
            (2.0 - WGS84::F) * WGS84::F / f64::powf(1.0 - WGS84::F, 2.0),
            epsilon = EPS_NORMAL
        );
    }

    #[test]
    fn normal_gravity_at_surface() {
        struct TestCase {
            name: &'static str,
            lat: f64,
            expected: f64,
            tol: f64,
        }
        let cases = [
            TestCase {
                name: "North Pole",
                lat: 90.0_f64.to_radians(),
                expected: WGS84::G_POLE,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "Equator",
                lat: 0.0,
                expected: WGS84::G_EQ,
                tol: EPS_SMALL,
            },
            // Expected value based on the implementation of WGS84 from AHRS Python package
            // https://ahrs.readthedocs.io/en/latest/wgs84.html#
            TestCase {
                name: "50°N",
                lat: 50.0_f64.to_radians(),
                expected: 9.810702135603085,
                tol: EPS_SMALL,
            },
        ];

        for c in cases {
            let got = WGS84::gravity_surface(c.lat);
            assert!(
                relative_eq!(got, c.expected, epsilon = c.tol),
                "test case: {}, got: {}, expected: {}, epsilon: {}",
                c.name,
                got,
                c.expected,
                c.tol
            );
        }
    }

    #[test]
    fn normal_gravity_at_0m() {
        struct TestCase {
            name: &'static str,
            lat: f64,
            alt: f64,
            expected: f64,
            tol: f64,
        }
        let cases = [
            TestCase {
                name: "North Pole",
                lat: 90.0_f64.to_radians(),
                alt: 0.0,
                expected: WGS84::G_POLE,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "Equator",
                lat: 0.0,
                alt: 0.0,
                expected: WGS84::G_EQ,
                tol: EPS_SMALL,
            },
        ];

        for c in cases {
            let got = WGS84::gravity_down(c.lat, c.alt);
            assert!(
                relative_eq!(got, c.expected, epsilon = c.tol),
                "test case: {}, got: {}, expected: {}, epsilon: {}",
                c.name,
                got,
                c.expected,
                c.tol
            );
        }
    }

    #[test]
    fn normal_gravity_alias() {
        assert_eq!(
            WGS84::gravity_down(45.0_f64.to_radians(), 0.0),
            WGS84::gravity_normal(45.0_f64.to_radians(), 0.0),
        );
    }

    #[test]
    fn normal_gravity_above_surface() {
        // Expected values based on the implementation of WGS84 from AHRS Python package
        // https://ahrs.readthedocs.io/en/latest/wgs84.html#

        struct TestCase {
            name: &'static str,
            lat: f64,
            alt: f64,
            expected: f64,
            tol: f64,
        }

        let cases = [
            TestCase {
                name: "50°N, 1000m",
                lat: 50.0_f64.to_radians(),
                alt: 1000.0,
                expected: 9.807617683884756,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "44.43771100°N, 80.41m",
                lat: 44.43771100_f64.to_radians(),
                alt: 80.41,
                expected: 9.805440780030928,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "52.37021570°N, 3.46m",
                lat: 52.37021570_f64.to_radians(),
                alt: 3.46,
                expected: 9.8127883621628,
                tol: EPS_SMALL,
            },
        ];

        for c in cases {
            let got = WGS84::gravity_down(c.lat, c.alt);
            assert!(
                relative_eq!(got, c.expected, epsilon = c.tol),
                "test case: {}, got: {}, expected: {}, epsilon: {}",
                c.name,
                got,
                c.expected,
                c.tol
            );
        }
    }

    #[test]
    fn gravity_north() {
        struct TestCase {
            name: &'static str,
            lat: f64,
            alt: f64,
            expected: f64,
            tol: f64,
        }
        let cases = [
            TestCase {
                name: "Equator, 0m",
                lat: 0.0,
                alt: 0.0,
                expected: 0.0,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "Equator, 1m",
                lat: 0.0,
                alt: 1.0,
                expected: 0.0,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "North Pole, 10m",
                lat: 90.0_f64.to_radians(),
                alt: 10.0,
                expected: 0.0,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "45°N, 80m",
                lat: 45.0_f64.to_radians(),
                alt: 80.0,
                expected: 80.0 * WGS84::G_NORTH_SCALING,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "45°N, 1/G_NORTH_SCALING",
                lat: 45.0_f64.to_radians(),
                alt: 1.0 / WGS84::G_NORTH_SCALING,
                expected: 1.0,
                tol: EPS_SMALL,
            },
        ];

        for c in cases {
            let got = WGS84::gravity_north(c.lat, c.alt);
            assert!(
                relative_eq!(got, c.expected, epsilon = c.tol),
                "test case: {}, got: {}, expected: {}, epsilon: {}",
                c.name,
                got,
                c.expected,
                c.tol
            );
        }
    }

    #[test]
    fn meridian_curvature_radius() {
        // Expected values based on the implementation of WGS84 from AHRS Python package
        // https://ahrs.readthedocs.io/en/latest/wgs84.html#

        struct TestCase {
            name: &'static str,
            lat: f64,
            expected: f64,
            tol: f64,
        }
        let cases = [
            TestCase {
                name: "Equator",
                lat: 0.0,
                expected: 6335439.327292819507420063018798828125,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "50°N",
                lat: 50.0_f64.to_radians(),
                expected: 6372955.92573519889265298843383789062500,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "North Pole",
                lat: 90.0_f64.to_radians(),
                expected: 6399593.62575849331915378570556640625,
                tol: EPS_SMALL,
            },
        ];

        for c in cases {
            let got = WGS84::radius_meridian(c.lat);
            assert!(
                relative_eq!(got, c.expected, epsilon = c.tol),
                "test case: {}, got: {}, expected: {}, epsilon: {}",
                c.name,
                got,
                c.expected,
                c.tol
            );
        }
    }

    #[test]
    fn normal_curvature_radius() {
        // Expected values based on the implementation of WGS84 from AHRS Python package
        // https://ahrs.readthedocs.io/en/latest/wgs84.html#

        struct TestCase {
            name: &'static str,
            lat: f64,
            expected: f64,
            tol: f64,
        }

        let cases = [
            TestCase {
                name: "Equator",
                lat: 0.0,
                expected: 6378137.0,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "50°N",
                lat: 50.0_f64.to_radians(),
                expected: 6390702.044194686226546764373779296875,
                tol: EPS_SMALL,
            },
            TestCase {
                name: "North Pole",
                lat: 90.0_f64.to_radians(),
                expected: 6399593.62575849331915378570556640625,
                tol: EPS_SMALL,
            },
        ];

        for c in cases {
            let got = WGS84::radius_normal(c.lat);
            assert!(
                relative_eq!(got, c.expected, epsilon = c.tol),
                "test case: {}, got: {}, expected: {}, epsilon: {}",
                c.name,
                got,
                c.expected,
                c.tol
            );
        }
    }

    #[test]
    fn earth_rotation_ned_norm() {
        let lat = 52.0;
        let earth_rotation = WGS84::earth_rotation_ned(lat).to_vec();
        let rotation_norm = (earth_rotation.iter().map(|x| x * x).sum::<f64>()).sqrt();
        assert_relative_eq!(rotation_norm, WGS84::W_IE);
    }

    #[test]
    fn gravity_ned_components() {
        let latitude = 50.0_f64.to_radians();
        let elevation = 1000.0;
        let ned = WGS84::gravity_ned(latitude, elevation);
        assert_relative_eq!(
            ned[0],
            WGS84::gravity_north(latitude, elevation),
            epsilon = EPS_SMALL
        );
        assert_relative_eq!(ned[1], 0.0, epsilon = EPS_SMALL);
        assert_relative_eq!(
            ned[2],
            WGS84::gravity_down(latitude, elevation),
            epsilon = EPS_SMALL
        );
    }
}
