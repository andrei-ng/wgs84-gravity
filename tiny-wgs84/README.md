# Tiny WGS84 Model

A tiny WGS84 model for calculating Gravity and Earth Radii using the WGS84 reference ellipsoid.

Implementation is based on the following references
- [Aided Navigation: GPS with High Rate Sensors (Jay A. Farrell)](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false)
- [Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems (GNSS Technology and Applications)](https://books.google.nl/books?id=t94fAgAAQBAJ&printsec=copyrigh)
- [International Gravity Formula (Walter D. Lambert)](https://earth.geology.yale.edu/~ajs/1945A/360.pdf)
- [AHRS: WGS84 documentation](https://ahrs.readthedocs.io/en/latest/wgs84.html)
- [Wikipedia: Theoretical gravity - Somigliana equation](https://en.wikipedia.org/wiki/Theoretical_gravity#Somigliana_equation)
- [AHRS Python implementation (wgs84.py)](https://github.com/Mayitzin/ahrs/blob/b179ad0449c6da5da4780533d7cc9bd522c3ef87/ahrs/utils/wgs84.py#L553)

 ## Features

 - **Default API**: Functions accept and return `f64` values in SI units: latitude and longitude in `[rad]`, altitude in `[m]`, and acceleration in `[m/s²]`
 - **Units API** (feature `units`): Functions accept and return strongly typed values: `Angle`, `Length` and `Acceleration` (based on `uom` library).

 See [uom](https://docs.rs/uom/latest/uom/) for more information about the units.

 Recommendation is to use the strongly typed, `units` API, with `uom` crate for type safety.

 ## Example

 ```rust
#[cfg(not(feature = "units"))]
{
    use tiny_wgs84::WGS84;

    // Default `raw` API 
    let lat_rad = 52.0_f64.to_radians();
    let alt_m = 1000.0; 
    let g = WGS84::gravity_down(lat_rad, alt_m);
}
 ```

```rust
// Strongly typed API: when compiled with --features units 
#[cfg(feature = "units")]
{
    use tiny_wgs84::WGS84;
    use uom::si::f64::{Angle, Length, Acceleration};
    use uom::si::angle::radian;
    use uom::si::length::meter;

    let lat = Angle::new::<radian>(52.0_f64.to_radians());
    let alt = Length::new::<meter>(1000.0);
    let g: Acceleration = WGS84::gravity_down(lat, alt);
    // `gravity_normal` is an alias  for `gravity_down`
    let g: Acceleration = WGS84::gravity_normal(lat, alt);

    let g_surface: Acceleration = WGS84::gravity_surface(lat);
    let g_ned: [Acceleration; 3] = WGS84::gravity_ned(lat, alt); // north, east, down [m/s^2]

    let r_m: Length = WGS84::radius_meridian(lat);
    let r_n: Length = WGS84::radius_normal(lat);
}
```

Docs: [docs.rs/strapdown-ins/tiny-wgs84](https://docs.rs/strapdown-ins/tiny-wgs84) 

License: MIT
