# wgs84-gravity

WGS84 gravity model constants and formulas (radii, normal gravity, Earth rotation).

All data inputs and outputs are in SI units.

Implementation is based on the following references
- [Aided Navigation: GPS with High Rate Sensors (Jay A. Farrell)](https://books.google.nl/books?id=yNujEvIMszYC&lpg=PP1&pg=PR3#v=onepage&q&f=false)
- [Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems (GNSS Technology and Applications)](https://books.google.nl/books?id=t94fAgAAQBAJ&printsec=copyrigh)
- [International Gravity Formula (Walter D. Lambert)](https://earth.geology.yale.edu/~ajs/1945A/360.pdf)
- [AHRS: WGS84 documentation](https://ahrs.readthedocs.io/en/latest/wgs84.html)
- [Wikipedia: Theoretical gravity - Somigliana equation](https://en.wikipedia.org/wiki/Theoretical_gravity#Somigliana_equation)
- [AHRS Python implementation (wgs84.py)](https://github.com/Mayitzin/ahrs/blob/b179ad0449c6da5da4780533d7cc9bd522c3ef87/ahrs/utils/wgs84.py#L553)

## Install
```toml
[dependencies]
wgs84-gravity = "0.1"
```

## Example
```rust
use wgs84_gravity::WGS84;

let lat = 52.0_f64.to_radians();
let alt = 100.0; // meters

// `gravity_normal` is and alias for `gravity_down`
let g = WGS84::gravity_down(lat, alt);
let g = WGS84::gravity_normal(lat, alt);

let g_surface = WGS84::gravity_surface(lat);
let g_ned = WGS84::gravity_ned(lat, alt); // north, east, down [m/s^2]

let r_m = WGS84::radius_meridian(lat);
let r_n = WGS84::radius_normal(lat);
let w_ie_ned = WGS84::earth_rotation_ned(lat);
```

Docs: [docs.rs/wgs84-gravity](https://docs.rs/wgs84-gravity) 

License: MIT
