use nalgebra::{self as na, UnitQuaternion};

pub fn skew_symmetric(vec: &na::Vector3<f64>) -> na::Matrix3<f64> {
    #[rustfmt::skip]
    let mat = na::Matrix3::new(
        0.0, -vec.z, vec.y,
      vec.z, 0.0, -vec.x,
     -vec.y, vec.x, 0.0,
    );
    mat
}

/// Computes quaternion representing orientation change from angular displacement vector.
///
/// Converts an angular displacement vector to a unit quaternion
/// using Rodrigues' formula with small-angle approximation for numerical stability.
pub fn quaternion_increment(angle_vec: na::Vector3<f64>) -> na::UnitQuaternion<f64> {
    const SMALL_ANGLE: f64 = 1e-8;
    let angle = angle_vec.norm();
    let q = if angle < SMALL_ANGLE {
        // first order approximation for small angles to avoid division by zero
        let q0 = 1.0 - 0.125 * angle * angle;
        let qv = 0.5 * angle_vec;
        na::Quaternion::from_parts(q0, qv)
    } else {
        // Rodrigues' formula for quaternion from angle axis vector
        let half = angle * 0.5;
        let q0 = f64::cos(half);
        let qv = angle_vec * f64::sin(half) / angle;
        na::Quaternion::from_parts(q0, qv)
    };
    UnitQuaternion::from_quaternion(q)
}

/// Compute a new a new quaternion orientation from a given orientation (as quaternion)and an angular displacement vector.
pub fn quaternion_rotation(
    q: na::UnitQuaternion<f64>,
    angle_vec: na::Vector3<f64>,
) -> na::UnitQuaternion<f64> {
    let q_from_angle_vec = quaternion_increment(angle_vec);
    UnitQuaternion::from_quaternion((q * q_from_angle_vec).normalize())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn skewsymmetric_matrix_from_vec() {
        assert_relative_eq!(
            skew_symmetric(&na::Vector3::zeros()),
            na::Matrix3::zeros(),
            epsilon = 1e-16
        );

        #[rustfmt::skip]
        assert_relative_eq!(
            skew_symmetric(&na::Vector3::new(1.0, 2.0, 3.0)),
            na::Matrix3::new(
                0.0, -3.0, 2.0,
                3.0, 0.0, -1.0,
               -2.0, 1.0, 0.0),
            epsilon = 1e-16
        );
    }
}
