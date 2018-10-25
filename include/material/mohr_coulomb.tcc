//! Constructor with id and material properties
template <unsigned Tdim>
mpm::MohrCoulomb<Tdim>::MohrCoulomb(unsigned id,
                                    const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties["density"].template get<double>();
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    friction_angle_ =
        material_properties["friction_angle"].template get<double>();
    dilation_angle_ =
        material_properties["dilation_angle"].template get<double>();
    cohesion_ = material_properties["cohesion"].template get<double>();
    tension_cutoff_ =
        material_properties["tension_cutoff"].template get<double>();
    porosity_ = material_properties["porosity"].template get<double>();
    properties_ = material_properties;
    // Calculate bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
    shear_modulus_ = youngs_modulus_ / (2.0 * (1 + poisson_ratio_));

    // Set elastic tensor
    this->compute_elastic_tensor();

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Return elastic tensor
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  // clang-format off
  // compute elasticityTensor
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on
  return true;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::MohrCoulomb<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr) {

  const double PI = std::atan(1.0) * 4.;

  // Friction and dilation in radians
  const double phi = friction_angle_ * PI / 180.;
  const double psi = dilation_angle_ * PI / 180.;

  // Elastic update of stresses
  Vector6d stress_update = stress + this->de_ * dstrain;

  const double fx = stress_update(0);
  const double fy = stress_update(1);
  const double fz = stress_update(2);
  const double fxy = stress_update(3);

  const double mohr_c = 0.5 * (fx + fy);
  const double mohr_r = 0.5 * sqrt((fx - fy) * (fx - fy) + 4.0 * fxy * fxy);

  const double f1 = mohr_c - mohr_r;
  const double f2 = mohr_c + mohr_r;

  // Stressess
  Eigen::Matrix<double, 3, 1> sigma;

  unsigned mohr_flag;
  if (f1 > fz) {
    sigma(0) = fz;
    sigma(1) = f1;
    sigma(2) = f2;
    mohr_flag = 2;
  } else if (f2 < fz) {
    sigma(0) = f1;
    sigma(1) = f2;
    sigma(2) = fz;
    mohr_flag = 3;
  } else {
    sigma(0) = f1;
    sigma(1) = fz;
    sigma(2) = f2;
    mohr_flag = 1;
  }

  // Mohr-Coulomb failure criteria
  const double nphi = (1.0 + sin(phi)) / (1.0 - sin(phi));
  const double npsi = (1.0 + sin(psi)) / (1.0 - sin(psi));
  const double f_s = sigma(0) - sigma(2) * nphi + 2.0 * cohesion_ * sqrt(nphi);
  const double f_t = tension_cutoff_ - sigma(2);
  const double alphap = sqrt(1.0 + nphi * nphi) + nphi;
  const double heichi = sigma(2) - tension_cutoff_ +
                        alphap * (sigma(0) - nphi * tension_cutoff_ +
                                  2.0 * cohesion_ * sqrt(nphi));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * shear_modulus_;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * shear_modulus_;


  double lamda_s, lamda_t;
  if (heichi < 0.) {
    if (f_s < 0.) {
      // tens + shear failure
      lamda_s = f_s / ((a1 - a2 * npsi) - (a2 - a1 * npsi) * nphi);
      sigma(0) -= lamda_s * (a1 - a2 * npsi);
      sigma(1) -= lamda_s * a2 * (1.0 - npsi);
      sigma(2) -= lamda_s * (a2 - a1 * npsi);
    }
  } else if (heichi > 0.) {
    if (f_t < 0.) {
      // tens failure
      lamda_t = f_t / a1;
      sigma(0) += lamda_t * a2;
      sigma(1) += lamda_t * a2;
      sigma(2) += lamda_t * a1;
      tension_cutoff_ = 0.0;
    }
  }

  double cs2, si2, dc2, dss;
  if (sigma(0) == sigma(2)) {
    cs2 = 1.0;
    si2 = 0.0;
  } else {
    cs2 = (fx - fy) / (f1 - f2);
    si2 = 2.0 * fxy / (f1 - f2);
  }

  if (mohr_flag == 1) {
    dc2 = (sigma(0) - sigma(2)) * cs2;
    dss = sigma(0) + sigma(2);
    stress_update(0) = 0.5 * (dss + dc2);
    stress_update(1) = 0.5 * (dss - dc2);
    stress_update(2) = sigma(1);
    stress_update(3) = 0.5 * (sigma(0) - sigma(2)) * si2;
    stress_update(4) = 0.0;
    stress_update(5) = 0.0;
  } else if (mohr_flag == 2) {
    dc2 = (sigma(1) - sigma(2)) * cs2;
    dss = sigma(1) + sigma(2);
    stress_update(0) = 0.5 * (dss + dc2);
    stress_update(1) = 0.5 * (dss - dc2);
    stress_update(2) = sigma(0);
    stress_update(3) = 0.5 * (sigma(1) - sigma(2)) * si2;
    stress_update(4) = 0.0;
    stress_update(5) = 0.0;
  } else {
    dc2 = (sigma(0) - sigma(1)) * cs2;
    dss = sigma(0) + sigma(1);
    stress_update(0) = 0.5 * (dss + dc2);
    stress_update(1) = 0.5 * (dss - dc2);
    stress_update(2) = sigma(2);
    stress_update(3) = 0.5 * (sigma(0) - sigma(1)) * si2;
    stress_update(4) = 0.0;
    stress_update(5) = 0.0;
  }

  return stress_update;
}
