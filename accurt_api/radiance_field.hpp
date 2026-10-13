#ifndef flick_radiance_field
#define flick_radiance_field

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>
#include "../numeric/rotation.hpp"
#include "../numeric/table.hpp"

namespace flick {
enum class irradiance_type {scalar, plane};

class radiance_field {
  pl_table t_;
  std::vector<double> mu_;
  std::vector<double> phi_;
  const double pi = std::numbers::pi;
public:
  radiance_field(const std::vector<double>& theta,
		 const std::vector<double>& phi,
		 const std::vector<std::vector<double>>& r) {
    ensure(r.size() == theta.size());
    ensure(std::is_sorted(theta.begin(), theta.end(), std::less<double>{}));
    ensure(std::is_sorted(phi.begin(), phi.end(), std::less<double>{}));
    ensure(theta.front() >= 0 && theta.back() <= pi);
    ensure(phi.front() >= -pi && phi.back() <= pi);

    mu_ = std::vector<double>(theta.rbegin(), theta.rend());
    for (auto& x : mu_)
      x = std::cos(x);
    phi_ = phi;
    const std::vector<std::vector<double>> reversed_r(r.rbegin(), r.rend());
    t_ = pl_table(mu_, phi_, reversed_r);
  }
  double value(double theta, double phi) const {
    return t_.value(cos(theta),phi);
  }
  auto& rotate_z_axis_to(const unit_vector& direction) {
    quaternion q = rotation_to(direction);
    const radiance_field original = *this;
    for (size_t i = 0; i < t_.row_values().size(); ++i) {
      for (size_t j = 0; j < t_.column_values().size(); ++j) {
	unit_vector u = normalize(rotate(unit_vector(mu_[i], phi_[j]), q));
	t_.set(i,j, original.value(u.mu(), u.phi()));
      }
    }   
    return *this;
  }
  auto& cosine_weight(const unit_vector& normal) {
    for (size_t i = 0; i < mu_.size(); ++i) {
      for (size_t j = 0; j < phi_.size(); ++j) {
	unit_vector r_hat(acos(mu_[i]), phi_[j]);
	const double cosine = std::max(0.0, dot(r_hat, normal));
	t_.set(i,j,t_.value_at(i,j)*cosine);
      }
    }   
    return *this;
  }
  double irradiance(double theta_low, double theta_high,
		    double phi_low, double phi_high,
		    irradiance_type type = irradiance_type::scalar) const {
    ensure(theta_low <= theta_high && phi_low <= phi_high);
    ensure(theta_low >= 0 && theta_high <= std::numbers::pi);
    ensure(phi_low >= -std::numbers::pi && phi_high <= std::numbers::pi);
    radiance_field r = *this;
    if (type == irradiance_type::plane)
      r.cosine_weight({0,0,1});
    return r.integral(theta_low, theta_high, phi_low, phi_high);
  }
private:
  static void ensure(bool b) {
    if (!b) {
      throw std::runtime_error("radiance_field");
    }
  }  
  double integral(double theta_low, double theta_high,
		  double phi_low, double phi_high) const {
    const double mu_low = std::cos(theta_high);
    const double mu_high = std::cos(theta_low);
    return t_.integral(mu_low, mu_high, phi_low, phi_high);
  }

 
};

}

#endif

