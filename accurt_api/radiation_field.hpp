#ifndef flick_radiation_field
#define flick_radiation_field

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>
#include "../numeric/rotation.hpp"
#include "../numeric/table.hpp"
#include "../numeric/vector.hpp"

namespace flick {
enum class irradiance_type {scalar, plane};

class radiation_field {
  pl_table radiance_; //[W/m2/m/sr]
  std::vector<double> mu_;
  std::vector<double> phi_;
  std::vector<flick::vector> dirac_delta_beam_; // [W/m2/m] normal to propagation
  static constexpr double pi = std::numbers::pi;
public:
  radiation_field(const std::vector<double>& theta,
		  const std::vector<double>& phi,
		  const std::vector<std::vector<double>>& r,
		  const std::vector<flick::vector>& dirac_delta_beam) {
    ensure(r.size() == theta.size());
    ensure(std::is_sorted(theta.begin(), theta.end(), std::greater<double>{}));
    ensure(std::is_sorted(phi.begin(), phi.end(), std::less<double>{}));
    mu_ = std::vector<double>(theta.begin(), theta.end());
    for (auto& x : mu_)
      x = std::clamp<double>(std::cos(x),-1,1);
    phi_ = phi;
    radiance_ = pl_table(mu_, phi_, r);
    dirac_delta_beam_ = dirac_delta_beam;
  }
  double radiance(double theta, double phi) const {
    return radiance_.value(cos(theta),phi);
  }
  const auto& dirac_delta_beam() const {
    return dirac_delta_beam_;
  }
  static bool beam_inside(const vector& beam,
                        double theta_low, double theta_high,
                        double phi_low, double phi_high) {
    if (norm(beam) == 0)
      return false;
    unit_vector direction = normalize(beam);
    double theta = direction.theta();
    double phi = std::remainder(direction.phi() - pi, 2*pi);
    bool theta_inside = theta >= theta_low && theta <= theta_high;
    bool phi_inside =
      (phi >= phi_low && phi <= phi_high) ||
      (phi + 2*pi >= phi_low && phi + 2*pi <= phi_high) ||
      (phi - 2*pi >= phi_low && phi - 2*pi <= phi_high);  
    return theta_inside && phi_inside;
  }
  auto& apply_cosine_weight(const unit_vector& normal) {
    for (size_t i = 0; i < mu_.size(); ++i) {
      for (size_t j = 0; j < phi_.size(); ++j) {
	unit_vector r_hat(acos(mu_[i]), phi_[j]+pi);
	const double cosine = std::max(0.0, dot(r_hat, normal));
	radiance_.set(i,j,radiance_.value_at(i,j)*cosine);
      }
    }
    for (auto& x : dirac_delta_beam_) {
      x *= std::max(0.0, dot(normalize(x), normal));
    }
    return *this;
  }
  double irradiance(double theta_low, double theta_high,
		    double phi_low, double phi_high,
		    irradiance_type type,
		    const unit_vector& normal = unit_vector(0,0)) const {
    radiation_field r = *this;
    if (type == irradiance_type::plane)
      r.apply_cosine_weight(normal);
    double b = 0;
    for (size_t i = 0; i < dirac_delta_beam_.size(); ++i) {
      if (beam_inside(dirac_delta_beam_[i],
		      theta_low, theta_high,
		      phi_low, phi_high)) {
	b += norm(r.dirac_delta_beam()[i]);
      }
    }
    return r.integral(theta_low, theta_high, phi_low, phi_high) + b;
  }
  auto& rotate_z_axis_to(const unit_vector& direction) {
    quaternion q = rotation_to(direction);
    const radiation_field original = *this;
    for (size_t i = 0; i < radiance_.row_values().size(); ++i) {
      for (size_t j = 0; j < radiance_.column_values().size(); ++j) {
	unit_vector u = normalize(rotate(unit_vector(acos(mu_[i]), phi_[j]+pi), q));
	radiance_.set(i,j, original.radiance(acos(u.mu()), u.phi()-pi));
      }
    }
    for (auto& x : dirac_delta_beam_) {
      x = rotate(x, inv(q));
    }
    return *this;
  }
  friend std::ostream& operator<<(std::ostream &os,
				    const radiation_field& r) {
    for (auto& x : r.dirac_delta_beam_) {
      os << x << "\n";
    }
    os << r.radiance_ << std::endl;
    return os;
  }
private:
  static void ensure(bool b) {
    if (!b) {
      throw std::runtime_error("radiation_field");
    }
  }  
  double integral(double theta_low, double theta_high,
		  double phi_low, double phi_high) const {
    ensure(theta_low <= theta_high && phi_low <= phi_high);
    ensure(theta_low >= 0 && theta_high <= pi);
    ensure(phi_low >= -pi && phi_high <= pi);
    const double mu_low = std::cos(theta_high);
    const double mu_high = std::cos(theta_low);
    return radiance_.integral(mu_low, mu_high, phi_low, phi_high);
  }
};
}

#endif

