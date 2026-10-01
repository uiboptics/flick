#ifndef flick_material_snow_impurity
#define flick_material_snow_impurity

#include "../material.hpp"
#include "../../numeric/function.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class snow_impurity : public base {
    pl_function a_p_;
    pl_function a_cdom_;
    double mass_concentration_; // kg/m3
    double sf_;
    static constexpr double to_nm_{1e9};
    static constexpr double g_to_kg_{1e-3};
  public:
    snow_impurity(const std::string& name, double concentration_scaling_factor=1)
      : sf_{concentration_scaling_factor} {
      std::string p = "/material/snow_impurity/iop_tables";
      a_p_ = read<pl_function>(name+"/a_particles.txt", p);
      a_cdom_ = read<pl_function>(name+"/a_cdom.txt", p);
      mass_concentration_ =
	read<matrix<double>>(name+"/concentration.txt", p).element(0,0)*g_to_kg_;
      a_p_.add_constant_extrapolation();
      a_cdom_.add_constant_extrapolation();
    }
    double mass_concentration() const {
      return mass_concentration_;
    }
    void concentration_scaling_factor(double sf) {
      sf_ = sf;
    }
    double absorption_coefficient() const override {
      double wl = wavelength()*to_nm_;
      return (a_p_.value(wl) + a_cdom_.value(wl)) * sf_;
    }
    double mass_scattering_coefficient() const
    /* [m^2/kg] at 550 nm */
    {
      return 1000;
    }
    double scattering_coefficient() const override {
      double slope = 0.5;
      return mass_scattering_coefficient() * mass_concentration_ * sf_ *
	pow(wavelength()/550e-9, -slope);
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      double g = 0.85;
      double cos_ang = std::clamp<double>(scattering_direction.mu(),-1,1);
      mueller m;
      return m.add(0,0,flick::henyey_greenstein(g).value(cos_ang));
    }
  };
}
}

#endif
