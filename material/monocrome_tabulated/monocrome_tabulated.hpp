#ifndef flick_material_monocrome_tabulated
#define flick_material_monocrome_tabulated

#include "../material.hpp"
#include "../../environment/input_output.hpp"
#include "../../numeric/function.hpp"
#include "../../polarization/mueller.hpp"

#include <string>
#include <numbers>

namespace flick {
namespace material {
  class monocrome_tabulated : public base {
  protected:
    double ac_;
    double sc_;
    double n_;
    tabulated_phase_function p_;
  public:
    monocrome_tabulated(double ac, double sc,
			const std::string& phase_function_name,
			double n=1)
      : ac_{ac}, sc_{sc}, n_{n} {
      auto f = read<pe_function>(phase_function_name + ".txt",
				 "/material/monocrome_tabulated/phase_functions");
      f.scale_x(std::numbers::pi/180);
      p_ = tabulated_phase_function(to_cos_x(f),4);
      p_.normalize();
      std::cout << p_.integral_4pi();
    }
    double absorption_coefficient() const override {
      return ac_;
    }
    double asymmetry_factor() {
      return p_.asymmetry_factor();
    }
    double scattering_coefficient() const override {
      return sc_;
    }
    double real_refractive_index() const override {
      return n_;
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      mueller m;
      return m.add(0,0, p_.value(cos(angle(scattering_direction))));
    }
  }; 
}
}

#endif
