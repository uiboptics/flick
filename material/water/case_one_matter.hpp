#ifndef flick_material_case_one_matter
#define flick_material_case_one_matter

#include "../material.hpp"
#include "cdom.hpp"
#include "phytoplankton.hpp"
#include <algorithm>
#include <cmath>

namespace flick {
namespace material {
  class case_one_matter : public base {
    phytoplankton phy_;
    cdom cdom_;
  public:
    explicit case_one_matter(double chl_concentration = 1e-6)
      : phy_{chl_concentration, 0.97}, cdom_{cdom_absorption_440(), 0.014} {
      set_wavelength(wavelength());
    }
    double absorption_coefficient() const override {
      return cdom_.absorption_coefficient() + phy_.absorption_coefficient();
    }
    double scattering_coefficient() const override {
      return phy_.scattering_coefficient();
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      return phy_.mueller_matrix(scattering_direction);
    }
    void set_wavelength(double wl) override {
      base::set_wavelength(wl);
      phy_.set_wavelength(wl);
      cdom_.set_wavelength(wl);
    }
    void set(const flick::pose& p) override {
      base::set(p);
      phy_.set(p);
      cdom_.set(p);
    }
  private:
    double cdom_absorption_440() {
      /*
      Morel & Maritorena (2001), Appendix B, equation 18.
      */
      phy_.set_wavelength(440e-9);
      return 0.2 * phy_.absorption_coefficient();
    }
  };
}
}

#endif

