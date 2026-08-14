#ifndef flick_parameterized_monodispersed_mie
#define flick_parameterized_monodispersed_mie

#include "basic_monodispersed_mie.hpp"
#include "../numeric/physics_function.hpp"
#include "../environment/input_output.hpp"

namespace flick { 
  class parameterized_monodispersed_mie : public basic_monodispersed_mie
  // Approximate Mie-code solutions for large spheres, see: Stamnes,
  // K., Hamre, B., Stamnes, J.J., Ryzhikov, G., Biryulina, M.,
  // Mahoney, R., Hauss, B. and Sei, A., 2011. Modeling of radiation
  // transport in coupled atmosphere-snow-ice-ocean systems. Journal
  // of Quantitative Spectroscopy and Radiative Transfer, 112(4),
  // pp.714-726. Some modifications to be consistent with full Mie
  // theory, which gives negative absorption for transparent particles
  // in an absorbing host medium - gives correctly less total
  // absorption when e.g. bubbles are added in water.
  {
    double Qa_{0};
    pl_function g0_ = read<pl_function>("mie/g_parameterized.txt");
    
    void update_efficiency() {
      Qa_ = parameterized_absorption_efficiency();
    }
    double asymmetry_factor() const {
      double n_r = real(relative_refractive_index());
      return pow(g0_.value(n_r), pow(1-Qa_, 0.6));
    }
  public:
    using basic_monodispersed_mie::basic_monodispersed_mie;
    void radius(double r) override {
      radius_ = r;
      update_efficiency();
    }
    double absorption_cross_section() const override {
      return Qa_ * geometrical_cross_section();
    }
    double scattering_cross_section() const override {
      return geometrical_extinction_cross_section() - absorption_cross_section();
    }
    stdvector scattering_matrix_element(size_t row, size_t col) const override
    // Note that integrating element 0,0 over all 4*pi solid angles gives
    // the scattering cross section, where we count from 0 instead of 1.
    {
      if (row == 0 && col == 0) {
	stdvector hg(angles_.size());
	double g = asymmetry_factor();
	for (size_t i=0; i<hg.size(); ++i) {
	  hg[i] = henyey_greenstein(g).phase_function(angles_[i]);
	}
	return hg * scattering_cross_section();
      }
      else
	return stdvector(angles_.size(),0);
    }    
  };
}

#endif
  
