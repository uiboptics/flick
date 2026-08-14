#ifndef flick_basic_monodispersed
#define flick_basic_monodispersed

#include "../numeric/constants.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/function.hpp"

namespace flick {
  class basic_monodispersed_mie {
  protected:
    const double pi_{constants::pi};

    stdcomplex m_host_;
    stdcomplex m_sphere_;
    double vacuum_wl_;
    double radius_{1e-6};
    stdvector angles_{0};
    stdcomplex wavenumber_in_host_{2*pi_*m_host_/vacuum_wl_};
    stdcomplex m_sphere_at_r0_;
    
    stdcomplex size_parameter_in_host() const {
      return wavenumber_in_host_ * radius_;
    }
    double size_parameter_in_vacuum() const {
      return 2*pi_*radius_/vacuum_wl_;
    }
    stdcomplex relative_refractive_index() const {
      return m_sphere_ / m_host_;
    }
    double geometrical_cross_section() const {
      return pi_ * pow(radius_,2);
    }
    double geometrical_extinction_cross_section() const {
      return 2 * geometrical_cross_section();
    }
    double parameterized_absorption_efficiency() const {
      stdcomplex n = relative_refractive_index();
      stdcomplex arg = 1./n * (pow(n,3) - pow(pow(n,2)-1., 3./2));
      double Qa0 = 8./3*imag(m_sphere_-m_host_) * real(size_parameter_in_host())
	* std::abs(arg);
      return 0.94 * tanh(Qa0/0.94);
    }
    double parameterized_absorption_cross_section() const {
      return parameterized_absorption_efficiency() * geometrical_cross_section();
    }
    void update_wavenumber_in_host() {
      wavenumber_in_host_ = 2*pi_*m_host_/vacuum_wl_;
    }
  public:
    basic_monodispersed_mie(const stdcomplex& m_host,
			     const stdcomplex& m_sphere,
			     double vacuum_wl)
      : m_host_{m_host}, m_sphere_{m_sphere}, vacuum_wl_{vacuum_wl},
	m_sphere_at_r0_{m_sphere} {
    }
    virtual ~basic_monodispersed_mie() = default;
    virtual double absorption_cross_section() const = 0;
    virtual double scattering_cross_section() const = 0;
    virtual stdvector scattering_matrix_element(size_t row,
						size_t col) const = 0;
    virtual void radius(double r) = 0;
    virtual void angles(const stdvector& a) {
      angles_ = a;
    }
    double radius() const {
      return radius_;
    }
    stdvector angles() const {
      return angles_;
    }
    void set_wavelength(double wl) {
      vacuum_wl_ = wl;
      update_wavenumber_in_host();
    }
    void set_refractive_indices(stdcomplex m_host, stdcomplex m_sphere) {
      m_host_ = m_host;
      m_sphere_ = m_sphere;
      m_sphere_at_r0_ = m_sphere;
      update_wavenumber_in_host();
    }
  };  
}

#endif
