#include "parameterized_monodispersed_mie.hpp"

namespace flick {  
  begin_test_case(p_mono_mie_test_A) {
    stdcomplex m_host{1.3,0};
    stdcomplex m_sphere{1.0,0};
    parameterized_monodispersed_mie mono_mie(m_host,
					     m_sphere,
					     500e-9);
    double r = 0.001;
    mono_mie.radius(r);
    check_close(mono_mie.scattering_cross_section(),2*constants::pi*pow(r,2),0.1);
    check_small(mono_mie.absorption_cross_section(),1e-12);
  } end_test_case()
  
  begin_test_case(p_mono_mie_test_B) {
    stdcomplex m_host{1.3,0.001};
    stdcomplex m_sphere{1.0,0};
    parameterized_monodispersed_mie mono_mie(m_host,
					     m_sphere,
					     500e-9);
    double r = 1e-6;
    mono_mie.radius(r);
  } end_test_case()
}
