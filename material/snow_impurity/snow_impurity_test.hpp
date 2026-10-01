#include "snow_impurity.hpp"

namespace flick {
  begin_test_case(snow_impurity_test) {
    material::snow_impurity si("EIK1");
    si.set_wavelength(550e-9);
    si.concentration_scaling_factor(2);
    check_close(si.mass_concentration(), 4.920*1e-3, 0.1_pct);
    check_close(si.absorption_coefficient(), (0.287872 -0.000192)*2, 0.1_pct);
    check_close(si.scattering_coefficient(),4.920e-3*2*1000,0.1_pct);    
  } end_test_case()
}
