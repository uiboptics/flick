#include "tabulated_iops.hpp"

namespace flick {
  begin_test_case(tabulated_iops_test) {
    double wl = 550e-9;
    double volume_fraction = 1e-4;
    std::string material_name = "small_marine_bubbles";
    material::tabulated_iops ti(material_name, volume_fraction);
    ti.set_wavelength(wl);
    check(ti.absorption_coefficient() < 0);
    
  } end_test_case()
}
