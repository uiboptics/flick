#include "single_layer_slab.hpp"

namespace flick {  
  begin_test_case(single_layer_slab_test) {
    using namespace flick;
    absorption_coefficient a(0);
    scattering_coefficient b(0);
    asymmetry_factor g(0);
    material::henyey_greenstein hg(a,b,g);
  
    incidence_angle theta_0(0);
    number_of_packages np(1e3);

    model::single_layer_slab sl_slab(thickness(1),hg,theta_0,np);
    sl_slab.set(bottom_albedo(0));
    check_small(sl_slab.hemispherical_reflectance(),1e-12);
    sl_slab.set(bottom_albedo(1));
    check_close(sl_slab.hemispherical_reflectance(),1,1e-12);
    sl_slab.set(bottom_albedo(0.5));
    check_close(sl_slab.hemispherical_reflectance(),0.5,10);
  
  } end_test_case()
}