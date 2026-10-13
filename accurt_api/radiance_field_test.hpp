#include "../environment/unit_test.hpp"
#include "radiance_field.hpp"

namespace flick {
  begin_test_case(radiance_field_test_A) {
    const double pi = std::numbers::pi;
    radiance_field rf({0,pi/2,pi},{-pi,0,pi},{{1,1,1},
					      {1,1,1},
					      {1,1,1}
      });
    check_close(rf.irradiance(0,pi,-pi,pi),4*pi, 0.1_pct);
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::scalar),4*pi, 1e-6_pct);
    rf.rotate_z_axis_to(unit_vector(pi/4,pi/3));
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::scalar),4*pi, 1e-6_pct);
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::plane),pi, 1e-6_pct);
    check_small(rf.irradiance(pi/2,pi,-pi,pi,irradiance_type::plane),1e-9);

  } end_test_case()
  
}
