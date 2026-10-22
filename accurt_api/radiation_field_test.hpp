#include "../environment/unit_test.hpp"
#include "radiation_field.hpp"

namespace flick {
  begin_test_case(radiation_field_test_A) {
    // Isotropic radiance, no beams
    const double pi = std::numbers::pi;
    radiation_field rf({pi,pi/2,0},{-pi,0,pi},{{1,1,1},
					       {1,1,1},
					       {1,1,1}}, {});
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::scalar),4*pi, 1e-6_pct);
    rf.rotate_z_axis_to(unit_vector(pi/4,pi/3));
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::scalar),4*pi, 1e-6_pct);
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::plane),pi, 1e-6_pct);
    check_close(rf.irradiance(pi/2,pi,-pi,pi,irradiance_type::plane,
			     unit_vector(0,0,-1)),pi,1e-6_pct);
    check_small(rf.irradiance(pi/2,pi,-pi,pi,irradiance_type::plane),1e-9);
  } end_test_case()
  
  begin_test_case(radiation_field_test_B) {
    // Isotropic radiance with beams
    const double pi = std::numbers::pi;
    radiation_field rf({pi,pi/2,0},{-pi,0,pi},{{1,1,1},
					       {1,1,1},
					       {1,1,1}}, {{0,0,-1},{0,0,1},{-1,0,0}});
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::plane),pi+1, 1e-6_pct);
    rf.rotate_z_axis_to(unit_vector{1,0,0});
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::plane),pi, 1e-6_pct);
  } end_test_case()
  
  begin_test_case(radiation_field_test_C) {
    // More radiance up than down, no beams
    const double pi = std::numbers::pi;
    radiation_field rf({pi,pi/2,0},{-pi,0,pi},{{0,0,0},
					       {1,1,1},
					       {1,1,1}}, {});
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::plane),pi, 1e-6_pct);
    check_small(rf.radiance(pi,0),1e-9);
    rf.rotate_z_axis_to(unit_vector{0,0,1}); // No rotation
    check_small(rf.radiance(pi,0),1e-9); 
    rf.rotate_z_axis_to(unit_vector{0,0,-1});
    check_small(rf.radiance(0,0),1e-9); 
    check_small(rf.irradiance(0,pi,-pi,pi,irradiance_type::plane), 1e-9);
  } end_test_case()
  
  begin_test_case(radiation_field_test_D) {
    /* Beam outside acceptance angle should not contribute */
    const double pi = std::numbers::pi;
    radiation_field rf({pi,pi/2,0},{-pi,0,pi},{{0,0,0},
					       {0,0,0},
					       {0,0,0}}, {{0,0,1}});
    check_small(rf.irradiance(pi/4,pi/2,0,pi,irradiance_type::scalar),1e-9);
  } end_test_case()

  begin_test_case(radiation_field_test_E) {
    // Ascending polar angles are accepted without changing the field.
    const double pi = std::numbers::pi;
    radiation_field rf({0,pi/2,pi},{-pi,0,pi},{{1,1,1},
					       {1,1,1},
					       {0,0,0}}, {});
    check_close(rf.irradiance(0,pi,-pi,pi,irradiance_type::plane),pi, 1e-6_pct);
    check_close(rf.radiance(0,0),1,1e-6_pct);
    check_small(rf.radiance(pi,0),1e-9);
  } end_test_case()
}
