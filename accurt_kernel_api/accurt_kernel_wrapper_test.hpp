#include "../environment/unit_test.hpp"
#include "accurt_kernel_wrapper.hpp"
#include <numbers>

namespace flick {
  begin_test_case(accurt_kernel_wrapper_test_A) {
    // Compare with van de Hulst 1980, vol 1, chapter 9, table 12,
    // p258, FLUX U and AVERAGE N for mu_0 = 1.0
    accurt_kernel_wrapper a;
    a.streams_upper_slab(6);
    a.insert_user_optical_depth(0);
    a.bottom_albedo(0);
    a.refractive_index(1.0);
    check(a.number_of_layers() == 0);
    a.insert_layer_at_bottom(1.0, 1.0, {1.0, 0.0, 0.0, 0.0, 0.0});
    check(a.number_of_layers() == 1);
    
    a.solve_rte();
    check_close(a.plane_irradiance_up()[0], 0.34133, 0.01_pct);
    check_close(a.scalar_irradiance_up()[0]/2, 0.37869, 0.01_pct);
  } end_test_case()

  begin_test_case(accurt_test_B) {
    // Compare with van de Hulst 1980, vol 2, chapter 13, table 35,
    // p431, FLUX
    double ssalb = 0.9;
    double od = 1;
    std::vector<double> p0_75 =
      {
	1.0,
	0.75,
	0.5625,
	0.421875,
	0.31640625,
	0.2373046875,
	0.177978515625
      };
    accurt_kernel_wrapper a;
    a.streams_upper_slab(6);
    a.insert_user_optical_depth(od);
    a.bottom_albedo(0);
    a.refractive_index(1.0);
    a.insert_layer_at_bottom(od, ssalb, p0_75);
    a.solve_rte();
    check_close(a.plane_irradiance_down()[0], 0.82792, 0.02_pct);
  } end_test_case()
    
  begin_test_case(accurt_test_C) {
    // Van de Hulst, chapter 13, page 418 and 419, mu_0=0.9
    double ssalb = 1;
    double od = 1;
    std::vector<double> p0_5 =
      {
	1.0,
	0.5,
	0.25,
	0.125,
	0.0625,
	0.03125,
	0.015625
      };
    accurt_kernel_wrapper a;
    double ang = std::numbers::pi-acos(0.9);
    a.source_zenith_angle(ang);
    a.streams_upper_slab(6);
    a.insert_user_optical_depth(0);
    a.bottom_albedo(0);
    a.refractive_index(1.0);
    a.insert_layer_at_bottom(od, ssalb, p0_5);
    a.solve_rte();
    double albedo = a.plane_irradiance_up()[0]/a.plane_irradiance_down()[0];
    check_close(albedo, 0.20197, 0.04_pct);
  } end_test_case()
  
  begin_test_case(accurt_test_D) {
    // Van de Hulst, chapter 13, page 418 and 419, mu_0=0.9
    accurt_kernel_wrapper a;
    a.insert_user_optical_depth(0);
    a.bottom_albedo(1);
    a.refractive_index(1.0);
    a.solve_rte();
    double albedo = a.plane_irradiance_up()[0]/a.plane_irradiance_down()[0];
    check_close(albedo, 1, 0.001_pct);
    check_close(a.scalar_irradiance_up()[0], 2.0, 1e-5_pct);
    check_close(a.scalar_irradiance_down()[0], 1.0, 1e-5_pct);
    check_close(a.beam_plane_irradiance_down()[0], 1, 1e-5_pct);
    check_small(a.beam_plane_irradiance_up()[0], 1e-9);
  } end_test_case()
    
  begin_test_case(accurt_test_E) {
    // Transparent refractive slab with Lambert white bottom
    // Test against AI-generated transmittance to bottom.
    accurt_kernel_wrapper a;
    double od1 = 1e-12;
    double od2 = od1;
    double ssa = 1;
    a.insert_layer_at_bottom(od1,ssa,{1,0});
    a.insert_layer_at_bottom(od2,ssa,{1,0});
    a.layers_upper_slab(1);
    a.source_zenith_angle(0);      
    a.insert_user_optical_depth(0);
    a.insert_user_optical_depth(od1+od2);
    a.streams_upper_slab(12);
    a.bottom_albedo(1);
    a.refractive_index(1.3);
    a.solve_rte();
    double T_bottom = 1.7694;
    check_close(a.plane_irradiance_up()[0], 1, 0.001_pct);
    check_close(a.plane_irradiance_down()[1], T_bottom, 0.001_pct);
  } end_test_case() 
}
