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
    a.insert_layer_at_bottom(1.0, 0, {1.0, 0.0, 0.0, 0.0, 0.0});
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
    a.insert_layer_at_bottom(od * ssalb, od * (1 - ssalb), p0_75);
    a.solve_rte();
    check_close(a.plane_irradiance_down()[0], 0.82792, 0.02_pct);
  } end_test_case()
    /*
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
    a.source_zenith_angle(std::numbers::pi-acos(0.9));
    a.streams_upper_slab(6);
    a.insert_user_optical_depth(0);
    a.bottom_albedo(0);
    a.refractive_index(1.0);
    a.insert_layer_at_bottom(od * ssalb, od * (1 - ssalb), p0_5);
    a.solve_rte();
    check_close(a.plane_irradiance_up()[0], 0.20197, 0.02_pct);
  } end_test_case()
    */
  begin_test_case(accurt_test_D) {
    /* ChatGPT solution of transparent refractive slab with Lambert white
       bottom. */
    double ssalb = 0;
    double od = 1e-6;
    accurt_kernel_wrapper a;
    a.source_zenith_angle(0);
    a.streams_upper_slab(6);
    a.insert_user_optical_depth(od);
    a.bottom_albedo(1);
    a.refractive_index(1.3);
    a.insert_layer_at_bottom(od * ssalb, od * (1 - ssalb), {1,0});
    a.solve_rte();
    double R_fresnel = 0.443566;
    double  T = 1/(1-R_fresnel);
    check_close(a.plane_irradiance_down()[0], T, 0.02_pct);
  } end_test_case()
}
