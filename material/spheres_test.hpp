#include "spheres.hpp"
#include "water/pure_water.hpp"
#include <numbers>

namespace flick {
  double compare_rayleigh(size_t row, size_t col, double angle) {
    double r = 1e-10;
    double sigma = 0.1;
    material::bubbles_in_ice<monodispersed_mie> bi(1,stdvector{angle},log(r),sigma);
    bi.percentage_accuracy(0.01);
    auto m1 = bi.mueller_matrix(unit_vector{angle,0});
    auto m2 = rayleigh_mueller(angle,0);
    return abs(1-m1.value(row,col)/m2.value(row,col))*100;
  }
  begin_test_case(spheres_test_A) {
    material::vacuum v;
    material::pure_water pw;
    double r = 10e-6;
    log_normal_distribution sd(log(r),0.001);
    double f = 0.1;  
    material::spheres<log_normal_distribution, material::vacuum,
		      material::pure_water,
		      parameterized_monodispersed_mie> s(f,stdvector{0},sd,v,pw);
    double p = 0.1;
    s.percentage_accuracy(p);
    s.set_wavelength(400e-9);
    check_close(s.scattering_coefficient(),3./2*f/r,p);

    double pi = std::numbers::pi;
    p = 0.01; //%
    check_small(compare_rayleigh(0,0,0),p);
    check_small(compare_rayleigh(2,2,0),p);
    check_small(compare_rayleigh(3,3,0),p);
    check_small(compare_rayleigh(0,1,pi/2),p);
    check_small(compare_rayleigh(1,0,pi/2),p);
  } end_test_case()
  
  begin_test_case(spheres_test_B) {
    // Check size-distribution integration convergence for small particles
    material::vacuum v;
    material::pure_water pw;
    double r = 1e-10;
    double sigma = 0.1;
    double angle = std::numbers::pi/2;
    material::bubbles_in_ice<monodispersed_mie> bi(1,stdvector{angle},log(r),sigma);
    bi.percentage_accuracy(0.1);
    bi.phase_function_only(true);
    bi.mueller_matrix(unit_vector{angle,0});
    check(true);
  } end_test_case()
}
