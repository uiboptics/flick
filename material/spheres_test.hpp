#include "spheres.hpp"
#include "water/pure_water.hpp"
#include <numbers>

namespace flick {
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
    material::bubbles_in_ice<monodispersed_mie> bi(1,stdvector{0,pi/2},log(1e-10),0.0001);
    p = 0.5; //%
    bi.percentage_accuracy(p);
    auto m1 = bi.mueller_matrix(unit_vector{0,0});
    auto m2 = rayleigh_mueller(0,0);
    check_close(m1.value(0,0),m2.value(0,0),p);
    check_close(m1.value(2,2),m2.value(2,2),p);
    check_close(m1.value(3,3),m2.value(3,3),p);

    p = 1;
    bi.percentage_accuracy(p);
    double theta = pi/2;
    m1 = bi.mueller_matrix(unit_vector{theta,0});
    m2 = rayleigh_mueller(theta,0);
    check_close(m1.value(0,0),m2.value(0,0),p);
    check_close(m1.value(1,1),m2.value(1,1),p);
    //check_close(m1.value(0,1),m2.value(0,1),p);
    //check_close(m1.value(1,0),m2.value(1,0),p);
    
  } end_test_case()

}
