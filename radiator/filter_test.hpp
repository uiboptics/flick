#include "filter.hpp"
#include "cie_d65.hpp"
#include "cie_a.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  begin_test_case(filter_test_A) {
    double p = 0.1;
    check_close(filter::erythema().transmittance(298e-9),1,p);
    check_close(filter::erythema().transmittance(310e-9),0.07447,p);
    check_close(filter::erythema().transmittance(350e-9),0.000707946,p);
    pl_function f = read<pl_function>("./toa_solar.txt");
    auto f2 = filter::cut_ends(280e-9,400e-9).transmission(f);
    check_close(uv_index(f2),uva_index(f2)+uvb_index(f2),0.3_pct);
    double wl0 = 0.5e-6;
    double fwhm = 5e-9;
    check_close(filter::gaussian(wl0,fwhm).weighted_average(f),
		filter::square(wl0,fwhm).weighted_average(f),1.0_pct);
  } end_test_case()
  
  begin_test_case(filter_test_B) {
    check(filter::sentinel3().closest_srf(300e-9)==0);
    check(filter::sentinel3().transmittance(400e-9) > 1e-5);
    check_small(filter::sentinel3().transmittance(300e-9));
    check_small(filter::sentinel3().transmittance(1500e-9));
    check(filter::sentinel3().closest_srf(1090e-9)==20);
    check(filter::sentinel3().closest_srf(938e-9)==19);
    check_close(filter::sentinel3().importance_sampled_wavelengths(1)[1], 412e-9,1_pct);

    filter::sentinel3 olci;
    olci.center_wavelength(762e-9);
    check(olci.closest_srf(762e-9)==12);
    check_close(olci.fwhm(),2.37e-9,1_pct);
    check(olci.transmittance(762e-9) > 1e7);
    check_small(olci.transmittance(754e-9));
 } end_test_case()
    
  begin_test_case(filter_test_C) {
    /* https://en.wikipedia.org/wiki/Standard_illuminant */
    size_t n = 500;
    pp_function E(range(380e-9,800e-9,n).linspace(),std::vector<double>(n,1));
    for (size_t i = 0; i<3; i++) {
      check_close(chromaticity(E)[i], 1./3, 0.03_pct);
    }
    pp_function A = radiator::cie_a().spectrum();
    check_close(chromaticity(A)[0], 0.44758, 0.03_pct);
    check_close(chromaticity(A)[1], 0.40745, 0.03_pct);

    pp_function D65 = radiator::cie_d65().spectrum();
    check_close(chromaticity(D65)[0], 0.31272, 0.03_pct);
    check_close(chromaticity(D65)[1], 0.32903, 0.03_pct);
    for (size_t i = 0; i<3; i++) {
      check_close(rgb(D65)[i], 1, 0.03_pct);
    }
  } end_test_case()
 
  begin_test_case(filter_test_D) {
    pl_function f = read<pl_function>("./toa_solar.txt");
    check_close(n_photons(f,400e-9,700e-9),1.4e21,5_pct);
  } end_test_case()
  
  begin_test_case(filter_test_E) {
    auto f = read<pl_function>("./toa_solar.txt");
    auto f2 = read<pp_function>("./toa_solar.txt");
    double wl_lo = 250e-9;
    double wl_hi = 4000e-9;
    double wl_0 = 400e-9;
    double p = 0.003_pct;
    double fwhm = 10e-9;
    double a = f.integral(wl_lo, wl_hi);
    double b1 = moving_average(f,filter::gaussian(wl_0,fwhm)).integral(wl_lo,wl_hi);
    double b2 = moving_average(f,filter::triangular(wl_0,fwhm)).integral(wl_lo,wl_hi);
    filter::square filter(fwhm);
    double b3 = moving_average(f,filter).integral(wl_lo,wl_hi);
    double b4 = moving_average(f2,filter::gaussian(fwhm)).integral();
    check_close(a,b1,p);
    check_close(a,b2,p);
    check_close(a,b3,p);
    check_close(f2.integral(),b4,p);
  } end_test_case()
}
