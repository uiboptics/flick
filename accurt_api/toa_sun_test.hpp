#include "../environment/unit_test.hpp"
#include "toa_sun.hpp"

namespace flick {
  begin_test_case(toa_sun_test_A) {
    /* Sun at zenith */
    std::string arg = "2026 3 20 12 0 0.0 -0.461033 2.041024 0";
    double ETR = 1373;
    toa_sun ts(arg);
    check_close(ts.spectrum().integral(200e-9,7000e-9),ETR,0.3_pct);

    arg.pop_back();
    arg += "10e-9";
    toa_sun ts_smoothed(arg);
    check_close(ts_smoothed.spectrum().integral(200e-9,7000e-9),ETR,0.3_pct);

  } end_test_case()
  
  begin_test_case(toa_sun_test_B) {
    /* October morning in Bergen */
    std::string arg = "2026 10 4 6 17 0.0 60.391 5.322 0";
    double ETR = 60;
    toa_sun ts(arg);
    
    check_close(ts.spectrum().integral(200e-9,7000e-9),ETR,0.3_pct);
  } end_test_case()
}
