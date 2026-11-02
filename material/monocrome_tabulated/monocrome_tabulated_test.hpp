#include "monocrome_tabulated.hpp"

namespace flick {
  begin_test_case(monocrome_tabulated_test) {
    material::monocrome_tabulated m(0, 1, "vsfPart1_bp0.1");
    double g = m.asymmetry_factor();
    tabulated_phase_function p = make_tabulated_phase_function(m);
    check_close(p.integral_4pi(), 1, 0.04_pct);
    check_close(p.asymmetry_factor(), g, 0.003_pct);
  } end_test_case()
}
