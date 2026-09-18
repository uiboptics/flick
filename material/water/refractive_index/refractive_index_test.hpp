#include "refractive_index.hpp"

namespace flick {
  begin_test_case(refractive_index_test) {
    material::water::refractive_index n(35,273+20);
    check_close(n.value(500e-9),1.34,1_pct);
    check_close(n.value(699.9e-9),n.value(700.1e-9),0.1_pct);
    n.set_S_T(0,273+20);
    check_close(n.value(500e-9),1.337,0.1_pct);
    check_close(n.value(2500e-9),1.26,0.15_pct);
  } end_test_case()
}
