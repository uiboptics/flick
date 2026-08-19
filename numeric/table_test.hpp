#include "../environment/input_output.hpp"
#include "table.hpp"

namespace flick {
  begin_test_case(table_test) {  
    pl_table t = read<pl_table>("numeric/table_test_data.txt");
    check_close(t.value(0,37.376e-12),0.7419,0.1_pct);
    check(t.value(0.97,1294.0e-12) < 0.7790 and t.value(0.97,1294.0e-12) > 0.7692);
  } end_test_case()
}
