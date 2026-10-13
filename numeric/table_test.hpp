#include "../environment/input_output.hpp"
#include "table.hpp"

namespace flick {
  begin_test_case(table_test_A) {  
    pl_table t = read<pl_table>("numeric/table_test_data.txt");
    check_close(t.value(0,37.376e-12),0.7419,0.1_pct);
    check(t.value(0.97,1294.0e-12) < 0.7790 and t.value(0.97,1294.0e-12) > 0.7692);
  } end_test_case()
  
  begin_test_case(table_test_B) {  
    pl_table t({-1,1},{-3,3},{{-7,5},
			      {-5,7}});
    check_close(t.value(0.1,1), 2.1, 1e-6_pct);
    check_small(t.integral(-1,1,-1,1), 1e-20);
    check_close(t.integral(-2,1,-1,1), -3, 1e-6_pct);
    check(t.set(0,0,-8).integral(-2,1,-1,1) < -5);
    
  } end_test_case()
}
