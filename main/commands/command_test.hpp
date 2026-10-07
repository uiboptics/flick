#include "radiator.hpp"

namespace flick {
  begin_test_case(command_test) {
    check(system("flick radiator planck 5800 1 > flick_test_tmp")==0);
    check(system("flick accurt -g flick_test_config > flick_test_tmp")==0);
    check(system("flick accurt flick_test_config > flick_test_tmp")==0);
  } end_test_case()
}

