#include "radiator.hpp"

namespace flick {
  begin_test_case(command_test) {
    check(system("flick radiator planck 5800 1")==0);
  } end_test_case()
}

