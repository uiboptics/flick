#include "../../../environment/unit_test.hpp"
#include "refractive_index_test.hpp"

int main() {
  using namespace flick;
  unit_test t("water");
  t.include<refractive_index_test>();
  t.run_test_cases();
  return 0;
}
