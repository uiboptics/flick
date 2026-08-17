#include "../../environment/unit_test.hpp"
#include "tabulated_iops_test.hpp"

int main() {
  using namespace flick;
  unit_test t("tabulated_iops");
  t.include<tabulated_iops_test>();
 
  t.run_test_cases();
  return 0;
}
