#include "../../environment/unit_test.hpp"
#include "snow_impurity_test.hpp"

int main() {
  using namespace flick;
  unit_test t("snow impurity");
  t.include<snow_impurity_test>();
 
  t.run_test_cases();
  return 0;
}
