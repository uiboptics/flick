#include "../../environment/unit_test.hpp"
#include "monocrome_tabulated_test.hpp"

int main() {
  using namespace flick;
  unit_test t("material");
  t.include<monocrome_tabulated_test>();
 
  t.run_test_cases();
  return 0;
}
