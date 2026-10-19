#include "accurt_kernel_wrapper_test.hpp"

int main() {
  using namespace flick;
  unit_test t("accurt kernel wrapper");
  
  t.include<accurt_kernel_wrapper_test_A>();
  t.include<accurt_test_B>();
  t.include<accurt_test_C>();
  
  t.run_test_cases();
  return 0;
} 
