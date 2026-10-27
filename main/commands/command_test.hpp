#include "radiator.hpp"

namespace flick {
  begin_test_case(command_test) {
    check(system("flick radiator planck 5800 1 > flick_test_tmp")==0);
    check(system("flick accurt -g flick_test_config > flick_test_tmp")==0);
    check(system("flick accurt flick_test_config > flick_test_tmp")==0);
    system("rm -f flick_test_tmp flick_test_config");
  } end_test_case()

  begin_test_case(command_accurt_parallel_test) {
    system("mkdir -p flick_test_parallel_tmp/parallel_000");
    system("mkdir -p flick_test_parallel_tmp/parallel_001");

    system("flick accurt -g test_config");
    system("flick text test_config set flick_tmp_directory_name flick_test_parallel_tmp/parallel_000 > test_config_tmp");
    system("mv test_config_tmp flick_test_parallel_tmp/parallel_000/test_config");

    system("flick accurt -g test_config");
    system("flick text test_config set flick_tmp_directory_name flick_test_parallel_tmp/parallel_001 > test_config_tmp");
    system("mv test_config_tmp flick_test_parallel_tmp/parallel_001/test_config");

    check(system("flick accurt --parallel flick_tmp > flick_test_tmp")==0);
    system("rm -fr flick_tmp");
    system("rm -fr flick_test_parallel_tmp");
    system("rm -f flick_test_tmp");
    system("rm -f flick_test_config");
    system("rm -f config");
    system("rm -f test_config");

  } end_test_case()
}

