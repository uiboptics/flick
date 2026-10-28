#include "radiator.hpp"

#include <cstdlib>
#include <filesystem>

namespace flick {
  begin_test_case(command_test) {
    check(system("flick radiator planck 5800 1 > flick_test_tmp")==0);
    check(system("flick accurt -g flick_test_config > flick_test_tmp")==0);
    check(system("flick accurt flick_test_config > flick_test_tmp")==0);
    namespace fs = std::filesystem;
    fs::remove("flick_test_tmp");
    fs::remove("flick_test_config");
  } end_test_case()
  
  begin_test_case(command_accurt_parallel_test) {
    namespace fs = std::filesystem;

    fs::create_directories(
			   "flick_test_parallel_tmp/parallel_000"
			   );
    fs::create_directories(
			   "flick_test_parallel_tmp/parallel_001"
			   );

    check(std::system("flick accurt -g test_config") == 0);
    check(std::system(
		      "flick text test_config set flick_tmp_directory_name "
		      "flick_test_parallel_tmp/parallel_000 > test_config_tmp"
		      ) == 0);
    
    fs::rename(
	       "test_config_tmp",
	       "flick_test_parallel_tmp/parallel_000/test_config"
	       );
    
    check(std::system("flick accurt -g test_config") == 0);
    check(std::system(
		      "flick text test_config set flick_tmp_directory_name "
		      "flick_test_parallel_tmp/parallel_001 > test_config_tmp"
		      ) == 0);
    
    fs::rename(
	       "test_config_tmp",
	       "flick_test_parallel_tmp/parallel_001/test_config"
	       );
    
    check(std::system(
		      "flick accurt --parallel flick_tmp > flick_test_tmp"
		      ) == 0);

    fs::remove_all("flick_tmp");
    fs::remove_all("flick_test_parallel_tmp");
    fs::remove("flick_test_tmp");
    fs::remove("flick_test_config");
    fs::remove("config");
    fs::remove("test_config");
    
  } end_test_case()  
}


