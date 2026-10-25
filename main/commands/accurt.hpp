#ifndef flick_command_accurt
#define flick_command_accurt

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "basic_command.hpp"
#include "../../material/atmosphere_ocean.hpp"
#include "../../environment/input_output.hpp"
#include "../../environment/configuration.hpp"
#include "../../accurt_api/accurt.hpp"

namespace flick {
  namespace command {
    class accurt : public basic_command {
    public:
      accurt() : basic_command("accurt") {};
      void run() {
	if (size() == 3 && a(1) == "-g" ) {
	  configuration_template::atmosphere_ocean().write(a(2));
	} else if (size() == 4 && a(1) == "-g" ) {
	  if (a(2) == "atmosphere") {
	    configuration_template::atmosphere().write(a(3));
	  }
	  else if (a(2) == "boa_transmittance") {
	    configuration_template::boa_transmittance().write(a(3));
	  }
	  else if (a(2) == "atmosphere_ocean") {
	    configuration_template::atmosphere_ocean().write(a(3));
	  }
	  else if (a(2) == "default") {
	    configuration_template::atmosphere_ocean().write(a(3));
	  }
	  else if (a(2) == "toa_reflectance") {
	    configuration_template::toa_reflectance().write(a(3));
	  }
	  else if (a(2) == "ocean_radiance") {
	    configuration_template::ocean_radiance().write(a(3));
	  }
	  else if (a(2) == "rs_reflectance") {
	    configuration_template::rs_reflectance().write(a(3));
	  }
	  else if (a(2) == "ramses_reflectance") {
	    configuration_template::ramses_reflectance().write(a(3));
	  }
	} else if (size() == 2 && opt(0) == "--parallel") {
	  namespace fs = std::filesystem;
	  std::vector<fs::path> configs;
	  for (const auto& entry : fs::directory_iterator(a(1))) {
	    const auto name = entry.path().filename().string();
	    if (entry.is_directory() && name.starts_with("parallel_")) {
	      std::vector<fs::path> files;
	      for (const auto& file : fs::directory_iterator(entry.path())) {
		if (file.is_regular_file()) {
		  files.push_back(file.path());
		}
	      }
	      if (files.size() > 1) {
		throw std::runtime_error(
		  "accurt command: more than one file found in '" +
		  entry.path().string() + "'");
	      }
	      if (files.size() == 1) {
		configs.push_back(files.front());
	      }
	    }
	  }
	  std::sort(configs.begin(), configs.end());
	  std::vector<pid_t> children;
	  std::vector<fs::path> output_files;
	  std::vector<fs::path> error_files;
	  for (size_t i=0; i<configs.size(); ++i) {
	    const auto id = std::to_string(static_cast<long long>(getpid())) +
	      "-" + std::to_string(i);
	    const auto output_file = fs::temp_directory_path() /
	      ("flick-accurt-" + id + ".out");
	    const auto error_file = fs::temp_directory_path() /
	      ("flick-accurt-" + id + ".err");
	    const pid_t child = fork();
	    if (child == 0) {
	      std::ofstream output(output_file);
	      std::cout.rdbuf(output.rdbuf());
	      try {
		std::cout << stream_radiation(configs.at(i).string());
		std::cout.flush();
		_exit(EXIT_SUCCESS);
	      }
	      catch (const std::exception& e) {
		std::ofstream error(error_file);
		error << e.what();
		error.close();
		_exit(EXIT_FAILURE);
	      }
	      catch (...) {
		std::ofstream error(error_file);
		error << "unknown child-process error";
		error.close();
		_exit(EXIT_FAILURE);
	      }
	    }
	    if (child < 0) {
	      for (const auto running_child : children) {
		waitpid(running_child, nullptr, 0);
	      }
	      throw std::runtime_error("accurt command: unable to create child process");
	    }
	    children.push_back(child);
	    output_files.push_back(output_file);
	    error_files.push_back(error_file);
	  }
	  std::string child_error;
	  for (size_t i=0; i<children.size(); ++i) {
	    int status = 0;
	    if (waitpid(children.at(i), &status, 0) < 0 ||
		!WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS) {
	      std::ifstream error(error_files.at(i));
	      std::ostringstream message;
	      message << error.rdbuf();
	      child_error += "accurt command: process for '" +
		configs.at(i).string() + "' failed";
	      if (!message.str().empty()) child_error += ": " + message.str();
	      child_error += "\n";
	    }
	  }
	  if (child_error.empty()) {
	    for (const auto& output_file : output_files) {
	      std::ifstream output(output_file);
	      std::cout << output.rdbuf();
	    }
	  }
	  for (const auto& output_file : output_files) fs::remove(output_file);
	  for (const auto& error_file : error_files) fs::remove(error_file);
	  if (!child_error.empty()) {
	    throw std::runtime_error(child_error);
	  }
	} else if (size() == 2 ) {
	  std::cout << stream_radiation(a(1));
	}
      }
    private:
	std::string stream_radiation(const std::string config) {
	auto c0 = read<configuration_template::basic_accurt>(config);
	std::string name = c0.get<std::string>("material_name");
	std::shared_ptr<material::base> m;
	if (name == "atmosphere") {
	  auto c = read<configuration_template::atmosphere>(config);
	  m = std::make_shared<material::atmosphere>(c);
	}
	else if(name == "atmosphere_ocean") {
	  auto c = read<configuration_template::atmosphere_ocean>(config);
	  m = std::make_shared<material::atmosphere_ocean>(c);
	}
	else {
	  throw std::runtime_error("accurt command: material '"+name+"' not found");
	}
	flick::accurt a(c0,m);
	std::ostringstream output;
	if (a.radiance_distribution_override() and not
	    a.conical_radiation_override()) {
	  const auto r = a.relative_distributed_radiance().first;
	  for (const auto& polar : r) {
	    for (const auto& azimuth : polar) {
	      for (const auto value : azimuth) {
		output << value << " ";
	      }
	      output << '\n';
	    }
	  }
	} else
	  output << a.relative_radiation();
	output << '\n';
	return output.str();
      }
    };    
  }
}

#endif
