#ifndef flick_input_output
#define flick_input_output

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <filesystem>

namespace flick {
  std::string path() {
    std::string path;
    char *p = getenv("FLICK_PATH");
    if (p==NULL)
      throw std::invalid_argument("FLICK_PATH not found. See Prerequisites");
    return std::string(p);
  }

  class columns {
    size_t n_cols_;
    std::vector<std::vector<double>> cols_;
  public:
    columns(size_t n_cols) : n_cols_{n_cols} {
      cols_.resize(n_cols_);
    }
    size_t size() const {
      return cols_.at(0).size();
    }
    const std::vector<double>& column(size_t n) const {
      return cols_.at(n);
    }
    friend std::istream& operator>>(std::istream &is, columns& c) {
      bool stop = false;
      while (!stop) {
	for (size_t i=0; i<c.n_cols_; ++i) {
	  double x;
	  if(not(is >> x)) {
	    stop = true;
	  } else {
	  c.cols_[i].push_back(x);
	  }
	}
      }
      return is;
    }
    friend std::ostream& operator<<(std::ostream &os, columns& c) {
      for (size_t i=0; i<c.column(0).size(); ++i) {
	for (size_t j=0; j<c.n_cols_; ++j) {
	  os << c.cols_[j][i] << " ";
	}
	os << "\n";
      }
      return os;
    }
  };
  
  struct two_columns : public columns {
    two_columns() : columns::columns(2) {}
  };
    
  class text {
    std::string t_;
    friend std::istream& operator>>(std::istream &is, text& t) {
      t.t_= std::string(std::istreambuf_iterator<char>(is),
	std::istreambuf_iterator<char>());
      return is;
    }
    friend std::ostream& operator<<(std::ostream &os,text& t) {
      os << t.t_;
      return os;
    }
  };
  
  template<typename T>
  T read(std::string file, const std::string& alternative_path="") {
    if (not std::filesystem::exists(file))
      file = alternative_path+"/"+file;
    if (not std::filesystem::exists(file))
      file = path()+file;
    std::ifstream ifs(file);
    if (!ifs)
      throw std::invalid_argument(file+" not found");
    T t;
    ifs >> t;
    ifs.close();
    return t;	
  }  

  template<typename T>
  void write(const T& t, const std::string& file, size_t precision=5) {
    std::ofstream ofs(file);
    if (!ofs)
      throw std::invalid_argument(file+" could not be opened");
    
    ofs << std::setprecision(precision) << t;
    ofs.close();
  }
}

#endif
  
