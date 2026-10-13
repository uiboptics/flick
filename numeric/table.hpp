#ifndef flick_table
#define flick_table
 
#include "function.hpp"

namespace flick {
  template<class Interpolation>
  class table
  {
    std::string header_;
    std::vector<double> row_values_;
    std::vector<double> col_values_;
    std::vector<function<Interpolation>> column_functions_;
  public:
    table() = default;
    table(const std::vector<double>& row_values,
	  const std::vector<double>& col_values,
	  const std::vector<std::vector<double>>& table_values)
      : row_values_{row_values}, col_values_{col_values} {
      column_functions_.resize(col_values_.size());
      for (size_t i=0; i<row_values_.size(); ++i) {
	for (size_t j=0; j<col_values_.size(); ++j) {
	  column_functions_[j].append({row_values_[i],table_values[i][j]});	  
	}
      }
      ensure(!column_functions_.empty());
      ensure(col_values_.size()==column_functions_.size());
      ensure(row_values_.size()==column_functions_.at(0).size());
    }
    std::string header() const {
      return header_;
    }
    const std::vector<double>& row_values() {
      return row_values_;
    }
    const std::vector<double>& column_values() {
      return col_values_;
    }
    function<Interpolation> row(size_t n) {
      function<Interpolation> f;
      for (size_t i=0; i < column_functions_.size(); ++i)
	f.append({col_values_[i],column_functions_[i].y().at(n)});
      return f;
    }
    function<Interpolation> column(size_t n) {
      return function<Interpolation>{column_functions_.at(n).x(),
	column_functions_.at(n).y()};
    }
    double value(double row_value, double column_value) const {
      function<Interpolation> f;
      for (size_t i=0; i < col_values_.size(); ++i)
	f.append({col_values_[i], column_functions_[i].value(row_value)});
      return f.value(column_value);
    }
    double value_at(size_t row, size_t column) {
      return column_functions_.at(column).y().at(row);
    }
    double integral(double row_low, double row_high,
		    double column_low, double column_high) const {
      function<Interpolation> f;
      for (size_t i=0; i < col_values_.size(); ++i)
	f.append({col_values_[i], column_functions_[i].integral(row_low, row_high)});
      return f.integral(column_low, column_high);
    }
    auto& set(size_t row, size_t column, double value) {
      column_functions_.at(column).set_y(row, value);
      return *this;
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const table<Interpolation>& f) {
      os << f.header();
      os << f.row_values_.size() << std::endl;
      os << f.col_values_.size() << std::endl;
      for (size_t i=0; i<f.row_values_.size(); ++i) {
	os << f.row_values_[i] << " ";
      }
      os << std::endl;
      for (size_t i=0; i<f.col_values_.size(); ++i) {
	os << f.col_values_[i] << " ";
      }
      os << std::endl;
      for (size_t i=0; i<f.row_values_.size(); ++i) {
	for (size_t j=0; j<f.col_values_.size(); ++j) {
	  os << f.column_functions_[j].value(f.row_values_[i]) << " ";	  
	}
	os << std::endl;
      }
      return os;
    }
    friend std::istream& operator>>(std::istream &is,
				    table<Interpolation>& f) {
      f.header_ = read_header(is);
      size_t n_rows, n_cols;
      is >> n_rows >> n_cols;
      f.row_values_.resize(n_rows);
      for (size_t i=0; i<n_rows; ++i) {
	is >> f.row_values_[i];
      }
      f.col_values_.resize(n_cols);
      for (size_t i=0; i<n_cols; ++i) {
	is >> f.col_values_[i];
      }
      f.column_functions_.resize(n_cols);
      double value;
      for (size_t i=0; i<n_rows; ++i) {
	for (size_t j=0; j<n_cols; ++j) {
	  is >> value;
	  f.column_functions_[j].append({f.row_values_[i],value});	  
	}
      }
      return is;
    }
    void ensure(bool b) const {
      if (!b)
	throw std::runtime_error("numeric table");
    }
  };

  using pl_table = table<piecewise_linear>;
  using pp_table = table<piecewise_power>;
  using pe_table = table<piecewise_exponential>;
}

#endif
