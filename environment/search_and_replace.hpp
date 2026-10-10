#ifndef flick_search_and_replace
#define flick_search_and_replace

#include "input_output.hpp"

namespace flick {
  class parameter_text {
    std::string s_;
    std::string begin_qualifier_ = "/*";
    std::string end_qualifier_ = "*/";
  public:
    parameter_text() = default;
    parameter_text(const std::string& s) : s_{s} {}
    void set_begin_qualifier(const std::string s) {
      begin_qualifier_ = s;
      end_qualifier_ = s;
    }
    void set_text_qualifiers(const std::string& begin,
			     const std::string& end) {
      begin_qualifier_ = begin;
      end_qualifier_ = end;
    }
    std::string get(const std::string& parameter) {
      size_t n = find_parameter(parameter);
      if (n == std::string::npos)
	throw std::runtime_error("\""+parameter+" = \""+" not found in parameter text");
      size_t n_first = n + parameter.length() + 2; 
      size_t n_last = s_.find(begin_qualifier_, n_first);
      return trim(s_.substr(n_first, n_last-n_first));
    }
    void set(const std::string& parameter, const std::string& value) {
      size_t n = find_parameter(parameter);
      if (n == std::string::npos)
	throw std::runtime_error("\""+parameter+" = \""+" not found in parameter text");
      std::string old_str = parameter + " = " + get(parameter);
      std::string new_str = parameter + " = " + value;
      s_.replace(n,old_str.length(), new_str);
    }
  private:
    size_t find_parameter(const std::string& parameter) const {
      if (parameter.empty())
	return std::string::npos;

      const std::string assignment = parameter + " = ";
      size_t position = 0;
      while (position < s_.size()) {
	size_t parameter_position = s_.find(assignment, position);
	size_t description_position = s_.find(begin_qualifier_, position);

	if (parameter_position == std::string::npos)
	  return std::string::npos;
	if (description_position == std::string::npos or
	    parameter_position < description_position)
	  return parameter_position;

	size_t description_end = s_.find(
	  end_qualifier_, description_position + begin_qualifier_.length());
	if (description_end == std::string::npos)
	  return std::string::npos;
	position = description_end + end_qualifier_.length();
	if (begin_qualifier_ == end_qualifier_) {
	  while (s_.compare(position, end_qualifier_.length(),
			    end_qualifier_) == 0)
	    position += end_qualifier_.length();
	}
      }
      return std::string::npos;
    }

    std::string trim(std::string s) const {
      const char* t = " \t\n\r\f\v";
      s.erase(0, s.find_first_not_of(t));
      s.erase(s.find_last_not_of(t) + 1);
      return s;
    }
    friend std::ostream& operator<<(std::ostream &os, const parameter_text& t) {
      os << t.s_;
      return os;
    }
    friend std::istream& operator>>(std::istream &is, parameter_text& t) {
      t.s_ = std::string(std::istreambuf_iterator<char>(is), {});
      return is;
    }
  };
}

#endif
