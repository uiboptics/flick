#ifndef flick_std_operators
#define flick_std_operators
#include <vector>
#include <complex>
#include <iostream>

namespace flick {
  using stdvector = std::vector<double>;
  using stdcomplex = std::complex<double>;
  using stdvectorc = std::vector<stdcomplex>;
  using namespace std::complex_literals;
  
  stdvector& operator+=(stdvector& v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] += v2[i];
    return v1;
  }
  stdvector operator+(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] += v2[i];
    return v1;
  }
  stdvectorc operator+(stdvectorc v1, const stdvectorc& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] += v2[i];
    return v1;
  }
  stdvector operator+(stdvector v, double k) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] += k;
    return v;
  }
  stdvector operator+(double k, stdvector v) {
    return v+k;
  }
  stdvector operator-(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] -= v2[i];
    return v1;
  }
  stdvectorc operator-(stdvectorc v1, const stdvectorc& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] -= v2[i];
    return v1;
  }
  stdvector operator-(stdvector v, double k) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] -= k;
    return v;
  }
  stdvector operator-(double k, stdvector v) {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = k-v[i];
    return v;
  }
   stdvector& operator*=(stdvector& v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] *= v2[i];
    return v1;
  }
  stdvector operator*(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] *= v2[i];
    return v1;
  }
  stdvectorc operator*(stdvectorc v1, const stdvectorc& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] *= v2[i];
    return v1;
  }
  stdvector operator*(stdvector v, double k) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] *= k;
    return v;
  }
  stdvector operator*(double k, stdvector v) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] *= k;
    return v;
  }
  stdvectorc operator*(stdcomplex k, stdvectorc v) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] *= k;
    return v;
  }
  stdvector operator/(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] /= v2[i];
    return v1;
  }
  stdvectorc operator/(stdvectorc v1, const stdvectorc& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] /= v2[i];
    return v1;
  }
  stdvector operator/(stdvector v, double k) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] /= k;
    return v;
  }
  stdvector operator/(double k, stdvector v) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] = k/v[i];
    return v;
  }
  stdvector operator^(stdvector v, double k) {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = std::pow(v[i],k);
    return v;
  }
  namespace vec {
    double sum(const stdvector& v) {
      double s=0;
      for (size_t i=0; i<v.size(); ++i)
	s += v[i];
      return s;
    }
    stdvector sin(stdvector v) {
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::sin(v[i]);
      return v;
    }
    stdvector cos(stdvector v) {
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::cos(v[i]);
      return v;
    }
    stdvector acos(stdvector v) {
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::acos(v[i]);
      return v;
    }
    stdvector exp(stdvector v) {
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::exp(v[i]);
      return v;
    }
    stdcomplex sum(const stdvectorc& v) {
      stdcomplex s{0,0};
      for (size_t i=0; i<v.size(); ++i)
	s += v[i];
      return s;
    }
    double mean(const stdvector& v) {
      return sum(v)/v.size();
    }
    double rms(const stdvector& v) {
      return sqrt(mean(v^2));
    }
    stdvector abs(stdvector v) {
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::abs(v[i]);
      return v;
    }
    stdvector abs(const stdvectorc& vc) {
      stdvector v(vc.size());
      for (size_t i=0; i<vc.size(); ++i)
	v[i] = std::abs(vc[i]);
      return v;
    }
    stdvectorc conj(stdvectorc v) {
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::conj(v[i]);
      return v;
    }
    stdvector real(const stdvectorc& vc) {
      stdvector v(vc.size());
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::real(vc[i]);
      return v;
    }
    stdvector imag(const stdvectorc& vc) {
      stdvector v(vc.size());
      for (size_t i=0; i<v.size(); ++i)
	v[i] = std::imag(vc[i]);
      return v;
    }
  }
  std::ostream& operator<<(std::ostream& out, const stdvector& v) {
    for (const auto& i: v)
      out << i << " ";
    return out;
  }
  std::ostream& operator<<(std::ostream& out, const stdvectorc& v) {
    for (const auto& i: v)
      out << i << " ";
    return out;
  }
}

#endif
