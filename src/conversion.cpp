#include "conversion.h"

std::vector<std::complex<double>> to_complex_vector(const Rcpp::ComplexVector& z) {
  std::vector<std::complex<double>> vec_z;
  vec_z.reserve(z.size());
  for (int i = 0; i < z.size(); i++) {
    vec_z.push_back(std::complex<double>(z[i].r, z[i].i));
  }
  return vec_z;
}

Rcpp::ComplexVector from_complex_vector(const std::vector<std::complex<double>>& vec) {
  Rcpp::ComplexVector res(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    res[i].r = vec[i].real();
    res[i].i = vec[i].imag();
  }
  return res;
}
