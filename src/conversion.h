#include <Rcpp.h>
#include <vector>
#include <complex>

std::vector<std::complex<double>> to_complex_vector(const Rcpp::ComplexVector& z);
Rcpp::ComplexVector from_complex_vector(const std::vector<std::complex<double>>& vec);
