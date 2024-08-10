#ifndef bessel_h
#define bessel_h

#include <Rcpp.h>
#include <complex>
#include <vector>
#include <limits>
#include <iostream>
#include <cmath>
extern "C" {
#include "zbsubs.h"
}

// Function to calculate BesselK for std::vector<std::complex<double>>
//
namespace bessel {
  std::vector<std::complex<double > > BesselK_complex_std(const std::vector<std::complex<double > >& z, double nu, bool expon_scaled, int verbose);
  std::vector<double> BesselK_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose);

  std::vector<double> BesselI_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose);
  std::vector<std::complex<double>> BesselI_complex_std(const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose);

  std::vector<double> BesselJ_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose);
  std::vector<std::complex<double>> BesselJ_complex_std(const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose);

  std::vector<double> BesselY_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose);
  std::vector<std::complex<double>> BesselY_complex_std(const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose);

  std::vector<std::complex<double>> BesselH_real_std(int m, const std::vector<double>& z, double nu, bool expon_scaled, int verbose);
  std::vector<std::complex<double>> BesselH_complex_std(int m, const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose);

  std::vector<double> AiryA_real_std(const std::vector<double>& z, int deriv, bool expon_scaled, int verbose);
  std::vector<std::complex<double>> AiryA_complex_std(const std::vector<std::complex<double>>& z, int deriv, bool expon_scaled, int verbose);

  std::vector<double> AiryB_real_std(const std::vector<double>& z, int deriv, bool expon_scaled, int verbose);
  std::vector<std::complex<double>> AiryB_complex_std(const std::vector<std::complex<double>>& z, int deriv, bool expon_scaled, int verbose);

}


#endif //
