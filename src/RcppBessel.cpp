#include <Rcpp.h>
#include "bessel.h"
#include "conversion.h"


SEXP BesselK_wrapper_real(const Rcpp::NumericVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<double> vec_z = Rcpp::as<std::vector<double>>(z);
  std::vector<double> result = bessel::BesselK_real_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(result);
}

SEXP BesselK_wrapper_complex(const Rcpp::ComplexVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<std::complex<double>> vec_z = to_complex_vector(z);
  std::vector<std::complex<double>> result = bessel::BesselK_complex_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}
SEXP BesselI_wrapper_real(const Rcpp::NumericVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<double> vec_z = Rcpp::as<std::vector<double>>(z);
  std::vector<double> result = bessel::BesselI_real_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(result);
}

SEXP BesselI_wrapper_complex(const Rcpp::ComplexVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<std::complex<double>> vec_z = to_complex_vector(z);
  std::vector<std::complex<double>> result = bessel::BesselI_complex_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}

SEXP BesselJ_wrapper_real(const Rcpp::NumericVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<double> vec_z = Rcpp::as<std::vector<double>>(z);
  std::vector<double> result = bessel::BesselJ_real_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(result);
}

SEXP BesselJ_wrapper_complex(const Rcpp::ComplexVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<std::complex<double>> vec_z = to_complex_vector(z);
  std::vector<std::complex<double>> result = bessel::BesselJ_complex_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}

// Wrapper functions for BesselY
SEXP BesselY_wrapper_real(const Rcpp::NumericVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<double> vec_z = Rcpp::as<std::vector<double>>(z);
  std::vector<double> result = bessel::BesselY_real_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(result);
}

SEXP BesselY_wrapper_complex(const Rcpp::ComplexVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<std::complex<double>> vec_z = to_complex_vector(z);
  std::vector<std::complex<double>> result = bessel::BesselY_complex_std(vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}

SEXP BesselH_wrapper_real(int m, const Rcpp::NumericVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<double> vec_z = Rcpp::as<std::vector<double>>(z);
  std::vector<std::complex<double>> result = bessel::BesselH_real_std(m, vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}

SEXP BesselH_wrapper_complex(int m, const Rcpp::ComplexVector& z, double nu, bool expon_scaled, int verbose) {
  std::vector<std::complex<double>> vec_z = to_complex_vector(z);
  std::vector<std::complex<double>> result = bessel::BesselH_complex_std(m, vec_z, nu, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}


SEXP AiryA_wrapper_real(const Rcpp::NumericVector& z, int deriv, bool expon_scaled, int verbose) {
  std::vector<double> vec_z = Rcpp::as<std::vector<double>>(z);
  std::vector<double> result = bessel::AiryA_real_std(vec_z, deriv, expon_scaled, verbose);
  return Rcpp::wrap(result);
}

SEXP AiryA_wrapper_complex(const Rcpp::ComplexVector& z, int deriv, bool expon_scaled, int verbose) {
  std::vector<std::complex<double>> vec_z = to_complex_vector(z);
  std::vector<std::complex<double>> result = bessel::AiryA_complex_std(vec_z, deriv, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}

SEXP AiryB_wrapper_real(const Rcpp::NumericVector& z, int deriv, bool expon_scaled, int verbose) {
  std::vector<double> vec_z = Rcpp::as<std::vector<double>>(z);
  std::vector<double> result = bessel::AiryB_real_std(vec_z, deriv, expon_scaled, verbose);
  return Rcpp::wrap(result);
}

SEXP AiryB_wrapper_complex(const Rcpp::ComplexVector& z, int deriv, bool expon_scaled, int verbose) {
  std::vector<std::complex<double>> vec_z = to_complex_vector(z);
  std::vector<std::complex<double>> result = bessel::AiryB_complex_std(vec_z, deriv, expon_scaled, verbose);
  return Rcpp::wrap(from_complex_vector(result));
}

//[[Rcpp::interfaces(cpp,r)]]
//[[Rcpp::export(bessel_k)]]
SEXP BesselK(SEXP z, double nu, bool expon_scaled = false, int verbose = 0) {
  if (Rf_isNumeric(z)) {
    Rcpp::NumericVector rz(z);
    bool has_negative = false;
    for (int i = 0; i < rz.size(); ++i) {
      if (rz[i] < 0) {
        has_negative = true;
        break;
      }
    }

    if (has_negative) {
      Rcpp::ComplexVector cz(rz.size());
      for (int i = 0; i < rz.size(); ++i) {
        Rcomplex c;
        c.r = rz[i];
        c.i = 0.0;
        cz[i] = c;
      }
      std::vector<std::complex<double>> vec_cz = Rcpp::as<std::vector<std::complex<double>>>(cz);
      std::vector<std::complex<double>> result = bessel::BesselK_complex_std(vec_cz, nu, expon_scaled, verbose);
      return Rcpp::wrap(result);
    } else {
      std::vector<double> vec_rz = Rcpp::as<std::vector<double>>(rz);
      std::vector<double> result = bessel::BesselK_real_std(vec_rz, nu, expon_scaled, verbose);
      return Rcpp::wrap(result); // Directly wrap the numeric result
    }
  } else if (Rf_isComplex(z)) {
    Rcpp::ComplexVector cz(z);
    std::vector<std::complex<double>> vec_cz = Rcpp::as<std::vector<std::complex<double>>>(cz);
    std::vector<std::complex<double>> result = bessel::BesselK_complex_std(vec_cz, nu, expon_scaled, verbose);
    return Rcpp::wrap(result);
  } else {
    Rcpp::stop("Unsupported input type");
  }
}

// BesselI
//[[Rcpp::interfaces(cpp,r)]]
//[[Rcpp::export(bessel_i)]]
SEXP BesselI(SEXP z, double nu, bool expon_scaled = false, int verbose = 0) {
  if (Rf_isNumeric(z)) {
    return BesselI_wrapper_real(Rcpp::as<Rcpp::NumericVector>(z), nu, expon_scaled, verbose);
  } else if (Rf_isComplex(z)) {
    return BesselI_wrapper_complex(Rcpp::as<Rcpp::ComplexVector>(z), nu, expon_scaled, verbose);
  } else {
    Rcpp::stop("Unsupported input type");
  }
}

// BesselJ
//[[Rcpp::interfaces(cpp,r)]]
//[[Rcpp::export(bessel_j)]]
SEXP BesselJ(SEXP z, double nu, bool expon_scaled = false, int verbose = 0) {
  if (Rf_isNumeric(z)) {
    return BesselJ_wrapper_real(Rcpp::as<Rcpp::NumericVector>(z), nu, expon_scaled, verbose);
  } else if (Rf_isComplex(z)) {
    return BesselJ_wrapper_complex(Rcpp::as<Rcpp::ComplexVector>(z), nu, expon_scaled, verbose);
  } else {
    Rcpp::stop("Unsupported input type");
  }
}

// BesselY
//[[Rcpp::interfaces(cpp,r)]]
//[[Rcpp::export(bessel_y)]]
SEXP BesselY(SEXP z, double nu, bool expon_scaled = false, int verbose = 0) {
  if (Rf_isNumeric(z)) {
    Rcpp::NumericVector rz(z);
    bool has_negative = false;
    for (int i = 0; i < rz.size(); ++i) {
      if (rz[i] < 0) {
        has_negative = true;
        break;
      }
    }
    if (has_negative) {
      Rcpp::ComplexVector cz(rz.size());
      for (int i = 0; i < rz.size(); ++i) {
        Rcomplex c;
        c.r = rz[i];
        c.i = 0.0;
        cz[i] = c;
      }
      std::vector<std::complex<double>> vec_cz = Rcpp::as<std::vector<std::complex<double>>>(cz);
      std::vector<std::complex<double>> result = bessel::BesselY_complex_std(vec_cz, nu, expon_scaled, verbose);
      return Rcpp::wrap(result);
    } else {
      std::vector<double> vec_rz = Rcpp::as<std::vector<double>>(rz);
      std::vector<double> result = bessel::BesselY_real_std(vec_rz, nu, expon_scaled, verbose);
      return Rcpp::wrap(result); // Directly wrap the numeric result
    }
  } else if (Rf_isComplex(z)) {
    Rcpp::ComplexVector cz(z);
    std::vector<std::complex<double>> vec_cz = Rcpp::as<std::vector<std::complex<double>>>(cz);
    std::vector<std::complex<double>> result = bessel::BesselY_complex_std(vec_cz, nu, expon_scaled, verbose);
    return Rcpp::wrap(result);
  } else {
    Rcpp::stop("Unsupported input type");
  }
}

// Main function to handle logic for different cases
//[[Rcpp::interfaces(cpp,r)]]
//[[Rcpp::export(bessel_h)]]
SEXP BesselH(int m, SEXP z, double nu, bool expon_scaled = false, int verbose = 0) {
  // Ensure m is either 1 or 2
  if (m != 1 && m != 2) {
    Rcpp::stop("Invalid value for m. It should be either 1 or 2.");
  }
  // Check if the input is numeric or complex
  if (Rf_isNumeric(z)) {
    return BesselH_wrapper_real(m, Rcpp::as<Rcpp::NumericVector>(z), nu, expon_scaled, verbose);
  } else if (Rf_isComplex(z)) {
    return BesselH_wrapper_complex(m, Rcpp::as<Rcpp::ComplexVector>(z), nu, expon_scaled, verbose);
  } else {
    Rcpp::stop("Unsupported input type");
  }
}


//[[Rcpp::interfaces(cpp,r)]]
//[[Rcpp::export(airy_a)]]
SEXP AiryA(SEXP z, int deriv = 0, bool expon_scaled = false, int verbose = 0) {
  if (Rf_isNumeric(z)) {
    return AiryA_wrapper_real(Rcpp::as<Rcpp::NumericVector>(z), deriv, expon_scaled, verbose);
  } else if (Rf_isComplex(z)) {
    return AiryA_wrapper_complex(Rcpp::as<Rcpp::ComplexVector>(z), deriv, expon_scaled, verbose);
  } else {
    Rcpp::stop("Unsupported input type");
  }
}

//[[Rcpp::interfaces(cpp,r)]]
//[[Rcpp::export(airy_b)]]
SEXP AiryB(SEXP z, int deriv = 0, bool expon_scaled = false, int verbose = 0) {
  if (Rf_isNumeric(z)) {
    return AiryB_wrapper_real(Rcpp::as<Rcpp::NumericVector>(z), deriv, expon_scaled, verbose);
  } else if (Rf_isComplex(z)) {
    return AiryB_wrapper_complex(Rcpp::as<Rcpp::ComplexVector>(z), deriv, expon_scaled, verbose);
  } else {
    Rcpp::stop("Unsupported input type");
  }
}
