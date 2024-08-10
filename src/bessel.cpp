#include "bessel.h"

namespace bessel {
  std::vector<std::complex<double>> BesselK_complex_std(const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose) {
    if (nu < 0) nu = -1.0 * nu;
    int n = z.size();
    int nSeq = 1;
    std::vector<std::complex<double>> result(n);
    for (int i = 0; i < n; i++) {
      double zr = std::real(z[i]);
      double zi = std::imag(z[i]);
      int ierr = verbose;
      int kode = expon_scaled ? 2 : 1;
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz;

      // Call to zbesk function from the included C code
      zbesk(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);

      if (ierr != 0) {
        std::string f_x = "'zbesk(" + std::to_string(zr) + (zi >= 0 ? " + " : " - ") + std::to_string(std::abs(zi)) + "i, nu=" + std::to_string(nu) + ")'";
        if (ierr == 3) {
          Rcpp::warning("%s large arguments -> precision loss (of at least half machine accuracy)", f_x);
        } else if (ierr == 2) {
          if (verbose) Rcpp::Rcout << f_x << "  -> overflow ; returning Inf\n";
          std::fill(cyr.begin(), cyr.end(), INFINITY);
          std::fill(cyi.begin(), cyi.end(), INFINITY);
        } else if (ierr == 4) {
          Rcpp::warning("%s  -> ierr=4: |z| or nu too large\n", f_x);
          std::fill(cyr.begin(), cyr.end(), NAN);
          std::fill(cyi.begin(), cyi.end(), NAN);
        } else {
          Rcpp::stop("%s unexpected error 'ierr = %d'", f_x, ierr);
        }
      }

      result[i] = std::complex<double>(cyr[0], cyi[0]);

      // Handle zero input after processing
      if (zr == 0.0 && zi == 0.0) {
        result[i] = std::complex<double>(std::numeric_limits<double>::infinity(), 0.0);
      }
    }
    return result;
  }


  // Function to calculate BesselK for std::vector<double>
  std::vector<double> BesselK_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose) {
    if (nu < 0) nu = -1.0 * nu;
    int nSeq = 1;
    int n = z.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
      double zr = z[i];
      double zi = 0.0;
      int ierr = verbose;
      int kode = expon_scaled ? 2 : 1;
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz;

      // Call to zbesk function from the included C code
      zbesk(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);

      if (ierr != 0) {
        std::string f_x = "'zbesk(" + std::to_string(zr) + (zi >= 0 ? " + " : " - ") + std::to_string(std::abs(zi)) + "i, nu=" + std::to_string(nu) + ")'";
        if (ierr == 3) {
          Rcpp::warning("%s large arguments -> precision loss (of at least half machine accuracy)", f_x);
        } else if (ierr == 2) {
          if (verbose) Rcpp::Rcout << f_x << "  -> overflow ; returning Inf\n";
          std::fill(cyr.begin(), cyr.end(), INFINITY);
          std::fill(cyi.begin(), cyi.end(), INFINITY);
        } else if (ierr == 4) {
          Rcpp::warning("%s  -> ierr=4: |z| or nu too large\n", f_x);
          std::fill(cyr.begin(), cyr.end(), NAN);
          std::fill(cyi.begin(), cyi.end(), NAN);
        } else {
          Rcpp::stop("%s unexpected error 'ierr = %d'", f_x, ierr);
        }
      }
      result[i] = cyr[0]; // Imaginary part should be zero
      if (zr == 0.0) {
        result[i] =std::numeric_limits<double>::infinity();
      }
    }
    return result;
  }

  // Function for real Bessel I
  std::vector<double> BesselI_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<double> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1;

    if (nu < 0) {
      if (nu == round(nu)) {
        return BesselI_real_std(z, -nu, expon_scaled, verbose);
      }
      std::vector<double> kf = BesselK_real_std(z, -nu, expon_scaled, verbose);
      if (expon_scaled) {
        for (size_t i = 0; i < z.size(); ++i) {
          kf[i] *= exp(-z[i] - std::abs(z[i]));
        }
      }
      std::vector<double> bi = BesselI_real_std(z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bi[i] + 2 / M_PI * sin(M_PI * nu) * kf[i];
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i];
      double zi = 0.0;
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesi(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselI for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::numeric_limits<double>::quiet_NaN();
      } else {
        result[i] = cyr[0];
      }
    }
    return result;
  }

  // Function for complex Bessel I
  std::vector<std::complex<double>> BesselI_complex_std(const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<std::complex<double>> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1;

    if (nu < 0) {
      if (nu == round(nu)) {
        return BesselI_complex_std(z, -nu, expon_scaled, verbose);
      }
      std::vector<std::complex<double>> kf = BesselK_complex_std(z, -nu, expon_scaled, verbose);
      if (expon_scaled) {
        for (size_t i = 0; i < z.size(); ++i) {
          kf[i] *= exp(-std::abs(z[i]));
        }
      }
      std::vector<std::complex<double>> bi = BesselI_complex_std(z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bi[i] + 2 / M_PI * sin(M_PI * nu) * kf[i];
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i].real();
      double zi = z[i].imag();
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesi(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselI for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
      } else {
        result[i] = std::complex<double>(cyr[0], cyi[0]);
      }
    }
    return result;
  }

  // Function for real Bessel J
  std::vector<double> BesselJ_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<double> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1;

    if (nu < 0) {
      if (expon_scaled) {
        Rcpp::Rcerr << "'expon.scaled=TRUE' not implemented for nu < 0" << std::endl;
        return std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN());
      }
      std::vector<double> nu_vec(z.size(), -nu);
      std::vector<double> bessel_j_neg_nu = BesselJ_real_std(z, -nu, expon_scaled, verbose);
      std::vector<double> bessel_y_neg_nu = BesselY_real_std(z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bessel_j_neg_nu[i] * std::cos(M_PI * nu) - bessel_y_neg_nu[i] * std::sin(M_PI * nu);
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i];
      double zi = 0.0;
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesj(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselJ for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::numeric_limits<double>::quiet_NaN();
      } else {
        result[i] = cyr[0];
      }
    }
    return result;
  }

  std::vector<std::complex<double>> BesselJ_complex_std(const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<std::complex<double>> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1;

    if (nu < 0) {
      if (expon_scaled) {
        Rcpp::Rcerr << "'expon.scaled=TRUE' not implemented for nu < 0" << std::endl;
        return std::vector<std::complex<double>>(z.size(), std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      }
      std::vector<std::complex<double>> bessel_j_neg_nu = BesselJ_complex_std(z, -nu, expon_scaled, verbose);
      std::vector<std::complex<double>> bessel_y_neg_nu = BesselY_complex_std(z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bessel_j_neg_nu[i] * std::cos(M_PI * nu) - bessel_y_neg_nu[i] * std::sin(M_PI * nu);
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i].real();
      double zi = z[i].imag();
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesj(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselJ for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
      } else {
        result[i] = std::complex<double>(cyr[0], cyi[0]);
      }
    }
    return result;
  }

  // Function for real Bessel Y
  std::vector<double> BesselY_real_std(const std::vector<double>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<double> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1;

    if (nu < 0) {
      if (expon_scaled) {
        Rcpp::Rcerr << "'expon.scaled=TRUE' not implemented for nu < 0" << std::endl;
        return std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN());
      }
      std::vector<double> nu_vec(z.size(), -nu);
      std::vector<double> bessel_y_neg_nu = BesselY_real_std(z, -nu, expon_scaled, verbose);
      std::vector<double> bessel_j_neg_nu = BesselJ_real_std(z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bessel_y_neg_nu[i] * std::cos(M_PI * nu) + bessel_j_neg_nu[i] * std::sin(M_PI * nu);
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      if (z[i] == 0) {
        result[i] = -std::numeric_limits<double>::infinity();
        continue;
      }
      double zr = z[i];
      double zi = 0.0;
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      std::vector<double> cwrkr(nSeq);
      std::vector<double> cwrki(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesy(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, cwrkr.data(), cwrki.data(), &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselY for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::numeric_limits<double>::quiet_NaN();
      } else {
        result[i] = cyr[0];
      }
    }
    return result;
  }

  // Function for complex Bessel Y
  std::vector<std::complex<double>> BesselY_complex_std(const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<std::complex<double>> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1;

    if (nu < 0) {
      if (expon_scaled) {
        Rcpp::Rcerr << "'expon.scaled=TRUE' not implemented for nu < 0" << std::endl;
        return std::vector<std::complex<double>>(z.size(), std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      }
      std::vector<std::complex<double>> bessel_y_neg_nu = BesselY_complex_std(z, -nu, expon_scaled, verbose);
      std::vector<std::complex<double>> bessel_j_neg_nu = BesselJ_complex_std(z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bessel_y_neg_nu[i] * std::cos(M_PI * nu) + bessel_j_neg_nu[i] * std::sin(M_PI * nu);
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      if (z[i] == std::complex<double>(0.0, 0.0)) {
        result[i] = std::complex<double>(-std::numeric_limits<double>::infinity(), 0.0);
        continue;
      }
      double zr = z[i].real();
      double zi = z[i].imag();
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      std::vector<double> cwrkr(nSeq);
      std::vector<double> cwrki(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesy(&zr, &zi, &nu, &kode, &nSeq, cyr.data(), cyi.data(), &nz, cwrkr.data(), cwrki.data(), &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselY for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
      } else {
        result[i] = std::complex<double>(cyr[0], cyi[0]);
      }
    }
    return result;
  }

  std::vector<std::complex<double>> BesselH_real_std(int m, const std::vector<double>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<std::complex<double>> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1;

    if (nu < 0) {
      if (expon_scaled) {
        Rcpp::Rcerr << "'expon.scaled=TRUE' not implemented for nu < 0" << std::endl;
        return std::vector<std::complex<double>>(z.size(), std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      }
      std::complex<double> pnu = ((m == 1) ? std::complex<double>(0, 1) : std::complex<double>(0, -1)) * M_PI * (-nu);
      std::vector<std::complex<double>> bessel_h_neg_nu = BesselH_real_std(m, z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bessel_h_neg_nu[i] * std::exp(pnu);
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i];
      double zi = 0.0;
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesh(&zr, &zi, &nu, &kode, &m, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselH for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
      } else {
        result[i] = std::complex<double>(cyr[0], cyi[0]);
      }
    }
    return result;
  }

  // Function for complex vectors
  std::vector<std::complex<double>> BesselH_complex_std(int m, const std::vector<std::complex<double>>& z, double nu, bool expon_scaled, int verbose) {
    std::vector<std::complex<double>> result(z.size());
    int kode = expon_scaled ? 2 : 1;
    int nSeq = 1; // Always 1 as per the assumption

    if (nu < 0) {
      if (expon_scaled) {
        Rcpp::Rcerr << "'expon.scaled=TRUE' not implemented for nu < 0" << std::endl;
        return std::vector<std::complex<double>>(z.size(), std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
      }
      std::complex<double> pnu = ((m == 1) ? std::complex<double>(0, 1) : std::complex<double>(0, -1)) * M_PI * (-nu);
      std::vector<std::complex<double>> bessel_h_neg_nu = BesselH_complex_std(m, z, -nu, expon_scaled, verbose);
      for (size_t i = 0; i < z.size(); ++i) {
        result[i] = bessel_h_neg_nu[i] * std::exp(pnu);
      }
      return result;
    }

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i].real();
      double zi = z[i].imag();
      std::vector<double> cyr(nSeq);
      std::vector<double> cyi(nSeq);
      int nz = 0;
      int ierr = 0;
      zbesh(&zr, &zi, &nu, &kode, &m, &nSeq, cyr.data(), cyi.data(), &nz, &ierr);
      if (ierr != 0) {
        if (verbose) {
          Rcpp::Rcerr << "Error computing BesselH for z[" << i << "]=" << z[i] << ": ierr=" << ierr << std::endl;
        }
        result[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
      } else {
        result[i] = std::complex<double>(cyr[0], cyi[0]);
      }
    }
    return result;
  }

  std::vector<double> AiryA_real_std(const std::vector<double>& z, int deriv, bool expon_scaled, int verbose) {
    if (deriv != 0 && deriv != 1) {
      Rcpp::Rcerr << "Invalid value for deriv. It should be either 0 or 1." << std::endl;
      return std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN());
    }
    std::vector<double> result(z.size());
    int kode = expon_scaled ? 2 : 1;

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i];
      double zi = 0.0;
      double air = 0.0;
      double aii = 0.0;
      int nz = 0;
      int ierr = 0;
      zairy(&zr, &zi, &deriv, &kode, &air, &aii, &nz, &ierr);
      if (ierr != 0) {
        std::string f_x = "zairy(" + std::to_string(zr) + " + " + std::to_string(zi) + "i, deriv=" + std::to_string(deriv) + ")";
        if (ierr == 3) {
          Rcpp::Rcerr << f_x << " large arguments -> precision loss (of at least half machine accuracy)" << std::endl;
        } else if (ierr == 2) {
          if (verbose) {
            Rcpp::Rcout << f_x << " -> overflow; returning Inf" << std::endl;
          }
          air = std::numeric_limits<double>::infinity();
          aii = std::numeric_limits<double>::infinity();
        } else if (ierr == 4) {
          Rcpp::Rcerr << f_x << " -> ierr=4: |z| too large" << std::endl;
          air = std::numeric_limits<double>::quiet_NaN();
          aii = (zi == 0.0) ? 0.0 : std::numeric_limits<double>::quiet_NaN();
        } else {
          Rcpp::Rcerr << f_x << " unexpected error 'ierr = " << ierr << "'" << std::endl;
          air = std::numeric_limits<double>::quiet_NaN();
          aii = std::numeric_limits<double>::quiet_NaN();
        }
      }
      result[i] = air;
    }
    return result;
  }

  // Function for complex vectors
  std::vector<std::complex<double>> AiryA_complex_std(const std::vector<std::complex<double>>& z, int deriv, bool expon_scaled, int verbose) {
    if (deriv != 0 && deriv != 1) {
      Rcpp::Rcerr << "Invalid value for deriv. It should be either 0 or 1." << std::endl;
      return std::vector<std::complex<double>>(z.size(), std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
    }
    std::vector<std::complex<double>> result(z.size());
    int kode = expon_scaled ? 2 : 1;

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i].real();
      double zi = z[i].imag();
      double air = 0.0;
      double aii = 0.0;
      int nz = 0;
      int ierr = 0;
      zairy(&zr, &zi, &deriv, &kode, &air, &aii, &nz, &ierr);
      if (ierr != 0) {
        std::string f_x = "zairy(" + std::to_string(zr) + " + " + std::to_string(zi) + "i, deriv=" + std::to_string(deriv) + ")";
        if (ierr == 3) {
          Rcpp::Rcerr << f_x << " large arguments -> precision loss (of at least half machine accuracy)" << std::endl;
        } else if (ierr == 2) {
          if (verbose) {
            Rcpp::Rcout << f_x << " -> overflow; returning Inf" << std::endl;
          }
          air = std::numeric_limits<double>::infinity();
          aii = std::numeric_limits<double>::infinity();
        } else if (ierr == 4) {
          Rcpp::Rcerr << f_x << " -> ierr=4: |z| too large" << std::endl;
          air = std::numeric_limits<double>::quiet_NaN();
          aii = std::numeric_limits<double>::quiet_NaN();
        } else {
          Rcpp::Rcerr << f_x << " unexpected error 'ierr = " << ierr << "'" << std::endl;
          air = std::numeric_limits<double>::quiet_NaN();
          aii = std::numeric_limits<double>::quiet_NaN();
        }
      }
      result[i] = std::complex<double>(air, aii);
    }
    return result;
  }

  std::vector<double> AiryB_real_std(const std::vector<double>& z, int deriv, bool expon_scaled, int verbose) {
    if (deriv != 0 && deriv != 1) {
      Rcpp::Rcerr << "Invalid value for deriv. It should be either 0 or 1." << std::endl;
      return std::vector<double>(z.size(), std::numeric_limits<double>::quiet_NaN());
    }

    std::vector<double> result(z.size());
    int kode = expon_scaled ? 2 : 1;

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i];
      double zi = 0.0;
      double bir = 0.0;
      double bii = 0.0;
      int ierr = 0;
      zbiry(&zr, &zi, &deriv, &kode, &bir, &bii, &ierr);
      if (ierr != 0) {
        std::string f_x = "zbiry(" + std::to_string(zr) + " + " + std::to_string(zi) + "i, deriv=" + std::to_string(deriv) + ")";
        if (ierr == 3) {
          Rcpp::Rcerr << f_x << " large arguments -> precision loss (of at least half machine accuracy)" << std::endl;
        } else if (ierr == 2) {
          if (verbose) {
            Rcpp::Rcerr << f_x << " -> overflow; returning Inf" << std::endl;
          }
          bir = std::numeric_limits<double>::infinity();
          bii = std::numeric_limits<double>::infinity();
        } else if (ierr == 4) {
          Rcpp::Rcerr << f_x << " -> ierr=4: |z| too large" << std::endl;
          bir = std::numeric_limits<double>::quiet_NaN();
          bii = 0.0;
        } else {
          Rcpp::Rcerr << f_x << " unexpected error 'ierr = " << ierr << "'" << std::endl;
          bir = std::numeric_limits<double>::quiet_NaN();
          bii = std::numeric_limits<double>::quiet_NaN();
        }
      }
      result[i] = bir;
    }
    return result;
  }

  // Function for complex vectors
  std::vector<std::complex<double>> AiryB_complex_std(const std::vector<std::complex<double>>& z, int deriv, bool expon_scaled, int verbose) {
    if (deriv != 0 && deriv != 1) {
      Rcpp::Rcerr << "Invalid value for deriv. It should be either 0 or 1." << std::endl;
      return std::vector<std::complex<double>>(z.size(), std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
    }

    std::vector<std::complex<double>> result(z.size());
    int kode = expon_scaled ? 2 : 1;

    for (size_t i = 0; i < z.size(); ++i) {
      double zr = z[i].real();
      double zi = z[i].imag();
      double bir = 0.0;
      double bii = 0.0;
      int ierr = 0;
      zbiry(&zr, &zi, &deriv, &kode, &bir, &bii, &ierr);
      if (ierr != 0) {
        std::string f_x = "zbiry(" + std::to_string(zr) + " + " + std::to_string(zi) + "i, deriv=" + std::to_string(deriv) + ")";
        if (ierr == 3) {
          Rcpp::Rcerr << f_x << " large arguments -> precision loss (of at least half machine accuracy)" << std::endl;
        } else if (ierr == 2) {
          if (verbose) {
            Rcpp::Rcerr << f_x << " -> overflow; returning Inf" << std::endl;
          }
          bir = std::numeric_limits<double>::infinity();
          bii = std::numeric_limits<double>::infinity();
        } else if (ierr == 4) {
          Rcpp::Rcerr << f_x << " -> ierr=4: |z| too large" << std::endl;
          bir = std::numeric_limits<double>::quiet_NaN();
          bii = std::numeric_limits<double>::quiet_NaN();
        } else {
          Rcpp::Rcerr << f_x << " unexpected error 'ierr = " << ierr << "'" << std::endl;
          bir = std::numeric_limits<double>::quiet_NaN();
          bii = std::numeric_limits<double>::quiet_NaN();
        }
      }
      result[i] = std::complex<double>(bir, bii);
    }
    return result;
  }
}

