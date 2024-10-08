// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/RcppBessel.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BesselK
SEXP BesselK(SEXP z, double nu, bool expon_scaled, int verbose);
static SEXP _RcppBessel_BesselK_try(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type expon_scaled(expon_scaledSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(BesselK(z, nu, expon_scaled, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _RcppBessel_BesselK(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_RcppBessel_BesselK_try(zSEXP, nuSEXP, expon_scaledSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// BesselI
SEXP BesselI(SEXP z, double nu, bool expon_scaled, int verbose);
static SEXP _RcppBessel_BesselI_try(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type expon_scaled(expon_scaledSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(BesselI(z, nu, expon_scaled, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _RcppBessel_BesselI(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_RcppBessel_BesselI_try(zSEXP, nuSEXP, expon_scaledSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// BesselJ
SEXP BesselJ(SEXP z, double nu, bool expon_scaled, int verbose);
static SEXP _RcppBessel_BesselJ_try(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type expon_scaled(expon_scaledSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(BesselJ(z, nu, expon_scaled, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _RcppBessel_BesselJ(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_RcppBessel_BesselJ_try(zSEXP, nuSEXP, expon_scaledSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// BesselY
SEXP BesselY(SEXP z, double nu, bool expon_scaled, int verbose);
static SEXP _RcppBessel_BesselY_try(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type expon_scaled(expon_scaledSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(BesselY(z, nu, expon_scaled, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _RcppBessel_BesselY(SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_RcppBessel_BesselY_try(zSEXP, nuSEXP, expon_scaledSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// BesselH
SEXP BesselH(int m, SEXP z, double nu, bool expon_scaled, int verbose);
static SEXP _RcppBessel_BesselH_try(SEXP mSEXP, SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type expon_scaled(expon_scaledSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(BesselH(m, z, nu, expon_scaled, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _RcppBessel_BesselH(SEXP mSEXP, SEXP zSEXP, SEXP nuSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_RcppBessel_BesselH_try(mSEXP, zSEXP, nuSEXP, expon_scaledSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// AiryA
SEXP AiryA(SEXP z, int deriv, bool expon_scaled, int verbose);
static SEXP _RcppBessel_AiryA_try(SEXP zSEXP, SEXP derivSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type deriv(derivSEXP);
    Rcpp::traits::input_parameter< bool >::type expon_scaled(expon_scaledSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(AiryA(z, deriv, expon_scaled, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _RcppBessel_AiryA(SEXP zSEXP, SEXP derivSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_RcppBessel_AiryA_try(zSEXP, derivSEXP, expon_scaledSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// AiryB
SEXP AiryB(SEXP z, int deriv, bool expon_scaled, int verbose);
static SEXP _RcppBessel_AiryB_try(SEXP zSEXP, SEXP derivSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type deriv(derivSEXP);
    Rcpp::traits::input_parameter< bool >::type expon_scaled(expon_scaledSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(AiryB(z, deriv, expon_scaled, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _RcppBessel_AiryB(SEXP zSEXP, SEXP derivSEXP, SEXP expon_scaledSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_RcppBessel_AiryB_try(zSEXP, derivSEXP, expon_scaledSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _RcppBessel_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("SEXP(*bessel_k)(SEXP,double,bool,int)");
        signatures.insert("SEXP(*bessel_i)(SEXP,double,bool,int)");
        signatures.insert("SEXP(*bessel_j)(SEXP,double,bool,int)");
        signatures.insert("SEXP(*bessel_y)(SEXP,double,bool,int)");
        signatures.insert("SEXP(*bessel_h)(int,SEXP,double,bool,int)");
        signatures.insert("SEXP(*airy_a)(SEXP,int,bool,int)");
        signatures.insert("SEXP(*airy_b)(SEXP,int,bool,int)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _RcppBessel_RcppExport_registerCCallable() { 
    R_RegisterCCallable("RcppBessel", "_RcppBessel_bessel_k", (DL_FUNC)_RcppBessel_BesselK_try);
    R_RegisterCCallable("RcppBessel", "_RcppBessel_bessel_i", (DL_FUNC)_RcppBessel_BesselI_try);
    R_RegisterCCallable("RcppBessel", "_RcppBessel_bessel_j", (DL_FUNC)_RcppBessel_BesselJ_try);
    R_RegisterCCallable("RcppBessel", "_RcppBessel_bessel_y", (DL_FUNC)_RcppBessel_BesselY_try);
    R_RegisterCCallable("RcppBessel", "_RcppBessel_bessel_h", (DL_FUNC)_RcppBessel_BesselH_try);
    R_RegisterCCallable("RcppBessel", "_RcppBessel_airy_a", (DL_FUNC)_RcppBessel_AiryA_try);
    R_RegisterCCallable("RcppBessel", "_RcppBessel_airy_b", (DL_FUNC)_RcppBessel_AiryB_try);
    R_RegisterCCallable("RcppBessel", "_RcppBessel_RcppExport_validate", (DL_FUNC)_RcppBessel_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppBessel_BesselK", (DL_FUNC) &_RcppBessel_BesselK, 4},
    {"_RcppBessel_BesselI", (DL_FUNC) &_RcppBessel_BesselI, 4},
    {"_RcppBessel_BesselJ", (DL_FUNC) &_RcppBessel_BesselJ, 4},
    {"_RcppBessel_BesselY", (DL_FUNC) &_RcppBessel_BesselY, 4},
    {"_RcppBessel_BesselH", (DL_FUNC) &_RcppBessel_BesselH, 5},
    {"_RcppBessel_AiryA", (DL_FUNC) &_RcppBessel_AiryA, 4},
    {"_RcppBessel_AiryB", (DL_FUNC) &_RcppBessel_AiryB, 4},
    {"_RcppBessel_RcppExport_registerCCallable", (DL_FUNC) &_RcppBessel_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppBessel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
