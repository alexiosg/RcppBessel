# RcppBessel

The `RcppBessel` package exports headers for the Bessel functions found in the 
`Bessel` package of Maechler by wrapping the functions in the C code (translated
and cleaned up from Amos' original FORTRAN code by Maechler) into C++. Not all 
functionality is exposed, in order to keep things simple, but may be expanded 
in the future to export more of the functions such as those for the asymptotic 
expansion of BesselK and BesselI for large nu and x. The package also exports 
the functions for use in R, as it was found to offer a speed up of up to 10x.

The initial motivation for writing this package was so as to allow the modified
Bessel function of the second kind to be used in C++ code which the author needed
in another package. The [Boost](https://www.boost.org/doc/libs/1_85_0/libs/math/doc/html/math_toolkit/bessel/mbessel.html) 
implementation does not allow for complex arguments
whilst the special math functions in [std library](https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_k)
is not available for Mac OS (via the shipped libc++).  Additionally, neither of those 
implementation has the option of exponential scaling.

# Install

Until the package makes its way to CRAN, you can install it as follows:

``` r
remotes::install_github("alexiosg/RcppBessel", dependencies = TRUE)
```
