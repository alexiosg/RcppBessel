#' The BesselK Function
#'
#' Computes the modified Bessel function of the second kind for real or complex inputs.
#' @usage bessel_k(z, nu, expon_scaled = FALSE, verbose = 0)
#' @param z A numeric or complex vector representing the input values at which to evaluate the Bessel function.
#' @param nu A double representing the order of the Bessel function.
#' @param expon_scaled A logical value indicating whether to use the exponentially scaled form of the Bessel function. Defaults to `FALSE`.
#' @param verbose An integer specifying the verbosity level for error messages. Defaults to `0`.
#'
#' @return A numeric or complex vector (depending on the input) containing the values of the \code{bessel_k} function evaluated at the points in \code{z}.
#'
#' @name bessel_k
#' @references {
#' \insertRef{Maechler2024}{RcppBessel}\cr
#' \insertRef{Amos1995}{RcppBessel}
#' }
#' @export
NULL

#' The BesselI Function
#'
#' Computes the modified Bessel function of the first kind for real or complex inputs.
#' @usage bessel_i(z, nu, expon_scaled = FALSE, verbose = 0)
#' @param z A numeric or complex vector representing the input values at which to evaluate the Bessel function.
#' @param nu A double representing the order of the Bessel function.
#' @param expon_scaled A logical value indicating whether to use the exponentially scaled form of the Bessel function. Defaults to `FALSE`.
#' @param verbose An integer specifying the verbosity level for error messages. Defaults to `0`.
#'
#' @return A numeric or complex vector (depending on the input) containing the values of the \code{bessel_i} function evaluated at the points in \code{z}.
#' @references {
#' \insertRef{Maechler2024}{RcppBessel}\cr
#' \insertRef{Amos1995}{RcppBessel}
#' }
#' @name bessel_i
#' @export
NULL

#' The BesselY Function
#'
#' Computes the Bessel function of the second kind (Neumann function) for real or complex inputs.
#' @usage bessel_y(z, nu, expon_scaled = FALSE, verbose = 0)
#' @param z A numeric or complex vector representing the input values at which to evaluate the Bessel function.
#' @param nu A double representing the order of the Bessel function.
#' @param expon_scaled A logical value indicating whether to use the exponentially scaled form of the Bessel function. Defaults to `FALSE`.
#' @param verbose An integer specifying the verbosity level for error messages. Defaults to `0`.
#'
#' @return A numeric or complex vector (depending on the input) containing the values of the \code{bessel_y} function evaluated at the points in \code{z}.
#' @references {
#' \insertRef{Maechler2024}{RcppBessel}\cr
#' \insertRef{Amos1995}{RcppBessel}
#' }
#' @name bessel_y
#' @export
NULL

#' The BesselJ Function
#'
#' Computes the Bessel function of the first kind for real or complex inputs.
#' @usage bessel_j(z, nu, expon_scaled = FALSE, verbose = 0)
#' @param z A numeric or complex vector representing the input values at which to evaluate the Bessel function.
#' @param nu A double representing the order of the Bessel function.
#' @param expon_scaled A logical value indicating whether to use the exponentially scaled form of the Bessel function. Defaults to `FALSE`.
#' @param verbose An integer specifying the verbosity level for error messages. Defaults to `0`.
#'
#' @return A numeric or complex vector (depending on the input) containing the values of the \code{bessel_j} function evaluated at the points in \code{z}.
#' @references {
#' \insertRef{Maechler2024}{RcppBessel}\cr
#' \insertRef{Amos1995}{RcppBessel}
#' }
#' @name bessel_j
#' @export
NULL

#' The BesselH Function
#'
#' Computes the Hankel function (Bessel function of the third kind) for real or complex inputs.
#' @usage bessel_h(m, z, nu, expon_scaled = FALSE, verbose = 0)
#' @param m An integer representing the type of Hankel function. It must be either `1` (for the first kind) or `2` (for the second kind).
#' @param z A numeric or complex vector representing the input values at which to evaluate the Hankel function.
#' @param nu A double representing the order of the Hankel function.
#' @param expon_scaled A logical value indicating whether to use the exponentially scaled form of the Hankel function. Defaults to `FALSE`.
#' @param verbose An integer specifying the verbosity level for error messages. Defaults to `0`.
#'
#' @return A complex vector containing the values of the \code{bessel_h} function evaluated at the points in \code{z}.
#' @references {
#' \insertRef{Maechler2024}{RcppBessel}\cr
#' \insertRef{Amos1995}{RcppBessel}
#' }
#' @name bessel_h
#' @export
NULL

#' The AiryA Function
#'
#' Computes the Airy function Ai for real or complex inputs.
#' @usage airy_a(z, deriv = 0, expon_scaled = FALSE, verbose = 0)
#' @param z A numeric or complex vector representing the input values at which to evaluate the Airy function.
#' @param deriv An integer indicating whether to compute the function (`0` for the function itself) or its first derivative (`1` for the first derivative). Defaults to `0`.
#' @param expon_scaled A logical value indicating whether to use the exponentially scaled form of the Airy function. Defaults to `FALSE`.
#' @param verbose An integer specifying the verbosity level for error messages. Defaults to `0`.
#'
#' @return A numeric or complex vector (depending on the input) containing the values of the \code{airy_a} function evaluated at the points in \code{z}.
#' @references {
#' \insertRef{Maechler2024}{RcppBessel}\cr
#' \insertRef{Amos1995}{RcppBessel}
#' }
#' @name airy_a
#' @export
NULL

#' The AiryB Function
#'
#' Computes the Airy function Bi for real or complex inputs.
#' @usage airy_b(z, deriv = 0, expon_scaled = FALSE, verbose = 0)
#' @param z A numeric or complex vector representing the input values at which to evaluate the Airy function.
#' @param deriv An integer indicating whether to compute the function (`0` for the function itself) or its first derivative (`1` for the first derivative). Defaults to `0`.
#' @param expon_scaled A logical value indicating whether to use the exponentially scaled form of the Airy function. Defaults to `FALSE`.
#' @param verbose An integer specifying the verbosity level for error messages. Defaults to `0`.
#'
#' @return A numeric or complex vector (depending on the input) containing the values of the \code{airy_b} function evaluated at the points in \code{z}.
#' @references {
#' \insertRef{Maechler2024}{RcppBessel}\cr
#' \insertRef{Amos1995}{RcppBessel}
#' }
#' @name airy_b
#' @export
NULL

