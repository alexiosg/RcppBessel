/* -*- mode: c; kept-new-versions: 25; kept-old-versions: 20 -*- */

/* MM: --- Renaming  zexp, zlog, zsqrt -- exist as *functions* in g77 & gfortran
 *** zbsubs.f
 -----------------------------------------------------------------
 >>>  zbsubs.for:  double precision subroutines
 ----------------------------------------------------------------- */

/* Produced by f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 *
 * zbsubs.f -- translated by f2c (version 20061008). http://www.netlib.org/f2c/
*/

#include <math.h>
#include <float.h> /* DBL_MIN etc */
#include <limits.h>
#include <stdlib.h> // abs()

/* possibly needed for 'verbose' and debugging */
#include <R_ext/Print.h>

#include <Rmath.h>

#include "zbsubs.h"

/* Table of constant values */

static int c__0 = 0;
static int c__1 = 1;
static int c__2 = 2;
static double pi  = 3.141592653589793238462643383; // = pi
static double hpi = 1.570796326794896619231321696; // = pi/2
static double thpi= 4.7123889803846898576939650749;// = 3/2 pi
static double rthpi=1.2533141373155002512078826424;// = sqrt(pi/2)
static double aic = 1.265512123484645396488945797; // = log( 2*sqrt(pi) ) = log(2) + log(pi)/2

static double ufl = DBL_MIN * 1e3; // = 2.225..e-305 ;  ufl := dexp(-elim)


void zbesh(double *zr, double *zi, double *fnu, int *kode, int *m, int *n,
	   double *cyr, double *cyi, int *nz, int *ierr)
{
/***begin prologue  zbesh
 ***date written   830501   (yymmdd)
 ***revision date  890801, 930101   (yymmdd)
 ***category no.  b5k
 ***keywords  h-bessel functions,bessel functions of complex argument,
             bessel functions of third kind,hankel functions
 ***author  amos, donald e., sandia national laboratories
 ***purpose  to compute the h-bessel functions of a complex argument
 ***description

                      ***a double precision routine***

       - on kode=1, zbesh computes an n member sequence of complex
         hankel (bessel) functions cy(j) = H(m,fnu+j-1,z) for kinds m = 1 or 2,
	 float, nonnegative orders fnu+j-1, j=1,...,n, and complex
         z != cmplx(0.,0.) in the cut plane -pi < arg(z) <= pi.

       - on kode=2, zbesh returns the scaled hankel functions

              cy(j) = exp(-mm*z*i) * H(m,fnu+j-1,z)
                mm := 3-2*m  (\in {-1, 1}),          and  i^2 = -1.

         which removes the exponential behavior in both the upper and
         lower half planes.  Definitions and notation are found in the
         NBS handbook of mathematical functions (ref. 1).


         input      zr,zi,fnu are double precision
           zr,zi  - z=cmplx(zr,zi), z != cmplx(0,0),
                    -pt < arg(z) <= pi
           fnu    - order of initial h function, fnu >= 0
           kode   - a parameter to indicate the scaling option
                    kode= 1  returns
                             cy(j)=H(m,fnu+j-1,z),   j=1,...,n
                        = 2  returns
                             cy(j)=H(m,fnu+j-1,z)*exp(-i*z*(3-2m))
                                  j=1,...,n  ,  i^2=-1
           m      - kind of hankel function, m=1 or 2
           n      - number of members in the sequence, n >= 1

         output     cyr,cyi are double precision
           cyr,cyi- double precision vectors whose first n components
                    contain float and imaginary parts for the sequence
                    cy(j)=H(m,fnu+j-1,z)  or
                    cy(j)=H(m,fnu+j-1,z)*exp(-i*z*(3-2m))  j=1,...,n
                    depending on kode, i^2=-1.
           nz     - number of components set to zero due to underflow,
                    nz= 0   , normal return
                    nz > 0 , first nz components of cy set to zero due
                              to underflow, cy(j)=cmplx(0,0)
                              j=1,...,nz when y > 0 and m=1 or
                              y < 0 and m=2. for the complmentary
                              half planes, nz states only the number
                              of underflows.
           ierr   - error flag
                    ierr=0, normal return - computation completed
                    ierr=1, input error   - no computation
                    ierr=2, overflow      - no computation, fnu too
                            large or cabs(z) too small or both
                    ierr=3, cabs(z) or fnu+n-1 large - computation done
                            but losses of signifcance by argument
                            reduction produce less than half of machine
                            accuracy
                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
                            tion because of complete losses of signifi-
                            cance by argument reduction
                    ierr=5, error              - no computation,
                            algorithm termination condition not met

 ***long description

         the computation is carried out by the relation

         H(m,fnu,z) = (1/mp)*exp(-mp*fnu)*K(fnu,z*exp(-mp))
             mp=mm*hpi*i,  mm=3-2*m,  hpi=pi/2,  i^2=-1

         for m=1 or 2 where the k bessel function is computed for the
         right half plane re(z) >= 0.  The k function is continued
         to the left half plane by the relation

         K(fnu,z*exp(mp)) = exp(-mp*fnu)*K(fnu,z) - mp*I(fnu,z)
         mp=mr*pi*i, mr=+1 or -1, re(z) > 0, i^2=-1

         where I(fnu,z) is the i bessel function.

         exponential decay of H(m,fnu,z) occurs in the upper half z
         plane for m=1 and the lower half z plane for m=2.  exponential
         growth occurs in the complementary half planes.  scaling
         by exp(-mm*z*i) removes the exponential behavior in the
         whole z plane for z to infinity.

         for negative orders,the formulae

               H(1,-fnu,z) = H(1,fnu,z)*cexp( pi*fnu*i)
               H(2,-fnu,z) = H(2,fnu,z)*cexp(-pi*fnu*i)     i^2 = -1

         can be used.

         in most complex variable computation, one must evaluate ele-
         mentary functions. when the magnitude of z or fnu+n-1 is
         large, losses of significance by argument reduction occur.
         consequently, if either one exceeds u1=sqrt(0.5/ur), then
         losses exceeding half precision are likely and an error flag
         ierr=3 is triggered where ur=dmax1(DBL_EPSILON,1.0d-18) is
         double precision unit roundoff limited to 18 digits precision.
         if either is larger than u2=0.5/ur, then all significance is
         lost and ierr=4. in order to use the int function, arguments
         must be further restricted not to exceed the largest machine
         int, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double precision
         arithmetic respectively. this makes u2 and u3 limiting in
         their respective arithmetics. this means that one can expect
         to retain, in the worst cases on 32 bit machines, no digits
         in single and only 7 digits in double precision arithmetic.
         similar considerations hold for other machines.

         the approximate relative error in the magnitude of a complex
         bessel function can be expressed by p*10**s where p=max(unit
         roundoff,1.0d-18) is the nominal precision and 10**s repre-
         sents the increase in error due to argument reduction in the
         elementary functions. here, s=max(1,abs(log10(cabs(z))),
         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
         have only absolute accuracy. this is most likely to occur when
         one component (in absolute value) is larger than the other by
         several orders of magnitude. if one component is 10**k larger
         than the other, then one can expect only max(abs(log10(p))-k,
         0) significant digits; or, stated another way, when k exceeds
         the exponent of p, no significant digits remain in the smaller
         component. however, the phase angle retains absolute accuracy
         because, in complex arithmetic with precision p, the smaller
         component will not (as a rule) decrease below p times the
         magnitude of the larger component. in these extreme cases,
         the principal phase angle is on the order of +p, -p, pi/2-p,
         or -pi/2+p.

 ***references  handbook of mathematical functions by m. abramowitz
                 and i. a. stegun, nbs ams series 55, u.s. dept. of
                 commerce, 1955.

               computation of bessel functions of complex argument
                 by d. e. amos, sand83-0083, may, 1983.

               computation of bessel functions of complex argument
                 and large order by d. e. amos, sand83-0643, may, 1983

               a subroutine package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, sand85-
                 1018, may, 1985

               a portable package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, acm
                 trans. math. software, vol. 12, no. 3, september 1986,
                 pp 265-273.

 ***routines called  zacon,zbknu,zbunk,zuoik,zabs,i1mach,d1mach
 ***end prologue  zbesh
 */

    /* Local variables */
    double arg, str, sti, atol;
    /* complex cy,z,zn,zt,csgn */

    int verbose = *ierr;

    /* ***first executable statement  zbesh */
    *nz = 0;
    if (*zr == 0. && *zi == 0.) {
	*zi = DBL_MIN;// FIXME -- better?  result strongly depends on  nu !!
    }
    if (*fnu < 0. || *m < 1 || *m > 2 || *kode < 1 || *kode > 2 || *n < 1) {
	*ierr = 1; return;
    } else
	*ierr = 0;
    int nn = *n, nw;
/* -----------------------------------------------------------------------
     set parameters related to machine constants.
     tol is the approximate unit roundoff limited to 1.0e-18.
     elim is the approximate exponential over- and underflow limit.
     exp(-elim) < exp(-alim)=exp(-elim)/tol    and
     exp(elim) > exp(alim)=exp(elim)*tol       are intervals near
     underflow and overflow limits where scaled arithmetic is done.
     rl is the lower boundary of the asymptotic expansion for large z.
     dig = number of base 10 digits in tol = 10**(-dig).
     fnul is the lower boundary of the asymptotic series for large fnu
 ----------------------------------------------------------------------- */
    double
	tol  = fmax2(DBL_EPSILON, 1e-18),
	r1m5 = M_LOG10_2,
	k    = fmin2(-DBL_MIN_EXP, DBL_MAX_EXP),
	elim = (k * r1m5 - 3.) * 2.303,
	aa = r1m5 * (DBL_MANT_DIG - 1),
	dig = fmin2(aa, 18.);
    aa *= 2.303;
    double
	alim = elim + fmax2(-aa, -41.45),
	fnul = (dig - 3.) * 6. + 10.,
	rl = dig * 1.2 + 3.,
	fn = *fnu + (double) (nn - 1);
    int mm = 3 - *m - *m;
    double fmm = (double) mm,
	znr = fmm * *zi,
	zni = -fmm * *zr,
/* -----------------------------------------------------------------------
     test for proper range
 ----------------------------------------------------------------------- */
	az = zabs(*zr, *zi),
	bb = ((double) INT_MAX) * .5;
    aa = fmin2(.5 / tol, bb);
    if (az > aa || fn > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa || fn > aa) {
	*ierr = 3;
    }
/* -----------------------------------------------------------------------
     overflow test on the last member of the sequence
 ----------------------------------------------------------------------- */
    if (az < ufl) {
	goto L230;
    }
    if (*fnu > fnul) { // L90:
	/* ---------------------------------------------------------------
	   uniform asymptotic expansions for fnu > fnul
	   --------------------------------------------------------------- */
	int mr = 0;
	if (!(znr >= 0. && (znr != 0. || zni >= 0. || *m != 2))) {
	    mr = -mm;
	    if (znr == 0. && zni < 0.) {
		znr = -znr;
		zni = -zni;
	    }
	}
	// L100:
	zbunk_(&znr, &zni, fnu, kode, &mr, &nn, cyr, cyi, &nw, &tol,
	       &elim, &alim);
	if (nw < 0) {
	    goto L240;
	}
	*nz += nw;

    } else { // *fnu <= fnul
	if (fn > 1.) {
	    if (fn <= 2.) {
		if (az <= tol) {
		    arg = az * .5;
		    double aln = -fn * log(arg);
		    if (aln > elim) {
			goto L230;
		    }
		}
		// goto L70;
	    } else { // fn > 2
		int nuf =
		    zuoik(&znr, &zni, fnu, *kode, 2, nn, cyr, cyi, tol, elim, alim);
		if (nuf < 0) {
		    goto L230;
		}
		*nz += nuf;
		nn -= nuf;
		/* --------------------------------------------------------------
		   here nn=n or nn=0 since nuf=0,nn, or -1 on return from cuoik
		   if nuf=nn, then cy(i)=czero for all i
		   -------------------------------------------------------------- */
		if (nn == 0) {
		    goto L140;
		}
	    }
	}
	// L70:
	if (znr < 0. || (znr == 0. && zni < 0. && *m == 2)) { // L80:
	    /* ---------------------------------------------------------------
	       left half plane computation
	       --------------------------------------------------------------- */
	    int mr = -mm;
	    zacon_(&znr, &zni, fnu, kode, &mr, &nn, cyr, cyi, &nw,
		   verbose, rl, fnul, tol, elim, alim);
	    if (nw < 0) {
		goto L240;
	    }
	    *nz = nw;
	} else {
	    /* ---------------------------------------------------------------
	       right half plane computation, xn >= 0. .and. (xn != 0. .or.
	       yn >= 0. .or. m=1)
	       --------------------------------------------------------------- */
	    *nz = zbknu(&znr, &zni, fnu, *kode, nn, verbose,
			cyr, cyi, tol, elim, alim);
	}
    }
    // L110:
/* -----------------------------------------------------------------------
     H(m,fnu,z) = -fmm*(i/hpi)*(zt**fnu)*K(fnu,-z*zt)

     zt=exp(-fmm*hpi*i) = cmplx(0.0,-fmm), fmm=3-2*m, m=1,2
 ----------------------------------------------------------------------- */
    double sgn = fsign(hpi, -fmm);
/* -----------------------------------------------------------------------
     calculate exp(fnu*hpi*i) to minimize losses of significance
     when fnu is large
 ----------------------------------------------------------------------- */
    int
	inu = (int) *fnu,
	inuh = inu / 2, // integer div: inuh = floor( nu/2 )
	ir = inu - (inuh << 1);
    arg = (*fnu - (double) ((float) (inu - ir))) * sgn;
    double
	rhpi = 1. / sgn,
	/*     zni = rhpi*dcos(arg)
	       znr = -rhpi*dsin(arg) */
	csgni = rhpi * cos(arg),
	csgnr = -rhpi * sin(arg);
    if (inuh % 2 == 1) { // nu/2 is odd  <==>  nu = 2(2k+1) (+1) = 4k + {2 or 3}
	/*     znr = -znr
	       zni = -zni */
	csgnr = -csgnr;
	csgni = -csgni;
    }
    // L120:
    double
	zti = -fmm,
	rtol = 1. / tol,
	ascle = ufl * rtol;
    for (int i = 0; i < nn; ++i) {
	/*
	  str = cyr(i)*znr - cyi(i)*zni
	  cyi(i) = cyr(i)*zni + cyi(i)*znr
	  cyr(i) = str
	  str = -zni*zti
	  zni = znr*zti
	  znr = str
	*/
	aa = cyr[i];
	bb = cyi[i];
	atol = 1.;
	if (fmax2(fabs(aa), fabs(bb)) <= ascle) {
	    aa *= rtol;
	    bb *= rtol;
	    atol = tol;
	}
	str = aa * csgnr - bb * csgni;
	sti = aa * csgni + bb * csgnr;
	cyr[i] = str * atol;
	cyi[i] = sti * atol;
	str = -csgni * zti;
	csgni = csgnr * zti;
	csgnr = str;
    }
    return;
L140:
    if (znr >= 0.)
	return;
L230:
    *nz = 0;	*ierr = 2;	return;
L240:
    if (nw == -1) {
	goto L230;
    }
    *nz = 0;	*ierr = 5;	return;
L260:
    *nz = 0;    *ierr = 4;	return;
} /* zbesh */

void zbesi(double *zr, double *zi, double *fnu, int *kode, int *n,
	   double *cyr, double *cyi, int *nz, int *ierr)
{
/***begin prologue  zbesi
 ***date written   830501   (yymmdd)
 ***revision date  890801, 930101   (yymmdd)
 ***category no.  b5k
 ***keywords  i-bessel function,complex bessel function,
             modified bessel function of the first kind
 ***author  amos, donald e., sandia national laboratories
 ***purpose  to compute i-bessel functions of complex argument
 ***description

                    ***a double precision routine***
         on kode=1, zbesi computes an n member sequence of complex
         bessel functions cy(j) = I(fnu+j-1,z) for float, nonnegative
         orders fnu+j-1, j=1,...,n and complex z in the cut plane
         -pi < arg(z) <= pi. on kode=2, zbesi returns the scaled
         functions

         cy(j) = exp(-abs(x)) * I(fnu+j-1,z)   j = 1,...,n , x=float(z)

         with the exponential growth removed in both the left and
         right half planes for z to infinity. definitions and notation
         are found in the nbs handbook of mathematical functions
         (ref. 1).

         input      zr,zi,fnu are double precision
           zr,zi  - z=cmplx(zr,zi),  -pi < arg(z) <= pi
           fnu    - order of initial I function, fnu >= 0
           kode   - a parameter to indicate the scaling option
              kode= 1  returns  cy(j) = I(fnu+j-1,z),              j=1,...,n
                  = 2  returns  cy(j) = I(fnu+j-1,z)*exp(-abs(x)), j=1,...,n

           n      - number of members of the sequence, n >= 1

         output     cyr,cyi are double precision
           cyr,cyi- double precision vectors whose first n components
                    contain float and imaginary parts for the sequence
                    cy(j) = I(fnu+j-1,z)  or
                    cy(j) = I(fnu+j-1,z)*exp(-abs(x))  j=1,...,n
                    depending on kode, x=float(z)
           nz     - number of components set to zero due to underflow,
                    nz= 0   , normal return
                    nz > 0 , last nz components of cy set to zero
                              to underflow, cy(j) = cmplx(0,0)
                              j = n-nz+1,...,n
           ierr   - error flag
                    ierr=0, normal return - computation completed
                    ierr=1, input error   - no computation
                    ierr=2, overflow      - no computation, float(z) too
                            large on kode=1
                    ierr=3, cabs(z) or fnu+n-1 large - computation done
                            but losses of signifcance by argument
                            reduction produce less than half of machine
                            accuracy
                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
                            tion because of complete losses of signifi-
                            cance by argument reduction
                    ierr=5, error              - no computation,
                            algorithm termination condition not met

 ***long description

         the computation is carried out by the power series for
         small cabs(z), the asymptotic expansion for large cabs(z),
         the miller algorithm normalized by the wronskian and a
         neumann series for imtermediate magnitudes, and the
         uniform asymptotic expansions for I(fnu,z) and J(fnu,z)
         for large orders. backward recurrence is used to generate
         sequences or reduce orders when necessary.

         the calculations above are done in the right half plane and
         continued into the left half plane by the formula

         I(fnu,z*exp(m*pi)) = exp(m*pi*fnu)*I(fnu,z)  float(z) > 0.0
                       m = +i or -i,  i**2=-1

         for negative orders,the formula

              I(-fnu,z) = I(fnu,z) + (2/pi)*sin(pi*fnu)*K(fnu,z)

         can be used. however,for large orders close to integers, the
         the function changes radically. when fnu is a large positive
         int,the magnitude of I(-fnu,z)=I(fnu,z) is a large
         negative power of ten. but when fnu is not an int,
         K(fnu,z) dominates in magnitude with a large positive power of
         ten and the most that the second term can be reduced is by
         unit roundoff from the coefficient. thus, wide changes can
         occur within unit roundoff of a large int for fnu. here,
         large means fnu > cabs(z).

         in most complex variable computation, one must evaluate ele-
         mentary functions. when the magnitude of z or fnu+n-1 is
         large, losses of significance by argument reduction occur.
         consequently, if either one exceeds u1=sqrt(0.5/ur), then
         losses exceeding half precision are likely and an error flag
         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
         double precision unit roundoff limited to 18 digits precision.
         if either is larger than u2=0.5/ur, then all significance is
         lost and ierr=4. in order to use the int function, arguments
         must be further restricted not to exceed the largest machine
         int, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double precision
         arithmetic respectively. this makes u2 and u3 limiting in
         their respective arithmetics. this means that one can expect
         to retain, in the worst cases on 32 bit machines, no digits
         in single and only 7 digits in double precision arithmetic.
         similar considerations hold for other machines.

         the approximate relative error in the magnitude of a complex
         bessel function can be expressed by p*10**s where p=max(unit
         roundoff,1.0e-18) is the nominal precision and 10**s repre-
         sents the increase in error due to argument reduction in the
         elementary functions. here, s=max(1,abs(log10(cabs(z))),
         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
         have only absolute accuracy. this is most likely to occur when
         one component (in absolute value) is larger than the other by
         several orders of magnitude. if one component is 10**k larger
         than the other, then one can expect only max(abs(log10(p))-k,
         0) significant digits; or, stated another way, when k exceeds
         the exponent of p, no significant digits remain in the smaller
         component. however, the phase angle retains absolute accuracy
         because, in complex arithmetic with precision p, the smaller
         component will not (as a rule) decrease below p times the
         magnitude of the larger component. in these extreme cases,
         the principal phase angle is on the order of +p, -p, pi/2-p,
         or -pi/2+p.

 ***references  handbook of mathematical functions by m. abramowitz
                 and i. a. stegun, nbs ams series 55, u.s. dept. of
                 commerce, 1955.

               computation of bessel functions of complex argument
                 by d. e. amos, sand83-0083, may, 1983.

               computation of bessel functions of complex argument
                 and large order by d. e. amos, sand83-0643, may, 1983

               a subroutine package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, sand85-
                 1018, may, 1985

               a portable package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, acm
                 trans. math. software, vol. 12, no. 3, september 1986,
                 pp 265-273.

 ***routines called  zbinu,zabs,i1mach,d1mach
 ***end prologue  zbesi
 */
    /* Initialized data */
    static double coner = 1.;
    static double conei = 0.;

    /* Local variables */
    static int nn, inu;
    static double fn, az, arg, sti,
	zni, str, znr, atol, rtol, ascle, csgni, csgnr;

    /* complex cone,csgn,cw,cy,czero,z,zn */

    /* ***first executable statement  zbesi */
    *ierr = 0;
    *nz = 0;
    if (*fnu < 0.) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*n < 1) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return;
    }
/* -----------------------------------------------------------------------
     set parameters related to machine constants.
     tol is the approximate unit roundoff limited to 1.0e-18.
     elim is the approximate exponential over- and underflow limit.
     exp(-elim) < exp(-alim)=exp(-elim)/tol    and
     exp(elim) > exp(alim)=exp(elim)*tol       are intervals near
     underflow and overflow limits where scaled arithmetic is done.
     rl is the lower boundary of the asymptotic expansion for large z.
     dig = number of base 10 digits in tol = 10**(-dig).
     fnul is the lower boundary of the asymptotic series for large fnu.
     ----------------------------------------------------------------------- */
    double
	tol = fmax2(DBL_EPSILON, 1e-18),
	r1m5 = M_LOG10_2,
	k    = fmin2(-DBL_MIN_EXP, DBL_MAX_EXP),
	elim = (k * r1m5 - 3.) * 2.303,
	aa = r1m5 * (DBL_MANT_DIG - 1),
	dig = fmin2(aa, 18.);
    aa *= 2.303;
    double
	alim = elim + fmax2(-aa, -41.45),
	rl   = dig * 1.2 + 3.,
	fnul = (dig - 3.) * 6. + 10.;
/* -----------------------------------------------------------------------
     test for proper range
 ----------------------------------------------------------------------- */
    az = zabs(*zr, *zi);
    fn = *fnu + (double) ((float) (*n - 1));
    aa = .5 / tol;
    double bb = (double) ((float) INT_MAX) * .5;
    aa = fmin2(aa,bb);
    if (az > aa) {
	goto L260;
    }
    if (fn > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    if (fn > aa) {
	*ierr = 3;
    }
    znr = *zr;
    zni = *zi;
    csgnr = coner;
    csgni = conei;
    if (*zr >= 0.) {
	goto L40;
    }
    znr = -(*zr);
    zni = -(*zi);
/* -----------------------------------------------------------------------
     calculate csgn=exp(fnu*pi*i) to minimize losses of significance
     when fnu is large
 ----------------------------------------------------------------------- */
    inu = (int) ((float) (*fnu));
    arg = (*fnu - (double) ((float) inu)) * pi;
    if (*zi < 0.) {
	arg = -arg;
    }
    csgnr = cos(arg);
    csgni = sin(arg);
    if (inu % 2 == 0) {
	goto L40;
    }
    csgnr = -csgnr;
    csgni = -csgni;
L40:
    *nz = zbinu(&znr, &zni, fnu, *kode, *n, cyr, cyi,
		rl, fnul, tol, elim, alim);
    if (*nz < 0) {
	goto L120;
    }
    if (*zr >= 0.) {
	return;
    }
/* -----------------------------------------------------------------------
     analytic continuation to the left half plane
 ----------------------------------------------------------------------- */
    nn = *n - *nz;
    if (nn == 0) {
	return;
    }
    rtol = 1. / tol;
    ascle = DBL_MIN * rtol * 1e3;
    for (int i = 0; i < nn; ++i) {
	/*
	  str = cyr(i)*csgnr - cyi(i)*csgni
	  cyi(i) = cyr(i)*csgni + cyi(i)*csgnr
	  cyr(i) = str */
	aa = cyr[i];
	bb = cyi[i];
	atol = 1.;
	if (fmax2(fabs(aa), fabs(bb)) <= ascle) {
	    aa *= rtol;
	    bb *= rtol;
	    atol = tol;
	}
	str = aa * csgnr - bb * csgni;
	sti = aa * csgni + bb * csgnr;
	cyr[i] = str * atol;
	cyi[i] = sti * atol;
	csgnr = -csgnr;
	csgni = -csgni;
    }
    return;
L120:
    if (*nz == -2) goto L130;
         *nz = 0;    *ierr = 2;    return;
L130:    *nz = 0;    *ierr = 5;    return;
L260:    *nz = 0;    *ierr = 4;    return;
} /* zbesi */

void zbesj(double *zr, double *zi, double *fnu, int *kode, int *n,
	   double *cyr, double *cyi, int *nz, int *ierr)
{
/***begin prologue  zbesj
 ***date written   830501   (yymmdd)
 ***revision date  890801, 930101   (yymmdd)
 ***category no.  b5k
 ***keywords  j-bessel function,bessel function of complex argument,
             bessel function of first kind
 ***author  amos, donald e., sandia national laboratories
 ***purpose  to compute the j-bessel function of a complex argument
 ***description

                      ***a double precision routine***
         on kode=1, zbesj computes an n member  sequence of complex
         bessel functions cy(j) = J(fnu+j-1,z) for float, nonnegative
         orders fnu+j-1, j=1,...,n and complex z in the cut plane
         -pi < arg(z) <= pi. on kode=2, zbesj returns the scaled
         functions

         cy(j)=exp(-abs(y))*J(fnu+j-1,z)   j = 1,...,n , y=aimag(z)

         which remove the exponential growth in both the upper and
         lower half planes for z to infinity. definitions and notation
         are found in the nbs handbook of mathematical functions
         (ref. 1).

         input      zr,zi,fnu are double precision
           zr,zi  - z=cmplx(zr,zi),  -pi < arg(z) <= pi
           fnu    - order of initial j function, fnu >= 0
           kode   - a parameter to indicate the scaling option
                    kode= 1  returns
                             cy(j)=J(fnu+j-1,z), j=1,...,n
                        = 2  returns
                             cy(j)=J(fnu+j-1,z)exp(-abs(y)), j=1,...,n
           n      - number of members of the sequence, n >= 1

         output     cyr,cyi are double precision
           cyr,cyi- double precision vectors whose first n components
                    contain float and imaginary parts for the sequence
                    cy(j)=J(fnu+j-1,z)  or
                    cy(j)=J(fnu+j-1,z)exp(-abs(y))  j=1,...,n
                    depending on kode, y=aimag(z).
           nz     - number of components set to zero due to underflow,
                    nz= 0   , normal return
                    nz > 0 , last nz components of cy set  zero due
                              to underflow, cy(j)=cmplx(0,0),
                              j = n-nz+1,...,n
           ierr   - error flag
                    ierr=0, normal return - computation completed
                    ierr=1, input error   - no computation
                    ierr=2, overflow      - no computation, aimag(z)
                            too large on kode=1
                    ierr=3, cabs(z) or fnu+n-1 large - computation done
                            but losses of signifcance by argument
                            reduction produce less than half of machine
                            accuracy
                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
                            tion because of complete losses of signifi-
                            cance by argument reduction
                    ierr=5, error              - no computation,
                            algorithm termination condition not met

 ***long description

         the computation is carried out by the formula

         J(fnu,z)=exp( fnu*pi*i/2)*I(fnu,-i*z)    aimag(z) >= 0.0

         J(fnu,z)=exp(-fnu*pi*i/2)*I(fnu, i*z)    aimag(z) < 0.0

         where i**2 = -1 and I(fnu,z) is the i bessel function.

         for negative orders,the formula

              J(-fnu,z) = J(fnu,z)*cos(pi*fnu) - Y(fnu,z)*sin(pi*fnu)

         can be used. however,for large orders close to integers, the
         the function changes radically. when fnu is a large positive
         int,the magnitude of J(-fnu,z)=J(fnu,z)*cos(pi*fnu) is a
         large negative power of ten. but when fnu is not an int,
         Y(fnu,z) dominates in magnitude with a large positive power of
         ten and the most that the second term can be reduced is by
         unit roundoff from the coefficient. thus, wide changes can
         occur within unit roundoff of a large int for fnu. here,
         large means fnu > cabs(z).

         in most complex variable computation, one must evaluate ele-
         mentary functions. when the magnitude of z or fnu+n-1 is
         large, losses of significance by argument reduction occur.
         consequently, if either one exceeds u1=sqrt(0.5/ur), then
         losses exceeding half precision are likely and an error flag
         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
         double precision unit roundoff limited to 18 digits precision.
         if either is larger than u2=0.5/ur, then all significance is
         lost and ierr=4. in order to use the int function, arguments
         must be further restricted not to exceed the largest machine
         int, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double precision
         arithmetic respectively. this makes u2 and u3 limiting in
         their respective arithmetics. this means that one can expect
         to retain, in the worst cases on 32 bit machines, no digits
         in single and only 7 digits in double precision arithmetic.
         similar considerations hold for other machines.

         the approximate relative error in the magnitude of a complex
         bessel function can be expressed by p*10**s where p=max(unit
         roundoff,1.0e-18) is the nominal precision and 10**s repre-
         sents the increase in error due to argument reduction in the
         elementary functions. here, s=max(1,abs(log10(cabs(z))),
         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
         have only absolute accuracy. this is most likely to occur when
         one component (in absolute value) is larger than the other by
         several orders of magnitude. if one component is 10**k larger
         than the other, then one can expect only max(abs(log10(p))-k,
         0) significant digits; or, stated another way, when k exceeds
         the exponent of p, no significant digits remain in the smaller
         component. however, the phase angle retains absolute accuracy
         because, in complex arithmetic with precision p, the smaller
         component will not (as a rule) decrease below p times the
         magnitude of the larger component. in these extreme cases,
         the principal phase angle is on the order of +p, -p, pi/2-p,
         or -pi/2+p.

 ***references  handbook of mathematical functions by m. abramowitz
                 and i. a. stegun, nbs ams series 55, u.s. dept. of
                 commerce, 1955.

               computation of bessel functions of complex argument
                 by d. e. amos, sand83-0083, may, 1983.

               computation of bessel functions of complex argument
                 and large order by d. e. amos, sand83-0643, may, 1983

               a subroutine package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, sand85-
                 1018, may, 1985

               a portable package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, acm
                 trans. math. software, vol. 12, no. 3, september 1986,
                 pp 265-273.

 ***routines called  zbinu, zabs, i1mach, d1mach
 ***end prologue  zbesj */

    *nz = 0;
    if (*fnu < 0. || *kode < 1 || *kode > 2 || *n < 1) {
	*ierr = 1;	return;
    }
    *ierr = 0;

/* -----------------------------------------------------------------------
     set parameters related to machine constants.
     tol is the approximate unit roundoff limited to 1.0e-18.
     elim is the approximate exponential over- and underflow limit.
     exp(-elim) < exp(-alim)=exp(-elim)/tol    and
     exp(elim) > exp(alim)=exp(elim)*tol       are intervals near
     underflow and overflow limits where scaled arithmetic is done.
     rl is the lower boundary of the asymptotic expansion for large z.
     dig = number of base 10 digits in tol = 10**(-dig).
     fnul is the lower boundary of the asymptotic series for large fnu.
     ----------------------------------------------------------------------- */
    double
	tol = fmax2(DBL_EPSILON, 1e-18),
	r1m5 = M_LOG10_2,
	k    = fmin2(-DBL_MIN_EXP, DBL_MAX_EXP),
	elim = (k * r1m5 - 3.) * 2.303,
	aa = r1m5 * (DBL_MANT_DIG - 1),
	dig = fmin2(aa, 18.);
    aa *= 2.303;
    double
	alim = elim + fmax2(-aa, -41.45),
	rl   = dig * 1.2 + 3.,
	fnul = (dig - 3.) * 6. + 10.;
/* -----------------------------------------------------------------------
     test for proper range
 ----------------------------------------------------------------------- */
    double az = zabs(*zr, *zi),
	fn = *fnu + (double) ((float) (*n - 1)),
	bb = (double) ((float) INT_MAX) * .5;
    aa = .5 / tol;
    aa = fmin2(aa,bb);
    if (az > aa || fn > aa) { // goto L260 -- L260:
	*nz = 0;	*ierr = 4;	return;
    }
    // else
    aa = sqrt(aa);
    if (az > aa || fn > aa) {
	*ierr = 3;
    }
/* -----------------------------------------------------------------------
     calculate csgn=exp(fnu*hpi*i) to minimize losses of significance
     when fnu is large
 ----------------------------------------------------------------------- */
    int inu = (int) ((float) (*fnu)),
	inuh = inu / 2,
	ir = inu - (inuh << 1);
    double cii = 1.,
	arg = (*fnu - (double) ((float) (inu - ir))) * hpi,
	csgnr = cos(arg),
	csgni = sin(arg);
    if (inuh % 2) { // L40:
	csgnr = -csgnr;
	csgni = -csgni;
    }
/* -----------------------------------------------------------------------
     zn is in the right half plane
 ----------------------------------------------------------------------- */
    double
	znr = *zi,
	zni = -(*zr);
    if (*zi >= 0.) {
	goto L50;
    }
    znr = -znr;
    zni = -zni;
    csgni = -csgni;
    cii = -cii;
L50:
    *nz = zbinu(&znr, &zni, fnu, *kode, *n, cyr, cyi,
		rl, fnul, tol, elim, alim);
    if (*nz >= 0) {
	int nl = *n - *nz;
	if (nl == 0)
	    return;
	// else :
	double rtol = 1. / tol,
	    ascle = DBL_MIN * rtol * 1e3;
	for (int i = 0; i < nl; ++i) {
	    /*
	      str = cyr(i)*csgnr - cyi(i)*csgni
	      cyi(i) = cyr(i)*csgni + cyi(i)*csgnr
	      cyr(i) = str */
	    aa = cyr[i];
	    bb = cyi[i];
	    double atol = 1.;
	    if (fmax2(fabs(aa), fabs(bb)) <= ascle) {
		aa *= rtol;
		bb *= rtol;
		atol = tol;
	    }
	    double
		str = aa * csgnr - bb * csgni,
		sti = aa * csgni + bb * csgnr;
	    cyr[i] = str * atol;
	    cyi[i] = sti * atol;
	    str = -csgni * cii;
	    csgni = csgnr * cii;
	    csgnr = str;
	}
    } else if (*nz == -2) {
	*nz = 0;    *ierr = 5;
    } else { // L140
	*nz = 0;    *ierr = 2;
    }
    return;

} /* zbesj ------------------------------------------------------------------------*/


void zbesk(double *zr, double *zi, double *fnu, int *kode, int *n,
	   double *cyr, double *cyi, int *nz, int *ierr)
{
/***begin prologue  zbesk
 ***date written   830501   (yymmdd)
 ***revision date  890801, 930101   (yymmdd)
 ***category no.  b5k
 ***keywords  k-bessel function,complex bessel function,
             modified bessel function of the second kind,
             bessel function of the third kind
 ***author  amos, donald e., sandia national laboratories
 ***purpose  to compute k-bessel functions of complex argument
 ***description

                      ***a double precision routine***

         on kode=1, zbesk computes an n member sequence of complex
         bessel functions cy(j)=K(fnu+j-1,z) for float, nonnegative
         orders fnu+j-1, j=1,...,n and complex z != cmplx(0.0,0.0)
         in the cut plane -pi < arg(z) <= pi. on kode=2, zbesk
         returns the scaled k functions,

         cy(j)=exp(z)*K(fnu+j-1,z) , j=1,...,n,

         which remove the exponential behavior in both the left and
         right half planes for z to infinity. definitions and
         notation are found in the nbs handbook of mathematical
         functions (ref. 1).

         input      zr,zi,fnu are double precision
           zr,zi  - z=cmplx(zr,zi), z != cmplx(0,0),
                    -pi < arg(z) <= pi
           fnu    - order of initial k function, fnu >= 0
           n      - number of members of the sequence, n >= 1
           kode   - a parameter to indicate the scaling option
                    kode= 1  returns
                             cy(j)=K(fnu+j-1,z), j=1,...,n
                        = 2  returns
                             cy(j)=K(fnu+j-1,z)*exp(z), j=1,...,n

         output     cyr,cyi are double precision
           cyr,cyi- double precision vectors whose first n components
                    contain float and imaginary parts for the sequence
                    cy(j)=K(fnu+j-1,z), j=1,...,n or
                    cy(j)=K(fnu+j-1,z)*exp(z), j=1,...,n
                    depending on kode
           nz     - number of components set to zero due to underflow.
                    nz= 0   , normal return
                    nz > 0 , first nz components of cy set to zero due
                              to underflow, cy(j)=cmplx(0,0),
                              j=1,...,n when x >= 0.0. when x < 0.0
                              nz states only the number of underflows
                              in the sequence.

           ierr   - error flag
                    ierr=0, normal return - computation completed
                    ierr=1, input error   - no computation
                    ierr=2, overflow      - no computation, fnu is
                            too large or cabs(z) is too small or both
                    ierr=3, cabs(z) or fnu+n-1 large - computation done
                            but losses of signifcance by argument
                            reduction produce less than half of machine
                            accuracy
                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
                            tion because of complete losses of signifi-
                            cance by argument reduction
                    ierr=5, error              - no computation,
                            algorithm termination condition not met

 ***long description

         equations of the reference are implemented for small orders
         dnu and dnu+1.0 in the right half plane x >= 0.0. forward
         recurrence generates higher orders. k is continued to the left
         half plane by the relation

         K(fnu,z*exp(mp)) = exp(-mp*fnu)*K(fnu,z)-mp*I(fnu,z)
         mp=mr*pi*i, mr=+1 or -1, re(z) > 0, i**2=-1

         where I(fnu,z) is the j bessel function.

         for large orders, fnu > fnul, the k function is computed
         by means of its uniform asymptotic expansions.

         for negative orders, the formula

                       K(-fnu,z) = K(fnu,z)

         can be used.

         zbesk assumes that a significant digit sinh(x) function is
         available.

         in most complex variable computation, one must evaluate ele-
         mentary functions. when the magnitude of z or fnu+n-1 is
         large, losses of significance by argument reduction occur.
         consequently, if either one exceeds u1=sqrt(0.5/ur), then
         losses exceeding half precision are likely and an error flag
         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
         double precision unit roundoff limited to 18 digits precision.
         if either is larger than u2=0.5/ur, then all significance is
         lost and ierr=4. in order to use the int function, arguments
         must be further restricted not to exceed the largest machine
         int, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double precision
         arithmetic respectively. this makes u2 and u3 limiting in
         their respective arithmetics. this means that one can expect
         to retain, in the worst cases on 32 bit machines, no digits
         in single and only 7 digits in double precision arithmetic.
         similar considerations hold for other machines.

         the approximate relative error in the magnitude of a complex
         bessel function can be expressed by p*10**s where p=max(unit
         roundoff,1.0e-18) is the nominal precision and 10**s repre-
         sents the increase in error due to argument reduction in the
         elementary functions. here, s=max(1,abs(log10(cabs(z))),
         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
         have only absolute accuracy. this is most likely to occur when
         one component (in absolute value) is larger than the other by
         several orders of magnitude. if one component is 10**k larger
         than the other, then one can expect only max(abs(log10(p))-k,
         0) significant digits; or, stated another way, when k exceeds
         the exponent of p, no significant digits remain in the smaller
         component. however, the phase angle retains absolute accuracy
         because, in complex arithmetic with precision p, the smaller
         component will not (as a rule) decrease below p times the
         magnitude of the larger component. in these extreme cases,
         the principal phase angle is on the order of +p, -p, pi/2-p,
         or -pi/2+p.

 ***references  handbook of mathematical functions by m. abramowitz
                 and i. a. stegun, nbs ams series 55, u.s. dept. of
                 commerce, 1955.

               computation of bessel functions of complex argument
                 by d. e. amos, sand83-0083, may, 1983.

               computation of bessel functions of complex argument
                 and large order by d. e. amos, sand83-0643, may, 1983.

               a subroutine package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, sand85-
                 1018, may, 1985

               a portable package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, acm
                 trans. math. software, vol. 12, no. 3, september 1986,
                 pp 265-273.

 ***routines called  zacon,zbknu,zbunk,zuoik,zabs,i1mach,d1mach
 ***end prologue  zbesk
 */

    int verbose = *ierr, nn = *n;

    /* Function Body */
    *nz = 0;
    if(*fnu < 0. || *kode < 1 || *kode > 2 || *n < 1) {
	*ierr = 1; return;
    }
    else
	*ierr = 0;
    if(*zi == 0. && *zr == 0.) { // K_nu(0) := 0 [ok for real axis]
	cyr[0] = cyi[0] = 0.;
	return;
    }


/* -----------------------------------------------------------------------
     set parameters related to machine constants.
     tol is the approximate unit roundoff limited to 1.0e-18.
     elim is the approximate exponential over- and underflow limit.
     exp(-elim) < exp(-alim)=exp(-elim)/tol    and
     exp(elim) > exp(alim)=exp(elim)*tol       are intervals near
     underflow and overflow limits where scaled arithmetic is done.
     rl is the lower boundary of the asymptotic expansion for large z.
     dig = number of base 10 digits in tol = 10**(-dig).
     fnul is the lower boundary of the asymptotic series for large fnu
 ----------------------------------------------------------------------- */
    double
	tol = fmax2(DBL_EPSILON, 1e-18),        // = DBL_EPSILON = 2.22..e-16
	r1m5 = M_LOG10_2,
	k    = fmin2(-DBL_MIN_EXP, DBL_MAX_EXP),// = 1021
	elim = (k * r1m5 - 3.) * 2.303,         // = 700.92179369444591
	aa = r1m5 * (DBL_MANT_DIG - 1),
	dig = fmin2(aa, 18.);                   // = 15.653559774527022
    aa *= 2.303;
    double
	alim = elim + fmax2(-aa, -41.45),	// = 664.87164553371019
	fnul = (dig - 3.) * 6. + 10.,		// = 85.921358647162123
	rl = dig * 1.2 + 3.;			// = 21.784271729432426
/* -----------------------------------------------------------------------------
     test for proper range
 ----------------------------------------------------------------------- */
    double
	az = zabs(*zr, *zi),
	fn = *fnu + (double) (nn - 1),
	bb = .5 * (double) INT_MAX;
    aa = fmin2(.5 / tol, bb);			// = 1073741823.5
    if (az > aa || fn > aa) { // L260:
	*nz = 0;	*ierr = 4;	return;
    }
    aa = sqrt(aa);				// = 32767.999992..
    if (az > aa || fn > aa) {
	*ierr = 3;
	if(verbose) Rprintf("zbesk():  |z| or nu+n-1 too large --> ierr=3\n");
    }
/* -----------------------------------------------------------------------
     overflow test on the last member of the sequence
   ----------------------------------------------------------------------- */
    if (az < ufl) {
	goto L_err_2;
    }
    if (*fnu > fnul) { // nu > 85.92 ..
	goto L80;
    }
    if (fn > 2.) { // L50:
	int nuf = zuoik(zr, zi, fnu, *kode, 2, nn, cyr, cyi, tol, elim, alim);
	if(verbose) Rprintf(" after zuoik(): nuf=%d, cy(r,i)[0]= (%g,%g)\n",
			    nuf, cyr[0], cyi[0]);
	if (nuf < 0) {
	    goto L_err_2;
	}
	*nz += nuf;
	nn -= nuf;
	/* -----------------------------------------------------------------------
	   here nn=n or nn=0 since nuf=0,nn, or -1 on return from cuoik
	   if nuf=nn, then cy(i)=czero for all i
	   ----------------------------------------------------------------------- */
	if (nn == 0)
	    goto L100;
    }
    else if (fn > 1 && az <= tol) {
	double aln = -fn * log(.5 * az);
	if (aln > elim) {
	    goto L_err_2;
	}
    }

    // L60:
    int nw;
    if (*zr >= 0.) {
	/* -----------------------------------------------------------------------
	   right half plane computation, float(z) >= 0.
	   ----------------------------------------------------------------------- */
	nw = zbknu(zr, zi, fnu, *kode, nn, verbose,
		   cyr, cyi, tol, elim, alim);
	if(verbose) Rprintf(" after zbknu(): nw=%d\n", nw);
	if (nw < 0) {
	    goto L200;
	}
	*nz = nw;
	return;
    }
    // L70:
    /* -----------------------------------------------------------------------
       left half plane computation
       pi/2 < arg(z) <= pi and -pi < arg(z) < -pi/2.
	   ----------------------------------------------------------------------- */
    if (*nz != 0) {
	goto L_err_2;
    }
    int mr = 1;
    if (*zi < 0.) {
	mr = -1;
    }
    zacon_(zr, zi, fnu, kode, &mr, &nn, cyr, cyi, &nw,
	   verbose, rl, fnul, tol, elim, alim);
    if(verbose) Rprintf(" after zacon(): nw=%d\n", nw);
    if (nw < 0) {
	goto L200;
    }
    *nz = nw;
    return;

L80:
/* -----------------------------------------------------------------------
     uniform asymptotic expansions for fnu > fnul
 ----------------------------------------------------------------------- */
    if (*zr >= 0.)
	mr = 0;
    else if (*zi >= 0.)
	mr = 1;
    else // (*zi < 0.)
	mr = -1;
    if(verbose) Rprintf(" fnu > fnul=%g: mr=%d;", fnul, mr);
    // L90:
    zbunk_(zr, zi, fnu, kode, &mr, &nn, cyr, cyi, &nw, &tol, &elim, &alim);
    if(verbose) Rprintf("  after zbunk(): nw=%d\n", nw);
    if (nw < 0) {
	goto L200;
    }
    *nz += nw;
    return;

L100:
    if (*zr < 0.) goto L_err_2;
    return;

L200:
    if (nw != -1) {
	*nz = 0;	*ierr = 5;	return;
    }
L_err_2: *nz = 0;	*ierr = 2;	return;

} /* zbesk */


void zbesy(double *zr, double *zi, double *fnu,
	   int *kode, int *n, double *cyr, double *cyi,
	   int *nz, double *cwrkr, double *cwrki, int *ierr)
{
/***begin prologue  zbesy
 ***date written   830501   (yymmdd)
 ***revision date  890801, 930101   (yymmdd)
 ***category no.  b5k
 ***keywords  y-bessel function,bessel function of complex argument,
             bessel function of second kind
 ***author  amos, donald e., sandia national laboratories
 ***purpose  to compute the y-bessel function of a complex argument
 ***description

                      ***a double precision routine***

         on kode=1, zbesy computes an n member sequence of complex
         bessel functions cy(j)=Y(fnu+j-1,z) for float, nonnegative
         orders fnu+j-1, j=1,...,n and complex z in the cut plane
         -pi < arg(z) <= pi. on kode=2, zbesy returns the scaled
         functions

         cy(j)=exp(-abs(y))*Y(fnu+j-1,z)   j = 1,...,n , y=aimag(z)

         which remove the exponential growth in both the upper and
         lower half planes for z to infinity. definitions and notation
         are found in the nbs handbook of mathematical functions
         (ref. 1).

         input      zr,zi,fnu are double precision
           zr,zi  - z=cmplx(zr,zi), z != cmplx(0,0),
                    -pi < arg(z) <= pi
           fnu    - order of initial y function, fnu >= 0
           kode   - a parameter to indicate the scaling option
                    kode= 1  returns cy(j) = Y(fnu+j-1,z),              j=1,...,n
                        = 2  returns cy(j) = Y(fnu+j-1,z)*exp(-abs(y)), j=1,...,n
                                             where y = aimag(z)
           n      - number of members of the sequence, n >= 1
           cwrkr, - double precision work vectors of dimension at
           cwrki    at least n

         output     cyr,cyi are double precision
           cyr,cyi- double precision vectors whose first n components
                    contain float and imaginary parts for the sequence
                    cy(j), j = 1..n, see above

           nz     - nz = 0 , a normal return
                    nz > 0 , nz components of cy set to zero due to
		             underflow (generally on kode=2)
           ierr   - error flag
                    ierr=0, normal return - computation completed
                    ierr=1, input error   - no computation
                    ierr=2, overflow      - no computation, fnu is
                            too large or cabs(z) is too small or both
                    ierr=3, cabs(z) or fnu+n-1 large - computation done
                            but losses of signifcance by argument
                            reduction produce less than half of machine
                            accuracy
                    ierr=4, cabs(z) or fnu+n-1 too large - no computa-
                            tion because of complete losses of signifi-
                            cance by argument reduction
                    ierr=5, error              - no computation,
                            algorithm termination condition not met

 ***long description

         the computation is carried out in terms of the I(fnu,z) and
         K(fnu,z) bessel functions in the right half plane by

             Y(fnu,z) = i*cc*I(fnu,arg) - (2/pi)*conjg(cc)*K(fnu,arg)

             Y(fnu,z) = conjg(Y(fnu,conjg(z)))

         for aimag(z) >= 0 and aimag(z) < 0 respectively, where
         cc=exp(i*pi*fnu/2), arg=z*exp(-i*pi/2) and i**2=-1.

         for negative orders,the formula

              Y(-fnu,z) = Y(fnu,z)*cos(pi*fnu) + J(fnu,z)*sin(pi*fnu)

         can be used. however,for large orders close to half odd
         integers the function changes radically. when fnu is a large
         positive half odd int,the magnitude of Y(-fnu,z)=J(fnu,z)*
         sin(pi*fnu) is a large negative power of ten. but when fnu is
         not a half odd int, Y(fnu,z) dominates in magnitude with a
         large positive power of ten and the most that the second term
         can be reduced is by unit roundoff from the coefficient. thus,
         wide changes can occur within unit roundoff of a large half
         odd int. here, large means fnu > cabs(z).

         in most complex variable computation, one must evaluate ele-
         mentary functions. when the magnitude of z or fnu+n-1 is
         large, losses of significance by argument reduction occur.
         consequently, if either one exceeds u1=sqrt(0.5/ur), then
         losses exceeding half precision are likely and an error flag
         ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
         double precision unit roundoff limited to 18 digits precision.
         if either is larger than u2=0.5/ur, then all significance is
         lost and ierr=4. in order to use the int function, arguments
         must be further restricted not to exceed the largest machine
         int, u3=i1mach(9). thus, the magnitude of z and fnu+n-1 is
         restricted by min(u2,u3). on 32 bit machines, u1,u2, and u3
         are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single precision
         arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double precision
         arithmetic respectively. this makes u2 and u3 limiting in
         their respective arithmetics. this means that one can expect
         to retain, in the worst cases on 32 bit machines, no digits
         in single and only 7 digits in double precision arithmetic.
         similar considerations hold for other machines.

         the approximate relative error in the magnitude of a complex
         bessel function can be expressed by p*10**s where p=max(unit
         roundoff,1.0e-18) is the nominal precision and 10**s repre-
         sents the increase in error due to argument reduction in the
         elementary functions. here, s=max(1,abs(log10(cabs(z))),
         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
         have only absolute accuracy. this is most likely to occur when
         one component (in absolute value) is larger than the other by
         several orders of magnitude. if one component is 10**k larger
         than the other, then one can expect only max(abs(log10(p))-k,
         0) significant digits; or, stated another way, when k exceeds
         the exponent of p, no significant digits remain in the smaller
         component. however, the phase angle retains absolute accuracy
         because, in complex arithmetic with precision p, the smaller
         component will not (as a rule) decrease below p times the
         magnitude of the larger component. in these extreme cases,
         the principal phase angle is on the order of +p, -p, pi/2-p,
         or -pi/2+p.

 ***references  handbook of mathematical functions by m. abramowitz
                 and i. a. stegun, nbs ams series 55, u.s. dept. of
                 commerce, 1955.

               computation of bessel functions of complex argument
                 by d. e. amos, sand83-0083, may, 1983.

               computation of bessel functions of complex argument
                 and large order by d. e. amos, sand83-0643, may, 1983

               a subroutine package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, sand85-
                 1018, may, 1985

               a portable package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, acm
                 trans. math. software, vol. 12, no. 3, september 1986,
                 pp 265-273.

 ***routines called  zbesi,zbesk,i1mach,d1mach
 ***end prologue  zbesy
 */

    /* Initialized data */
    static double cipr[4] = { 1.,0.,-1.,0. };
    static double cipi[4] = { 0.,1.,0.,-1. };

    /* Local variables */
    int i, nz1, nz2, ifnu;
    double arg, exi, exr, sti, zni, zui, str, znr,
	zvi, zzi, zur, zvr, zzr, ffnu, atol, rhpi,
	rtol, ascle, csgni, csgnr, cspni, cspnr;

    /* complex cwrk,cy,c1,c2,ex,hci,z,zu,zv
     */

    int verbose = *ierr;

    *nz = 0;
    if (*fnu < 0. || *kode < 1 || *kode > 2 || *n < 1) {
	*ierr = 1; return;
    }
    else
	*ierr = 0;
    if(*zi == 0. && *zr == 0.) { // Y_nu(0) := 0 [ok for real axis]
	*cyr = *cyi = 0.; return;
    }

    zzr = *zr;
    zzi = *zi;
    if (*zi < 0.) {
	zzi = -zzi;
    }
    znr = zzi;
    zni = -zzr;
    *ierr = verbose; // verbose for zbesi()
    zbesi(&znr, &zni, fnu, kode, n, cyr, cyi, &nz1, ierr);
    if (*ierr != 0 && *ierr != 3) {
	/* L90: */ *nz = 0; return;
    }
    *ierr = verbose;
    zbesk(&znr, &zni, fnu, kode, n, cwrkr, cwrki, &nz2, ierr);
    if (*ierr != 0 && *ierr != 3) {
	/* L90: */ *nz = 0; return;
    }
    *nz = imin2(nz1,nz2);
    ifnu = (int) *fnu;
    ffnu = *fnu - (double) ifnu;
    arg = hpi * ffnu;
    csgnr = cos(arg);
    csgni = sin(arg);
    int i4 = ifnu % 4;
    str   = csgnr * cipr[i4] - csgni * cipi[i4];
    csgni = csgnr * cipi[i4] + csgni * cipr[i4];
    csgnr = str;
    rhpi = 1. / hpi;
    cspnr = csgnr * rhpi;
    cspni = -csgni * rhpi;
    str = -csgni;
    csgni = csgnr;
    csgnr = str;
    if (*kode == 2) {
	goto L60;
    }
    for (i = 0; i < *n; ++i) {
	/* cy(i) = csgn*cy(i)-cspn*cwrk(i) */
	str = csgnr * cyr[i] - csgni * cyi[i];
	str -= cspnr * cwrkr[i] - cspni * cwrki[i];
	sti = csgnr * cyi[i] + csgni * cyr[i];
	sti -= cspnr * cwrki[i] + cspni * cwrkr[i];
	cyr[i] = str;
	cyi[i] = sti;
	str = -csgni; csgni =  csgnr; csgnr = str;
	str =  cspni; cspni = -cspnr; cspnr = str;
    }
    if (*zi < 0.) {
	for (int i = 0; i < *n; ++i) {
	    cyi[i] = -cyi[i];
	}
    }
    return;

L60:
    exr = cos(*zr);
    exi = sin(*zr);
    double
	tol  = fmax2(DBL_EPSILON, 1e-18),
	d1m5 = M_LOG10_2,
	k    = fmin2(-DBL_MIN_EXP, DBL_MAX_EXP),
	/* -----------------------------------------------------------------------
	   elim is the approximate exponential under- and overflow limit
	   ----------------------------------------------------------------------- */
	elim = (double) (k * d1m5 - 3.) * 2.303,
	ey   = 0.,
	tay  = fabs(*zi + *zi);
    if (tay < elim) {
	ey = exp(-tay);
    }
    str   = (exr * cspnr - exi * cspni) * ey;
    cspni = (exr * cspni + exi * cspnr) * ey;
    cspnr = str;
    *nz = 0;
    rtol = 1. / tol;
    ascle = DBL_MIN * rtol * 1e3;
    for (int i = 0; i < *n; ++i) {
	/*----------------------------------------------------------------------
	  cy(i) = csgn*cy(i)-cspn*cwrk(i): products are computed in
	  scaled mode if cy(i) or cwrk(i) are close to underflow to
	  prevent underflow in an intermediate computation.
	  ---------------------------------------------------------------------- */
	zvr = cwrkr[i];
	zvi = cwrki[i];
	atol = 1.;
	if (fmax2(fabs(zvr), fabs(zvi)) <= ascle) {
	    zvr *= rtol;
	    zvi *= rtol;
	    atol = tol;
	}
	str = (zvr * cspnr - zvi * cspni) * atol;
	zvi = (zvr * cspni + zvi * cspnr) * atol;
	zvr = str;
	zur = cyr[i];
	zui = cyi[i];
	atol = 1.;
	if (fmax2(fabs(zur), fabs(zui)) <= ascle) {
	    zur *= rtol;
	    zui *= rtol;
	    atol = tol;
	}
	str = (zur * csgnr - zui * csgni) * atol;
	zui = (zur * csgni + zui * csgnr) * atol;
	zur = str;
	cyr[i] = zur - zvr;
	cyi[i] = zui - zvi;
	if (*zi < 0.) {
	    cyi[i] = -cyi[i];
	}
	if (cyr[i] == 0. && cyi[i] == 0. && ey == 0.) {
	    ++(*nz);
	}
	str = -csgni; csgni =  csgnr; csgnr = str;
	str =  cspni; cspni = -cspnr; cspnr = str;
    }
    return;

} /* zbesy */

void zairy(double *zr, double *zi, int *id, int *kode,
	   double *air, double *aii, int *nz, int *ierr)
{
/***begin prologue  zairy
 ***date written   830501   (yymmdd)
 ***revision date  890801, 930101   (yymmdd)
 ***category no.  b5k
 ***keywords  airy function,bessel functions of order one third
 ***author  amos, donald e., sandia national laboratories
 ***purpose  to compute airy functions ai(z) and dai(z) for complex z
 ***description

 ***a double precision routine***
 on kode=1, zairy computes the complex airy function ai(z) or
 its derivative dai(z)/dz on id=0 or id=1 respectively. on
 kode=2, a scaling option cexp(zta)*ai(z) or cexp(zta)*
 dai(z)/dz is provided to remove the exponential decay in
 -pi/3 < arg(z) < pi/3 and the exponential growth in
 pi/3 < abs(arg(z)) < pi where zta=(2/3)*z*csqrt(z).

 while the airy functions ai(z) and dai(z)/dz are analytic in
 the whole z plane, the corresponding scaled functions defined
 for kode=2 have a cut along the negative float axis.
 definitions and notation are found in the nbs handbook of
 mathematical functions (ref. 1).

 input      zr,zi are double precision
 zr,zi  - z=cmplx(zr,zi)
 id     - order of derivative, id=0 or id=1
 kode   - a parameter to indicate the scaling option
 kode= 1  returns
 ai=ai(z)                on id=0 or
 ai=dai(z)/dz            on id=1
 = 2  returns
 ai=cexp(zta)*ai(z)       on id=0 or
 ai=cexp(zta)*dai(z)/dz   on id=1 where
 zta=(2/3)*z*csqrt(z)

 output     air,aii are double precision
 air,aii- complex answer depending on the choices for id and
 kode
 nz     - underflow indicator
 nz= 0   , normal return
 nz= 1   , ai=cmplx(0,0) due to underflow in
 -pi/3 < arg(z) < pi/3 on kode=1
 ierr   - error flag
 ierr=0, normal return - computation completed
 ierr=1, input error   - no computation
 ierr=2, overflow      - no computation, float(zta)
 too large on kode=1
 ierr=3, cabs(z) large      - computation completed
 losses of signifcance by argument reduction
 produce less than half of machine accuracy
 ierr=4, cabs(z) too large  - no computation
 complete loss of accuracy by argument
 reduction
 ierr=5, error              - no computation,
 algorithm termination condition not met

 ***long description

 ai and dai are computed for cabs(z) > 1.0 from the k bessel
 functions by

 ai(z)=c*sqrt(z)*K(1/3,zta) , dai(z)=-c*z*K(2/3,zta)
 c=1.0/(pi*sqrt(3.0))
 zta=(2/3)*z**(3/2)

 with the power series for cabs(z) <= 1.0.

 in most complex variable computation, one must evaluate ele-
 mentary functions. when the magnitude of z is large, losses
 of significance by argument reduction occur. consequently, if
 the magnitude of zeta=(2/3)*z**1.5 exceeds u1=sqrt(0.5/ur),
 then losses exceeding half precision are likely and an error
 flag ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
 double precision unit roundoff limited to 18 digits precision.
 also, if the magnitude of zeta is larger than u2=0.5/ur, then
 all significance is lost and ierr=4. in order to use the int
 function, zeta must be further restricted not to exceed the
 largest int, u3=i1mach(9). thus, the magnitude of zeta
 must be restricted by min(u2,u3). on 32 bit machines, u1,u2,
 and u3 are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single
 precision arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double
 precision arithmetic respectively. this makes u2 and u3 limit-
 ing in their respective arithmetics. this means that the mag-
 nitude of z cannot exceed 3.1e+4 in single and 2.1e+6 in
 double precision arithmetic. this also means that one can
 expect to retain, in the worst cases on 32 bit machines,
 no digits in single precision and only 7 digits in double
 precision arithmetic. similar considerations hold for other
 machines.

 the approximate relative error in the magnitude of a complex
 bessel function can be expressed by p*10**s where p=max(unit
 roundoff,1.0e-18) is the nominal precision and 10**s repre-
 sents the increase in error due to argument reduction in the
 elementary functions. here, s=max(1,abs(log10(cabs(z))),
 abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
 cabs(z),abs(exponent of fnu)) ). however, the phase angle may
 have only absolute accuracy. this is most likely to occur when
 one component (in absolute value) is larger than the other by
 several orders of magnitude. if one component is 10**k larger
 than the other, then one can expect only max(abs(log10(p))-k,
 0) significant digits; or, stated another way, when k exceeds
 the exponent of p, no significant digits remain in the smaller
 component. however, the phase angle retains absolute accuracy
 because, in complex arithmetic with precision p, the smaller
 component will not (as a rule) decrease below p times the
 magnitude of the larger component. in these extreme cases,
 the principal phase angle is on the order of +p, -p, pi/2-p,
 or -pi/2+p.

 ***references  handbook of mathematical functions by m. abramowitz
 and i. a. stegun, nbs ams series 55, u.s. dept. of
 commerce, 1955.

 computation of bessel functions of complex argument
 and large order by d. e. amos, sand83-0643, may, 1983

 a subroutine package for bessel functions of a complex
 argument and nonnegative order by d. e. amos, sand85-
 1018, may, 1985

 a portable package for bessel functions of a complex
 argument and nonnegative order by d. e. amos, acm
 trans. math. software, vol. 12, no. 3, september 1986,
 pp 265-273.

 ***routines called  zacai,zbknu,zexp_sub,zsqrt_sub,zabs,i1mach,d1mach
 ***end prologue  zairy
 */
    /* Initialized data */
    static double tth = .666666666666666667;
    static double c1 = .35502805388781724;
    static double c2 = .258819403792806799;
    static double coef = .183776298473930683;
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;

    /* Local variables */
    static int nn, mr, iflag;
    static double bb, cc, fnu, cyi[1], cyr[1], ptr, sfac,
	csqi, ztai, csqr, ztar;
    /* complex ai,cone,csq,cy,s1,s2,trm1,trm2,z,zta,z3 */

    int verbose = *ierr;

    *nz = 0;
    if (*id < 0 || *id > 1 || *kode < 1 || *kode > 2) {
	*ierr = 1;	return;
    }
    *ierr = 0;
    double
	az = zabs(*zr, *zi),
	tol = fmax2(DBL_EPSILON, 1e-18),
	fid = (double) *id;

    if(verbose) Rprintf("zairy(z, id=%d, kode=%d, *), |z|=%g\n", *id, *kode, az);
    if (az > 1.) {
	goto L70;
    }

    if (az < tol) { // L170:
	double
	    aa = DBL_MIN * 1e3,
	    s1r = zeror,
	    s1i = zeroi;
	if (*id != 1) {
	    if (az > aa) {
		s1r = c2 * *zr;
		s1i = c2 * *zi;
	    }
	    // L180:
	    *air = c1 - s1r;
	    *aii = -s1i;
	} else { // *id == 1  L190:
	    *air = -c2;
	    *aii = 0.;
	    aa = sqrt(aa);
	    if (az > aa) {
		s1r = (*zr * *zr - *zi * *zi) * .5;
		s1i = *zr * *zi;
	    }
	    // L200:
	    *air += c1 * s1r;
	    *aii += c1 * s1i;
	}
	return;
    }

    /* -----------------------------------------------------------------------
       power series for  tol <= |z| <= 1.
       -----------------------------------------------------------------------*/
    double
	s1r = coner,
	s1i = conei,
	s2r = coner,
	s2i = conei,
	aa = az * az;
    if (aa >= tol / az) {
	double
	    trm1r = coner,
	    trm1i = conei,
	    trm2r = coner,
	    trm2i = conei,
	    atrm = 1.,
	    str = *zr * *zr - *zi * *zi,
	    sti = *zr * *zi + *zi * *zr,
	    z3r = str * *zr - sti * *zi,
	    z3i = str * *zi + sti * *zr,
	    az3 = az * aa,
	    ak = fid + 2.,
	    bk = 3. - fid - fid,
	    ck = 4. - fid,
	    dk = fid + 3. + fid,
	    d1 = ak * dk,
	    d2 = bk * ck,
	    ad = fmin2(d1,d2);
	ak = fid * 9. + 24.;
	bk = 30. - fid * 9.;
	for (int k = 1; k <= 25; ++k) {
	    str   = (trm1r * z3r - trm1i * z3i) / d1;
	    trm1i = (trm1r * z3i + trm1i * z3r) / d1;
	    trm1r = str;
	    s1r += trm1r;
	    s1i += trm1i;
	    str   = (trm2r * z3r - trm2i * z3i) / d2;
	    trm2i = (trm2r * z3i + trm2i * z3r) / d2;
	    trm2r = str;
	    s2r += trm2r;
	    s2i += trm2i;
	    atrm = atrm * az3 / ad;
	    d1 += ak;
	    d2 += bk;
	    ad = fmin2(d1,d2);
	    if (atrm < tol * ad)
		break;

	    ak += 18.;
	    bk += 18.;
	}
    }
    // L40:
    if (*id != 1) {
	*air = s1r * c1 - c2 * (*zr * s2r - *zi * s2i);
	*aii = s1i * c1 - c2 * (*zr * s2i + *zi * s2r);
	if (*kode != 1) {
	    double str, sti;
	    zsqrt_sub__(zr, zi, &str, &sti);
	    ztar = tth * (*zr * str - *zi * sti);
	    ztai = tth * (*zr * sti + *zi * str);
	    zexp_sub__(&ztar, &ztai, &str, &sti);
	    ptr = *air * str - *aii * sti;
	    *aii = *air * sti + *aii * str;
	    *air = ptr;
	}
	return;
    }
    // L50:  id == 1
    *air = -s2r * c2;
    *aii = -s2i * c2;
    if (az > tol) {
	double
	    str = *zr * s1r - *zi * s1i,
	    sti = *zr * s1i + *zi * s1r;
	cc = c1 / (fid + 1.);
	*air += cc * (str * *zr - sti * *zi);
	*aii += cc * (str * *zi + sti * *zr);
    }
    // L60:
    if (*kode != 1) {
	double str, sti;
	zsqrt_sub__(zr, zi, &str, &sti);
	ztar = tth * (*zr * str - *zi * sti);
	ztai = tth * (*zr * sti + *zi * str);
	zexp_sub__(&ztar, &ztai, &str, &sti);
	ptr = str * *air - sti * *aii;
	*aii = str * *aii + sti * *air;
	*air = ptr;
    }
    return;

L70:
    /* -----------------------------------------------------------------------
       case for cabs(z) > 1.0
       ----------------------------------------------------------------------- */
    fnu = (fid + 1.) / 3.;
    /* -----------------------------------------------------------------------
       set parameters related to machine constants.
       tol is the approximate unit roundoff limited to 1.0d-18.
       elim is the approximate exponential over- and underflow limit.
       exp(-elim) < exp(-alim)=exp(-elim)/tol    and
       exp(elim) > exp(alim)=exp(elim)*tol       are intervals near
       underflow and overflow limits where scaled arithmetic is done.
       rl is the lower boundary of the asymptotic expansion for large z.
       dig = number of base 10 digits in tol = 10**(-dig).
       ----------------------------------------------------------------------- */
    double
	r1m5 = M_LOG10_2,
	k = fmin2(-DBL_MIN_EXP, DBL_MAX_EXP),
	elim = (k * r1m5 - 3.) * 2.303;
    aa = r1m5 * (double) (DBL_MANT_DIG - 1);
    double dig = fmin2(aa, 18.);
    aa *= 2.303;
    double
	alim = elim + fmax2(-aa, -41.45),
	rl = dig * 1.2 + 3.,
	alaz = log(az);
/* --------------------------------------------------------------------------
   test for proper range
   ----------------------------------------------------------------------- */
    aa = .5 / tol;
    bb = (double) ((float) INT_MAX) * .5;
    aa = fmin2(aa, bb);
    aa = pow(aa, tth);
    if (az > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    zsqrt_sub__(zr, zi, &csqr, &csqi);
    ztar = tth * (*zr * csqr - *zi * csqi);
    ztai = tth * (*zr * csqi + *zi * csqr);
/* -----------------------------------------------------------------------
   re(zta) <= 0 when re(z) < 0, especially when im(z) is small
   ----------------------------------------------------------------------- */
    iflag = 0;
    sfac = 1.;
    double ak = ztai;
    if (*zr < 0.) {
	ztar = -fabs(ztar);
	ztai = ak;
    }
    // L80:
    if(*zi == 0. && *zr <= 0.) {
	ztar = 0.;
	ztai = ak;
    }
    // L90:
    aa = ztar;
    if (aa >= 0. && *zr > 0.) {
	goto L110;
    }
    if (*kode == 2) {
	goto L100;
    }
/* -----------------------------------------------------------------------
   overflow test
   ----------------------------------------------------------------------- */
    if (aa > -alim) {
	goto L100;
    }
    aa = -aa + alaz * .25;
    iflag = 1;
    sfac = tol;
    if (aa > elim) {
	goto L_err_2;
    }

L100:
/* -----------------------------------------------------------------------
   cbknu and cacon return exp(zta)*K(fnu,zta) on kode=2
   ----------------------------------------------------------------------- */
    mr = 1;
    if (*zi < 0.) {
	mr = -1;
    }
    zacai_(&ztar, &ztai, &fnu, kode, &mr, &c__1, cyr, cyi, &nn, &rl, &tol, &elim, &alim);
    if (nn < 0) {
	goto L280;
    }
    *nz += nn;
    goto L130;

L110:
    if (*kode == 2) {
	goto L120;
    }
/* -----------------------------------------------------------------------
   underflow test
   ----------------------------------------------------------------------- */
    if (aa < alim) {
	goto L120;
    }
    aa = -aa - alaz * .25;
    iflag = 2;
    sfac = 1. / tol;
    if (aa < -elim) {
	goto L210;
    }
L120:
    *nz = zbknu(&ztar, &ztai, &fnu, *kode, 1, verbose,
		cyr, cyi, tol, elim, alim);
L130:
    s1r = cyr[0] * coef;
    s1i = cyi[0] * coef;
    if (iflag == 0) {
	if (*id != 1) {
	    *air = csqr * s1r - csqi * s1i;
	    *aii = csqr * s1i + csqi * s1r;
	} else { // id == 1  L140:
	    *air = -(*zr * s1r - *zi * s1i);
	    *aii = -(*zr * s1i + *zi * s1r);
	}
	return;
    }
    // else iflag != 0  L150:
    s1r *= sfac;
    s1i *= sfac;
    if (*id != 1) {
	double
	    t = s1r * csqr - s1i * csqi;
	s1i   = s1r * csqi + s1i * csqr;
	s1r = t;
	*air = s1r / sfac;
	*aii = s1i / sfac;
    } else { // id == 1  L160:
	double
	    t = -(s1r * *zr - s1i * *zi);
	s1i   = -(s1r * *zi + s1i * *zr);
	s1r = t;
	*air = s1r / sfac;
	*aii = s1i / sfac;
    }
    return;


L210: *nz = 1;    *air = zeror;    *aii = zeroi;    return;

L_err_2: *nz = 0;    *ierr = 2;    return;
L280: if (nn == -1) goto L_err_2;
    // else
      *nz = 0;    *ierr = 5;    return;
L260: *nz = 0;    *ierr = 4;    return;
} /* zairy */

void zbiry(double *zr, double *zi, int *id, int *kode,
	   double *bir, double *bii, int *ierr)
{
/***begin prologue  zbiry
 ***date written   830501   (yymmdd)
 ***revision date  890801, 930101   (yymmdd)
 ***category no.  b5k
 ***keywords  airy function,bessel functions of order one third
 ***author  amos, donald e., sandia national laboratories
 ***purpose  to compute airy functions bi(z) and dbi(z) for complex z
 ***description

                      ***a double precision routine***
         on kode=1, cbiry computes the complex airy function bi(z) or
         its derivative dbi(z)/dz on id=0 or id=1 respectively. on
         kode=2, a scaling option cexp(-axzta)*bi(z) or cexp(-axzta)*
         dbi(z)/dz is provided to remove the exponential behavior in
         both the left and right half planes where
         zta=(2/3)*z*csqrt(z)=cmplx(xzta,yzta) and axzta=abs(xzta).
         definitions and notation are found in the nbs handbook of
         mathematical functions (ref. 1).

         input      zr,zi are double precision
           zr,zi  - z=cmplx(zr,zi)
           id     - order of derivative, id=0 or id=1
           kode   - a parameter to indicate the scaling option
                    kode= 1  returns
                             bi=bi(z)                 on id=0 or
                             bi=dbi(z)/dz             on id=1
                        = 2  returns
                             bi=cexp(-axzta)*bi(z)     on id=0 or
                             bi=cexp(-axzta)*dbi(z)/dz on id=1 where
                             zta=(2/3)*z*csqrt(z)=cmplx(xzta,yzta)
                             and axzta=abs(xzta)

         output     bir,bii are double precision
           bir,bii- complex answer depending on the choices for id and
                    kode
           ierr   - error flag
                    ierr=0, normal return - computation completed
                    ierr=1, input error   - no computation
                    ierr=2, overflow      - no computation, float(z)
                            too large on kode=1
                    ierr=3, cabs(z) large      - computation completed
                            losses of signifcance by argument reduction
                            produce less than half of machine accuracy
                    ierr=4, cabs(z) too large  - no computation
                            complete loss of accuracy by argument
                            reduction
                    ierr=5, error              - no computation,
                            algorithm termination condition not met

 ***long description

         bi and dbi are computed for cabs(z) > 1.0 from the I bessel
         functions by

                bi(z)=c*sqrt(z)*( I(-1/3,zta) + I(1/3,zta) )
               dbi(z)=c *  z  * ( I(-2/3,zta) + I(2/3,zta) )
                               c=1.0/sqrt(3.0)
                             zta=(2/3)*z**(3/2)

         with the power series for cabs(z) <= 1.0.

         in most complex variable computation, one must evaluate ele-
         mentary functions. when the magnitude of z is large, losses
         of significance by argument reduction occur. consequently, if
         the magnitude of zeta=(2/3)*z**1.5 exceeds u1=sqrt(0.5/ur),
         then losses exceeding half precision are likely and an error
         flag ierr=3 is triggered where ur=dmax1(d1mach(4),1.0d-18) is
         double precision unit roundoff limited to 18 digits precision.
         also, if the magnitude of zeta is larger than u2=0.5/ur, then
         all significance is lost and ierr=4. in order to use the int
         function, zeta must be further restricted not to exceed the
         largest int, u3=i1mach(9). thus, the magnitude of zeta
         must be restricted by min(u2,u3). on 32 bit machines, u1,u2,
         and u3 are approximately 2.0e+3, 4.2e+6, 2.1e+9 in single
         precision arithmetic and 1.3e+8, 1.8e+16, 2.1e+9 in double
         precision arithmetic respectively. this makes u2 and u3 limit-
         ing in their respective arithmetics. this means that the mag-
         nitude of z cannot exceed 3.1e+4 in single and 2.1e+6 in
         double precision arithmetic. this also means that one can
         expect to retain, in the worst cases on 32 bit machines,
         no digits in single precision and only 7 digits in double
         precision arithmetic. similar considerations hold for other
         machines.

         the approximate relative error in the magnitude of a complex
         bessel function can be expressed by p*10**s where p=max(unit
         roundoff,1.0e-18) is the nominal precision and 10**s repre-
         sents the increase in error due to argument reduction in the
         elementary functions. here, s=max(1,abs(log10(cabs(z))),
         abs(log10(fnu))) approximately (i.e. s=max(1,abs(exponent of
         cabs(z),abs(exponent of fnu)) ). however, the phase angle may
         have only absolute accuracy. this is most likely to occur when
         one component (in absolute value) is larger than the other by
         several orders of magnitude. if one component is 10**k larger
         than the other, then one can expect only max(abs(log10(p))-k,
         0) significant digits; or, stated another way, when k exceeds
         the exponent of p, no significant digits remain in the smaller
         component. however, the phase angle retains absolute accuracy
         because, in complex arithmetic with precision p, the smaller
         component will not (as a rule) decrease below p times the
         magnitude of the larger component. in these extreme cases,
         the principal phase angle is on the order of +p, -p, pi/2-p,
         or -pi/2+p.

 ***references  handbook of mathematical functions by m. abramowitz
                 and i. a. stegun, nbs ams series 55, u.s. dept. of
                 commerce, 1955.

               computation of bessel functions of complex argument
                 and large order by d. e. amos, sand83-0643, may, 1983

               a subroutine package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, sand85-
                 1018, may, 1985

               a portable package for bessel functions of a complex
                 argument and nonnegative order by d. e. amos, acm
                 trans. math. software, vol. 12, no. 3, september 1986,
                 pp 265-273.

 ***routines called  zbinu,zabs,zdiv,zsqrt_sub,d1mach,i1mach
 ***end prologue  zbiry
     complex bi,cone,csq,cy,s1,s2,trm1,trm2,z,zta,z3
 ***first executable statement  zbiry */

    if (*id < 0 || *id > 1 || *kode < 1 || *kode > 2) {
	*ierr = 1;
	return;
    }

    /* Initialized data */

    static double tth = .666666666666666667;
    static double c1 = .614926627446000736;
    static double c2 = .448288357353826359;
    static double coef = .577350269189625765;
    static double coner = 1.;
    static double conei = 0.;

    int verbose = *ierr;

    double bb, cc, rl,
	s1i, s2i, s1r, s2r, eaa, fnu, dig, cyi[2], cyr[2], sfac, alim,
	ztar, csqi, fnul, ztai, csqr;
    int nz = 0;
    *ierr = 0;
    double az = zabs(*zr, *zi),
	tol = fmax2(DBL_EPSILON, 1e-18),
	fid = (double) *id;

    if(verbose) Rprintf("zbiry(z, id=%d, kode=%d, *), |z|=%g\n", *id, *kode, az);

    if (az > 1.f) {
	goto L70;
    }
/* -----------------------------------------------------------------------
     power series for cabs(z) <= 1.
 ----------------------------------------------------------------------- */
    s1r = coner;
    s1i = conei;
    s2r = coner;
    s2i = conei;
    if (az < tol) {
	goto L130;
    }
    double aa = az * az;
    if (aa >= tol / az) {
	double
	    trm1r = coner,
	    trm1i = conei,
	    trm2r = coner,
	    trm2i = conei,
	    atrm = 1.,
	    str = *zr * *zr - *zi * *zi,
	    sti = *zr * *zi + *zi * *zr,
	    z3r = str * *zr - sti * *zi,
	    z3i = str * *zi + sti * *zr,
	    az3 = az * aa,
	    ak = fid + 2.,
	    bk = 3. - fid - fid,
	    ck = 4. - fid,
	    dk = fid + 3. + fid,
	    d1 = ak * dk,
	    d2 = bk * ck,
	    ad = fmin2(d1,d2);
	ak = fid * 9. + 24.;
	bk = 30. - fid * 9.;
	for (int k = 1; k <= 25; ++k) {
	    str = (trm1r * z3r - trm1i * z3i) / d1;
	    trm1i = (trm1r * z3i + trm1i * z3r) / d1;
	    trm1r = str;
	    s1r += trm1r;
	    s1i += trm1i;
	    str = (trm2r * z3r - trm2i * z3i) / d2;
	    trm2i = (trm2r * z3i + trm2i * z3r) / d2;
	    trm2r = str;
	    s2r += trm2r;
	    s2i += trm2i;
	    atrm = atrm * az3 / ad;
	    d1 += ak;
	    d2 += bk;
	    ad = fmin2(d1,d2);
	    if (atrm < tol * ad)
		break;

	    ak += 18.;
	    bk += 18.;
	}
    }
    // L40:
    if (*id != 1) {
	*bir = c1 * s1r + c2 * (*zr * s2r - *zi * s2i);
	*bii = c1 * s1i + c2 * (*zr * s2i + *zi * s2r);
	if (*kode != 1) {
	    double str, sti;
	    zsqrt_sub__(zr, zi, &str, &sti);
	    ztar = tth * (*zr * str - *zi * sti);
	    ztai = tth * (*zr * sti + *zi * str);
	    aa = ztar;
	    aa = -fabs(aa);
	    eaa = exp(aa);
	    *bir *= eaa;
	    *bii *= eaa;
	}
	return;
    }

    // L50:  *id == 1
    *bir = s2r * c2;
    *bii = s2i * c2;
    if (az > tol) {
	double
	    str = s1r * *zr - s1i * *zi,
	    sti = s1r * *zi + s1i * *zr;
	cc = c1 / (fid + 1.);
	*bir += cc * (str * *zr - sti * *zi);
	*bii += cc * (str * *zi + sti * *zr);
    }
    // L60:
    if (*kode != 1) {
	double str, sti;
	zsqrt_sub__(zr, zi, &str, &sti);
	ztar = tth * (*zr * str - *zi * sti);
	ztai = tth * (*zr * sti + *zi * str);
	aa = ztar;
	aa = -fabs(aa);
	eaa = exp(aa);
	*bir *= eaa;
	*bii *= eaa;
    }
    return;

/* -----------------------------------------------------------------------
     case for cabs(z) > 1.0
 ----------------------------------------------------------------------- */
L70:
    fnu = (fid + 1.) / 3.;
/* -----------------------------------------------------------------------
     set parameters related to machine constants.
     tol is the approximate unit roundoff limited to 1.0e-18.
     elim is the approximate exponential over- and underflow limit.
     exp(-elim) < exp(-alim)=exp(-elim)/tol    and
     exp(elim) > exp(alim)=exp(elim)*tol       are intervals near
     underflow and overflow limits where scaled arithmetic is done.
     rl is the lower boundary of the asymptotic expansion for large z.
     dig = number of base 10 digits in tol = 10**(-dig).
     fnul is the lower boundary of the asymptotic series for large fnu.
 ----------------------------------------------------------------------- */
    double
        k    = fmin2(-DBL_MIN_EXP, DBL_MAX_EXP),
        r1m5 = M_LOG10_2,
        elim = (k * r1m5 - 3.) * 2.303;

    aa = r1m5 * (DBL_MANT_DIG - 1);
    dig = fmin2(aa, 18.);
    aa *= 2.303;
    alim = elim + fmax2(-aa, -41.45);
    rl = dig * 1.2 + 3.;
    fnul = (dig - 3.) * 6. + 10.;
/* -----------------------------------------------------------------------
     test for range
 ----------------------------------------------------------------------- */
    aa = .5 / tol;
    bb = (double) ((float) INT_MAX) * .5;
    aa = fmin2(aa,bb);
    aa = pow(aa, tth);
    if (az > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    zsqrt_sub__(zr, zi, &csqr, &csqi);
    ztar = tth * (*zr * csqr - *zi * csqi);
    ztai = tth * (*zr * csqi + *zi * csqr);
/* -----------------------------------------------------------------------
     re(zta) <= 0 when re(z) < 0, especially when im(z) is small
 ----------------------------------------------------------------------- */
    sfac = 1.;
    double ak = ztai;
    if (*zr < 0.) {
	ztar = -fabs(ztar);
	ztai = ak;
    }
    // L80:
    if (*zi == 0. && *zr <= 0.) {
	ztar = 0.;
	ztai = ak;
    }
    // L90:
    aa = ztar;
    if (*kode != 2) {
	/* ------------------------------------------------------------------
	   overflow test
	   ------------------------------------------------------------------ */
	bb = fabs(aa);
	if (bb >= alim) {
	    bb += log(az) * .25;
	    sfac = tol;
	    if (bb > elim) {
		goto L190;
	    }
	}
    }
    // L100:
    double fmr = 0.;
    if (aa < 0. || *zr <= 0.) {
	if (*zi < 0.)
	    fmr = -pi;
	else
	    fmr = pi;
	ztar = -ztar;
	ztai = -ztai;
    }
    // L110:
/* -----------------------------------------------------------------------
     aa=factor for analytic continuation of I(fnu,zta)
     kode=2 returns exp(-abs(xzta))*I(fnu,zta) from zbesi
 ----------------------------------------------------------------------- */
    nz = zbinu(&ztar, &ztai, &fnu, *kode, 1, cyr, cyi, rl, fnul, tol, elim, alim);
    if (nz < 0) {
	goto L200;
    }
    aa = fmr * fnu;
    double
	str = cos(aa),
	sti = sin(aa);
    s1r = (str * cyr[0] - sti * cyi[0]) * sfac;
    s1i = (str * cyi[0] + sti * cyr[0]) * sfac;
    fnu = (2. - fid) / 3.;
    nz = zbinu(&ztar, &ztai, &fnu, *kode, 2, cyr, cyi, rl, fnul, tol, elim, alim);
    cyr[0] *= sfac;
    cyi[0] *= sfac;
    cyr[1] *= sfac;
    cyi[1] *= sfac;
/* -----------------------------------------------------------------------
     backward recur one step for orders -1/3 or -2/3
 ----------------------------------------------------------------------- */
    zdiv(*cyr, *cyi, ztar, ztai, &str, &sti);
    s2r = (fnu + fnu) * str + cyr[1];
    s2i = (fnu + fnu) * sti + cyi[1];
    aa = fmr * (fnu - 1.);
    str = cos(aa);
    sti = sin(aa);
    s1r = coef * (s1r + s2r * str - s2i * sti);
    s1i = coef * (s1i + s2r * sti + s2i * str);
    if (*id != 1) {
	str = csqr * s1r - csqi * s1i;
	s1i = csqr * s1i + csqi * s1r;
	s1r = str;
	*bir = s1r / sfac;
	*bii = s1i / sfac;
    } else { // L120:  *id == 1
	str = *zr * s1r - *zi * s1i;
	s1i = *zr * s1i + *zi * s1r;
	s1r = str;
	*bir = s1r / sfac;
	*bii = s1i / sfac;
    }
    return;

L130:
    aa = c1 * (1. - fid) + fid * c2;
    *bir = aa;
    *bii = 0.;
    return;

L190: nz = 0;	*ierr = 2;	return;

L200:
    if (nz == -1) goto L190;
    nz = 0;     *ierr = 5;	return;

L260: nz = 0;	*ierr = 4;	return;
} /* zbiry */

//-----------AUXILIARIES-------------------------------------------------------

void
zmlt(double ar, double ai, double br, double bi,
     double *cr, double *ci)
{
/* double precision complex multiply,  c = a*b  */
    *cr = ar * br - ai * bi;
    *ci = ar * bi + ai * br;
    return;
} /* zmlt */

void
zdiv(double ar, double ai, double br, double bi,
     double *cr, double *ci)
{
/* double precision complex divide  c = a/b */

    double
	bm = 1. / zabs(br, bi),
	cc = br * bm,
	cd = bi * bm;
    *cr = (ar * cc + ai * cd) * bm;
    *ci = (ai * cc - ar * cd) * bm;
    return;
} /* zdiv */

/* Subroutine */ int
zsqrt_sub__(double *ar, double *ai, double *br, double *bi)
{
/*  double precision complex square root,  b = csqrt(a)  */

    /* Initialized data */
    static double drt = .7071067811865475244008443621;

    /* Local variables */
    static double zm;
    static double dtheta;

    zm = zabs(*ar, *ai);
    zm = sqrt(zm);
    if (*ar == 0.) {
	goto L10;
    }
    if (*ai == 0.) {
	goto L20;
    }
    dtheta = atan(*ai / *ar);
    if (dtheta <= 0.) {
	goto L40;
    }
    if (*ar < 0.) {
	dtheta -= pi;
    }
    goto L50;
L10:
    if (*ai > 0.) {
	goto L60;
    }
    if (*ai < 0.) {
	goto L70;
    }
    *br = 0.;
    *bi = 0.;
    return 0;
L20:
    if (*ar > 0.) {
	goto L30;
    }
    *br = 0.;
    *bi = sqrt(fabs(*ar));
    return 0;
L30:
    *br = sqrt(*ar);
    *bi = 0.;
    return 0;
L40:
    if (*ar < 0.) {
	dtheta += pi;
    }
L50:
    dtheta *= .5;
    *br = zm * cos(dtheta);
    *bi = zm * sin(dtheta);
    return 0;
L60:
    *br = zm * drt;
    *bi = zm * drt;
    return 0;
L70:
    *br = zm * drt;
    *bi = -zm * drt;
    return 0;
} /* zsqrt_sub__ */

/* Subroutine */ int
zexp_sub__(double *ar, double *ai, double *br, double *bi)
{
/*  double precision complex exponential function  b = exp(a) */

    /* Builtin functions */
    double exp(double), cos(double), sin(double);

    /* Local variables */
    static double ca, cb, zm;

    zm = exp(*ar);
    ca = zm * cos(*ai);
    cb = zm * sin(*ai);
    *br = ca;
    *bi = cb;
    return 0;
} /* zexp_sub__ */

/* Subroutine */ int
zlog_sub__(double *ar, double *ai, double *br, double *bi, int *ierr)
{
/*   double precision complex logarithm  b = clog(a)
 *
 *   ierr=0,normal return      ierr=1, z=cmplx(0.0,0.0)
 */

    /* Initialized data */

    /* Builtin functions */
    double atan(double), log(double);

    /* Local variables */
    static double zm;
    static double dtheta;

    *ierr = 0;
    if (*ar == 0.) {
	goto L10;
    }
    if (*ai == 0.) {
	goto L20;
    }
    dtheta = atan(*ai / *ar);
    if (dtheta <= 0.) {
	goto L40;
    }
    if (*ar < 0.) {
	dtheta -= pi;
    }
    goto L50;
L10:
    if (*ai == 0.) {
	goto L60;
    }
    *bi = hpi;
    *br = log(fabs(*ai));
    if (*ai < 0.) {
	*bi = -(*bi);
    }
    return 0;
L20:
    if (*ar > 0.) {
	goto L30;
    }
    *br = log(fabs(*ar));
    *bi = pi;
    return 0;
L30:
    *br = log(*ar);
    *bi = 0.;
    return 0;
L40:
    if (*ar < 0.) {
	dtheta += pi;
    }
L50:
    zm = zabs(*ar, *ai);
    *br = log(zm);
    *bi = dtheta;
    return 0;
L60:
    *ierr = 1;
    return 0;
} /* zlog_sub__ */

double zabs(double zr, double zi)
{
/* zabs computes the absolute value or magnitude of a double
   precision complex variable cmplx(zr,zi)
*/

    /* Builtin functions */
    double sqrt(double);

    double
	u = fabs(zr),
	v = fabs(zi),
	s = u + v;
/* -----------------------------------------------------------------------
     s*1.0d0 makes an unnormalized underflow on cdc machines into a
     true floating zero
 ----------------------------------------------------------------------- */
    s *= 1.;
    if (s == 0.) { // L20:
	return 0.;
    } else if (u > v) { // L10:
	double q = v / u;
	return u * sqrt(q * q + 1.);
    } else {
	double q = u / v;
	return v * sqrt(q * q + 1.);
    }
} /* zabs */


int zbknu(double *zr, double *zi, double *fnu,
	  int kode, int n, int verbose,
	  double *yr, double *yi,
	  double tol, double elim, double alim)
{
/*
 *** refers to  zbesi,zbesk,zairy,zbesh

     zbknu computes the k bessel function in the right half z plane.

 *** routines called  dgamln,i1mach,d1mach,zkscl,zshch,zuchk,zabs,zdiv,
                    zexp_sub,zlog_sub,zmlt,zsqrt_sub
 */
    int nz = 0; // the return value

    /* Initialized data */
    static int kmax = 30;
    static double czeror = 0.;
    static double czeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;
    static double ctwor = 2.;
    static double r1 = 2.;
    static double spi = 1.90985931710274403; // = 6/pi
    static double fpi = 1.89769999331517738;
    static double tth = .666666666666666666;


    /* Local variables */
    int i, k, kk, inub, idum, kflag;
    double cyr[2], cyi[2],
	a1, a2, g1, g2, t1, t2, aa, bb, fc, ak, bk,
	fi, fk, as, fr, p_i,p_r, qi,qr, tm,
	p1i, p2i, s1i, s2i, p2m, p1r, p2r, s1r, s2r, cbi, cbr,
	cki, ckr, fhs, fks, rak, czi, czr, csi, csr,
	pti, ptr, zdi, zdr,
	cchi, cchr, alas, cshi,
	cshr, fmui, fmur,
	smui, smur, coefi, ascle, coefr;
    /*
     complex z,y,a,b,rz,smu,fu,fmu,f,flrz,cz,s1,s2,csh,cch
     complex ck,p,q,coef,p1,p2,cbk,pt,czero,cone,ctwo,st,ez,cs,dk
    */

    double
	caz = zabs(*zr, *zi),
	csclr = 1. / tol,
	cssr[3] = {csclr, 1., tol},
	csrr[3] = {tol, 1., csclr},
	by0 = DBL_MIN * 1e3 / tol,
	bry[3] = {by0, 1./by0, DBL_MAX};
    int
	iflag = 0,
	inu = (int) (*fnu + .5);// = "round"(fnu)
    double // fnu == inu + dnu  == "round"(fnu) + dnu;   dnu in  [ -1/2, 1/2 )
	dnu = *fnu - (double) inu,
	rcaz = 1. / caz,
	str = *zr * rcaz,
	sti = -(*zi) * rcaz,
	rzr = (str + str) * rcaz,
	rzi = (sti + sti) * rcaz;

    if(verbose) Rprintf("zbknu(z, nu, *): |z|=%g, inu=%d, dnu=%g \n", caz, inu, dnu);
    double dnu2 = 0.;
    if (fabs(dnu) == .5) {
	goto L110;
    }
    if (fabs(dnu) > tol) {
	dnu2 = dnu * dnu;
    }
    if (caz > r1) {
	goto L110;
    }
/* -----------------------------------------------------------------------
     series for cabs(z) <= r1
 ----------------------------------------------------------------------- */
    if(verbose) Rprintf(" .. series for |z| <= r1=%g:  ", r1);
    zlog_sub__(&rzr, &rzi, &smur, &smui, &idum);
    fmur = smur * dnu;
    fmui = smui * dnu;
    zshch_(&fmur, &fmui, &cshr, &cshi, &cchr, &cchi);
    if (dnu != 0.) {
	fc = dnu * pi;
	fc /= sin(fc);
	smur = cshr / dnu;
	smui = cshi / dnu;
    } else
	fc = 1.;
    if(verbose) Rprintf(" fc=%g, ", fc);
/* -----------------------------------------------------------------------
     gam(1-z)*gam(1+z)=pi*z/sin(pi*z), t1=1/gam(1-dnu), t2=1/gam(1+dnu)
 ----------------------------------------------------------------------- */
    t2 = exp(-dgamln(dnu + 1., &idum));
    t1 = 1. / (t2 * fc);
    if(verbose) Rprintf(" t2=exp(-log Gamma(.))=%g, t1 := 1/(t2*fc) = %g\n", t2, t1);
    if (fabs(dnu) <= .1) {
	/* --------------------------------------------------------------
	   series for f0 to resolve indeterminacy for small abs(dnu)
	   -------------------------------------------------------------- */
	static double cc[8] = { .577215664901532861,-.0420026350340952355,
				-.0421977345555443367,.00721894324666309954,
				-2.15241674114950973e-4,-2.01348547807882387e-5,
				1.13302723198169588e-6,6.11609510448141582e-9 };
	double ak = 1.,
	    s = cc[0];
	for (int k = 2; k <= 8; ++k) {
	    ak *= dnu2;
	    tm = cc[k - 1] * ak;
	    s += tm;
	    if (fabs(tm) < tol)
		break;
	}
	/* L30: */ g1 = -s;
    } else { // |dnu| > 0.1  (L40)
	g1 = (t1 - t2) / (dnu + dnu);
    }
    // L50:
    g2 = (t1 + t2) * .5;
    fr = fc * (cchr * g1 + smur * g2);
    fi = fc * (cchi * g1 + smui * g2);
    zexp_sub__(&fmur, &fmui, &str, &sti);
    p_r = str * .5 / t2;
    p_i = sti * .5 / t2;
    zdiv(0.5, 0., str, sti, &ptr, &pti);
    qr = ptr / t1;
    qi = pti / t1;
    s1r = fr;
    s1i = fi;
    s2r = p_r;
    s2i = p_i;
    ak = 1.;
    a1 = 1.;
    ckr = coner;
    cki = conei;
    bk = 1. - dnu2;
    if (inu <= 0 && n == 1) {
	/* ----------------------------------------------------------------
	   generate K(fnu,z), 0  <=  fnu  <  0.5  and  n=1
	   ---------------------------------------------------------------- */
	if(verbose) Rprintf(" .. zbknu(*, n=1, *): 0 <= fnu=%g < 0.5 (caz=%g >? tol=%g)\n",
			    *fnu, caz, tol);
	if (caz >= tol) {
	    zmlt(*zr, *zi, *zr, *zi, &czr, &czi);
	    czr *= .25;
	    czi *= .25;
	    t1 = caz * .25 * caz;
	    do { // L60:
		fr = (fr * ak + p_r + qr) / bk;
		fi = (fi * ak + p_i + qi) / bk;
		str = 1. / (ak - dnu);
		p_r *= str;
		p_i *= str;
		str = 1. / (ak + dnu);
		qr *= str;
		qi *= str;
		str = ckr * czr - cki * czi;
		rak = 1. / ak;
		cki = (ckr * czi + cki * czr) * rak;
		ckr = str * rak;
		s1r = ckr * fr - cki * fi + s1r;
		s1i = ckr * fi + cki * fr + s1i;
		a1 = a1 * t1 * rak;
		bk = bk + ak + ak + 1.;
		ak += 1.;
	    } while (a1 > tol);
	}
	// L70:
	if (kode == 1) {
	    yr[0] = s1r;
	    yi[0] = s1i;
	} else { // expon.scaled
	    zexp_sub__(zr, zi, &str, &sti);
	    zmlt(s1r, s1i, str, sti, yr, yi);
	}
	return nz;
    }
    // else: n > 1 or inu >= 1

    if(verbose)
	Rprintf(" .. 'L80': generate K(.) for fwd ref.(caz=%g, tol=%g, ak=%0f.)\n",
		caz, tol, ak);
/* L80:
   -----------------------------------------------------------------------
   generate K(dnu,z) and K(dnu+1,z) for forward recurrence
   ----------------------------------------------------------------------- */
    if (caz >= tol) {
	zmlt(*zr, *zi, *zr, *zi, &czr, &czi);
	czr *= .25;
	czi *= .25;
	t1 = caz * .25 * caz;
	do { // L90:
	    // MM: very parallel to L60 -- (FIXME? DRY!) -- additionally s2[ri]
	    fr = (fr * ak + p_r + qr) / bk;
	    fi = (fi * ak + p_i + qi) / bk;
	    str = 1. / (ak - dnu);
	    p_r *= str;
	    p_i *= str;
	    str = 1. / (ak + dnu);
	    qr *= str;
	    qi *= str;
	    str = ckr * czr - cki * czi;
	    rak = 1. / ak;
	    cki = (ckr * czi + cki * czr) * rak;
	    ckr = str * rak;
	    s1r = ckr * fr - cki * fi + s1r;
	    s1i = ckr * fi + cki * fr + s1i;
	    str = p_r - fr * ak;
	    sti = p_i - fi * ak;
	    s2r = ckr * str - cki * sti + s2r;
	    s2i = ckr * sti + cki * str + s2i;
	    a1 = a1 * t1 * rak;
	    bk = bk + ak + ak + 1.;
	    ak += 1.;
	} while (a1 > tol);
    }
    // L100:
    if(verbose) Rprintf("  .. 'L100': final ak=%.0f", ak);
    kflag = 2;
    a1 = *fnu + 1.;
    ak = a1 * fabs(smur);
    if (ak > alim) {
	kflag = 3;
    }
    if(verbose) Rprintf("  ak=%g >/<= alim=%g --> kflag=%d", ak, alim, kflag);

    str = cssr[kflag - 1];
    p2r = s2r * str;
    p2i = s2i * str;
    zmlt(p2r, p2i, rzr, rzi, &s2r, &s2i);
    s1r *= str;
    s1i *= str;
    if (kode != 1) {
	zexp_sub__(zr, zi, &fr, &fi);
	zmlt(s1r, s1i, fr, fi, &s1r, &s1i);
	zmlt(s2r, s2i, fr, fi, &s2r, &s2i);
    }
    goto Loop_fwd_rec;

/* -----------------------------------------------------------------------
     iflag=0 means no underflow occurred
     iflag=1 means an underflow occurred- computation proceeds with
     kode=2 and a test for on scale values is made during forward
     recursion
 ----------------------------------------------------------------------- */
L110:
    if(verbose) Rprintf("   .. 'L110': ");
    zsqrt_sub__(zr, zi, &str, &sti);
    zdiv(rthpi, czeroi, str, sti, &coefr, &coefi);
    kflag = 2;
    if (kode != 2) {
	if (*zr > alim) { // L290:
	    if(verbose) Rprintf(" zr > alim ==> kode = 2, iflag = 1;");
	    /* ---------------------------------------------------------
	       scale by dexp(z), iflag = 1 cases
	       --------------------------------------------------------- */
	    kode = 2;
	    iflag = 1; // kflag = 2;
	} else {
	    if(verbose) Rprintf(" normal: zr <= alim;");
	    str = exp(-(*zr)) * cssr[kflag - 1];
	    sti = -str * sin(*zi);
	    str *= cos(*zi);
	    zmlt(coefr, coefi, str, sti, &coefr, &coefi);
	}
    }
    // L120:
    if (fabs(dnu) == .5) {
	goto L300;
    }
/* -----------------------------------------------------------------------
     miller algorithm for cabs(z) > r1
 ----------------------------------------------------------------------- */
    ak = cos(pi * dnu);
    ak = fabs(ak);
    if (ak == czeror) {
	goto L300;
    }
    fhs = fabs(.25 - dnu2);
    if (fhs == czeror) {
	goto L300;
    }
/* -----------------------------------------------------------------------
     compute r2=f(e). if cabs(z) >= r2, use forward recurrence to
     determine the backward index k. r2=f(e) is a straight line on
     12 <= e <= 60. e is computed from 2**(-e)=b**(1-i1mach(14))=
     tol where b is the base of the arithmetic.
 ----------------------------------------------------------------------- */
    t1 = (double) (DBL_MANT_DIG - 1);
    t1 = t1 * M_LOG10_2 * 3.321928094;
    t1 = fmax2(t1, 12.);
    t1 = fmin2(t1, 60.);
    t2 = tth * t1 - 6.;
    if (*zr == 0.) {
	t1 = hpi;
    } else {
	t1 = fabs(atan(*zi / *zr));
    }
    if(verbose) Rprintf("   .. t1=%g, t2=%g\n", t1, t2);
    if (t2 <= caz) {
	/* ----------------------------------------------------------------
	   forward recurrence loop when cabs(z) >= r2
	   ---------------------------------------------------------------- */
	double etest = ak / (pi * caz * tol);
	if(verbose) Rprintf(" => |z| >= t2: etest=%g <?> coner=%g ", etest, coner);
	fk = coner;
	if (etest < coner) {
	    goto L180;
	}
	fks = ctwor;
	ckr = caz + caz + ctwor;
	p1r = czeror;
	p2r = coner;
	for (i = 1; i <= kmax; ++i) {
	    ak = fhs / fks;
	    cbr = ckr / (fk + coner);
	    ptr = p2r;
	    p2r = cbr * p2r - p1r * ak;
	    p1r = ptr;
	    ckr += ctwor;
	    fks = fks + fk + fk + ctwor;
	    fhs = fhs + fk + fk;
	    fk += coner;
	    str = fabs(p2r) * fk;
	    if (etest < str) {
		// L160 :
		fk += spi * t1 * sqrt(t2 / caz);
		fhs = fabs(.25 - dnu2);
		goto L180;
	    }
	}
        // L310:
	return -2;

    } else { // L170:
	if(verbose) Rprintf(" => |z| < t2 ");
	/* ----------------------------------------------------------------
	   compute backward index k for cabs(z) < t2
	   ---------------------------------------------------------------- */
	a2 = sqrt(caz);
	ak = fpi * ak / (tol * sqrt(a2));
	aa = t1 * 3. / (caz + 1.);
	bb = t1 * 14.7 / (caz + 28.);
	ak = (log(ak) + caz * cos(aa) / (caz * .008 + 1.)) / cos(bb);
	fk = ak * .12125 * ak / caz + 1.5;
    }
L180:
/* -----------------------------------------------------------------------
     backward recurrence loop for miller algorithm
 ----------------------------------------------------------------------- */
    k = (int) fk;
    if(verbose) Rprintf("  L180: bkw rec. for Miller algo, k=%d:\n", k);
    fk = (double) k;
    fks = fk * fk;
    p1r = czeror;
    p1i = czeroi;
    p2r = tol;
    p2i = czeroi;
    csr = p2r;
    csi = p2i;
    for (i = 1; i <= k; ++i) {
	a1 = fks - fk;
	ak = (fks + fk) / (a1 + fhs);
	rak = 2. / (fk + coner);
	cbr = (fk + *zr) * rak;
	cbi = *zi * rak;
	ptr = p2r;
	pti = p2i;
	p2r = (ptr * cbr - pti * cbi - p1r) * ak;
	p2i = (pti * cbr + ptr * cbi - p1i) * ak;
	p1r = ptr;
	p1i = pti;
	csr += p2r;
	csi += p2i;
	fks = a1 - fk + coner;
	fk -= coner;
    }
/* -----------------------------------------------------------------------
     compute (p2/cs)=(p2/cabs(cs))*(conjg(cs)/cabs(cs)) for better
     scaling
 ----------------------------------------------------------------------- */
    tm = zabs(csr, csi);
    ptr = 1. / tm;
    s1r = p2r * ptr;
    s1i = p2i * ptr;
    csr *= ptr;
    csi = -csi * ptr;
    zmlt(coefr, coefi, s1r, s1i, &str, &sti);
    zmlt(str, sti, csr, csi, &s1r, &s1i);
    if (inu <= 0 && n == 1) {
	zdr = *zr;
	zdi = *zi;
	if (iflag == 1) {
	    goto L270;
	}
	goto L240;
    }

    // else L200: inu > 0 ||  n > 1
/* -----------------------------------------------------------------------
     compute p1/p2=(p1/cabs(p2)*conjg(p2)/cabs(p2) for scaling
 ----------------------------------------------------------------------- */
    tm = zabs(p2r, p2i);
    ptr = 1. / tm;
    p1r *= ptr;
    p1i *= ptr;
    p2r *= ptr;
    p2i = -p2i * ptr;
    zmlt(p1r, p1i, p2r, p2i, &ptr, &pti);
    str = dnu + .5 - ptr;
    sti = -pti;
    zdiv(str, sti, *zr, *zi, &str, &sti);
    str += 1.;
    zmlt(str, sti, s1r, s1i, &s2r, &s2i);

Loop_fwd_rec: // L210
/* -----------------------------------------------------------------------
     forward recursion on the three term recursion with relation with
     scaling near exponent extremes on kflag=1 or kflag=3
 ----------------------------------------------------------------------- */
    str = dnu + 1.;
    ckr = str * rzr;
    cki = str * rzi;
    if (n == 1) {
	--inu;
    }
    if(verbose) Rprintf(" 'Loop_fwd_rec': ck = %g + i*%g, inu=%d (<= 0 ?), iflag=%d\n",
			ckr, cki, inu, iflag);
    if (inu <= 0) {
	if (n <= 1) {
	    s1r = s2r;
	    s1i = s2i;
	}
	// L215:
	zdr = *zr;
	zdi = *zi;
	if (iflag == 1)
	    goto L270;
	// else goto L240;

    } else { // inu > 0
        // L220:
	inub = 1;
	if (iflag == 1) {
	    goto L261;
	}
    L225:
	p1r  = csrr[kflag - 1];
	ascle = bry[kflag - 1];
	for (i = inub; i <= inu; ++i) {
	    str = s2r;
	    sti = s2i;
	    s2r = ckr * str - cki * sti + s1r;
	    s2i = ckr * sti + cki * str + s1i;
	    s1r = str;
	    s1i = sti;
	    ckr += rzr;
	    cki += rzi;
	    if (kflag >= 3)
		continue; // for(i ..) loop

	    p2r = s2r * p1r;
	    p2i = s2i * p1r;
	    str = fabs(p2r);
	    sti = fabs(p2i);
	    p2m = fmax2(str,sti);
	    if (p2m <= ascle)
		continue; // for(i ..) loop

	    ++kflag;
	    ascle = bry[kflag - 1];
	    s1r *= p1r;
	    s1i *= p1r;
	    s2r = p2r;
	    s2i = p2i;
	    str = cssr[kflag - 1];
	    s1r *= str;
	    s1i *= str;
	    s2r *= str;
	    s2i *= str;
	    p1r = csrr[kflag - 1];
	    // L230: ;
	}
	if (n == 1) {
	    s1r = s2r;
	    s1i = s2i;
	}
    } // else
L240:
    str = csrr[kflag - 1];
    yr[0] = s1r * str;
    yi[0] = s1i * str;
    if (n == 1) {
	return nz;
    }
    yr[1] = s2r * str;
    yi[1] = s2i * str;
    if (n == 2) {
	return nz;
    }
    kk = 2;
L250: //---  for(kk=3, kk <= n; ++kk)  -----------
    ++kk;
    if (kk > n) {
	return nz;
    }
    p1r  = csrr[kflag - 1];
    ascle = bry[kflag - 1];
    for (i = kk-1; i < n; ++i) {
	p2r = s2r;
	p2i = s2i;
	s2r = ckr * p2r - cki * p2i + s1r;
	s2i = cki * p2r + ckr * p2i + s1i;
	s1r = p2r;
	s1i = p2i;
	ckr += rzr;
	cki += rzi;
	p2r = s2r * p1r;
	p2i = s2i * p1r;
	yr[i] = p2r;
	yi[i] = p2i;
	if (kflag >= 3)
	    continue; // for(i ..) loop

	str = fabs(p2r);
	sti = fabs(p2i);
	p2m = fmax2(str,sti);
	if (p2m <= ascle)
	    continue; // for(i ..) loop

	++kflag;
	ascle = bry[kflag - 1];
	s1r *= p1r;
	s1i *= p1r;
	s2r = p2r;
	s2i = p2i;
	str = cssr[kflag - 1];
	s1r *= str;
	s1i *= str;
	s2r *= str;
	s2i *= str;
	p1r = csrr[kflag - 1];
    }
    return nz;

/* -----------------------------------------------------------------------
     iflag=1 cases, forward recurrence on scaled values on underflow
 ----------------------------------------------------------------------- */
L261:
  {
    double helim = elim * .5,
	elm = exp(-elim),
	celmr = elm;
    int ic = -1,
	j = 2;
    if(verbose) Rprintf("  L261: elm=%g, inu=%d\n", elm, inu);
    ascle = bry[0];
    zdr = *zr;
    zdi = *zi;
    for (i = 1; i <= inu; ++i) {
	str = s2r;
	sti = s2i;
	s2r = str * ckr - sti * cki + s1r;
	s2i = sti * ckr + str * cki + s1i;
	s1r = str;
	s1i = sti;
	ckr += rzr;
	cki += rzi;
	as = zabs(s2r, s2i);
	alas = log(as);
	p2r = -zdr + alas;
	if (p2r < -elim) {
	    goto L263;
	}
	zlog_sub__(&s2r, &s2i, &str, &sti, &idum);
	p2r = -zdr + str;
	p2i = -zdi + sti;
	p2m = exp(p2r) / tol;
	p1r = p2m * cos(p2i);
	p1i = p2m * sin(p2i);
	if (zuchk(p1r, p1i, ascle, tol) == 0) {
	    j = 3 - j; // 1 <--> 2 (flip)
	    cyr[j - 1] = p1r;
	    cyi[j - 1] = p1i;
	    if (ic == i - 1) {
		// L264:
		kflag = 1;
		inub = i + 1;
		s2r = cyr[j - 1];
		s2i = cyi[j - 1];
		j = 3 - j;
		s1r = cyr[j - 1];
		s1i = cyi[j - 1];
		if (inub <= inu) {
		    goto L225;
		}
		if (n == 1) {
		    s1r = s2r;
		    s1i = s2i;
		}
		goto L240;
	    }
	    ic = i;

	} else {
 L263:
	    if (alas >= helim) {
		zdr -= elim;
		s1r *= celmr;
		s1i *= celmr;
		s2r *= celmr;
		s2i *= celmr;
	    }
	}
	// L262:
    } // for(i ..)

    if (n == 1) {
	s1r = s2r;
	s1i = s2i;
    }

L270:
    yr[0] = s1r;  yi[0] = s1i;
    if (n > 1) {
	yr[1] = s2r;  yi[1] = s2i;
    }
    // L280:
    ascle = bry[0];
    zkscl_(&zdr, &zdi, fnu, &n, yr, yi, &nz, &rzr, &rzi, &ascle, &tol, &elim);
    inu = n - nz;
    if(verbose) Rprintf("  L270, after zkscl(): nz=%d, new inu=%d\n", nz, inu);
    if (inu <= 0) {
	return nz;
    }
    kk = nz;
    s1r = yr[kk];
    s1i = yi[kk];
    yr[kk] = s1r * csrr[0];
    yi[kk] = s1i * csrr[0];
    if (inu == 1) {
	return nz;
    }
    kk = nz + 1;
    s2r = yr[kk];
    s2i = yi[kk];
    yr[kk] = s2r * csrr[0];
    yi[kk] = s2i * csrr[0];
    if (inu == 2) {
	return nz;
    }
    t2 = *fnu + (double) kk;
    ckr = t2 * rzr;
    cki = t2 * rzi;
    kflag = 1;
    goto L250;
  }

/* -----------------------------------------------------------------------
     fnu=half odd int case, dnu=-0.5
 ----------------------------------------------------------------------- */
L300:
    if(verbose) Rprintf(" 'L300': coef = %g + i*%g\n", coefr, coefi);
    s1r = coefr;
    s1i = coefi;
    s2r = coefr;
    s2i = coefi;
    goto Loop_fwd_rec;

} /* zbknu */

/* Subroutine */ int
zkscl_(double *zrr, double *zri, double *fnu,
       int *n, double *yr, double *yi, int *nz, double * rzr, double *rzi,
       double *ascle, double *tol, double *elim)
{
/*
 ***refer to  zbesk

     set k functions to zero on underflow, continue recurrence
     on scaled functions until two members come on scale, then
     return with min(nz+2,n) values scaled by 1/tol.

 ***routines called  zuchk,zabs,zlog_sub
 */
    /* Initialized data */

    static double zeror = 0.;
    static double zeroi = 0.;

    /* Local variables */
    int i, ic, kk, nn, idum;
    double as, fn, s1i, s2i, s1r, s2r, acs, cki, elm, ckr,
	cyi[2], cyr[2], csi, csr, zdi, zdr, str,
	alas, helim, celmr;
    /* complex ck,cs,cy,czero,rz,s1,s2,y,zr,zd,celm

     Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    ic = 0;
    nn = imin2(2,*n);
    for (i = 1; i <= nn; ++i) {
	s1r = yr[i];
	s1i = yi[i];
	cyr[i - 1] = s1r;
	cyi[i - 1] = s1i;
	as = zabs(s1r, s1i);
	acs = -(*zrr) + log(as);
	++(*nz);
	yr[i] = zeror;
	yi[i] = zeroi;
	if (acs >= -(*elim)) {
	    zlog_sub__(&s1r, &s1i, &csr, &csi, &idum);
	    csr -= *zrr;
	    csi -= *zri;
	    str = exp(csr) / *tol;
	    csr = str * cos(csi);
	    csi = str * sin(csi);
	    if (zuchk(csr, csi, *ascle, *tol) == 0) {
		yr[i] = csr;
		yi[i] = csi;
		ic = i;
		--(*nz);
	    }
	}
    }
    if (*n == 1) {
	return 0;
    }
    if (ic <= 1) {
	yr[1] = zeror;
	yi[1] = zeroi;
	*nz = 2;
    }
    // L20:
    if (*n == 2 || *nz == 0) {
	return 0;
    }
    fn = *fnu + 1.;
    ckr = fn * *rzr;
    cki = fn * *rzi;
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    helim = *elim * .5;
    elm = exp(-(*elim));
    celmr = elm;
    zdr = *zrr;
    zdi = *zri;

/*     find two consecutive y values on scale. scale recurrence if
     s2 gets larger than exp(elim/2) */

    for (i = 3; i <= *n; ++i) {
	kk = i;
	csr = s2r;
	csi = s2i;
	s2r = ckr * csr - cki * csi + s1r;
	s2i = cki * csr + ckr * csi + s1i;
	s1r = csr;
	s1i = csi;
	ckr += *rzr;
	cki += *rzi;
	as = zabs(s2r, s2i);
	alas = log(as);
	acs = -zdr + alas;
	++(*nz);
	yr[i] = zeror;
	yi[i] = zeroi;
	if (acs >= -(*elim)) {
	    zlog_sub__(&s2r, &s2i, &csr, &csi, &idum);
	    csr -= zdr;
	    csi -= zdi;
	    str = exp(csr) / *tol;
	    csr = str * cos(csi);
	    csi = str * sin(csi);
	    if (zuchk(csr, csi, *ascle, *tol) == 0) {
		yr[i] = csr;
		yi[i] = csi;
		--(*nz);
		if (ic == kk - 1) {
		    // L40:
		    *nz = kk - 2;
		    goto L45;
		}
		ic = kk;
		continue; // goto L30
	    }
	}
	// L25:
	if (alas >= helim) {
	    zdr -= *elim;
	    s1r *= celmr;
	    s1i *= celmr;
	    s2r *= celmr;
	    s2i *= celmr;
	}
    } // for(i ..)

    *nz = (ic == *n) ? *n : (*n - 1);
L45:
    for (i = 1; i <= *nz; ++i) {
	yr[i] = zeror;
	yi[i] = zeroi;
    }
    return 0;
} /* zkscl_ */

/* Subroutine */ int
zshch_(double *zr, double *zi, double *cshr, double *cshi,
       double *cchr, double *cchi)
{
/***begin prologue  zshch
 ***refer to  zbesk,zbesh

     zshch computes the complex hyperbolic functions csh=sinh(x+i*y)
     and cch=cosh(x+i*y), where i**2=-1.

 ***end prologue  zshch
*/

    /* Builtin functions */
    double sinh(double), cosh(double), sin(double), cos(double);

    /* Local variables */
    double
	sh = sinh(*zr),
	ch = cosh(*zr),
	sn = sin(*zi),
	cn = cos(*zi);
    *cshr = sh * cn;
    *cshi = ch * sn;
    *cchr = ch * cn;
    *cchi = sh * sn;
    return 0;
} /* zshch_ */

void
zrati(double zr, double zi, double fnu, int n, double tol,
      double *cyr, double *cyi)
{
/*
 ***refer to  zbesi,zbesk,zbesh

 zrati computes ratios of i bessel functions by backward
 recurrence.  the starting index is determined by forward
 recurrence as described in j. res. of nat. bur. of standards-b,
 mathematical sciences, vol 77b, p111-114, september, 1973,
 bessel functions i and j of complex argument and int order,
 by d. j. sookne.

 ***routines called  zabs,zdiv

 */
    /* Initialized data */

    static double czeror = 0.;
    static double czeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;
    static double rt2 = 1.41421356237309505;
    /* complex z,cy(1),cone,czero,p1,p2,t1,rz,pt,cdfnu */

    --cyi;
    --cyr;

    /* Function Body */
    int
	inu  = (int) fnu,
	idnu = inu + n - 1;
    double
	az = zabs(zr, zi);
    int magz = (int) az;
    double fnup = fmax((double)(magz + 1), (double) idnu);
    int id = idnu - magz - 1,
	itime = 1,
	k = 1;
    double
	ptr = 1. / az,
	rzr =  ptr * (zr + zr) * ptr,
	rzi = -ptr * (zi + zi) * ptr,
	t1r = rzr * fnup,
	t1i = rzi * fnup,
	p2r = -t1r,
	p2i = -t1i,
	p1r = coner,
	p1i = conei;
    t1r += rzr;
    t1i += rzi;
    if (id > 0)
	id = 0;
    double
	ap2 = zabs(p2r, p2i),
	ap1 = zabs(p1r, p1i),
/* -----------------------------------------------------------------------
   the overflow test on K(fnu+i-1,z) before the call to cbknu
   guarantees that p2 is on scale.  Scale test1 and all subsequent
   p2 values by ap1 to ensure that an overflow does not occur prematurely.
   ----------------------------------------------------------------------- */
	arg = (ap2 + ap2) / (ap1 * tol),
	test1 = sqrt(arg),
	test = test1,
	rap1 = 1. / ap1;
    p1r *= rap1;
    p1i *= rap1;
    p2r *= rap1;
    p2i *= rap1;
    ap2 *= rap1;

Loop:
    ++k;
    ap1 = ap2;
    ptr = p2r;
    double pti = p2i;
    p2r = p1r - (t1r * ptr - t1i * pti);
    p2i = p1i - (t1r * pti + t1i * ptr);
    p1r = ptr;
    p1i = pti;
    t1r += rzr;
    t1i += rzi;
    ap2 = zabs(p2r, p2i);
    if (ap1 <= test) {
	goto Loop;
    }
    if (itime == 2) {
	goto L20;
    }
    double
	ak = zabs(t1r, t1i) * .5,
	flam = ak + sqrt(ak * ak - 1.),
	rho = fmin2(ap2 / ap1, flam);
    test = test1 * sqrt(rho / (rho * rho - 1.));
    itime = 2;
    goto Loop;
    //-------------
    int kk;
L20:
    kk = k + 1 - id;
    ak = (double) kk;
    t1r = ak;
    t1i = czeroi;
    double dfnu = fnu + (double) (n - 1);
    p1r = 1. / ap2;
    p1i = czeroi;
    p2r = czeror;
    p2i = czeroi;
    for (int i = 1; i <= kk; ++i) {
	ptr = p1r;
	pti = p1i;
	rap1 = dfnu + t1r;
	double
	    ttr = rzr * rap1,
	    tti = rzi * rap1;
	p1r = ptr * ttr - pti * tti + p2r;
	p1i = ptr * tti + pti * ttr + p2i;
	p2r = ptr;
	p2i = pti;
	t1r -= coner;
    }
    if (p1r == czeror && p1i == czeroi) {
	p1r = tol;
	p1i = tol;
    }
    zdiv(p2r, p2i, p1r, p1i, &cyr[n], &cyi[n]);
    if (n == 1) {
	return;
    }
    k = n - 1;
    ak = (double) k;
    t1r = ak;
    t1i = czeroi;
    double
	cdfnur = fnu * rzr,
	cdfnui = fnu * rzi;
    for (int i = 2; i <= n; ++i) {
	ptr = cdfnur + (t1r * rzr - t1i * rzi) + cyr[k + 1];
	pti = cdfnui + (t1r * rzi + t1i * rzr) + cyi[k + 1];
	ak = zabs(ptr, pti);
	if (ak == czeror) {
	    ptr = tol;
	    pti = tol;
	    ak = tol * rt2;
	}
	double rak = coner / ak;
	cyr[k] =  rak * ptr * rak;
	cyi[k] = -rak * pti * rak;
	t1r -= coner;
	--k;
    }
    return;
} /* zrati */

void
zs1s2_(double *zrr, double *zri, double *s1r,
       double *s1i, double *s2r, double *s2i, int *nz,
       double *ascle, double *alim, int *iuf)
{
/*
 ***refer to  zbesk,zairy

     zs1s2 tests for a possible underflow resulting from the
     addition of the i and k functions in the analytic con-
     tinuation formula where s1=k function and s2=i function.
     on kode=1 the i and k functions are different orders of
     magnitude, but for kode=2 they can be of the same order
     of magnitude and the maximum must be at least one
     precision above the underflow limit.

 ***routines called  zabs,zexp_sub,zlog_sub
 */

    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;

    *nz = 0;
    double // complex s1
	as1 = zabs(*s1r, *s1i),
	as2 = zabs(*s2r, *s2i);
    if ((*s1r == 0. && *s1i == 0.) || as1 == 0.) {
	goto L10;
    }
    double
	aln = -(*zrr) - *zrr + log(as1),
	s1dr = *s1r,
	s1di = *s1i;
    *s1r = zeror;
    *s1i = zeroi;
    as1 = zeror;
    if (aln >= -(*alim)) {
	int idum;
	double c1i, c1r;
	zlog_sub__(&s1dr, &s1di, &c1r, &c1i, &idum);
	c1r = c1r - *zrr - *zrr;
	c1i = c1i - *zri - *zri;
	zexp_sub__(&c1r, &c1i, s1r, s1i);
	as1 = zabs(*s1r, *s1i);
	++(*iuf);
    }
L10:
    if (fmax2(as1,as2) <= *ascle) {
	*s1r = zeror;
	*s1i = zeroi;
	*s2r = zeror;
	*s2i = zeroi;
	*nz = 1;
	*iuf = 0;
    }
    return;
} /* zs1s2_ */

/* Subroutine */ int
zbunk_(double *zr, double *zi, double *fnu,
       int *kode, int *mr, int *n, double *yr, double *yi,
       int *nz, double *tol, double *elim, double *alim)
{
/***begin prologue  zbunk
 ***refer to  zbesk,zbesh

     zbunk computes the k bessel function for fnu > fnul.
     according to the uniform asymptotic expansion for K(fnu,z)
     in zunk1 and the expansion for H(2,fnu,z) in zunk2

 ***routines called  zunk1,zunk2
 ***end prologue  zbunk
 */
    static double ax, ay;

    /* complex y,z */

    /* Function Body */
    *nz = 0;
    ax = fabs(*zr) * 1.7321;
    ay = fabs(*zi);
    if (ay <= ax) {
	/* -----------------------------------------------------------------------
	   asymptotic expansion for K(fnu,z) for large fnu applied in
	   -pi/3 <= arg(z) <= pi/3
	   ----------------------------------------------------------------------- */
	zunk1_(zr, zi, fnu, kode, mr, n, yr,yi, nz, tol, elim, alim);
    } else {
	/* -----------------------------------------------------------------------
	   asymptotic expansion for H(2,fnu,z*exp(m*hpi)) for large fnu
	   applied in pi/3 < abs(arg(z)) <= pi/2 where m=+i or -i
	   and hpi=pi/2
	   ----------------------------------------------------------------------- */
	zunk2_(zr, zi, fnu, kode, mr, n, yr,yi, nz, tol, elim, alim);
    }
    return 0;
} /* zbunk_ */

int zmlri(double *zr, double *zi, double *fnu, int kode, int n,
	  double *yr, double *yi, double tol)
{
/*
 ***refer to  zbesi,zbesk

     zmlri computes the i bessel function for re(z) >= 0.0 by the
     miller algorithm normalized by a neumann series.

 ***routines called  dgamln,d1mach,zabs,zexp_sub,zlog_sub,zmlt
 */

    /* Initialized data */

    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;

    /* Local variables */
    int i, k, m, iaz, inu, kk, km, idum, ifnu, itime;
    double ak, bk, ap, at, az, p1i, p2i, p1r,p2r, ack, cki,ckr, fnf, fkk,
	rho, pti,ptr, raz, sti,str, rzi,rzr,
	rho2, flam, fkap, scle, tst, tfnf, sumi,sumr, cnormi,cnormr;

    /* complex ck,cnorm,cone,ctwo,czero,pt,p1,p2,rz,sum,y,z
     Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
    scle = DBL_MIN / tol;
    az = zabs(*zr, *zi);
    iaz = (int) az;
    ifnu = (int) *fnu;
    inu = ifnu + n - 1;
    at = ((double) iaz) + 1.;
    raz = 1. / az;
    str =   *zr  * raz;
    sti = -(*zi) * raz;
    ckr = str * at * raz;
    cki = sti * at * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    p1r = zeror;
    p1i = zeroi;
    p2r = coner;
    p2i = conei;
    ack = (at + 1.) * raz;
    rho = ack + sqrt(ack * ack - 1.);
    rho2 = rho * rho;
    tst = (rho2 + rho2) / ((rho2 - 1.) * (rho - 1.));
    tst /= tol;
/* -----------------------------------------------------------------------
     compute relative truncation error index for series
 ----------------------------------------------------------------------- */
    ak = at;
    for (i = 1; i <= 80; ++i) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r - (ckr * ptr - cki * pti);
	p2i = p1i - (cki * ptr + ckr * pti);
	p1r = ptr;
	p1i = pti;
	ckr += rzr;
	cki += rzi;
	ap = zabs(p2r, p2i);
	if (ap > tst * ak * ak) {
	    goto L20;
	}
	ak += 1.;
    }
    return -2;
L20:
    ++i;
    k = 0;
    if (inu < iaz) {
	goto L40;
    }
/* -----------------------------------------------------------------------
     compute relative truncation error for ratios
 ----------------------------------------------------------------------- */
    p1r = zeror;
    p1i = zeroi;
    p2r = coner;
    p2i = conei;
    at = (double) (inu + 1);
    str = *zr * raz;
    sti = -(*zi) * raz;
    ckr = str * at * raz;
    cki = sti * at * raz;
    ack = at * raz;
    tst = sqrt(ack / tol);
    itime = 1;
    for (k = 1; k <= 80; ++k) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r - (ckr * ptr - cki * pti);
	p2i = p1i - (ckr * pti + cki * ptr);
	p1r = ptr;
	p1i = pti;
	ckr += rzr;
	cki += rzi;
	ap = zabs(p2r, p2i);
	if (ap < tst)
	    continue; // for(k ..)
	if (itime == 2) {
	    goto L40;
	}
	ack = zabs(ckr, cki);
	flam = ack + sqrt(ack * ack - 1.);
	fkap = ap / zabs(p1r, p1i);
	rho = fmin2(flam,fkap);
	tst *= sqrt(rho / (rho * rho - 1.));
	itime = 2; // leave next time
    }
    return -2; // not converged in 80 steps
L40:
/* -----------------------------------------------------------------------
     backward recurrence and sum normalizing relation
 ----------------------------------------------------------------------- */
    ++k;
    kk = fmax2(i + iaz, k + inu);
    fkk = (double) kk;
    p1r = zeror;
    p1i = zeroi;
/* -----------------------------------------------------------------------
     scale p2 and sum by scle
 ----------------------------------------------------------------------- */
    p2r = scle;
    p2i = zeroi;
    fnf = *fnu - (double) ifnu;
    tfnf = fnf + fnf;
    bk = dgamln(fkk + tfnf + 1., &idum)
	- dgamln(fkk + 1., &idum)
	- dgamln(tfnf + 1., &idum);
    bk = exp(bk);
    sumr = zeror;
    sumi = zeroi;
    km = kk - inu;
    for (i = 1; i <= km; ++i) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
	p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
	p1r = ptr;
	p1i = pti;
	ak = 1. - tfnf / (fkk + tfnf);
	ack = bk * ak;
	sumr += (ack + bk) * p1r;
	sumi += (ack + bk) * p1i;
	bk = ack;
	fkk += -1.;
    }
    yr[n] = p2r;
    yi[n] = p2i;
    if (n == 1) {
	goto L70;
    }
    for (i = 2; i <= n; ++i) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
	p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
	p1r = ptr;
	p1i = pti;
	ak = 1. - tfnf / (fkk + tfnf);
	ack = bk * ak;
	sumr += (ack + bk) * p1r;
	sumi += (ack + bk) * p1i;
	bk = ack;
	fkk += -1.;
	m = n - i + 1;
	yr[m] = p2r;
	yi[m] = p2i;
    }
L70:
    if (ifnu <= 0) {
	goto L90;
    }
    for (i = 1; i <= ifnu; ++i) {
	ptr = p2r;
	pti = p2i;
	p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
	p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
	p1r = ptr;
	p1i = pti;
	ak = 1. - tfnf / (fkk + tfnf);
	ack = bk * ak;
	sumr += (ack + bk) * p1r;
	sumi += (ack + bk) * p1i;
	bk = ack;
	fkk += -1.;
    }
L90:
    ptr = *zr;
    pti = *zi;
    if (kode == 2) {
	ptr = zeror;
    }
    zlog_sub__(&rzr, &rzi, &str, &sti, &idum);
    p1r = -fnf * str + ptr;
    p1i = -fnf * sti + pti;
    ap = dgamln(fnf + 1., &idum);
    ptr = p1r - ap;
    pti = p1i;
/* -----------------------------------------------------------------------
     the division cexp(pt)/(sum+p2) is altered to avoid overflow
     in the denominator by squaring large quantities
 ----------------------------------------------------------------------- */
    p2r += sumr;
    p2i += sumi;
    ap = zabs(p2r, p2i);
    p1r = 1. / ap;
    zexp_sub__(&ptr, &pti, &str, &sti);
    ckr = str * p1r;
    cki = sti * p1r;
    ptr = p2r * p1r;
    pti = -p2i * p1r;
    zmlt(ckr, cki, ptr, pti, &cnormr, &cnormi);
    for (i = 1; i <= n; ++i) {
	str = yr[i] * cnormr - yi[i] * cnormi;
	yi[i] = yr[i] * cnormi + yi[i] * cnormr;
	yr[i] = str;
    }
    return 0;
} // zmlri

// TODO: int function() returning current '*nz'
void
zwrsk_(double *zrr, double *zri, double *fnu,
       int *kode, int *n, double *yr, double *yi, int *nz,
       double *cwr, double *cwi, double *tol, double *elim, double *alim)
{
/*
 ***refer to  zbesi,zbesk

 zwrsk computes the i bessel function for re(z) >= 0.0 by
 normalizing the i function ratios from zrati by the wronskian

 ***routines called  d1mach,zbknu,zrati,zabs
 */

    int i, nw;
    double ascle, csclr, cinui, cinur, c1i, c2i, c1r, c2r, act, acw, cti, ctr, pti, sti, ptr, str, ract;
    /* complex cinu,cscl,ct,cw,c1,c2,rct,st,y,zr */

    int verbose = 0; // FIXME -- argument

/*----------------------------------------------------------------------
     I(fnu+i-1,z) by backward recurrence for ratios
     Y(i)=I(fnu+i,z)/I(fnu+i-1,z) from crati normalized by the
     wronskian with K(fnu,z) and K(fnu+1,z) from cbknu.
 -----------------------------------------------------------------------
*/
    /* Parameter adjustments */
    --yi;
    --yr;
    --cwr;
    --cwi;

    nw = zbknu(zrr, zri, fnu, *kode, 2, verbose, &cwr[1], &cwi[1], *tol, *elim, *alim);
    if (nw != 0) { // L50:
	*nz = (nw == -2) ? -2 : -1; return;
    }
    *nz = 0;
    zrati(*zrr, *zri, *fnu, *n, *tol, &yr[1], &yi[1]);
/* -----------------------------------------------------------------------
     recur forward on I(fnu+1,z) = r(fnu,z)*I(fnu,z),
     r(fnu+j-1,z)=Y(j),  j=1,...,n
 ----------------------------------------------------------------------- */
    if (*kode == 1) {
	cinur = 1.;
	cinui = 0.;
    } else {
	cinur = cos(*zri);
	cinui = sin(*zri);
    }

/* -----------------------------------------------------------------------
     on low exponent machines the k functions can be close to both
     the under and overflow limits and the normalization must be
     scaled to prevent over or underflow. cuoik has determined that
     the result is on scale.
 ----------------------------------------------------------------------- */
    acw = zabs(cwr[2], cwi[2]);
    ascle = DBL_MIN * 1e3 / *tol;
    if (acw <= ascle) {
	csclr = 1. / *tol;
    } else { // L20:
	csclr = (acw >= (1. / ascle)) ? *tol : 1.;
    }
    // L30:
    c1r = cwr[1] * csclr;
    c1i = cwi[1] * csclr;
    c2r = cwr[2] * csclr;
    c2i = cwi[2] * csclr;
    str = yr[1];
    sti = yi[1];
/* -----------------------------------------------------------------------
     cinu=cinu*(conjg(ct)/cabs(ct))*(1.0d0/cabs(ct) prevents
     under- or overflow prematurely by squaring cabs(ct)
 ----------------------------------------------------------------------- */
    ptr = str * c1r - sti * c1i;
    pti = str * c1i + sti * c1r;
    ptr += c2r;
    pti += c2i;
    ctr = *zrr * ptr - *zri * pti;
    cti = *zrr * pti + *zri * ptr;
    act = zabs(ctr, cti);
    ract = 1. / act;
    ctr *=  ract;
    cti *= -ract;
    ptr = cinur * ract;
    pti = cinui * ract;
    cinur = ptr * ctr - pti * cti;
    cinui = ptr * cti + pti * ctr;
    yr[1] = cinur * csclr;
    yi[1] = cinui * csclr;
    for (i = 2; i <= *n; ++i) {
	ptr   = str * cinur - sti * cinui;
	cinui = str * cinui + sti * cinur;
	cinur = ptr;
	str = yr[i];
	sti = yi[i];
	yr[i] = cinur * csclr;
	yi[i] = cinui * csclr;
    }
    return;
} /* zwrsk_ */

// TODO: int function() returning current '*nz'
void
zseri_(double *zr, double *zi, double *fnu,
       int *kode, int *n, double *yr, double *yi,
       int *nz, double *tol, double *elim, double *alim)
{
/***begin prologue  zseri
 ***refer to  zbesi,zbesk

     zseri computes the i bessel function for float(z) >= 0.0 by
     means of the power series for large cabs(z) in the
     region cabs(z) <= 2*sqrt(fnu+1). nz=0 is a normal return.
     nz > 0 means that the last nz components were set to zero
     due to underflow. nz < 0 means underflow occurred, but the
     condition cabs(z) <= 2*sqrt(fnu+1) was violated and the
     computation must be completed in another routine with n=n-abs(nz).

 ***routines called  dgamln,d1mach,zuchk,zabs,zdiv,zlog_sub,zmlt
 ***end prologue  zseri
 */
    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;

    /* Local variables */
    static int i, k, l, m;
    static double s, aa;
    static int ib;
    static double ak;
    static int il;
    static double az;
    static int nn;
    static double wi[2], rs, ss;
    static double wr[2], s1i, s2i, s1r, s2r, cki, acz, arm, ckr, czi, hzi,
	     raz, czr, sti, hzr, rzi, str, rzr, ak1i, ak1r, rtr1, dfnu;
    static int idum;
    static double atol, fnup;
    static int iflag;
    static double coefi, ascle, coefr, crscr;

    /* complex ak1,ck,coef,cone,crsc,cscl,cz,czero,hz,rz,s1,s2,y,z
     Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    az = zabs(*zr, *zi);
    if (az == 0.) {
	goto L160;
    }
    arm = DBL_MIN * 1e3;
    rtr1 = sqrt(arm);
    crscr = 1.;
    iflag = 0;
    if (az < arm) {
	goto L150;
    }
    hzr = *zr * .5;
    hzi = *zi * .5;
    czr = zeror;
    czi = zeroi;
    if (az <= rtr1) {
	goto L10;
    }
    zmlt(hzr, hzi, hzr, hzi, &czr, &czi);
L10:
    acz = zabs(czr, czi);
    nn = *n;
    zlog_sub__(&hzr, &hzi, &ckr, &cki, &idum);
L20:
    dfnu = *fnu + (double) ((float) (nn - 1));
    fnup = dfnu + 1.;
/* -----------------------------------------------------------------------
     underflow test
 ----------------------------------------------------------------------- */
    ak1r = ckr * dfnu;
    ak1i = cki * dfnu;
    ak = dgamln(fnup, &idum);
    ak1r -= ak;
    if (*kode == 2) {
	ak1r -= *zr;
    }
    if (ak1r > -(*elim)) {
	goto L40;
    }
L30:
    ++(*nz);
    yr[nn] = zeror;
    yi[nn] = zeroi;
    if (acz > dfnu) {
	goto L190;
    }
    --nn;
    if (nn == 0) {
	return;
    }
    goto L20;
L40:
    if (ak1r > -(*alim)) {
	goto L50;
    }
    iflag = 1;
    ss = 1. / *tol;
    crscr = *tol;
    ascle = arm * ss;
L50:
    aa = exp(ak1r);
    if (iflag == 1) {
	aa *= ss;
    }
    coefr = aa * cos(ak1i);
    coefi = aa * sin(ak1i);
    atol = *tol * acz / fnup;
    il = imin2(2,nn);
    for (i = 1; i <= il; ++i) {
	dfnu = *fnu + (double) ((float) (nn - i));
	fnup = dfnu + 1.;
	s1r = coner;
	s1i = conei;
	if (acz < *tol * fnup) {
	    goto L70;
	}
	ak1r = coner;
	ak1i = conei;
	ak = fnup + 2.;
	s = fnup;
	aa = 2.;
L60:
	rs = 1. / s;
	str = ak1r * czr - ak1i * czi;
	sti = ak1r * czi + ak1i * czr;
	ak1r = str * rs;
	ak1i = sti * rs;
	s1r += ak1r;
	s1i += ak1i;
	s += ak;
	ak += 2.;
	aa = aa * acz * rs;
	if (aa > atol) {
	    goto L60;
	}
L70:
	s2r = s1r * coefr - s1i * coefi;
	s2i = s1r * coefi + s1i * coefr;
	wr[i - 1] = s2r;
	wi[i - 1] = s2i;
	if (iflag == 0) {
	    goto L80;
	}
	if (zuchk(s2r, s2i, ascle, *tol) != 0) {
	    goto L30;
	}
L80:
	m = nn - i + 1;
	yr[m] = s2r * crscr;
	yi[m] = s2i * crscr;
	if (i == il) {
	    goto L90;
	}
	zdiv(coefr, coefi, hzr, hzi, &str, &sti);
	coefr = str * dfnu;
	coefi = sti * dfnu;
L90:
	;
    }
    if (nn <= 2) {
	return;
    }
    k = nn - 2;
    ak = (double) ((float) k);
    raz = 1. / az;
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    if (iflag == 1) {
	goto L120;
    }
    ib = 3;
L100:
    for (i = ib; i <= nn; ++i) {
	yr[k] = (ak + *fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
	yi[k] = (ak + *fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
	ak += -1.;
	--k;
    }
    return;
/* -----------------------------------------------------------------------
     recur backward with scaled values
 ----------------------------------------------------------------------- */
L120:
/* -----------------------------------------------------------------------
     exp(-alim)=exp(-elim)/tol=approx. one precision above the
     underflow limit = ascle = d1mach(1)*ss*1.0d+3
 ----------------------------------------------------------------------- */
    s1r = wr[0];
    s1i = wi[0];
    s2r = wr[1];
    s2i = wi[1];
    for (l = 3; l <= nn; ++l) {
	ckr = s2r;
	cki = s2i;
	s2r = s1r + (ak + *fnu) * (rzr * ckr - rzi * cki);
	s2i = s1i + (ak + *fnu) * (rzr * cki + rzi * ckr);
	s1r = ckr;
	s1i = cki;
	ckr = s2r * crscr;
	cki = s2i * crscr;
	yr[k] = ckr;
	yi[k] = cki;
	ak += -1.;
	--k;
	if (zabs(ckr, cki) > ascle) {
	    goto L140;
	}
    }
    return;
L140:
    ib = l + 1;
    if (ib > nn) {
	return;
    }
    goto L100;
L150:
    *nz = *n;
    if (*fnu == 0.) {
	--(*nz);
    }
L160:
    yr[1] = zeror;
    yi[1] = zeroi;
    if (*fnu != 0.) {
	goto L170;
    }
    yr[1] = coner;
    yi[1] = conei;
L170:
    if (*n == 1) {
	return;
    }
    for (i = 2; i <= *n; ++i) {
	yr[i] = zeror;
	yi[i] = zeroi;
    }
    return;
/* -----------------------------------------------------------------------
     return with nz < 0 if cabs(z*z/4) > fnu+n-nz-1 complete
     the calculation in cbinu with n=n-iabs(nz)
 ----------------------------------------------------------------------- */
L190:
    *nz = -(*nz);
    return;
} /* zseri_ */

int zasyi(double *zr, double *zi, double *fnu,
	  int kode, int n, double *yr, double *yi,
	  double rl, double tol, double elim, double alim)
{
/***begin prologue  zasyi
 ***refer to  zbesi,zbesk

     zasyi computes the i bessel function for float(z) >= 0.0 by
     means of the asymptotic expansion for large cabs(z) in the
     region cabs(z) > max(rl,fnu*fnu/2).

     returns 'nz'
     nz = 0 is a normal return.
     nz < 0 indicates an overflow on kode=1.

 ***routines called  d1mach,zabs,zdiv,zexp_sub,zmlt,zsqrt_sub
 */

    /* Initialized data */
    static double rtpi = .159154943091895336; // = 1/(2 pi)
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;

    /* Local variables */
    static int i, j, k, m, ib, il, jl, nn, inu, koded;
    static double s, aa, bb, ak, bk, az,
	p1i, s2i, p1r, s2r, cki, dki, fdn, arg, aez, arm, ckr,
	dkr, czi, ezi, sgn,
	raz, czr, ezr, sqk, sti, rzi, str, rzr, ak1i,
	ak1r, cs1i, cs2i, cs1r, cs2r, dnu2, rtr1, dfnu, atol;
    /* complex ak1,ck,cone,cs1,cs2,cz,czero,dk,ez,p1,rz,s2,y,z */

    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    az = zabs(*zr, *zi);
    arm = DBL_MIN * 1e3;
    rtr1 = sqrt(arm);
    il = imin2(2, n);
    dfnu = *fnu + (double) ((float) (n - il));
/* -----------------------------------------------------------------------
     overflow test
 ----------------------------------------------------------------------- */
    raz = 1. / az;
    str = *zr * raz;
    sti = -(*zi) * raz;
    ak1r = rtpi * str * raz;
    ak1i = rtpi * sti * raz;
    zsqrt_sub__(&ak1r, &ak1i, &ak1r, &ak1i);
    czr = *zr;
    czi = *zi;
    if (kode == 2) {
	czr = zeror;
	czi = *zi;
    }
    if (fabs(czr) > elim) {
	// L100:
	return -1;
    }
    dnu2 = dfnu + dfnu;
    koded = 1;
    if (fabs(czr) <= alim || n <= 2) {
	koded = 0;
	zexp_sub__(&czr, &czi, &str, &sti);
	zmlt(ak1r, ak1i, str, sti, &ak1r, &ak1i);
    }
    // L20:
    fdn = 0.;
    if (dnu2 > rtr1) {
	fdn = dnu2 * dnu2;
    }
    ezr = *zr * 8.;
    ezi = *zi * 8.;
/* -----------------------------------------------------------------------
     when z is imaginary, the error test must be made relative to the
     first reciprocal power since this is the leading term of the
     expansion for the imaginary part.
 ----------------------------------------------------------------------- */
    aez = az * 8.;
    s = tol / aez;
    jl = (int) ((float) (rl + rl)) + 2;
    p1r = zeror;
    p1i = zeroi;
    if (*zi != 0.) {
	/* -----------------------------------------------------------------------
	   calculate exp(pi*(0.5+fnu+n-il)*i) to minimize losses of
	   significance when fnu or n is large
	   ----------------------------------------------------------------------- */
	inu = (int) ((float) (*fnu));
	arg = (*fnu - (double) ((float) inu)) * pi;
	inu = inu + n - il;
	ak = -sin(arg);
	bk = cos(arg);
	if (*zi < 0.) {
	    bk = -bk;
	}
	p1r = ak;
	p1i = bk;
	if (inu % 2) {
	    p1r = -p1r;
	    p1i = -p1i;
	}
    }
    // L30:
    for (k = 1; k <= il; ++k) {
	sqk = fdn - 1.;
	atol = s * fabs(sqk);
	sgn = 1.;
	cs1r = coner;
	cs1i = conei;
	cs2r = coner;
	cs2i = conei;
	ckr = coner;
	cki = conei;
	ak = 0.;
	aa = 1.;
	bb = aez;
	dkr = ezr;
	dki = ezi;
	for (j = 1; j <= jl; ++j) {
	    zdiv(ckr, cki, dkr, dki, &str, &sti);
	    ckr = str * sqk;
	    cki = sti * sqk;
	    cs2r += ckr;
	    cs2i += cki;
	    sgn = -sgn;
	    cs1r += ckr * sgn;
	    cs1i += cki * sgn;
	    dkr += ezr;
	    dki += ezi;
	    aa = aa * fabs(sqk) / bb;
	    bb += aez;
	    ak += 8.;
	    sqk -= ak;
	    if (aa <= atol) // converged:
		goto L50;
	}
	// L110:
	return -2;

L50:
	s2r = cs1r;
	s2i = cs1i;
	if (*zr + *zr < elim) {
	    double
		tzr = *zr + *zr,
		tzi = *zi + *zi,
		d__1 = -tzr,
		d__2 = -tzi;
	    zexp_sub__(&d__1, &d__2,   &str, &sti);
	    zmlt(str, sti,  p1r,  p1i, &str, &sti);
	    zmlt(str, sti, cs2r, cs2i, &str, &sti);
	    s2r += str;
	    s2i += sti;
	}
	// L60:
	fdn = fdn + dfnu * 8. + 4.;
	p1r = -p1r;
	p1i = -p1i;
	m = n - il + k;
	yr[m] = s2r * ak1r - s2i * ak1i;
	yi[m] = s2r * ak1i + s2i * ak1r;
    }
    if (n <= 2)
	return 0;
    // n >= 3 : ---
    nn = n;
    k = nn - 2;
    ak = (double) k;
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    ib = 3;
    for (i = ib; i <= nn; ++i) {
	yr[k] = (ak + *fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
	yi[k] = (ak + *fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
	ak += -1.;
	--k;
    }
    if(koded) {
	zexp_sub__(&czr, &czi, &ckr, &cki);
	for (i = 1; i <= nn; ++i) {
	    str = yr[i] * ckr - yi[i] * cki;
	    yi[i] = yr[i] * cki + yi[i] * ckr;
	    yr[i] = str;
	}
    }
    return 0;
} /* zasyi_ */

int zuoik(double *zr, double *zi, double *fnu,
	  int kode, int ikflg, int n, double *yr, double *yi,
	  double tol, double elim, double alim)
{
/***begin prologue  zuoik
 ***refer to  zbesi,zbesk,zbesh

     zuoik computes the leading terms of the uniform asymptotic
     expansions for the i and k functions and compares them
     (in logarithmic form) to alim and elim for over and underflow
     where alim < elim. if the magnitude, based on the leading
     exponential, is less than alim or greater than -alim, then
     the result is on scale. if not, then a refined test using other
     multipliers (in logarithmic form) is made based on elim. here
     exp(-elim)=smallest machine number*1.0e+3 and exp(-alim)=
     exp(-elim)/tol

     ikflg=1 means the i sequence is tested
          =2 means the k sequence is tested
     nuf = 0 means the last member of the sequence is on scale
         =-1 means an overflow would occur
     ikflg=1 and nuf > 0 means the last nuf y values were set to zero
             the first n-nuf values must be set by another routine
     ikflg=2 and nuf = n means all y values were set to zero
     ikflg=2 and 0 < nuf < n not considered. y must be set by
             another routine

 ***routines called  zuchk,zunhj,zunik,d1mach,zabs,zlog_sub
 ***end prologue  zuoik
 */
    int nuf = 0; // return value

    /* Initialized data */

    static double zeror = 0.;
    static double zeroi = 0.;

    /* Local variables */
    int i, idum, iform, nn;
    double cwrki[16], cwrkr[16];
    double aarg, aphi, argi, argr, ascle, ax, ay,
	fnn, gnn, gnu, rcz,
	czi, czr, phii, phir, sti, str, zbi, zbr, zni, znr, zri, zrr,
	asumi, asumr, bsumi, bsumr, sumi, sumr,
        zeta1i, zeta1r,	zeta2i, zeta2r;

    /*  complex arg,asum,bsum,cwrk,cz,czero,phi,sum,y,z,zb,zeta1,zeta2,zn, zr */

    nn = n;
    zrr = *zr;
    zri = *zi;
    if (*zr < 0.) {
	zrr = -(*zr);
	zri = -(*zi);
    }
    // L10:
    zbr = zrr;
    zbi = zri;
    ax = fabs(*zr) * 1.7321;
    ay = fabs(*zi);
    iform = 1;
    if (ay > ax) {
	iform = 2;
    }
    gnu = fmax2(*fnu, 1.);
    if (ikflg != 1) {
	fnn = (double) nn;
	gnn = *fnu + fnn - 1.;
	gnu = fmax2(gnn,fnn);
    }
    // L20:
/* -----------------------------------------------------------------------
     only the magnitude of arg and phi are needed along with the
     float parts of zeta1, zeta2 and zb. no attempt is made to get
     the sign of the imaginary part correct.
 ----------------------------------------------------------------------- */
    if (iform != 2) {
	int init = 0;
	zunik_(&zrr, &zri, &gnu, &ikflg, &c__1, &tol, &init,
	       &phir, &phii, &zeta1r, &zeta1i, &zeta2r, &zeta2i,
	       &sumr, &sumi, cwrkr, cwrki);
	czr = -zeta1r + zeta2r;
	czi = -zeta1i + zeta2i;
	// goto L50;
    } else { // L30:  iform == 2
	znr = zri;
	zni = -zrr;
	if (*zi <= 0.)
	    znr = -znr;
	// L40:
	zunhj_(&znr, &zni, &gnu, &c__1, &tol, &phir, &phii,
	       &argr, &argi, &zeta1r, &zeta1i, &zeta2r, &zeta2i,
	       &asumr, &asumi, &bsumr, &bsumi);
	czr = -zeta1r + zeta2r;
	czi = -zeta1i + zeta2i;
	aarg = zabs(argr, argi);
    }
    // L50:
    if (kode != 1) {
	czr -= zbr;
	czi -= zbi;
    }
    // L60:
    if (ikflg != 1) {
	czr = -czr;
	czi = -czi;
    }
    // L70:
    aphi = zabs(phir, phii);
    rcz = czr;
/* -----------------------------------------------------------------------
     overflow test
 ----------------------------------------------------------------------- */
    if (rcz > elim) {
	goto L210;
    }
    if (rcz < alim) {
	goto L80;
    }
    rcz += log(aphi);
    if (iform == 2) {
	rcz = rcz - log(aarg) * .25 - aic;
    }
    if (rcz > elim) {
	goto L210;
    }
    goto L130;

L80:
/* -----------------------------------------------------------------------
     underflow test
 ----------------------------------------------------------------------- */
    if (rcz < -elim) {
	goto L90;
    }
    if (rcz > -alim) {
	goto L130;
    }
    rcz += log(aphi);
    if (iform == 2) {
	rcz = rcz - log(aarg) * .25 - aic;
    }
    if (rcz > -elim) {
	goto L110;
    }
L90:
    for (i = 0; i < nn; ++i) {
	yr[i] = zeror;
	yi[i] = zeroi;
    }
    nuf = nn;
    return nuf;

L110:
    ascle = DBL_MIN * 1e3 / tol;
    zlog_sub__(&phir, &phii, &str, &sti, &idum);
    czr += str;
    czi += sti;
    if (iform != 1) {
	zlog_sub__(&argr, &argi, &str, &sti, &idum);
	czr = czr - str * .25 - aic;
	czi -= sti * .25;
    }
    // L120:
    ax = exp(rcz) / tol;
    ay = czi;
    czr = ax * cos(ay);
    czi = ax * sin(ay);
    if (zuchk(czr, czi, ascle, tol) != 0) {
	goto L90;
    }
L130:
    if (ikflg == 2 || n == 1) {
	return nuf;
    }
/* -----------------------------------------------------------------------
     set underflows on i sequence
 ----------------------------------------------------------------------- */
L140:
    gnu = *fnu + (double) ((float) (nn - 1));
    if (iform != 2) {
	int init = 0;
	zunik_(&zrr, &zri, &gnu, &ikflg, &c__1, &tol, &init, &phir, &phii,
	       &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
	czr = -zeta1r + zeta2r;
	czi = -zeta1i + zeta2i;
    } else { // iform == 2 : L150 :
	zunhj_(&znr, &zni, &gnu, &c__1, &tol, &phir, &phii, &argr, &argi,
	       &zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr, &bsumi);
	czr = -zeta1r + zeta2r;
	czi = -zeta1i + zeta2i;
	aarg = zabs(argr, argi);
    }
    // L160:
    if (kode != 1) {
	czr -= zbr;
	czi -= zbi;
    }
    // L170:
    aphi = zabs(phir, phii);
    rcz = czr;
    if (rcz < -elim) {
	goto L180;
    }
    if (rcz > -alim) {
	return nuf;
    }
    rcz += log(aphi);
    if (iform == 2) {
	rcz = rcz - log(aarg) * .25 - aic;
    }
    if (rcz > -elim) {
	goto L190;
    }
L180:
    --nn; // 0-indexing:
    yr[nn] = zeror;
    yi[nn] = zeroi;
    ++nuf;
    if (nn == 0) {
	return nuf;
    }
    goto L140;

L190:
    ascle = DBL_MIN * 1e3 / tol;
    zlog_sub__(&phir, &phii, &str, &sti, &idum);
    czr += str;
    czi += sti;
    if (iform != 1) {
	zlog_sub__(&argr, &argi, &str, &sti, &idum);
	czr = czr - str * .25 - aic;
	czi -= sti * .25;
    }
    // L200:
    ax = exp(rcz) / tol;
    ay = czi;
    czr = ax * cos(ay);
    czi = ax * sin(ay);
    if (zuchk(czr, czi, ascle, tol) != 0) {
	goto L180;
    }
    return nuf;

L210: return -1;

} /* zuoik_ */

void
zacon_(double *zr, double *zi, double *fnu,
       int *kode, int *mr, int *n, double *yr, double *yi,
       int *nz,
       int verbose,
       double rl, double fnul, double tol, double elim, double alim)
{
/***begin prologue  zacon
 ***refer to  zbesk,zbesh

     zacon applies the analytic continuation formula

         K(fnu,zn*exp(mp))=K(fnu,zn)*exp(-mp*fnu) - mp*I(fnu,zn)
                 mp=pi*mr*cmplx(0.0,1.0)

     to continue the k function from the right half to the left
     half z plane

 ***routines called  zbinu,zbknu,zs1s2,d1mach,zabs,zmlt
 ***end prologue  zacon
 */
    /* Initialized data */

    static double zeror = 0.;
    static double coner = 1.;

    /* Local variables */
    int i, iuf, inu, kflag, nn, nw;
    double fn, yy, c1i, c1r, c2i, c2r, c1m, as2, s1i, s1r, s2i, s2r,
	arg, cki, ckr, cpn,
	fmr, csr, azn, sgn,
	pti, ptr, spn, sti, str, zni, znr,
	rzi, rzr, sc1i, sc1r, sc2i, sc2r, cscl, cscr,
	razn,
	ascle, bscle,
	csgni, csgnr, cspni, cspnr;
    double bry[3], cyi[2], cyr[2], csrr[3], cssr[3];

    /* complex ck,cone,cscl,cscr,csgn,cspn,cy,czero,c1,c2,rz,sc1,sc2,st,
     * s1,s2,y,z,zn */

    *nz = 0;
    znr = -(*zr);
    zni = -(*zi);
    nn = *n;
    nw = zbinu(&znr, &zni, fnu, *kode, nn, yr, yi,
	       rl, fnul, tol, elim, alim);
    if (nw < 0) {
	*nz = (nw == -2) ? -2 : -1; return; // L90
    }

/* -----------------------------------------------------------------------
     analytic continuation to the left half plane for the k function
 ----------------------------------------------------------------------- */
    nn = imin2(2,*n);
    nw = zbknu(&znr, &zni, fnu, *kode, nn, verbose, cyr, cyi, tol, elim, alim);
    if (nw != 0) {
	*nz = (nw == -2) ? -2 : -1; return; // L90
    }
    s1r = cyr[0];
    s1i = cyi[0];
    fmr = (double) (*mr);
    sgn = -fsign(pi, fmr);
    csgnr = zeror;
    csgni = sgn;
    if (*kode != 1) {
	yy = -zni;
	cpn = cos(yy);
	spn = sin(yy);
	zmlt(csgnr, csgni, cpn, spn, &csgnr, &csgni);
    }

/* -----------------------------------------------------------------------
     calculate cspn=exp(fnu*pi*i) to minimize losses of significance
     when fnu is large
 ----------------------------------------------------------------------- */
    // L10:
    inu = (int) ((float) (*fnu));
    arg = (*fnu - (double) ((float) inu)) * sgn;
    cpn = cos(arg);
    spn = sin(arg);
    cspnr = cpn;
    cspni = spn;
    if (inu % 2 == 1) {
	cspnr = -cspnr;
	cspni = -cspni;
    }
    // L20:
    iuf = 0;
    c1r = s1r;
    c1i = s1i;
    c2r = yr[0];
    c2i = yi[0];
    ascle = DBL_MIN * 1e3 / tol;
    if (*kode != 1) {
	zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, &alim, &iuf);
	*nz += nw;
	sc1r = c1r;
	sc1i = c1i;
    }
    // L30:
    zmlt(cspnr, cspni, c1r, c1i, &str, &sti);
    zmlt(csgnr, csgni, c2r, c2i, &ptr, &pti);
    yr[0] = str + ptr;
    yi[0] = sti + pti;
    if (*n == 1) {
	return;
    }
    cspnr = -cspnr;
    cspni = -cspni;
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = s2r;
    c1i = s2i;
    c2r = yr[1];
    c2i = yi[1];
    if (*kode != 1) {
	zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, &alim, &iuf);
	*nz += nw;
	sc2r = c1r;
	sc2i = c1i;
    } else { // -Wall:
	sc2r = sc2i = 0.;
    }

    // L40:
    zmlt(cspnr, cspni, c1r, c1i, &str, &sti);
    zmlt(csgnr, csgni, c2r, c2i, &ptr, &pti);
    yr[1] = str + ptr;
    yi[1] = sti + pti;
    if (*n == 2) {
	return;
    }
    cspnr = -cspnr;
    cspni = -cspni;
    azn = zabs(znr, zni);
    razn = 1. / azn;
    str = znr * razn;
    sti = -zni * razn;
    rzr = (str + str) * razn;
    rzi = (sti + sti) * razn;
    fn = *fnu + 1.;
    ckr = fn * rzr;
    cki = fn * rzi;
/* -----------------------------------------------------------------------
     scale near exponent extremes during recurrence on k functions
 ----------------------------------------------------------------------- */
    cscl = 1. / tol;
    cscr = tol;
    cssr[0] = cscl;    cssr[1] = coner;    cssr[2] = cscr;
    csrr[0] = cscr;    csrr[1] = coner;    csrr[2] = cscl;
    bry [0] = ascle;   bry [1] = 1./ascle; bry [2] = DBL_MAX;
    as2 = zabs(s2r, s2i);
    if (as2 <= bry[0]) {
	kflag = 1;
    } else { // L50:
	if (as2 >= bry[1])
	    kflag = 3;
	else
	    kflag = 2;
    }
    // L60:
    bscle = bry[kflag - 1];
    s1r *= cssr[kflag - 1];
    s1i *= cssr[kflag - 1];
    s2r *= cssr[kflag - 1];
    s2i *= cssr[kflag - 1];
    csr  = csrr[kflag - 1];
    for (i = 2; i < *n; ++i) { // 0-indexing yr[], yi[]
	str = s2r;
	sti = s2i;
	s2r = ckr * str - cki * sti + s1r;
	s2i = ckr * sti + cki * str + s1i;
	s1r = str;
	s1i = sti;
	c1r = s2r * csr;
	c1i = s2i * csr;
	str = c1r;
	sti = c1i;
	c2r = yr[i];
	c2i = yi[i];
	if (*kode != 1 && iuf >= 0) {
	    zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, &alim, &iuf);
	    *nz += nw;
	    sc1r = sc2r;
	    sc1i = sc2i;
	    sc2r = c1r;
	    sc2i = c1i;
	    if (iuf == 3) {
		iuf = -4;
		s1r = sc1r * cssr[kflag - 1];
		s1i = sc1i * cssr[kflag - 1];
		s2r = sc2r * cssr[kflag - 1];
		s2i = sc2i * cssr[kflag - 1];
		str = sc2r;
		sti = sc2i;
	    }
	}
	// L70:
	ptr = cspnr * c1r - cspni * c1i;
	pti = cspnr * c1i + cspni * c1r;
	yr[i] = ptr + csgnr * c2r - csgni * c2i;
	yi[i] = pti + csgnr * c2i + csgni * c2r;
	ckr += rzr;
	cki += rzi;
	cspnr = -cspnr;
	cspni = -cspni;
	if (kflag >= 3) {
	    continue; // for(i )
	}
	ptr = fabs(c1r);
	pti = fabs(c1i);
	c1m = fmax2(ptr,pti);
	if (c1m <= bscle) {
	    continue; // for(i )
	}
	++kflag;
	bscle = bry[kflag - 1];
	s1r *= csr;
	s1i *= csr;
	s2r = str;
	s2i = sti;
	s1r *= cssr[kflag - 1];
	s1i *= cssr[kflag - 1];
	s2r *= cssr[kflag - 1];
	s2i *= cssr[kflag - 1];
	csr = csrr[kflag - 1];
    }  // for (i = 2; i < *n; ...

    return;
} /* zacon_ */

int zbinu(double *zr, double *zi, double *fnu,
	  int kode, int n, double *cyr, double *cyi,
	  double rl, double fnul, double tol, double elim, double alim)
{
/*
 ***refer to  zbesh,zbesi,zbesj,zbesk,zairy,zbiry

     zbinu computes the i function in the right half z plane

 ***routines called  zabs,zasyi,zbuni,zmlri,zseri,zuoik,zwrsk
*/
    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;

    int nz = 0; // = return value
    int i, nw, nui, inw, nlast;
    double dfnu, cwi[2], cwr[2];

    double az = zabs(*zr, *zi);
    int nn = n;
    if (az > 2.) {
	dfnu = *fnu + (double) (n - 1);
	if (az * az * .25 > dfnu + 1.) {
	    goto L20;
	}
    }
    // L10:
/* -----------------------------------------------------------------------
     power series
 ----------------------------------------------------------------------- */
    zseri_(zr, zi, fnu, &kode, &nn, cyr, cyi, &nw, &tol, &elim, &alim);
    inw = abs(nw);
    nz += inw;
    nn -= inw;
    if (nn == 0 || nw >= 0) {
	return nz;
    }
    dfnu = *fnu + (double) (nn - 1);
L20:
    if (az < rl) {
	goto L40;
    }
    if (dfnu > 1.) {
	if (az + az < dfnu * dfnu) {
	    goto L50;
	}
    }
/* -----------------------------------------------------------------------
     asymptotic expansion for large z
 ----------------------------------------------------------------------- */
    // L30:
    nw = zasyi(zr, zi, fnu, kode, nn, cyr, cyi, rl, tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    return nz;
L40:
    if (dfnu <= 1.) {
	goto L70;
    }
L50:
/* -----------------------------------------------------------------------
     overflow and underflow test on i sequence for miller algorithm
 ----------------------------------------------------------------------- */
    nw = zuoik(zr, zi, fnu, kode, 1, nn, cyr, cyi, tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    nz += nw;
    nn -= nw;
    if (nn == 0) {
	return nz;
    }
    dfnu = *fnu + (double) (nn - 1);
    if (dfnu > fnul || az > fnul) {
	goto L110;
    }

Loop_60: //=====================================================================
    if (az > rl) {
	goto L80;
    }
L70:
/* -----------------------------------------------------------------------
     miller algorithm normalized by the series
 ----------------------------------------------------------------------- */
    nw = zmlri(zr, zi, fnu, kode, nn, cyr, cyi, tol);
    if (nw < 0) {
	goto L130;
    }
    return nz;
L80:
/* -----------------------------------------------------------------------
     miller algorithm normalized by the wronskian
 -----------------------------------------------------------------------
 -----------------------------------------------------------------------
     overflow test on k functions used in wronskian
 ----------------------------------------------------------------------- */
    nw = zuoik(zr, zi, fnu, kode, 2, 2, cwr, cwi, tol, elim, alim);
    if (nw < 0) {
	nz = nn;
	for (i = 0; i < nn; ++i) {
	    cyr[i] = zeror;
	    cyi[i] = zeroi;
	}
	return nz;
    }
    // L100:
    if (nw > 0) {
	goto L130;
    }
    zwrsk_(zr, zi, fnu, &kode, &nn, cyr, cyi, &nw, cwr, cwi, &tol, &elim, &alim);
    if (nw < 0) {
	goto L130;
    }
    return nz;

L110:
/* -----------------------------------------------------------------------
     increment fnu+nn-1 up to fnul, compute and recur backward
 ----------------------------------------------------------------------- */
    nui = imax2(0, (int) (fnul - dfnu) + 1);
    nw = zbuni(zr, zi, fnu, kode, nn, nui,
	       cyr, cyi, &nlast,
	       fnul, tol, elim, alim);
    if (nw < 0) {
	goto L130;
    }
    nz += nw;
    if (nlast == 0) {
	return nz;
    }
    nn = nlast;
    goto Loop_60; //============================================================

L130: return ((nw == -2) ? -2 : -1);

} /* zbinu */


double dgamln(double z, int *ierr)
{
/***begin prologue  dgamln
***date written   830501   (yymmdd)
***revision date  830501   (yymmdd)
***category no.  b5f
***keywords  gamma function,logarithm of gamma function
***author  amos, donald e., sandia national laboratories
***purpose  to compute the logarithm of the gamma function
***description

               **** a double precision routine ****
         dgamln computes the natural log of the gamma function for
         z > 0.  the asymptotic expansion is used to generate values
         greater than zmin which are adjusted by the recursion
         g(z+1)=z*g(z) for z <= zmin.  the function was made as
         portable as possible by computimg zmin from the number of base
         10 digits in a word, rln=amax1(-alog10(r1mach(4)),0.5e-18)
         limited to 18 digits of (relative) accuracy.

         since int arguments are common, a table look up on 100
         values is used for speed of execution.

     description of arguments

         input      z is d0uble precision
           z      - argument, z > 0

         output      dgamln is double precision
           dgamln  - natural log of the gamma function at z != 0
           ierr    - error flag
                     ierr=0, normal return, computation completed
                     ierr=1, z <= 0,    no computation


*** references  computation of bessel functions of complex argument
                 by d. e. amos, sand83-0083, may, 1983.

***end prologue  dgamln
*/

    /* Initialized data */

    /* ln(gamma(n)), n=1,100 : */
    static const double gln[100] = { 0.,0.,.693147180559945309,
	    1.791759469228055,3.17805383034794562,4.78749174278204599,
	    6.579251212010101,8.5251613610654143,10.6046029027452502,
	    12.8018274800814696,15.1044125730755153,17.5023078458738858,
	    19.9872144956618861,22.5521638531234229,25.1912211827386815,
	    27.8992713838408916,30.6718601060806728,33.5050734501368889,
	    36.3954452080330536,39.339884187199494,42.335616460753485,
	    45.380138898476908,48.4711813518352239,51.6066755677643736,
	    54.7847293981123192,58.0036052229805199,61.261701761002002,
	    64.5575386270063311,67.889743137181535,71.257038967168009,
	    74.6582363488301644,78.0922235533153106,81.5579594561150372,
	    85.0544670175815174,88.5808275421976788,92.1361756036870925,
	    95.7196945421432025,99.3306124547874269,102.968198614513813,
	    106.631760260643459,110.320639714757395,114.034211781461703,
	    117.771881399745072,121.533081515438634,125.317271149356895,
	    129.123933639127215,132.95257503561631,136.802722637326368,
	    140.673923648234259,144.565743946344886,148.477766951773032,
	    152.409592584497358,156.360836303078785,160.331128216630907,
	    164.320112263195181,168.327445448427652,172.352797139162802,
	    176.395848406997352,180.456291417543771,184.533828861449491,
	    188.628173423671591,192.739047287844902,196.866181672889994,
	    201.009316399281527,205.168199482641199,209.342586752536836,
	    213.532241494563261,217.736934113954227,221.956441819130334,
	    226.190548323727593,230.439043565776952,234.701723442818268,
	    238.978389561834323,243.268849002982714,247.572914096186884,
	    251.890402209723194,256.221135550009525,260.564940971863209,
	    264.921649798552801,269.291097651019823,273.673124285693704,
	    278.067573440366143,282.474292687630396,286.893133295426994,
	    291.323950094270308,295.766601350760624,300.220948647014132,
	    304.686856765668715,309.164193580146922,313.652829949879062,
	    318.152639620209327,322.663499126726177,327.185287703775217,
	    331.717887196928473,336.261181979198477,340.815058870799018,
	    345.379407062266854,349.954118040770237,354.539085519440809,
	    359.134205369575399 };
    /* coefficients of asymptotic expansion : */
    static const double cf[22] = { .0833333333333333333,-.00277777777777777778,
	    7.93650793650793651e-4,-5.95238095238095238e-4,
	    8.41750841750841751e-4,-.00191752691752691753,
	    .00641025641025641026,-.0295506535947712418,.179644372368830573,
	    -1.39243221690590112,13.402864044168392,-156.848284626002017,
	    2193.10333333333333,-36108.7712537249894,691472.268851313067,
	    -15238221.5394074162,382900751.391414141,-10882266035.7843911,
	    347320283765.002252,-12369602142269.2745,488788064793079.335,
	    -21320333960919373.9 };

    static const double con = 1.83787706640934548; /* = ln(2*pi) */

    /* Local variables */
    int k, nz = -9999/*Wall*/;
    double s, t1, fz, zp, tlg, trm, tst, zsq, zinc, zdmy;

    *ierr = 0;
    if (z <= 0.) { *ierr = 1; return 0.; }

    if (z <= 101.) {
	nz = (int) (z);
	fz = z - (double) nz;
	if (fz <= 0. && nz <= 100) {
	    return gln[nz - 1];
	}
    }
    /* L10: */
    double
	wdtol = fmax2(DBL_EPSILON, 5e-19),
	rln = M_LOG10_2 * (double) DBL_MANT_DIG,
	fln = fmin2(rln, 20.);
    fln = fmax2(fln, 3.) - 3.; /* => fln in [0, 17] */
    double zm = fln * .3875 + 1.8;   /*  zm in [1.8, 8.3875] */
    int mz = (int) zm + 1; /* in {2, 3,.., 9} */
    double zmin = (double) mz;
    if (z < zmin) { /* z < zmin <= 9  --> nz was defined above */
	zinc = zmin - (double) nz;
	zdmy = z + zinc;
    } else {
	zinc = 0.;
	zdmy = z;
    }
    zp = 1. / zdmy;
    t1 = cf[0] * zp;
    s = t1;
    if (zp >= wdtol) {
	zsq = zp * zp;
	tst = t1 * wdtol;
	for (k = 2; k <= 22; ++k) {
	    zp *= zsq;
	    trm = cf[k - 1] * zp;
	    if (fabs(trm) < tst)
		break;
	    s += trm;
	}
    }
    /* L40: */
    if (zinc == 0.) {
	tlg = log(z);
	return z * (tlg - 1.) + (con - tlg) * .5 + s;
    } else { /* L50: */
	zp = 1.;
	for (int i = 0; i < (int) zinc; ++i)
	    zp *= z + (double) (i);

	tlg = log(zdmy);
	return zdmy * (tlg - 1.) - log(zp) + (con - tlg) * .5 + s;
    }
} /* dgamln */

void
zacai_(double *zr, double *zi, double *fnu,
       int *kode, int *mr, int *n, double *yr, double *yi,
       int *nz, double *rl, double *tol, double *elim, double *alim)
{
/***begin prologue  zacai
 ***refer to  zairy

     zacai applies the analytic continuation formula

         K(fnu,zn*exp(mp))=K(fnu,zn)*exp(-mp*fnu) - mp*I(fnu,zn)
                 mp=pi*mr*cmplx(0.0,1.0)

     to continue the k function from the right half to the left
     half z plane for use with zairy where fnu=1/3 or 2/3 and n=1.
     zacai is the same as zacon with the parts for larger orders and
     recurrence removed. a recursive call to zacon can result if zacon
     is called from zairy.

 ***routines called  zasyi,zbknu,zmlri,zseri,zs1s2,d1mach,zabs
 ***end prologue  zacai
 */

    int verbose = 0; // FIXME -- argument

    /* Local variables */
    int inu, nn, nw;
    double arg, az, fmr, sgn, yy, c1i, c1r, c2i, c2r,
	cyi[2], cyr[2],
	zni, znr, dfnu, ascle, csgni, csgnr, cspni, cspnr;

    /*  complex csgn,cspn,c1,c2,y,z,zn,cy
     */

    /* Function Body */
    *nz = 0;
    znr = -(*zr);
    zni = -(*zi);
    az = zabs(*zr, *zi);
    nn = *n;
    dfnu = *fnu + (double) ((float) (*n - 1));
    if (az <= 2.) { // L10:
	/* -----------------------------------------------------------------------
	   power series for the i function
	   ----------------------------------------------------------------------- */
	zseri_(&znr, &zni, fnu, kode, &nn, yr, yi, &nw, tol, elim, alim);

    } else if ((az * az * .25 > dfnu + 1.)  && (az >= *rl)) {
	/* -----------------------------------------------------------------------
	   asymptotic expansion for large z for the i function
	   ----------------------------------------------------------------------- */
	nw = zasyi(&znr, &zni, fnu, *kode, nn, yr, yi, *rl, *tol, *elim, *alim);
	if (nw < 0) {
	    goto L80;
	}

    } else { // L30:
	/* -----------------------------------------------------------------------
	   miller algorithm normalized by the series for the i function
	   ----------------------------------------------------------------------- */
	nw = zmlri(&znr, &zni, fnu, *kode, nn, yr, yi, *tol);
	if (nw < 0) {
	    goto L80;
	}
    }

    // L40:
/* -----------------------------------------------------------------------
     analytic continuation to the left half plane for the k function
 ----------------------------------------------------------------------- */
    nw = zbknu(&znr, &zni, fnu, *kode, 1, verbose, cyr, cyi, *tol, *elim, *alim);
    if (nw != 0) {
	goto L80;
    }
    fmr = (double) (*mr);
    sgn = -fsign(pi, fmr);
    csgnr = 0.;
    csgni = sgn;
    if (*kode != 1) {
	yy = -zni;
	csgnr = -csgni * sin(yy);
	csgni *= cos(yy);
    }
    // L50:
/* -----------------------------------------------------------------------
     calculate cspn=exp(fnu*pi*i) to minimize losses of significance
     when fnu is large
 ----------------------------------------------------------------------- */
    inu = (int) *fnu;
    arg = (*fnu - (double) inu) * sgn;
    cspnr = cos(arg);
    cspni = sin(arg);
    if (inu % 2 == 1) {
	cspnr = -cspnr;
	cspni = -cspni;
    }
    // L60:
    c1r = cyr[0];
    c1i = cyi[0];
    c2r = yr[0];
    c2i = yi[0];
    if (*kode != 1) {
	int iuf = 0;
	ascle = DBL_MIN * 1e3 / *tol;
	zs1s2_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim, &iuf);
	*nz += nw;
    }
    // L70:
    yr[0] = cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i;
    yi[0] = cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r;
    return;

L80: *nz = (nw == -2) ? -2 : -1; return;

} /* zacai_ */

int zuchk(double yr, double yi, double ascle, double tol)
{
/* ***refer to zseri,zuoik,zunk1,zunk2,zuni1,zuni2,zkscl

      y (complex, = yr + i * yi) enters as a scaled quantity whose magnitude
        is greater than exp(-alim) =: ascle = 1.0e+3*d1mach(1)/tol.

   The test is made to see
   if the magnitude of the float or imaginary part would underflow
   when y is scaled (by tol) to its proper value.  y is accepted
   if the underflow is at least one precision below the magnitude
   of the largest component; otherwise the phase angle does not have
   absolute accuracy and an underflow is assumed.
 */
    double
	wr = fabs(yr),
	wi = fabs(yi),
	st = fmin2(wr,wi);
    if (st > ascle) {
	return 0;
    } else {
	double ss = fmax2(wr,wi);
	st /= tol;
	return (ss < st) ? 1 : 0;
    }
} /* zuchk */

/* Subroutine */ int
zunik_(double *zrr, double *zri, double *fnu,
       int *ikflg, int *ipmtr, double *tol, int *init,
       double *phir, double *phii, double *zeta1r, double *zeta1i,
       double *zeta2r, double *zeta2i, double *sumr, double *sumi,
       double *cwrkr, double *cwrki)
{
/***begin prologue  zunik
 ***refer to  zbesi,zbesk

        zunik computes parameters for the uniform asymptotic
        expansions of the i and k functions on ikflg= 1 or 2
        respectively by

        w(fnu,zr) = phi*exp(zeta)*sum

        where       zeta=-zeta1 + zeta2       or
                          zeta1 - zeta2

        the first call must have init=0. subsequent calls with the
        same zr and fnu will return the i or k function on ikflg=
        1 or 2 with no change in init. cwrk is a complex work
        array. ipmtr=0 computes all parameters. ipmtr=1 computes phi,
        zeta1,zeta2.

 ***routines called  zdiv,zlog_sub,zsqrt_sub,d1mach
 ***end prologue  zunik
 */
    /* Initialized data */

    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;
    static double con[2] = { .398942280401432678,1.25331413731550025 };
    static double c__[120] = { 1.,-.208333333333333333,.125,
	    .334201388888888889,-.401041666666666667,.0703125,
	    -1.02581259645061728,1.84646267361111111,-.8912109375,.0732421875,
	    4.66958442342624743,-11.2070026162229938,8.78912353515625,
	    -2.3640869140625,.112152099609375,-28.2120725582002449,
	    84.6362176746007346,-91.8182415432400174,42.5349987453884549,
	    -7.3687943594796317,.227108001708984375,212.570130039217123,
	    -765.252468141181642,1059.99045252799988,-699.579627376132541,
	    218.19051174421159,-26.4914304869515555,.572501420974731445,
	    -1919.457662318407,8061.72218173730938,-13586.5500064341374,
	    11655.3933368645332,-5305.64697861340311,1200.90291321635246,
	    -108.090919788394656,1.7277275025844574,20204.2913309661486,
	    -96980.5983886375135,192547.001232531532,-203400.177280415534,
	    122200.46498301746,-41192.6549688975513,7109.51430248936372,
	    -493.915304773088012,6.07404200127348304,-242919.187900551333,
	    1311763.6146629772,-2998015.91853810675,3763271.297656404,
	    -2813563.22658653411,1268365.27332162478,-331645.172484563578,
	    45218.7689813627263,-2499.83048181120962,24.3805296995560639,
	    3284469.85307203782,-19706819.1184322269,50952602.4926646422,
	    -74105148.2115326577,66344512.2747290267,-37567176.6607633513,
	    13288767.1664218183,-2785618.12808645469,308186.404612662398,
	    -13886.0897537170405,110.017140269246738,-49329253.664509962,
	    325573074.185765749,-939462359.681578403,1553596899.57058006,
	    -1621080552.10833708,1106842816.82301447,-495889784.275030309,
	    142062907.797533095,-24474062.7257387285,2243768.17792244943,
	    -84005.4336030240853,551.335896122020586,814789096.118312115,
	    -5866481492.05184723,18688207509.2958249,-34632043388.1587779,
	    41280185579.753974,-33026599749.8007231,17954213731.1556001,
	    -6563293792.61928433,1559279864.87925751,-225105661.889415278,
	    17395107.5539781645,-549842.327572288687,3038.09051092238427,
	    -14679261247.6956167,114498237732.02581,-399096175224.466498,
	    819218669548.577329,-1098375156081.22331,1008158106865.38209,
	    -645364869245.376503,287900649906.150589,-87867072178.0232657,
	    17634730606.8349694,-2167164983.22379509,143157876.718888981,
	    -3871833.44257261262,18257.7554742931747,286464035717.679043,
	    -2406297900028.50396,9109341185239.89896,-20516899410934.4374,
	    30565125519935.3206,-31667088584785.1584,23348364044581.8409,
	    -12320491305598.2872,4612725780849.13197,-1196552880196.1816,
	    205914503232.410016,-21822927757.5292237,1247009293.51271032,
	    -29188388.1222208134,118838.426256783253 };

    /* Local variables */
    static int i, j, k, l;
    static double ac;
    static double si, ti, sr, tr, t2i, t2r, rfn, sri, sti, zni, srr, str,
	    znr;
    static int idum;
    static double test, crfni, crfnr;
    /* complex cfn,con,cone,crfn,cwrk,czero,phi,s,sr,sum,t,t2,zeta1,
     * zeta2,zn,zr */

    /* Parameter adjustments */
    --cwrki;
    --cwrkr;

    /* Function Body */

    if (*init != 0) {
	goto L40;
    }
/* -----------------------------------------------------------------------
     initialize all variables
 ----------------------------------------------------------------------- */
    rfn = 1. / *fnu;
/* -----------------------------------------------------------------------
     overflow test (zr/fnu too small)
 ----------------------------------------------------------------------- */
    test = DBL_MIN * 1e3;
    ac = *fnu * test;
    if (fabs(*zrr) <= ac && fabs(*zri) <= ac) {
	*zeta1r = fabs(log(test)) * 2. + *fnu;
	*zeta1i = 0.;
	*zeta2r = *fnu; *zeta2i = 0.;
	*phir = 1.;	*phii = 0.;
	return 0;
    }
    // L15:
    tr = *zrr * rfn;
    ti = *zri * rfn;
    sr = coner + (tr * tr - ti * ti);
    si = conei + (tr * ti + ti * tr);
    zsqrt_sub__(&sr, &si, &srr, &sri);
    str = coner + srr;
    sti = conei + sri;
    zdiv(str, sti, tr, ti, &znr, &zni);
    zlog_sub__(&znr, &zni, &str, &sti, &idum);
    *zeta1r = *fnu * str;
    *zeta1i = *fnu * sti;
    *zeta2r = *fnu * srr;
    *zeta2i = *fnu * sri;
    zdiv(coner, conei, srr, sri, &tr, &ti);
    srr = tr * rfn;
    sri = ti * rfn;
    zsqrt_sub__(&srr, &sri, &cwrkr[16], &cwrki[16]);
    *phir = cwrkr[16] * con[*ikflg - 1];
    *phii = cwrki[16] * con[*ikflg - 1];
    if (*ipmtr != 0) {
	return 0;
    }
    zdiv(coner, conei, sr, si, &t2r, &t2i);
    cwrkr[1] = coner;
    cwrki[1] = conei;
    crfnr = coner;
    crfni = conei;
    ac = 1.;
    l = 1;
    for (k = 2; k <= 15; ++k) {
	sr = zeror;
	si = zeroi;
	for (j = 1; j <= k; ++j) {
	    ++l;
	    str = sr * t2r - si * t2i + c__[l - 1];
	    si = sr * t2i + si * t2r;
	    sr = str;
	}
	str = crfnr * srr - crfni * sri;
	crfni = crfnr * sri + crfni * srr;
	crfnr = str;
	cwrkr[k] = crfnr * sr - crfni * si;
	cwrki[k] = crfnr * si + crfni * sr;
	ac *= rfn;
	test = fabs(cwrkr[k]) + fabs(cwrki[k]);
	if (ac < *tol && test < *tol) {
	    goto L30;
	}
    }
    k = 15;
L30:
    *init = k;
L40:
    if (*ikflg == 2) {
	goto L60;
    }
/* -----------------------------------------------------------------------
     compute sum for the i function
 ----------------------------------------------------------------------- */
    sr = zeror;
    si = zeroi;
    for (i = 1; i <= *init; ++i) {
	sr += cwrkr[i];
	si += cwrki[i];
    }
    *sumr = sr;
    *sumi = si;
    *phir = cwrkr[16] * con[0];
    *phii = cwrki[16] * con[0];
    return 0;
L60:
/* -----------------------------------------------------------------------
     compute sum for the k function
 ----------------------------------------------------------------------- */
    sr = zeror;
    si = zeroi;
    tr = coner;
    for (i = 1; i <= *init; ++i) {
	sr += tr * cwrkr[i];
	si += tr * cwrki[i];
	tr = -tr;
    }
    *sumr = sr;
    *sumi = si;
    *phir = cwrkr[16] * con[1];
    *phii = cwrki[16] * con[1];
    return 0;
} /* zunik_ */

/* Subroutine */ int
zunhj_(double *zr, double *zi, double *fnu,
       int *ipmtr, double *tol, double *phir, double *phii,
       double *argr, double *argi, double *zeta1r, double *zeta1i,
       double *zeta2r, double *zeta2i, double *asumr,
       double *asumi, double *bsumr, double *bsumi)
{
/***begin prologue  zunhj
 ***refer to  zbesi,zbesk

     references
         handbook of mathematical functions by m. abramowitz and i.a.
         stegun, ams55, national bureau of standards, 1965, chapter 9.

         asymptotics and special functions by f.w.j. olver, academic
         press, n.y., 1974, page 420

     abstract
         zunhj computes parameters for bessel functions c(fnu,z) =
         J(fnu,z), Y(fnu,z) or H(i,fnu,z) i=1,2 for large orders fnu
         by means of the uniform asymptotic expansion

         c(fnu,z)=c1*phi*( asum*airy(arg) + c2*bsum*dairy(arg) )

         for proper choices of c1, c2, airy and dairy where airy is
         an airy function and dairy is its derivative.

               (2/3)*fnu*zeta**1.5 = zeta1-zeta2,

         zeta1=0.5*fnu*clog((1+w)/(1-w)), zeta2=fnu*w for scaling
         purposes in airy functions from cairy or cbiry.

         mconj=sign of aimag(z), but is ambiguous when z is float and
         must be specified. ipmtr=0 returns all parameters. ipmtr=
         1 computes all except asum and bsum.

 ***routines called  zabs,zdiv,zlog_sub,zsqrt_sub,d1mach
 ***end prologue  zunhj
 */
    /* Initialized data */

    static double ar[14] = { 1.,.104166666666666667,.0835503472222222222,
	    .12822657455632716,.291849026464140464,.881627267443757652,
	    3.32140828186276754,14.9957629868625547,78.9230130115865181,
	    474.451538868264323,3207.49009089066193,24086.5496408740049,
	    198923.119169509794,1791902.00777534383 };
    static double br[14] = { 1.,-.145833333333333333,
	    -.0987413194444444444,-.143312053915895062,-.317227202678413548,
	    -.942429147957120249,-3.51120304082635426,-15.7272636203680451,
	    -82.2814390971859444,-492.355370523670524,-3316.21856854797251,
	    -24827.6742452085896,-204526.587315129788,-1838444.9170682099 };
    static double c__[105] = { 1.,-.208333333333333333,.125,
	    .334201388888888889,-.401041666666666667,.0703125,
	    -1.02581259645061728,1.84646267361111111,-.8912109375,.0732421875,
	    4.66958442342624743,-11.2070026162229938,8.78912353515625,
	    -2.3640869140625,.112152099609375,-28.2120725582002449,
	    84.6362176746007346,-91.8182415432400174,42.5349987453884549,
	    -7.3687943594796317,.227108001708984375,212.570130039217123,
	    -765.252468141181642,1059.99045252799988,-699.579627376132541,
	    218.19051174421159,-26.4914304869515555,.572501420974731445,
	    -1919.457662318407,8061.72218173730938,-13586.5500064341374,
	    11655.3933368645332,-5305.64697861340311,1200.90291321635246,
	    -108.090919788394656,1.7277275025844574,20204.2913309661486,
	    -96980.5983886375135,192547.001232531532,-203400.177280415534,
	    122200.46498301746,-41192.6549688975513,7109.51430248936372,
	    -493.915304773088012,6.07404200127348304,-242919.187900551333,
	    1311763.6146629772,-2998015.91853810675,3763271.297656404,
	    -2813563.22658653411,1268365.27332162478,-331645.172484563578,
	    45218.7689813627263,-2499.83048181120962,24.3805296995560639,
	    3284469.85307203782,-19706819.1184322269,50952602.4926646422,
	    -74105148.2115326577,66344512.2747290267,-37567176.6607633513,
	    13288767.1664218183,-2785618.12808645469,308186.404612662398,
	    -13886.0897537170405,110.017140269246738,-49329253.664509962,
	    325573074.185765749,-939462359.681578403,1553596899.57058006,
	    -1621080552.10833708,1106842816.82301447,-495889784.275030309,
	    142062907.797533095,-24474062.7257387285,2243768.17792244943,
	    -84005.4336030240853,551.335896122020586,814789096.118312115,
	    -5866481492.05184723,18688207509.2958249,-34632043388.1587779,
	    41280185579.753974,-33026599749.8007231,17954213731.1556001,
	    -6563293792.61928433,1559279864.87925751,-225105661.889415278,
	    17395107.5539781645,-549842.327572288687,3038.09051092238427,
	    -14679261247.6956167,114498237732.02581,-399096175224.466498,
	    819218669548.577329,-1098375156081.22331,1008158106865.38209,
	    -645364869245.376503,287900649906.150589,-87867072178.0232657,
	    17634730606.8349694,-2167164983.22379509,143157876.718888981,
	    -3871833.44257261262,18257.7554742931747 };
    static double alfa[180] = { -.00444444444444444444,
	    -9.22077922077922078e-4,-8.84892884892884893e-5,
	    1.65927687832449737e-4,2.4669137274179291e-4,
	    2.6599558934625478e-4,2.61824297061500945e-4,
	    2.48730437344655609e-4,2.32721040083232098e-4,
	    2.16362485712365082e-4,2.00738858762752355e-4,
	    1.86267636637545172e-4,1.73060775917876493e-4,
	    1.61091705929015752e-4,1.50274774160908134e-4,
	    1.40503497391269794e-4,1.31668816545922806e-4,
	    1.23667445598253261e-4,1.16405271474737902e-4,
	    1.09798298372713369e-4,1.03772410422992823e-4,
	    9.82626078369363448e-5,9.32120517249503256e-5,
	    8.85710852478711718e-5,8.42963105715700223e-5,
	    8.03497548407791151e-5,7.66981345359207388e-5,
	    7.33122157481777809e-5,7.01662625163141333e-5,
	    6.72375633790160292e-5,6.93735541354588974e-4,
	    2.32241745182921654e-4,-1.41986273556691197e-5,
	    -1.1644493167204864e-4,-1.50803558053048762e-4,
	    -1.55121924918096223e-4,-1.46809756646465549e-4,
	    -1.33815503867491367e-4,-1.19744975684254051e-4,
	    -1.0618431920797402e-4,-9.37699549891194492e-5,
	    -8.26923045588193274e-5,-7.29374348155221211e-5,
	    -6.44042357721016283e-5,-5.69611566009369048e-5,
	    -5.04731044303561628e-5,-4.48134868008882786e-5,
	    -3.98688727717598864e-5,-3.55400532972042498e-5,
	    -3.1741425660902248e-5,-2.83996793904174811e-5,
	    -2.54522720634870566e-5,-2.28459297164724555e-5,
	    -2.05352753106480604e-5,-1.84816217627666085e-5,
	    -1.66519330021393806e-5,-1.50179412980119482e-5,
	    -1.35554031379040526e-5,-1.22434746473858131e-5,
	    -1.10641884811308169e-5,-3.54211971457743841e-4,
	    -1.56161263945159416e-4,3.0446550359493641e-5,
	    1.30198655773242693e-4,1.67471106699712269e-4,
	    1.70222587683592569e-4,1.56501427608594704e-4,
	    1.3633917097744512e-4,1.14886692029825128e-4,
	    9.45869093034688111e-5,7.64498419250898258e-5,
	    6.07570334965197354e-5,4.74394299290508799e-5,
	    3.62757512005344297e-5,2.69939714979224901e-5,
	    1.93210938247939253e-5,1.30056674793963203e-5,
	    7.82620866744496661e-6,3.59257485819351583e-6,
	    1.44040049814251817e-7,-2.65396769697939116e-6,
	    -4.9134686709848591e-6,-6.72739296091248287e-6,
	    -8.17269379678657923e-6,-9.31304715093561232e-6,
	    -1.02011418798016441e-5,-1.0880596251059288e-5,
	    -1.13875481509603555e-5,-1.17519675674556414e-5,
	    -1.19987364870944141e-5,3.78194199201772914e-4,
	    2.02471952761816167e-4,-6.37938506318862408e-5,
	    -2.38598230603005903e-4,-3.10916256027361568e-4,
	    -3.13680115247576316e-4,-2.78950273791323387e-4,
	    -2.28564082619141374e-4,-1.75245280340846749e-4,
	    -1.25544063060690348e-4,-8.22982872820208365e-5,
	    -4.62860730588116458e-5,-1.72334302366962267e-5,
	    5.60690482304602267e-6,2.313954431482868e-5,
	    3.62642745856793957e-5,4.58006124490188752e-5,
	    5.2459529495911405e-5,5.68396208545815266e-5,
	    5.94349820393104052e-5,6.06478527578421742e-5,
	    6.08023907788436497e-5,6.01577894539460388e-5,
	    5.891996573446985e-5,5.72515823777593053e-5,
	    5.52804375585852577e-5,5.3106377380288017e-5,
	    5.08069302012325706e-5,4.84418647620094842e-5,
	    4.6056858160747537e-5,-6.91141397288294174e-4,
	    -4.29976633058871912e-4,1.83067735980039018e-4,
	    6.60088147542014144e-4,8.75964969951185931e-4,
	    8.77335235958235514e-4,7.49369585378990637e-4,
	    5.63832329756980918e-4,3.68059319971443156e-4,
	    1.88464535514455599e-4,3.70663057664904149e-5,
	    -8.28520220232137023e-5,-1.72751952869172998e-4,
	    -2.36314873605872983e-4,-2.77966150694906658e-4,
	    -3.02079514155456919e-4,-3.12594712643820127e-4,
	    -3.12872558758067163e-4,-3.05678038466324377e-4,
	    -2.93226470614557331e-4,-2.77255655582934777e-4,
	    -2.59103928467031709e-4,-2.39784014396480342e-4,
	    -2.20048260045422848e-4,-2.00443911094971498e-4,
	    -1.81358692210970687e-4,-1.63057674478657464e-4,
	    -1.45712672175205844e-4,-1.29425421983924587e-4,
	    -1.14245691942445952e-4,.00192821964248775885,
	    .00135592576302022234,-7.17858090421302995e-4,
	    -.00258084802575270346,-.00349271130826168475,
	    -.00346986299340960628,-.00282285233351310182,
	    -.00188103076404891354,-8.895317183839476e-4,
	    3.87912102631035228e-6,7.28688540119691412e-4,
	    .00126566373053457758,.00162518158372674427,.00183203153216373172,
	    .00191588388990527909,.00190588846755546138,.00182798982421825727,
	    .0017038950642112153,.00155097127171097686,.00138261421852276159,
	    .00120881424230064774,.00103676532638344962,
	    8.71437918068619115e-4,7.16080155297701002e-4,
	    5.72637002558129372e-4,4.42089819465802277e-4,
	    3.24724948503090564e-4,2.20342042730246599e-4,
	    1.28412898401353882e-4,4.82005924552095464e-5 };
    static double beta[210] = { .0179988721413553309,
	    .00559964911064388073,.00288501402231132779,.00180096606761053941,
	    .00124753110589199202,9.22878876572938311e-4,
	    7.14430421727287357e-4,5.71787281789704872e-4,
	    4.69431007606481533e-4,3.93232835462916638e-4,
	    3.34818889318297664e-4,2.88952148495751517e-4,
	    2.52211615549573284e-4,2.22280580798883327e-4,
	    1.97541838033062524e-4,1.76836855019718004e-4,
	    1.59316899661821081e-4,1.44347930197333986e-4,
	    1.31448068119965379e-4,1.20245444949302884e-4,
	    1.10449144504599392e-4,1.01828770740567258e-4,
	    9.41998224204237509e-5,8.74130545753834437e-5,
	    8.13466262162801467e-5,7.59002269646219339e-5,
	    7.09906300634153481e-5,6.65482874842468183e-5,
	    6.25146958969275078e-5,5.88403394426251749e-5,
	    -.00149282953213429172,-8.78204709546389328e-4,
	    -5.02916549572034614e-4,-2.94822138512746025e-4,
	    -1.75463996970782828e-4,-1.04008550460816434e-4,
	    -5.96141953046457895e-5,-3.1203892907609834e-5,
	    -1.26089735980230047e-5,-2.42892608575730389e-7,
	    8.05996165414273571e-6,1.36507009262147391e-5,
	    1.73964125472926261e-5,1.9867297884213378e-5,
	    2.14463263790822639e-5,2.23954659232456514e-5,
	    2.28967783814712629e-5,2.30785389811177817e-5,
	    2.30321976080909144e-5,2.28236073720348722e-5,
	    2.25005881105292418e-5,2.20981015361991429e-5,
	    2.16418427448103905e-5,2.11507649256220843e-5,
	    2.06388749782170737e-5,2.01165241997081666e-5,
	    1.95913450141179244e-5,1.9068936791043674e-5,
	    1.85533719641636667e-5,1.80475722259674218e-5,
	    5.5221307672129279e-4,4.47932581552384646e-4,
	    2.79520653992020589e-4,1.52468156198446602e-4,
	    6.93271105657043598e-5,1.76258683069991397e-5,
	    -1.35744996343269136e-5,-3.17972413350427135e-5,
	    -4.18861861696693365e-5,-4.69004889379141029e-5,
	    -4.87665447413787352e-5,-4.87010031186735069e-5,
	    -4.74755620890086638e-5,-4.55813058138628452e-5,
	    -4.33309644511266036e-5,-4.09230193157750364e-5,
	    -3.84822638603221274e-5,-3.60857167535410501e-5,
	    -3.37793306123367417e-5,-3.15888560772109621e-5,
	    -2.95269561750807315e-5,-2.75978914828335759e-5,
	    -2.58006174666883713e-5,-2.413083567612802e-5,
	    -2.25823509518346033e-5,-2.11479656768912971e-5,
	    -1.98200638885294927e-5,-1.85909870801065077e-5,
	    -1.74532699844210224e-5,-1.63997823854497997e-5,
	    -4.74617796559959808e-4,-4.77864567147321487e-4,
	    -3.20390228067037603e-4,-1.61105016119962282e-4,
	    -4.25778101285435204e-5,3.44571294294967503e-5,
	    7.97092684075674924e-5,1.031382367082722e-4,
	    1.12466775262204158e-4,1.13103642108481389e-4,
	    1.08651634848774268e-4,1.01437951597661973e-4,
	    9.29298396593363896e-5,8.40293133016089978e-5,
	    7.52727991349134062e-5,6.69632521975730872e-5,
	    5.92564547323194704e-5,5.22169308826975567e-5,
	    4.58539485165360646e-5,4.01445513891486808e-5,
	    3.50481730031328081e-5,3.05157995034346659e-5,
	    2.64956119950516039e-5,2.29363633690998152e-5,
	    1.97893056664021636e-5,1.70091984636412623e-5,
	    1.45547428261524004e-5,1.23886640995878413e-5,
	    1.04775876076583236e-5,8.79179954978479373e-6,
	    7.36465810572578444e-4,8.72790805146193976e-4,
	    6.22614862573135066e-4,2.85998154194304147e-4,
	    3.84737672879366102e-6,-1.87906003636971558e-4,
	    -2.97603646594554535e-4,-3.45998126832656348e-4,
	    -3.53382470916037712e-4,-3.35715635775048757e-4,
	    -3.04321124789039809e-4,-2.66722723047612821e-4,
	    -2.27654214122819527e-4,-1.89922611854562356e-4,
	    -1.5505891859909387e-4,-1.2377824076187363e-4,
	    -9.62926147717644187e-5,-7.25178327714425337e-5,
	    -5.22070028895633801e-5,-3.50347750511900522e-5,
	    -2.06489761035551757e-5,-8.70106096849767054e-6,
	    1.1369868667510029e-6,9.16426474122778849e-6,
	    1.5647778542887262e-5,2.08223629482466847e-5,
	    2.48923381004595156e-5,2.80340509574146325e-5,
	    3.03987774629861915e-5,3.21156731406700616e-5,
	    -.00180182191963885708,-.00243402962938042533,
	    -.00183422663549856802,-7.62204596354009765e-4,
	    2.39079475256927218e-4,9.49266117176881141e-4,
	    .00134467449701540359,.00148457495259449178,.00144732339830617591,
	    .00130268261285657186,.00110351597375642682,
	    8.86047440419791759e-4,6.73073208165665473e-4,
	    4.77603872856582378e-4,3.05991926358789362e-4,
	    1.6031569459472163e-4,4.00749555270613286e-5,
	    -5.66607461635251611e-5,-1.32506186772982638e-4,
	    -1.90296187989614057e-4,-2.32811450376937408e-4,
	    -2.62628811464668841e-4,-2.82050469867598672e-4,
	    -2.93081563192861167e-4,-2.97435962176316616e-4,
	    -2.96557334239348078e-4,-2.91647363312090861e-4,
	    -2.83696203837734166e-4,-2.73512317095673346e-4,
	    -2.6175015580676858e-4,.00638585891212050914,
	    .00962374215806377941,.00761878061207001043,.00283219055545628054,
	    -.0020984135201272009,-.00573826764216626498,
	    -.0077080424449541462,-.00821011692264844401,
	    -.00765824520346905413,-.00647209729391045177,
	    -.00499132412004966473,-.0034561228971313328,
	    -.00201785580014170775,-7.59430686781961401e-4,
	    2.84173631523859138e-4,.00110891667586337403,
	    .00172901493872728771,.00216812590802684701,.00245357710494539735,
	    .00261281821058334862,.00267141039656276912,.0026520307339598043,
	    .00257411652877287315,.00245389126236094427,.00230460058071795494,
	    .00213684837686712662,.00195896528478870911,.00177737008679454412,
	    .00159690280765839059,.00142111975664438546 };
    static double gama[30] = { .629960524947436582,.251984209978974633,
	    .154790300415655846,.110713062416159013,.0857309395527394825,
	    .0697161316958684292,.0586085671893713576,.0504698873536310685,
	    .0442600580689154809,.0393720661543509966,.0354283195924455368,
	    .0321818857502098231,.0294646240791157679,.0271581677112934479,
	    .0251768272973861779,.0234570755306078891,.0219508390134907203,
	    .020621082823564624,.0194388240897880846,.0183810633800683158,
	    .0174293213231963172,.0165685837786612353,.0157865285987918445,
	    .0150729501494095594,.0144193250839954639,.0138184805735341786,
	    .0132643378994276568,.0127517121970498651,.0122761545318762767,
	    .0118338262398482403 };
    static double ex1 = .333333333333333333;
    static double ex2 = .666666666666666667;
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;


    /* Local variables */
    static int j, k, l, m, l1, l2, is, jr, ks, ju,
	lr, kp1, ias, ibs, lrp1, idum, kmax;
    static double ap[30], p_i[30], p_r[30];
    static double cri[14],crr[14], dri[14],drr[14], upi[14],upr[14];
    static double ac, aw2, pp, wi,wr, t2i,t2r, w2i,w2r, ang, fn13, fn23,
	raw, raw2, sti,str, zai,zar, zbi,zbr, zci,zcr,
	rfn13, atol, btol, tfni,tfnr, azth, rfnu, test,
	zthi,zthr, tzai,tzar, rfnu2, zetai,zetar, ptfni,ptfnr,
	sumai,sumar, sumbi,sumbr, razth, rzthi,rzthr,
	rtzti,rtztr, przthi,przthr;

    /* complex arg,asum,bsum,cfnu,cone,cr,czero,dr,p,phi,przth,ptfn,
    * rfn13,rtzta,rzth,suma,sumb,tfn,t2,up,w,w2,z,za,zb,zc,zeta,zeta1,
    * zeta2,zth */

    rfnu = 1. / *fnu;
/* -----------------------------------------------------------------------
     overflow test (z/fnu too small)
 ----------------------------------------------------------------------- */
    test = DBL_MIN * 1e3;
    ac = *fnu * test;
    if (fabs(*zr) <= ac && fabs(*zi) <= ac) {
	*zeta1r = fabs(log(test)) * 2. + *fnu;
	*zeta1i = 0.;
	*zeta2r = *fnu;
	*zeta2i = 0.;
	*phir = 1.; *phii = 0.;
	*argr = 1.; *argi = 0.;
	return 0;
    }
    // L15:
    zbr = *zr * rfnu;
    zbi = *zi * rfnu;
    rfnu2 = rfnu * rfnu;
/* -----------------------------------------------------------------------
     compute in the fourth quadrant
 ----------------------------------------------------------------------- */
    fn13 = pow(*fnu, ex1);
    fn23 = fn13 * fn13;
    rfn13 = 1. / fn13;
    w2r = coner - zbr * zbr + zbi * zbi;
    w2i = conei - zbr * zbi - zbr * zbi;
    aw2 = zabs(w2r, w2i);
    if (aw2 > .25) {
	goto L130;
    }
/* -----------------------------------------------------------------------
     power series for cabs(w2) <= 0.25d0
 ----------------------------------------------------------------------- */
    k = 1;
    p_r[0] = coner;
    p_i[0] = conei;
    sumar = gama[0];
    sumai = zeroi;
    ap[0] = 1.;
    if (aw2 < *tol) {
	goto L20;
    }
    for (k = 2; k <= 30; ++k) {
	p_r[k - 1] = p_r[k - 2] * w2r - p_i[k - 2] * w2i;
	p_i[k - 1] = p_r[k - 2] * w2i + p_i[k - 2] * w2r;
	sumar += p_r[k - 1] * gama[k - 1];
	sumai += p_i[k - 1] * gama[k - 1];
	ap[k - 1] = ap[k - 2] * aw2;
	if (ap[k - 1] < *tol) {
	    goto L20;
	}
    }
    k = 30;
L20:
    kmax = k;
    zetar = w2r * sumar - w2i * sumai;
    zetai = w2r * sumai + w2i * sumar;
    *argr = zetar * fn23;
    *argi = zetai * fn23;
    zsqrt_sub__(&sumar, &sumai, &zar, &zai);
    zsqrt_sub__(&w2r, &w2i, &str, &sti);
    *zeta2r = str * *fnu;
    *zeta2i = sti * *fnu;
    str = coner + ex2 * (zetar * zar - zetai * zai);
    sti = conei + ex2 * (zetar * zai + zetai * zar);
    *zeta1r = str * *zeta2r - sti * *zeta2i;
    *zeta1i = str * *zeta2i + sti * *zeta2r;
    zar += zar;
    zai += zai;
    zsqrt_sub__(&zar, &zai, &str, &sti);
    *phir = str * rfn13;
    *phii = sti * rfn13;
    if (*ipmtr == 1) {
	goto L120;
    }
/* -----------------------------------------------------------------------
     sum series for asum and bsum
 ----------------------------------------------------------------------- */
    sumbr = zeror;
    sumbi = zeroi;
    for (k = 1; k <= kmax; ++k) {
	sumbr += p_r[k - 1] * beta[k - 1];
	sumbi += p_i[k - 1] * beta[k - 1];
    }
    *asumr = zeror;
    *asumi = zeroi;
    *bsumr = sumbr;
    *bsumi = sumbi;
    l1 = 0;
    l2 = 30;
    btol = *tol * (fabs(*bsumr) + fabs(*bsumi));
    atol = *tol;
    pp = 1.;
    ias = 0;
    ibs = 0;
    if (rfnu2 < *tol) {
	goto L110;
    }
    for (is = 2; is <= 7; ++is) {
	atol /= rfnu2;
	pp *= rfnu2;
	if (ias == 1) {
	    goto L60;
	}
	sumar = zeror;
	sumai = zeroi;
	for (k = 1; k <= kmax; ++k) {
	    m = l1 + k;
	    sumar += p_r[k - 1] * alfa[m - 1];
	    sumai += p_i[k - 1] * alfa[m - 1];
	    if (ap[k - 1] < atol) {
		goto L50;
	    }
	}
L50:
	*asumr += sumar * pp;
	*asumi += sumai * pp;
	if (pp < *tol) {
	    ias = 1;
	}
L60:
	if (ibs == 1) {
	    goto L90;
	}
	sumbr = zeror;
	sumbi = zeroi;
	for (k = 1; k <= kmax; ++k) {
	    m = l2 + k;
	    sumbr += p_r[k - 1] * beta[m - 1];
	    sumbi += p_i[k - 1] * beta[m - 1];
	    if (ap[k - 1] < atol)
		break;
	}
	*bsumr += sumbr * pp;
	*bsumi += sumbi * pp;
	if (pp < btol) {
	    ibs = 1;
	}
L90:
	if (ias == 1 && ibs == 1) {
	    goto L110;
	}
	l1 += 30;
	l2 += 30;
    }
L110:
    *asumr += coner;
    pp = rfnu * rfn13;
    *bsumr *= pp;
    *bsumi *= pp;
L120:
    return 0;
/* -----------------------------------------------------------------------
     cabs(w2) > 0.25d0
 ----------------------------------------------------------------------- */
L130:
    zsqrt_sub__(&w2r, &w2i, &wr, &wi);
    if (wr < 0.) {
	wr = 0.;
    }
    if (wi < 0.) {
	wi = 0.;
    }
    str = coner + wr;
    sti = wi;
    zdiv(str, sti, zbr, zbi, &zar, &zai);
    zlog_sub__(&zar, &zai, &zcr, &zci, &idum);
    if (zci < 0.) {
	zci = 0.;
    }
    if (zci > hpi) {
	zci = hpi;
    }
    if (zcr < 0.) {
	zcr = 0.;
    }
    zthr = (zcr - wr) * 1.5;
    zthi = (zci - wi) * 1.5;
    *zeta1r = zcr * *fnu;
    *zeta1i = zci * *fnu;
    *zeta2r = wr * *fnu;
    *zeta2i = wi * *fnu;
    azth = zabs(zthr, zthi);
    ang = thpi;
    if (zthr >= 0. && zthi < 0.) {
	goto L140;
    }
    ang = hpi;
    if (zthr == 0.) {
	goto L140;
    }
    ang = atan(zthi / zthr);
    if (zthr < 0.) {
	ang += pi;
    }
L140:
    pp = pow(azth, ex2);
    ang *= ex2;
    zetar = pp * cos(ang);
    zetai = pp * sin(ang);
    if (zetai < 0.) {
	zetai = 0.;
    }
    *argr = zetar * fn23;
    *argi = zetai * fn23;
    zdiv(zthr, zthi, zetar, zetai, &rtztr, &rtzti);
    zdiv(rtztr, rtzti, wr, wi, &zar, &zai);
    tzar = zar + zar;
    tzai = zai + zai;
    zsqrt_sub__(&tzar, &tzai, &str, &sti);
    *phir = str * rfn13;
    *phii = sti * rfn13;
    if (*ipmtr == 1) {
	goto L120;
    }
    raw = 1. / sqrt(aw2);
    str = wr * raw;
    sti = -wi * raw;
    tfnr = str * rfnu * raw;
    tfni = sti * rfnu * raw;
    razth = 1. / azth;
    str = zthr * razth;
    sti = -zthi * razth;
    rzthr = str * razth * rfnu;
    rzthi = sti * razth * rfnu;
    zcr = rzthr * ar[1];
    zci = rzthi * ar[1];
    raw2 = 1. / aw2;
    str = w2r * raw2;
    sti = -w2i * raw2;
    t2r = str * raw2;
    t2i = sti * raw2;
    str = t2r * c__[1] + c__[2];
    sti = t2i * c__[1];
    upr[1] = str * tfnr - sti * tfni;
    upi[1] = str * tfni + sti * tfnr;
    *bsumr = upr[1] + zcr;
    *bsumi = upi[1] + zci;
    *asumr = zeror;
    *asumi = zeroi;
    if (rfnu < *tol) {
	goto L220;
    }
    przthr = rzthr;
    przthi = rzthi;
    ptfnr = tfnr;
    ptfni = tfni;
    upr[0] = coner;
    upi[0] = conei;
    pp = 1.;
    btol = *tol * (fabs(*bsumr) + fabs(*bsumi));
    ks = 0;
    kp1 = 2;
    l = 3;
    ias = 0;
    ibs = 0;
    for (lr = 2; lr <= 12; lr += 2) {
	lrp1 = lr + 1;
/* -----------------------------------------------------------------------
     compute two additional cr, dr, and up for two more terms in
     next suma and sumb
 ----------------------------------------------------------------------- */
	for (k = lr; k <= lrp1; ++k) {
	    ++ks;
	    ++kp1;
	    ++l;
	    zar = c__[l - 1];
	    zai = zeroi;
	    for (j = 2; j <= kp1; ++j) {
		++l;
		str = zar * t2r - t2i * zai + c__[l - 1];
		zai = zar * t2i + zai * t2r;
		zar = str;
	    }
	    str = ptfnr * tfnr - ptfni * tfni;
	    ptfni = ptfnr * tfni + ptfni * tfnr;
	    ptfnr = str;
	    upr[kp1 - 1] = ptfnr * zar - ptfni * zai;
	    upi[kp1 - 1] = ptfni * zar + ptfnr * zai;
	    crr[ks - 1] = przthr * br[ks];
	    cri[ks - 1] = przthi * br[ks];
	    str = przthr * rzthr - przthi * rzthi;
	    przthi = przthr * rzthi + przthi * rzthr;
	    przthr = str;
	    drr[ks - 1] = przthr * ar[ks + 1];
	    dri[ks - 1] = przthi * ar[ks + 1];
	}
	pp *= rfnu2;
	if (ias == 1) {
	    goto L180;
	}
	sumar = upr[lrp1 - 1];
	sumai = upi[lrp1 - 1];
	ju = lrp1;
	for (jr = 1; jr <= lr; ++jr) {
	    --ju;
	    sumar += crr[jr - 1] * upr[ju - 1] - cri[jr - 1] * upi[ju - 1];
	    sumai += crr[jr - 1] * upi[ju - 1] + cri[jr - 1] * upr[ju - 1];
	}
	*asumr += sumar;
	*asumi += sumai;
	test = fabs(sumar) + fabs(sumai);
	if (pp < *tol && test < *tol) {
	    ias = 1;
	}
L180:
	if (ibs == 1) {
	    goto L200;
	}
	sumbr = upr[lr + 1] + upr[lrp1 - 1] * zcr - upi[lrp1 - 1] * zci;
	sumbi = upi[lr + 1] + upr[lrp1 - 1] * zci + upi[lrp1 - 1] * zcr;
	ju = lrp1;
	for (jr = 1; jr <= lr; ++jr) {
	    --ju;
	    sumbr = sumbr + drr[jr - 1] * upr[ju - 1] - dri[jr - 1] * upi[ju - 1];
	    sumbi = sumbi + drr[jr - 1] * upi[ju - 1] + dri[jr - 1] * upr[ju - 1];
	}
	*bsumr += sumbr;
	*bsumi += sumbi;
	test = fabs(sumbr) + fabs(sumbi);
	if (pp < btol && test < btol) {
	    ibs = 1;
	}
L200:
	if (ias == 1 && ibs == 1) {
	    goto L220;
	}
    }
L220:
    *asumr += coner;
    str = -(*bsumr) * rfn13;
    sti = -(*bsumi) * rfn13;
    zdiv(str, sti, rtztr, rtzti, bsumr, bsumi);
    goto L120;
} /* zunhj_ */

void zunk1_(double *zr, double *zi, double *fnu,
	    int *kode, int *mr, int *n, double *yr, double *yi,
	    int *nz, double *tol, double *elim, double *alim)
{
/***begin prologue  zunk1
 ***refer to  zbesk

     zunk1 computes K(fnu,z) and its analytic continuation from the
     right half plane to the left half plane by means of the
     uniform asymptotic expansion.
     mr indicates the direction of rotation for analytic continuation.
     nz=-1 means an overflow will occur

 ***routines called  zkscl,zs1s2,zuchk,zunik,d1mach,zabs
 ***end prologue  zunk1
 */

    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;

    /* Local variables */
    static int i, j, k, m, ib, ic, il, kk, nw;
    static int ifn, iuf, inu, init[2];
    static int iflag, kflag, kdflg, ipard, initd;


    /* complex cfn,ck,cone,crsc,cs,cscl,csgn,cspn,csr,css,cwrk,cy,czero,
     *    c1,c2,phi,phid,rz,sum,sumd,s1,s2,y,z,zeta1,zeta1d,zeta2,zeta2d,zr */

    static double fn;
    static double c1i, c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, ang,
	    asc, cki, fnf;
    static double ckr;
    static double cyi[2], fmr, csr, sgn;
    static double bry[3], cyr[2], sti, rzi, zri, str, rzr, zrr, aphi,
	    cscl, phii[2], crsc;
    static double phir[2];
    static double csrr[3], cssr[3], rast, sumi[2], razr;
    static double sumr[2];
    static double ascle;
    static double phidi;
    static double csgni, phidr;
    static double cspni, cwrki[48]	/* was [16][3] */, sumdi;
    static double cspnr, cwrkr[48]	/* was [16][3] */, sumdr;
    static double zeta1i[2], zeta2i[2], zet1di, zet2di,
	zeta1r[2], zeta2r[2], zet1dr, zet2dr;

    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    kdflg = 1;
    *nz = 0;
/* -----------------------------------------------------------------------
     exp(-alim)=exp(-elim)/tol=approx. one precision greater than
     the underflow limit
 ----------------------------------------------------------------------- */
    cscl = 1. / *tol;
    crsc = *tol;
    cssr[0] = cscl;      csrr[0] = crsc;
    cssr[1] = coner;     csrr[1] = coner;
    cssr[2] = crsc;      csrr[2] = cscl;

    bry[0] = DBL_MIN * 1e3 / *tol;
    bry[1] = 1. / bry[0];
    bry[2] = DBL_MAX;
    zrr = *zr;
    zri = *zi;
    if (*zr < 0.) {
	zrr = -(*zr);
	zri = -(*zi);
    }
    j = 2;
    for (i = 1; i <= *n; ++i) {
/*----------------------------------------------------------------------
     j flip flops between 1 and 2 in j = 3 - j
 ----------------------------------------------------------------------- */
	j = 3 - j;
	fn = *fnu + (double) ((float) (i - 1));
	init[j - 1] = 0;
	zunik_(&zrr, &zri, &fn, &c__2, &c__0, tol, &init[j - 1],
	       &phir[j - 1], &phii[j - 1],
	       &zeta1r[j - 1], &zeta1i[j - 1], &zeta2r[j - 1], &zeta2i[j - 1],
	       &sumr[j - 1], &sumi[j - 1],
	       &cwrkr[(j << 4) - 16], &cwrki[(j << 4) - 16]);
	if (*kode == 1) {
	    s1r = zeta1r[j - 1] - zeta2r[j - 1];
	    s1i = zeta1i[j - 1] - zeta2i[j - 1];
	} else {
	    str = zrr + zeta2r[j - 1];
	    sti = zri + zeta2i[j - 1];
	    rast = fn / zabs(str, sti);
	    str = str * rast * rast;
	    sti = -sti * rast * rast;
	    s1r = zeta1r[j - 1] - str;
	    s1i = zeta1i[j - 1] - sti;
	}
	rs1 = s1r;
/* -----------------------------------------------------------------------
     test for underflow and overflow
 ----------------------------------------------------------------------- */
	if (fabs(rs1) > *elim) {
	    goto L60;
	}
	if (kdflg == 1) {
	    kflag = 2;
	}
	if (fabs(rs1) < *alim) {
	    goto L40;
	}
/* -----------------------------------------------------------------------
     refine  test and scale
 ----------------------------------------------------------------------- */
	aphi = zabs(phir[j - 1], phii[j - 1]);
	rs1 += log(aphi);
	if (fabs(rs1) > *elim) {
	    goto L60;
	}
	if (kdflg == 1) {
	    kflag = 1;
	}
	if (rs1 < 0.) {
	    goto L40;
	}
	if (kdflg == 1) {
	    kflag = 3;
	}
L40:
/* -----------------------------------------------------------------------
     scale s1 to keep intermediate arithmetic on scale near
     exponent extremes
 ----------------------------------------------------------------------- */
	s2r = phir[j - 1] * sumr[j - 1] - phii[j - 1] * sumi[j - 1];
	s2i = phir[j - 1] * sumi[j - 1] + phii[j - 1] * sumr[j - 1];
	str = exp(s1r) * cssr[kflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s1r * s2i + s2r * s1i;
	s2r = str;
	if (kflag != 1) {
	    goto L50;
	}
	nw = zuchk(s2r, s2i, *bry, *tol);
	if (nw != 0) {
	    goto L60;
	}
L50:
	cyr[kdflg - 1] = s2r;
	cyi[kdflg - 1] = s2i;
	yr[i] = s2r * csrr[kflag - 1];
	yi[i] = s2i * csrr[kflag - 1];
	if (kdflg == 2) {
	    goto L75;
	}
	kdflg = 2;
	goto L70;
L60:
	if (rs1 > 0.) {
	    goto L300;
	}
/* -----------------------------------------------------------------------
     for zr < 0.0, the i function to be added will overflow
 ----------------------------------------------------------------------- */
	if (*zr < 0.) {
	    goto L300;
	}
	kdflg = 1;
	yr[i] = zeror;
	yi[i] = zeroi;
	++(*nz);
	if (i == 1) {
	    goto L70;
	}
	if (yr[i - 1] == zeror && yi[i - 1] == zeroi) {
	    goto L70;
	}
	yr[i - 1] = zeror;
	yi[i - 1] = zeroi;
	++(*nz);
L70:
	;
    } /* end for(i in 1:*n) */

    i = *n;
L75:
    razr = 1. / zabs(zrr, zri);
    str = zrr * razr;
    sti = -zri * razr;
    rzr = (str + str) * razr;
    rzi = (sti + sti) * razr;
    ckr = fn * rzr;
    cki = fn * rzi;
    ib = i + 1;
    if (*n < ib) {
	goto L160;
    }
/* -----------------------------------------------------------------------
     test last member for underflow and overflow. set sequence to zero
     on underflow.
 ----------------------------------------------------------------------- */
    fn = *fnu + (double) ((float) (*n - 1));
    ipard = 1;
    if (*mr != 0) {
	ipard = 0;
    }
    initd = 0;
    zunik_(&zrr, &zri, &fn, &c__2, &ipard, tol, &initd, &phidr, &phidi,
	   &zet1dr, &zet1di, &zet2dr, &zet2di, &sumdr, &sumdi,
	   &cwrkr[32], &cwrki[32]);
    if (*kode == 1) {
	s1r = zet1dr - zet2dr;
	s1i = zet1di - zet2di;
    } else {
	str = zrr + zet2dr;
	sti = zri + zet2di;
	rast = fn / zabs(str, sti);
	rast *= rast;
	str *= rast;
	sti *= -rast;
	s1r = zet1dr - str;
	s1i = zet1di - sti;
    }
    /* L90; */
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
	goto L95;
    }
    if (fabs(rs1) < *alim) {
	goto L100;
    }
/* ----------------------------------------------------------------------------
     refine estimate and test
 ------------------------------------------------------------------------- */
    aphi = zabs(phidr, phidi);
    rs1 += log(aphi);
    if (fabs(rs1) < *elim) {
	goto L100;
    }
L95:
    if (fabs(rs1) > 0.) {
	goto L300;
    }
/* -----------------------------------------------------------------------
     for zr < 0.0, the i function to be added will overflow
 ----------------------------------------------------------------------- */
    if (*zr < 0.) {
	goto L300;
    }
    *nz = *n;
    for (int i = 1; i <= *n; ++i) {
	yr[i] = zeror;
	yi[i] = zeroi;
    }
    return;
/* ---------------------------------------------------------------------------
     forward recur for remainder of the sequence
 ---------------------------------------------------------------------------- */
L100:
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[kflag - 1];
    ascle = bry[kflag - 1];
    for (int i = ib; i <= *n; ++i) {
	c2r = s2r;
	c2i = s2i;
	s2r = ckr * c2r - cki * c2i + s1r;
	s2i = ckr * c2i + cki * c2r + s1i;
	s1r = c2r;
	s1i = c2i;
	ckr += rzr;
	cki += rzi;
	c2r = s2r * c1r;
	c2i = s2i * c1r;
	yr[i] = c2r;
	yi[i] = c2i;
	if (kflag >= 3) {
	    goto L120;
	}
	str = fabs(c2r);
	sti = fabs(c2i);
	c2m = fmax2(str,sti);
	if (c2m <= ascle) {
	    goto L120;
	}
	++kflag;
	ascle = bry[kflag - 1];
	s1r *= c1r;
	s1i *= c1r;
	s2r = c2r;
	s2i = c2i;
	s1r *= cssr[kflag - 1];
	s1i *= cssr[kflag - 1];
	s2r *= cssr[kflag - 1];
	s2i *= cssr[kflag - 1];
	c1r = csrr[kflag - 1];
L120:
	;
    }
L160:
    if (*mr == 0) {
	return;
    }
/* -----------------------------------------------------------------------
     analytic continuation for re(z) < 0
 ----------------------------------------------------------------------- */
    *nz = 0;
    fmr = (double) (*mr);
    sgn = -fsign(pi, fmr);
/* -----------------------------------------------------------------------
     cspn and csgn are coeff of k and i functions resp.
 ----------------------------------------------------------------------- */
    csgni = sgn;
    inu = (int) ((float) (*fnu));
    fnf = *fnu - (double) ((float) inu);
    ifn = inu + *n - 1;
    ang = fnf * sgn;
    cspnr = cos(ang);
    cspni = sin(ang);
    if (ifn % 2 == 0) {
	goto L170;
    }
    cspnr = -cspnr;
    cspni = -cspni;
L170:
    asc = bry[0];
    iuf = 0;
    kk = *n;
    kdflg = 1;
    --ib;
    ic = ib - 1;
    for (k = 1; k <= *n; ++k) {
	fn = *fnu + (double) ((float) (kk - 1));
/* -----------------------------------------------------------------------
     logic to sort out cases whose parameters were set for the k
     function above
 ----------------------------------------------------------------------- */
	m = 3;
	if (*n > 2) {
	    goto L175;
	}
L172:
	initd = init[j - 1];
	phidr = phir[j - 1];
	phidi = phii[j - 1];
	zet1dr = zeta1r[j - 1];
	zet1di = zeta1i[j - 1];
	zet2dr = zeta2r[j - 1];
	zet2di = zeta2i[j - 1];
	sumdr = sumr[j - 1];
	sumdi = sumi[j - 1];
	m = j;
	j = 3 - j;
	goto L180;
L175:
	if (kk == *n && ib < *n) {
	    goto L180;
	}
	if (kk == ib || kk == ic) {
	    goto L172;
	}
	initd = 0;
L180:
	zunik_(&zrr, &zri, &fn, &c__1, &c__0, tol, &initd, &phidr, &phidi, &
		zet1dr, &zet1di, &zet2dr, &zet2di, &sumdr, &sumdi, &cwrkr[(m

		<< 4) - 16], &cwrki[(m << 4) - 16]);
	if (*kode == 1) {
	    s1r = -zet1dr + zet2dr;
	    s1i = -zet1di + zet2di;
	} else {
	    str = zrr + zet2dr;
	    sti = zri + zet2di;
	    rast = fn / zabs(str, sti);
	    rast *= rast;
	    str *= rast;
	    sti *= -rast;
	    s1r = -zet1dr + str;
	    s1i = -zet1di + sti;
	}
/* L210: */
/* -----------------------------------------------------------------------
     test for underflow and overflow
 ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (fabs(rs1) > *elim) {
	    goto L260;
	}
	if (kdflg == 1) {
	    iflag = 2;
	}
	if (fabs(rs1) < *alim) {
	    goto L220;
	}
/* -----------------------------------------------------------------------
     refine  test and scale
 ----------------------------------------------------------------------- */
	aphi = zabs(phidr, phidi);
	rs1 += log(aphi);
	if (fabs(rs1) > *elim) {
	    goto L260;
	}
	if (kdflg == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L220;
	}
	if (kdflg == 1) {
	    iflag = 3;
	}
L220:
	str = phidr * sumdr - phidi * sumdi;
	sti = phidr * sumdi + phidi * sumdr;
	s2r = -csgni * sti;
	s2i = csgni * str;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L230;
	}
	nw = zuchk(s2r, s2i, *bry, *tol);
	if (nw == 0) {
	    goto L230;
	}
	s2r = zeror;
	s2i = zeroi;
L230:
	cyr[kdflg - 1] = s2r;
	cyi[kdflg - 1] = s2i;
	c2r = s2r;
	c2i = s2i;
	s2r *= csrr[iflag - 1];
	s2i *= csrr[iflag - 1];
/* -----------------------------------------------------------------------
     add i and k functions, k sequence in Y(i), i=1,n
 ----------------------------------------------------------------------- */
	s1r = yr[kk];
	s1i = yi[kk];
	if (*kode == 2) {
	    zs1s2_(&zrr, &zri, &s1r, &s1i, &s2r, &s2i, &nw, &asc, alim, &iuf);
	    *nz += nw;
	}
	yr[kk] = s1r * cspnr - s1i * cspni + s2r;
	yi[kk] = cspnr * s1i + cspni * s1r + s2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	if (c2r != 0. || c2i != 0.) {
	    goto L255;
	}
	kdflg = 1;
	goto L270;
L255:
	if (kdflg == 2) {
	    goto L275;
	}
	kdflg = 2;
	goto L270;
L260:
	if (rs1 > 0.) {
	    goto L300;
	}
	s2r = zeror;
	s2i = zeroi;
	goto L230;
L270:
	;
    }
    k = *n;
L275:
    il = *n - k;
    if (il == 0) {
	return;
    }
/* -----------------------------------------------------------------------
     recur backward for remainder of i sequence and add in the
     k functions, scaling the i sequence during recurrence to keep
     intermediate arithmetic on scale near exponent extremes.
 ----------------------------------------------------------------------- */
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    csr = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    fn = (double) ((float) (inu + il));
    for (i = 1; i <= il; ++i) {
	c2r = s2r;
	c2i = s2i;
	s2r = s1r + (fn + fnf) * (rzr * c2r - rzi * c2i);
	s2i = s1i + (fn + fnf) * (rzr * c2i + rzi * c2r);
	s1r = c2r;
	s1i = c2i;
	fn += -1.;
	c2r = s2r * csr;
	c2i = s2i * csr;
	ckr = c2r;
	cki = c2i;
	c1r = yr[kk];
	c1i = yi[kk];
	if (*kode == 2) {
	    zs1s2_(&zrr, &zri, &c1r, &c1i, &c2r, &c2i, &nw, &asc, alim, &iuf);
	    *nz += nw;
	}
	/* L280: */
	yr[kk] = c1r * cspnr - c1i * cspni + c2r;
	yi[kk] = c1r * cspni + c1i * cspnr + c2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	if (iflag >= 3) {
	    goto L290;
	}
	c2r = fabs(ckr);
	c2i = fabs(cki);
	c2m = fmax2(c2r,c2i);
	if (c2m <= ascle) {
	    goto L290;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= csr;
	s1i *= csr;
	s2r = ckr;
	s2i = cki;
	s1r *= cssr[iflag - 1];
	s1i *= cssr[iflag - 1];
	s2r *= cssr[iflag - 1];
	s2i *= cssr[iflag - 1];
	csr = csrr[iflag - 1];
L290:
	;
    }
    return;
L300:
    *nz = -1;
    return;
} /* zunk1_ */

void zunk2_(double *zr, double *zi, double *fnu,
	    int *kode, int *mr, int *n, double *yr, double * yi,
	    int *nz, double *tol, double *elim, double *alim)
{
/***begin prologue  zunk2
 ***refer to  zbesk

     zunk2 computes K(fnu,z) and its analytic continuation from the
     right half plane to the left half plane by means of the
     uniform asymptotic expansions for H(kind,fnu,zn) and J(fnu,zn)
     where zn is in the right half plane, kind=(3-mr)/2, mr=+1 or
     -1. here zn=zr*i or -zr*i where zr=z if z is in the right
     half plane or zr=-z if z is in the left half plane. mr indic-
     ates the direction of rotation for analytic continuation.
     nz=-1 means an overflow will occur

 ***routines called  zairy,zkscl,zs1s2,zuchk,zunhj,d1mach,zabs
 ***end prologue  zunk2
 */

    /* Initialized data */

    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double cr1r = 1.;
    static double cr1i = 1.73205080756887729;
    static double cr2r = -.5;
    static double cr2i = -.866025403784438647;
    static double cipr[4] = { 1.,0.,-1.,0. };
    static double cipi[4] = { 0.,-1.,0.,1. };


    /* Local variables */
    static int i, j, k, ib, ic, il, kk, in, nw;
    static int nai, ifn, iuf, inu, ndai, idum, iflag, kflag, kdflg, ipard;
    /*
     complex ai,arg,argd,asum,asumd,bsum,bsumd,cfn,ci,cip,ck,cone,crsc,
    * cr1,cr2,cs,cscl,csgn,cspn,csr,css,cy,czero,c1,c2,dai,phi,phid,rz,
    * s1,s2,y,z,zb,zeta1,zeta1d,zeta2,zeta2d,zn,zr
    */
    static double fn;
    static double yy, c1i, c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r,
	    ang, asc, car, cki, fnf;
    static double aii, air;
    static double csi, ckr;
    static double cyi[2], fmr, sar, csr, sgn, zbi;
    static double bry[3], cyr[2], pti, sti, zbr, zni, rzi, ptr, zri, str,
	    znr, rzr, zrr, daii, aarg;
    static double dair, aphi, argi[2], cscl, phii[2], crsc, argr[2];
    static double phir[2], csrr[3], cssr[3], rast, razr;
    static double argdi, ascle;
    static double phidi, argdr;
    static double csgni, phidr, cspni, asumi[2], bsumi[2];
    static double cspnr, asumr[2], bsumr[2];
    static double zeta1i[2], zeta2i[2], zet1di, zet2di, zeta1r[2], zeta2r[2],
	zet1dr, zet2dr, asumdi, bsumdi, asumdr, bsumdr;

    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    kdflg = 1;
    *nz = 0;
/* -----------------------------------------------------------------------
     exp(-alim)=exp(-elim)/tol=approx. one precision greater than
     the underflow limit
 ----------------------------------------------------------------------- */
    cscl = 1. / *tol;
    crsc = *tol;
    cssr[0] = cscl;
    cssr[1] = coner;
    cssr[2] = crsc;
    csrr[0] = crsc;
    csrr[1] = coner;
    csrr[2] = cscl;
    bry[0] = DBL_MIN * 1e3 / *tol;
    bry[1] = 1. / bry[0];
    bry[2] = DBL_MAX;
    zrr = *zr;
    zri = *zi;
    if (*zr < 0.) {
	zrr = -(*zr);
	zri = -(*zi);
    }
    yy = zri;
    znr = zri;
    zni = -zrr;
    zbr = zrr;
    zbi = zri;
    inu = (int) ((float) (*fnu));
    fnf = *fnu - (double) inu;
    ang = -hpi * fnf;
    car = cos(ang);
    sar = sin(ang);
    c2r = hpi * sar;
    c2i = -hpi * car;
    kk = inu % 4 + 1;
    str = c2r * cipr[kk - 1] - c2i * cipi[kk - 1];
    sti = c2r * cipi[kk - 1] + c2i * cipr[kk - 1];
    csr = cr1r * str - cr1i * sti;
    csi = cr1r * sti + cr1i * str;
    if (yy <= 0.) {
	znr = -znr;
	zbi = -zbi;
    }
/* L20 :*/
/* -----------------------------------------------------------------------
     K(fnu,z) is computed from H(2,fnu,-i*z) where z is in the first
     quadrant. fourth quadrant values (yy <= 0.0e0) are computed by
     conjugation since the k function is float on the positive float axis
 ----------------------------------------------------------------------- */
    j = 2;
    for (i = 1; i <= *n; ++i) {
	/*---------------------------------------------------------------
	  j flip flops between 1 and 2 in j = 3 - j
	  --------------------------------------------------------------- */
	j = 3 - j;
	fn = *fnu + (double) ((float) (i - 1));
	zunhj_(&znr, &zni, &fn, &c__0, tol, &phir[j - 1], &phii[j - 1],
	       &argr[j - 1], &argi[j - 1], &zeta1r[j - 1], &zeta1i[j - 1],
	       &zeta2r[j - 1], &zeta2i[j - 1],
	       &asumr[j - 1], &asumi[j - 1], &bsumr[j - 1], &bsumi[j - 1]);
	if (*kode == 1) {
	    s1r = zeta1r[j - 1] - zeta2r[j - 1];
	    s1i = zeta1i[j - 1] - zeta2i[j - 1];
	} else {
	    str = zbr + zeta2r[j - 1];
	    sti = zbi + zeta2i[j - 1];
	    rast = fn / zabs(str, sti);
	    str = str * rast * rast;
	    sti = -sti * rast * rast;
	    s1r = zeta1r[j - 1] - str;
	    s1i = zeta1i[j - 1] - sti;
	}
/* -----------------------------------------------------------------------
     test for underflow and overflow
 ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (fabs(rs1) > *elim) {
	    goto L70;
	}
	if (kdflg == 1) {
	    kflag = 2;
	}
	if (fabs(rs1) < *alim) {
	    goto L50;
	}
/* -----------------------------------------------------------------------
     refine  test and scale
 ----------------------------------------------------------------------- */
	aphi = zabs(phir[j - 1], phii[j - 1]);
	aarg = zabs(argr[j - 1], argi[j - 1]);
	rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
	if (fabs(rs1) > *elim) {
	    goto L70;
	}
	if (kdflg == 1) {
	    kflag = 1;
	}
	if (rs1 < 0.) {
	    goto L50;
	}
	if (kdflg == 1) {
	    kflag = 3;
	}
L50:
/* -----------------------------------------------------------------------
     scale s1 to keep intermediate arithmetic on scale near
     exponent extremes
 ----------------------------------------------------------------------- */
	c2r = argr[j - 1] * cr2r - argi[j - 1] * cr2i;
	c2i = argr[j - 1] * cr2i + argi[j - 1] * cr2r;
	zairy(&c2r, &c2i, &c__0, &c__2,  &air,  &aii,  &nai, &idum);
	zairy(&c2r, &c2i, &c__1, &c__2, &dair, &daii, &ndai, &idum);
	str = dair * bsumr[j - 1] - daii * bsumi[j - 1];
	sti = dair * bsumi[j - 1] + daii * bsumr[j - 1];
	ptr = str * cr2r - sti * cr2i;
	pti = str * cr2i + sti * cr2r;
	str = ptr + (air * asumr[j - 1] - aii * asumi[j - 1]);
	sti = pti + (air * asumi[j - 1] + aii * asumr[j - 1]);
	ptr = str * phir[j - 1] - sti * phii[j - 1];
	pti = str * phii[j - 1] + sti * phir[j - 1];
	s2r = ptr * csr - pti * csi;
	s2i = ptr * csi + pti * csr;
	str = exp(s1r) * cssr[kflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s1r * s2i + s2r * s1i;
	s2r = str;
	if (kflag != 1) {
	    goto L60;
	}
	nw = zuchk(s2r, s2i, *bry, *tol);
	if (nw != 0) {
	    goto L70;
	}
L60:
	if (yy <= 0.) {
	    s2i = -s2i;
	}
	cyr[kdflg - 1] = s2r;
	cyi[kdflg - 1] = s2i;
	yr[i] = s2r * csrr[kflag - 1];
	yi[i] = s2i * csrr[kflag - 1];
	str = csi;
	csi = -csr;
	csr = str;
	if (kdflg == 2) {
	    goto L85;
	}
	kdflg = 2;
	goto L80;
L70:
	if (rs1 > 0.) {
	    goto L320;
	}
/* -----------------------------------------------------------------------
     for zr < 0.0, the i function to be added will overflow
 ----------------------------------------------------------------------- */
	if (*zr < 0.) {
	    goto L320;
	}
	kdflg = 1;
	yr[i] = zeror;
	yi[i] = zeroi;
	++(*nz);
	str = csi;
	csi = -csr;
	csr = str;
	if (i > 1 && (yr[i - 1] != zeror || yi[i - 1] != zeroi)) {
	    yr[i - 1] = zeror;
	    yi[i - 1] = zeroi;
	    ++(*nz);
	}
L80:
	;
    }
    i = *n;
L85:
    razr = 1. / zabs(zrr, zri);
    str = zrr * razr;
    sti = -zri * razr;
    rzr = (str + str) * razr;
    rzi = (sti + sti) * razr;
    ckr = fn * rzr;
    cki = fn * rzi;
    ib = i + 1;
    if (*n < ib) {
	goto L180;
    }
/* -----------------------------------------------------------------------
     test last member for underflow and overflow. set sequence to zero
     on underflow.
 ----------------------------------------------------------------------- */
    fn = *fnu + (double) ((float) (*n - 1));
    ipard = 1;
    if (*mr != 0) {
	ipard = 0;
    }
    zunhj_(&znr, &zni, &fn, &ipard, tol, &phidr, &phidi, &argdr, &argdi,
	   &zet1dr, &zet1di, &zet2dr, &zet2di,
	   &asumdr, &asumdi, &bsumdr, &bsumdi);
    if (*kode == 1) {
	s1r = zet1dr - zet2dr;
	s1i = zet1di - zet2di;
    } else { /* kode == 2 */
	str = zbr + zet2dr;
	sti = zbi + zet2di;
	rast = fn / zabs(str, sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = zet1dr - str;
	s1i = zet1di - sti;
    }
    /* L100: */
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
	goto L105;
    }
    if (fabs(rs1) < *alim) {
	goto L120;
    }
/* ----------------------------------------------------------------------------
     refine estimate and test
 ------------------------------------------------------------------------- */
    aphi = zabs(phidr, phidi);
    rs1 += log(aphi);
    if (fabs(rs1) < *elim) {
	goto L120;
    }
L105:
    if (rs1 > 0.) {
	goto L320;
    }
/* -----------------------------------------------------------------------
     for zr < 0.0, the i function to be added will overflow
 ----------------------------------------------------------------------- */
    if (*zr < 0.) {
	goto L320;
    }
    *nz = *n;
    for (i = 1; i <= *n; ++i) {
	yr[i] = zeror;
	yi[i] = zeroi;
    }
    return;
L120:
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[kflag - 1];
    ascle = bry[kflag - 1];
    for (i = ib; i <= *n; ++i) {
	c2r = s2r;
	c2i = s2i;
	s2r = ckr * c2r - cki * c2i + s1r;
	s2i = ckr * c2i + cki * c2r + s1i;
	s1r = c2r;
	s1i = c2i;
	ckr += rzr;
	cki += rzi;
	c2r = s2r * c1r;
	c2i = s2i * c1r;
	yr[i] = c2r;
	yi[i] = c2i;
	if (kflag >= 3) {
	    goto L130;
	}
	str = fabs(c2r);
	sti = fabs(c2i);
	c2m = fmax2(str,sti);
	if (c2m <= ascle) {
	    goto L130;
	}
	++kflag;
	ascle = bry[kflag - 1];
	s1r *= c1r;
	s1i *= c1r;
	s2r = c2r;
	s2i = c2i;
	s1r *= cssr[kflag - 1];
	s1i *= cssr[kflag - 1];
	s2r *= cssr[kflag - 1];
	s2i *= cssr[kflag - 1];
	c1r = csrr[kflag - 1];
L130:
	;
    }
L180:
    if (*mr == 0) {
	return;
    }
/* -----------------------------------------------------------------------
     analytic continuation for re(z) < 0
 ----------------------------------------------------------------------- */
    *nz = 0;
    fmr = (double) (*mr);
    sgn = -fsign(pi, fmr);
/* -----------------------------------------------------------------------
     cspn and csgn are coeff of k and i funcions resp.
 ----------------------------------------------------------------------- */
    csgni = sgn;
    if (yy <= 0.) {
	csgni = -csgni;
    }
    ifn = inu + *n - 1;
    ang = fnf * sgn;
    cspnr = cos(ang);
    cspni = sin(ang);
    if (ifn % 2 == 0) {
	goto L190;
    }
    cspnr = -cspnr;
    cspni = -cspni;
L190:
/* -----------------------------------------------------------------------
     cs=coeff of the j function to get the i function. I(fnu,z) is
     computed from exp(i*fnu*hpi)*J(fnu,-i*z) where z is in the first
     quadrant. fourth quadrant values (yy <= 0.0e0) are computed by
     conjugation since the i function is float on the positive float axis
 ----------------------------------------------------------------------- */
    csr = sar * csgni;
    csi = car * csgni;
    in = ifn % 4 + 1;
    c2r = cipr[in - 1];
    c2i = cipi[in - 1];
    str = csr * c2r + csi * c2i;
    csi = -csr * c2i + csi * c2r;
    csr = str;
    asc = bry[0];
    iuf = 0;
    kk = *n;
    kdflg = 1;
    --ib;
    ic = ib - 1;
    for (k = 1; k <= *n; ++k) {
	fn = *fnu + (double) ((float) (kk - 1));
/* -----------------------------------------------------------------------
     logic to sort out cases whose parameters were set for the k
     function above
 ----------------------------------------------------------------------- */
	if (*n > 2) {
	    goto L175;
	}
/* Loop : */
L172:
	phidr = phir[j - 1];
	phidi = phii[j - 1];
	argdr = argr[j - 1];
	argdi = argi[j - 1];
	zet1dr = zeta1r[j - 1];
	zet1di = zeta1i[j - 1];
	zet2dr = zeta2r[j - 1];
	zet2di = zeta2i[j - 1];
	asumdr = asumr[j - 1];
	asumdi = asumi[j - 1];
	bsumdr = bsumr[j - 1];
	bsumdi = bsumi[j - 1];
	j = 3 - j;
	goto L210;
L175:
	if (kk == *n && ib < *n) {
	    goto L210;
	}
	if (kk == ib || kk == ic) {
	    goto L172;
	}
	zunhj_(&znr, &zni, &fn, &c__0, tol, &phidr, &phidi, &argdr, &argdi,
	       &zet1dr, &zet1di, &zet2dr, &zet2di,
	       &asumdr, &asumdi, &bsumdr, &bsumdi);
L210:
	if (*kode == 1) {
	    s1r = -zet1dr + zet2dr;
	    s1i = -zet1di + zet2di;
	} else {
	    str = zbr + zet2dr;
	    sti = zbi + zet2di;
	    rast = fn / zabs(str, sti);
	    str = str * rast * rast;
	    sti = -sti * rast * rast;
	    s1r = -zet1dr + str;
	    s1i = -zet1di + sti;
	}
	/* L230: */
/* -----------------------------------------------------------------------
     test for underflow and overflow
 ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (fabs(rs1) > *elim) {
	    goto L280;
	}
	if (kdflg == 1) {
	    iflag = 2;
	}
	if (fabs(rs1) < *alim) {
	    goto L240;
	}
/* -----------------------------------------------------------------------
     refine  test and scale
 ----------------------------------------------------------------------- */
	aphi = zabs(phidr, phidi);
	aarg = zabs(argdr, argdi);
	rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
	if (fabs(rs1) > *elim) {
	    goto L280;
	}
	if (kdflg == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L240;
	}
	if (kdflg == 1) {
	    iflag = 3;
	}
L240:
	zairy(&argdr, &argdi, &c__0, &c__2,  &air,  &aii,  &nai, &idum);
	zairy(&argdr, &argdi, &c__1, &c__2, &dair, &daii, &ndai, &idum);
	str = dair * bsumdr - daii * bsumdi;
	sti = dair * bsumdi + daii * bsumdr;
	str += air * asumdr - aii * asumdi;
	sti += air * asumdi + aii * asumdr;
	ptr = str * phidr - sti * phidi;
	pti = str * phidi + sti * phidr;
	s2r = ptr * csr - pti * csi;
	s2i = ptr * csi + pti * csr;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L250;
	}
	nw = zuchk(s2r, s2i, *bry, *tol);
	if (nw == 0) {
	    goto L250;
	}
	s2r = zeror;
	s2i = zeroi;
L250:
	if (yy <= 0.) {
	    s2i = -s2i;
	}
	cyr[kdflg - 1] = s2r;
	cyi[kdflg - 1] = s2i;
	c2r = s2r;
	c2i = s2i;
	s2r *= csrr[iflag - 1];
	s2i *= csrr[iflag - 1];
/* -----------------------------------------------------------------------
     add i and k functions, k sequence in Y(i), i=1,n
 ----------------------------------------------------------------------- */
	s1r = yr[kk];
	s1i = yi[kk];
	if (*kode == 2) {
	    zs1s2_(&zrr, &zri, &s1r, &s1i, &s2r, &s2i, &nw, &asc, alim, &iuf);
	    *nz += nw;
	}
        /* L270: */
	yr[kk] = s1r * cspnr - s1i * cspni + s2r;
	yi[kk] = s1r * cspni + s1i * cspnr + s2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	str = csi;
	csi = -csr;
	csr = str;
	if (c2r != 0. || c2i != 0.) {
	    goto L255;
	}
	kdflg = 1;
	goto L290;
L255:
	if (kdflg == 2) {
	    goto L295;
	}
	kdflg = 2;
	goto L290;
L280:
	if (rs1 > 0.) {
	    goto L320;
	}
	s2r = zeror;
	s2i = zeroi;
	goto L250;
L290:
	;
    }
    k = *n;
L295:
    il = *n - k;
    if (il == 0) {
	return;
    }
/* -----------------------------------------------------------------------
     recur backward for remainder of i sequence and add in the
     k functions, scaling the i sequence during recurrence to keep
     intermediate arithmetic on scale near exponent extremes.
 ----------------------------------------------------------------------- */
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    csr = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    fn = (double) ((float) (inu + il));
    for (i = 1; i <= il; ++i) {
	c2r = s2r;
	c2i = s2i;
	s2r = s1r + (fn + fnf) * (rzr * c2r - rzi * c2i);
	s2i = s1i + (fn + fnf) * (rzr * c2i + rzi * c2r);
	s1r = c2r;
	s1i = c2i;
	fn += -1.;
	c2r = s2r * csr;
	c2i = s2i * csr;
	ckr = c2r;
	cki = c2i;
	c1r = yr[kk];
	c1i = yi[kk];
	if (*kode == 2) {
	    zs1s2_(&zrr, &zri, &c1r, &c1i, &c2r, &c2i, &nw, &asc, alim, &iuf);
	    *nz += nw;
	}
	/* L300: */
	yr[kk] = c1r * cspnr - c1i * cspni + c2r;
	yi[kk] = c1r * cspni + c1i * cspnr + c2i;
	--kk;
	cspnr = -cspnr;
	cspni = -cspni;
	if (iflag >= 3) {
	    goto L310;
	}
	c2r = fabs(ckr);
	c2i = fabs(cki);
	c2m = fmax2(c2r,c2i);
	if (c2m <= ascle) {
	    goto L310;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= csr;
	s1i *= csr;
	s2r = ckr;
	s2i = cki;
	s1r *= cssr[iflag - 1];
	s1i *= cssr[iflag - 1];
	s2r *= cssr[iflag - 1];
	s2i *= cssr[iflag - 1];
	csr = csrr[iflag - 1];
L310:
	;
    }
    return;
L320:
    *nz = -1;
    return;
} /* zunk2_ */

int zbuni(double *zr, double *zi, double *fnu,
	  int kode, int n, int nui,
	  double *yr, double *yi, int *nlast,
	  double fnul, double tol, double elim, double alim)
{
/***begin prologue  zbuni
 ***refer to  zbesi,zbesk

     zbuni computes the i bessel function for large cabs(z) >
     fnul and fnu+n-1 < fnul. the order is increased from
     fnu+n-1 greater than fnul by adding nui and computing
     according to the uniform asymptotic expansion for I(fnu,z)
     on iform=1 and the expansion for J(fnu,z) on iform=2

 ***routines called  zuni1,zuni2,zabs,d1mach
 ***end prologue  zbuni
 */
    int nz = 0; // return value

    /* Local variables */
    int i, k,  nl, nw, iflag, iform;
    double ax, ay, ascle, csclr, cscrr,
	c1i, c1m, c1r, s1i, s2i, s1r, s2r, cyi[2], gnu, raz,
	cyr[2], sti, bry[3], rzi, str, rzr, dfnu, fnui;

    /* complex cscl,cscr,cy,rz,st,s1,s2,y,z
     Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
    ax = fabs(*zr) * 1.7321;
    ay = fabs(*zi);
    iform = 1;
    if (ay > ax) {
	iform = 2;
    }
    if (nui == 0) {
	goto L60;
    }
    fnui = (double) ((float) (nui));
    dfnu = *fnu + (double) ((float) (n - 1));
    gnu = dfnu + fnui;
    if (iform == 2) { /* L10: */
	/* -----------------------------------------------------------------------
	   asymptotic expansion for J(fnu,z*exp(m*hpi)) for large fnu
	   applied in pi/3 < abs(arg(z)) <= pi/2 where m=+i or -i
	   and hpi=pi/2
	   ----------------------------------------------------------------------- */
	zuni2_(zr, zi, &gnu, &kode, &c__2, cyr, cyi, &nw, nlast,
	       &fnul, &tol, &elim, &alim);
    } else { /* iform == 1 */
	/* -----------------------------------------------------------------------
	   asymptotic expansion for I(fnu,z) for large fnu applied in
	   -pi/3 <= arg(z) <= pi/3
	   ----------------------------------------------------------------------- */
	zuni1_(zr, zi, &gnu, &kode, &c__2, cyr, cyi, &nw, nlast,
	       &fnul, &tol, &elim, &alim);
    }

    if (nw < 0) {
	goto L50;
    }
    if (nw != 0) {
	goto L90;
    }
    str = zabs(*cyr, *cyi);
/* ----------------------------------------------------------------------
     scale backward recurrence, bry(3) is defined but never used
 ---------------------------------------------------------------------- */
    bry[0] = DBL_MIN * 1e3 / tol;
    bry[1] = 1. / bry[0];
    bry[2] = bry[1];
    iflag = 2;
    ascle = bry[1];
    csclr = 1.;
    if (str > bry[0]) {
	goto L21;
    }
    iflag = 1;
    ascle = bry[0];
    csclr = 1. / tol;
    goto L25;
L21:
    if (str < bry[1]) {
	goto L25;
    }
    iflag = 3;
    ascle = bry[2];
    csclr = tol;
L25:
    cscrr = 1. / csclr;
    s1r = cyr[1] * csclr;
    s1i = cyi[1] * csclr;
    s2r = cyr[0] * csclr;
    s2i = cyi[0] * csclr;
    raz = 1. / zabs(*zr, *zi);
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    for (i = 1; i <= nui; ++i) {
	str = s2r;
	sti = s2i;
	s2r = (dfnu + fnui) * (rzr * str - rzi * sti) + s1r;
	s2i = (dfnu + fnui) * (rzr * sti + rzi * str) + s1i;
	s1r = str;
	s1i = sti;
	fnui += -1.;
	if (iflag >= 3) {
	    goto L30;
	}
	str = s2r * cscrr;
	sti = s2i * cscrr;
	c1r = fabs(str);
	c1i = fabs(sti);
	c1m = fmax2(c1r,c1i);
	if (c1m <= ascle) {
	    goto L30;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= cscrr;
	s1i *= cscrr;
	s2r = str;
	s2i = sti;
	csclr *= tol;
	cscrr = 1. / csclr;
	s1r *= csclr;
	s1i *= csclr;
	s2r *= csclr;
	s2i *= csclr;
L30:
	;
    }
    yr[n] = s2r * cscrr;
    yi[n] = s2i * cscrr;
    if (n == 1) {
	return nz;
    }
    nl = n - 1;
    fnui = (double) ((float) nl);
    k = nl;
    for (i = 1; i <= nl; ++i) {
	str = s2r;
	sti = s2i;
	s2r = (*fnu + fnui) * (rzr * str - rzi * sti) + s1r;
	s2i = (*fnu + fnui) * (rzr * sti + rzi * str) + s1i;
	s1r = str;
	s1i = sti;
	str = s2r * cscrr;
	sti = s2i * cscrr;
	yr[k] = str;
	yi[k] = sti;
	fnui += -1.;
	--k;
	if (iflag >= 3) {
	    continue;
	}
	c1r = fabs(str);
	c1i = fabs(sti);
	c1m = fmax2(c1r,c1i);
	if (c1m <= ascle) {
	    continue;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= cscrr;
	s1i *= cscrr;
	s2r = str;
	s2i = sti;
	csclr *= tol;
	cscrr = 1. / csclr;
	s1r *= csclr;
	s1i *= csclr;
	s2r *= csclr;
	s2i *= csclr;
	// L40:
    }
    return nz;

L50: return (nw == -2) ? -2 : -1;


L60:
    if (iform == 1) {
	/* -----------------------------------------------------------------------
	   asymptotic expansion for I(fnu,z) for large fnu applied in
	   -pi/3 <= arg(z) <= pi/3
	   ----------------------------------------------------------------------- */
	zuni1_(zr, zi, fnu, &kode, &n, &yr[1], &yi[1], &nw, nlast,
	       &fnul, &tol, &elim, &alim);
    } else { /* L70: */
	/* -----------------------------------------------------------------------
	   asymptotic expansion for J(fnu,z*exp(m*hpi)) for large fnu
	   applied in pi/3 < abs(arg(z)) <= pi/2 where m=+i or -i
	   and hpi=pi/2
	   ----------------------------------------------------------------------- */
	zuni2_(zr, zi, fnu, &kode, &n, &yr[1], &yi[1], &nw, nlast,
	       &fnul, &tol, &elim, &alim);
    }
    if (nw < 0) {
	goto L50;
    }
    return nw;

L90:
    *nlast = n;
    return nz;
} /* zbuni_ */

void zuni1_(double *zr, double *zi, double *fnu, int *kode, int *n,
	    double *yr, double *yi, int *nz, int *nlast,
	    double *fnul, double *tol, double *elim, double *alim)
{
/*
***begin prologue  zuni1
***refer to  zbesi,zbesk

     zuni1 computes I(fnu,z)  by means of the uniform asymptotic
     expansion for I(fnu,z) in -pi/3 <= arg z <= pi/3.

     fnul is the smallest order permitted for the asymptotic
     expansion. nlast=0 means all of the y values were set.
     nlast != 0 is the number left to be computed by another
     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1 < fnul.
     Y(i)=czero for i=nlast+1,n

 ***routines called  zuchk,zunik,zuoik,d1mach,zabs
***end prologue  zuni1 */

    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;

    /* Local variables */
    static int i, k, m, nd, nn, nw, nuf, init, iflag;
/* complex cfn,cone,crsc,cscl,csr,css,cwrk,czero,c1,c2,phi,rz,sum,s1,
 *  s2,y,z,zeta1,zeta2 */
    static double fn, c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r,
	bry[3], cyr[2], cyi[2], rzi, rzr, sti, str, aphi,
	cscl, crsc, phii, phir, csrr[3], cssr[3], rast, sumi, sumr,
	ascle, cwrki[16], cwrkr[16], zeta1i, zeta2i, zeta1r, zeta2r;

/*     Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    nd = *n;
    *nlast = 0;
/* -----------------------------------------------------------------------
     computed values with exponents between alim and elim in mag-
     nitude are scaled to keep intermediate arithmetic on scale,
     exp(alim)=exp(elim)*tol
 ----------------------------------------------------------------------- */
    cscl = 1. / *tol;
    crsc = *tol;
    cssr[0] = cscl;
    cssr[1] = coner;
    cssr[2] = crsc;
    csrr[0] = crsc;
    csrr[1] = coner;
    csrr[2] = cscl;
    bry[0] = DBL_MIN * 1e3 / *tol;
/* -----------------------------------------------------------------------
     check for underflow and overflow on first member
 ----------------------------------------------------------------------- */
    fn = fmax2(*fnu,1.);
    init = 0;
    zunik_(zr, zi, &fn, &c__1, &c__1, tol, &init, &phir, &phii, &zeta1r, &
	    zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
    if (*kode == 1) {
	goto L10;
    }
    str = *zr + zeta2r;
    sti = *zi + zeta2i;
    rast = fn / zabs(str, sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zeta1r + str;
    s1i = -zeta1i + sti;
    goto L20;
L10:
    s1r = -zeta1r + zeta2r;
    s1i = -zeta1i + zeta2i;
L20:
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
	goto L130;
    }
L30:
    nn = imin2(2,nd);
    for (i = 1; i <= nn; ++i) {
	fn = *fnu + (double) ((float) (nd - i));
	init = 0;
	zunik_(zr, zi, &fn, &c__1, &c__0, tol, &init, &phir, &phii, &zeta1r, &
		zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr, cwrki);
	if (*kode == 1) {
	    goto L40;
	}
	str = *zr + zeta2r;
	sti = *zi + zeta2i;
	rast = fn / zabs(str, sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = -zeta1r + str;
	s1i = -zeta1i + sti + *zi;
	goto L50;
L40:
	s1r = -zeta1r + zeta2r;
	s1i = -zeta1i + zeta2i;
L50:
/* -----------------------------------------------------------------------
     test for underflow and overflow
 ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (fabs(rs1) > *elim) {
	    goto L110;
	}
	if (i == 1) {
	    iflag = 2;
	}
	if (fabs(rs1) < *alim) {
	    goto L60;
	}
/* -----------------------------------------------------------------------
     refine  test and scale
 ----------------------------------------------------------------------- */
	aphi = zabs(phir, phii);
	rs1 += log(aphi);
	if (fabs(rs1) > *elim) {
	    goto L110;
	}
	if (i == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L60;
	}
	if (i == 1) {
	    iflag = 3;
	}
L60:
/* -----------------------------------------------------------------------
     scale s1 if cabs(s1) < ascle
 ----------------------------------------------------------------------- */
	s2r = phir * sumr - phii * sumi;
	s2i = phir * sumi + phii * sumr;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L70;
	}
	nw = zuchk(s2r, s2i, *bry, *tol);
	if (nw != 0) {
	    goto L110;
	}
L70:
	cyr[i - 1] = s2r;
	cyi[i - 1] = s2i;
	m = nd - i + 1;
	yr[m] = s2r * csrr[iflag - 1];
	yi[m] = s2i * csrr[iflag - 1];
    }
    if (nd <= 2) {
	goto L100;
    }
    rast = 1. / zabs(*zr, *zi);
    str = *zr * rast;
    sti = -(*zi) * rast;
    rzr = (str + str) * rast;
    rzi = (sti + sti) * rast;
    bry[1] = 1. / bry[0];
    bry[2] = DBL_MAX;
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    k = nd - 2;
    fn = (double) ((float) k);
    for (i = 3; i <= nd; ++i) {
	c2r = s2r;
	c2i = s2i;
	s2r = s1r + (*fnu + fn) * (rzr * c2r - rzi * c2i);
	s2i = s1i + (*fnu + fn) * (rzr * c2i + rzi * c2r);
	s1r = c2r;
	s1i = c2i;
	c2r = s2r * c1r;
	c2i = s2i * c1r;
	yr[k] = c2r;
	yi[k] = c2i;
	--k;
	fn += -1.;
	if (iflag >= 3) {
	    goto L90;
	}
	str = fabs(c2r);
	sti = fabs(c2i);
	c2m = fmax2(str,sti);
	if (c2m <= ascle) {
	    goto L90;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= c1r;
	s1i *= c1r;
	s2r = c2r;
	s2i = c2i;
	s1r *= cssr[iflag - 1];
	s1i *= cssr[iflag - 1];
	s2r *= cssr[iflag - 1];
	s2i *= cssr[iflag - 1];
	c1r = csrr[iflag - 1];
L90:
	;
    }
L100:
    return;
/* -----------------------------------------------------------------------
     set underflow and update parameters
 ----------------------------------------------------------------------- */
L110:
    if (rs1 > 0.) {
	goto L120;
    }
    yr[nd] = zeror;
    yi[nd] = zeroi;
    ++(*nz);
    --nd;
    if (nd == 0) {
	goto L100;
    }
    nuf = zuoik(zr, zi, fnu, *kode, 1, nd, &yr[1], &yi[1], *tol, *elim, *alim);
    if (nuf < 0) {
	goto L120;
    }
    nd -= nuf;
    *nz += nuf;
    if (nd == 0) {
	goto L100;
    }
    fn = *fnu + (double) ((float) (nd - 1));
    if (fn >= *fnul) {
	goto L30;
    }
    *nlast = nd;
    return;
L120:
    *nz = -1;
    return;
L130:
    if (rs1 > 0.) {
	goto L120;
    }
    *nz = *n;
    for (i = 1; i <= *n; ++i) {
	yr[i] = zeror;
	yi[i] = zeroi;
    }
    return;
} /* zuni1_ */

void zuni2_(double *zr, double *zi, double *fnu,
	    int *kode, int *n, double *yr, double *yi, int *nz,
	    int *nlast, double *fnul, double *tol, double *elim, double *alim)
{
/***begin prologue  zuni2
 ***refer to  zbesi,zbesk

     zuni2 computes I(fnu,z) in the right half plane by means of
     uniform asymptotic expansion for J(fnu,zn) where zn is z*i
     or -z*i and zn is in the right half plane also.

     fnul is the smallest order permitted for the asymptotic
     expansion. nlast=0 means all of the y values were set.
     nlast != 0 is the number left to be computed by another
     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1 < fnul.
     Y(i)=czero for i=nlast+1,n

 ***routines called  zairy,zuchk,zunhj,zuoik,d1mach,zabs
 ***end prologue  zuni2
 */

    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double cipr[4] = { 1.,0.,-1.,0. };
    static double cipi[4] = { 0.,1.,0.,-1. };

    /* Local variables */
    static int j, k, nd, in, nn, nw, nai, nuf, inu, ndai, idum, iflag;

    /* complex ai,arg,asum,bsum,cfn,ci,cid,cip,cone,crsc,cscl,csr,css,
       czero,c1,c2,dai,phi,rz,s1,s2,y,z,zb,zeta1,zeta2,zn : */
    static double aii, air, fn;
    static double c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, ang, car;
    static double cyi[2], cyr[2], sar;
    static double bry[3], raz, sti, str, zbi, zbr, zni, znr, rzi, rzr,
	    cidi, aarg;
    static double daii, dair, aphi, argi, argr, cscl, phii, crsc;
    static double phir, csrr[3], cssr[3], rast;
    static double ascle, asumi, asumr, bsumi, bsumr;
    static double zeta1i, zeta2i, zeta1r, zeta2r;

     /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */

    *nz = 0;
    nd = *n;
    *nlast = 0;
/* -----------------------------------------------------------------------
     computed values with exponents between alim and elim in mag-
     nitude are scaled to keep intermediate arithmetic on scale,
     exp(alim)=exp(elim)*tol
 ----------------------------------------------------------------------- */
    cscl = 1. / *tol;
    crsc = *tol;
    cssr[0] = cscl;
    cssr[1] = coner;
    cssr[2] = crsc;
    csrr[0] = crsc;
    csrr[1] = coner;
    csrr[2] = cscl;
    bry[0] = DBL_MIN * 1e3 / *tol;
/* -----------------------------------------------------------------------
     zn is in the right half plane after rotation by ci or -ci
 ----------------------------------------------------------------------- */
    znr = *zi;
    zni = -(*zr);
    zbr = *zr;
    zbi = *zi;
    cidi = -coner;
    inu = (int) ((float) (*fnu));
    ang = hpi * (*fnu - (double) ((float) inu));
    c2r = cos(ang);
    c2i = sin(ang);
    car = c2r;
    sar = c2i;
    in = inu + *n - 1;
    in = in % 4 + 1;
    str = c2r * cipr[in - 1] - c2i * cipi[in - 1];
    c2i = c2r * cipi[in - 1] + c2i * cipr[in - 1];
    c2r = str;
    if (*zi > 0.) {
	goto L10;
    }
    znr = -znr;
    zbi = -zbi;
    cidi = -cidi;
    c2i = -c2i;
L10:
/* -----------------------------------------------------------------------
     check for underflow and overflow on first member
 ----------------------------------------------------------------------- */
    fn = fmax2(*fnu,1.);
    zunhj_(&znr, &zni, &fn, &c__1, tol, &phir, &phii, &argr, &argi, &zeta1r, &
	    zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr, &bsumi);
    if (*kode == 1) {
	goto L20;
    }
    str = zbr + zeta2r;
    sti = zbi + zeta2i;
    rast = fn / zabs(str, sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zeta1r + str;
    s1i = -zeta1i + sti;
    goto L30;
L20:
    s1r = -zeta1r + zeta2r;
    s1i = -zeta1i + zeta2i;
L30:
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
	goto L150;
    }
L40:
    nn = imin2(2,nd);
    for (int i = 1; i <= nn; ++i) {
	fn = *fnu + (double) ((float) (nd - i));
	zunhj_(&znr, &zni, &fn, &c__0, tol, &phir, &phii, &argr, &argi, &
		zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr, &
		bsumi);
	if (*kode == 1) {
	    goto L50;
	}
	str = zbr + zeta2r;
	sti = zbi + zeta2i;
	rast = fn / zabs(str, sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = -zeta1r + str;
	s1i = -zeta1i + sti + fabs(*zi);
	goto L60;
L50:
	s1r = -zeta1r + zeta2r;
	s1i = -zeta1i + zeta2i;
L60:
/* -----------------------------------------------------------------------
     test for underflow and overflow
 ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (fabs(rs1) > *elim) {
	    goto L120;
	}
	if (i == 1) {
	    iflag = 2;
	}
	if (fabs(rs1) < *alim) {
	    goto L70;
	}
/* -----------------------------------------------------------------------
     refine  test and scale
 -----------------------------------------------------------------------
 ----------------------------------------------------------------------- */
	aphi = zabs(phir, phii);
	aarg = zabs(argr, argi);
	rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
	if (fabs(rs1) > *elim) {
	    goto L120;
	}
	if (i == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L70;
	}
	if (i == 1) {
	    iflag = 3;
	}
L70:
/* -----------------------------------------------------------------------
     scale s1 to keep intermediate arithmetic on scale near
     exponent extremes
 ----------------------------------------------------------------------- */
	zairy(&argr, &argi, &c__0, &c__2,  &air,  &aii,  &nai, &idum);
	zairy(&argr, &argi, &c__1, &c__2, &dair, &daii, &ndai, &idum);
	str = dair * bsumr - daii * bsumi;
	sti = dair * bsumi + daii * bsumr;
	str += air * asumr - aii * asumi;
	sti += air * asumi + aii * asumr;
	s2r = phir * str - phii * sti;
	s2i = phir * sti + phii * str;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag == 1) {
	    nw = zuchk(s2r, s2i, *bry, *tol);
	    if (nw != 0) {
		goto L120;
	    }
	}
	/* L80: */
	if (*zi <= 0.) {
	    s2i = -s2i;
	}
	str = s2r * c2r - s2i * c2i;
	s2i = s2r * c2i + s2i * c2r;
	s2r = str;
	cyr[i - 1] = s2r;
	cyi[i - 1] = s2i;
	j = nd - i + 1;
	yr[j] = s2r * csrr[iflag - 1];
	yi[j] = s2i * csrr[iflag - 1];
	str = -c2i * cidi;
	c2i = c2r * cidi;
	c2r = str;
    }
    if (nd <= 2) {
	goto L110;
    }
    raz = 1. / zabs(*zr, *zi);
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    bry[1] = 1. / bry[0];
    bry[2] = DBL_MAX;
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    k = nd - 2;
    fn = (double) ((float) k);
    for (int i = 3; i <= nd; ++i) {
	c2r = s2r;
	c2i = s2i;
	s2r = s1r + (*fnu + fn) * (rzr * c2r - rzi * c2i);
	s2i = s1i + (*fnu + fn) * (rzr * c2i + rzi * c2r);
	s1r = c2r;
	s1i = c2i;
	c2r = s2r * c1r;
	c2i = s2i * c1r;
	yr[k] = c2r;
	yi[k] = c2i;
	--k;
	fn += -1.;
	if (iflag >= 3) {
	    goto L100;
	}
	str = fabs(c2r);
	sti = fabs(c2i);
	c2m = fmax2(str,sti);
	if (c2m <= ascle) {
	    goto L100;
	}
	++iflag;
	ascle = bry[iflag - 1];
	s1r *= c1r;
	s1i *= c1r;
	s2r = c2r;
	s2i = c2i;
	s1r *= cssr[iflag - 1];
	s1i *= cssr[iflag - 1];
	s2r *= cssr[iflag - 1];
	s2i *= cssr[iflag - 1];
	c1r = csrr[iflag - 1];
L100:
	;
    }
L110:
    return;
L120:
    if (rs1 > 0.) {
	goto L140;
    }
/* -----------------------------------------------------------------------
     set underflow and update parameters
 ----------------------------------------------------------------------- */
    yr[nd] = zeror;
    yi[nd] = zeroi;
    ++(*nz);
    --nd;
    if (nd == 0) {
	goto L110;
    }
    nuf = zuoik(zr, zi, fnu, *kode, 1, nd, &yr[1], &yi[1], *tol, *elim, *alim);
    if (nuf < 0) {
	goto L140;
    }
    nd -= nuf;
    *nz += nuf;
    if (nd == 0) {
	goto L110;
    }
    fn = *fnu + (double) ((float) (nd - 1));
    if (fn < *fnul) {
	goto L130;
    }
/*      fn = cidi
      j = nuf + 1
      k = mod(j,4) + 1
      s1r = cipr(k)
      s1i = cipi(k)
      if (fn < 0) s1i = -s1i
      str = c2r*s1r - c2i*s1i
      c2i = c2r*s1i + c2i*s1r
      c2r = str */
    in = inu + nd - 1;
    in = in % 4 + 1;
    c2r = car * cipr[in - 1] - sar * cipi[in - 1];
    c2i = car * cipi[in - 1] + sar * cipr[in - 1];
    if (*zi <= 0.) {
	c2i = -c2i;
    }
    goto L40;

L130: *nlast = nd; return;

L140: *nz = -1; return;

L150:
    if (rs1 > 0.) {
	goto L140;
    }
    *nz = *n;
    for (int i = 1; i <= *n; ++i) {
	yr[i] = zeror;
	yi[i] = zeroi;
    }
    return;
} /* zuni2_ */
