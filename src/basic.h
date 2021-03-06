#ifndef GSL_VS_R_H
#define GSL_VS_R_H 1

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <errno.h>
#define NDEBUG 1
#include <assert.h>

#define EXTERN extern "C"
#define PRINTF Rprintf
#define RF_NAN NA_REAL 
#define RF_NEGINF R_NegInf
#define RF_INF R_PosInf
#define RF_ISNA ISNAN
#define GAUSS_RANDOM(SIGMA) rnorm(0.0, SIGMA)
#define UNIFORM_RANDOM unif_rand()
#define POISSON_RANDOM(x) rpois(x)
#define SQRT2 M_SQRT2
#define SQRTPI M_SQRT_PI
#define INVPI M_1_PI
#define RF_M_SQRT_3 M_SQRT_3
//#define TWOPI M_2PI
#define PIHALF M_PI_2 
#define T_PI M_2_PI
#define TWOPI 6.283185307179586476925286766559
#define INVLOG2 1.442695040888963
#define INVSQRTTWO 0.70710678118654752440084436210
#define INVSQRTTWOPI 0.39894228040143270286
#define SQRTTWOPI 2.5066282746310002416
#define MINUSINVLOG005 0.3338082006953340674649
#define SQRTINVLOG005 0.5777613700268771079749
#define LOG05 -0.69314718055994528623
#define INFTY 1e9

#endif /* GSL_VS_R_H */


