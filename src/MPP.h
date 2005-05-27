#ifndef MPP_H
#define MPP_H 1
#include "GSLvsR.h"
#include "MarkedPointProcess.h"

#define NUMBER_L_NORMS 5 /* L_infty, L_2, L_1, L_robust L_anitrobust1 */
#define NUMBERWEIGHTS  7 /* id, 1/(cumsum # in previous bins), 1/sqrt(cumsum), 
			    #, sqrt(), var(bin), ??? */
#define NUMBERTESTS 38   /* NUMBER_L_NORMS * NUMBERWEIGHTS + 3 */
#define INT_EPSILON 0.00000001


#define NN_FACTOR 0
#define NN_MAX 1 /* last number plus 1 */

#define RC_NR 0
#define RC_RADIUS 1
#define RC_HEIGHT 2
#define RC_MAX 3 /* last number plus 1 */

#define MPPERR_NODATA 801 /* ncol of data matrix is zero */
#define MPPERR_PARAMS 802 /* number of parameters do not match required one */
#define MPPERR_POINTS 803 /* number of locations is too small -- should be dealt otherwise in future */
#define MPPERR_MEMALLOC 804
#define MPPERR_MODELNR 805 /* incorrect model number */
#define MPPERR_COORD 806   /* coordinates are NULL vector */
#define MPPERR_ 80
#define NOERROR 0

#define MPP_MAXCHAR 18
#define MPP_MODELS  3   /* maximal total number of models */
typedef void (*mpp_fct)(double *, double *, int*, double*, double*, int*,
			double*, double*);
typedef struct mpp_type {
  char name[MPP_MAXCHAR];
  int nparam, factorpos;
  mpp_fct fct;
}  mpp_type;


extern int currentmppNr;
extern mpp_type mpp_model[MPP_MODELS];
extern bool MPP_BARNARDBESAG;
extern int MPP_PRINTLEVEL;


void InitMPPModelList();
void MPPErrorMessage(int error);

#endif /* MPP_H */









