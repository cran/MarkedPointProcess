/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for simulation and analysis of marked point processes:
    analysis of mpp

 Copyright (C) 2004 -- 2006 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/



/*  TO DO:  **************************
future problem: *theory* on what is an good/optimal test statistic 
                in the MCMC test?
*******************
 
*/



/* 

STRUCTURE:

   input data:
       x : matrix of ``row'' individuals 
       data : [ (multivariate data) columns ] repeated rep times; 
              rep chosen arbitrarily
       row : number of indiviuals of the species

    n species of row[i] individuals containing
       in x[i]: d rows of coordinates (dimension d)
       in data[i]  col rows of (multivariate) data
   
   Eres = n species conditioned on n other species of col data partioned into
          lb bins
   SDres= n species conditioned on n other species of the lower triangular
          matrix (including diagonal) of the multivariate data partioned 
	  into lb bins
   GAMres=lower triangular matrix of (n*col) multivariate data partioned into 
          lb bins
   The number of elements in the E-,SD-,and GAM-BIN have the same structure as
   the corresponding result arrays

   Test= same structure as the result arrays except that the lb bins are 
         replaced by a certain number (NUMBERTESTS) of test statistics
*/

/*
What is left to do:
  * SD raushauen oder auf Korrelation im multivariaten erweitern.
  * Isham's correlation function programmieren
  * KMM ist noch durch E(0)^2 zu dividieren
*/


#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
// #include "RFsimu.h"
// never uncomment the following line
#include "MPP.h"
#include "auxiliary.h"

#define SLD sizeof(long double)

/* weight # raushauen . */
#define XXtolerance 1E-13
#define epsilonvariance 1E-7
#define epsilonvariance2 1E-14

#define MPP_INF 1E50;


typedef enum TestName { ETest, VTest, STest, NoTest } TestName;
TestName testname;
int MPP_PRINTLEVEL;

//  test(E,EBIN,&lb,&n,&col,&rep,p,standarddev,ETEST);
//  test(VAR,VARBIN,&lb,&n,&collowtri,&rep,p,varstdev,VARTEST);
void test(double *E, int *N, long double *sigma, int *lb, int *species, 
	  int *col, int *rep, double *p, long double *sd, double *statistic)
  /*
    E  : E, VAR, or SQ
    N  : number of points in a bins
    sigma   : the (unnormed) besag-diggle weights
    lb : number of bins
    species : number of species
    col     : number of colums; col-variate marks
    rep     : times the data are repeatedly measured or simulated
    p  : a probability parameter for the robust estimation
    sd : the standard deviation of a single mark
         of the whole population of a single species (for each repetition)
    statistic (out) : the test statistics for the currently 7 x 9 (+1) tests
  */
{
  int basiczaehler, zaehler; ///
  double *X, *X2, max, min, *Weight, sumweight, dummy, q, twoq, q2, robustX, 
    sumsqrtweight, sumcumsqrtweight, sumcumweight, sumsqrtcumweight, sumsigma;
  int i, spec, condspec, c, nrstatistic, nW;
  long segment, segtest, LBspec2COL, LBspec2COLrep, NTspec2COL, LBcol, 
    LBcolSPEC, segsd, 
    eseg, e0seg, sseg, zerosegment, statsegment, repsegment, weightNr, nrweight, 
    colSpecies, repsdseg;
  double a1, //a2, a3, a4, a5, 
    b1, //b2, b3, b4, b5, 
    antirob1//, antirob2, antirob3, antirob4, antirob5
    ;
  int NUMBERWEIGHTSxNUMBER_L_NORMS = NUMBERWEIGHTS * NUMBER_L_NORMS;
  double invlbM1 = 1.0/ ((double) *lb -1.0); 
  X = (double*) malloc(sizeof(double) * *lb);
  X2 = (double*) malloc(sizeof(double) * *lb);
  Weight = (double*) malloc(sizeof(double) * (*lb-1) * NUMBERWEIGHTS);
  LBspec2COL = *col * *species * *species * *lb; 
  LBspec2COLrep = *rep * LBspec2COL;
  NTspec2COL = *col * *species * *species * NUMBERTESTS;
  colSpecies = *col * *species;
  LBcol = *lb * *col;
  LBcolSPEC = LBcol * (*species +1);
  segment         /* pointer to Ebin_ij[0] and E_ij[0] */ 
    = zerosegment /* pointer to E_ii(0) */
    = segtest     /* pointer to tests for E_ij */ 
    = segsd       /* pointer to estimated standarddeviation of E_ii(0) */ 
    = 0;
  basiczaehler=0; ///
  for(spec=0; spec<*species; spec++, zerosegment+=LBcolSPEC, segsd+=*col) {      
    for(condspec=0; condspec<*species; condspec++, zerosegment-=LBcol) {
      for (c=0; c<*col; 
	   c++, segment+=*lb, segtest+=NUMBERTESTS, zerosegment+=*lb) {
	/* bins might be empty because of missing values. Otherwise the following
	   code could have been placed at the very beginning 
           however, treating of empty bins is not programmed yet
	*/       
	basiczaehler++; ///
	sumsqrtweight=sqrt(sumweight = (double) N[zerosegment]); 
	// Ebin(0) as starting value ...
	sumsqrtcumweight=sumcumweight = sumcumsqrtweight =0.0;
	nrweight = 0;
	for (i=1; i< *lb; i++) {// 1 not 0!; it is not necessary to consider E(0)
	  /* 
	     0 : constant; weights are 1 / (#bins - 1)
	     1 : 1.0 / sum^{i-1} #
	     2 : sqrt ( 1.0 / sum^{i-1} # ) 
	     3 : 1.0 / sum^{i-1} sqrt #
	     4 : #
	     5 : sqrt #
	     6 : sd (besag-diggle); first constant weights
	  */
	  Weight[nrweight++]= invlbM1;  // weight[0] , nrweight =1 afterwards
	  sumcumweight     += (Weight[nrweight]=1.0/sumweight); 
	  //weight[1]=1, nrweight = 1
          sumsqrtcumweight += (Weight[nrweight+1]=sqrt(Weight[nrweight]));//2, 1
	  nrweight+=2;              //  , 3
	  sumcumsqrtweight += (Weight[nrweight++]=1.0/sumsqrtweight);   // 3, 4
          sumweight        += (Weight[nrweight]=(double) N[segment+i]); // 4, 4
	  sumsqrtweight += (Weight[nrweight+1]=sqrt(Weight[nrweight])); // 5, 4
	  nrweight += 3; //  , 7 (( +2 to catch up; +1 to jump the besag weights)
	}
	sumweight     -= (double) N[zerosegment];
	sumsqrtweight -= sqrt((double) N[zerosegment]);
	sumcumweight     = 1.0/sumcumweight;
	sumsqrtcumweight = 1.0/sumsqrtcumweight;
	sumcumsqrtweight = 1.0/sumcumsqrtweight;
	sumweight        = 1.0/sumweight;
	sumsqrtweight    = 1.0/sumsqrtweight;
	
	
	nrweight = 0; 
	for (i=1;i<*lb;i++) { 
	  nrweight++; /* jump over weight==1.0 */
	  Weight[nrweight++] *= sumcumweight; 	 
	  Weight[nrweight++] *= sumsqrtcumweight; 
	  Weight[nrweight++] *= sumcumsqrtweight; 
	  Weight[nrweight++] *= sumweight;	 
	  Weight[nrweight++] *= sumsqrtweight;   
	  nrweight++; // sumsigma, besag-diggle
	}

	zaehler = basiczaehler;
	for (repsdseg=statsegment=repsegment=0; repsegment<LBspec2COLrep;
	     repsegment+=LBspec2COL, statsegment+=NTspec2COL, 
	       repsdseg+=colSpecies
	       ,zaehler += *col * *species * *species///
	     ) {
	  //  coefficients for antirobust estimators
	  dummy = sd[c + segsd + repsdseg];
	  if (dummy<0) {PRINTF("dummy ***** =%e\n",dummy); assert(false);}
          a1 = 0.10 * dummy;    b1 = 0.25 / a1;
	  //a2 = 0.20 * dummy;    b2 = 0.25 / a2;
          //a3 = 0.30 * dummy;    b3 = 0.25 / a3;
          //a4 = 0.40 * dummy;    b4 = 0.25 / a4;
          //a5 = 0.50 * dummy;    b5 = 0.25 / a5;
				    
	  eseg  = repsegment + segment;
	  e0seg = repsegment + zerosegment;
	  sumsigma = 0.0; // actually sum of inverse of sigma
	  for (i=1; i < *lb; i++) { // 1 not 0!
	    if (RF_ISNA(E[eseg+i]) || RF_ISNA(E[e0seg])) X2[i] = X[i] = RF_NAN;
	    else { 
	      X[i] = fabs(E[eseg+i]-E[e0seg]); //!! see next lines  
	      /* 
		 E[e0seg] is E_ii(0)
		 therefore zerosegment is reset after c-loop
		 and zerosegment is increased by  factor * (spec+ONE)
	      */
	      X2[i] = X[i] * X[i];
	      if (!RF_ISNA(sigma[eseg+i])) sumsigma += 1.0 / sigma[eseg+i];
	    }
	  }
	  X[0] = 0.0;

	  { // besag weights
	    int segment;
	    segment = NUMBERWEIGHTS-1;
	    sumsigma = 1 / sumsigma;
	    for (i=1; i < *lb; i++) { // 1 not 0!
	      Weight[segment] = 
		(RF_ISNA(sigma[eseg+i])) ? 0 : sumsigma / sigma[eseg+i];
		/*
		  berechnung des zaehlers:
		 */
	      segment += NUMBERWEIGHTS;
	    }
	  }

	  q  = quantile(X,*lb,*p);

	  q2 = q * q;
	  twoq = 2.0 * q;
	  
	  sseg = 
	    statsegment /* rep segment of tests */ 
	    + segtest /* pointer to test of E_ij */;
	  for (nW=0; nW<NUMBERWEIGHTSxNUMBER_L_NORMS; nW+=NUMBER_L_NORMS) 
	    statistic[sseg+nW] = 0.0;
          if (RF_ISNA(X[0])) {PRINTF("NAN: %f",X[0]);} 
	  if (RF_ISNA(E[eseg])) { max = RF_NEGINF; min = RF_INF;} 
	  else { max=min=E[eseg]; } 
	  for (weightNr=0, i=1; i<*lb; i++) { // run through all bins
	    if (!(RF_ISNA(X[i]))) {
	      if (E[eseg+i]<min) min=E[eseg+i];
	      if (E[eseg+i]>max) max=E[eseg+i]; 
	      nrstatistic = sseg;	    
	      robustX = X[i] < q ? X2[i] : twoq * X[i] - q2;
	      antirob1= X[i] < a1 ? X[i] : b1 * (dummy=X[i]+a1) * dummy;
	      //antirob2= X[i] < a2 ? X[i] : b2 * (dummy=X[i]+a2) * dummy;
	      //antirob3= X[i] < a3 ? X[i] : b3 * (dummy=X[i]+a3) * dummy;
	      //antirob4= X[i] < a4 ? X[i] : b4 * (dummy=X[i]+a4) * dummy;
	      //antirob5= X[i] < a5 ? X[i] : b5 * (dummy=X[i]+a5) * dummy;
	      
	      for (nW=0; nW < NUMBERWEIGHTS; nW++) {	      
		double wx,wx2,aW; 
		if (RF_ISNA(wx  =  X[i] * Weight[weightNr])) wx=0;
		if (RF_ISNA(wx2 = X2[i] * Weight[weightNr])) wx2=0;
		if (wx>statistic[nrstatistic]) {statistic[nrstatistic]=wx;} // 1
		nrstatistic++;
		statistic[nrstatistic++] += wx2;  // 2
		statistic[nrstatistic++] += wx;   // 3
		if (RF_ISNA(aW=robustX * Weight[weightNr])) aW=0;
		statistic[nrstatistic++] += aW;   // 4		
		if (RF_ISNA(aW=antirob1 * Weight[weightNr])) aW=0;
		statistic[nrstatistic++] += aW;   // 5
		//if (RF_ISNA(aW=antirob2 * Weight[weightNr])) aW=0;
		//statistic[nrstatistic++] += aW; 
		//if (RF_ISNA(aW=antirob3 * Weight[weightNr])) aW=0;
		//statistic[nrstatistic++] += aW; 
		//if (RF_ISNA(aW=antirob4 * Weight[weightNr])) aW=0;
		//statistic[nrstatistic++] += aW;
		//if (RF_ISNA(aW=antirob5 * Weight[weightNr])) aW=0;
		//statistic[nrstatistic++] += aW; 
		weightNr++;
	      }
	    } else {
	      weightNr += NUMBERWEIGHTS;
	    }
	  }
	  statistic[sseg + NUMBERWEIGHTSxNUMBER_L_NORMS] = max - min; 
	  // last but third one !!

	  // if ((testname==ETest) && (zaehler==-1230000)) 

	} /* rep */
      } /* col */
    }  /* conditional species */
  } /* spec */
  free(X);
  free(X2);
  free(Weight);
}

int mcf_internal(double *E, double *ETEST, int *EBIN,
		 double *VAR, double *VARTEST,
		 double *SQ, double *SQTEST, int *VARBIN,
		 double *KMM, int *KMMBIN,
		 double *GAM, int *GAMBIN,
		 double *p, double *bin, int lb, int dim, int n,
		 int col, int rep, double **X, double **DATA, int *ROW,
		 bool staticchoice
		 )

     /* NOTE : for multivariate data, the test statistics by the anti-robust 
	functions are not correctly calculated as varstdev is not calculated 
	for covariances !

     double *E,*ETEST;     results: E function and the test statistic 
			   ETEST : Matrix of all the ETEST statistics 
     int    *EBIN;         vector with the number of points in each bin 
     double *VAR,*VARTEST;   V- E^2 
     int    *VARBIN;     
     double *KMM; // *KMMTEST; stoyan's kmm 
     int    *KMMBIN;    
     double *GAM; // *GAMTEST; mark variogram 
     int    *GAMBIN;
     
     double  *SQ,*SQTEST;   sqrt(V - E^2) 
     
     double *p,            p-value for the "robust" functions in the test 
            *bin;          boundaries for the bins 
     int lb,dim,n,col,rep;
        lb  : length(bin)
	dim : dimension of the space, usually 2
	n   : number of species
	col : number of variables for multivariate marks
	rep : number of repeated measurement of the data
     
     double **X,**DATA;  lists of coordinates ([ROW x 2]-matrix) and data, 
                         respectively 
     int *ROW;           list with the number of marked points 
     bool staticchoice : good idea for debugging
     */
{
  int spec[2], row[2], collowtri, Xsumcol, sumcol, d, low, halflb, swap,
    antiswap, intsqrtn, member, *classes, *sort, ind, intertotal, error, end;
  long endfor, k, kk, i, j, endfor2, intdummy, Etestseg[2], Etotal, 
    VARtotal, Datatotal, Datatotal1, Etesttotal,
    GAMtotal, meantotal, meantotalseg, meantotalrep, 
    vartotal, vartotalseg,
    LBsumcolN, LBsumcolNrep,nnLBcollowtri,  LBcollowtriN2rep, LBsumcoltri, 
    LBsumcoltriRep,
    NTsumcolN,  NTsumcolNrep, NTcollowtriN2rep, 
    *ROWcol,
    GAMresseg[2], Eresseg[2], VARresseg[2], seg[2], seg0, seg1, l, var_l,
    var0, var1, LBcolN, LBcollowtriN, LBcol2, 
    NTcolN;
  long double dummy, *mean, *standarddev, *varstdev, *varExp, *sdstdev, *mom2,
    *mom22, *mom21, *mom12, *intersum1, *intersum2;
  double *pos;

  classes = sort = NULL;
  ROWcol = NULL;
  mean = standarddev = varstdev = varExp = sdstdev = mom2 = mom22 = mom21 = 
    mom12 =  intersum1 = intersum2 = NULL;
  pos =NULL;
  
  sumcol = n * col;
  collowtri= (col * (col+1)) / 2 ;
  halflb =  lb/2;
  LBcolN = lb * col * n;
  NTcolN = NUMBERTESTS * col * n;
  LBcollowtriN = lb * collowtri * n;
  LBcol2 = lb * col * col;
  LBsumcolN = n * sumcol * lb;
  LBsumcolNrep =LBsumcolN * rep;
  nnLBcollowtri = LBcollowtriN * n;
  LBcollowtriN2rep = nnLBcollowtri * rep;
  LBsumcoltri = (sumcol * (sumcol+1) * lb)/2; 
  LBsumcoltriRep =  LBsumcoltri * rep;
  NTsumcolN = n * sumcol * NUMBERTESTS; 
  NTsumcolNrep = NTsumcolN * rep;
  NTcollowtriN2rep = NUMBERTESTS * collowtri * n * n * rep;
  Xsumcol = (col * sumcol - (col*(col-1))/2) * lb; 

  meantotalseg = n * col; 
  vartotalseg = n * collowtri;
  meantotalrep = rep * meantotalseg;
  if (((ROWcol = (long *) malloc(sizeof(long) * n)) == NULL) ||
      ((mom2 = (long double *) malloc(SLD * LBsumcolNrep)) == NULL) || 
      // for variation of E
      ((mom22 = (long double *) malloc(SLD * LBcollowtriN2rep)) == NULL) ||
      // for variation of VAR
      ((mom21 = (long double *) malloc(SLD * LBcollowtriN2rep)) == NULL) ||
      // for variation of VAR
      ((mom12 = (long double *) malloc(SLD * LBcollowtriN2rep)) == NULL) ||
      // for variation of VAR
      ((mean  = (long double *) calloc(meantotalrep, SLD )) == NULL) ||
      ((standarddev = (long double *) calloc(meantotalrep, SLD))== NULL) ||
      // standarddev for mean
      ((varstdev = (long double *) calloc(rep * vartotalseg, SLD)) == NULL) || 
      // standarddev for V(ar), array is correct
      ((sdstdev  = (long double *) calloc(rep * vartotalseg, SLD)) == NULL) || 
      // standarddev for SD
      ((varExp   = (long double *) calloc(meantotalrep, SLD)) == NULL)
      // expectation of V(ar)
      ) { error=MPPERR_MEMALLOC; goto ErrorHandling; }

  /* 
     we have to estimate the standarddeviation for (mean), Var, SD;
     the standard deviation is used to calculate the anti-robust tests, 
     and there the "breakpoint" between linear and quadratic punishment
  
     To this end, a variant of the jackknife technique is used: the data are 
     partioned and put into z classes, where z is approximately sqrt(ROW). Then 
     for each of the classes the mean, variance, and standard deviation is 
     calculated; the variation among the classes, devided by approximately 
     sqrt(number of members in each class), gives the estimated variation of
     mean, V, SD
  */
  for (i=0; i<n; i++) {ROWcol[i] =  ROW[i] * col;}
  for (end=col, var_l=l=0, i=0; i<n; i++, end += col) {
    intsqrtn= (int) (sqrt(ROW[i])+INT_EPSILON);// sqrt(n) classes
    member = ROW[i] / intsqrtn;               // with about sqrt(n) elements each
    if ( ((classes = (int*) malloc(sizeof(int) * (intsqrtn+1))) == NULL) ||
	 ((pos = (double*) malloc(sizeof(double) * ROW[i])) == NULL) ||
	 ((sort = (int*) malloc(sizeof(int) * ROW[i])) == NULL) ||
	 ((intersum1 = (long double*) malloc(SLD *  rep)) == NULL) ||
	 ((intersum2 = (long double*) malloc(SLD *  rep)) == NULL)
	 ) { error=MPPERR_MEMALLOC; goto ErrorHandling; }
    kk = ROW[i] % intsqrtn;  // the remaining k elements are put to the first k
    //                          classes             
    classes[0] = 0;
    for (j=1; j<=intsqrtn; j++) {
      classes[j] = classes[j-1] + member;
      if (kk>0) { classes[j]++; kk--;}
    }
    
    for (j=0; j<ROW[i]; j++) {  
      pos[j]= (staticchoice) ? ((double) j) : (UNIFORM_RANDOM);   
      sort[j]=j; 
    }
    // missing values einarbeiten ! also letztendlich die Klasseneinteilung 
    // besser waehlen !
    // Im folgenden Aufruf wird angenommen, dass keine missing values auftreten!!
    
    orderdouble(pos, sort, 0,ROW[i]-1);
   
    for (j=0; l<end; j+=ROW[i], var_l+=end-l, l++) { //col !
      // here somehow the second loop for calculating covariances is missing 
      for (intdummy=0, k=0; k<intsqrtn; k++) { //classes
	for (intertotal=0; intertotal<rep; intertotal++) {
	  intersum1[intertotal] = intersum2[intertotal] = 0.0;
	}
	for (ind=classes[k]; ind<classes[k+1]; ind++){// within each classes 
	  register double dummy;
	  if (!RF_ISNA(DATA[i][j+sort[ind]])) {// sort[ind]==random permutation
	    for (Datatotal=meantotal=intertotal=0; meantotal<meantotalrep;
		 meantotal+=meantotalseg, Datatotal+=ROWcol[i], intertotal++) {
	      mean[l+meantotal] += (dummy=DATA[i][j+sort[ind]+Datatotal]);  
	      intersum1[intertotal] += dummy;
	      standarddev[l+meantotal] += (dummy *= dummy);
	      intersum2[intertotal] += dummy;	      
	    }
	    intdummy++;
	  }
	}  
	for (vartotal=meantotal=intertotal=0; intertotal<rep; 
	     intertotal++, meantotal+=meantotalseg, vartotal+=vartotalseg){
	  dummy = (intersum2[intertotal] - 
		   intersum1[intertotal] * intersum1[intertotal] 
		   / (classes[k+1]-classes[k])) / (classes[k+1]-classes[k]-1.0);
	  if (dummy<0.0) dummy=0.0;
	  varExp[l + meantotal]   += dummy;
	  varstdev[var_l + vartotal] += dummy * dummy; 
	  sdstdev [var_l + vartotal] += sqrt(dummy);	
	}
      } // k (classes)
      for (vartotal=meantotal=0; meantotal<meantotalrep; 
	   meantotal+=meantotalseg, vartotal+=vartotalseg) {
	double doubledummy, doublesqrtn, doublesqrtnSq;
	doubledummy = (double) intdummy;
	doublesqrtn = (double) intsqrtn;
	doublesqrtnSq = (doublesqrtn - 1.0) * doublesqrtn;
	mean[l + meantotal] /= doubledummy;  

	standarddev[l + meantotal] -= 
	  doubledummy * mean[l + meantotal] * mean[l + meantotal];
	standarddev[l + meantotal] = (standarddev[l + meantotal]>0) ? 
	  sqrt(standarddev[l + meantotal] / (doubledummy - 1.0)) : 0.0;	
        varstdev[var_l + vartotal] -= 
	  varExp[l + meantotal] * varExp[l + meantotal] / doublesqrtn;
	varstdev[var_l + vartotal] = (varstdev[var_l + vartotal]>0) ? 
	  sqrt(varstdev[var_l + vartotal] / doublesqrtnSq) : 0.0 ;	
	sdstdev[var_l + vartotal] = varExp[l + meantotal] - 
	  sdstdev[var_l + vartotal] * sdstdev[var_l + vartotal] / doublesqrtn;
        sdstdev[var_l + vartotal] = (sdstdev[var_l + vartotal]>0) ?
	  sqrt( sdstdev[var_l + vartotal] / doublesqrtnSq) : 0.0;	
	assert(standarddev[l+meantotal]>=0);
      }
    } // for j (col)
    free(classes); free(pos); free(sort); free (intersum1); free(intersum2);
    classes =  sort = NULL;
    intersum1 = intersum2 = NULL;
    pos = NULL;
  }
  free(varExp);
  varExp = NULL;
  /* end estimation of the standard deviations */
 

  for (i=0; i<LBsumcolNrep; i++) {mom2[i]=E[i]=0.0;}
  for (i=0; i<LBsumcolN; i++) {EBIN[i]=0;}
  for (i=0; i<LBcollowtriN2rep; i++) {mom22[i]=mom21[i]=mom12[i]=VAR[i]=0.0;}
  for (i=0; i<nnLBcollowtri; i++) {VARBIN[i]=0;}
  for (i=0; i<LBsumcoltriRep; i++) {KMM[i]=GAM[i]=0.0;}
  for (i=0; i<LBsumcoltri; i++) {KMMBIN[i]=GAMBIN[i]=0;}  
  for (i=0; i<NTsumcolNrep; i++) {ETEST[i]=0.0; }
  for (i=0; i<NTcollowtriN2rep; i++) {SQTEST[i]=VARTEST[i]=0.0;}
 
  // ****** summing up ******
  Eresseg[0] = 0;
  VARresseg[0] = 0;
  GAMresseg[0] = 0;
  Etestseg[0] = NUMBERTESTS-2; /* there are 2 tests on the nonbinned E fct */ 
  for (spec[0]=0;spec[0]<n;spec[0]++) { /* species 1 */
    if (MPP_PRINTLEVEL>4) PRINTF(" spec=%d",spec[0]);
    row[0] = ROW[spec[0]];
    Eresseg[1] = Eresseg[0];
    VARresseg[1] = VARresseg[0];
    Etestseg[1] = Etestseg[0];
    /// GAMresseg[1] = 0; not needed.
    for (spec[1]=spec[0];spec[1]<n;spec[1]++) { /* species 2 */
      row[1]=ROW[spec[1]];
      for (i=0; i<row[0]; i++) { /* locations of species 1 */
	for (j=0; j<row[1]; j++) { /* locations of species 2 */
// printf("%d %d\n", i,j);
	  double dist;
 	  dist=0;
          seg[0] = i;
          seg[1] = j;
	  for (d=0;d<dim;d++) {
	    dist += (dummy=X[spec[0]][seg[0]]-X[spec[1]][seg[1]]) * dummy;
	    seg[0] += row[0]; 
	    seg[1] += row[1];
	  }
          dist = sqrt(dist);
	  if (dist<bin[lb]) {
	    { 
	      // search which bin u[i] is in
	      register int up, cur;
	      low=0; up= lb-1; cur= halflb;
	      while (low!=up) {
//  printf("%d %d %d %e %e\n", low, up, cur, bin[cur], dist);
		if (dist> bin[cur]) {low=cur;} else {up=cur-1;} /* ( * ; * ] */
		cur=(up+low+1)/2;
	      }
	    } 	
    
	    /* calculate E location i given i,j, and j given j,i*/
	    seg[0] = i;
	    seg[1] = j;
	    endfor = (int) (spec[0]!=spec[1]);
	    for (swap = 0;swap<=endfor;swap++) {
	      antiswap = 1 -swap;
	      seg0 = seg[swap];
	      k=0;
	      for (var0=0; var0<col; var0++) {
//   printf("swap, var %d %d %d \n", swap, endfor, var0, col);
		if (!RF_ISNA(DATA[spec[swap]][seg0])) {
		  EBIN[intdummy=low+(var0+spec[antiswap]*col)*lb+Eresseg[swap]]++;
		  for (Datatotal=Etotal=0; Etotal<LBsumcolNrep; // for all repets
		       Etotal+=LBsumcolN, Datatotal+=ROWcol[spec[swap]]) {
		    E[intdummy + Etotal] += DATA[spec[swap]][seg0+Datatotal];
		    // the following line could be programmed more 
		    // efficiently, using the calculations of AR, finally
		    // by just copying the relevant data
		    mom2[intdummy + Etotal] += 
			  (long double)DATA[spec[swap]][seg0+Datatotal] *
			(long double)DATA[spec[swap]][seg0+Datatotal];
		  }
		
// von hier bis
		  intdummy=(var0+spec[antiswap]*col)*NUMBERTESTS+Etestseg[swap];
		  for (meantotal=Datatotal=Etesttotal=0; Etesttotal<NTsumcolNrep;
		       Etesttotal += NTsumcolN, Datatotal += ROWcol[spec[swap]],
			   meantotal += meantotalseg){
		    ETEST[intdummy+1+Etesttotal] +=
			(dummy = fabs(DATA[spec[swap]][seg0+Datatotal] -
				      mean[spec[swap]*col + var0 + meantotal])); 
		    ETEST[intdummy+Etesttotal] += dummy *dummy;
#ifdef DEBUG
		    ETEST[intdummy+1+Etesttotal]=ETEST[intdummy+Etesttotal]=RF_NAN;
#endif
		  }
// ...hier die 2 extra schaetzer, die sich letztendlich sehr schlecht verhalten
		}

		/* calculate VAR */
	        seg1 = seg0;
		for (var1=var0; var1<col; var1++) {
		  intdummy=low+(k+spec[antiswap]*collowtri)*lb+VARresseg[swap];
  		  if (!RF_ISNA(DATA[spec[swap]][seg0]) && 
		      !RF_ISNA(DATA[spec[swap]][seg1])) {
		    VARBIN[intdummy]++;
		    for (Datatotal=VARtotal=0;VARtotal<LBcollowtriN2rep;
		         VARtotal+=nnLBcollowtri, Datatotal+=ROWcol[spec[swap]]){
		      register long double dummy;
		      dummy = (long double) DATA[spec[swap]][seg0+Datatotal] * 
			  (long double) DATA[spec[swap]][seg1+Datatotal];
		      //  !!!!!!!!! not checked yet if it seg0 and seg1 correspond
		      // to the indices later on in the final calculations
		      mom22[intdummy+VARtotal] += dummy * dummy;
		      mom21[intdummy+VARtotal] += 
			  dummy * (long double) DATA[spec[swap]][seg0+Datatotal];
		      mom12[intdummy+VARtotal] += 
			  dummy * (long double) DATA[spec[swap]][seg1+Datatotal];
		      VAR[intdummy+VARtotal]+= dummy;
		    }
		  }// RF_ISNA
		  k++;
		  seg1 += row[swap]; // not antiswap !
		} // for var1=var0
		seg0 += row[swap];
	      }
	    }    
	    seg0 = i; seg1 = j;
	    /* calculate Variogram gamma*/
            /* ordering of the points :

	       |v1=v0 v1=0      v1=0  |v1=v0        v1=0  ||v1=v0 v1=0 
               | to    to        to   | to           to   || to    to  
               |col-1 col-1     col-1 |col-1        col-1 ||col-1 col-1
               |-----|----|....|------|-----|....|--------||-----|-----|....
                s0=0                                        s0=1  
                v0=0                   v0=1       v1=col-1  v0=0
                s1=0  s1=1      s1=n-1 s1=1       s1=n-1    s1=1  s1=2
		
	     Therefore, the first block s1=1 contains
	        sumcol+(sumcol-1)+...+(sumcol-ncol+1)=ncol*sumcol-ncol*(ncol-1)/2
	      elements; block s1=i:
                ncol*(n+1-i)*ncol - ncol*(ncol-1)/2 = 
                                     ncol*sumcol - (i-1)*ncol^2 - ncol*(ncol-1)/2

              (( total : 
	      n^2*ncol^2 - ncol^2*n*(n-1)/2  - n*ncol*(ncol-1)/2 
	        = n^2*ncol^2 - ncol^2*n*n/2 +  ncol^2*n/2 - n*ncol^2/2 + n*ncol/2
		= n^2 ncol^2/2 + n*ncol/2 
		= n * ncol * (n * ncol +1) /2
	      ))  
	    */
	    for (var0=0; var0<col; var0++) { 
	      int Var0Segment,startfor; // can be speeded up
	      //Var0Segment = var0*(n-spec[0])*col - (var0*(var0-1))/2;
	      /* [col+(n-s1-1)*col] + [col-1+(n-s1-1)*col] + .. 
		                                     + [col-(var-1)+(n-s1-1col)]
		 = var0*(n-s1)*col - var0*(var0-1)/2 
	      */
	      //SpecSegment=(spec[1]-spec[0])*col-var0; 
	      //  (col-var0) + (s1-s9-1)*col;
	      if (spec[0]!=spec[1]){seg[0]=i;    seg[1]=j;    startfor=0;} 
	      // col*col fcts m_i*m_j to average over
	      else                 {seg[0]=seg0; seg[1]=seg1; startfor=var0;}
	      // col(col+1)/2 fcts only because of symmetry 
	      Var0Segment = (spec[1]-spec[0]+ var0 * (n-spec[0])) * col 
		- (var0*(var0+1))/2;
	      for (var1=startfor; var1<col; var1++) { //!
		if (!RF_ISNA(DATA[spec[0]][seg0]) && 
		    !RF_ISNA(DATA[spec[1]][seg[1]])) { 
		  KMMBIN[intdummy=low+GAMresseg[0]+(Var0Segment+var1)*lb]++;
		  for (Datatotal1=Datatotal=GAMtotal=0; GAMtotal<LBsumcoltriRep;
		       GAMtotal+=LBsumcoltri,
			 Datatotal+=ROWcol[spec[0]],Datatotal1+=ROWcol[spec[1]]){
		    KMM[intdummy+GAMtotal]+=DATA[spec[0]][seg0+Datatotal]*
			DATA[spec[1]][seg[1]+Datatotal1];
		  }
		  if (!RF_ISNA(DATA[spec[1]][seg1]) && 
		      !RF_ISNA(DATA[spec[0]][seg[0]])) { 
	            GAMBIN[intdummy]++;    
		    for (Datatotal1=Datatotal=GAMtotal=0;
			 GAMtotal<LBsumcoltriRep;
			 GAMtotal+=LBsumcoltri, Datatotal+=ROWcol[spec[0]],
			   Datatotal1+=ROWcol[spec[1]]){ 
		      GAM[intdummy+GAMtotal]+=
			  (DATA[spec[0]][seg0+Datatotal] - 
			   DATA[spec[1]][seg1+Datatotal1]) *
			  (DATA[spec[0]][seg[0]+Datatotal] - 
			   DATA[spec[1]][seg[1]+Datatotal1]);
		    }
		  } // if !RF_ISNA, add, gamma
		}// if !RF_ISNA, kmm
		seg[0] += row[0]; seg[1] += row[1];
	      }
	      seg0 += row[0]; seg1 += row[1];
	    }
	  } /* disr < bin[lb] */
	} /* j */
      } /* i */
      Eresseg[1] +=   LBcolN;
      Etestseg[1] +=  NTcolN;
      VARresseg[1] +=  LBcollowtriN;
    } /* spec[1] */
    Eresseg[0] += LBcolN;
    Etestseg[0] +=  NTcolN;
    VARresseg[0] +=  LBcollowtriN;
    GAMresseg[0] += Xsumcol; 
    Xsumcol -= LBcol2;  
  } /* spec[0] */

//  printf("final calculations\n");

  /* final calculations  */ 
  /* E */
  endfor= lb * sumcol * n; 
  for (i=0; i<endfor; i++) { 
    long double doubleEBINi;
    doubleEBINi = (long double) EBIN[i];
    for (Etotal = 0; Etotal < LBsumcolNrep; Etotal += LBsumcolN){ 
      // sqrt of mom2 will be calculated later on !!	  
      if (EBIN[i]<=1) {
	// mom2[i+Etotal] = RF_INF;  // changed back 31.12.01
        if (EBIN[i]==1) {mom2[i+Etotal] = RF_INF;} // E remains the same
	else { E[i+Etotal] = mom2[i+Etotal] = RF_NAN; } 
      } else {
	E[i+Etotal] /= doubleEBINi; 
	mom2[i+Etotal] = (mom2[i+Etotal] - doubleEBINi * E[i+Etotal] * 
			  E[i+Etotal]) / (doubleEBINi - 1.0);
	if (mom2 [i+Etotal] < 0)  {
	  if (mom2[i+Etotal]<= -XXtolerance) {
	    PRINTF("\n ********* %f %f ***\n\n",doubleEBINi, E[i+Etotal]);
	    PRINTF(" %d: %d  ",i, Etotal);
	    PRINTF(" %e %e %e\n",mom2[i+Etotal], E[i+Etotal], doubleEBINi);
	    assert(false);
	  }
	  if ((MPP_PRINTLEVEL>5) )       
	    PRINTF("\n *********  %e %e %f %f ***", 
		   (double) mom2[i+Etotal], doubleEBINi, E[i+Etotal]);
	  mom2[i+Etotal]=0; 
	}	
      }
    }
  }

  /* VAR */
  Eresseg[0] = 0;
  kk = 0;      
  endfor = n * n;
  for (i=0; i<endfor;i++) {
    for(j=0; j<col; j++) {
      endfor2 = (col-j)*lb;
      for (k=0; k<endfor2;k+=lb) {
	for(l=0; l<lb; l++) {
	  for (Etotal = VARtotal = 0; VARtotal < LBcollowtriN2rep; 
	       Etotal += LBsumcolN, VARtotal += nnLBcollowtri){

	    register long double dummy,ExE, doubleVARBINkk;
	    
	    doubleVARBINkk = (double) VARBIN[kk];
	    ExE = E[Eresseg[0]+l+Etotal] * E[Eresseg[0]+k+l+Etotal];
	    dummy = (VAR[kk+VARtotal] - doubleVARBINkk * ExE) / 
	      (doubleVARBINkk-1.0);

	    // check if multivariate one is ok.!!
	    if (VARBIN[kk]>=3) {
	      long double XX;
	      XX = (mom22[kk+VARtotal] 
		    - 2.0 * mom21[kk+VARtotal] * 
		    (long double) E[Eresseg[0]+k+l+Etotal] 
		    - 2.0 * mom12[kk+VARtotal] * 
		    (long double) E[Eresseg[0]+l+Etotal]
		    + 4.0 *  (long double) VAR[kk+VARtotal] * ExE 
		    + (doubleVARBINkk - 1.0) * mom2[Eresseg[0]+k+l+Etotal] *
		    (long double)  E[Eresseg[0]+l+Etotal] * 
		    (long double) E[Eresseg[0]+l+Etotal]
		    + (doubleVARBINkk - 1.0) * mom2[Eresseg[0]+l+Etotal] * 
		    (long double) E[Eresseg[0]+k+l+Etotal] * 
		    (long double)E[Eresseg[0]+k+l+Etotal] 
		    - doubleVARBINkk * ExE * ExE) / (doubleVARBINkk-2.0) - 
		dummy*dummy;
	      if ((XX<0)) {
		if ((XX < - XXtolerance * dummy * dummy) 
		    &&
		    ((fabs(XX) > XXtolerance)||(fabs(dummy) > XXtolerance))) {
		  PRINTF("bin %d\n",VARBIN[kk]);
		  PRINTF(" %e\n", 
			 (double) ((XX + dummy * dummy) * (doubleVARBINkk-2)));
		  PRINTF(" %4.3f(%4.3f)\n",(double) mom21[kk+VARtotal], 
			 (double) mom21[kk+VARtotal]);
		  PRINTF(" [%e %e, %d, %e %e; %e %e]\n", 
			 (double)mom22[kk+VARtotal], 
			 (double)VAR[kk+VARtotal], VARBIN[kk], (double) ExE,
			 (double)mom2[Eresseg[0]+k+l+Etotal],(double)dummy,
			 sqrt((VAR[kk+VARtotal]- (double) VARBIN[kk]* ExE ) / 
			      (doubleVARBINkk-1.0)));
		  PRINTF("--> %e > %e [%e %e]\n",
			 (double) XX, (double) (- XXtolerance * dummy * dummy),
			 (double)XXtolerance, (double) dummy);
		  assert(false);
		}
		XX=0;
	      }
	      mom22[kk+VARtotal] = sqrt(XX);    
	    } else { 
	      mom22[kk+VARtotal] = RF_NAN;
	    }
	    if (VARBIN[kk]>=2) { // must be after above if statement !!
	      VAR[kk+VARtotal]=dummy;
	      SQ[kk+VARtotal] = sqrt(fabs(dummy));
	      if(dummy<0) {SQ[kk+VARtotal] = -SQ[kk+VARtotal];}
	    } else {
	      VAR[kk+VARtotal] = SQ[kk+VARtotal] = RF_NAN;
	    }
	  }
	  kk++;
	}             
      }
      Eresseg[0]+=lb;
    }
  }

  /* GAM */
  endfor=(sumcol * (sumcol+1) * lb)/2; 
  for (i=0;i<endfor;i++){
    dummy = 1 / (double) GAMBIN[i];
    for (GAMtotal=0; GAMtotal<LBsumcoltriRep; GAMtotal+=LBsumcoltri){
      KMM[i+GAMtotal] *= dummy;
      GAM[i+GAMtotal] *= (dummy * 0.5);
    }
  }

  /*          ********  TESTS ********            */
  if (MPP_PRINTLEVEL>3) PRINTF("\nE TEST ");
   testname = ETest;

  for (i=0;i<LBsumcolNrep;i++) // Standard-Abweichung
    mom2[i] = (mom2[i]<epsilonvariance2) ? epsilonvariance : sqrt(mom2[i]);
  test(E, EBIN, mom2, &lb, &n, &col, &rep, p, standarddev, ETEST);

  if (MPP_PRINTLEVEL>3)  PRINTF("\nVAR TEST ");
  testname = VTest;
  for (i=0; i<LBcollowtriN2rep; i++) // Standard-Abweichung
    if (mom22[i]<epsilonvariance) mom22[i]=epsilonvariance;
  test(VAR, VARBIN, mom22, &lb, &n, &collowtri, &rep, p, varstdev, VARTEST);
  
  if (MPP_PRINTLEVEL>3) PRINTF("\nSQ TEST ");
  for (i=0;i<LBcollowtriN2rep;i++) { 
    assert( (mom22[i]>=0) || RF_ISNA(mom22[i])); 
    if (!RF_ISNA(mom22[i])) mom22[i] = sqrt(mom22[i]); 
  }
  testname = STest;
  test(SQ, VARBIN, mom22, &lb, &n, &collowtri, &rep, p, sdstdev, SQTEST);
 
  free(mean);  
  free(standarddev); 
  free(varstdev);
  free(sdstdev); 
  free(mom2);  
  free(mom22);       
  free(mom21);    
  free(mom12); 
  free(ROWcol);
  return NOERROR;

ErrorHandling:
  if (classes!=NULL) free(classes);
  if (sort!=NULL)    free(sort);
  if (ROWcol!=NULL)  free(ROWcol);
  if (mean!=NULL)    free(mean);
  if (standarddev!=NULL) free(standarddev); 
  if (varstdev!=NULL)    free(varstdev);
  if (varExp!=NULL)  free(varExp);
  if (sdstdev!=NULL) free(sdstdev);
  if (mom2!=NULL)    free(mom2);
  if (mom22!=NULL)   free(mom22);
  if (mom21!=NULL)   free(mom21);
  if (mom12!=NULL)   free(mom12);
  if (pos!=NULL)     free(pos);
  if (intersum1!=NULL)   free(intersum1);
  if (intersum2!=NULL)   free(intersum2);
  return error;
}


void store(int value, int *whereto, int *intestbin, int nestbin) {	  
  int low;
  if ((value>=intestbin[0]) && (value<intestbin[nestbin])) {
    register int up, cur;
    low=0; up= nestbin-1; cur= nestbin/2;
    while (low!=up) {
      if (value >= intestbin[cur]) {low=cur;} else {up=cur-1;} // [ * ; * ) 
      cur=(up+low+1)/2;
    }
    (whereto[low])++;
  } else {
    if (MPP_PRINTLEVEL>0) 
      PRINTF("value (%d) outside inttestbin [%d, %d] \n",
	     value,intestbin[0],intestbin[nestbin]); 
  }
}

void MCtest(int *repet, double *coord, double *data, int *npoints,
	    int *dim, double *simu,  int *PrintLevel,
	    double *bin, int *nbin, double *estbin, int *nestbin,
	    int *Etestposition, int *VARtestposition, int *SQtestposition, 
	    int *Maxtests, int *nmaxtests, int *MAXtestposition,
	    int *error, int *additive, int *staticchoice)
{
    /* repet  : MC repetitions. is usually 99
       coord  : are the fixed coordinates of the mpp (npoints coordinates)
       data   : vector of (univariate) marks
       npoints : number of points
       dim    : spatial dimension
       simu   : gives the simulation results (Barnard type and Gaussrf)
       bin    : the margins of the bins for the E and the VAR function
       nbin   : number of bins! i.e. length(bin)-1
       estbin : usually 0.00,0.01,...,1.00;  margins for the ranks in percent
       nestbin : number of Etest bins, so length(estbin)-1
       nestbin : number of estbins, i.,e length(estbin)-1
           NOTE: if ncol==1 then Etestposition and VARtestposition returns the 
	   rank in Etestposition[0..NUMBERTESTS-1] and 
	                                      VARtestposition[0..NUMBERTESTS-1]
	   the actual values of estbin and nestbin are ignored!
	   Etestposition, VARtestposition
       Etestposition : 
       VarTesttposition :
       SQtestposition   :
       MaxTests (input var) gives the nmaxtests combinations of  E and V tests
           currently unused by R
       MaxTestpos contains the results,
       additive : usually, if data measured only once, MCtest returns the 
                  Etest values is able to overwrite this, i.e., a Etest 
                  matrix is always (expected and) returned;
                  furthermore the entries of the matrix are not set to zero
       staticchoice: if 1 then intervals in mcf_internal are chosen 
             deterministically; 
     */

  double *E,*ETest; 
  int    *Ebin;
  double *VAR,*VARTest; 
  int    *VARbin;
  double *KMM; 
  int    *KMMbin;
  double *GAM; 
  int    *GAMbin;
  double *SQ,*SQTest; ///

  double **X, **DATA, *DUMMY;

   int nn, i, j, k,
    *intestbin, 
    /* the values in estbin are given in quantiles, here the equivalent,
       absolute values are stored. Example: repet=99, estbin=0.0,..,1
       then intestbin 0,1,..,99
    */
    Epos[NUMBERTESTS], VARpos[NUMBERTESTS],SQpos[NUMBERTESTS],*MAXpos; 
  /* they store the actual rank of the data-"E"-statistic in the MC test */
    
  double p;
  int n, col, endrepet;
  long simresN; simresN=0;

 /////////////////
  MPP_PRINTLEVEL = *PrintLevel;
  p = 0.7;// p-value for "robust" function in test 
  n = 1;  // number of species
  col =1; // univariate data only
  ////////////////
  E=ETest=VAR=VARTest=KMM=GAM=SQ=SQTest=NULL;
  X=DATA=NULL;
  DUMMY=NULL;
  Ebin=VARbin=KMMbin=GAMbin=intestbin=MAXpos=NULL;

  if ((intestbin = (int*) malloc(sizeof(int) * (*nestbin+1)))==NULL) {
    *error=MPPERR_MEMALLOC; goto ErrorHandling;}
  nn = *repet + 1;

  GetRNGstate();

  if (*additive) {
    /*
      store Etest result in large matrix, so that position in matrix is 
      equivalent to value otherwise return the Etest value directly
    
      if additive, then the user of R has to ensure that Etest,etc are set
      to zero at the very beginning
      
      additive is necessary, if data are simulated for different coordinates, 
      but where the results in the Etest table should be added up
    */
    for(i=0; i<=*nestbin; i++) {
      if (estbin[i]>1e-10) 
	intestbin[i] = ((int) (estbin[i] * (double) nn - 1e-10)) + 1;        
      else intestbin[i] = 0;
    }
  }

  {
    int double_nn_nbin,  int_nn_nbin, double_nn_NT;
    double_nn_nbin = sizeof(double) * nn * *nbin; 
    int_nn_nbin = sizeof(int) * nn * *nbin; 
    double_nn_NT = sizeof(double) * nn * NUMBERTESTS;    
    if (
	((E      = (double*) malloc(double_nn_nbin))==NULL) ||
	((VAR    = (double*) malloc(double_nn_nbin))==NULL) ||
	((SQ     = (double*) malloc(double_nn_nbin))==NULL) ||
	((KMM    = (double*) malloc(double_nn_nbin))==NULL) ||
	((GAM    = (double*) malloc(double_nn_nbin))==NULL) ||
	((ETest  = (double*) malloc(double_nn_NT))==NULL) ||
	((VARTest= (double*) malloc(double_nn_NT))==NULL) ||
	((SQTest = (double*) malloc(double_nn_NT))==NULL) ||
	((Ebin   = (int*) malloc(int_nn_nbin))==NULL) ||
	((VARbin = (int*) malloc(int_nn_nbin))==NULL) ||
	((KMMbin = (int*) malloc(int_nn_nbin))==NULL) ||
	((GAMbin = (int*) malloc(int_nn_nbin))==NULL) ||
	((MAXpos = (int*) malloc(sizeof(int) * *nmaxtests))==NULL) ||
	((X      = (double **) malloc(sizeof(double*) * n))==NULL) ||
	((DATA   = (double **) malloc(sizeof(double*) * n))==NULL) ||
	((DATA[0]= (double*) malloc(sizeof(double) * *npoints * nn))==NULL) 
	) { *error=MPPERR_MEMALLOC; goto ErrorHandling; }
  }
  if ((X[0] = coord)==0) {*error=MPPERR_COORD; goto ErrorHandling;}
  
  // first put the data into DATA, then the simulations ...    
  for (i=0; i<*npoints; i++) { DATA[0][i] = data[i]; }   
  endrepet = *npoints * *repet;
  for (i=0; i<endrepet; i++) DATA[0][*npoints+i] = simu[i];
 

  /* calculate values and teststatistics */
  if ((*error = 
       mcf_internal(E, ETest, Ebin, VAR, VARTest, SQ, SQTest, VARbin,
		    KMM, KMMbin, GAM, GAMbin, &p, bin, *nbin, *dim, n, 
		    col, nn, X, DATA, npoints, (bool) *staticchoice))!=0) {
    goto ErrorHandling;
  }
  
  // calculate rank of test statistic 
  k = NUMBERTESTS;
  for (i=0;i<NUMBERTESTS;i++){ Epos[i]=VARpos[i]=SQpos[i]=0; }
  for (i=1; i<nn; i++){//nn : number of total sets (true data + simulated ones)
    for (j=0; j<NUMBERTESTS; j++) { // j indexes the results for the true data
      if (ETest[j]>ETest[k]) {Epos[j]++;}
      if (VARTest[j]>VARTest[k]) {VARpos[j]++;}
      if (SQTest[j]>SQTest[k]) {SQpos[j]++;}
      k++;
      // k run from NUMBERTESTS to NUMBERTESTS * nn,
      // so through all comparison data sets
      }   
  }

  for (j=i=0; j<*nmaxtests; j++) {
    // unused in R -- work to be done
    MAXpos[j] = (Epos[Maxtests[i]] > SQpos[Maxtests[i+1]] ?
		 Epos[Maxtests[i]] : SQpos[Maxtests[i+1]]);
    i += 2;
  }
  
  if (MPP_PRINTLEVEL > 2) {
    PRINTF("  E");
    for (i=0;i<NUMBERTESTS;i++){ PRINTF(" %d",Epos[i]);}
    PRINTF("  V");
    for (i=0;i<NUMBERTESTS;i++){ PRINTF(" %d",VARpos[i]); }
    PRINTF(" SD");
    for (i=0;i<NUMBERTESTS;i++){ PRINTF(" %d",SQpos[i]); }
    PRINTF("\n");
  }
  
  if (*additive) {  
    // store ranks in matrix 
    for (i=0,k=0;i<NUMBERTESTS;i++,k+=*nestbin){
      store(Epos[i], &(Etestposition[k]), intestbin, *nestbin);
      store(VARpos[i], &(VARtestposition[k]), intestbin, *nestbin);
      store(SQpos[i], &(SQtestposition[k]), intestbin, *nestbin);
      }
    for (i=0,k=0;i<*nmaxtests;i++,k+=*nestbin)
      store(MAXpos[i], &(MAXtestposition[k]), intestbin, *nestbin);
  } else { // !*additive
    // store rank value 
    if (MPP_PRINTLEVEL>4) PRINTF("storing ranks directly...\n");
    for (i=0;i<NUMBERTESTS;i++){
      Etestposition[i]=Epos[i]; 
      VARtestposition[i]=VARpos[i]; 
      SQtestposition[i]=SQpos[i];	
    }
    for (i=0; i<*nmaxtests; i++) MAXtestposition[i]=MAXpos[i];
  }

  PutRNGstate();

  free(intestbin);
  free(E);    free(ETest);   free(Ebin);
  free(VAR);  free(VARTest);
  free(SQ);   free(SQTest);  free(VARbin);
  free(KMM);  free(KMMbin);
  free(GAM);  free(GAMbin);
  free(X);
  free(DATA[0]);   free(DATA); 
  free(MAXpos);
  return;
  
 ErrorHandling:   
  if (MPP_PRINTLEVEL>0) MPPErrorMessage(*error);
  if (intestbin!=NULL) free(intestbin);
  if (E!=NULL) free(E);    
  if (ETest!=NULL) free(ETest);   
  if (Ebin!=NULL) free(Ebin);
  if (VAR!=NULL) free(VAR);  
  if (VARTest!=NULL) free(VARTest);
  if (SQ!=NULL) free(SQ);   
  if (SQTest!=NULL) free(SQTest); 
  if (VARbin!=NULL) free(VARbin);
  if (KMM!=NULL) free(KMM);  
  if (KMMbin!=NULL) free(KMMbin);
  if (GAM!=NULL) free(GAM); 
  if (GAMbin!=NULL) free(GAMbin);
  if (X!=NULL) free(X);
  if (DATA[0]!=NULL) free(DATA[0]);  
  if (DATA!=NULL) free(DATA); 
  if (MAXpos!=NULL) free(MAXpos);
}



void mcf(double *E, /* size: rep * n^2 * col * #bins */
	 double *ETEST,/* size: rep * n^2*col * #Tests */
	 int *EBIN,/* size: rep *  n^2 * #bins */
	 double *VAR, /* rep * n^2 * c (c+1)/2 * #bins
		       (lower triangle including diagonal) */
	 double *VARTEST,
	 double *SQ,  /* rep * n^2 * c (c+1)/2 * #bins
		       (lower triangle including diagonal) */
	 double *SQTEST, 
	 int *VARBIN,
	 double *KMM,/* rep * [sum_i c_i] * [1+sum c_i]/2 * #bins */ 
	 int *KMMBIN, /* [sum_i c_i] * [1+sum c_i]/2 * #bins */
	 double *GAM,/* rep * [sum_i c_i] * [1+sum c_i]/2 * #bins */
	 int *GAMBIN,/* [sum_i c_i] * [1+sum c_i]/2 * #bins */
	 int *error, int *PrintLevel,
	 double *p,
	 double *bin, /*bin[0] MUST contain exactly the origin, bin=(-1,0,....)*/
	 int *lb, int *dim, int *n, int *col,
	 int *staticchoice, int *rep,
	 ...)
{ 
  double **X,   // coordinates of the individuals of a species
    **DATA;   // their marks 
  int *ROW;   // the number of individuals of a species
    //  ... = X[species1], DATA[species1], ROW[species1], X[species2], 
         
  long i;

  MPP_PRINTLEVEL = *PrintLevel;
  X = DATA = NULL;
  ROW = NULL; 

  va_list ap;
  va_start(ap,rep); 
  if ( ((X = (double **) malloc(sizeof(double*) * *n)) == NULL) ||
       ((DATA = (double **) malloc(sizeof(double*) * *n)) == NULL) ||
       ((ROW = (int *) malloc(sizeof(int) * *n)) == NULL) 
       ) { *error=MPPERR_MEMALLOC; goto ErrorHandling; }
  for (i=0; i<*n; i++) {
    X[i]    = va_arg(ap, double *);    
    DATA[i] = va_arg(ap, double *);    
    ROW[i]  = *(va_arg(ap,int *));
  }
  va_end(ap);

  if ((*error=mcf_internal(E, ETEST, EBIN, VAR, VARTEST, SQ, SQTEST, VARBIN, KMM,
			  KMMBIN, GAM, GAMBIN, p, bin, *lb, *dim, *n, *col, *rep,
			   X, DATA, ROW, (bool) *staticchoice)) !=0)
    { goto ErrorHandling;}
 
  free(ROW);
  free(X);
  free(DATA);
  *error = 0;
  return;
  
 ErrorHandling:
  if (MPP_PRINTLEVEL>0) MPPErrorMessage(*error);
  if (ROW!=NULL) free(ROW);
  if (DATA!=NULL) free(DATA);
  if (X!=NULL) free(X);
}



