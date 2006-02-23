#ifndef MarkedPointProcesses_H
#define MarkedPointProcesses_H 1

EXTERN void nearestneighbour(double *x, double *y, int *lx, double *px, 
			     double *py, int *lp, double *param, double *nnd); 

EXTERN void randomcoins(double *x, double *y, int *lx, double *px, double *py, 
			int *lp, double *param, double *res);
EXTERN void randomvariance(double *x, double *y, int *lx, double *px, 
			   double *py, int *lp, double *param, double *res);

EXTERN void GetNrMPPParameters(int *nr, int *n, int *np);
EXTERN void GetMPPModelName(int *nr,char **name);
EXTERN void GetMPPModelNr(char** name, int *n, int *nr);

EXTERN void GetmppParameters(int *numberLnorms, int *numberweights, 
			     int *numbertests, int *mppmaxchar,
			     int *modelnr);
EXTERN void SetmppParameters(int *action, int *BarnardBesag);

EXTERN void GenerateMPPData(double *coord, 
			    int *nrow, 
			    int *coordmethod, 
			    double *edgecorrection,
			    double *window,
			    double *lambda,
			    double *data,
			    int *ncol,
			    int *mppnr, int *nmpp, double *mppPList, int *nPList,
			    int *PrintLevel, int *error);

EXTERN void MCtest(int *repet, double *coord, double *data, int *npoints,
		   int *dim, double *simu, int *PrintLevel,
		   double *bin, int *nbin, 
		   int *Etestposition, int *VARtestposition, int *SQtestposition,
		   int *Maxtests, int *nmaxtests, int *MAXtestposition,
		   int *error, int *additive, int *copyANDstatic);

EXTERN void mcf(double *E, /* size: rep * n^2 * col * #bins */
	 double *ETEST,/* size: rep * n^2*col * #Tests */
	 int *EBIN,/* size: rep *  n^2 * #bins */
	 double *VAR,/* rep * n^2 * c (c+1)/2 * #bins
		      (lower triangle including diagonal) */
	 double *VARTEST,
	 double *SQ,  /* rep * n^2 * c (c+1)/2 * #bins 
		       (lower triangle including diagonal) */
	 double *SQTEST, 
	 int *VARBIN,
	 double *KMM,/* rep * [sum_i c_i] * [1+sum c_i]/2 * #bins */ 
	 int *KMMBIN, /* rep * [sum_i c_i] * [1+sum c_i]/2 * #bins */
	 double *GAM,/* rep * [sum_i c_i] * [1+sum c_i]/2 * #bins */
	 int *GAMBIN,/* rep * [sum_i c_i] * [1+sum c_i]/2 * #bins */
	 int *error, int *PrintLevel,
	 double *p,
	 double *bin, /* bin[0] MUST contain exactly the origin, bin=(-1,0,...)*/
	 int *lb, int *dim, int *n, int *col, int *staticchoice, int *rep, ...);

#endif /* MarkedPointProcesses_H */









