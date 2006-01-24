
/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for simulation and analysis of marked point processes:
    simulation of mpp

 Copyright (C) 2004 -- 2006  Martin Schlather, 

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

// Version 27.6.00
//#define DEBUG 1
#define CHECK 1

#include <math.h>
#include <assert.h>
#include "MPP.h"

#define PRECISION 1E-11

void nearestneighbour(double *x, double *y, int *lx, double *px, double *py,
		      int *lp, double *param, double *nnd) 
       /* x,y : coord. of the points where the nearest neighbour distances   
		are to be calculated (excluding the point itself)
	  lx  : number of points
	  px,py:coordinates of points of the point process
	  lp  : number of points of the point process
	  nnd : nearest neighbour distance; nnd is a vector of length lx
	  nnd == -1 if there is no nearest neighbour
       */   
{ 
  int k;
  {
    register int i, j;
    register double dist, mindist;
    for (i=0;  i<*lx; i++){
      mindist = RF_INF;
      for (j=0; j<*lp; j++){
	if (((dist=(dist=x[i]-px[j]) * dist+(dist=y[i]-py[j]) * dist)<mindist)
	    && (dist>0)) {mindist=dist;};
      };
      nnd[i]=mindist; 
    };
  }
  for (k=0;  k<*lx; k++) {nnd[k] = sqrt(nnd[k]) * param[NN_FACTOR];}
} 


/* the function additiveboolean calculates the convolution of 
   a fct (e.g. sphere, triangle, below) with the points of a point process
   at some given locations
   if fct is "sphere", see switch below, and if the locations are chosen to 
   be the points of the point process, then the marks of the process given 
   in Waelder and Stoyan, example 1 are obtained
*/

void randomcoins(double *x, double *y, int *lx, double *px, double *py, int *lp,
		     double *param, double *res) 
/* x,y   : coordinates of the locations of interest;
   lx    : # of locations of interest;
   px,py : coordinates of the points of the point process
   lp    : # of points of the point process
   param : fct   : chooses the boolean function, see below
           r     : is the radius of the sphere; 
   res   : the calculated marks for the point(s) (x,y)
*/
{
  int fct; double r;
  register double xx,yy,rsq;
  register int k,i; 

  fct = (int) (param[RC_NR] + 0.5);
  r = param[RC_RADIUS];

  for (k=0; k<*lx; res[k++]=0.0);
  switch (fct) {
  case 1 :  /* disk of height 1 and fixed radius r*/
    rsq = r * r; 
    for (k=0; k<*lx; k++){
      for(i=0;i<*lp;i++){
	if ((xx=x[k]-px[i]) * xx + (xx=y[k]-py[i])*xx<rsq) {res[k]+=1.0;}
      }
    }
    break;
  case 2 :  /* cone of height one and fixed radius r*/
    return;
    for (k=0; k<*lx; k++){
      for(i=0;i<*lp;i++){
	xx = 1.0-fabs(x[k]-px[i]) * r; 
	yy = (1.0-fabs(y[k]-py[i])) * r; 
	if ((xx>0.0) && (yy>0.0)) {res[k]+=xx * yy;}
      }
    }
    break;
  default: assert(false);
  }
  
  register double height;
  height = param[RC_HEIGHT];
  for (k=0; k<*lx; k++) { res[k] *= height; }; 
}


void randomvariance(double *x, double *y, int *lx, double *px, double *py, 
		    int *lp, double *param, double *res) 
{ 
  int k;
  randomcoins(x, y, lx, px, py, lp, param, res);
  for (k=0; k<*lx; k++) { res[k] = GAUSS_RANDOM(res[k]); } 
}



// coordparameter positions for method 1 and 2
#define XMIN     0
#define XMAX     1
#define YMIN     2
#define YMAX     3


/* up to now only 2-dimensional data is simulated */
void GenerateMPPData(double *coord, 
		     int *npoints, 
		     int *coordmodel,
		     double *edgecorrection,
		     double *window,
		     double *lambda,
		     double *data,
		     int *ncol,
		     int *mppnr, int *nmpp, double *mppPList, int *nPList,     
		     int *PrintLevel, int *error)
{
  /*    
    coord : matrix of *npoints x 2 entries (which might be empty ones)
    npoints  : is passed back giving the actrow !
    coordmodel :
        0 : fixed coordinates
        1 : npoints coordinates scattered uniforlmy in 
            [window[0],window[1]]  x  [window[2], window[3]]
        2 : actrow=min { Poisson[window[4]], *npoints } 
            npoints coordinates scattered uniforlmy in 
            [window[0],window[1]]  x  [window[2], window[3]],		    
     window : see #define above
     lambda: intensity for the Poisson processes
     data : ncol independent set of generated marks, 
            the data are stored in the first actrow rows of each of the ncol
            columns
     ncol : number of independent sets of marks to be generated given
            the locations
     mppnr: method numbers of the traditional methods to generate marks
            (0: nearest neighbour; 1:random coins
     nmpp : number of traditional methods
     mppPList : parameters for the methods, all in a sequence
     nPlist   : number of total parameters; control parameter;
                the required number of parameters is calculated and compared 
                with nPList	    
   */
  
  int actrow, np, curparam;
  long i,k, npointsncol, segment;
  double diffx,diffy,minx,maxx,miny,maxy;
  double *G, *px, *py;
  bool marksinitialised;

  MPP_PRINTLEVEL = *PrintLevel;
  px = py = G = NULL;
  npointsncol = *npoints * *ncol;

  if (*ncol<1) {*error=MPPERR_NODATA; goto ErrorHandling;} 
  
  GetRNGstate();

  if (currentmppNr==-1) InitMPPModelList();

  *error = 0;
   
  px =  coord;
  py =  &coord[*npoints];
  np =  *npoints;
  marksinitialised = true; /* future models will probably create locations 
			      and marks at once */

  switch (*coordmodel) {
  case 0 : actrow=*npoints; /* coordinated fixed and passed by coord */
    marksinitialised = false;
    break;
  case 1 : /* uniformly distributed */   
    actrow=*npoints; 
    minx=window[XMIN];
    miny=window[YMIN];
    diffx = window[XMAX]-window[XMIN];
    diffy = window[YMAX]-window[YMIN];

    for (i=0;i<actrow;i++) {
      coord[i]=UNIFORM_RANDOM * diffx + minx;
      coord[i+ *npoints] = UNIFORM_RANDOM * diffy + miny;
    }     
    marksinitialised = false;
   break;
  case 2: /* Poisson process with intensity lambda */
    // add edgecorrected points to a point list (assume Poisson distributed), 
    // see below at the moment this is not done very efficiently.
    minx = window[XMIN];
    miny = window[YMIN]; 
    diffx = window[XMAX] - window[XMIN];
    diffy = window[YMAX] - window[YMIN];
    for (;;) {
      actrow=(int) (INT_EPSILON + POISSON_RANDOM(*lambda * diffx * diffy));
      if (actrow<=*npoints) break;
      else // the chance be get there is extremely small
	if (MPP_PRINTLEVEL>1) 
	PRINTF("actrow (%d) > *npoints=%d, limited to *npoints \n",
	       actrow, *npoints);
    }
    for (i=0;i<actrow;i++) {
      coord[i] = UNIFORM_RANDOM * diffx + minx;
      coord[i+ *npoints] = UNIFORM_RANDOM * diffy + miny;
    }   
    marksinitialised = false;
    break;
  default : assert(false);
  }

  if (actrow<2) { *error=MPPERR_POINTS; goto ErrorHandling;} 
    
  if (!marksinitialised) 
    for (segment=0; segment<npointsncol; segment+=*npoints) {      
      k = segment+actrow;
      for (i=segment; i<k; i++) { data[i] = 0.0; }
    }

  // edge correction of points, recommended for traditional
  // methods of mark generation
  if (*edgecorrection>0.0) {
    assert(*coordmodel>=0 && *coordmodel<=2); 
    minx = window[XMIN] - *edgecorrection;
    miny = window[YMIN] - *edgecorrection;
    maxx = window[XMAX] + *edgecorrection;
    maxy = window[YMAX] + *edgecorrection;
    diffx = maxx - minx;
    diffy = maxy - miny;
    np = (int) (INT_EPSILON + POISSON_RANDOM(*lambda * diffx * diffy));
    if ( ((px = (double*) malloc(sizeof(double) * (np + actrow)))==NULL) ||
	 ((py = (double*) malloc(sizeof(double) * (np + actrow)))==NULL) ) {
      *error=MPPERR_MEMALLOC; goto ErrorHandling;}
    k = 0; 
    for (i=0; i<np; i++) { 
      px[k]= UNIFORM_RANDOM * diffx + minx;
      py[k]= UNIFORM_RANDOM * diffy + miny;
      if ((px[k]<window[XMIN]) // i.e., outside the sampling area...
	  ||(px[k]>window[XMAX]) ||(py[k]<window[YMIN]) || (py[k]>window[YMAX])){
	k++; 
      }
    }
    for (i=0; i<actrow; i++) { 
      px[i+k] = coord[i];
      py[i+k] = coord[i+*npoints];
    }
    np = actrow+k;
  }

  // generate the marks of the mpp models
  // first check if the number of passed parameters is correct
  for (curparam=i=0; i<*nmpp; i++) {
    if ((mppnr[i]<0) || (mppnr[i]>=currentmppNr)) 
      {*error=MPPERR_MODELNR; goto ErrorHandling;}
    curparam += mpp_model[mppnr[i]].nparam;
  }
  if (curparam != *nPList) { *error=MPPERR_PARAMS; goto ErrorHandling; }

  // generation of the marks
  curparam = 0;
  if ((G = (double *) malloc(sizeof(double) * actrow))==NULL) {
    *error = MPPERR_MEMALLOC; goto ErrorHandling;}
  for (i=0; i<*nmpp; i++) {
    mpp_model[mppnr[i]].fct(coord,&coord[*npoints],&actrow,px,py,&np, 
			    &(mppPList[curparam]),G);
    for (k=0;k<actrow;k++)
      for (segment=k; segment<npointsncol; segment+=*npoints) 
	data[segment] += G[k];
    curparam += mpp_model[mppnr[i]].nparam;
  }
  free(G);
  
  PutRNGstate();

  if (*edgecorrection) { 
    if (px!=NULL) free(px);
    if (py!=NULL) free(py); 
  }
  *npoints = actrow; // !! actrow is passed back
  return;
  
 ErrorHandling:
  if (G!=NULL) free(G);
  if (*edgecorrection) { 
    if (px!=NULL) free(px);
    if (py!=NULL) free(py); 
  }
  if (MPP_PRINTLEVEL>0) MPPErrorMessage(*error);
}

