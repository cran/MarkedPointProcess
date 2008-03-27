//#define DEBUG 1


/* 
 Authors
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Collection of auxiliary functions

 Copyright (C) 2001 -- 2006 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <math.h>
#include <unistd.h>
#include <assert.h>
#include "auxiliary.h"

/*
void xyInchToPixels(double *xy, double *Pixels) {
  DevDesc* dd;
  dd = CurrentDevice();
  Pixels[0] = GConvertXUnits(xy[0], INCHES, DEVICE, dd);
  Pixels[1] = GConvertYUnits(xy[1], INCHES, DEVICE, dd);
}
*/


void orderdouble(double *d, int *pos, int start, int end)
     /* quicksort algorithm, slightly modified:
        does not sort the data, but d[pos] will be ordered 
	NOTE: pos must have the values 0,1,2,...,start-end !
	(orderdouble is a kind of sorting pos according to
	the variable d)
     */
{
  int randpos, pivot, left, right, pivotpos, swap;
  double Dpivot;

  if( start < end )
  {   
    //GetRNGstate();randpos = start + (int) (UNIFORM_RANDOM * (end-start+1)); PutRNGstate();
    randpos = (int) (0.5 * (start + end));
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    Dpivot = d[pivot];
   
    pivotpos=start; 
    left = start;
    right=end+1;   

    while (left < right) {      
      while (++left < right && d[pos[left]] < Dpivot) pivotpos++;
      while (--right > left && d[pos[right]] > Dpivot) ;      
      if (left < right) {
        swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
        pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    orderdouble( d, pos, start, pivotpos-1);
    orderdouble( d, pos, pivotpos + 1, end);
  }
}


void quicksortdouble(double *d, int start, int end) {
  int i=start, j=end;
  double h, x=d[(start+end)/2];
  do {    
    while (d[i]<x) i++; 
    while (d[j]>x) j--;
    if (i<=j)
      {
	h=d[i]; d[i]=d[j]; d[j]=h;
	i++; j--;
      }
  } while (i<=j);
  if (start<j) quicksortdouble(d, start, j);
  if (i<end) quicksortdouble(d, i, end);
}


// the next function should be obsolete (soon)!
#define EPSILON_QUANTILE 0.00000001
double quantile(double *X, int lb, double p) 
{
  int *pos,i,j; 
  double result,*Y;
  assert( (pos = (int *) malloc(sizeof(int) * lb)) != NULL);
  assert( (Y = (double *) malloc(sizeof(double) * lb)) != NULL);
  for (j=i=0; i<lb; i++) { 
    if (!(RF_ISNA(X[i]))) {
      pos[j] = j; /* needed for orderdouble */  
      Y[j] = X[i];
      j++;
    }
  }
  if (j==0) {free(pos); free(Y); return RF_NAN;}
  orderdouble(Y, pos, 0, j-1);  
  
  result = Y[pos[(int) (p * (double) j + EPSILON_QUANTILE)]];
  free(pos);
  free(Y);
  return result;
} 
void Rquantile(double *X,int *lb,double *p, double *res) { *res = quantile(X,*lb,*p); }



void vectordist(double *v, int *Dim, double *dist, int *diag){
  int m, n, d, l, dim, r, lr, dr, add;
  add = (*diag==0) ? 1 : 0;
  l = Dim[0];
  dim = Dim[1] * l;
  lr = (l * (l + 1 - 2 * add)) / 2;
  for (r=0, m=0; m<l; m++) { // if add==1 loop is one to large, but
      // but doesn't matter
    for (n=m+add; n<l; n++, r++) {
      for (d=0, dr=0; d<dim; d+=l, dr+=lr) {
	  dist[r + dr] = v[n + d] - v[m + d];
      }
    }
  }
} 

static int ORDERDIM;
static double *ORDERD;
bool smaller(int i, int j)
{
  double *x, *y;
  int d;
  x = ORDERD + i * ORDERDIM;
  y = ORDERD + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
     if (x[d] != y[d]) return x[d] < y[d];
  return false;
}

bool greater(int i, int j)
{
  double *x, *y;
  int d;
  x = ORDERD + i * ORDERDIM;
  y = ORDERD + j * ORDERDIM;
  for(d=0; d<ORDERDIM; d++)
    if (x[d] != y[d]) return x[d] > y[d];
  return false;
}

void order(int *pos, int start, int end) {
  int randpos, pivot, left, right, pivotpos, swap;
  
  if( start < end ) {   
    //GetRNGstate();randpos = start + (int) (UNIFORM_RANDOM * (end-start+1)); PutRNGstate();
    randpos = (int) (0.5 * (start + end));
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    
    pivotpos=start; 
    left = start;
    right=end+1;   
    while (left < right) {      
      while (++left < right && smaller(pos[left], pivot)) pivotpos++;
      while (--right > left && greater(pos[right], pivot));      
      if (left < right) {
	swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
	pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    order(pos, start, pivotpos-1 );
    order(pos, pivotpos + 1, end );
  }
}

void ordering(double *d, int len, int dim, int *pos) 
{
  /* quicksort algorithm, slightly modified:
     does not sort the data, but d[pos] will be ordered 
     NOTE: pos must have the values 0,1,2,...,start-end !
     (orderdouble is a kind of sorting pos according to
     the variable d)
  */  
  int i;
  for (i=0; i<len;i++) pos[i]=i;
  ORDERD = d;
  ORDERDIM = dim;
  order(pos, 0, len-1);
}


void Ordering(double *d, int *len, int *dim, int *pos) 
{
  int i;
  for (i=0; i<*len; i++) pos[i]=i;
  ORDERD = d;
  ORDERDIM = *dim;
  order(pos, 0, *len-1 );
}

