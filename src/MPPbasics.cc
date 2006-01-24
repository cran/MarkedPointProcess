/* 
 Authors
 Martin Schlather, schlath@hsu-hh.de

 library for simulation and analysis of marked point processes:
    basic input/output

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


#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "MPP.h"


int currentmppNr=-1;
mpp_type mpp_model[MPP_MODELS];


void MPPErrorMessage(int error) {
  char EM[50];
  if ((error<800) || (error>=900)) {
    if (MPP_PRINTLEVEL>8) PRINTF("error returned to calling function in MPP");
    return;
  }
  switch (error) {
  case MPPERR_NODATA  : strcpy(EM,"No mpp data to be simulated? Hmm"); break;
  case MPPERR_PARAMS  : strcpy(EM,"Number of mpp parameter not correct"); break;
  case MPPERR_POINTS  : strcpy(EM,"Not enough points"); break;  
  case MPPERR_MEMALLOC: strcpy(EM,"Memory allocation error"); break;
  case MPPERR_MODELNR : strcpy(EM,"specified number of model not allowed");
    break;
  case MPPERR_ : strcpy(EM,""); break;
  default : PRINTF("** %d\n", error); assert(false);
  } 
  PRINTF(" mpp error: %s.\n",EM);
}


void InitMPPModelList()
{

  if (currentmppNr==-1) {
    currentmppNr=0;
    assert(currentmppNr<MPP_MODELS); 

    strncpy(mpp_model[currentmppNr].name, "nearest neighbour", MPP_MAXCHAR-1);
    mpp_model[currentmppNr].nparam = NN_MAX; 
    mpp_model[currentmppNr].fct = nearestneighbour;
    assert(currentmppNr<MPP_MODELS); 
    currentmppNr++;

    strncpy(mpp_model[currentmppNr].name, "random coins", MPP_MAXCHAR-1);
    mpp_model[currentmppNr].nparam = RC_MAX; 
    mpp_model[currentmppNr].fct = randomcoins;
    assert(currentmppNr<MPP_MODELS); 
    currentmppNr++;

    strncpy(mpp_model[currentmppNr].name, "variance by coins", MPP_MAXCHAR-1);
    mpp_model[currentmppNr].nparam = RC_MAX; 
    mpp_model[currentmppNr].fct = randomvariance;
    assert(currentmppNr<MPP_MODELS); 
    currentmppNr++;

  } else assert(currentmppNr==3);
}

void GetNrMPPParameters(int *nr, int *n, int *np){
  int i;
  for (i=0; i<*n; i++) np[i]=mpp_model[nr[i]].nparam;
}

void GetMPPModelName(int *nr,char **name){
  if (currentmppNr==-1) InitMPPModelList();
  if ((*nr<0) ||(*nr>=currentmppNr)) {strncpy(*name,"",MPP_MAXCHAR); return;}
  strncpy(*name,mpp_model[*nr].name,MPP_MAXCHAR);
}


void GetMPPModelNr(char **name, int *n, int *nr) 
{
  unsigned int ln,v;
  unsigned int nn;
  // == -1 if no matching function is found
  // == -2 if multiple matching fnts are found, without one matching exactly
  // if more than one match exactly, the last one is taken (enables overwriting 
  // standard functions)
  if (currentmppNr==-1) InitMPPModelList();
  nn = (unsigned int) *n;

  for (v=0; v<nn; v++) {
    nr[v]=0;
    ln=strlen(name[v]);
    while ( (nr[v]<currentmppNr) && strncmp(name[v],mpp_model[nr[v]].name,ln)) {
      (nr[v])++;
    }
    if (nr[v]<currentmppNr) { 
      // a matching function is found. Are there other functions that match?
      int j; 
      bool exactmatching, multiplematching;
      exactmatching=(ln==strlen(mpp_model[nr[v]].name));
      multiplematching=false;
      j=nr[v]+1; // if two or more covariance functions have the same name 
      //            the last one is taken 
      while (j<currentmppNr) {
	while ( (j<currentmppNr) && strncmp(name[v],mpp_model[j].name,ln)) j++;
	if (j<currentmppNr) {
	  if (ln==strlen(mpp_model[j].name)) {nr[v]=j; exactmatching=true;} 
	  else {multiplematching=true;}
	}
	j++;
      } 
      if (!exactmatching && multiplematching) {nr[v]=-2;}
    } else nr[v]=-1;
  }
}

void GetmppParameters(int *numberLnorms, int *numberweights, 
		      int *numbertests, int *mppmaxchar, int *modelnr){
  // achtung!! auch in tests/CHECK.R verwendet!!
  if (currentmppNr==-1) InitMPPModelList();
  *numberLnorms = NUMBER_L_NORMS;
  *numberweights= NUMBERWEIGHTS;
  *numbertests  = NUMBERTESTS;
  *mppmaxchar   = MPP_MAXCHAR;
  *modelnr      = currentmppNr;  
}

