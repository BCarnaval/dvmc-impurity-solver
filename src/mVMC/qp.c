/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Quantum Projection
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

//#include "qp.h"
#ifndef _QP_SRC
#define _QP_SRC
#include "global.h"
#include <complex.h>
//#include "legendrepoly.h"
#include "workspace.h"
#include "gauleg.c"


void UpdateQPWeight() {
  int i,j,offset;
  double complex tmp; //TBC

  if(FlagOptTrans>0) {
    for(i=0;i<NOptTrans;i++) {
      offset = i*NQPFix;
      tmp = OptTrans[i];
      for(j=0;j<NQPFix;j++) {
        QPFullWeight[offset+j] = tmp * QPFixWeight[j];
      }
    }
  } else {
    for(j=0;j<NQPFix;j++) {
      QPFullWeight[j] = QPFixWeight[j];
    }
  }

  return;
}


/* calculate SPGLCos, SPGLSin and QPFullWeight */
void InitQPWeight() {
  int i,j,idx;
  double *beta;
  double *weight;
  double w;

  RequestWorkSpaceDouble(2*NSPGaussLeg);
  beta   = GetWorkSpaceDouble(NSPGaussLeg);
  weight = GetWorkSpaceDouble(NSPGaussLeg);

  if(NSPGaussLeg==1) {
    beta[0] = 0.0;
    weight[0] = 1.0;
    SPGLCos[0] = 1.0;
    SPGLSin[0] = 0.0;
    SPGLCosSin[0] = 0.0;
    SPGLCosCos[0] = 1.0;
    SPGLSinSin[0] = 0.0;

    for(j=0;j<NMPTrans;j++) {
      //printf("XDEBUG: j=%d %lf %lf\n",j, creal(ParaQPTrans[j]), cimag(ParaQPTrans[j]));
      QPFixWeight[j] = ParaQPTrans[j];
    }

  } else {
    GaussLeg(0, M_PI, beta, weight, NSPGaussLeg);

    #pragma omp parallel for default(shared) private(i,j,w,idx)
    for(i=0;i<NSPGaussLeg;i++) {
      SPGLCos[i]    = cos(0.5 * beta[i])+0*I;    //TBC
      SPGLSin[i]    = sin(0.5 * beta[i])+0*I;    //TBC
      SPGLCosSin[i] = SPGLCos[i]*SPGLSin[i]+0*I; //TBC
      SPGLCosCos[i] = SPGLCos[i]*SPGLCos[i]+0*I; //TBC
      SPGLSinSin[i] = SPGLSin[i]*SPGLSin[i]+0*I; //TBC
      
      w = 0.5*sin(beta[i])*weight[i]*LegendrePoly(cos(beta[i]), NSPStot);
      
      for(j=0;j<NMPTrans;j++) {
        idx = i + j*NSPGaussLeg;
        QPFixWeight[idx] = w * ParaQPTrans[j];
      }
    }
  }

  UpdateQPWeight();

  ReleaseWorkSpaceDouble();
  return;
}




#endif
