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
 * Stochastic Reconfiguration method by PDPOSV
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

//#include "stcopt_pdposv.h"
#ifndef _SRC_STCOPT_PDPOSV
#define _SRC_STCOPT_PDPOSV

/* calculate the parameter change r[nSmat] from SOpt.
   The result is gathered in rank 0. */
int stcOptMain(double *r, const int nSmat, const int *smatToParaIdx, MPI_Comm comm) {
  /* global vector */
  double *w; /* workspace */
  /* distributed matrix (row x col) */
  double *s; /* overlap matrix (nSmat x nSmat) */
  double *g; /* energy gradient and parameter change (nSmat x 1) */

  int sSize,gSize;
  int wSize=nSmat;

  /* for MPI */
  int rank,size;
  int dims[2]={0,0};
  MPI_Comm comm_col;

  /* index table */
  int *irToSmatIdx; /* table for local indx to Smat index */
  int *irToParaIdx, *icToParaIdx; /* tables for local index to Para index */

  /* for BLACS */
  int ictxt,nprow,npcol,myprow,mypcol;
  char procOrder='R';

  /* for array descriptor of SCALAPACK */
  int m,n;
  int mlld,mlocr,mlocc; /* for matrix */
  int vlld,vlocr,vlocc; /* for vector */
  int mb=64, nb=64; /* blocking factor */
  int irsrc=0, icsrc=0;
  int descs[9],descg[9];

  /* for PDSYEVD */
  char uplo;
  int nrhs, is, js, ig, jg;

  int info;
  const double ratioDiag = 1.0+DSROptStaDel;
  int si,pi,pj,idx;
  int ir,ic;

  const int srOptSize = SROptSize;//TBC
  const double dSROptStepDt = DSROptStepDt;
  const double srOptHO_0 = creal(SROptHO[0]);
 // const double complex srOptHO_0 = SROptHO[0];
  double complex *srOptOO=SROptOO;
  double complex *srOptHO=SROptHO;

  StartTimer(55);

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  MPI_Dims_create(size,2,dims);

  /* initialize the BLACS context */
  nprow=dims[0]; npcol=dims[1];
  ictxt = Csys2blacs_handle(comm);
  Cblacs_gridinit(&ictxt, &procOrder, nprow, npcol);
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myprow, &mypcol);

  /* initialize array descriptors for distributed matrix s */
  m=n=nSmat; /* matrix size */
  mlocr = M_NUMROC(&m, &mb, &myprow, &irsrc, &nprow);
  mlocc = M_NUMROC(&n, &nb, &mypcol, &icsrc, &npcol);
  mlld = (mlocr>0) ? mlocr : 1;
  M_DESCINIT(descs, &m, &n, &mb, &nb, &irsrc, &icsrc, &ictxt, &mlld, &info);
  sSize = (mlocr*mlocc>0) ? mlocr*mlocc : 1;

  /* initialize array descriptors for distributed vector g */
  m=nSmat; n=1; /* vector */
  vlocr = mlocr; /* M_NUMROC(&m, &mb, &myprow, &irsrc, &nprow); */
  vlocc = M_NUMROC(&n, &nb, &mypcol, &icsrc, &npcol);
  vlld = (vlocr>0) ? vlocr : 1;
  M_DESCINIT(descg, &m, &n, &mb, &nb, &irsrc, &icsrc, &ictxt, &vlld, &info);
  gSize = (vlocr*vlocc>0) ? vlocr*vlocc : 1;

  //printf("DEBUG: n=%d nSmat=%d gSize=%d  sSize=%d wSize=%d: srOptSize=%d \n",n,nSmat,gSize,sSize,wSize,srOptSize);
  //printf("DEBUG: vlocc=%d mlocc=%d \n",vlocc,mlocc);
  /* allocate memory */
  RequestWorkSpaceDouble(sSize+gSize+wSize);
  s = GetWorkSpaceDouble(sSize);
  g = GetWorkSpaceDouble(gSize);
  w = GetWorkSpaceDouble(wSize);

  /* Para indices of the distributed vector and calculate them */
  RequestWorkSpaceInt(2*mlocr+mlocc);
  irToSmatIdx = GetWorkSpaceInt(mlocr);
  irToParaIdx = GetWorkSpaceInt(mlocr);
  icToParaIdx = GetWorkSpaceInt(mlocc);

  #pragma omp parallel for default(shared) private(ir,si)
  for(ir=0;ir<mlocr;ir++) {
    si = (ir/mb)*nprow*mb + myprow*mb + (ir%mb);
    irToSmatIdx[ir] = si;
    irToParaIdx[ir] = smatToParaIdx[si];
  }
  #pragma omp parallel for default(shared) private(ic)
  for(ic=0;ic<mlocc;ic++) {
    icToParaIdx[ic] = smatToParaIdx[(ic/nb)*npcol*nb + mypcol*nb + (ic%nb)];
  }

  StopTimer(55);
  StartTimer(56);
  /* calculate the overlap matrix S */
  //printf("YDEBUG: %d %d %lf \n",mlocc,mlocr,ratioDiag);
  #pragma omp parallel for default(shared) private(ic,ir,pi,pj,idx)
  #pragma loop noalias
  for(ic=0;ic<mlocc;ic++) {
    pj = icToParaIdx[ic]; /* Para index (global) */
    for(ir=0;ir<mlocr;ir++) {
      pi = irToParaIdx[ir]; /* Para index (global) */
      idx = ir + ic*mlocr; /* local index (row major) */

      /* S[i][j] = xOO[i+1][j+1] - xOO[0][i+1] * xOO[0][j+1]; */
      s[idx] = creal(srOptOO[(pi+2)*(2*srOptSize)+(pj+2)]) - creal(srOptOO[pi+2]) * creal(srOptOO[pj+2]);
      /* modify diagonal elements */
      //printf("DEBUG: idx=%d %d %d s[]=%lf \n",idx,pi,pj,s[idx]);
      //printf("XDEBUG %d %d %lf \n",ic,ir,s[idx]);
      if(pi==pj) s[idx] *= ratioDiag; // TBC
    }
  }

  /* calculate the energy gradient g and multiply (-dt) */
  if(vlocc>0) {
    #pragma omp parallel for default(shared) private(ir,pi)
    #pragma loop noalias
    #pragma loop norecurrence irToParaIdx
    for(ir=0;ir<vlocr;ir++) {
      pi = irToParaIdx[ir]; /* Para index (global) */
      
      /* energy gradient = 2.0*( xHO[i+1] - xHO[0] * xOO[0][i+1]) */
      /* g[i] = -dt * (energy gradient) */
      g[ir] = -dSROptStepDt*2.0*(creal(srOptHO[pi+2]) - srOptHO_0 * creal(srOptOO[pi+2]));
      //printf("ZDEBUG: %d %lf \n",ir,g[ir]);
    }
  }

  StopTimer(56);

#ifdef _DEBUG_STCOPT
  fprintf(stderr, "g:\n");
  for(ir=0; ir<vlocr; ++ir){
    fprintf(stderr, "%lg\n", g[ir]);
  }
  fprintf(stderr, "S:\n");
  for(ic=0;ic<mlocc;ic++) {
    for(ir=0;ir<mlocr;ir++) {
      idx = ir + ic*mlocr; /* local index (row major) */
      fprintf(stderr, "%lg ", s[idx]);
    }
    fprintf(stderr, "\n");
  }
#endif

  StartTimer(57);

  /***** solve the linear equation S*r=g by PDPOSV *****/
  uplo='U'; n=nSmat; nrhs=1; is=1; js=1; ig=1; jg=1;
  M_PDPOSV(&uplo, &n, &nrhs, s, &is, &js, descs, 
           g, &ig, &jg, descg, &info); /* g is overwritten with the solution. */
  /* error handle */
  if(info!=0) {
    if(rank==0) fprintf(stderr,"error: PDPOSV info=%d\n",info);
    return info;
  }
  /* end of diagonalization */

  StopTimer(57);
  StartTimer(58);

  /* create a communicator whose process has the same value of mypcol */
  MPI_Comm_split(comm,mypcol,myprow,&comm_col); 

  /* clear workspace */
  for(si=0;si<nSmat;si++) w[si] = 0.0;

  if(mypcol==0) {
    /* copy the solution to workspace */
    #pragma omp parallel for default(shared) private(ir,si)
    #pragma loop noalias
    #pragma loop norecurrence irToSmatIdx
    for(ir=0;ir<vlocr;ir++) {
      si = irToSmatIdx[ir]; /* Smat index (global) */
      w[si] = g[ir];
    }

    /* gather the solution to r on rank=0 process */
    SafeMpiReduce(w,r,nSmat,comm_col);
  }

  StopTimer(58);

#ifdef _DEBUG_STCOPT
  fprintf(stderr, "r:\n");
  for(si=0; si<nSmat; ++si){
    fprintf(stderr, "%lg\n", r[si]);
  }
#endif

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  MPI_Comm_free(&comm_col);
  Cblacs_gridexit(ictxt);
  return info;
}

#endif
