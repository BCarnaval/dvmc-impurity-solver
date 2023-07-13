/*
mVMC - A numerical solver package for a wide range of quantum lattice models
based on many-variable Variational Monte Carlo method Copyright (C) 2016 The
University of Tokyo, All rights reserved.

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
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#ifndef _SRC_VMCCAL
#define _SRC_VMCCAL

// #define _DEBUG_VMCCAL
// #define _DEBUG_VMCCAL_DETAIL

void clearPhysQuantity();

void calculateOO_real(double *srOptOO, double *srOptHO, const double *srOptO,
                      const double w, const double e, const int srOptSize);
void calculateOO_Store_real(double *srOptOO_real, double *srOptHO_real,
                            double *srOptO_real, const double w, const double e,
                            int srOptSize, int sampleSize);

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
  int val = (int)(percentage * 100);
  int lpad = (int)(percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\rComplete %3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

void VMCMainCal(MPI_Comm comm) {
  int *eleIdx, *eleCfg, *eleNum, *eleProjCnt;
  double complex e, ip;
  double w;
  double sqrtw;
  double complex we;
  double iprogress_d;

  const int qpStart = 0;
  const int qpEnd = NQPFull;
  int sample, sampleStart, sampleEnd, sampleSize;
  int i, info, tmp_i;

  /* optimazation for Kei */
  const int nProj = NProj;
  double complex *srOptO = SROptO;
  double *srOptO_real = SROptO_real;

  int rank, size, int_i;
  MPI_Comm_size(comm, &size);

  MPI_Comm_rank(comm, &rank);
#ifdef _DEBUG_VMCCAL
  printf("  Debug: SplitLoop\n");
#endif
  SplitLoop(&sampleStart, &sampleEnd, NVMCSample, rank, size);
  // printf(" %d / %d   %d %d \n",rank,size,sampleStart,sampleEnd);

  /* initialization */
  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);
  for (sample = sampleStart; sample < sampleEnd; sample++) {

    if (NVMCCalMode > 0 && PrintProgress > 0) {
      if (rank == 0 && size == 1) {
        iprogress_d = (1.0 * sample) / (sampleEnd - 1);
        printProgress(iprogress_d);
      }
    }
    // if(rank==0) printf("  Debug: sample=%d: CalculateMAll \n",sample);
#ifdef _DEBUG_VMCCAL
    int sample_to_print = 10000;
#endif

    eleIdx = EleIdx + sample * Nsize;
    eleCfg = EleCfg + sample * Nsite2;
    eleNum = EleNum + sample * Nsite2;
    eleProjCnt = EleProjCnt + sample * NProj;

    StartTimer(40);
#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: CalculateMAll \n", sample);
#endif
    assert(AllComplexFlag == 0);
    info = CalculateMAll_real(eleIdx, qpStart,
                              qpEnd); // InvM_real,PfM_real will change
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0; tmp_i < NQPFull * (Nsize * Nsize + 1); tmp_i++)
      InvM[tmp_i] = InvM_real[tmp_i]; // InvM will be used in SlaterElmDiff_fcmp
    StopTimer(40);

    if (info != 0) {
      fprintf(stderr,
              "warning: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n",
              rank, sample, info);
      continue;
    }
#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: CalculateIP \n", sample);
#endif
    assert(AllComplexFlag == 0);
    ip = CalculateIP_real(PfM_real, qpStart, qpEnd, MPI_COMM_SELF);

#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: LogProjVal \n", sample);
#endif
    /* calculate reweight */
    // w = exp(2.0*(log(fabs(ip))+x) - logSqPfFullSlater[sample]);
    w = 1.0;
#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: isfinite \n", sample);
#endif
    if (!isfinite(w)) {
      fprintf(stderr, "warning: VMCMainCal rank:%d sample:%d w=%e\n", rank,
              sample, w);
      continue;
    }

    // double u,t;

    StartTimer(41);
    /* calculate energy */
#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: calculateHam \n", sample);
#endif
    assert(AllComplexFlag == 0);
#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: calculateHam_real \n", sample);
#endif
    // e = CalculateHamiltonian_real(creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt);
    e = CalculateU_T_real(sample, creal(ip), eleIdx, eleCfg, eleNum,
                          eleProjCnt); //,&u,&t);
    // printf("  %4.3f %4.3f  %4.3f\n\n", e,test_u,test_t);
    // printf("MDEBUG: %lf %lf \n",creal(e),cimag(e));
    StopTimer(41);

#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: e = %lf %lf \n", sample, creal(e), cimag(e));
#endif
    if (!isfinite(creal(e) + cimag(e))) {
      fprintf(stderr, "warning: VMCMainCal rank:%d sample:%d e=%e\n", rank,
              sample, creal(e)); // TBC
      continue;
    }

    Wc += w;
    Etot += w * e;
    Etot2 += w * conj(e) * e;

    // print each energy measurement to file
    // if(sample>=sampleEnd-1000)
    // if(NVMCCalMode==0){

    /*
    if(PrintEnergy) {
      //int ii;
      fprintf(File_E,"%d %16.12f %16.12f\n",sample, w, e);
      if(NVMCCalMode>1){
        if(PrintConfig) {
          int jj,si;
          fprintf(File_Config,"|");
          for(si=0;si<2;si++){
            for(jj=0;jj<Nsite;jj++){
            //fprintf(File_Config,"%d %d %d ↑ %d ↓ \n",sample, jj, eleNum[jj],
    eleNum[jj+Nsite]); fprintf(File_Config,"%d",eleNum[jj+si*Nsite]);
            }
            if(si==0) fprintf(File_Config,",");
          }
          fprintf(File_Config,">  ");
          fprintf(File_Config,"%5d %16.12f   %4.3f %4.3f\n",sample, e, test_t,
    test_u);
        }
      }



      //if(fabs(creal(e)-creal(Etot/Wc))/(creal(Etot/Wc))>0.2) {
      //printf("Outlier energy measurement: %16.12f\n", e);
      //printf("Electron Configuration: \n");
      //for(ii=0;ii<Nsite;ii++){
      //  printf("%d %d ↑ %d ↓ \n", ii, eleIdx[ii], eleIdx[ii+Nsite]);
      //}
      //}
    }
    */
    //}

#ifdef _DEBUG_VMCCAL
    if (sample % sample_to_print == 0)
      printf("  Debug: sample=%d: calculateOpt \n", sample);
#endif
    if (NVMCCalMode == 0) {
      /* Calculate O for correlation fauctors */
      srOptO[0] = 1.0 + 0.0 * I; //   real
      srOptO[1] = 0.0 + 0.0 * I; //   real
#pragma loop noalias
      for (i = 0; i < nProj; i++) {
        srOptO[(i + 1) * 2] = (double)(eleProjCnt[i]); // even real
        srOptO[(i + 1) * 2 + 1] = 0.0 + 0.0 * I;       // odd  comp
      }

      StartTimer(42);
      /* SlaterElmDiff */
      SlaterElmDiff_fcmp(SROptO + 2 * NProj + 2, ip,
                         eleIdx); // TBC: using InvM not InvM_real
      StopTimer(42);

      //[s] this part will be used for real varaibles
      assert(AllComplexFlag == 0);
#pragma loop noalias
      for (i = 0; i < SROptSize; i++) {
        srOptO_real[i] = creal(srOptO[2 * i]);
      }
      //[e]

      StartTimer(43);
      /* Calculate OO and HO */
      if (NSRCG == 0 && NStoreO == 0) {
        assert(AllComplexFlag == 0);
        calculateOO_real(SROptOO_real, SROptHO_real, SROptO_real, w, creal(e),
                         SROptSize);
      } else {
        we = w * e;
        sqrtw = sqrt(w);
        assert(AllComplexFlag == 0);
#pragma omp parallel for default(shared) private(int_i)
        for (int_i = 0; int_i < SROptSize; int_i++) {
          // SROptO_Store for fortran
          SROptO_Store_real[int_i + sample * SROptSize] =
              sqrtw * SROptO_real[int_i];
          SROptHO_real[int_i] += creal(we) * SROptO_real[int_i];
        }
      }
      StopTimer(43);

    } else if (NVMCCalMode == 1) {
      fprintf(stdout, "OUPS: not implemented anymore.\n");
      exit(-1);
      //      StartTimer(42);
      //      /* Calculate Green Function */
      // #ifdef _DEBUG_VMCCAL
      //      fprintf(stdout, "Debug: Start: CalcGreenFunc\n");
      // #endif
      //      CalculateGreenFunc_real(w,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
      // #ifdef _DEBUG_VMCCAL
      //      fprintf(stdout, "Debug: End: CalcGreenFunc\n");
      // #endif
      //      StopTimer(42);

    } else if (NVMCCalMode == 2) {
      StartTimer(42);
#ifdef _DEBUG_VMCCAL
      if (sample % sample_to_print == 0)
        fprintf(stdout, "Debug: Start: CalculateStaticQuantities_real\n");
      fflush(stdout);
#endif

      int ii, jj;
      // printf("%16.12f %16.12f
      // %16.12f\n",creal(e),creal(Etot/Wc),fabs(creal(e)-creal(Etot/Wc))/(creal(Etot/Wc)));
      // if(fabs(creal(e)-creal(Etot/Wc))/(fabs(creal(Etot/Wc)))>0.2) {
      // if(){
      // printf("Outlier energy measurement: %16.12f\n", e);
      // printf("Electron Configuration: \n");
      for (ii = 0; ii < Nsite; ii++) {
        // printf("%d %d ↑ %d ↓ \n", ii, eleNum[ii], eleNum[ii+Nsite]);
        // if(eleNum[ii]+eleNum[ii+Nsite]>1){
        if (0) {
          printf("Outlier energy measurement: %16.12d\n", e);
          printf("Electron Configuration: \n");
          for (jj = 0; jj < Nsite; jj++) {
            printf("%d %d ↑ %d ↓ \n", jj, eleNum[jj], eleNum[jj + Nsite]);
          }
        }
      }

      CalculateStaticQuantities_real(w, eleIdx, eleCfg, eleNum, eleProjCnt);
#ifdef _DEBUG_VMCCAL
      if (sample % sample_to_print == 0)
        fprintf(stdout, "Debug: End: CalculateStaticQuantities_real\n");
      fflush(stdout);
#endif
      StopTimer(42);

    } else if (NVMCCalMode == 3) {
      StartTimer(42);
      // Calculate Dynamical Green Function
#ifdef _DEBUG_VMCCAL
      if (sample % sample_to_print == 0)
        fprintf(stdout, "Debug: Start: CalculateDynamicalGreenFunc_real\n");
      fflush(stdout);
#endif

      // printf("Calling CalculateDynamicalGreenFunc_real\n");
      CalculateDynamicalGreenFunc_real(creal(w), creal(ip), eleIdx, eleCfg,
                                       eleNum, eleProjCnt,
                                       sample % sampleChunk);

#ifdef _DEBUG_VMCCAL
      if (sample % sample_to_print == 0)
        fprintf(stdout, "Debug: End: CalculateDynamicalGreenFunc_real\n");
      fflush(stdout);
#endif
      StopTimer(42);
    }
  } // end of for(sample)
  // exit(0);

  // calculate OO and HO at NVMCCalMode==0
  if (NVMCCalMode == 0) {
    if (NSRCG != 0 || NStoreO != 0) {
      sampleSize = sampleEnd - sampleStart;
      assert(AllComplexFlag == 0);
      StartTimer(45);
      calculateOO_Store_real(SROptOO_real, SROptHO_real, SROptO_Store_real,
                             creal(w), creal(e), SROptSize, sampleSize);
      StopTimer(45);
    }
  }

  return;
}

void clearPhysQuantity() {
  int i, n;
  double complex *vec;
  double *vec_real;
  //[s] MERGE BY TM
  Wc = Etot = Etot2 = Sztot = Sztot2 = 0.0; // fsz
  // Wc = Etot = Etot2 = 0.0;
  Dbtot = Dbtot2 = 0.0;
  //[e] MERGE BY TM
  if (NVMCCalMode == 0) {
    /* SROptOO, SROptHO, SROptO */
    if (NSRCG != 0) {
      n = (2 * SROptSize) * 4; // TBC
    } else {
      n = (2 * SROptSize) * (2 * SROptSize + 2); // TBC
    }
    vec = SROptOO;
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < n; i++)
      vec[i] = 0.0 + 0.0 * I;
    // only for real variables
    if (NSRCG != 0) {
      n = (SROptSize)*4; // TBC
    } else {
      n = (SROptSize) * (SROptSize + 2); // TBC
    }
    vec_real = SROptOO_real;
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < n; i++)
      vec_real[i] = 0.0;
  } else if (NVMCCalMode == 1) {
    /* CisAjs, CisAjsCktAlt, CisAjsCktAltDC */
    n = NCisAjs + NCisAjsCktAlt + NCisAjsCktAltDC;
    vec = PhysCisAjs;
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < n; i++)
      vec[i] = 0.0 + 0.0 * I;

    assert(NLanczosMode <= 0);
  } else if (NVMCCalMode == 2) {
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < Nsite * ONE_SITE_PHYS_QTY; i++)
      PhysN1[i] = 0.0;
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < NCisAjs * TWO_SITES_PHYS_QTY; i++)
      PhysN2[i] = 0.0;
  } else if (NVMCCalMode == 3) {
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < Nsite * Nsite * NExcitation * NExcitation; i++) {
      Phys_nCHAm_averaged[i] = 0.0;
      Phys_nAHCm_averaged[i] = 0.0;
      Phys_nCAm_averaged[i] = 0.0;
      Phys_nACm_averaged[i] = 0.0;
    }
  }
  return;
}

void calculateOO_Store_real(double *srOptOO_real, double *srOptHO_real,
                            double *srOptO_Store_real, const double w,
                            const double e, int srOptSize, int sampleSize) {

  int i, j;
  char jobz, uplo;
  double alpha, beta, o;

  alpha = 1.0;
  beta = 0.0;

  jobz = 'N';
  uplo = 'T';
  if (NSRCG == 0) {
    M_DGEMM(&jobz, &uplo, &srOptSize, &srOptSize, &sampleSize, &alpha,
            srOptO_Store_real, &srOptSize, srOptO_Store_real, &srOptSize, &beta,
            srOptOO_real, &srOptSize);
  } else {
#pragma omp parallel for default(shared) private(i)
#pragma loop noalias
    for (i = 0; i < srOptSize; ++i) {
      srOptOO_real[i] = 0.0;
      srOptOO_real[i + srOptSize] = 0.0;
    }
    for (j = 0; j < sampleSize; ++j) {
#pragma omp parallel for default(shared) private(i, o)
#pragma loop noalias
      for (i = 0; i < srOptSize; ++i) {
        o = srOptO_Store_real[i + j * srOptSize];
        srOptOO_real[i] += o;
        srOptOO_real[i + srOptSize] += o * o;
      }
    }
  }

  return;
}

void calculateOO_real(double *srOptOO, double *srOptHO, const double *srOptO,
                      const double w, const double e, const int srOptSize) {
  double we = w * e;
  int m, n, incx, incy, lda;
  m = n = lda = srOptSize;
  incx = incy = 1;

  /* OO[i][j] += w*O[i]*O[j] */
  M_DGER(&m, &n, &w, srOptO, &incx, srOptO, &incy, srOptOO, &lda);

  /* HO[i] += w*e*O[i] */
  M_DAXPY(&n, &we, srOptO, &incx, srOptHO, &incy);

  return;
}

#endif
