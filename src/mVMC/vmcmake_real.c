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
 * make sample
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
// #include "vmcmake_real.h"
#ifndef _SRC_VMCMAKE_REAL
#define _SRC_VMCMAKE_REAL

#include "global.h"
#include "slater.h"
// #include "matrix.h"
// #include "pfupdate_real.h"
// #include "qp_real.h"
#include "splitloop.h"
#include "vmcmake.h"

int makeInitialSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                      const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  int flag = 1, flagRdc, loop = 0;
  int ri, mi, si, msi, rsi;
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  do {
/* initialize */
#pragma omp parallel for default(shared) private(msi)
    for (msi = 0; msi < nsize; msi++)
      eleIdx[msi] = -1;
#pragma omp parallel for default(shared) private(rsi)
    for (rsi = 0; rsi < nsite2; rsi++)
      eleCfg[rsi] = -1;

    /* local spin */
    for (ri = 0; ri < Nsite; ri++) {
      if (LocSpn[ri] == 1) {
        do {
          mi = gen_rand32() % Ne;
          si = (genrand_real2() < 0.5) ? 0 : 1;
        } while (eleIdx[mi + si * Ne] != -1);
        eleCfg[ri + si * Nsite] = mi;
        eleIdx[mi + si * Ne] = ri;
      }
    }

    /* itinerant electron */
    for (si = 0; si < 2; si++) {
      for (mi = 0; mi < Ne; mi++) {
        if (eleIdx[mi + si * Ne] == -1) {
          do {
            ri = gen_rand32() % Nsite;
          } while (eleCfg[ri + si * Nsite] != -1 || LocSpn[ri] == 1);
          eleCfg[ri + si * Nsite] = mi;
          eleIdx[mi + si * Ne] = ri;
        }
      }
    }

/* EleNum */
#pragma omp parallel for default(shared) private(rsi)
#pragma loop noalias
    for (rsi = 0; rsi < nsite2; rsi++) {
      eleNum[rsi] = (eleCfg[rsi] < 0) ? 0 : 1;
    }

    MakeProjCnt(eleProjCnt, eleNum);

    flag = CalculateMAll_fcmp(eleIdx, qpStart, qpEnd);
    // printf("DEBUG: maker4: PfM=%lf\n",creal(PfM[0]));
    if (size > 1) {
      MPI_Allreduce(&flag, &flagRdc, 1, MPI_INT, MPI_MAX, comm);
      flag = flagRdc;
    }

    loop++;
    if (loop > 100) {
      if (rank == 0)
        fprintf(stderr, "error: makeInitialSample: Too many loops\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  } while (flag > 0);

  return 0;
}

void copyFromBurnSample(int *eleIdx, int *eleCfg, int *eleNum,
                        int *eleProjCnt) {
  int i, n;
  const int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2 * Nsite + 2 * Nsite + NProj;
#pragma loop noalias
  for (i = 0; i < n; i++)
    eleIdx[i] = burnEleIdx[i];
  return;
}

void copyToBurnSample(const int *eleIdx, const int *eleCfg, const int *eleNum,
                      const int *eleProjCnt) {
  int i, n;
  int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2 * Nsite + 2 * Nsite + NProj;
#pragma loop noalias
  for (i = 0; i < n; i++)
    burnEleIdx[i] = eleIdx[i];
  return;
}

void saveEleConfig(const int sample, const double complex logIp,
                   const int *eleIdx, const int *eleCfg, const int *eleNum,
                   const int *eleProjCnt) {
  int i, offset;
  double x;
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  const int nProj = NProj;

  offset = sample * nsize;
#pragma loop noalias
  for (i = 0; i < nsize; i++)
    EleIdx[offset + i] = eleIdx[i];
  offset = sample * nsite2;
#pragma loop noalias
  for (i = 0; i < nsite2; i++)
    EleCfg[offset + i] = eleCfg[i];
#pragma loop noalias
  for (i = 0; i < nsite2; i++)
    EleNum[offset + i] = eleNum[i];
  offset = sample * nProj;
#pragma loop noalias
  for (i = 0; i < nProj; i++)
    EleProjCnt[offset + i] = eleProjCnt[i];

  x = LogProjVal(eleProjCnt);
  logSqPfFullSlater[sample] = 2.0 * (x + creal(logIp)); // TBC

  return;
}

void ReduceCounter(MPI_Comm comm) {
#ifdef _mpi_use
  int n = Counter_max;
  int recv[n];
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  MPI_Allreduce(Counter, recv, n, MPI_INT, MPI_SUM, comm);
  if (rank == 0) {
    int i;
    for (i = 0; i < n; i++)
      Counter[i] = recv[i];
  }
#endif
  return;
}

/* The mi-th electron with spin s hops to site rj */
void makeCandidate_hopping(int *mi_, int *ri_, int *rj_, int *s_,
                           int *rejectFlag_, const int *eleIdx,
                           const int *eleCfg) {
  const int icnt_max = Nsite * Nsite;
  int icnt;
  int mi, ri, rj, s, flag;

  flag = 0; // FALSE
  do {
    mi = gen_rand32() % Ne;
    s = (genrand_real2() < 0.5) ? 0 : 1;
    ri = eleIdx[mi + s * Ne];
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32() % Nsite;
    if (icnt > icnt_max) {
      flag = 1; // TRUE
      break;
    }
    icnt += 1;
  } while (eleCfg[rj + s * Nsite] != -1 || LocSpn[rj] == 1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_ = s;
  *rejectFlag_ = flag;

  return;
}

/* The mi-th electron with spin s exchanges with the electron on site rj with
 * spin 1-s */
void makeCandidate_exchange(int *mi_, int *ri_, int *rj_, int *s_,
                            int *rejectFlag_, const int *eleIdx,
                            const int *eleCfg, const int *eleNum) {
  int mi, mj, ri, rj, s, t, flag;

  flag = 1; // TRUE
  for (ri = 0; ri < Nsite; ri++) {
    if ((eleNum[ri] + eleNum[ri + Nsite]) == 1) {
      flag = 0; // FALSE
      break;
    }
  }
  if (flag) {
    *rejectFlag_ = flag;
    return;
  }

  do {
    mi = gen_rand32() % Ne;
    s = (genrand_real2() < 0.5) ? 0 : 1;
    ri = eleIdx[mi + s * Ne];
  } while (eleCfg[ri + (1 - s) * Nsite] != -1);
  do {
    mj = gen_rand32() % Ne;
    t = 1 - s;
    rj = eleIdx[mj + t * Ne];
  } while (eleCfg[rj + (1 - t) * Nsite] != -1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_ = s;
  *rejectFlag_ = flag;

  return;
}

/* The mi-th electron with spin s hops to site rj */
void updateEleConfig(int mi, int ri, int rj, int s, int *eleIdx, int *eleCfg,
                     int *eleNum) {
  eleIdx[mi + s * Ne] = rj;
  eleCfg[ri + s * Nsite] = -1;
  eleCfg[rj + s * Nsite] = mi;
  eleNum[ri + s * Nsite] = 0;
  eleNum[rj + s * Nsite] = 1;
  return;
}

void revertEleConfig(int mi, int ri, int rj, int s, int *eleIdx, int *eleCfg,
                     int *eleNum) {
  eleIdx[mi + s * Ne] = ri;
  eleCfg[ri + s * Nsite] = mi;
  eleCfg[rj + s * Nsite] = -1;
  eleNum[ri + s * Nsite] = 1;
  eleNum[rj + s * Nsite] = 0;
  return;
}

UpdateType getUpdateType(int path) {
  if (path == 0) {
    return HOPPING;
  } else if (path == 1) {
    return (genrand_real2() < 0.5) ? EXCHANGE
                                   : HOPPING; /* exchange or hopping */
  } else if (path == 2) {
    if (iFlgOrbitalGeneral == 0) {
      return EXCHANGE;
    } else {
      if (TwoSz == -1) { // Sz is not conserved
        return (genrand_real2() < 0.5) ? EXCHANGE : LOCALSPINFLIP;
            /* exchange or localspinflip */ // fsz
      } else {
        return EXCHANGE; /* exchange */
      }
    }
  } else if (path == 3) {        // for KondoGC
    if (genrand_real2() < 0.5) { // for conduction electrons
      return HOPPING;            /* hopping */
    } else { /* Exchange for conductions and local spins, localspinflip for
                local spins　*/
      return (genrand_real2() < 0.5)
                 ? EXCHANGE
                 : LOCALSPINFLIP; /* exchange or localspinflip */
    }
  }
  return NONE;
}

void VMCMakeSample_real(MPI_Comm comm) {
  int outStep, nOutStep;
  int inStep, nInStep;
  UpdateType updateType;
  int mi, mj, ri, rj, s, t, i;
  int ii, jj, si;
  int nAccept = 0;
  int sample;
  double rnd, rnd_p, w_save;

  double logIpOld, logIpNew;
      /* logarithm of inner product <phi|L|x> */ // is this ok ? TBC
  int projCntNew[NProj];
  double pfMNew_real[NQPFull];
  double x, w; // TBC x will be complex number

  int qpStart, qpEnd;
  int rejectFlag;
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  SplitLoop(&qpStart, &qpEnd, NQPFull, rank, size);

  StartTimer(30);
  if (BurnFlag == 0) {
    makeInitialSample(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, qpStart,
                      qpEnd, comm);
  } else {
    copyFromBurnSample(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt);
  }

  CalculateMAll_real(TmpEleIdx, qpStart, qpEnd);
  // printf("DEBUG: maker1: PfM=%lf\n",creal(PfM[0]));
  logIpOld = CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm);
  if (!isfinite(logIpOld)) {
    if (rank == 0)
      fprintf(stderr, "waring: VMCMakeSample remakeSample logIpOld=%e\n",
              creal(logIpOld)); // TBC
    makeInitialSample(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, qpStart,
                      qpEnd, comm);
    CalculateMAll_real(TmpEleIdx, qpStart, qpEnd);
    // printf("DEBUG: maker2: PfM=%lf\n",creal(PfM[0]));
    logIpOld = CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm);
    BurnFlag = 0;
  }
  StopTimer(30);

  nOutStep = (BurnFlag == 0) ? NVMCWarmUp + NVMCSample : NVMCSample + 1;
  nInStep = NVMCInterval * Nsite;

  for (i = 0; i < Counter_max; i++)
    Counter[i] = 0; /* reset counter */

  for (outStep = 0; outStep < nOutStep; outStep++) {
    for (inStep = 0; inStep < nInStep; inStep++) {

      updateType = getUpdateType(NExUpdatePath);

      if (updateType == HOPPING) { /* hopping */
        Counter[0]++;

        StartTimer(31);
        makeCandidate_hopping(&mi, &ri, &rj, &s, &rejectFlag, TmpEleIdx,
                              TmpEleCfg);
        StopTimer(31);

        if (rejectFlag)
          continue;

        StartTimer(32);
        StartTimer(60);
        /* The mi-th electron with spin s hops to site rj */
        updateEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        UpdateProjCnt(ri, rj, s, projCntNew, TmpEleProjCnt, TmpEleNum);
        StopTimer(60);
        StartTimer(61);
        // CalculateNewPfM2(mi,s,pfMNew,TmpEleIdx,qpStart,qpEnd);
        CalculateNewPfM2_real(mi, s, pfMNew_real, TmpEleIdx, qpStart, qpEnd);
        // printf("DEBUG: out %d in %d pfMNew=%lf
        // \n",outStep,inStep,creal(pfMNew[0]));
        StopTimer(61);

        StartTimer(62);
        /* calculate inner product <phi|L|x> */
        // logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        logIpNew = CalculateLogIP_real(pfMNew_real, qpStart, qpEnd, comm);
        StopTimer(62);

        /* Metroplis */
        x = LogProjRatio(projCntNew, TmpEleProjCnt);
        w = exp(2.0 * (x + (logIpNew - logIpOld)));
        if (!isfinite(w))
          w = -1.0; /* should be rejected */

        rnd = genrand_real2();

        /*
        printf("%d %d %8.7f %8.7f %8.7f %8.7f %8.7f", outStep, inStep, rnd, w,
        x, logIpNew, logIpOld); printf("|"); for(si=0;si<2;si++){
          for(jj=0;jj<Nsite;jj++){
            printf("%d",TmpEleNum[jj+si*Nsite]);
          }
          if(si==0) printf(",");
        }
        printf(">  \n");
        */
        // if (w > genrand_real2()) { /* accept */
        if (w > rnd) {
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          StartTimer(63);
          /*
          printf("%d  %8.7f %8.7f ", outStep, rnd, w);
          printf("|");
          for(si=0;si<2;si++){
            for(jj=0;jj<Nsite;jj++){
              printf("%d",TmpEleNum[jj+si*Nsite]);
            }
            if(si==0) printf(",");
          }
          printf(">  \n");
          */
          rnd_p = rnd;
          w_save = w;
          UpdateMAll_real(mi, s, TmpEleIdx, qpStart, qpEnd);
          //            UpdateMAll(mi,s,TmpEleIdx,qpStart,qpEnd);
          StopTimer(63);

          for (i = 0; i < NProj; i++)
            TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[1]++;
        } else { /* reject */
          revertEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        }
        StopTimer(32);

      } else if (updateType == EXCHANGE) { /* exchange */
        Counter[2]++;

        StartTimer(31);
        makeCandidate_exchange(&mi, &ri, &rj, &s, &rejectFlag, TmpEleIdx,
                               TmpEleCfg, TmpEleNum);
        StopTimer(31);

        if (rejectFlag)
          continue;

        StartTimer(33);
        StartTimer(65);

        /* The mi-th electron with spin s exchanges with the electron on site rj
         * with spin 1-s */
        t = 1 - s;
        mj = TmpEleCfg[rj + t * Nsite];

        /* The mi-th electron with spin s hops to rj */
        updateEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        UpdateProjCnt(ri, rj, s, projCntNew, TmpEleProjCnt, TmpEleNum);
        /* The mj-th electron with spin t hops to ri */
        updateEleConfig(mj, rj, ri, t, TmpEleIdx, TmpEleCfg, TmpEleNum);
        UpdateProjCnt(rj, ri, t, projCntNew, projCntNew, TmpEleNum);

        StopTimer(65);
        StartTimer(66);

        CalculateNewPfMTwo2_real(mi, s, mj, t, pfMNew_real, TmpEleIdx, qpStart,
                                 qpEnd);
        StopTimer(66);
        StartTimer(67);

        /* calculate inner product <phi|L|x> */
        logIpNew = CalculateLogIP_real(pfMNew_real, qpStart, qpEnd, comm);

        StopTimer(67);

        /* Metroplis */
        x = LogProjRatio(projCntNew, TmpEleProjCnt);
        w = exp(2.0 * (x + (logIpNew - logIpOld))); // TBC
        if (!isfinite(w))
          w = -1.0; /* should be rejected */

        rnd = genrand_real2();

        /*
        printf("%d  %d %8.7f %8.7f %8.7f %8.7f %8.7f", outStep, inStep, rnd, w,
        x, logIpNew, logIpOld); printf("|"); for(si=0;si<2;si++){
          for(jj=0;jj<Nsite;jj++){
            printf("%d",TmpEleNum[jj+si*Nsite]);
          }
          if(si==0) printf(",");
        }
        printf(">  \n");
        */

        // if (w > genrand_real2()) { /* accept */
        if (w > rnd) { /* accept */
          /*
          printf("%d %8.7f %8.7f ", outStep, rnd, w);
          printf("|");
          for(si=0;si<2;si++){
            for(jj=0;jj<Nsite;jj++){
              printf("%d",TmpEleNum[jj+si*Nsite]);
            }
            if(si==0) printf(",");
          }
          printf(">  \n");
          */
          rnd_p = rnd;
          w_save = w;
          StartTimer(68);
          UpdateMAllTwo_real(mi, s, mj, t, ri, rj, TmpEleIdx, qpStart, qpEnd);
          StopTimer(68);

          for (i = 0; i < NProj; i++)
            TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[3]++;
        } else { /* reject */
          revertEleConfig(mj, rj, ri, t, TmpEleIdx, TmpEleCfg, TmpEleNum);
          revertEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        }
        StopTimer(33);
      }

      if (nAccept > Nsite) {
        StartTimer(34);
        /* recal PfM and InvM */
        CalculateMAll_real(TmpEleIdx, qpStart, qpEnd);
        // printf("DEBUG: maker3: PfM=%lf\n",creal(PfM[0]));
        logIpOld = CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm);
        StopTimer(34);
        nAccept = 0;
      }
    } /* end of instep */

    StartTimer(35);
    /* save Electron Configuration */
    if (outStep >= nOutStep - NVMCSample) {
      sample = outStep - (nOutStep - NVMCSample);
      /*
      printf("%d %8.7f %8.7f ", sample, rnd_p, w_save);
      printf("|");
      for(si=0;si<2;si++){
        for(jj=0;jj<Nsite;jj++){
          printf("%d",TmpEleNum[jj+si*Nsite]);
        }
        if(si==0) printf(",");
      }
      printf(">  \n");
      */
      saveEleConfig(sample, logIpOld, TmpEleIdx, TmpEleCfg, TmpEleNum,
                    TmpEleProjCnt);
    }
    StopTimer(35);

  } /* end of outstep */

  copyToBurnSample(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt);
  BurnFlag = 1;
  return;
}

void MakeExactSample(MPI_Comm comm) {
  // int mi,mj,ri,rj,t,i,j;
  int sample, i, j;
  int icount = 0, jcount = 0;
  int isample = 0;
  int iel, jel;
  int NAllCfg = pow(2, Nsite);
  int NExactSample = pow(2, Nsite2);
  int EleExactCfg[2 * Nsite * NExactSample]; // the number of samples //
  int EleExactIdx[2 * Ne * NExactSample];    // the number of samples //
  int EleExactNum[2 * Nsite * NExactSample]; // the number of samples //

  // double complex logIpOld=0.0; // logarithm of inner product <phi|L|x> //
  int TempNumUp[Nsite], TempNumDn[Nsite];
  // int
  // TempEleCfg[Nsite],TempEleIdx[Ne],TempEleNum[Nsite],TempEleProjCnt[Nsite];
  // double x,w;
  // char filename[20];

  int icfg, jcfg, sampleStart, sampleEnd;
  int rank, size, tmp;

  int *moto = TmpEleIdx;

  // FILE *cfp;
  // int accepttime;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  SplitLoop(&sampleStart, &sampleEnd, NExactSample, rank, size);
  for (icfg = 0; icfg < NAllCfg; icfg++) {
    icount = 0;
    tmp = icfg;
    for (i = 0; i < Nsite; i++) {
      TempNumUp[i] = tmp % 2;
      if (TempNumUp[i] == 1) {
        icount += 1;
      }
      tmp = tmp / 2;
    }
    for (jcfg = 0; jcfg < NAllCfg; jcfg++) {
      jcount = 0;
      tmp = jcfg;
      for (j = 0; j < Nsite; j++) {
        TempNumDn[j] = tmp % 2;
        if (TempNumDn[j] == 1) {
          jcount += 1;
        }
        tmp = tmp / 2;
      }
      // printf("icfg=%d,icount=%d\n",icfg,icount);
      if (icount == Ne && jcount == Ne) {
        iel = 0;
        jel = 0;
        // printf("isample=%d,AllSample=%d\n",isample,sampleEnd-sampleStart);
        for (i = 0; i < Nsite; i++) {
          EleExactNum[isample * 2 * Nsite + i] = TempNumUp[i];
          EleExactCfg[isample * 2 * Nsite + i] = -1;
          EleExactNum[isample * 2 * Nsite + i + Nsite] = TempNumDn[i];
          EleExactCfg[isample * 2 * Nsite + i + Nsite] = -1;
          if (TempNumUp[i] == 1) {
            EleExactIdx[isample * 2 * Ne + iel] = i;
            EleExactCfg[isample * 2 * Nsite + i] = iel;
            iel++;
          }
          if (TempNumDn[i] == 1) {
            EleExactIdx[isample * 2 * Ne + jel + Ne] = i;
            EleExactCfg[isample * 2 * Nsite + i + Nsite] = jel;
            jel++;
          }
        }
        isample += 1;
      }
    }
  }

  NVMCSample = NExactSample = isample;
  SplitLoop(&sampleStart, &sampleEnd, NExactSample, rank, size);
  for (sample = sampleStart; sample < sampleEnd; sample++) {
    // printf("sample=%d\n",sample);
    //  save Electron Configuration //
    TmpEleIdx = EleExactIdx + sample * 2 * Ne;
    // printf("Tmp=%d,moto=%d\n",TmpEleIdx[0],EleExactIdx[sample*Ne]);
    TmpEleCfg = EleExactCfg + sample * 2 * Nsite;
    // printf("Tmp=%d,moto=%d\n",TmpEleCfg[0],EleExactCfg[sample*Nsite]);
    TmpEleNum = EleExactNum + sample * 2 * Nsite;
    // printf("Tmp=%d,moto=%d\n",TmpEleNum[0],EleExactNum[sample*Nsite]);
    MakeProjCnt(TmpEleProjCnt, TmpEleNum);
    // printf("Tmp=%d\n",TmpEleProjCnt[0]);
    // saveEleConfig(sample,logIpOld,0.0,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
    saveEleConfig(sample, 0.0, TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt);

    // printf("%d  ",sample);
    // for(i=0;i<Nsite;i++){
    //   printf("%d",TmpEleNum[i]);
    //   printf("%d",TmpEleNum[i+Nsite]);
    // }
    // printf("\n");
    /*int si,jj;
    printf("|");
    for(si=0;si<2;si++){
      for(jj=0;jj<Nsite;jj++){
        printf("%d",TmpEleNum[jj+si*Nsite]);
      }
      if(si==0) printf(",");
    }
    printf(">  \n");
    */
  } // end of outstep //
  TmpEleIdx = moto;

  return;
}

#endif
