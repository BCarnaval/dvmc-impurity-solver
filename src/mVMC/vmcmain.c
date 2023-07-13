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
 * main program
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
/* #include "fjcoll.h" */
#include "vmcmain.h"
// #include "physcal_lanczos.h"

// #define _DEBUG
// #define _DEBUG_DUMP_SROPTO_STORE
// #define _DEBUG_DUMP_SROPTOO
// #define _DEBUG_DUMP_PARA

int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1,
               MPI_Comm comm_child2);
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1,
               MPI_Comm comm_child2);
void outputData();
void StdFace_main(char *fname);

/*main program*/
int main(int argc, char *argv[]) {
  /* input file name */
  char fileDefList[256];
  char fileInitPara[256];

  int info = 0;

  int flagReadInitPara = 0;

  /* for Standard mode (-s option)*/
  // int flagStandard = 0;
  /* for getopt() */
  // int option;
  extern char *optarg;
  extern int optind, opterr, optopt;
  /* for strtol() */
  extern int errno;
  // char *endptr;
  // long num;

  /* for MPI */
  int rank0 = 0, size0 = 1;
  int group1 = 0; //,rank1,rank2,size1,size2;
  MPI_Comm comm0 = 0, comm1 = 0, comm2 = 0;

  MPI_Init(&argc, &argv);
  NThread = omp_get_max_threads();

  printf("NThread: %d\n", NThread);

  InitTimer();
  StartTimer(0);
  StartTimer(1);
  StartTimer(10);

  /* set input filename */
  strcpy(fileDefList, argv[optind]);
  if (argc - optind > 1) {
    flagReadInitPara = 1;
    strcpy(fileInitPara, argv[optind + 1]);
  }

  MPI_Comm_dup(MPI_COMM_WORLD, &comm0);

  MPI_Comm_rank(comm0, &rank0);
  MPI_Comm_size(comm0, &size0);
  StopTimer(10);
  /*
   Standard mode: generating input files
  */
  MPI_Barrier(comm0);

  StartTimer(11);
  if (rank0 == 0)
    fprintf(stdout, "Start: Read *def files.\n");
  ReadDefFileNInt(fileDefList, comm0);
  if (rank0 == 0)
    fprintf(stdout, "End  : Read *def files.\n");
  StopTimer(11);

  StartTimer(12);
  SetMemoryDef();
  StopTimer(12);

  StartTimer(11);
  if (rank0 == 0)
    fprintf(stdout, "Start: Read parameters from *def files.\n");
  ReadDefFileIdxPara(fileDefList, comm0);
  if (rank0 == 0)
    fprintf(stdout, "End  : Read parameters from *def files.\n");
  StopTimer(11);

  StartTimer(12);
  if (rank0 == 0)
    fprintf(stdout, "Start: Set memories.\n");
  SetMemory();
  if (rank0 == 0)
    fprintf(stdout, "End  : Set memories.\n");
  StopTimer(12);

  /* split MPI coummunicator */
#ifdef _mpi_use
  StartTimer(10);
  group1 = rank0 / NSplitSize;
  int rank1, rank2, size1, size2;
  MPI_Comm_split(comm0, group1, rank0, &comm1);
  MPI_Comm_size(comm1, &size1);
  MPI_Comm_rank(comm1, &rank1);
  int group2 = rank1;
  MPI_Comm_split(comm0, group2, rank0, &comm2);
  MPI_Comm_size(comm2, &size2);
  MPI_Comm_rank(comm2, &rank2);

  if (size0 % NSplitSize != 0 && rank0 == 0) {
    fprintf(stderr, "warning: load imbalance. MPI_size0=%d NSplitSize=%d\n",
            size0, NSplitSize);
  }
  /*   printf("rank=%d group1=%d rank1=%d rank2=%d size1=%d size2=%d\n", */
  /*      rank,group1,rank1,rank2,size1,size2); */
  StopTimer(10);
#endif

  /* initialize Mersenne Twister */
  init_gen_rand(RndSeed + group1);
  /* get the size of work space for LAPACK and PFAPACK */
  LapackLWork = getLWork_fcmp(); // TBC

  StartTimer(13);
  /* initialize variational parameters */
  if (rank0 == 0)
    fprintf(stdout, "Start: Initialize parameters.\n");
  InitParameter(); /* Run parallelly for synchronization of random generator */
  if (flagReadInitPara > 0 && rank0 == 0)
    ReadInitParameter(fileInitPara);
  //[s] add read parameters respectively
  if (rank0 == 0) {
    if (!ReadInputParameters(fileDefList, comm0) == 0) {
      //[ToDo]: Add Error procedure
      info = 1;
    }
    fprintf(stdout, "End  : Initialize parameters.\n");
  }
  //[e] add read parameters respectively

  SyncModifiedParameter(comm0);
  StopTimer(13);

  /* initialize variables for quantum projection */
  if (rank0 == 0)
    fprintf(stdout, "Start: Initialize variables for quantum projection.\n");
  InitQPWeight();
  if (rank0 == 0)
    fprintf(stdout, "End  : Initialize variables for quantum projection.\n");
  /* initialize output files */
  if (rank0 == 0)
    InitFile(fileDefList, rank0);

  StopTimer(1);

  if (NVMCCalMode == 0) {
    StartTimer(2);
    /*-- VMC Parameter Optimization --*/
    if (rank0 == 0)
      fprintf(stdout, "Start: Optimize VMC parameters.\n");
    info = VMCParaOpt(comm0, comm1, comm2);
    if (rank0 == 0)
      fprintf(stdout, "End  : Optimize VMC parameters.\n");
    StopTimer(2);
  } else if ((NVMCCalMode == 1) || (NVMCCalMode == 2) || (NVMCCalMode == 3)) {
    StartTimer(2);
    /*-- VMC Physical Quantity Calculation --*/
    if (rank0 == 0)
      fprintf(stdout, "Start: Calculate VMC physical quantities.\n");
    // if(NVMCCalMode==3)   read_StdFace_L_W();
    info = VMCPhysCal(comm0, comm1, comm2);
    if (rank0 == 0)
      fprintf(stdout, "End  : Calculate VMC physical quantities.\n");
    StopTimer(2);
  } else {
    info = 1;
    if (rank0 == 0)
      fprintf(stderr, "error: NVMCCalMode must be 0 1 2 or 3.\n");
  }

  StopTimer(0);
  if (rank0 == 0) {
    if (NVMCCalMode == 0) {
      OutputTimerParaOpt();
    } else if (NVMCCalMode == 1) {
      OutputTimerPhysCal();
    }

    /* close output files */
    CloseFile(rank0);
    fprintf(stdout, "Start: Free Memory.\n");
  }
  FreeMemory();
  FreeMemoryDef();
  if (rank0 == 0)
    fprintf(stdout, "End: Free Memory.\n");

  MPI_Finalize();
  if (rank0 == 0)
    fprintf(stdout, "Finish calculation.\n");

  return info;
}

/*-- VMC Parameter Optimization --*/
int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1,
               MPI_Comm comm_child2) {
  int step;
  int info;
  int rank;
  int tmp_i; // DEBUG
  int iprogress;
  MPI_Comm_rank(comm_parent, &rank);

  for (step = 0; step < NSROptItrStep; step++) {
    // printf("0 DUBUG make:step=%d TwoSz=%d\n",step,TwoSz);
    if (rank == 0) {
      OutputTime(step);
      if (NSROptItrStep < 20) {
        iprogress = (int)(100.0 * step / NSROptItrStep);
        printf("Progress of Optimization: %d %%.\n", iprogress);
      } else if (step % (NSROptItrStep / 20) == 0) {
        iprogress = (int)(100.0 * step / NSROptItrStep);
        printf("Progress of Optimization: %d %%.\n", iprogress);
      }
    }

    StartTimer(20);
    // printf("1 DUBUG make:step=%d \n",step);
    assert(iFlgOrbitalGeneral == 0); // sz is conserved
    UpdateSlaterElm_fcmp();
    // printf("2 DUBUG make:step=%d \n",step);
    UpdateQPWeight();
    StopTimer(20);
    StartTimer(3);
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, MakeSample.\n", step);
#endif
    assert(AllComplexFlag == 0 && iFlgOrbitalGeneral == 0); // real & sz=0
    // only for real TBC
    StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0; tmp_i < NQPFull * (2 * Nsite) * (2 * Nsite); tmp_i++)
      SlaterElm_real[tmp_i] = creal(SlaterElm[tmp_i]);
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0; tmp_i < NQPFull * (Nsize * Nsize + 1); tmp_i++)
      InvM_real[tmp_i] = creal(InvM[tmp_i]);
    StopTimer(69);
    assert(NProjBF == 0);
    // SlaterElm_real will be used in CalculateMAll, note that SlaterElm will
    // not change before SR
    if (ExactSample) {
      MakeExactSample(comm_child1);
    } else {
      // printf("Generating Sample via MC\n");
      VMCMakeSample_real(comm_child1);
    }
    // only for real TBC
    StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0; tmp_i < NQPFull * (Nsize * Nsize + 1); tmp_i++)
      InvM[tmp_i] = InvM_real[tmp_i] + 0.0 * I;
    StopTimer(69);
    // only for real TBC

    StopTimer(3);
    StartTimer(4);
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, MainCal.\n", step);
#endif
    assert(NProjBF == 0);
    assert(iFlgOrbitalGeneral == 0); // sz is conserved
    VMCMainCal(comm_child1);
    StopTimer(4);
    StartTimer(21);
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, AverageWE.\n", step);
#endif
    WeightAverageWE(comm_parent);
    StartTimer(25); // DEBUG
#ifdef _DEBUG_DETAIL
    printf("Debug: step %d, SROpt.\n", step);
#endif
    assert(AllComplexFlag == 0 && iFlgOrbitalGeneral == 0); // real & sz =0
    WeightAverageSROpt_real(comm_parent);
    StopTimer(25);
    ReduceCounter(comm_child2);
    StopTimer(21);
    StartTimer(22);
    /* output zvo_out and zvo_var */
    if (rank == 0)
      outputData();
    StopTimer(22);

#ifdef _DEBUG_DUMP_SROPTO_STORE
    if (rank == 0) {
      assert(AllComplexFlag == 0 && iFlgOrbitalGeneral == 0); // real & sz=0
      for (i = 0; i < SROptSize * NVMCSample; i++) {
        fprintf(stderr, "DEBUG: SROptO_Store_real[%d]=%lf +I*%lf\n", i,
                creal(SROptO_Store_real[i]), cimag(SROptO_Store_real[i]));
      }
    }
#endif

#ifdef _DEBUG_DUMP_SROPTOO
    if (rank == 0) {
      assert(AllComplexFlag == 0 && iFlgOrbitalGeneral == 0); // real & sz=0
      for (i = 0; i < (NSRCG == 0 ? SROptSize * SROptSize : SROptSize * 2);
           i++) {
        fprintf(stderr, "DEBUG: SROptOO_real[%d]=%lf +I*%lf\n", i,
                creal(SROptOO_real[i]), cimag(SROptOO_real[i]));
      }
      for (i = 0; i < SROptSize; i++) {
        fprintf(stderr, "DEBUG: SROptHO_real[%d]=%lf +I*%lf\n", i,
                creal(SROptHO_real[i]), cimag(SROptHO_real[i]));
      }
      for (i = 0; i < SROptSize; i++) {
        fprintf(stderr, "DEBUG: SROptO_real[%d]=%lf +I*%lf\n", i,
                creal(SROptO_real[i]), cimag(SROptO_real[i]));
      }
    }
#endif

    StartTimer(5);
    if (NSRCG != 0) {
      info = StochasticOptCG(comm_parent);
    } else {
      info = StochasticOpt(comm_parent);
    }
    // info = StochasticOptDiag(comm_parent);
    StopTimer(5);

#ifdef _DEBUG_DUMP_PARA
    for (int i = 0; i < NPara; ++i) {
      fprintf(stderr, "DEBUG: Para[%d] = %lf %lf\n", i, creal(Para[i]),
              cimag(Para[i]));
    }
#endif

    // DEBUG
    // abort();

    if (info != 0) {
      if (rank == 0)
        fprintf(stderr, "Error: StcOpt info=%d step=%d\n", info, step);
      return info;
    }

    StartTimer(23);
    SyncModifiedParameter(comm_parent);
    StopTimer(23);

    if (step >= NSROptItrStep - NSROptItrSmp) {
      StoreOptData(step - (NSROptItrStep - NSROptItrSmp));
    }

    FlushFile(step, rank);
  }

  /* output zqp_opt */
  if (rank == 0) {
    OutputTime(NSROptItrStep);
    fprintf(stdout, "Start: Output opt params.\n");
    OutputOptData();
    fprintf(stdout, "End: Output opt params.\n");
  }

  return 0;
}

/*-- VMC Physical Quantity Calculation --*/
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1,
               MPI_Comm comm_child2) {
  int ismp, tmp_i;
  int rank;
  MPI_Comm_rank(comm_parent, &rank);

  if (rank == 0)
    fprintf(stdout, "Start: UpdateSlaterElm.\n");
  StartTimer(20);
  assert(iFlgOrbitalGeneral == 0); // sz is conserved
  UpdateSlaterElm_fcmp();
  StopTimer(20);
  if (rank == 0) {
    fprintf(stdout, "End  : UpdateSlaterElm.\n");
    fprintf(stdout, "Start: Sampling.\n");
  }
  for (ismp = 0; ismp < NDataQtySmp; ismp++) {
    if (rank == 0) {
      printf("measurement %d\n", ismp + 1);
    }
    if (rank == 0)
      OutputTime(ismp);
    FlushFile(0, rank);
    InitFilePhysCal(ismp, rank);
    StartTimer(3);
    assert(NProjBF == 0);
    assert(AllComplexFlag == 0 && iFlgOrbitalGeneral == 0); // real & sz=0
    // only for real TBC
    StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0; tmp_i < NQPFull * (2 * Nsite) * (2 * Nsite); tmp_i++)
      SlaterElm_real[tmp_i] = creal(SlaterElm[tmp_i]);
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0; tmp_i < NQPFull * (Nsize * Nsize + 1); tmp_i++)
      InvM_real[tmp_i] = creal(InvM[tmp_i]);
    StopTimer(69);
    // SlaterElm_real will be used in CalculateMAll, note that SlaterElm will
    // not change before SR
    if (ExactSample) {
      MakeExactSample(comm_child1);
    } else {
      VMCMakeSample_real(comm_child1);
    }
    // only for real TBC
    StartTimer(69);
#pragma omp parallel for default(shared) private(tmp_i)
    for (tmp_i = 0; tmp_i < NQPFull * (Nsize * Nsize + 1); tmp_i++)
      InvM[tmp_i] = InvM_real[tmp_i] + 0.0 * I;
    StopTimer(69);
    // only for real TBC

    StopTimer(3);
    StartTimer(4);
    if (rank == 0) {
      fprintf(stdout, "End  : Sampling.\n");
      fprintf(stdout, "Start: Main calculation.\n");
    }
    assert(NProjBF == 0);
    assert(iFlgOrbitalGeneral == 0);
    VMCMainCal(comm_child1);
    if (PrintProgress) {
      if (rank == 0)
        fprintf(stdout, "\nEnd  : Main calculation.\n");
    } else {
      if (rank == 0)
        fprintf(stdout, "End  : Main calculation.\n");
    }

    StopTimer(4);

    if (NVMCCalMode == 1) {
      StartTimer(21);
      WeightAverageWE(comm_parent);
      WeightAverageGreenFunc(comm_parent);
      ReduceCounter(comm_child2);

      StopTimer(21);
      StartTimer(22);
      /* output zvo_out and green functions */
      if (rank == 0)
        outputData();
      CloseFilePhysCal(rank);

      StopTimer(22);
      // printf("hey5\n");
    } else if (NVMCCalMode == 2) {
      // printf("before averaging.\n"); fflush(stdout);
      WeightAverageWE(comm_parent);
      WeightAverageStaticQuantities_real(comm_parent);
      ReduceCounter(comm_child2);
      if (rank == 0) {
        outputData();
        fclose(FileEn);
        fclose(FileN1);
        fclose(FileN2);
      }
    } else if (NVMCCalMode == 3) {
      // printf("before averaging.\n"); fflush(stdout);
      WeightAverageWE(comm_parent);
      WeightAverageDynamicalGreenFunc_real(comm_parent);
      ReduceCounter(comm_child2);
      if (rank == 0) {
        outputData();
        // fclose(File_nCHAm_bin);
      }
    }
    StopTimer(5);
  }

  if (rank == 0)
    OutputTime(NDataQtySmp);

  return 0;
}

void outputData() {
  int i, j;

  // zvo_out.dat
  //[s] MERGE BY TM
  // fprintf(FileOut, "% .18e % .18e % .18e \n", Etot, Etot2, (Etot2 -
  // Etot*Etot)/(Etot*Etot));
  //   fprintf(FileOut, "% .18e % .18e  % .18e % .18e \n",
  //   creal(Etot),cimag(Etot), creal(Etot2), creal((Etot2 -
  //   Etot*Etot)/(Etot*Etot)));
  if (NVMCCalMode < 2) {
    fprintf(FileOut, "% .18e % .18e  % .18e % .18e %.18e %.18e\n", creal(Etot),
            cimag(Etot), creal(Etot2),
            creal((Etot2 - Etot * Etot) / (Etot * Etot)), creal(Sztot),
            creal(Sztot2));
    // fprintf(FileOut, "% .18e % .18e % .18e \n", Etot, Etot2, (Etot2 -
    // Etot*Etot)/(Etot*Etot));
    // fprintf(FileOut, "% .18e % .18e  % .18e % .18e \n", creal(Etot),
    // cimag(Etot), creal(Etot2),
    //         creal((Etot2 - Etot * Etot) / (Etot * Etot)));
    //[e] MERGE BY TM

    // zvo_var.dat
    if (FlagBinary == 0) { // formatted output
      fprintf(FileVar, "% .18e % .18e 0.0 % .18e % .18e 0.0 ", creal(Etot),
              cimag(Etot), creal(Etot2), cimag(Etot2));
      for (i = 0; i < NPara; i++)
        fprintf(FileVar, "% .18e % .18e 0.0 ", creal(Para[i]), cimag(Para[i]));
      fprintf(FileVar, "\n");
      // for(i=0;i<NPara;i++)  printf("DEBUG:i=%d: % .18e % .18e  \n",i,
      // creal(Para[i]),cimag(Para[i]));
    } else { // binary output
      fwrite(Para, sizeof(double), NPara, FileVar);
    }
  }

  if (NVMCCalMode == 1) {
    // zvo_cisajs.dat
    if (NCisAjs > 0) {
      if (NLanczosMode < 2) {
        for (i = 0; i < NCisAjs; i++) {
          fprintf(FileCisAjs, "%d %d %d %d % .18e  % .18e \n", CisAjsIdx[i][0],
                  CisAjsIdx[i][1], CisAjsIdx[i][2], CisAjsIdx[i][3],
                  creal(PhysCisAjs[i]), cimag(PhysCisAjs[i]));
        }
      } else {
        for (i = 0; i < NCisAjsLz; i++) {
          int idx = iOneBodyGIdx[CisAjsLzIdx[i][0] + CisAjsLzIdx[i][1] * Nsite]
                                [CisAjsLzIdx[i][2] + CisAjsLzIdx[i][3] * Nsite];
          // fprintf(stdout, "Debug: idx= %d value= % .18e % .18e\n", idx,
          // creal(PhysCisAjs[idx]), cimag(PhysCisAjs[idx]));
          fprintf(FileCisAjs, "%d %d %d %d % .18e % .18e \n",
                  CisAjsLzIdx[idx][0], CisAjsLzIdx[idx][1], CisAjsLzIdx[idx][2],
                  CisAjsLzIdx[idx][3], creal(PhysCisAjs[idx]),
                  cimag(PhysCisAjs[idx]));
        }
      }
      fprintf(FileCisAjs, "\n");
    }
    // zvo_cisajscktalt.dat
    if (NCisAjsCktAlt > 0) {
      for (i = 0; i < NCisAjsCktAlt; i++)
        fprintf(FileCisAjsCktAlt, "% .18e  % .18e ", creal(PhysCisAjsCktAlt[i]),
                cimag(PhysCisAjsCktAlt[i]));
      fprintf(FileCisAjsCktAlt, "\n");
    }

    // zvo_cisajscktaltdc.dat
    if (NCisAjsCktAltDC > 0) {
      for (i = 0; i < NCisAjsCktAltDC; i++) {
        fprintf(FileCisAjsCktAltDC, "%d %d %d %d %d %d %d %d % .18e % .18e\n",
                CisAjsCktAltDCIdx[i][0], CisAjsCktAltDCIdx[i][1],
                CisAjsCktAltDCIdx[i][2], CisAjsCktAltDCIdx[i][3],
                CisAjsCktAltDCIdx[i][4], CisAjsCktAltDCIdx[i][5],
                CisAjsCktAltDCIdx[i][6], CisAjsCktAltDCIdx[i][7],
                creal(PhysCisAjsCktAltDC[i]), cimag(PhysCisAjsCktAltDC[i]));
      }
      fprintf(FileCisAjsCktAltDC, "\n");
    }
    assert(NLanczosMode <= 0);
  }

  else if (NVMCCalMode == 2) {
    fprintf(FileEn, "% .18e % .18e  % .18e % .18e %.18e %.18e\n", creal(Etot),
            cimag(Etot), creal(Etot2),
            creal((Etot2 - Etot * Etot) / (Etot * Etot)), creal(Sztot),
            creal(Sztot2));
    // printf("trying to print files.\n"); fflush(stdout);
    if (NCisAjs > 0) {
      // for (i = 0; i < NCisAjs; i++) {
      //     fprintf(FileCisAjs, "%d %d %d %d % .18e  % .18e \n",
      //     CisAjsIdx[i][0], CisAjsIdx[i][1], CisAjsIdx[i][2],
      //             CisAjsIdx[i][3], creal(PhysCisAjs[i]),
      //             cimag(PhysCisAjs[i]));
      // }
      for (i = 0; i < NCisAjs; i++) {
        fprintf(FileN2, "%d %d %d %d ", CisAjsIdx[i][0], CisAjsIdx[i][1],
                CisAjsIdx[i][2], CisAjsIdx[i][3]);
        for (j = 0; j < TWO_SITES_PHYS_QTY; j++)
          fprintf(FileN2, "% 0.6e   ", PhysN2[i + NCisAjs * j]);
        //          fprintf(FileN, "% .8e  % .8e   ",
        //          creal(PhysN2[i+NCisAjs*j]), cimag(PhysN2[i+NCisAjs*j]));
        fprintf(FileN2, "\n");
      }
    }
    for (i = 0; i < Nsite; i++) {
      fprintf(FileN1, "%d ", i);
      for (j = 0; j < ONE_SITE_PHYS_QTY; j++)
        fprintf(FileN1, "% 0.2e   ", PhysN1[i + Nsite * j]);
      fprintf(FileN1, "\n");
    }
    fprintf(FileN1, "\n");
    printf("ending of print files.\n");
    fflush(stdout);
  } else if (NVMCCalMode == 3) {
    // print binary output:
    long int totalSize = NExcitation * NExcitation * Nsite * Nsite;
    fwrite(&NExcitation, sizeof(int), 1, File_nCHAm_bin);
    fwrite(&Excitation_L, sizeof(int), 1, File_nCHAm_bin);
    fwrite(&Excitation_W, sizeof(int), 1, File_nCHAm_bin);

    convert_double2float_fwrite(Phys_nCAm_averaged, totalSize, File_nCHAm_bin);
    convert_double2float_fwrite(Phys_nACm_averaged, totalSize, File_nCHAm_bin);
    convert_double2float_fwrite(Phys_nCHAm_averaged, totalSize, File_nCHAm_bin);
    convert_double2float_fwrite(Phys_nAHCm_averaged, totalSize, File_nCHAm_bin);
  }
  return;
}
