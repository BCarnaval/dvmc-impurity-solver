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
 * initialization of files
 *-------------------------------------------------------------
 * by Satoshi Morita 
 *-------------------------------------------------------------*/

#ifndef _INITFILE_SRC
#define _INITFILE_SRC

void InitFile(char *xNameListFile, int rank) {
  char fileName[D_FileNameMaxL];

  if(rank!=0) return;

  //sprintf(fileName, "%s_cfg_%03d.dat", CDataFileHead, NDataIdxStart);
  //writeConfig(xNameListFile, fileName);

  sprintf(fileName, "%s_time_%03d.dat", CDataFileHead, NDataIdxStart);
  FileTime = fopen(fileName, "w");

  if(NVMCCalMode==0) {
    sprintf(fileName, "%s_SRinfo.dat", CDataFileHead);
    FileSRinfo = fopen(fileName, "w");
    if(SRFlag == 0){
      fprintf(FileSRinfo,
            "#Npara Msize optCut diagCut sDiagMax  sDiagMin    absRmax       imax\n");
    }else{
      fprintf(FileSRinfo,
            "#Npara Msize optCut diagCut sEigenMax  sEigenMin    absRmax       imax\n");
    }

    sprintf(fileName, "%s_out_%03d.dat", CDataFileHead, NDataIdxStart);
    FileOut = fopen(fileName, "w");

    if(FlagBinary==0) {
      sprintf(fileName, "%s_var_%03d.dat", CDataFileHead, NDataIdxStart);
      FileVar = fopen(fileName, "w");
    } else {
      sprintf(fileName, "%s_varbin_%03d.dat", CDataFileHead, NDataIdxStart);
      FileVar = fopen(fileName, "wb");
      fwrite(&NPara,sizeof(int),1,FileVar);
      fwrite(&NSROptItrStep,sizeof(int),1,FileVar);
    }

    if(PrintEnergy){
      sprintf(fileName,"Energy.dat");
      File_E = fopen(fileName, "w");
    }

  }

  return;
}

void InitFilePhysCal(int i, int rank) {
  char fileName[D_FileNameMaxL];
  int idx = i+NDataIdxStart;
  int one = 1;

  if(rank!=0) return;

  if(NVMCCalMode<2){
    sprintf(fileName, "%s_out_%03d.dat", CDataFileHead, idx);
    FileOut = fopen(fileName, "w");
   
    if(FlagBinary==0) {
      sprintf(fileName, "%s_var_%03d.dat", CDataFileHead, idx);
      FileVar = fopen(fileName, "w");
    } else {
      sprintf(fileName, "%s_varbin_%03d.dat", CDataFileHead, idx);
      FileVar = fopen(fileName, "wb");
      fwrite(&NPara,sizeof(int),1,FileVar);
      fwrite(&one,sizeof(int),1,FileVar);
    }
  }

  if(NVMCCalMode==1){
  /* Green function */
    if(NCisAjs>0){
      sprintf(fileName, "%s_cisajs_%03d.dat", CDataFileHead, idx);
      FileCisAjs = fopen(fileName, "w");
    }
   

    if(NCisAjsCktAlt>0){
      sprintf(fileName, "%s_cisajscktaltex_%03d.dat", CDataFileHead, idx);
      FileCisAjsCktAlt = fopen(fileName, "w");
    }

    if(NCisAjsCktAltDC>0){
      sprintf(fileName, "%s_cisajscktalt_%03d.dat", CDataFileHead, idx);
      FileCisAjsCktAltDC = fopen(fileName, "w");
    }
    
    if(NLanczosMode>0){
      sprintf(fileName, "%s_ls_out_%03d.dat", CDataFileHead, idx);
      FileLS = fopen(fileName, "w");

      sprintf(fileName, "%s_ls_qqqq_%03d.dat", CDataFileHead, idx);
      FileLSQQQQ = fopen(fileName, "w");
      
      if(NLanczosMode>1){
  #ifdef _DEBUG
        sprintf(fileName, "%s_ls_qcisajsq_%03d.dat",
                CDataFileHead, idx);
        FileLSQCisAjsQ = fopen(fileName, "w");
        sprintf(fileName, "%s_ls_qcisajscktaltq_%03d.dat",
                CDataFileHead, idx);
        FileLSQCisAjsCktAltQ = fopen(fileName, "w");

   #endif
        sprintf(fileName, "%s_ls_cisajs_%03d.dat",
                CDataFileHead, idx);
        FileLSCisAjs = fopen(fileName, "w");

        sprintf(fileName, "%s_ls_cisajscktalt_%03d.dat",
                CDataFileHead, idx);
        FileLSCisAjsCktAlt = fopen(fileName, "w");
      }
    }
  
  }
  else if(NVMCCalMode==2){
    //if(NCisAjs>0){
    //  sprintf(fileName, "%s_cisajs_%03d.dat", CDataFileHead, idx);
    //  FileCisAjs = fopen(fileName, "w");
    //}

    if(PrintEnergy){
      sprintf(fileName,"Energy_%03d.dat",idx);
      File_E = fopen(fileName, "w");
    }

    if(PrintConfig){
      sprintf(fileName,"Config_%03d.dat",idx);
      File_Config = fopen(fileName, "w");
    }
    
    sprintf(fileName, "%s_energy_%03d.dat", CDataFileHead, idx);
    FileEn = fopen(fileName, "w");
    
    sprintf(fileName, "%s_physN1_%03d.dat", CDataFileHead, idx);
    FileN1 = fopen(fileName, "w");
    
    sprintf(fileName, "%s_physN2_%03d.dat", CDataFileHead, idx);
    FileN2 = fopen(fileName, "w");
/*
    sprintf(fileName, "%s_physAC_%03d.dat", CDataFileHead, idx);
    File_AC = fopen(fileName, "w");
    
    sprintf(fileName, "%s_physACN_%03d.dat", CDataFileHead, idx);
    File_ACN = fopen(fileName, "w");
    
    sprintf(fileName, "%s_physNACN_%03d.dat", CDataFileHead, idx);
    File_NACN = fopen(fileName, "w");

    sprintf(fileName, "%s_phys_nCHAm_%03d.dat", CDataFileHead, idx);
    File_nCHAm = fopen(fileName, "w");*/
  }
  else if(NVMCCalMode==3){
    //if(NCisAjs>0){
    //  sprintf(fileName, "%s_cisajs_%03d.dat", CDataFileHead, idx);
    //  FileCisAjs = fopen(fileName, "w");
    //}

    //sprintf(fileName, "%s_nCHAm_nAHCm_%03d.dat", CDataFileHead, idx);
    //File_nCHAm = fopen(fileName, "w");

    sprintf(fileName, "%s_nCHAm_nAHCm_%03d.bin", CDataFileHead, idx);
    File_nCHAm_bin = fopen(fileName, "wb");

    //sprintf(fileName, "%s_nCAm_nACm_nCHAm_nAHCm_%03d.dat", CDataFileHead, idx);
    //File_nCHAm_compact = fopen(fileName, "w");
  }
  return;
}

void CloseFile(int rank) {
  if(rank!=0) return;

  fclose(FileTime);

  if(NVMCCalMode==0) {
    fclose(FileSRinfo);
    fclose(FileOut);
    fclose(FileVar);
    if(PrintEnergy)fclose(File_E);
  }

  return;
}

void CloseFilePhysCal(int rank) {
  if(rank!=0) return;

  if(NVMCCalMode<2) {
    fclose(FileOut);
    fclose(FileVar);
  }

  if(PrintEnergy) fclose(File_E);
  if(PrintConfig) fclose(File_Config);
  
  if(NCisAjs>0){
    fclose(FileCisAjs);
  }
  if(NCisAjsCktAlt>0){
    fclose(FileCisAjsCktAlt);
  }
  if(NCisAjsCktAltDC>0){
    fclose(FileCisAjsCktAltDC);
  }
  
  if(NLanczosMode>0){
    fclose(FileLS);
    fclose(FileLSQQQQ);
    
    if(NLanczosMode>1){
#ifdef _DEBUG
      fclose(FileLSQCisAjsQ);
      fclose(FileLSQCisAjsCktAltQ);
#endif
      fclose(FileLSCisAjs);
      fclose(FileLSCisAjsCktAlt);
    }
  }

  return;
}

void FlushFile(int step, int rank) {
  if(rank!=0) return;

  if(step%NFileFlushInterval==0) {
    fflush(FileTime);
    if(NVMCCalMode==0) {
      fflush(FileSRinfo);
      fflush(FileOut);
      fflush(FileVar);
    }
  }
  return;
}


#endif
