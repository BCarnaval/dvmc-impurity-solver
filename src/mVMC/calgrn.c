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
 * Cauculate Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#ifndef _CALGRN_SRC
#define _CALGRN_SRC



void CalculateStaticQuantities_real(const double w, 
                               int *eleIdx, int *eleCfg,
                               int *eleNum, int *eleProjCnt) {

  int idx;
  int ri,rj,s;
  int *myEleNum;
  
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  {
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      rj = CisAjsIdx[idx][2];
      s  = CisAjsIdx[idx][3];
      
      //Doublon-Holon TJS
      PhysN2[idx+NCisAjs*0] += w*myEleNum[ri+s*Nsite]*myEleNum[ri+(1-s)*Nsite]
                                   *(1.0-myEleNum[rj+s*Nsite])*(1.0-myEleNum[rj+(1-s)*Nsite]);
                                   
      //Doublon-Doublon TJS
      PhysN2[idx+NCisAjs*1] += w*myEleNum[ri+s*Nsite]*myEleNum[ri+(1-s)*Nsite]
                                      *myEleNum[rj+s*Nsite]*myEleNum[rj+(1-s)*Nsite];

      //Charge-Doublon TJS
      PhysN2[idx+NCisAjs*2] += w*myEleNum[ri+s*Nsite] *myEleNum[rj+s*Nsite]*myEleNum[rj+(1-s)*Nsite];
      
      //n_sigma (1-n_sigma) MC
      PhysN2[idx+NCisAjs*3] += w*myEleNum[ri+s*Nsite] *(1.0-myEleNum[rj+s*Nsite]);

      //S^z_i S^z_j
      PhysN2[idx+NCisAjs*4] += w*(myEleNum[ri]-myEleNum[ri+Nsite])*(myEleNum[rj]-myEleNum[rj+Nsite]);
    }

    for(ri=0;ri<Nsite;ri++) {
      s=0;
      //Doublon MC
      PhysN1[ri+Nsite*0] += w*myEleNum[ri]*myEleNum[ri+Nsite];

      //Holon MC
      PhysN1[ri+Nsite*1] += w*(1.0-myEleNum[ri])*(1.0-myEleNum[ri+Nsite]);
            
      //Density_up (s==0) MC
      PhysN1[ri+Nsite*2] += w*myEleNum[ri];

      //Density_down (s==1) MC
      PhysN1[ri+Nsite*3] += w*myEleNum[ri+Nsite];

    }
  }
  ReleaseWorkSpaceThreadInt();
  return;
}




int del(int i,int j){
  return ((i==j)? 1:0);
}


// This function defines every type of excitation (0 to 7) and its 
// anticommutation with different combination of creation operator (C)
// and annihilation operator (A)

// Charlebois and Imada 2019, arxiv v1 (equation B5)
int Commute_Nat_(commuting_with commuting, int ra, int rb, int t, int ri, int rj, int s, int rm, int rn, int u, int *eleNum) {
  int sign;
  if((commuting==with_CisCmuAnuAjs) || (commuting==with_AisCmuAnuCjs)) {
    if(commuting==with_CisCmuAnuAjs) {sign = 1;}
    else if(commuting==with_AisCmuAnuCjs) {sign = -1;}
    
    if(t==0){ 
      return 1;
    }
    else if(t==1){ 
      return (eleNum[ra+(1-s)*Nsite] + del((1-s),u) * (del(ra,rm)-del(ra,rn)));
    }
    else if(t==2){
      if(ra==rb){
        return (eleNum[ra +s *Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,s)    + sign * (del(ra,ri)-del(ra,rj)));      
      }
      return (eleNum[ra+   s *Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,s)    + sign * (del(ra,ri)-del(ra,rj)) )
            *(eleNum[rb+   s *Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,s)    + sign * (del(rb,ri)-del(rb,rj)) );
    }
    else if(t==3){
      return (eleNum[ra+(1-s)*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,1-s))
            *(eleNum[rb+(1-s)*Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,1-s));
    }
    else if(t==4){
      return (eleNum[ra+(1-s)*Nsite] + (del(ra,rm)-del(ra,rn)) * del(u,1-s))
            *(eleNum[rb+   s *Nsite] + (del(rb,rm)-del(rb,rn)) * del(u,s)    + sign * (del(rb,ri)-del(rb,rj)) )   ;
    }
    else{
      printf("oups, error %d\n", t);
      exit(1);
      return -1;
    }
  
  }
  else if((commuting==with_CisAjs) || (commuting==with_AisCjs)) {
  
    if(commuting==with_CisAjs) {sign = 1;}
    else if(commuting==with_AisCjs) {sign = -1;}
    
    if(t==0){ //no charge
      return 1;
    }
    else if(t==1){ // electron reverse-spin
      return (eleNum[ra+(1-s)*Nsite]) ;
    }
    else if(t==2){ 
      if(ra==rb) {
        return (eleNum[ra+    s*Nsite] + sign * (del(ra,ri)-del(ra,rj) ) );
      }
        return (eleNum[ra+    s*Nsite] + sign * (del(ra,ri)-del(ra,rj) ) )
              *(eleNum[rb+    s*Nsite] + sign * (del(rb,ri)-del(rb,rj) ) );      
    }
    else if(t==3){ 
      return (eleNum[ra+(1-s)*Nsite] )
            *(eleNum[rb+(1-s)*Nsite] );
    }
    else if(t==4){ 
      return (eleNum[ra+(1-s)*Nsite] )
            *(eleNum[rb+    s*Nsite] + sign * (del(rb,ri)-del(rb,rj) ) );
    }
    else{
      printf("oups, error %d\n", t);
      exit(1);
      return -1;
    }
  }
  else if(commuting==with_nothing) {
    if(t==0){ //no charge
      return 1;
    }    
    else if(t==1){ //1 electron reverse-spin
      return (eleNum[ra+(1-s)*Nsite]);
    }
    else if(t==2){
      return (eleNum[ra+s*Nsite]*eleNum[rb+s*Nsite]);
    }
    else if(t==3){
      return (eleNum[ra+(1-s)*Nsite]*eleNum[rb+(1-s)*Nsite]);
    }
    else if(t==4){
      return (eleNum[ra+(1-s)*Nsite]*eleNum[rb+s*Nsite]);
    }
    else{
      printf("oups, error %d\n", t);
      exit(1);
      return -1;
    }
  }
  exit(1);
  return -2; 
}


//this function is useful to reduce memory by a factor of two.
//we could argue that we are loosing precision. But in a Monte Carlo
//sampling, it is better to keep double for the whole calculation until 
//the end, in order to minimize the discretization noise addeed at each 
//sample. But the result at the end still have quite a lot of noise, so
//it is ok to reduce the final result to a float when saving in a binary file.
//Maxime Charlebois
void convert_double2float_fwrite(double *array, long int totalSize, FILE *fp) {
  long int ii;
  for(ii=0;ii<totalSize;ii++) data_float_convert[ii] = (float) array[ii];
  fwrite(data_float_convert, sizeof(float), totalSize, fp);
  return;
}


void CalculateDynamicalGreenFunc_real(const double w, const double ip, 
                                      int *eleIdx, int *eleCfg,
                                      int *eleNum, int *eleProjCnt, int sample) {
  int idx;
  int ri,rj,s;
  int *myEleIdx, *myEleNum, *myProjCntNew; //, *myEleCfg
  double *myBuffer_real;
  
  //printf("Begin CalculateDynamicalGreenFunc_real\n");

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadDouble(NQPFull+2*Nsize);
  
#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer_real)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew   = GetWorkSpaceThreadInt(NProj);
    myBuffer_real  = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize; idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
    
    int rm, rn, u;
    int idx_int, idx_trans;

    //printf("\n");
    //printf("\n");
    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) 
    for(ri=0;ri<Nsite;ri++) {
      //printf("%d ",ri);
     for(rj=0;rj<Nsite;rj++) {
      for(s=0;s<2;s++) {
       // Charlebois and Imada 2019, arxiv v1 (equation B1)         
       Local_CA[ri+Nsite*rj+Nsite*Nsite*s] = GreenFunc1_real(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                                                              myProjCntNew,myBuffer_real); 
       //if(Local_CA[ri+Nsite*rj+Nsite*Nsite*s]>100){
       // printf("ri : %d, rj : %d, Local_CA : %f \n", ri, rj, Local_CA[ri+Nsite*rj+Nsite*Nsite*s]);
       //}
      }
     }
    }

    //printf("\n");
    #pragma omp for private(idx,ri,rj,s,idx_trans,rm,rn,u) schedule(dynamic) 
    for(ri=0;ri<Nsite;ri++) {
      //for(ri=0;ri<Nsite;ri++){
     //printf("%d ",ri);
     for(rj=0;rj<Nsite;rj++) {
       //for(rj=0;rj<Nsite;rj++) {
      s=0;

      for(idx_trans=0;idx_trans<NTransfer;idx_trans++) {
      
        rm = Transfer[idx_trans][0];
        rn = Transfer[idx_trans][2];
        u  = Transfer[idx_trans][3];

        // Charlebois and Imada 2019, arxiv v1 (equation B2)
        Local_CisAjsCmuAnu[idx_trans + NTransfer*(ri+Nsite*rj)] = GreenFunc2_real(ri,rj,rm,rn,s,u,ip,
                        myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer_real);
	//Local_CisAjsCmuAnu[idx_trans + NTransfer*(ri+Nsite*rj)] = GreenFunc2_real(ri,rj,rm,rn,s,u,ip,
        //                myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer_real);
      }
     }
    }
    
    //printf("\n");
    #pragma omp for private(idx,ri,rj,s,idx_trans,rm,rn,u) schedule(dynamic) 
    for(ri=0;ri<Nsite;ri++) {
     //printf("%d ",ri);
     for(rj=0;rj<Nsite;rj++) {
       int idx1 = ri+Nsite*rj;
       int idx2 = rj+Nsite*ri;
        
       int idx1_full = ri+Nsite*rj;
       int idx2_full = rj+Nsite*ri;

     
       s=0; // just doing spin up
       int idx_exc_mm, idx_exc_nn;
      
       for(idx_exc_mm=0;idx_exc_mm<NExcitation;idx_exc_mm++){        
        int file_line_idx_i = idx_exc_mm+NExcitation*ri; // BEWARE: order matters in the file
        int t   = ChargeExcitationIdx[file_line_idx_i][0];  // type
        int rii = ChargeExcitationIdx[file_line_idx_i][1];
        int ra1 = ChargeExcitationIdx[file_line_idx_i][2];
        int ra2 = ChargeExcitationIdx[file_line_idx_i][3];
                
        if(ri != rii) printf("something wrong in excitations? %d == %d \n", ri, rii);
                
        int idx_vector1 = idx_exc_mm + NExcitation*(sample + sampleChunk*idx1);
        int idx_vector2 = idx_exc_mm + NExcitation*(sample + sampleChunk*idx2);
        
        //printf("%d %d %d %d %d\n", ChargeExcitationIdx[idx_exc][0],ChargeExcitationIdx[idx_exc][1],ChargeExcitationIdx[idx_exc][2],ChargeExcitationIdx[idx_exc][3],ChargeExcitationIdx[idx_exc][4]);
        //printf("%d %d \n", Excitation_L, Excitation_W); fflush(stdout);
        //printf("%d %d %d %d \n", ra1,ra1,rb1,rb2);        
        //printf("%d %d %d \n", ChargeExcitationIdx[idx_exc][0],ChargeExcitationIdx[idx_exc][1],ChargeExcitationIdx[idx_exc][2]);
        
        // <phi|ca|x> / <phi|x>
        // Charlebois and Imada 2019, arxiv v1 (equation B7)
        double CA_tmp = Local_CA[idx1_full] * Commute_Nat_(with_CisAjs,  ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);        
        double O_CA = CA_tmp;        
        //O_CA_vec2[idx_vector2] = CA_tmp;//conj(CA_tmp);  // (equation B10)      

        // <phi|ac|x> / <phi|x> = delta_{ri,rj} * <phi|x> / <phi|x> - <phi|ca|x> / <phi|x>           <-- need to reverse indices
        // Charlebois and Imada 2019, arxiv v1 (equation B6)
        double AC_tmp  = del(ri,rj);
        AC_tmp -= Local_CA[idx2_full];
        AC_tmp *= Commute_Nat_(with_AisCjs, ra1, ra2, t, ri, rj, s, 0,0,0, myEleNum);
        double O_AC = AC_tmp;
        //O_AC_vec2[idx_vector2] = AC_tmp;//conj(AC_tmp);  // (equation B10)
        
        //
        // <phi|H_U|x> 
        double tmp_int_AHC=0.0;
        double tmp_int_CHA=0.0;
                
        for(idx_int=0;idx_int<NCoulombIntra;idx_int++) {
          rm = CoulombIntra[idx_int];
          double factor = ParaCoulombIntra[idx_int] *
                   ( 1.*del(rj,rm) + myEleNum[rm+s*Nsite])*myEleNum[rm+(1-s)*Nsite];
          tmp_int_AHC += factor;        

          factor = ParaCoulombIntra[idx_int] *
                   ( -1.*del(rj,rm) + myEleNum[rm+s*Nsite])*myEleNum[rm+(1-s)*Nsite];
          tmp_int_CHA += factor;
        }
        
        double H_AC = tmp_int_AHC * AC_tmp ;
        double H_CA = tmp_int_CHA * CA_tmp ;
        
        // <phi|H_T|x> / <phi|x>
        for(idx_trans=0;idx_trans<NTransfer;idx_trans++) {
          //*
          
          rm = Transfer[idx_trans][0];
          rn = Transfer[idx_trans][2];
          u  = Transfer[idx_trans][3];
          
          int idx_green0 = ri+Nsite*rn;
          // Charlebois and Imada 2019, arxiv v1 (equation B9)
          double tmp = -1.0 * ParaTransfer[idx_trans] 
                     * (Local_CisAjsCmuAnu[idx_trans + NTransfer*idx1] 
                        - del(rm,rj) * del(s,u) * Local_CA[idx_green0] ) 
                     * Commute_Nat_(with_CisCmuAnuAjs, ra1, ra2, t, ri, rj, s, rm, rn, u, myEleNum) ;
          H_CA += tmp;
          //H_CA_vec2[idx_vector2] += tmp;//conj(tmp); // (equation B11)
          
          
          int idx_green1 = rj+Nsite*ri;
          int idx_green1_full = rj+Nsite*ri;
          int idx_green2 = rm+Nsite*rn+Nsite*Nsite*u;
          int idx_green3 = rm+Nsite*ri;
          // Charlebois and Imada 2019, arxiv v1 (equation B8)
          tmp = -1.0 * ParaTransfer[idx_trans]
                     * ( - Local_CisAjsCmuAnu[idx_trans + NTransfer*idx_green1] 
                         + del(ri,rj) * Local_CA[idx_green2] 
                         + del(rn,rj) * del(s,u) * ( del(rm,ri) - Local_CA[idx_green3] ))
                     * Commute_Nat_(with_AisCmuAnuCjs, ra1, ra2, t, ri, rj, s, rm, rn, u, myEleNum) ;
          H_AC += tmp;
          //H_AC_vec2[idx_vector2] += tmp;//conj(tmp); // (equation B11)
          
          
        }        
        
        for(idx_exc_nn=0;idx_exc_nn<NExcitation;idx_exc_nn++){                
        // <phi|x>
          int file_line_idx_j = idx_exc_nn+NExcitation*rj; // BEWARE: order matters in the file
          int t2  = ChargeExcitationIdx[file_line_idx_j][0];  // type
          int rjj = ChargeExcitationIdx[file_line_idx_j][1];
          int rb1 = ChargeExcitationIdx[file_line_idx_j][2];
          int rb2 = ChargeExcitationIdx[file_line_idx_j][3];
          if(rj != rjj) printf("something wrong in excitations? %d == %d \n", rj, rjj);
          
          double O0 = w * ((double) (Commute_Nat_(with_nothing, rb1, rb2, t2, 0,0,s, 0,0,0, myEleNum)));

          //phys_nCAm_averaged[nn + NExcitation*(mm + NExcitation*(ii + Nsite*jj))]; 
          //S_CA[ii + Nsite*(mm + NExcitation*(jj + Nsite*nn))]
          
          int final_idx = ri + Nsite*(idx_exc_mm + NExcitation*(rj + Nsite*idx_exc_nn));
          
          Phys_nACm_averaged[final_idx] += O_AC * O0 ;
          Phys_nCAm_averaged[final_idx] += O_CA * O0 ;

          Phys_nAHCm_averaged[final_idx] += H_AC * O0 ;
          Phys_nCHAm_averaged[final_idx] += H_CA * O0 ;

        }
      }
     }
    }//*/

  } //


  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();  
  return;
}








#endif
