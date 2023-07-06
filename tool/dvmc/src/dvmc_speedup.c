/*
*/
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>

//
void interrupt_handler(int sig) {exit(0);}
//


//simpler, faster than function above (require binary output (*.bin)
int mergeOutputBin(int n_file, const char ** string_list, int *n_exc, int Nsite, int dummy, 
                   float *S_CA, float *S_AC, float *H_CA, float *H_AC, int verbose) {
           
  // to be able to use ctrl-c
  signal(SIGINT, interrupt_handler); 
  //
  
  double factor = 1./n_file;
  int ii;
  
  for(ii = 0; ii<n_file; ii++){
    printf("%s\n",string_list[ii]);
  }

  FILE *fp;
  fp = fopen(string_list[0], "rb");
  if (fp == NULL) {
    fprintf(stdout, "error: no '%s' found.\n",string_list[ii]); 
    fclose(fp);
    exit(-1);
  }
  
  int NExcitation, Excitation_L, Excitation_W;
  if(fread(&NExcitation,  sizeof(int), 1, fp) ==0){
    fprintf(stdout, "error: reading NExcitation in file '%s'.\n",string_list[ii]); 
    exit(-1);
  }
  if(fread(&Excitation_L,  sizeof(int), 1, fp) ==0){
    fprintf(stdout, "error: reading Excitation_L in file '%s'.\n",string_list[ii]); 
    exit(-1);
  }
  if(fread(&Excitation_W,  sizeof(int), 1, fp) ==0){
    fprintf(stdout, "error: reading Excitation_W in file '%s'.\n",string_list[ii]); 
    exit(-1);
  }
  
  fclose(fp);
  
  //int Nsite = Excitation_L*Excitation_W;
  long int jj,size = NExcitation*NExcitation*Nsite*Nsite;
  float * data_read  = (float *)calloc((4*size),sizeof(float));
  
  //fread(data_read, sizeof(double), 4*size, fp);
  
  if(verbose) fprintf(stdout, "\n\n%d sites to read per file\n\n",Nsite);
  for(ii = 0; ii<n_file; ii++){
    if(verbose) fprintf(stdout, "reading file '%s'\n",string_list[ii]);
    fp = fopen(string_list[ii], "rb");
    if (fp == NULL) {
      fprintf(stdout, "error: no '%s' found.\n",string_list[ii]); 
      fclose(fp);
      exit(-1);
    }
    
    int NExcitation_tmp, Excitation_L_tmp, Excitation_W_tmp;
    if(fread(&NExcitation_tmp,  sizeof(int), 1, fp) ==0){
      fprintf(stdout, "error: reading NExcitation in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }
    if(fread(&Excitation_L_tmp,  sizeof(int), 1, fp) ==0){
      fprintf(stdout, "error: reading Excitation_L in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }
    if(fread(&Excitation_W_tmp,  sizeof(int), 1, fp) ==0){
      fprintf(stdout, "error: reading Excitation_W in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }
    
    if(NExcitation_tmp!=NExcitation || Excitation_L_tmp!=Excitation_L || Excitation_W_tmp!=Excitation_W){
      fprintf(stdout, "error: incoherent header in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }

    if(fread(data_read, sizeof(float), 4*size, fp) ==0){
      fprintf(stdout, "error: reading Excitation_W in file '%s'.\n",string_list[ii]); 
      exit(-1);
    }

    for(jj=0;jj<size;jj++){
      S_CA[jj] += factor*data_read[jj+size*0]; 
      S_AC[jj] += factor*data_read[jj+size*1];
      H_CA[jj] += factor*data_read[jj+size*2]; 
      H_AC[jj] += factor*data_read[jj+size*3]; 
      }
    fclose(fp);    
  }

  free(data_read);

  if(verbose) {
    printf("\nData read:\n");
    for(ii=0;ii<10;ii++) printf("% 4.5f % 4.5f % 4.5f % 4.5f \n",S_CA[ii], S_AC[ii], H_CA[ii], H_AC[ii]);
    printf("...\n");
  }
  
  return 0;
}



int greenFrom_e_U_Uinv_S(int sign, int Nw, double *w, int dim, double complex *e, double complex *U, double complex *Uinv_S, double Omega, double mu, double eta, double complex *g) //)H, double *S,)
{
  int ii,nn;
   
  for(ii=0;ii<Nw;ii++) {
    g[ii]=0.0;
    double complex z = w[ii] + I*eta + sign*Omega + mu;
    for(nn=0;nn<dim;nn++){
      g[ii] += U[nn]*Uinv_S[dim*nn] /(z-sign*e[nn]);
    }
  }

  return 0;
}

int greenFrom_e_U_Uinv_S_general(int sign, int ll, int pp, int Ns, int Nw, double *w, int Nex, double complex *e, double complex *U, double complex *Uinv_S, double Omega, double mu, double eta, double complex *g) //)H, double *S,)
{
  int ii,nn,ll2,pp2;
  int ind1, ind2;

  for(ii=0;ii<Ns*Ns*Nw;ii++) g[ii]=0.0;

  for(ii=0;ii<Nw;ii++) {
    double complex z = w[ii] + I*eta + sign*Omega + mu;

    for(pp2=0;pp2<Ns;pp2++){
      for(nn=0;nn<Nex;nn++){
        ind1 = pp2 + Ns*(nn+Nex*(ll+Ns*0));
        ind2 = pp + Ns*(0+Nex*(pp2+Ns*nn));
        g[ii] += U[ind1]*Uinv_S[ind2]/(z-sign*e[nn]);
	    }
	  }
  }
  return 0;
}

int fourier_Green(int Lx, int Ly, double complex *gr, double complex *gk)
{

  int ii, jj, mm, nn, ll, pp;
  int Ns = Lx*Ly;
  int N[2] ={Lx,Ly} ;

  double kx[Ns], ky[Ns], x[Ns], y[Ns];
  double pi = acos(-1.0);

  int ntemp, den;

  // get coordinates for lattice sites
  for(ii=0;ii<Ns;ii++){
    ntemp=ii;
    for(jj=1;jj>=0;jj--){
      den=1;
      for(ll=0;ll<jj;ll++){
	den=den*N[ll];
      }
      if(jj==0) {
	x[ii]=ntemp/den;
	ntemp=ntemp-x[ii]*den;
	x[ii]=x[ii];//+1;
      }
      if(jj==1) {
	y[ii]=ntemp/den;
	ntemp=ntemp-y[ii]*den;
	y[ii]=y[ii];//+1;
      }
    }
    //printf("%d\t%f\t%f\n",ii,x[ii],y[ii]);
  }

  for(ii=0;ii<Ns;ii++){
    kx[ii] = (2*pi*x[ii])/Lx;
    ky[ii] = (2*pi*y[ii])/Ly;
    //printf("%d\t%f\t%f\t%f\t%f\n",ii,(2*pi*x[ii]),(2*pi*y[ii]),kx[ii],ky[ii]);
  }

  for(ii=0;ii<Ns*Ns;ii++) gk[ii] = 0.0;

  /*
  for(ii=00;ii<Ns;ii++){
    for(jj=00;jj<Ns;jj++){
      printf("%d\t%d\t%f\t%f\n",ii,jj,creal(gr[Ns*ii+jj]),cimag(gr[Ns*ii+jj]));
    }
  }
  */

  // fourier transform G_{k,kp} = 1/Ns * sum_{ij} e^ik.r_i*e^ikp.r_j * G_{ij}
  for(mm=0;mm<Ns;mm++){
    for(nn=0;nn<Ns;nn++){
      for(ii=0;ii<Ns;ii++){
	for(jj=0;jj<Ns;jj++){
	  double complex phase = cexp(I*(kx[mm]*x[ii]+ky[mm]*y[ii]))*cexp(-1.0*I*(kx[nn]*x[jj]+ky[nn]*y[jj]));
	  //double complex phase = I*(kx[mm]*x[ii]+ky[mm]*y[ii]);
	  //if(mm==nn) printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",kx[mm],ky[mm],x[ii],y[ii],creal(phase),cimag(phase),creal(gr[Ns*ii+jj]),cimag(gr[Ns*ii+jj]));
	  gk[Ns*mm+nn] += (1.0/Ns)*phase*gr[Ns*ii+jj];
	}
      }
      //if(mm==nn) 
      //printf("%d\t%d\t%d\t%f\t%f\n",mm,nn,Ns*mm+nn,creal(gk[Ns*mm+nn]),cimag(gk[Ns*mm+nn]));
    }
  }

  return 0;

}

