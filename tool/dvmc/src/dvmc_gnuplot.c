/*
zlib license:

Copyright (c) 2017,2019 Maxime Charlebois

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.

A large part of this code has been imported from
onebody.tar.gz

on the webpage:
https://www.physique.usherbrooke.ca/source_code/
*/

#include <stdlib.h>
#include <stdio.h>
//----------------
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <termios.h>
//----------------
#include <sys/stat.h>

#ifdef __APPLE__
    #define _TERMINAL "wxt"
#else
    #define _TERMINAL "wxt"
#endif



void prepareTerminalInputs()
{
    /////////simple trick to get the input from terminal:
    struct termios oldSettings, newSettings;    
    tcgetattr( fileno( stdin ), &oldSettings ); // get settings
    newSettings = oldSettings;
    newSettings.c_lflag &= (~ICANON & ~ECHO);  // turn off the ECHO and ICANON option
    // ECHO: controls whether input is immediately re-echoed as output.
    // ICANON: the terminal buffers a line at a time, and enables line editing. 
    // Turning both off input is made available to programs immediately
    // This is also known as “cbreak” mode.
    tcsetattr( fileno( stdin ), TCSANOW, &newSettings ); // set settings
}   


void positionPlot(FILE *pipe, float lm,float rm,float tm,float bm){
  fprintf (pipe,"set lmargin at screen %4.3f; set rmargin at screen %4.3f;\n", lm,rm);
  fprintf (pipe,"set tmargin at screen %4.3f; set bmargin at screen %4.3f;\n", tm,bm);
}

FILE *gpc_init_image ()
{
	FILE *pipe;
	
  //int error = system("gnuplot tmp.gp");  
  pipe = popen ("gnuplot > /dev/null 2>&1", "w");                  // Open pipe to Gnuplot and check for error
  if (pipe == NULL)
  {
    printf ("Can not find the required Gnuplot executable.\n");
    printf ("\nGraph creation failure\n");
    exit (1);
  }

  fflush (pipe);                                    // flush the pipe
  return (pipe);
}

int gpc_plot_image (FILE *pipe,float *pData, int freq, int Nx, int Ny, float zMin, float zMax, int Nw)
{
  int i, j;
  char title[80];
  
  double chem_y1 = (Ny-0.);
  double chem_y2 = (Ny+0.5);
  double chem_x = (Nx+1.)*freq / ((double) Nw);
  
  //fprintf (pipe, "set arrow 1 from %f,%f to %f,%f  lw 2 lc 'black' front\n",chem_x,chem_y1,chem_x,chem_y2);
  
  fprintf (pipe, "set palette defined (0 'white',0.016 'white', 0.04 '#fff879',0.16 '#fe7e00',0.27 '#d81800', 0.4 '#780000', 0.6 'black')\n");
  fprintf (pipe, "set pm3d map\nset size square\nset cbtics ticspace\nset cbrange [%1.3f:%1.3f]\n", zMin, zMax);
  
  fprintf (pipe, "Nx=%d\nNy=%d\nset xrange [-0.5:Nx-0.5]\nset yrange [-0.5:Ny-0.5]\n",Nx+1,Ny+1);

  //positionPlot(pipe, 0.1,0.9,0.6,0.1);
  fprintf (pipe, "plot '-' matrix with image\n"); // Set plot format

  for (j = 0; j < (Ny+1); j++)                 // For every row
  {
    for (i = 0; i < (Nx+1); i++)               // For every pixel in the row
    {
      fprintf (pipe, "%f ", pData[freq*Nx*Ny+ (((i+Nx/2)%Nx)* Ny) + ((j+Ny/2)%Ny)]);
    }
    fprintf (pipe, "\n");                           // End of isoline scan
  }
  fprintf (pipe, "\ne\ne\n");                       // End of spectrogram dataset
  //fprintf (pipe, "unset arrow 1 \n");
  fflush (pipe);                                    // Flush the pipe
  //fprintf (pipe, "unset multiplot \n");
  return (0);
}

void gpc_close (FILE *pipe)
{
  fprintf (pipe, "exit\n");                         // Close GNUPlot
  pclose (pipe);                                    // Close the pipe to Gnuplot
}




int main(int argc, char* argv[]) {
    if((argc <5) || (argc >6))
    {
        printf("usage example (Nx=8,Ny=8,Nw=1500,NFermi=645,cbar_max=0.6):\n$ vmc_gnuplot 8 8 1500 645 0.6\n");
        exit(0);
    }
    int Nx = atoi(argv[1]);
    int Ny = atoi(argv[2]);
    int Nw = atoi(argv[3]);
    int freq = Nw/2;
    float cbar_max = 0.6;
    if(argc >=5) freq = atoi(argv[4]);
    if(argc >=6) cbar_max = atof(argv[5]);
    
    printf("w=%d    \n",freq);
    
    printf("Nx=%d,Ny=%d,Nw=%d\n",Nx,Ny,Nw);
    printf("First give focus to the command line and second press and hold z or x to change the frame\n");
        
    FILE *hImage;
    
    hImage = gpc_init_image ();
    
    FILE *dataAkw = fopen("output/Akw_all.dat","r");
    if (dataAkw == NULL)
    {
        fprintf(stderr, "Failed to load file output/Akw_all.dat\n");
        exit(1);
    }
    
    int i, j;
    float * dataCube = (float*) malloc(Nx*Ny*Nw*sizeof(float));
    
    for (i = 0; i < Nw; i++)
    {
        for (j = 0; j < Nx*Ny; j++)
        {
            if (fscanf(dataAkw,"%f ", &dataCube[Nx*Ny*i+j]) != 1)
            {
                fprintf(stderr, "Failed to read dataCube at i=%d, j=%d\n", i, j);
                exit(1);
            }
            //printf("%4.5f",dataCube[i+Nw*j]);
        }
    }
    
    
    fclose(dataAkw);
    gpc_plot_image(hImage,dataCube,freq,Nx,Ny,0.0,cbar_max,Nw); 
    
    //exit(0);

    prepareTerminalInputs();
    
    while ( 1 ) // yes, 1
    {
        fd_set set;
        struct timeval tv;
        
        tv.tv_sec = 10;
        tv.tv_usec = 0;
        
        FD_ZERO( &set );
        FD_SET( fileno( stdin ), &set );
        
        int res = select( fileno( stdin )+1, &set, NULL, NULL, &tv );
        if( res > 0 )
        {
            char c;
            char ch[2561];
            if (read( fileno( stdin ), &ch, 2561 )==1){
                c=ch[0];
                int step=1;
                if(c=='x') {
                 if(freq<Nw-step) {
                   freq+=step; printf("            \rw=%d    ",freq); fflush(stdout); gpc_plot_image(hImage,dataCube,freq,Nx,Ny,0.0,cbar_max,Nw);}
                 }
                else if(c=='z') { 
                 if(freq>=step) {
                   freq-=step; printf("            \rw=%d    ",freq); fflush(stdout); gpc_plot_image(hImage,dataCube,freq,Nx,Ny,0.0,cbar_max,Nw);}
                 }
                
                else if(c=='p') {
                  FILE *tmpSlice = fopen("output/tmpSlice.dat","w");
                  for (j = 0; j < (Ny+1); j++)                 // For every row
                  {
                    for (i = 0; i < (Nx+1); i++)               // For every pixel in the row
                    {
                      fprintf (tmpSlice,"%f ", dataCube[freq*Nx*Ny+ (((i+Nx/2)%Nx)* Ny) + ((j+Ny/2)%Ny)]);
                    }
                    fprintf (tmpSlice,"\n");
                  }
                  fclose(tmpSlice);
                }
                if(c==' ') {printf("\n"); fflush(stdout); gpc_plot_image(hImage,dataCube,freq,Nx,Ny,0.0,cbar_max,Nw);  }
            }
        }
        else if( res < 0 )
        {
            perror( "select error\n" );
            break;
        }
    }
    
    //next two lines never used because we ctrl-c out of the program.
    //tcsetattr( fileno( stdin ), TCSANOW, &oldSettings );
    //gpc_close (hImage);
    //free(dataCube);
    exit(0);
}



