#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "io.h"

//====================================================================
//
//    Program lower_resolution:     lower the resoltion of an N-body 
//        initial conditions file by a factor of 8 (2 per dimension).
//        It only works for single files and dark matter. Indexes are
//        assumed to be correct, if not the particle's position will 
//        be messed up.
//
//
//    Compile:
//       gcc -lm lower_resolution.c -o lower_resolution
//
//    Comments:
//       file io routines based on read_snap.c,and io.c (Volker).
//
//    Written by:
//       Miguel Angel Aragon Calvo (miguel@astro.rug.nl) 

//
//====================================================================


struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  //--- Fills to 256 Bytes 
} header1, header2;


//--- Global variables
int     NumPart, NumPartLow, Ngas;
int     n_low, n_high;
int     *Id, *Id_low;                                               //--- Pointer to particle's ID
double  Time, Redshift;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P1, *P2;                                                //--- Pointer to structure



//====================================================
//================  MAIN FINCTION  ===================
//====================================================
int main(int argc, char **argv)
{
  char path[200], input_fname[200], basename[200];
  int  type, snapshot_number, files;


  if (argc != 5)
    {
      printf("<-----------------------------------------------------------> \n");
      printf("Usage: \n");
      printf("       lower_resolution SNAP128_000 128 SNAP064_000 64 \n");
      printf("<-----------------------------------------------------------> \n");
      exit(0);
    }

  //--- Number of files per snapshot
  files=1;                               

  //--- TEMPORAL
  n_high = atoi(argv[2]);
  n_low  = atoi(argv[4]);

  printf("Reading particles from file...\n");
  load_snapshot(argv[1], files);
  
  printf("Merge particles\n");
  merge_particles();

  printf("Write particles to file\n");
  savepositions_ioformat1(argv[3]);

}



//===============================================================
//        Here the particle data is at your disposal
//===============================================================
int merge_particles(void)
{
  int i, j, k, cont;
  int s1, s2, s3, s4, s5, s6, s7, s8;
  int ind; 

  
  //--- Address particle array with offset 1 !!!
  cont   = 1;
  for(i=0; i<=n_low-1; i++)
    {
      for(j=0; j<=n_low-1; j++)
	{
	  for(k=0; k<=n_low-1; k++)
	    {
	      
	      ind = 2*i + 2*n_high*j + 2*n_high*n_high*k;
	      
	      s1 =  ind+1;
	      s2 =  ind+1+1;
	      s3 =  ind+n_high+1;
	      s4 =  ind+n_high+1+1;
	      s5 =  ind+n_high*n_high+1;
	      s6 =  ind+n_high*n_high+1+1;
	      s7 =  ind+n_high*n_high+n_high+1;
	      s8 =  ind+n_high*n_high+n_high+1+1;
	     
	      //--- Average position
	      P2[cont].Pos[0] = (P1[s1].Pos[0]+P1[s2].Pos[0]+P1[s3].Pos[0]+P1[s4].Pos[0]+P1[s5].Pos[0]+P1[s6].Pos[0]+P1[s7].Pos[0]+P1[s8].Pos[0])/8.0;
	      P2[cont].Pos[1] = (P1[s1].Pos[1]+P1[s2].Pos[1]+P1[s3].Pos[1]+P1[s4].Pos[1]+P1[s5].Pos[1]+P1[s6].Pos[1]+P1[s7].Pos[1]+P1[s8].Pos[1])/8.0;
	      P2[cont].Pos[2] = (P1[s1].Pos[2]+P1[s2].Pos[2]+P1[s3].Pos[2]+P1[s4].Pos[2]+P1[s5].Pos[2]+P1[s6].Pos[2]+P1[s7].Pos[2]+P1[s8].Pos[2])/8.0;
	      
	      //--- Average velocities
	      P2[cont].Vel[0] = (P1[s1].Vel[0]+P1[s2].Vel[0]+P1[s3].Vel[0]+P1[s4].Vel[0]+P1[s5].Vel[0]+P1[s6].Vel[0]+P1[s7].Vel[0]+P1[s8].Vel[0])/8.0;
	      P2[cont].Vel[1] = (P1[s1].Vel[1]+P1[s2].Vel[1]+P1[s3].Vel[1]+P1[s4].Vel[1]+P1[s5].Vel[1]+P1[s6].Vel[1]+P1[s7].Vel[1]+P1[s8].Vel[1])/8.0;
	      P2[cont].Vel[2] = (P1[s1].Vel[2]+P1[s2].Vel[2]+P1[s3].Vel[2]+P1[s4].Vel[2]+P1[s5].Vel[2]+P1[s6].Vel[2]+P1[s7].Vel[2]+P1[s8].Vel[2])/8.0;

	      cont++;
	      
	    }
	}
    }
}











