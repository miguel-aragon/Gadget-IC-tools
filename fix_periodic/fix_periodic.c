#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "io.h"

//====================================================================
//
//    Program fix_periodic
//       Fix the periodic boundary conditions of a gadget file.
//       Zeldovich displacementes can move the particles out
//       of the box
//
//
//    Comments:
//       file io routines based on read_snap.c,and io.c (Volker).
//
//    Written by:
//       Miguel Angel Aragon Calvo (miguel@pha.jhu.edu) 
//
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
} header1;


//--- Global variables
int     NumPart, NumPartLow, Ngas;
float   n_high;
int     *Id;                                               //--- Pointer to particle's ID
double  Time, Redshift;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P1;                                                //--- Pointer to structure


//====================================================
//================  MAIN FINCTION  ===================
//====================================================
int main(int argc, char **argv)
{
  char path[200], input_fname[200], basename[200];
  int  type, snapshot_number, files;


  if (argc != 3)
    {
      printf("<-----------------------------------------------------------> \n");
      printf("Usage: \n");
      printf("       fix_periodic FILE_NO_PERIODIC.GAD FILE_PERIODIC.GAD \n");
      printf("<-----------------------------------------------------------> \n");
      exit(0);
    }

  //--- Number of files per snapshot
  files=1;                               

  printf("Reading particles from 64 bit fortran file...\n");
  load_snapshot(argv[1], files);
  
  printf("Fixing periodic boundaries...\n");
  fix_periodic();

  printf("Writing particles to   32 bit fortran file\n");
  savepositions_ioformat1(argv[2]);

}









