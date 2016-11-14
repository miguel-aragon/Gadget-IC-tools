#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "io.h"
#include "fix_periodic.h"

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
int     NumPart, Ngas;
float   n_high;
int     *Id;                                         //--- Pointer to particle's ID
double  Time, Redshift;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P1;                                            //--- Pointer to structure


//=================================================================
/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
//=================================================================
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k;
  int dummy;
  int    ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file '%s'\n",buf);
	  exit(0);
	}

      printf("Reading file '%s' ...\n",buf); 

      //--- Read header
      fread(&dummy, sizeof(dummy), 1, fd);     
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      //--- Number of particles
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	NumPart+= header1.npart[k];
      
      //--- Gas particles
      Ngas= header1.npart[0];

      //--- Particles with mass
      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      //--- Allocate memory for particles
      if(i==0)
	allocate_memory();

      
      printf("   Read positions...\n"); 
      //--- Read Particle's postitions (START AAT 1!!!)
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P1[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      printf("   Read velocities...\n");
      //--- Read Particle's velocities (START AAT 1!!!)
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P1[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
    
      printf("   Read ID's...\n");
      //--- Read Particle's ID (START AT 1!!!)
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      //--- Read masses (if specified)
      if(ntot_withmasses>0)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P1[pc_new].Type=k;

	      if(header1.mass[k]==0)
		fread(&P1[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P1[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
      

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P1[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P1[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P1[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P1[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }

  printf("Ready reading particles...\n");

  Time= header1.time;
  Redshift= header1.time;
}




//=================================================================
/* this routine allocates the memory for the 
 * particle data.
 */
//=================================================================
int allocate_memory(void)
{
  printf("allocating memory...\n");

  //--- High resolution structure
  if(!(P1=malloc(NumPart*sizeof(struct particle_data) )))
    {
      fprintf(stderr,"failed to allocate memory HIGH RES.\n");
      exit(0);
    }
  
  P1--;  //--- start with offset 1

  //--- High resolution ID
  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;       //--- start with offset 1 

  printf("allocating memory...done\n");
}


//=================================================================
//
//=================================================================
void fix_periodic(void){

  int   n;
  float BoxSize = header1.BoxSize;

  printf("Fixing %d particles inside %f", NumPart, BoxSize);

  //--- Loop over particles and fix the ones that are out of the box
  for(n=1;n<=NumPart;n++)
    {

      //--- X
      if (P1[n].Pos[0] > BoxSize) P1[n].Pos[0] = P1[n].Pos[0]-BoxSize;
      if (P1[n].Pos[0] <    0   ) P1[n].Pos[0] = P1[n].Pos[0]+BoxSize;

      //--- Y
      if (P1[n].Pos[1] > BoxSize) P1[n].Pos[1] = P1[n].Pos[1]-BoxSize;
      if (P1[n].Pos[1] <    0   ) P1[n].Pos[1] = P1[n].Pos[1]+BoxSize;

      //--- Z
      if (P1[n].Pos[2] > BoxSize) P1[n].Pos[2] = P1[n].Pos[2]-BoxSize;
      if (P1[n].Pos[2] <    0   ) P1[n].Pos[2] = P1[n].Pos[2]+BoxSize;

    }

}


//=================================================================
//
//=================================================================
int savepositions_ioformat1(char *fname)
{
  FILE *fd_out;
  char buf[100];
  float dummy[3];
  int i,k;
  int   blklen,masscount;

#define BLKLEN my_fwrite(&blklen, sizeof(blklen), 1, fd_out);

   //--- Output file
  sprintf(buf,"%s",fname);

  if((fd_out=fopen(buf,"w")))
    {

      printf("Start writting file '%s' ...\n",buf);

      //--- Write Header in one pass
      blklen=sizeof(header1);
      BLKLEN;
      my_fwrite(&header1, sizeof(header1), 1, fd_out);
      BLKLEN;

     
      blklen=NumPart*3*sizeof(float);   //--- FORTRAN record (particle's array)
      //--- Write Postitons
      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P1[i].Pos[k];
	  my_fwrite(dummy,sizeof(float),3,fd_out);
	}
      BLKLEN;


      //--- Write Velocities     
      BLKLEN;
      for(i=1;i<=NumPart;i++)
	{
	  for(k=0;k<3;k++)
	    dummy[k]=P1[i].Vel[k];
	  my_fwrite(dummy,sizeof(float),3,fd_out);
	}
      BLKLEN;
      
      
      //--- Write ID's
      blklen=NumPart*sizeof(int);   //--- FORTRAN record (particle's array)
      BLKLEN;
      for(i=0;i<=NumPart-1;i++)
      {
	my_fwrite(&i,sizeof(int),1,fd_out);

      }
      BLKLEN;

      fclose(fd_out);
      printf("Done  writting file '%s'\n",buf);

    }
  else
    {
      fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
    }
}



//------------------------------------------------------
//                     my_fwrite
//------------------------------------------------------
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

  if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("I/O error (fwrite) on has occured.\n");
      fflush(stdout);
    }
  return nwritten;
}


//------------------------------------------------------
//                     my_fread
//------------------------------------------------------
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;

  if((nread=fread(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("I/O error (fread) has occured.\n");
      fflush(stdout);
    }
  return nread;
}
