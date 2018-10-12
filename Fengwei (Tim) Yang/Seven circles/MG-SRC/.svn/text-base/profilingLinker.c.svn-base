#include <stdio.h>
#include <stdlib.h>
#ifdef MPE
#include "mpe.h"
#endif

void mpe_init_log_();
void mpe_finish_log_();
#ifdef CRAYCOMP
void coutput_(long int *grids, long int *noprocs, long int *nblocks);
void coutput2_(long int *grids, long int *noprocs, long int *nblocks);
#else
void coutput_(int *grids, int *noprocs, int *nblocks);
void coutput2_(int *grids, int *noprocs, int *nblocks);
#endif


void mpe_init_log_()
{
  printf("src/profilingLinker.c : MPE_Init_log");
#ifdef MPE
  MPE_Init_log();
#endif
}

void mpe_finish_log_()
{
  printf("src/profilingLinker.c : MPE_Finish_log");
#ifdef MPE
  MPE_Finish_log("prof-out-MPE");
#endif
}

#ifdef CRAYCOMP
void coutput_(long int *grids, long int *noprocs, long int *nblocks)
#else
void coutput_(int *grids, int *noprocs, int *nblocks)
#endif
{
  int i, p;
  printf("\ncoutput_\n");
  for (i=0; i<*grids; i++)
  {
  	printf("%d ", i); 
  	for (p=0; p<*noprocs; p++)
/*		printf("%5d (%d)", nblocks[i+2+p*(*grids+2)], i+2+p*(*grids+2));*/
		printf("%5d ", nblocks[i+2+p*(*grids+2)]);
	printf("\n");
  }
  return;
}

#ifdef CRAYCOMP
void coutput2_(long int *grids, long int *noprocs, long int *nblocks)
#else
void coutput2_(int *grids, int *noprocs, int *nblocks)
#endif
{
  int i, p;
  printf("\ncoutput2_\n");
  for (i=0; i<*grids; i++)
  {
  	printf("%d ", i+1); 
  	for (p=0; p<*noprocs; p++)
		printf("%5d ", nblocks[i+p*(*grids)]);
	printf("\n");
  }
  return;
}
