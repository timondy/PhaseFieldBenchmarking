#include <stdio.h>
#include <stdlib.h>
#include "mpe.h"

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

void coutput_(int *grids, int *noprocs, int *nblocks)
{
  int i, p;
  printf("\ncoutput_\n");
  for (i=0; i<*grids; i++)
  {
  	printf("%d ", i); 
  	for (p=0; p<*noprocs; p++)
		printf("%5d ", nblocks[i+2+p*(*grids+2)]);
	printf("\n");
  }
  return;
}
void coutput2_(int *grids, int *noprocs, int *nblocks)
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
