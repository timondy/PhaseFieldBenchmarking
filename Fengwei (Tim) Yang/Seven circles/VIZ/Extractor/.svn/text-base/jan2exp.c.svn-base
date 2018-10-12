#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  long int N, i;
  double *data, *coords, d1, d2;
  char fname[16];
  FILE *fptr, *lattice, *pyr;
  
  if (argc==1) {printf("Jan2Exp usage : ./Jan2Exp <filenumber>\n"); return -1;}
  sprintf(fname, "/tmp/%s.dat", argv[1]);
  fptr = fopen(fname, "r");

  fscanf(fptr, "%lf %lf %ld %ld", &d1, &d2, &i, &N);

  printf("N=%ld\n", N);

  coords = malloc(2*N*sizeof(double));
  data = malloc(2*N*sizeof(double));

  for (i=0; i<N; i++)
	fscanf(fptr, "%lf %lf %lf %lf", &coords[2*i], &coords[2*i+1], &data[2*i], &data[2*i+1]);

  fclose(fptr);

  sprintf(fname, "/tmp/%snd.lat", argv[1]);
  lattice = fopen(fname, "w");

  fprintf(lattice,"#!/usr/explorer/bin/explorer cxLattice plain 1.0\n");   /* HEADER */
  fprintf(lattice," 1\n");                                                 /* DIMENSIONS */
  fprintf(lattice," %d\n", N);                                 /* NUMBER OF DATA POINTS PER DIMENSION*/
  fprintf(lattice," 2\n");                                                 /* NUMBER OF DATA VALUES PER NODE */
  fprintf(lattice," 4\n");                                                 /* DATA FORMAT: 3=float, 4=double */
  fprintf(lattice," 2\n");                                                 /* COORDINATE TYPE: 2=curvilinear */
  fprintf(lattice," 1\n");                                                 /* NUMBER OF DATA SETS PRESENT */
  fprintf(lattice," 2\n");                                                 /* COORDINATE DIMENSIONS */
 
  for (i=0; i<N; i++)
	fprintf(lattice, "%f %f\n", coords[2*i], coords[2*i+1]); 

  for (i=0; i<N; i++)
	fprintf(lattice, "%f %f\n", data[2*i], data[2*i+1]); 

  fclose(lattice);

#ifdef PYRAMID
  pyr = fopen("/tmp/1000.pyr", "w");

  fprintf(pyr,"#!/usr/explorer/bin/explorer cxPyramid plain 1.0\n");   /* HEADER */
  fprintf(pyr,"include %s\n", "/tmp/1000nd.lat");                               /* INCLUDE LATTICE COORDINATES AND DATA */
  fprintf(pyr," 1\n");                                                 /* NUMBER OF PYRAMID LAYERS */
  fprintf(pyr," 1\n");                                                 /* COMPRESSION TYPE: 1=unique */
  fprintf(pyr," 0\n");                                                 /* COMPRESSION TYPE: 4=tetrahedra */
  fprintf(pyr," 0 0\n");                                               /* LAYER 1: COMPRESSED */
  fprintf(pyr," 0 0\n");                                               /* LAYER 2: COMPRESSED */
  fprintf(pyr," 0 0\n");                                               /* LAYER 3: # TETS, # CONNECTIONS */

  fclose(pyr);
#endif

  return 0;
}
