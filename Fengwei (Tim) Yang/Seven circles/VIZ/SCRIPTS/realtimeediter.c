#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX(a,b) ((a)>(b) ? (a):(b) )
#define MIN(a,b) ((a)<(b) ? (a):(b) )


#define AVERAGEOVER 201

void shiftTenset(double *tenset, double time, double loc, double rad);
int testFopen(char *fname, FILE *fptr);

int main(int argc, char **argv)
{
  int i, it;
  char* fileExhausted;
  double old=0.0, new, vold=0.0;
  char fname[128], *lineRadOld, *lineRadNew, *lineLocOld, *lineLocNew, *tmp;
  FILE /* *fptr, *gptr,*/ *hptr, *iptr, *jptr, *kptr;
  double *tenset, oldRad, oldLoc, tmpD, dt1, dt10, dt200, rad10, rad200, vel1, vel10, vel200, avt;

/* Initialise line strings */
  lineRadOld=malloc(1024*sizeof(char));
  lineRadNew=malloc(1024*sizeof(char));
  lineLocOld=malloc(1024*sizeof(char));
  lineLocNew=malloc(1024*sizeof(char));

  if (argc<=1)
  {
	printf("REALSTRIPPER: \nRemoves duplicate timesteps and calculates average speeds and radii");
	printf("Usage:\n        ./REALSTRIPPER 0    : Processes tip_loc.txt and tip_rad.txt in current directory\n");
	printf("        ./REALSTRIPPER n1 [n2 [n3 [...]]]    : Processes tip_loc.<N?> and tip_rad.<N?> in /tmp/\n");
	return 0;
  }


  /* Loop over supplied input files */
  for (i=1; i<argc; i++)
  {

  	/* Open files */
  	if (strcmp(argv[i],"0")==0) 
  	{
		/*sprintf(fname, "real_time.txt");
		fptr=fopen(fname, "r");*/
		sprintf(fname, "tip_loc.txt");
		hptr=fopen(fname, "r");
		if (-1==testFopen(fname, hptr)) continue;
		sprintf(fname, "tip_rad.txt");
		jptr=fopen(fname, "r");
		if (-1==testFopen(fname, jptr)) continue;
  	} else {
  		/*sprintf(fname, "/tmp/real_time.%s", argv[i]);
  		fptr=fopen(fname, "r");*/
  		sprintf(fname, "/tmp/tip_loc.%s", argv[i]);
  		hptr=fopen(fname, "r");
		if (-1==testFopen(fname, hptr)) continue;
  		sprintf(fname, "/tmp/tip_rad.%s", argv[i]);
  		jptr=fopen(fname, "r");
		if (-1==testFopen(fname, jptr)) continue;
  	}
  	/*sprintf(fname, "/tmp/real_time_ed.%s", argv[i]);
  	gptr=fopen(fname, "w");*/
  	sprintf(fname, "/tmp/tip_loc_ed.%s", argv[i]);
  	iptr=fopen(fname, "w");
	if (-1==testFopen(fname, iptr)) continue;
  	sprintf(fname, "/tmp/tip_rad_ed.%s", argv[i]);
  	kptr=fopen(fname, "w");
	if (-1==testFopen(fname, kptr)) continue;

	printf("Processing %s...\n", argv[i]);
  	tenset = calloc(3*AVERAGEOVER,sizeof(double));

  	/* Read first line of files */
  	fileExhausted = fgets(lineRadOld, 1023, jptr);
  	sscanf(lineRadOld, "%le %le", &old, &oldLoc);
  	fileExhausted = fgets(lineLocOld, 1023, hptr);
  	/*printf("Old loc was %f\n", oldLoc);*/
  	for (i=0; i<AVERAGEOVER; i++)
  	{
		tenset[3*i] = 0.0;
		tenset[3*i+1] = oldLoc;
  	}

  	/* Set initial times to dt */
  	it=0;

  	while (fileExhausted != NULL)
  	{
		/* Read next Lines */
	  	fileExhausted = fgets(lineRadNew, 1023, jptr);
	  	sscanf(lineRadNew, "%le", &new);
	  	fileExhausted = fgets(lineLocNew, 1023, hptr);
	  	sscanf(lineLocNew, "%le", &new);
	
		if (new>old) 
		{

	  		sscanf(lineRadOld, "%le %le", &tmpD, &oldRad);
	  		sscanf(lineLocOld, "%le %le", &tmpD, &oldLoc);

			shiftTenset(tenset, old, oldLoc, oldRad);		it++;
	
			/* Write out previous lines */
			lineRadOld[strlen(lineRadOld)-1] = ' ';
			lineLocOld[strlen(lineLocOld)-1] = ' ';

			dt1 = tenset[0]-tenset[3];	dt10 = tenset[0]-tenset[3*10];	dt200 = tenset[0]-tenset[3*200];
		 
			vel1 = (tenset[1]-tenset[4])/dt1;	
			vel10= (tenset[1]-tenset[31])/dt10;	
			vel200=(tenset[1]-tenset[601])/dt200;
			rad10 = 0.0;	rad200=0.0;	avt=0.0;
			for (i=0; i<MIN(10,it); i++) rad10 += tenset[3*i+2];		rad10/=MIN(10,it); 		
			for (i=0; i<MIN(200,it); i++) rad200 += tenset[3*i+2]; 		rad200/=MIN(200,it);
			for (i=0; i<MIN(200,it); i++) avt += tenset[3*i];		avt/=MIN(200,it);
			if (it==1) {vel1=0.0;vel10=0.0;vel200=0.0;}

/*			fprintf(iptr, "%s %e %e %e %e %e\n", lineLocOld, dt1, vel1, vel10, vel200, tenset[300]);*/  /* Last time is time with 100 
											before and after of the 200 averages to match Andy*/

/*			fprintf(kptr, "%s %e %e %e %e\n", lineRadOld, dt1, rad10, rad200, tenset[300]);*/ /* Radius similar */

				fprintf(iptr, "%s %e %e %e %e %e\n", lineLocOld, dt1, vel1, vel10, vel200, avt);  /* Last time is average of the 200  times*/
				fprintf(kptr, "%s %e %e %e %e\n", lineRadOld, dt1, rad10, rad200, avt); /* Radius similar */

			/* Swap pointers, etc, for next valid step */
			vold=old;

			old=new;
		}

		tmp = lineLocOld;
		lineLocOld = lineLocNew;
		lineLocNew = tmp;
	
		tmp = lineRadOld;
		lineRadOld = lineRadNew;
		lineRadNew = tmp;
  	}
  	/* Last line not printed at minute as don't know it was successful */

  	/* Close files */

	/*  fclose(fptr);*/
	/*  fclose(gptr);*/
	  fclose(hptr);
	  fclose(iptr);
	  fclose(jptr);
	  fclose(kptr);


	  free(tenset);
  }
  free (lineLocNew);
  free (lineLocOld);
  free (lineRadNew);
  free (lineRadOld);

  return 0;
}

void shiftTenset(double *tenset, double time, double loc, double rad)
{
  int i;

  for (i=MAX(AVERAGEOVER-2,198); i>=0; i--)
  {	
	tenset[3*(i+1)] = tenset[3*i];
	tenset[3*(i+1)+1] = tenset[3*i+1];
	tenset[3*(i+1)+2] = tenset[3*i+2];
  }
  tenset[0]=time;
  tenset[1]=loc;
  tenset[2]=rad;

/*  for (i=0; i<AVERAGEOVER; i++)
	printf("%e ", tenset[3*i]);
  printf("\n");
  for (i=0; i<AVERAGEOVER; i++)
	printf("%e ", tenset[3*i+1]);
  printf("\n-----------------\n");*/


  return;
}

int testFopen(char *fname, FILE *fptr)
{
  if (fptr==NULL)
  {
	printf("ERROR: File %s failed to open\n", fname);
	return -1; 
  } else {
	return 0;
  }
}
