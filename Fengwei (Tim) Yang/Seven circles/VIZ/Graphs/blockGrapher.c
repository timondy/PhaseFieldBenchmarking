#include <stdio.h>
#include <stdlib.h>

/****************************************************************************
 ***  blockGrapher                                                        ***
 ***                                                                      ***
 ***  Chris Goodyer, May 2012                                             ***
 ***                                                                      ***
 ***  Produces directed graphs from output data of Paramesh applications  ***
 ***  for visualization using GraphViz.                                   ***
 ****************************************************************************/

void digraph_output(int **data, int *N, int NP, int step);

int main(int argc, char **argv)
{
  int *N, NP=-1, **data, i, p, tmp1, tmp2;
  char fname[128], Scase[20], root[32];
  FILE *fptr;
  int args=1, step=-1, process=0;

  /* Read in number of processors */
  sprintf(root,"/tmp/");

  /* Process command line inputs */
  if (argc>0)
  {
   	while (args < argc)
        {
                if (strcmp(argv[args],"--noprocs") == 0 || strcmp(argv[args],"-n") == 0) {
                        /* How many processors were used? */
                       	args++;
                       	if (argv[args])
                               	sprintf(Scase,"%s", argv[args]);
                       	NP=atof(Scase);
                } else if (strcmp(argv[args],"--step") == 0 || strcmp(argv[args],"-s") == 0) {
                        /* Which timestep output to process */
                        args++;
                        if (argv[args])
                                sprintf(Scase,"%s", argv[args]);
                        step=atoi(Scase);
                } else if (strcmp(argv[args],"--root") == 0 || strcmp(argv[args],"-r") == 0) {
                        /* Setting root directory for files */
                        args++;
                        if (argv[args])
                                sprintf(root,"%s", argv[args]);
                } else if (strcmp(argv[args],"--png") == 0 ) {
			process=1;
                } else if (strcmp(argv[args],"--ps") == 0 ) {
			process=2;
                } else if (strcmp(argv[args],"--help") == 0 ) {
                        /* SCC case */
			printf("%s: Turns output data from Phase Field neigh() block structure into a directed graph using GraphViz.\n", argv[0]);
			printf("Command to produce the output: 'neato -Tpng /tmp/digraph.gv > /tmp/digraph.png'\n");
			printf("Usage:\n\t--noprocs <np>|-n <np>      : How many processors were used in the run\n");
			printf(        "\t--step <number>|-s <number> : Which timestep do you want to visualize\n");
			printf(        "\t[--root <dir>|-r <dir>]     : Directory containing files [default /tmp/]\n");
			printf(        "\t[--help|-h]                 : Display this help\n");
			return -1;
		} else {
			printf("WARNING: %s is an unknown argument.  Use --help to see available options.\n", argv[args]);
		}
		args++;
	}
  }
  
  if (NP<1||step<1) 
  {
        printf("\033[31;1;41m");
	printf("ERROR : NP=%d step=%d\n",NP, step);
        printf("\033[0m\n");
	return -1;
  }

  data = malloc(NP*sizeof(int*));
  N = malloc(NP*sizeof(int));

  /* Loop over processors */
  for (p=0; p<NP; p++)
  {
	sprintf(fname, "%s/cgblocks_%.2d_%.5d.%.2d", root, NP, step, p);
  	fptr = fopen(fname,"r");
	if (fptr==NULL) 
	{
	        printf("\033[31;1;41m");
	        printf("Unable to open %s for reading.  Quitting...\n", fname);
	        printf("\033[0m\n");
	        return -1;
	}

  	/* Read in number of points followed by the data */
  	fscanf(fptr, "%d", &N[p]);

  	data[p] = malloc(2*N[p]*6*sizeof(int));
  	for (i=0; i<2*N[p]; i++)
	{
		fscanf(fptr, "%d %d %d %d %d %d", &data[p][6*i+0], &data[p][6*i+1], &data[p][6*i+2], &data[p][6*i+3], &data[p][6*i+4], &data[p][6*i+5]);
		
	}

	fclose(fptr);
  }

  digraph_output(data, N, NP, step);

  switch (process)
  {
	case 1 :
		sprintf(fname, "neato -Tpng /tmp/digraph.gv > /tmp/digraph_%d.png",step);
		system(fname);
		break;
	case 2 :
		sprintf(fname, "neato -Tps /tmp/digraph.gv > /tmp/digraph.ps");
		system(fname);
		break;
	default:
		break;
  }

  return 0;
}

void digraph_output(int **data, int *N, int NP, int step)
{
  FILE *fptr;
  int p, i, colstyle=1;
  char colourstring[20];
  
  fptr = fopen("/tmp/digraph.gv", "w");

  fprintf(fptr, "##Command to produce the output: 'neato -Tpng /tmp/digraph.gv > /tmp/digraph.png'\n");
  fprintf(fptr, "digraph CoarseBlocks {\n\n");
  if (colstyle==1) fprintf(fptr, "node [colorscheme=dark28]\n");

  for (p = 0; p<NP; p++ )
  {
	for (i = 0; i<N[p]; i++ )
	{
  		/* List connections */
		/* Colour red for links between processors */
		/* Fixing lengths means we approximate a 3-d grid structure rather than just a seemingly random graph */
		if (p==data[p][12*i+6]) sprintf(colourstring, "");
			else sprintf(colourstring, ",color=red");
		if (data[p][12*i  ]>0)
			fprintf(fptr, "bk_%d_%d -> bk_%d_%d [len=2.0%s] \n", (i+1),p, data[p][12*i  ],data[p][12*i+6],colourstring);
		if (p==data[p][12*i+7]) sprintf(colourstring, "");
			else sprintf(colourstring, ",color=red");
		if (data[p][12*i+1]>0)
			fprintf(fptr, "bk_%d_%d -> bk_%d_%d [len=2.0%s] \n", (i+1),p, data[p][12*i+1],data[p][12*i+7],colourstring);
		if (p==data[p][12*i+8]) sprintf(colourstring, "");
			else sprintf(colourstring, ",color=red");
		if (data[p][12*i+2]>0)
			fprintf(fptr, "bk_%d_%d -> bk_%d_%d [len=2.5%s] \n", (i+1),p, data[p][12*i+2],data[p][12*i+8],colourstring);
		if (p==data[p][12*i+9]) sprintf(colourstring, "");
			else sprintf(colourstring, ",color=red");
		if (data[p][12*i+3]>0)
			fprintf(fptr, "bk_%d_%d -> bk_%d_%d [len=2.5%s] \n", (i+1),p, data[p][12*i+3],data[p][12*i+9],colourstring);
		if (p==data[p][12*i+10]) sprintf(colourstring, "");
			else sprintf(colourstring, ",color=red");
		if (data[p][12*i+4]>0)
			fprintf(fptr, "bk_%d_%d -> bk_%d_%d [len=3.0,style=dotted%s] \n", (i+1),p, data[p][12*i+4],data[p][12*i+10],colourstring);
		if (p==data[p][12*i+11]) sprintf(colourstring, "");
			else sprintf(colourstring, ",color=red");
		if (data[p][12*i+5]>0)
			fprintf(fptr, "bk_%d_%d -> bk_%d_%d [len=3.0,style=dotted%s] \n", (i+1),p, data[p][12*i+5],data[p][12*i+11],colourstring);

		/* Finally set up the style for this node itself - colour using chosen style by processor number */

		if (colstyle==1)
		{
			/* This option selects from the chosen colorscheme set based on p */
			fprintf(fptr, "bk_%d_%d  [color=%d,style=filled,fontcolor=white,label=\"%d\\n   p%d\"]; \n", (i+1),p,p+1,i+1,p);
		} else {
			/* This loop for manually assigned colours */
		        switch (p) 
			{
				case 0 :
						sprintf(colourstring, "navyblue");
						break;
				case 1 :
						sprintf(colourstring, "tomato");
						break;
				case 2 :
						sprintf(colourstring, "darkgreen");
						break;
				case 3 :
						sprintf(colourstring, "blueviolet");
						break;
				default :
						sprintf(colourstring, "black");
						break;
			}
			fprintf(fptr, "bk_%d_%d  [color=%s,style=filled,fontcolor=white,label=\"%d\\n   p%d\"]; \n", (i+1),p,colourstring,i+1,p);
		}
	}
  }

  /* Set the global options for the graph */

  fprintf(fptr, "overlap=false\n");		/* Prefer lines not to cross */
/*  fprintf(fptr, "pack=true\n");*/
  fprintf(fptr, "rankdir=BT\n");		/* Put the root node at the bottom (generally) */
  fprintf(fptr, "label=\"Top level block distribution\\nTimestep %d, noprocs=%d\"\n", step, NP);	/* Title */
/*  fprintf(fptr, "fontsize=12\n");*/
  fprintf(fptr, "size=\"7.27,11.69\"\n");	/* A4 */
  fprintf(fptr, "}\n");

  fprintf(fptr, "\n");
  fclose(fptr);
  return;
}
