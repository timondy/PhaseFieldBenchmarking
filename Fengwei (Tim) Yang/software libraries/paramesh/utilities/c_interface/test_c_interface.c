/*----------------------------------------------------------------------
! Paramesh - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------*/
#include "underscore.h"
/*#define LAHEY*/

#ifdef REAL8
#define float double
#endif

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "math.h"

#ifdef UNDERSCORE
#define c_amr_initialize c_amr_initialize_
#define c_amr_1blk_guardcell_reset c_amr_1blk_guardcell_reset_
#define c_amr_refine_derefine c_amr_refine_derefine_
#define c_mpi_morton_bnd_prolong c_mpi_morton_bnd_prolong_
#define c_amr_prolong c_amr_prolong_
#define c_amr_flux_conserve c_amr_flux_conserve_
#define c_amr_guardcell c_amr_guardcell_
#define c_amr_restrict c_amr_restrict_
#define c_amr_1blk_copy_soln c_amr_1blk_copy_soln_
#define c_mpi_amr_comm_setup c_mpi_amr_comm_setup_
#define c_amr_1blk_guardcell c_amr_1blk_guardcell_
#define c_amr_1blk_to_perm c_amr_1blk_to_perm_
#define c_amr_plotfile_chombo c_amr_plotfile_chombo_
#define c_amr_checkpoint_re c_amr_checkpoint_re_
#define c_amr_checkpoint_wr c_amr_checkpoint_wr_
#define c_amr_close c_amr_close_
#endif

#ifdef DOUBLE_UNDERSCORE
#define c_amr_initialize c_amr_initialize__
#define c_amr_1blk_guardcell_reset c_amr_1blk_guardcell_reset__
#define c_amr_refine_derefine c_amr_refine_derefine__
#define c_mpi_morton_bnd_prolong c_mpi_morton_bnd_prolong__
#define c_amr_prolong c_amr_prolong__
#define c_amr_flux_conserve c_amr_flux_conserve__
#define c_amr_guardcell c_amr_guardcell__
#define c_amr_restrict c_amr_restrict__
#define c_amr_1blk_copy_soln c_amr_1blk_copy_soln__
#define c_mpi_amr_comm_setup c_mpi_amr_comm_setup__
#define c_amr_1blk_guardcell c_amr_1blk_guardcell__
#define c_amr_1blk_to_perm c_amr_1blk_to_perm__
#define c_amr_plotfile_chombo c_amr_plotfile_chombo__
#define c_amr_checkpoint_re c_amr_checkpoint_re__
#define c_amr_checkpoint_wr c_amr_checkpoint_wr__
#define c_amr_close c_amr_close__
#endif

extern void c_amr_initialize(void);
extern void c_amr_1blk_guardcell_reset(void);
extern void c_amr_refine_derefine(void);
extern void c_mpi_morton_bnd_prolong(int*, int*, int*);
extern void c_amr_prolong(int*, int*, int*);
extern void c_amr_flux_conserve(int*, int*, int*);
extern void c_amr_guardcell(int*, int*, int*, int*, int*, int*);
extern void c_amr_restrict(int*, int*, int*, int*);
extern void c_amr_1blk_copy_soln(int*);
extern void c_mpi_amr_comm_setup(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern void c_amr_1blk_guardcell(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern void c_amr_1blk_to_perm(int*, int*, int*, int*, int*, int*, int*);
extern void c_amr_checkpoint_wr(int*,int*,char*,float *,float *,float *,float *,float *);
extern void c_amr_checkpoint_re(int*,int*,char*,float *,float *,float *,float *,float *);
extern void c_amr_plotfile_chombo(int*);
extern void c_amr_close(void);

#ifdef LAHEY
MAIN__(iargc, iargv) 
#else
main(iargc, iargv) 
#endif
     int iargc;
     char** iargv;
{

  FILE* test_log;
  char file_name[80];
  int mype; int i; int ii; int iii; int iiii; int j; int k; int lb; 
  int nnvar; int ivar; int loop_count; int iunit1; int iopt; int nlayers;
  int tag_offset; int nprocs;
  int lguard; int lprolong; int lflux; int ledge; int lrestrict; int lfulltree;
  int lcc; int lfc; int lec; int lnc; 
  int nlayersx; int nlayersy; int nlayersz;
  int flux_dir; int l_srl_only; int ldiag; int icoord; int idest;
  int nguard0; int nguard_work0;
  float dx; float dy; float dz; float x0; float y0; float z0;
  float ax; float ay; float az; float value;
  float eps;
  int total_errors;
  int nsub;

  /* initialize paramesh */

  MPI_Init(&iargc, &iargv);
  c_amr_initialize();
  eps = 1.e10;

  /* NOTE THESE INCLUDES MUST GO AFTER c_amr_initialize and they must be      in this order !!!*/
#include "paramesh_dimensions.h"
#include "tree.h"
#include "physicaldata.h"
#include "workspace.h"
#include "io.h"

  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  total_errors = 0;

  nguard0 = *nguard**npgs;
  nguard_work0 = *nguard_work**npgs;

  ax = 1.;
  ay = 10.;
  az = 100.;

  printf("mdim %d\n",*mdim);
  printf("mfaces, nfaces = %d %d\n",*mfaces,*nfaces);
  printf("mchild, nchild = %d %d\n",*mchild,*nchild);
  printf("nguard, npgs, nguard0 = %d %d %d\n",*nguard,*npgs,nguard0);


  /* initialize tree data structure */

  for (lb = 0; lb < *maxblocks_tr; lb++) {
    for (i = 0; i < *mdim; i++) {
      (*bsize)[lb][i] = -1.;
    }
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*lrefine)[lb] = -1;
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*nodetype)[lb] = -1;
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*stay)[lb] = 0;
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*refine)[lb] = 0;
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*derefine)[lb] = 0;
  }


  for (lb = 0; lb < *maxblocks_tr; lb++) {
    for (i = 0; i < 2; i++) {
      (*parent)[lb][i] = -1;
    }
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    for (j = 0; j < *nchild; j++) {
      for (i = 0; i < 2; i++) {
	(*child)[lb][j][i] = -1;
	j = j + 1;
      }
    }
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*which_child)[lb] = -1;
  }
  for (lb = 0; lb < *maxblocks; lb++) {
    for (i = 0; i < *ndim; i++) {
      (*coord)[lb][i] = -1.;
      j = j + 1;
    }
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    for (j = 0; j < *ndim; j++) {
      for (i = 0; i < 2; i++) {
	(*bnd_box)[lb][j][i] = -1.;
	j = j + 1;
      }
    }
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    for (j = 0; j < *mfaces; j++) {
      for (i = 0; i < 2; i++) {
	(*neigh)[lb][j][i] = -1;
      }
    }
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*empty)[lb] = -1;
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    for (i = 0; i < *mflags; i++) {
      (*bflags)[lb][i] = -1;
    }
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    (*work_block)[lb] = 0.;
  }
  for (lb = 0; lb < *maxblocks_tr; lb++) {
    for (iiii = 0; iiii < 1+ (2* *k3d); iiii++) {
      for (iii = 0; iii < 1+ (2* *k2d); iii++) {
	for (ii = 0; ii < 3; ii++) {
	  (*surr_blks)[lb][iiii][iii][ii] = -1;
	}
      }
    }
  }

  if (*nvar > 0) {

  for (lb = 0; lb < *maxblocks; lb++) {
    for (k = 0; k < *ku_bnd; k++) {
      for (j = 0; j < *ju_bnd; j++) {
	for (i = 0; i < *iu_bnd; i++) {
	  for (nnvar = 0; nnvar < *nvar; nnvar++) {
	    (*unk)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }

  }  /* end if (*nvar) */

  if (*nvarcorn > 0) {

  for (lb = 0; lb < *maxblocksn; lb++) {
    for (k = 0; k < *ku_bnd + *k3d; k++) {
      for (j = 0; j < *ju_bnd + *k2d; j++) {
	for (i = 0; i < *iu_bnd + 1; i++) {
	  for (nnvar = 0; nnvar < *nbndvarc; nnvar++) {
	    (*unk_n)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }

  }  /* end if (*nvarcorn) */

  if (*nfacevar > 0) {

  for (lb = 0; lb < *maxblocksf; lb++) {
    for (k = 0; k < *ku_bnd; k++) {
      for (j = 0; j < *ju_bnd; j++) {
	for (i = 0; i < *iu_bnd+1; i++) {
	  for (nnvar = 0; nnvar < *nbndvar; nnvar++) {
	    (*facevarx)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksf; lb++) {
    for (k = 0; k < *ku_bnd; k++) {
      for (j = 0; j < *ju_bnd+*k2d; j++) {
	for (i = 0; i < *iu_bnd; i++) {
	  for (nnvar = 0; nnvar < *nbndvar; nnvar++) {
	    (*facevary)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksf; lb++) {
    for (k = 0; k < *ku_bnd+*k3d; k++) {
      for (j = 0; j < *ju_bnd; j++) {
	for (i = 0; i < *iu_bnd; i++) {
	  for (nnvar = 0; nnvar < *nbndvar; nnvar++) {
	    (*facevarz)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }

  }   /* end if (*nfacvar > 0) */

  for (lb = 0; lb < *maxblocksfl; lb++) {
    for (k = 0; k < *ku_bndi-nguard0**k3d; k++) {
      for (j = 0;  j < *ju_bndi-nguard0; j++) {
	for (i = 0; i < 2; i++) {
	  for (nnvar = 0; nnvar < *nfluxes; nnvar++) {
	    (*flux_x)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksfl; lb++) {
    for (k = 0;  k < *ku_bndi-nguard0**k3d;  k++) {
      for (j = 0; j < 2; j++) {
	for (i = 0; i < *iu_bndi-nguard0; i++) {
	  for (nnvar = 0; nnvar < *nfluxes; nnvar++) {
	    (*flux_y)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksfl; lb++) {
    for (k = 0; k < 2; k++) {
      for (j = 0; j < *ju_bndi-nguard0; j++) {
	for (i = 0; i < *iu_bndi-nguard0; i++) {
	  for (nnvar = 0; nnvar < *nfluxes; nnvar++) {
	    (*flux_z)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }

  if (*nvaredge > 0) {

  for (lb = 0; lb < *maxblocksue; lb++) {
    for (k = 0; k < *ku_bnd+*k3d; k++) {
      for (j = 0; j < *ju_bnd+*k2d; j++) {
	for (i = 0; i < *iu_bnd; i++) {
	  for (nnvar = 0; nnvar < *nbndvare; nnvar++) {
	    (*unk_e_x)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksue; lb++) {
    for (k = 0; k < *ku_bnd+*k3d; k++) {
      for (j = 0; j < *ju_bnd; j++) {
	for (i = 0; i < *iu_bnd+1; i++) {
	  for (nnvar = 0; nnvar < *nbndvare; nnvar++) {
	    (*unk_e_y)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksue; lb++) {
    for (k = 0; k < *ku_bnd; k++) {
      for (j = 0; j < *ju_bnd+*k2d; j++) {
	for (i = 0; i < *iu_bnd+1; i++) {
	  for (nnvar = 0; nnvar < *nbndvare; nnvar++) {
	    (*unk_e_z)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }

  }  /* end if (*nvaredge) */

  for (lb = 0; lb < *maxblocksfl; lb++) {
    for (k = 0; k < *ku_bndi-nguard0**k3d; k++) {
      for (j = 0;  j < *ju_bndi-nguard0; j++) {
	for (i = 0; i < 2; i++) {
	  for (nnvar = 0; nnvar < *nfluxes; nnvar++) {
	    (*flux_x)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksfl; lb++) {
    for (k = 0;  k < *ku_bndi-nguard0**k3d;  k++) {
      for (j = 0; j < 2; j++) {
	for (i = 0; i < *iu_bndi-nguard0; i++) {
	  for (nnvar = 0; nnvar < *nfluxes; nnvar++) {
	    (*flux_y)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }
  for (lb = 0; lb < *maxblocksfl; lb++) {
    for (k = 0; k < 2; k++) {
      for (j = 0; j < *ju_bndi-nguard0; j++) {
	for (i = 0; i < *iu_bndi-nguard0; i++) {
	  for (nnvar = 0; nnvar < *nfluxes; nnvar++) {
	    (*flux_z)[lb][k][j][i][nnvar] = 0.;
	  }
	}
      }
    }
  }


  c_amr_1blk_guardcell_reset();

  *lrefine_max = 5;
  *lrefine_min = 1;

  /* setup initial grid state */

  *lnblocks = 0;

  if (mype == 0) {
    *lnblocks = 1;
    for (i = 0; i < *mdim; i++) {
      (*bsize)[0][i] = 1.;
    }
    for (i = 0; i < *mdim; i++) {
      (*coord)[0][i] = .5;
    }
    for (i = 0; i < *mdim; i++) {
      (*bnd_box)[0][i][0] = 0.;
      (*bnd_box)[0][i][1] = 1.;
    }
    (*nodetype)[0] = 1;
    (*lrefine)[0] = 1;
    
    for (j = 0; j < *mfaces; j++) {
      for (i = 0; i < 2; i++) {
	(*neigh)[0][j][i] = -21;
      }
    }

    *refine[0] = 1;
  }

  for (i = 0; i < *mboundaries; i++) {
    (*boundary_index)[i] = -21;
  }
  for (k = 0; k < 2; k++) {
    for (j = 1; j < 3; j++) {
      (*boundary_box)[k][j][0] = -1.e10;
      (*boundary_box)[k][j][1] =  1.e10;
    }
  }
  (*boundary_box)[0][0][0] = -1.e10;
  (*boundary_box)[0][0][1] = 0.;
  (*boundary_box)[1][0][0] = 1.;
  (*boundary_box)[0][0][1] = 1.e10;
  if (*ndim >= 2) {
    (*boundary_box)[2][0][0] = -1.e10;
    (*boundary_box)[3][0][0] = -1.e10;
    (*boundary_box)[2][0][1] =  1.e10;
    (*boundary_box)[3][0][1] =  1.e10;
    (*boundary_box)[2][2][0] = -1.e10;
    (*boundary_box)[3][2][0] = -1.e10;
    (*boundary_box)[2][2][1] =  1.e10;
    (*boundary_box)[3][2][1] =  1.e10;
    (*boundary_box)[2][1][0] = -1.e10;
    (*boundary_box)[2][1][1] = 0.;
    (*boundary_box)[3][1][0] = 1.;
    (*boundary_box)[3][1][1] = 1.e10;
  }
  if (*ndim == 3) {
    (*boundary_box)[4][0][0] = -1.e10;
    (*boundary_box)[4][1][0] = -1.e10;
    (*boundary_box)[5][0][0] = -1.e10;
    (*boundary_box)[5][1][0] = -1.e10;
    (*boundary_box)[4][0][1] =  1.e10;
    (*boundary_box)[4][1][1] =  1.e10;
    (*boundary_box)[5][0][1] =  1.e10;
    (*boundary_box)[5][1][1] =  1.e10;
    (*boundary_box)[4][2][0] = -1.e10;
    (*boundary_box)[4][2][1] = 0.;
    (*boundary_box)[5][2][0] = 1.;
    (*boundary_box)[5][2][1] = 1.e10;
  }

  /* START THE TEST */

  /* Set the solution array to be the grid points x,y, or z */

  for (lb = 0; lb < *lnblocks; lb++) {
    if ((*nodetype)[lb] == 1 || *advance_all_levels == 1) {

      if (*ndim == 3) {dz = (*bsize)[lb][2]/((float) *nzb);}
      dy = (*bsize)[lb][1]/((float) *nyb);
      dx = (*bsize)[lb][0]/((float) *nxb);

      for (k = nguard0**k3d; k < *nzb + nguard0**k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0; i++) {
	x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-nguard0);
	y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nvar; ivar++) {
	  (*unk)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}
      }
      }
      }
	
      if (*nfacevar > 0) {

      for (k = nguard0**k3d; k < *nzb + nguard0**k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0 + 1; i++) {
	x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-nguard0);
	y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nfacevar; ivar++) {
	  (*facevarx)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }

      if (*ndim >= 2) {
      for (k = nguard0**k3d; k < *nzb + nguard0**k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d + *k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0; i++) {
	x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-nguard0);
	y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nfacevar; ivar++) {
	  (*facevary)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }
      }

      if (*ndim == 3) {
      for (k = nguard0**k3d; k < *nzb + nguard0**k3d + *k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0; i++) {
	x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-nguard0);
	y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nfacevar; ivar++) {
	  (*facevarz)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }
      }

      } /* end nfacevar */

      if (*nvaredge > 0) {

      for (k = nguard0**k3d; k < *nzb + nguard0**k3d + *k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d + *k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0; i++) {
	x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-nguard0);
	y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nvaredge; ivar++) {
	  (*unk_e_x)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }

      if (*ndim >= 2) {
      for (k = nguard0**k3d; k < *nzb + nguard0**k3d + *k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0 + 1; i++) {
	x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-nguard0);
	y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nvaredge; ivar++) {
	  (*unk_e_y)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }
      }

      if (*ndim == 3) {
      for (k = nguard0**k3d; k < *nzb + nguard0**k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d + *k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0 + 1; i++) {
	x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-nguard0);
	y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nvaredge; ivar++) {
	  (*unk_e_z)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }
      }

      } /* end nvaredge */

      if (*nvarcorn > 0) {

      for (k = nguard0**k3d; k < *nzb + nguard0**k3d + *k3d; k++) {
      for (j = nguard0**k2d; j < *nyb + nguard0**k2d + *k2d; j++) {
      for (i = nguard0; i < *nxb + nguard0 + 1; i++) {
	x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-nguard0);
	y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-nguard0);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-nguard0);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nvarcorn; ivar++) {
	  (*unk_n)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }

      } /* end nvarcorn */

    }
  }


  /* REFINE THE MESH SEVERAL TIMES */

  loop_count = 0;
  while(loop_count < 3) {
    for (lb = 0; lb < *lnblocks; lb++) {
      (*refine)[lb] = 0;
    }
    if (loop_count < 2) {
      for (lb = 0; lb < *lnblocks; lb++) {
	(*refine)[lb] = 1;
      }
    } else {
      if (loop_count == 2) {
	if (*ndim == 3) {
	  for (lb = 0; lb < *lnblocks; lb++) {
	    if ((*coord)[lb][0] == 0.125 && (*coord)[lb][1] == 0.125 && (*coord)[lb][2] == 0.125) {(*refine)[lb] = 1;}
	    if ((*coord)[lb][0] == 0.375 && (*coord)[lb][1] == 0.375 && (*coord)[lb][2] == 0.375) {(*refine)[lb] = 1;}
	    if ((*coord)[lb][0] == 0.625 && (*coord)[lb][1] == 0.875 && (*coord)[lb][2] == 0.875) {(*refine)[lb] = 1;}
	  }
	} else {
	  for (lb = 0; lb < *lnblocks; lb++) {
	    if ((*coord)[lb][0] == 0.125 && (*coord)[lb][1] == 0.125) {(*refine)[lb] = 1;}
	    if ((*coord)[lb][0] == 0.375 && (*coord)[lb][1] == 0.375) {(*refine)[lb] = 1;}
	    if ((*coord)[lb][0] == 0.625 && (*coord)[lb][1] == 0.875) {(*refine)[lb] = 1;}
	  }	    
	}
      }
    }

    c_amr_refine_derefine();

    /* prolong the data to the new mesh */

    tag_offset = 100;
    c_mpi_morton_bnd_prolong(&mype,&nprocs,&tag_offset);

    iopt = 1;
    nlayers = *nguard;
    c_amr_prolong(&mype, &iopt, &nlayers);

    /* test fluxes */
    if (*nfluxvar > 0) {

      for (lb = 0; lb < *lnblocks; lb++) {

      if (*ndim == 3) {dz = (*bsize)[lb][2]/((float) *nzb);}
      dy = (*bsize)[lb][1]/((float) *nyb);
      dx = (*bsize)[lb][0]/((float) *nxb);

      for (k = 0; k < *ku_bndi-nguard0**k3d; k++) {
      for (j = 0; j < *ju_bndi-nguard0; j++) {
      for (i = 0; i < 2;        i++) {
	x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + (*bsize)[lb][0]*(float)i;
	y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)(j+1);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)(k+1);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nfluxvar; ivar++) {
	  (*flux_x)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }

      if (*ndim >= 2) {
      for (k = 0; k < *ku_bndi-nguard0**k3d; k++) {
      for (j = 0; j < 2;        j++) {
      for (i = 0; i < *iu_bndi-nguard0; i++) {
	x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)(i+1);
	y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + (*bsize)[lb][1]*(float)j;
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)(k+1);
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nfluxvar; ivar++) {
	  (*flux_y)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }
      }

      if (*ndim == 3) {
      for (k = 0;          k < 2;        k++) {
      for (j = 0; j < *ju_bndi-nguard0; j++) {
      for (i = 0; i < *iu_bndi-nguard0; i++) {
	x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)(i+1);
	y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)(j+1);
	if (*ndim == 3) {
	  z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + (*bsize)[lb][2]*(float)k;
	} else {
	  z0 = 0.;
	}
	value = ax*x0 + ay*y0 + az*z0;
	for (ivar = 0; ivar < *nfluxvar; ivar++) {
	  (*flux_z)[lb][k][j][i][ivar] = value * (float)(ivar+1);
	}

      }
      }
      }
      }

      }

    } /* end nfluxvar */

    nsub = 0;
    flux_dir = 0;
    c_amr_flux_conserve(&mype, &nsub, &flux_dir);

    if (*no_permanent_guardcells == 0) {
      iopt = 1;
      nlayers = *nguard;
      c_amr_guardcell(&mype, &iopt, &nlayers, 
		      &nlayers, &nlayers, &nlayers);
    } else {

      iopt = -1;
      if (*advance_all_levels == 0) {
	i = 1;
	j = 0;
	k = 1;
	c_amr_restrict(&mype,&i,&j,&k);
      } else {
	iopt = -1;
	c_amr_1blk_copy_soln(&iopt);
      }

      iopt = 1;
      lcc = 0;
      lfc = 0;
      lec = 0;
      lnc = 0;
      if (*nvar > 0) {lcc = 1;}
      if (*nfacevar > 0) {lfc = 1;}
      if (*nvaredge > 0) {lec = 1;}
      if (*nvarcorn > 0) {lnc = 1;}
      tag_offset = 100;
      lguard = 1;
      lprolong = 0;
      lflux = 0;
      ledge = 0;
      lrestrict = 0;
      lfulltree = 0;
      nlayersx = *nguard;
      nlayersy = *nguard;
      nlayersz = *nguard;
      flux_dir = 0;

      c_mpi_amr_comm_setup(&mype, &nprocs, &lguard, &lprolong, &lflux, 
			   &ledge, &lrestrict, &lfulltree, &iopt, 
			   &lcc, &lfc, &lec, &lnc, &tag_offset, 
			   &nlayersx, &nlayersy, &nlayersz, &flux_dir);
    }

    loop_count = loop_count + 1;
  }
      

  /* test of unk communications */

  for (ii = 0; ii < nprocs; ii++) {

    if (mype == ii) {
      
      for (lb = 0; lb < *lnblocks; lb++) {

	if (*no_permanent_guardcells == 1) {

	  iopt = 1;
	  nlayers = *nguard;
	  nlayersx = *nguard;
	  nlayersy = *nguard;
	  nlayersz = *nguard;
	  lcc = 0;
	  lfc = 0;
	  lec = 0;
	  lnc = 0;
	  if (*nvar > 0) {lcc = 1;}
	  if (*nfacevar > 0) {lfc = 1;}
	  if (*nvaredge > 0) {lec = 1;}
	  if (*nvarcorn > 0) {lnc = 1;}
	  l_srl_only = 0;
	  icoord = 0;
	  ldiag = 0;
	  if (*diagonals == 1) {ldiag = 1;}
	  if ((*nodetype)[lb] == 1 || *advance_all_levels == 1) {

	    c_amr_1blk_guardcell(&mype, &iopt, &nlayers, &lb, &mype, 
				 &lcc, &lfc, &lec, &lnc, &l_srl_only, 
				 &icoord, &ldiag, 
				 &nlayersx, &nlayersy, &nlayersz); 

	    idest = 0; 
	    c_amr_1blk_to_perm(&lcc, &lfc, &lec, &lnc, &lb, &iopt, &idest);

	  } 
	  
	} else {

	  for(k = 0; k < *ku_bnd1; k++) {
	    for (j = 0; j < *ju_bnd1; j++) {
	      for (i = 0; i < *iu_bnd1; i++) {
		for (ivar = 0; ivar < *nvar; ivar++) {
		  (*unk1)[0][k][j][i][ivar] = (*unk)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }

	  if (*nfacevar > 0) {

	  for(k = 0; k < *ku_bnd1; k++) {
	    for (j = 0; j < *ju_bnd1; j++) {
	      for (i = 0; i < *iu_bnd1+1; i++) {
		for (ivar = 0; ivar < *nfacevar; ivar++) {
		  (*facevarx1)[0][k][j][i][ivar] = (*facevarx)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }

	  if (*ndim >= 2) {
	  for(k = 0; k < *ku_bnd1; k++) {
	    for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	      for (i = 0; i < *iu_bnd1; i++) {
		for (ivar = 0; ivar < *nfacevar; ivar++) {
		  (*facevary1)[0][k][j][i][ivar] = (*facevary)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }
	  }

	  if (*ndim == 3) {
	  for(k = 0; k < *ku_bnd1 + *k3d; k++) {
	    for (j = 0; j < *ju_bnd1; j++) {
	      for (i = 0; i < *iu_bnd1; i++) {
		for (ivar = 0; ivar < *nfacevar; ivar++) {
		  (*facevarz1)[0][k][j][i][ivar] = (*facevarz)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }
	  }

	  }  /* end if nfacevar */
	    
	  if (*nvaredge > 0) {

	  for(k = 0; k < *ku_bnd1 + *k3d; k++) {
	    for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	      for (i = 0; i < *iu_bnd1; i++) {
		for (ivar = 0; ivar < *nedgevar; ivar++) {
		  (*unk_e_x1)[0][k][j][i][ivar] = (*unk_e_x)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }

	  if (*ndim >= 2) {
	  for(k = 0; k < *ku_bnd1 + *k3d; k++) {
	    for (j = 0; j < *ju_bnd1; j++) {
	      for (i = 0; i < *iu_bnd1 + 1; i++) {
		for (ivar = 0; ivar < *nvaredge; ivar++) {
		  (*unk_e_y1)[0][k][j][i][ivar] = (*unk_e_y)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }
	  }

	  if (*ndim == 3) {
	  for(k = 0; k < *ku_bnd1; k++) {
	    for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	      for (i = 0; i < *iu_bnd1 + 1; i++) {
		for (ivar = 0; ivar < *nvaredge; ivar++) {
		  (*unk_e_z1)[0][k][j][i][ivar] = (*unk_e_z)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }
	  }

	  }  /* end if nvaredge */

	  if (*nvarcorn > 0) {

	  for(k = 0; k < *ku_bnd1 + *k3d; k++) {
	    for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	      for (i = 0; i < *iu_bnd1 + 1; i++) {
		for (ivar = 0; ivar < *nvarcorn; ivar++) {
		  (*unk_n1)[0][k][j][i][ivar] = (*unk_n)[lb][k][j][i][ivar];
		}
	      }
	    }
	  }

	  } /* end if nvarcorn */

 	}

	/* Perform check of solution everywhere, 
	   including in the guardcells */

	if ((*nodetype)[lb] == 1 || *advance_all_levels == 1) {

	  if (*ndim == 3) dz = (*bsize)[lb][2]/(float)*nzb;
	  dy = (*bsize)[lb][1]/(float)*nyb;
	  dx = (*bsize)[lb][0]/(float)*nxb;

	  for (k = 0; k < *ku_bnd1; k++) {	
	  for (j = 0; j < *ju_bnd1; j++) {
	  for (i = 0; i < *iu_bnd1; i++) { 

	    x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-*nguard);
	    y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nvar; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*unk1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in unk at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("unk1, value = %f %f\n",
			(*unk1)[0][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }

	  if (*nfacevar > 0) {

	    /* facevarx */

	  for (k = 0; k < *ku_bnd1; k++) {	
	  for (j = 0; j < *ju_bnd1; j++) {
	  for (i = 0; i < *iu_bnd1+1; i++) { 

	    x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-*nguard);
	    y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nfacevar; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*facevarx1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in facevarx at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("facevarx1, value = %f %f\n",
			(*facevarx1)[0][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }

	  /* facevary */

	  if (*ndim >= 2) {
	  for (k = 0; k < *ku_bnd1; k++) {	
	  for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	  for (i = 0; i < *iu_bnd1; i++) { 

	    x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-*nguard);
	    y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nfacevar; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*facevary1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in facevary at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("facevary1, value = %f %f\n",
			(*facevary1)[0][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }
	  }

	  /* facevarz */

	  if (*ndim == 3) {
	  for (k = 0; k < *ku_bnd1 + *k3d; k++) {	
	  for (j = 0; j < *ju_bnd1; j++) {
	  for (i = 0; i < *iu_bnd1; i++) { 

	    x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-*nguard);
	    y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nfacevar; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*facevarz1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in facevarz at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("facevarz1, value = %f %f\n",
			(*facevarz1)[0][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }
	  }
	
	  } /* end nfacevar */

	  if (*nfluxvar > 0) {

	  /* flux_x */

	  for (k = 0; k < *ku_bndi-nguard0**k3d; k++) {	
	  for (j = 0; j < *ju_bndi-nguard0; j++) {
	  for (i = 0; i < 2;                i++) { 

	    x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + (*bsize)[lb][0]*(float)i;
	    y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)(j+1);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)(k+1); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nfluxvar; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*flux_x)[lb][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in flux_x at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("flux_x, value = %f %f\n",
			(*flux_x)[lb][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }

	  /* flux_y */

	  if (*ndim >= 2) {
	  for (k = 0; k < *ku_bndi-nguard0**k3d; k++) {	
	  for (j = 0; j < 2;                j++) {
	  for (i = 0; i < *iu_bndi-nguard0; i++) { 

	    x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)(i+1);
	    y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + (*bsize)[lb][1]*(float)j;
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)(k+1); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nfluxvar; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*flux_y)[lb][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in flux_y at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("flux_y, value = %f %f\n",
			(*flux_y)[lb][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }
	  }

	  /* flux_z */

	  if (*ndim == 3) {
	  for (k = 0; k < 2;                k++) {	
	  for (j = 0; j < *ju_bndi-nguard0; j++) {
	  for (i = 0; i < *iu_bndi-nguard0; i++) { 

	    x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)(i+1);
	    y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)(j+1);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + (*bsize)[lb][2]*(float)k; 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nfluxvar; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*flux_z)[lb][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in flux_z at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("flux_z, value = %f %f\n",
			(*flux_z)[lb][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }
	  }

	  } /* end nfluxvar */

	  if (*nvaredge > 0) {

	    /* unk_e_x */

	  for (k = 0; k < *ku_bnd1 + *k3d; k++) {	
	  for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	  for (i = 0; i < *iu_bnd1; i++) { 

	    x0 = (*coord)[lb][0] - .5*((*bsize)[lb][0]+dx) + dx*(float)((i+1)-*nguard);
	    y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nvaredge; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*unk_e_x1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in unk_e_x at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("unk_e_x1, value = %f %f\n",
			(*unk_e_x1)[0][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }

	  /* unk_e_y */

	  if (*ndim >= 2) {
	  for (k = 0; k < *ku_bnd1 + *k3d; k++) {	
	  for (j = 0; j < *ju_bnd1; j++) {
	  for (i = 0; i < *iu_bnd1 + 1; i++) { 

	    x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-*nguard);
	    y0 = (*coord)[lb][1] - .5*((*bsize)[lb][1]+dy) + dy*(float)((j+1)-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nvaredge; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*unk_e_y1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in unk_e_y at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("unk_e_y1, value = %f %f\n",
			(*unk_e_y1)[0][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }
	  }

	  /* unk_e_z */

	  if (*ndim == 3) {
	  for (k = 0; k < *ku_bnd1; k++) {	
	  for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	  for (i = 0; i < *iu_bnd1 + 1; i++) { 

	    x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-*nguard);
	    y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*((*bsize)[lb][2]+dz) + dz*(float)((k+1)-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nvaredge; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*unk_e_z1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in unk_e_z at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("unk_e_z1, value = %f %f\n",
			(*unk_e_z1)[0][k][j][i][ivar], 
			value);
	      }

	    }
	  }
	  }
	  }
	  }

	  } /* end nvaredge */

	  if (nvarcorn > 0) {

	  for (k = 0; k < *ku_bnd1 + *k3d; k++) {	
	  for (j = 0; j < *ju_bnd1 + *k2d; j++) {
	  for (i = 0; i < *iu_bnd1 + 1; i++) { 

	    x0 = (*coord)[lb][0] - .5*(*bsize)[lb][0] + dx*(float)(i-*nguard);
	    y0 = (*coord)[lb][1] - .5*(*bsize)[lb][1] + dy*(float)(j-*nguard);
	    if (*ndim == 3) {
	      z0 = (*coord)[lb][2] - .5*(*bsize)[lb][2] + dz*(float)(k-*nguard); 
	    } else {
	      z0 = 0.;
	    }

	    for (ivar = 0; ivar < *nvarcorn; ivar++) {

	      value = ax*x0 + ay*y0 + az*z0;
	      value = value*(float)(ivar+1);

	      if (fabs(value - (*unk_n1)[0][k][j][i][ivar]) > eps) {
		total_errors = total_errors + 1;
		printf ("ERROR in unk_n at ivar,i,j,k,lb %d %d %d %d %d\n",
			ivar,i,j,k,lb);
		printf ("unk_n1, value = %f %f\n",
			(*unk_n1)[0][k][j][i][ivar], 
			value);
	      }
	    }
	  }
	  }
	  }

	  } /* end if nvarcorn */

	} /* end if nodetype */ 
	} /* end loop over blocks */
      } /* end if (mype */
      MPI_Barrier(MPI_COMM_WORLD);
    } /* end for (ii */


  iunit1 = 96;
  // Causes problems on thunderhead because it is parallel io
  //c_amr_plotfile_chombo(&iunit1);

  printf("Wrote chombovis file\n");

  char check_format[80];
  for (i = 0; i < 80; i++) {
    check_format[i] = 0;
  }
  strncpy(check_format,"default",strlen("default"));

  int l_with_guardcells;
  float user1;
  float user2;
  float user3;
  float user4;
  float user5;
  l_with_guardcells = 0;
  user1 = 0.;
  user2 = 0.;
  user3 = 0.;
  user4 = 0.;
  user5 = 0.;

  // Had to add this print to make things work on thunderhead with the
  // NAG compiler.  Don't know why !
  printf("check_format = %s\n",check_format);

  c_amr_checkpoint_wr(&iunit1,&l_with_guardcells,\
		      check_format,&user1,&user2,&user3,&user4,&user5);

  if (mype == 0) {
    printf("checkpoint file written\n");
  }

  // Had to add this print to make things work on thunderhead with the
  // NAG compiler.  Don't know why !
  printf("check_format = %s\n",check_format);

  c_amr_checkpoint_re(&iunit1,&l_with_guardcells,\
		      check_format,&user1,&user2,&user3,&user4,&user5);

  if (mype == 0) {
    printf("checkpoint file read back in\n");
  }

  if (mype == 0) {
    for (i = 0; i < 80; i++) {
      file_name[i] = 0;
    }
    int iloop;
    // strip off fortran junk from output_dir
    for(iloop=0; iloop<=strlen(output_dir)-1; iloop++) {
      if (output_dir[iloop] == ' ') { output_dir[iloop]=0; break; }
    }
    strncpy(file_name,output_dir,strlen(output_dir));
    strncat(file_name,"test.log",strlen("test.log"));
    test_log = fopen(file_name,"w+");
    fprintf(test_log,"test_c_interface completed succesfully\n");
    fprintf(test_log,"total errors = %d\n",total_errors);
    printf("test_c_interface completed succesfully\n");
    printf("total errors = %d\n",total_errors);
    if (total_errors > 0) {
      fprintf(test_log,"test_c_interface failed\n");
    }
    fclose(test_log);
  }
  
  c_amr_close();
  
}
