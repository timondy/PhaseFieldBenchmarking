#include "mpi.h"
#include <ctype.h>

#include "underscore.h"
#define float double
#define MPI_FLOAT MPI_DOUBLE


void get_new_type(int nxb, int nyb, int nzb, int nvar, int iorder[],
		  int nx, int ny, int nz,
                  MPI_Datatype* new_type) {

  MPI_Datatype type1;
  MPI_Datatype type2;

  int iii; int udim[4]; int udim_tot[4];
  for (iii = 0; iii < 4; iii++) {
    if (iorder[iii] == 1) {
      udim[iii] = 1;
      udim_tot[iii] = nvar;
    }
    if (iorder[iii] == 2) {
      udim[iii] = nxb;
      udim_tot[iii] = nx;
    }
    if (iorder[iii] == 3) {
      udim[iii] = nyb;
      udim_tot[iii] = ny;
    }
    if (iorder[iii] == 4) {
      udim[iii] = nzb;
      udim_tot[iii] = nz;
    }
  }

  MPI_Type_hvector(udim[1],
		   udim[0],
		   udim_tot[0] * sizeof(float),
		   MPI_FLOAT,
		   &type1);
  MPI_Type_commit(&type1);
  MPI_Type_hvector(udim[2],
		   1,
		   udim_tot[0] * udim_tot[1] * sizeof(float),
		   type1,
		   &type2);
  MPI_Type_commit(&type2);
  MPI_Type_hvector(udim[3],
		   1,
		   udim_tot[0] * udim_tot[1] * udim_tot[2] * sizeof(float),
		   type2,
		   new_type);
  MPI_Type_commit(new_type);

  MPI_Type_free(&type1);
  MPI_Type_free(&type2);

}

void get_udim_mpiio(int nvar, int nx, int ny, int nz, 
	      int ivar, int istart, int jstart, int kstart,
	      int iorder[4], int udim[4], int it[4]) {

  int iii;
  for (iii = 0; iii < 4; iii++) {
    if (iorder[iii] == 1) {
      udim[iii] = nvar;
      it[iii] = ivar;
    }
    if (iorder[iii] == 2) {
      udim[iii] = nx;
      it[iii] = istart;
    }
    if (iorder[iii] == 3) {
      udim[iii] = ny;
      it[iii] = jstart;
    }
    if (iorder[iii] == 4) {
      udim[iii] = nz;
      it[iii] = kstart;
    }
  }
}  

#ifdef UNDERSCORE
void write_blocks_mpiio_r8_ (char* file_name_in, 
#endif
#ifdef DOUBLE_UNDERSCORE
void write_blocks_mpiio_r8__ (char* file_name_in, 
#endif
#ifndef UNDERSCORE 
#ifndef DOUBLE_UNDERSCORE
void write_blocks_mpiio_r8 (char* file_name_in, 
#endif
#endif
			   int*  tot_blocks,
			   int*  tot_blocks_wr,
			   int*  max_lnblocks,
			   int*  lnblocks,
			   int*  n_to_left,
			   int*  mdim,
			   int*  ndim,
			   int*  ngid,
			   int*  mflags,
			   int   lrefine[], 
			   int   nodetype[], 
			   int   which_child[], 
			   int   gid[], 
			   int   bflags[],
			   float coord[],
			   float bnd_box[],
			   float work_block[],
			   float unk[],
			   int*  nvar,
			   int*  nvar_chk_cc,
			   int   checkp_on_cc[],
			   float facevarx[], float facevary[], float facevarz[],
			   int*  nbndvar,
			   int*  nvar_chk_fc,
			   int   checkp_on_fc[],
			   float unk_e_x[], float unk_e_y[], float unk_e_z[],
			   int*  nbndvare,
			   int*  nvar_chk_ec,
			   int   checkp_on_ec[],
			   float unk_n[],
			   int*  nbndvarc,
			   int*  nvar_chk_nc,
			   int   checkp_on_nc[],
			   int*  nx,  int* ny, int* nz,
			   int*  il0, int* iu0,
			   int*  jl0, int* ju0, 
			   int*  kl0, int* ku0,
			   float* user_attr_1,
			   float* user_attr_2,
			   float* user_attr_3,
			   float* user_attr_4,
			   float* user_attr_5,
			   int iorder[])
{

  int i;
  int ii;
  char file_name[80];
  MPI_File file_handle;
  MPI_Status status;
  MPI_Offset disp;
  MPI_Offset disp_old;
  MPI_Datatype etype;
  MPI_Datatype filetype;
  int mype;
  int nprocs;
  int istart;
  int jstart;
  int kstart;
  int iend;
  int jend;
  int kend;
  int lb;
  int k;
  int j;
  int ivar;
  int k2d, k3d;

  MPI_Datatype unk_type;
  MPI_Datatype unk_type1;
  MPI_Datatype unk_type2;
  int unk_size;

  MPI_Datatype facevarx_type;
  int facevarx_size;
  MPI_Datatype facevary_type;
  int facevary_size;
  MPI_Datatype facevarz_type;
  int facevarz_size;

  MPI_Datatype unk_e_x_type;
  int unk_e_x_size;
  MPI_Datatype unk_e_y_type;
  int unk_e_y_size;
  MPI_Datatype unk_e_z_type;
  int unk_e_z_size;

  MPI_Datatype unk_n_type;
  int unk_n_size;


  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /*
     Don't know why I had to add this code, but it appears that some
     non-numeric, non-chararcter data gets inserted into this file name
     when called from the c_interface.  This code will remove any weird
     characters
  */
  for(i=0; i<80; i++) {
    if (!isalnum(file_name_in[i]) && !ispunct(file_name_in[i])) {
      for (ii=i;ii<79;ii++) {
        file_name_in[ii]=file_name_in[ii+1];
      }
    }
  }

  /* strip off fortran junk from the file_name */
  strncpy(file_name, file_name_in, 80);
  for(i=79; i>=0; i--) {
    if (file_name[i] != ' ') { file_name[i+1]=0; break; }
  }


  /* OPEN THE FILE */

  MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handle);



  /* WRITE THE USER ATTRIBUTES */
  disp = (MPI_Offset)0;

  if (mype == 0) {
    MPI_File_write_at(file_handle, disp, tot_blocks, 1, MPI_INT, &status);
  }
  disp = disp + (MPI_Offset)sizeof(int);

  if (mype == 0) {
    MPI_File_write_at(file_handle, disp, user_attr_1, 1, MPI_FLOAT, &status);
    disp = disp + (MPI_Offset)sizeof(float);
    MPI_File_write_at(file_handle, disp, user_attr_2, 1, MPI_FLOAT, &status);
    disp = disp + (MPI_Offset)sizeof(float);
    MPI_File_write_at(file_handle, disp, user_attr_3, 1, MPI_FLOAT, &status);
    disp = disp + (MPI_Offset)sizeof(float);
    MPI_File_write_at(file_handle, disp, user_attr_4, 1, MPI_FLOAT, &status);
    disp = disp + (MPI_Offset)sizeof(float);
    MPI_File_write_at(file_handle, disp, user_attr_5, 1, MPI_FLOAT, &status);
  }

  disp_old = (MPI_Offset)(sizeof(int) + 5*sizeof(float));

  if (*lnblocks > 0) {

  /* lrefine */
  disp = disp_old + (MPI_Offset)(*n_to_left*sizeof(int));
  MPI_File_write_at(file_handle, disp, lrefine, *lnblocks, MPI_INT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*sizeof(int));

  /* nodetype */
  disp = disp_old + (MPI_Offset)(*n_to_left*sizeof(int));
  MPI_File_write_at(file_handle, disp, nodetype, *lnblocks, MPI_INT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*sizeof(int));

  /* which_child */
  disp = disp_old + (MPI_Offset)(*n_to_left*sizeof(int));
  MPI_File_write_at(file_handle, disp, which_child, *lnblocks, MPI_INT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*sizeof(int));

  /* gid */
  disp = disp_old + (MPI_Offset)(*n_to_left*(*ngid)*sizeof(int));
  MPI_File_write_at(file_handle, disp, gid, *lnblocks*(*ngid), MPI_INT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*(*ngid)*sizeof(int));

  /* bflags */
  disp = disp_old + (MPI_Offset)(*n_to_left*(*mflags)*sizeof(int));
  MPI_File_write_at(file_handle, disp, bflags, *lnblocks*(*mflags), MPI_INT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*(*mflags)*sizeof(int));

  /* coord */
  disp = disp_old + (MPI_Offset)(*n_to_left*(*mdim)*sizeof(float));
  MPI_File_write_at(file_handle, disp, coord, *lnblocks*(*mdim), MPI_FLOAT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*(*mdim)*sizeof(float));

  /* bnd_box */
  disp = disp_old + (MPI_Offset)(*n_to_left*(2 * *mdim)*sizeof(float));
  MPI_File_write_at(file_handle, disp, bnd_box, *lnblocks*(2 * *mdim), MPI_FLOAT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*(2 * *mdim)*sizeof(float));

  /* work_block */
  disp = disp_old + (MPI_Offset)(*n_to_left*sizeof(float));
  MPI_File_write_at(file_handle, disp, work_block, *lnblocks, MPI_FLOAT, &status);
  disp_old = disp_old + (MPI_Offset)(*tot_blocks*sizeof(float));

  }

  /* BLOCK DATA (unk, facevar's, unk_e's and unk_n */

  kstart = *kl0;
  kend   = *ku0;
  jstart = *jl0;
  jend   = *ju0;
  istart = *il0;
  iend   = *iu0;
  if (*ndim < 3) {
    kstart = 0;
    kend   = 1;
  }
  if (*ndim < 2) {
    kstart = 0;
    kend   = 1;
    jstart = 0;
    jend   = 1;
  }

  k2d = 0;
  if (*ndim >= 2) k2d = 1;
  k3d = 0;
  if (*ndim == 3) k3d = 1;

  disp = disp_old + 
    (MPI_Offset)((*n_to_left)*
		 (*nvar_chk_cc)*
		 (kend-kstart)*
		 (jend-jstart)*
		 (iend-istart)*
		 sizeof(float))
    + (MPI_Offset)((*n_to_left)*
                   (*nvar_chk_fc)*
                   (kend-kstart)*
                   (jend-jstart)*
                   (iend-istart+1)*
                   sizeof(float))
    + (MPI_Offset)((*n_to_left)*
                   (*nvar_chk_fc)*
                   (kend-kstart)*
                   (jend-jstart+k2d)*
                   (iend-istart)*
                   sizeof(float))
    + (MPI_Offset)((*n_to_left)*
                   (*nvar_chk_fc)*
                   (kend-kstart+k3d)*
                   (jend-jstart)*
                   (iend-istart)*
                   sizeof(float))
    + (MPI_Offset)((*n_to_left)*
                   (*nvar_chk_ec)*
                   (kend-kstart+k3d)*
                   (jend-jstart+k2d)*
                   (iend-istart)*
                   sizeof(float))
    + (MPI_Offset)((*n_to_left)*
                   (*nvar_chk_ec)*
                   (kend-kstart+k3d)*
                   (jend-jstart)*
                   (iend-istart+1)*
                   sizeof(float))
    + (MPI_Offset)((*n_to_left)*
                   (*nvar_chk_ec)*
                   (kend-kstart)*
                   (jend-jstart+k2d)*
                   (iend-istart+1)*
                   sizeof(float))
    + (MPI_Offset)((*n_to_left)*
                   (*nvar_chk_nc)*
                   (kend-kstart+k3d)*
                   (jend-jstart+k2d)*
                   (iend-istart+1)*
                   sizeof(float))
    ;


  get_new_type(iend-istart, jend-jstart, kend-kstart, *nvar,
	       iorder,
	       *nx, *ny, *nz, &unk_type);
  MPI_Type_size(unk_type, &unk_size);

  get_new_type(iend-istart+1, jend-jstart, kend-kstart, *nbndvar,
	       iorder,
	       *nx+1, *ny, *nz, &facevarx_type);
  MPI_Type_size(facevarx_type, &facevarx_size);

  get_new_type(iend-istart, jend-jstart+k2d, kend-kstart, *nbndvar,
	       iorder,
	       *nx, *ny+k2d, *nz, &facevary_type);
  MPI_Type_size(facevary_type, &facevary_size);

  get_new_type(iend-istart, jend-jstart, kend-kstart+k3d, *nbndvar,
	       iorder,
	       *nx, *ny, *nz+k3d, &facevarz_type);
  MPI_Type_size(facevarz_type, &facevarz_size);

  get_new_type(iend-istart, jend-jstart+k2d, kend-kstart+k3d, *nbndvare,
	       iorder,
	       *nx, *ny+k2d, *nz+k3d, &unk_e_x_type);
  MPI_Type_size(unk_e_x_type, &unk_e_x_size);

  get_new_type(iend-istart+1, jend-jstart, kend-kstart+k3d, *nbndvare,
	       iorder,
	       *nx+1, *ny, *nz+k3d, &unk_e_y_type);
  MPI_Type_size(unk_e_y_type, &unk_e_y_size);

  get_new_type(iend-istart+1, jend-jstart+k2d, kend-kstart, *nbndvare,
	       iorder,
	       *nx+1, *ny+k2d, *nz, &unk_e_z_type);
  MPI_Type_size(unk_e_z_type, &unk_e_z_size);

  get_new_type(iend-istart+1, jend-jstart+k2d, kend-kstart+k3d, *nbndvarc,
	       iorder,
	       *nx+1, *ny+k2d, *nz+k3d, &unk_n_type);
  MPI_Type_size(unk_n_type, &unk_n_size);

  int udim[4]; int it[4];
  for (lb = 0; lb < *lnblocks; lb++) {

    /* unk data */

    if (*nvar_chk_cc > 0) {

      for (ivar=0; ivar < *nvar; ivar++) {

	get_udim_mpiio(*nvar, *nx, *ny, *nz, ivar, istart, jstart, kstart, iorder, udim, it);

        if (checkp_on_cc[ivar] == 1) {

          MPI_File_write_at(file_handle, disp,
			   &unk[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				(udim[2] * udim[1] * udim[0] * it[3]) +
				(udim[1] * udim[0] * it[2]) +
				(udim[0] * it[1]) +
				it[0]] ,
			   1, unk_type, &status);

          disp = disp + (MPI_Offset)unk_size;
        }

      } /* end ivar loop */

    } /* end if (nvar_chk_cc */

    if (*nvar_chk_fc > 0) {

      for (ivar = 0; ivar < *nbndvar; ivar++) {
	
	get_udim_mpiio(*nbndvar, *nx+1, *ny, *nz, ivar, istart, jstart, kstart, iorder, udim, it);

	if (checkp_on_fc[ivar + 2 * ivar] == 1) {
	  MPI_File_write_at(file_handle, disp,
			   &facevarx[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				     (udim[2] * udim[1] * udim[0] * it[3]) +
				     (udim[1] * udim[0] * it[2]) +
				     (udim[0] * it[1]) +
				     it[0]] ,
			   1, facevarx_type, &status);
	  
	  disp = disp + (MPI_Offset)facevarx_size;
	}

      } /* end ivar loop */

      /* facevary data */

      for (ivar = 0; ivar < *nbndvar; ivar++) {
	
	get_udim_mpiio(*nbndvar, *nx, *ny+k2d, *nz, ivar, istart, jstart, kstart, iorder, udim, it);

	if (checkp_on_fc[ivar + 1 + 2 * ivar] == 1) {
	  MPI_File_write_at(file_handle, disp,
			   &facevary[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				     (udim[2] * udim[1] * udim[0] * it[3]) +
				     (udim[1] * udim[0] * it[2]) +
				     (udim[0] * it[1]) +
				     it[0]] ,
			   1, facevary_type, &status);
	  
	  disp = disp + (MPI_Offset)facevary_size;
	}

      } /* end ivar loop */

      /* facevarz data */

      for (ivar = 0; ivar < *nbndvar; ivar++) {
	
	get_udim_mpiio(*nbndvar, *nx, *ny, *nz+k3d, ivar, istart, jstart, kstart, iorder, udim, it);

	if (checkp_on_fc[ivar + 2 + 2 * ivar] == 1) {
	  MPI_File_write_at(file_handle, disp,
			   &facevarz[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				     (udim[2] * udim[1] * udim[0] * it[3]) +
				     (udim[1] * udim[0] * it[2]) +
				     (udim[0] * it[1]) +
				     it[0]] ,
			   1, facevarz_type, &status);
	  
	  disp = disp + (MPI_Offset)facevarz_size;
	}
	
      } /* end ivar loop */

    } /* end if (nvar_chk_fc */

    /* unk_e_x data */

    if (*nvar_chk_ec > 0) {

      for (ivar = 0; ivar < *nbndvare; ivar++) {

	get_udim_mpiio(*nbndvare, *nx, *ny+k2d, *nz+k3d, ivar, istart, jstart, kstart, iorder, udim, it);

	if (checkp_on_ec[ivar + 2 * ivar] == 1) {
	  MPI_File_write_at(file_handle, disp,
			   &unk_e_x[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				    (udim[2] * udim[1] * udim[0] * it[3]) +
				    (udim[1] * udim[0] * it[2]) +
				    (udim[0] * it[1]) +
				    it[0]] ,
			   1, unk_e_x_type, &status);
		
	  disp = disp + (MPI_Offset)unk_e_x_size;
	}

      } /* end ivar loop */

      /* unk_e_y data */

      for (ivar = 0; ivar < *nbndvare; ivar++) {

	get_udim_mpiio(*nbndvare, *nx+1, *ny, *nz+k3d, ivar, istart, jstart, kstart, iorder, udim, it);

	if (checkp_on_ec[ivar + 1 + 2 * ivar] == 1) {
	  MPI_File_write_at(file_handle, disp,
			   &unk_e_y[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				    (udim[2] * udim[1] * udim[0] * it[3]) +
				    (udim[1] * udim[0] * it[2]) +
				    (udim[0] * it[1]) +
				    it[0]] ,
			   1, unk_e_y_type, &status);
	  
	  disp = disp + (MPI_Offset)unk_e_y_size;
	}

      } /* end ivar loop */

      /* unk_e_z data */

      for (ivar = 0; ivar < *nbndvare; ivar++) {

	get_udim_mpiio(*nbndvare, *nx+1, *ny+k2d, *nz, ivar, istart, jstart, kstart, iorder, udim, it);

	if (checkp_on_ec[ivar + 2 + 2 * ivar] == 1) {
	  MPI_File_write_at(file_handle, disp,
			   &unk_e_z[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				    (udim[2] * udim[1] * udim[0] * it[3]) +
				    (udim[1] * udim[0] * it[2]) +
				    (udim[0] * it[1]) +
				    it[0]] ,
			   1, unk_e_z_type, &status);
	  
	  disp = disp + (MPI_Offset)unk_e_z_size;
	}

      } /* end ivar loop */

    } /* end if (nvar_chk_ec */

    if (*nvar_chk_nc > 0) {

      for (ivar=0; ivar < *nbndvarc; ivar++) {
	
      get_udim_mpiio(*nbndvarc, *nx+1, *ny+k2d, *nz+k3d, ivar, istart, jstart, kstart, iorder, udim, it);
	if (checkp_on_nc[ivar] == 1) {
	  MPI_File_write_at(file_handle, disp,
			   &unk_n[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
				  (udim[2] * udim[1] * udim[0] * it[3]) +
				  (udim[1] * udim[0] * it[2]) +
				  (udim[0] * it[1]) +
				  it[0]] ,
			   1, unk_n_type, &status);
	  
	  disp = disp + (MPI_Offset)unk_n_size;
	}
	      
      } /* end ivar loop */
      
    } /* end if (nvar_chk_nc */

  } /* end loop over blocks */

  MPI_Type_free(&unk_type);
  MPI_Type_free(&facevarx_type);
  MPI_Type_free(&facevary_type);
  MPI_Type_free(&facevarz_type);
  MPI_Type_free(&unk_e_x_type);
  MPI_Type_free(&unk_e_y_type);
  MPI_Type_free(&unk_e_z_type);
  MPI_Type_free(&unk_n_type);

  /* close the file */

  MPI_File_close(&file_handle);
  
}

     
