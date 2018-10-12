#include "hdf5.h"
#include "mpi.h"
#include <stdio.h>
#include <ctype.h>

#include "underscore.h"
#define float double
#define H5T_NATIVE_FLOAT H5T_NATIVE_DOUBLE

void CREATE_Block_data_type_hdf5_r8(hid_t*,int*,int*,int*,int*,int*,int*,
				    int*,int*,int*,int*,int*,int*,int*,int*);

void OPEN_File_hdf5_r8 (char* file_name_in, hid_t* file_id) {

  char file_name[80];
  int i; int ii;

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

  /* strip off fortran junk from the filename */
  strncpy(file_name, file_name_in, 80);
  for(i=79; i>=0; i--) {
    if (file_name[i] != ' ') { file_name[i+1]=0; break; }
  }

  *file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
}


void READ_Header_hdf5_r8(hid_t* file_id, int* tot_blocks,
			 float* user_attr_1,
			 float* user_attr_2,
			 float* user_attr_3,
			 float* user_attr_4,
			 float* user_attr_5) {

  hid_t    dset_header, group_id, dataspace_header, attribute_id;
  herr_t   status;
  int      temp[1], mpi_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  dset_header = H5Dopen(*file_id, "Global Data");
  /* all procs read from the same location in the same file, in parallel */
  /* broadcast not needed */
  status = H5Dread(dset_header, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
		   H5P_DEFAULT, temp);
  H5Dclose(dset_header);

  *tot_blocks = temp[0];

  group_id = H5Gopen(*file_id,"/");

  /* attribute 1 */
  attribute_id = H5Aopen_name(group_id, "user attribute 1");
  if (mpi_rank == 0) {
    status = H5Aread(attribute_id, H5T_NATIVE_FLOAT, user_attr_1);
  }
  H5Aclose(attribute_id);

  /* attribute 2 */
  attribute_id = H5Aopen_name(group_id, "user attribute 2");
  if (mpi_rank == 0) {
    status = H5Aread(attribute_id, H5T_NATIVE_FLOAT, user_attr_2);
  }
  H5Aclose(attribute_id);

  /* attribute 3 */
  attribute_id = H5Aopen_name(group_id, "user attribute 3");
  if (mpi_rank == 0) {
    status = H5Aread(attribute_id, H5T_NATIVE_FLOAT, user_attr_3);
  }
  H5Aclose(attribute_id);

  /* attribute 4 */
  attribute_id = H5Aopen_name(group_id, "user attribute 4");
  if (mpi_rank == 0) {
    status = H5Aread(attribute_id, H5T_NATIVE_FLOAT, user_attr_4);
  }
  H5Aclose(attribute_id);

  /* attribute 5 */
  attribute_id = H5Aopen_name(group_id, "user attribute 5");
  if (mpi_rank == 0) {
    status = H5Aread(attribute_id, H5T_NATIVE_FLOAT, user_attr_5);
  }
  H5Aclose(attribute_id);

  H5Gclose(group_id);
}


void OPEN_Block_data_set_hdf5_r8(hid_t* dset_block, 
				 hid_t* file_id)
{

  *dset_block = H5Dopen(*file_id, "Block Data");

}


void var_to_block_pointer(float unk[], int iorder[], int nvar, int nx, int ny, int nz,
			  int kstart, int kend,
			  int jstart, int jend,
			  int istart, int iend,
			  int checkp_on[], int off1, int off2, int lb,
			  int *off, float *block_pointer)
{

  int udim[4]; int iii;
  for (iii = 0; iii < 4; iii++) {
    if (iorder[iii] == 1) {
      udim[iii] = nvar;
	  }
    if (iorder[iii] == 2) {
      udim[iii] = nx;
    }
    if (iorder[iii] == 3) {
      udim[iii] = ny;
    }
    if (iorder[iii] == 4) {
      udim[iii] = nz;
    }
  }
  
  int i; int j; int k; int ivar;
  for (k = kstart; k < kend; k++) {
    for (j = jstart; j < jend; j++) {
      for (i = istart; i < iend; i++) {
	for (ivar=0; ivar < nvar; ivar++) {
	  
	  int it[4];
	  for (iii = 0; iii < 4; iii++) {
	    if (iorder[iii] == 1) {
	      it[iii] = ivar;
	    }
	    if (iorder[iii] == 2) {
	      it[iii] = i;
	    }
	    if (iorder[iii] == 3) {
	      it[iii] = j;
	    }
	    if (iorder[iii] == 4) {
	      it[iii] = k;
	    }
	  }
	  
	  if (checkp_on[ivar + off1*ivar + off2] == 1) {
	    unk[(udim[3] * udim[2] * udim[1] * udim[0] * lb) +
		(udim[2] * udim[1] * udim[0] * it[3]) +
		(udim[1] * udim[0] * it[2]) +
		(udim[0] * it[1]) +
		it[0]] =
	      *(block_pointer+*off);
	    *off = *off + 1;
	  }
	  
	} /* end ivar loop */
      } /* end i loop */
    } /* end j loop */
  } /* end k loop */ 
  
}

void READ_Blocks_hdf5_r8(hid_t* file_id, 
			 hid_t* dset_block, 
			 hid_t* block_type_id, 
			 int* tot_blocks, 
			 int* lnblocks, 
			 int* mdim,
			 int* ndim,
			 int* ngid,
			 int* mflags,
			 int  lrefine[],
			 int  nodetype[],
			 int  which_child[],
			 int  gid[],
			 int  bflags[],
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
                         int   iorder[])
{
  hid_t  memspace, filespace;
  hsize_t  count[1];
  hssize_t offset[1];
  hid_t plist_id;
  int   ivar, i, j, k, lb, i2, j2, k2, nx2, ny2, ivar2;
  int   istart, iend, jstart, jend, kstart, kend;
  herr_t status;
  int mpi_rank, mpi_size;
  int n_to_left;
  int icount, proc;
  int k2d, k3d;
  int size_block_type, off;
  float * block_pointer;

  k2d = 0;
  if (*ndim >= 2) k2d = 1;
  k3d = 0;
  if (*ndim == 3) k3d = 1;

  size_block_type = (sizeof(float)*3) +
    (sizeof(float) * *ngid) * + (sizeof(float) * *mflags) +
    (sizeof(float) * *mdim) + (sizeof(float) * 2 * *mdim) +
    (sizeof(float)) +
    (sizeof(float) * *nvar_chk_cc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_fc * (*iu0-*il0) * (*ju0-*jl0) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0) * (*ju0-*jl0+k2d) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0) * (*ku0-*kl0+k3d)) +
    (sizeof(float) * *nvar_chk_ec * (*iu0-*il0+1) * (*ju0-*jl0+k2d) * (*ku0-*kl0)) +
    (sizeof(float) * *nvar_chk_nc * (*iu0-*il0+1) * (*ju0-*jl0+k2d) * (*ku0-*kl0+k3d));

  block_pointer = (float *) malloc(size_block_type);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
  icount = 0;
  *lnblocks = 0;
  n_to_left = 0;
  proc = 0;
  while (icount < *tot_blocks) {
    if (proc == mpi_rank) {
      *lnblocks = *lnblocks + 1;
    } 
    if (proc < mpi_rank) { 
      n_to_left = n_to_left + 1; 
    } 
    icount = icount + 1;
    proc = proc + 1;
    if (proc > mpi_size-1) { proc = 0; }
  }


  nx2 = *iu0 - *il0;
  ny2 = *ju0 - *jl0;
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

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  /*H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);*/
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  for (lb = 0; lb < *lnblocks; lb++) {

    count[0] = 1;
    memspace = H5Screate_simple(1, count, NULL);
    
    offset[0] = n_to_left + lb; 

    filespace = H5Dget_space(*dset_block);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			NULL, count, NULL);
    

    status = H5Dread(*dset_block, *block_type_id, memspace, filespace, 
		     plist_id, block_pointer);
    
    off = 0;
    lrefine[lb] = (int)*(block_pointer+off);
    off = off + 1;
    nodetype[lb] = (int)*(block_pointer+off);
    off = off + 1;
    which_child[lb] = (int)*(block_pointer+off);
    off = off + 1;
    for (i = 0; i < *ngid; i++) {
      gid[(*ngid * lb) + i] = (int)*(block_pointer+off);
      off = off + 1;
    }
    for (i = 0; i < *mflags; i++) {
      bflags[(*mflags * lb) + i] = (int)*(block_pointer+off);
      off = off + 1;
    }
    for (i = 0; i < *mdim; i++) {
      coord[(*mdim * lb) + i] = *(block_pointer+off);
      off = off + 1;
    }
    for (i = 0; i < 2 * *mdim; i++) {
      bnd_box[(2 * *mdim * lb) + i] = *(block_pointer+off);
      off = off + 1;
    }
    work_block[lb] = *(block_pointer+off);
    off = off + 1;

    /* unk data */

    if (*nvar_chk_cc > 0) {

      var_to_block_pointer(unk, iorder, *nvar, *nx, *ny, *nz,
			   kstart, kend,
			   jstart, jend,
			   istart, iend,
			   checkp_on_cc, 0, 0, lb,
			   &off, block_pointer);
  
    } /* end if nvar_chk_cc */

    /* facevarx data */
    
    if (*nvar_chk_fc > 0) {

      var_to_block_pointer(facevarx, iorder, *nbndvar, *nx+1, *ny, *nz,
			   kstart, kend,
			   jstart, jend,
			   istart, iend+1,
			   checkp_on_fc, 2, 0, lb,
			   &off, block_pointer);

    /* facevary data */

      var_to_block_pointer(facevary, iorder, *nbndvar, *nx, *ny+k2d, *nz,
			   kstart, kend,
			   jstart, jend+k2d,
			   istart, iend,
			   checkp_on_fc, 2, 1, lb,
			   &off, block_pointer);

    /* facevarz data */

      var_to_block_pointer(facevarz, iorder, *nbndvar, *nx, *ny, *nz+k3d,
			   kstart, kend+k3d,
			   jstart, jend,
			   istart, iend,
			   checkp_on_fc, 2, 2, lb,
			   &off, block_pointer);

    } /* end if (nvar_chk_fc */
    
    /* unk_e_x data */
    
    if (*nvar_chk_ec > 0) {

      var_to_block_pointer(unk_e_x, iorder, *nbndvare, *nx, *ny+k2d, *nz+k3d,
			   kstart, kend+k3d,
			   jstart, jend+k2d,
			   istart, iend,
			   checkp_on_ec, 2, 0, lb,
			   &off, block_pointer);

    /* unk_e_y data */

      var_to_block_pointer(unk_e_y, iorder, *nbndvare, *nx+1, *ny, *nz+k3d,
			   kstart, kend+k3d,
			   jstart, jend,
			   istart, iend+1,
			   checkp_on_ec, 2, 1, lb,
			   &off, block_pointer);
    
    /* unk_e_z data */

      var_to_block_pointer(unk_e_z, iorder, *nbndvare, *nx+1, *ny+k2d, *nz,
			   kstart, kend,
			   jstart, jend+k2d,
			   istart, iend+1,
			   checkp_on_ec, 2, 2, lb,
			   &off, block_pointer);
    
    } /* end if (nvar_chk_ec */

    /* unk_n data */

    if (*nvar_chk_nc > 0) {

      var_to_block_pointer(unk_n, iorder, *nbndvarc, *nx+1, *ny+k2d, *nz+k3d,
			   kstart, kend+k3d,
			   jstart, jend+k2d,
			   istart, iend+1,
			   checkp_on_nc, 0, 0, lb,
			   &off, block_pointer);

    } /* end if nvar_chk_nc */
    
    H5Sclose(filespace);
    H5Sclose(memspace);
    
  }

  H5Pclose(plist_id);
  free(block_pointer);
  
}


#ifdef UNDERSCORE
void read_blocks_hdf5_r8_ (char* file_name_in, 
#endif
#ifdef DOUBLE_UNDERSCORE
void read_blocks_hdf5_r8__ (char* file_name_in, 
#endif
#ifndef UNDERSCORE 
#ifndef DOUBLE_UNDERSCORE
void read_blocks_hdf5_r8 (char* file_name_in, 
#endif
#endif
			  int*  tot_blocks,
			  int*  lnblocks,
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
			  int*  nx, int* ny, int* nz,
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
  hid_t  file_id, dset_block;  /* file and dataset ids */
  hid_t  block_type_id;
  int ivar;

  OPEN_File_hdf5_r8(file_name_in, &file_id);

  READ_Header_hdf5_r8(&file_id, tot_blocks,
		      user_attr_1,
		      user_attr_2,
		      user_attr_3,
		      user_attr_4,
		      user_attr_5);
  
  CREATE_Block_data_type_hdf5_r8(&block_type_id,
				 mdim,
				 ndim,
				 ngid,
				 mflags,
				 nvar_chk_cc,
				 nvar_chk_fc,
				 nvar_chk_ec,
				 nvar_chk_nc,
				 il0, iu0,
				 jl0, ju0,
				 kl0, ku0);

  OPEN_Block_data_set_hdf5_r8(&dset_block, 
			      &file_id);
  
  READ_Blocks_hdf5_r8(&file_id, 
		      &dset_block, 
		      &block_type_id, 
		      tot_blocks, 
		      lnblocks, 
		      mdim,
		      ndim,
		      ngid,
		      mflags,
		      lrefine,
		      nodetype,
		      which_child,
		      gid,
		      bflags,
		      coord,
		      bnd_box,
		      work_block,
		      unk,
		      nvar,
		      nvar_chk_cc,
		      checkp_on_cc,
		      facevarx, facevary, facevarz,
		      nbndvar,
		      nvar_chk_fc,
		      checkp_on_fc,
		      unk_e_x, unk_e_y, unk_e_z,
		      nbndvare,
		      nvar_chk_ec,
		      checkp_on_ec,
		      unk_n,
		      nbndvarc,
		      nvar_chk_nc,
		      checkp_on_nc,
		      nx,  ny, nz,
		      il0, iu0,
		      jl0, ju0,
		      kl0, ku0,
		      iorder);
  
  /* Close resources */
  H5Dclose(dset_block);
  H5Fclose(file_id);
  H5Tclose(block_type_id);

}
