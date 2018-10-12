#include "hdf5.h"
#include "mpi.h"
#include <math.h>
#include <ctype.h>

#include "underscore.h"

void get_udim_chombo_r4(int nvar, int nx, int ny, int nz,
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

void CREATE_File4 (char* file_name_in, hid_t* file_id) {

  hid_t plist_id;
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

  /* Set up file access property list with parllel I/O access */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  /* Create a new file collectively */
  *file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

}



void WRITE_Header4(hid_t*  file_id, 
		   int*    num_components,
		   int*    num_levels,
		   char*   compNames,
		   int*    ndim,
		   int*    mdim,
		   int*    nxb,
		   int*    nyb,
		   int*    nzb,
                   int*    lnblocks,
		   float   bnd_box[],
		   int     no_at_level[]) {

  hsize_t  dims[1];
  hid_t    attribute_id, dataspace_id, group_id, group_id2;
  hid_t    string_type, prob_domain_type_id;
  hid_t    dset_id;
  herr_t   status;
  int      mpi_rank;
  int      temp[1];
  float    rtemp[1];
  char     ctemp[20];
  char     name[20];
  int      i, ii;
  float    dx, dy, dz;
  float xleft, yleft, zleft;
  float xright, yright, zright;
  float xl, yl, zl;
  int lb;

  typedef struct prob_domain_type_3d {
    int lo_i;
    int lo_j;
    int lo_k;
    int hi_i;
    int hi_j;
    int hi_k;
  } prob_domain_type_3d;
  prob_domain_type_3d prob_domain_3d[1];

  typedef struct prob_domain_type_2d {
    int lo_i;
    int lo_j;
    int hi_i;
    int hi_j;
  } prob_domain_type_2d;
  prob_domain_type_2d prob_domain_2d[1];

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  group_id = H5Gopen(*file_id, "/");


  /* num_components */
  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "num_components", 
			   H5T_NATIVE_INT, dataspace_id, 
			   H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, num_components);
  }
  H5Sclose(dataspace_id);
  H5Aclose(attribute_id);


  /* iteration, this is a dummy */
  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "iteration", 
			   H5T_NATIVE_INT, dataspace_id, 
			   H5P_DEFAULT);
  temp[0] = 0;
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, temp);
  }
  H5Sclose(dataspace_id);
  H5Aclose(attribute_id);


  /* num_levels */
  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "num_levels", 
			   H5T_NATIVE_INT, dataspace_id, 
			   H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, num_levels);
  }
  H5Sclose(dataspace_id);
  H5Aclose(attribute_id);

  /* anisotropic stretch's x, y, z */
  typedef struct floatvect_id_type {
    float x;
    float y;
    float z;
  }  floatvect_id_type;
  floatvect_id_type floatvect[1];
  hid_t floatvect_id;
  floatvect_id = H5Tcreate (H5T_COMPOUND, sizeof(floatvect_id_type));

  H5Tinsert(floatvect_id, "x",  
            HOFFSET(floatvect_id_type, x), H5T_NATIVE_FLOAT);
  H5Tinsert(floatvect_id, "y",  
            HOFFSET(floatvect_id_type, y), H5T_NATIVE_FLOAT);
  H5Tinsert(floatvect_id, "z",  
            HOFFSET(floatvect_id_type, z), H5T_NATIVE_FLOAT);

  /* Find the maximum (right) sides of the domain */

  xl = -1.e30;
  yl = -1.e30;
  zl = -1.e30;
  for (lb = 0; lb < *lnblocks; lb++) {
    if (bnd_box[lb* *mdim*2 + 2*0 + 1] > xl) {
      xl = bnd_box[lb* *mdim*2 + 2*0 + 1];
    }
    if (bnd_box[lb* *mdim*2 + 2*1 + 1] > yl) {
      yl = bnd_box[lb* *mdim*2 + 2*1 + 1];
    }
    if (bnd_box[lb* *mdim*2 + 2*2 + 1] > zl) {
      zl = bnd_box[lb* *mdim*2 + 2*2 + 1];
    }
  }

  status = MPI_Allreduce(&xl, &xright, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  status = MPI_Allreduce(&yl, &yright, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  status = MPI_Allreduce(&zl, &zright, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

  /* Find the minimum (left) sides of the domain */

  xl = 1.e30;
  yl = 1.e30;
  zl = 1.e30;
  for (lb = 0; lb < *lnblocks; lb++) {
    if (bnd_box[lb* *mdim*2 + 2*0 + 0] < xl) {
      xl = bnd_box[lb* *mdim*2 + 2*0 + 0];
    }
    if (bnd_box[lb* *mdim*2 + 2*1 + 0] < yl) {
      yl = bnd_box[lb* *mdim*2 + 2*1 + 0];
    }
    if (bnd_box[lb* *mdim*2 + 2*2 + 0] < zl) {
      zl = bnd_box[lb* *mdim*2 + 2*2 + 0];
    }
  }

  status = MPI_Allreduce(&xl, &xleft, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  status = MPI_Allreduce(&yl, &yleft, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  status = MPI_Allreduce(&zl, &zleft, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

  floatvect[0].x = 1.;
  floatvect[0].y = (yright-yleft) / (xright-xleft);
  floatvect[0].y = floatvect[0].y * *nxb / *nyb;
  floatvect[0].z = (zright-zleft) / (xright-xleft);
  floatvect[0].z = floatvect[0].z * *nxb / *nzb;

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "anisotropic", 
			   floatvect_id, dataspace_id, 
			   H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, floatvect_id, floatvect);
  }
  H5Sclose(dataspace_id);
  H5Aclose(attribute_id);


  /* time, this is a dummy */
  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "time", 
			   H5T_NATIVE_FLOAT, dataspace_id, 
			   H5P_DEFAULT);
  rtemp[0] = 0;
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, rtemp);
  }
  H5Sclose(dataspace_id);
  H5Aclose(attribute_id);


  /* compNames */
  for (i = 0; i < *num_components; i++) {
    dataspace_id = H5Screate(H5S_SCALAR);
    string_type  = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, 20);
    sprintf(name, "component_%i",i);
    attribute_id = H5Acreate(group_id, name, 
			     string_type, dataspace_id, 
			     H5P_DEFAULT);
    for (ii = 0; ii < 20; ii++) {
      ctemp[ii] = compNames[i*20+ii];
    }
    if (mpi_rank == 0) {
      status = H5Awrite(attribute_id, string_type, ctemp);
    }
    H5Sclose(dataspace_id);
    H5Aclose(attribute_id);
    H5Tclose(string_type);
  }

  /* close group "/" */
  H5Gclose(group_id);

  /* Group Chombo global */
  group_id = H5Gcreate(*file_id, "/Chombo_global", 0);

  /* Add attribute "SpaceDim" to group "/Chombo_global" */
  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "SpaceDim", 
			   H5T_NATIVE_INT, dataspace_id, 
			   H5P_DEFAULT);
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, ndim);
  }
  H5Sclose(dataspace_id);
  H5Aclose(attribute_id);

  /* Add attribute "testReal" to group "/Chombo_global" */
  /* This is a dummy */
  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "testReal", 
			   H5T_NATIVE_FLOAT, dataspace_id, 
			   H5P_DEFAULT);
  rtemp[0] = 0;
  if (mpi_rank == 0) {
    status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, rtemp);
  }
  H5Sclose(dataspace_id);
  H5Aclose(attribute_id);

  /* close group "/Chombo_global" */
  H5Gclose(group_id);

  /* Create Groups for each level */
  for (i = 0; i < *num_levels; i++) {

    sprintf(name, "/level_%i",i);
    group_id = H5Gcreate(*file_id, name, 0);
    
    /* Add attribute for ref_ratio */
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate(group_id, "ref_ratio", 
			     H5T_NATIVE_INT, dataspace_id, 
			     H5P_DEFAULT);
    temp[0] = 2;
    if (mpi_rank == 0) {
      status = H5Awrite(attribute_id, H5T_NATIVE_INT, temp);
    }
    H5Sclose(dataspace_id);
    H5Aclose(attribute_id);
    
    /* Add attribute for "dt" */
    /* This is a dummy */
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate(group_id, "dt", 
			     H5T_NATIVE_FLOAT, dataspace_id, 
			     H5P_DEFAULT);
    rtemp[0] = 0.;
    if (mpi_rank == 0) {
      status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, rtemp);
    }
    H5Sclose(dataspace_id);
    H5Aclose(attribute_id);
    
    /* Add attribute for "dx" */
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate(group_id, "dx", 
			     H5T_NATIVE_FLOAT, dataspace_id, 
			     H5P_DEFAULT);
    /* This is not really right, but things should look OK 
       as long as the cell sizes change by only a factor of 2
       from one level to the next
       This will change the total domain size 
       What we really want is a 1 x 1 x 1 box in all cases, but relative
       sizes should be OK .
       Maybe its best just to assume that the top level is just 1 block ??? */
    dx = (xright-xleft)/(pow(2,i)*(float) *nxb);
    if (mpi_rank == 0) {
      status = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &dx);
    }
    H5Sclose(dataspace_id);
    H5Aclose(attribute_id);

    /* Add attribute for "prob_domain" */

    if (*ndim == 3) {

    prob_domain_type_id = H5Tcreate (H5T_COMPOUND, sizeof(prob_domain_type_3d));

    H5Tinsert(prob_domain_type_id, "lo_i",  
	      HOFFSET(prob_domain_type_3d, lo_i), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "lo_j",  
	      HOFFSET(prob_domain_type_3d, lo_j), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "lo_k",  
	      HOFFSET(prob_domain_type_3d, lo_k), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "hi_i",  
	      HOFFSET(prob_domain_type_3d, hi_i), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "hi_j",  
	      HOFFSET(prob_domain_type_3d, hi_j), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "hi_k",  
	      HOFFSET(prob_domain_type_3d, hi_k), H5T_NATIVE_INT);

    dims[0] = 1;
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate(group_id, "prob_domain", 
			     prob_domain_type_id, dataspace_id, 
			     H5P_DEFAULT);

    prob_domain_3d[0].lo_i = 0;
    prob_domain_3d[0].lo_j = 0;
    prob_domain_3d[0].lo_k = 0;
    prob_domain_3d[0].hi_i = (pow(2,i) * *nxb) - 1;
    prob_domain_3d[0].hi_j = (pow(2,i) * *nyb) - 1;
    prob_domain_3d[0].hi_k = (pow(2,i) * *nzb) - 1;
    
    if (mpi_rank == 0) {
      status = H5Awrite(attribute_id, prob_domain_type_id, prob_domain_3d);
    }
    H5Sclose(dataspace_id);
    H5Aclose(attribute_id);

    } else {

    prob_domain_type_id = H5Tcreate (H5T_COMPOUND, sizeof(prob_domain_type_2d));

    H5Tinsert(prob_domain_type_id, "lo_i",  
	      HOFFSET(prob_domain_type_2d, lo_i), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "lo_j",  
	      HOFFSET(prob_domain_type_2d, lo_j), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "hi_i",  
	      HOFFSET(prob_domain_type_2d, hi_i), H5T_NATIVE_INT);
    H5Tinsert(prob_domain_type_id, "hi_j",  
	      HOFFSET(prob_domain_type_2d, hi_j), H5T_NATIVE_INT);

    dims[0] = 1;
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate(group_id, "prob_domain", 
			     prob_domain_type_id, dataspace_id, 
			     H5P_DEFAULT);

    prob_domain_2d[0].lo_i = 0;
    prob_domain_2d[0].lo_j = 0;
    prob_domain_2d[0].hi_i = (pow(2,i) * *nxb) - 1;
    prob_domain_2d[0].hi_j = (pow(2,i) * *nyb) - 1;
    
    if (mpi_rank == 0) {
      status = H5Awrite(attribute_id, prob_domain_type_id, prob_domain_2d);
    }
    H5Sclose(dataspace_id);
    H5Aclose(attribute_id);

    }

    /* Add "boxes" dataset */
    dims[0] = no_at_level[i];
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate(group_id, "boxes", prob_domain_type_id, 
			dataspace_id, H5P_DEFAULT);
    H5Dclose(dset_id);
    H5Sclose(dataspace_id);

    H5Tclose(prob_domain_type_id);

    /* Add "data:datatype=0" dataset */
    dims[0] = *num_components * *nxb * *nyb * *nzb * no_at_level[i];
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate(group_id, "data:datatype=0", H5T_NATIVE_FLOAT, 
			dataspace_id, H5P_DEFAULT);
    H5Dclose(dset_id);
    H5Sclose(dataspace_id);

    /* Group "data_attributes" */
    group_id2 = H5Gcreate(group_id, "data_attributes", 0);
    
    /* Add attribute "comps" to group "data_attributes" */
    dataspace_id = H5Screate(H5S_SCALAR);
    attribute_id = H5Acreate(group_id2, "comps", 
			     H5T_NATIVE_INT, dataspace_id, 
			     H5P_DEFAULT);
    if (mpi_rank == 0) {
      status = H5Awrite(attribute_id, H5T_NATIVE_INT, num_components);
    }
    H5Sclose(dataspace_id);
    H5Aclose(attribute_id);

    /* close group "data_attributes" */
    H5Gclose(group_id2);

    /* close the group for level i */
    H5Gclose(group_id);

  } /* end loop over levels */

}


void WRITE_Blocks4(hid_t* file_id, 
		   int* num_levels,
		   int* num_components,
		   int* nxb,
		   int* nyb,
		   int* nzb,
		   int* max_lnblocks, 
		   int* n_to_left, 
		   int  n_to_left_level[], 
		   int* lnblocks, 
		   int* mdim,
		   int* ndim,
		   int* mflags,
		   int lrefine[], 
		   int nodetype[], 
		   int which_child[], 
		   int bflags[],
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
  int rank_array;
  hsize_t  dims[1], dims_array[1];
  hssize_t offset[1];
  hid_t plist_id;
  hid_t dset_id, dataspace_id, group_id, prob_domain_type_id;
  int ivar, i, j, k, lb, i2, j2, k2, nx2, ny2, ivar2;
  int istart, iend, jstart, jend, kstart, kend, k3d;
  herr_t status;
  char name[20];
  int  n_at_level[*num_levels];
  int  block_no_at_level[*lnblocks];
  float xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz;
  int lb2;
  float unk_t[*nxb * *nyb * *nzb];
  int mpi_rank;
  int ii, jj, kk;
  float xl, yl, zl;
  float xleft, yleft, zleft;


  typedef struct prob_domain_type_3d {
    int lo_i;
    int lo_j;
    int lo_k;
    int hi_i;
    int hi_j;
    int hi_k;
  } prob_domain_type_3d;
  prob_domain_type_3d boxes_3d[1];
  
  typedef struct prob_domain_type_2d {
    int lo_i;
    int lo_j;
    int hi_i;
    int hi_j;
  } prob_domain_type_2d;
  prob_domain_type_2d boxes_2d[1];
  

  /* Create property list for collective dataset write */
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  /*H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);*/
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  /* Create "prob_domain_id" */

  if (*ndim == 3) {

  prob_domain_type_id = H5Tcreate (H5T_COMPOUND, sizeof(prob_domain_type_3d));
  
  H5Tinsert(prob_domain_type_id, "lo_i",  
	    HOFFSET(prob_domain_type_3d, lo_i), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "lo_j",  
	    HOFFSET(prob_domain_type_3d, lo_j), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "lo_k",  
	    HOFFSET(prob_domain_type_3d, lo_k), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "hi_i",  
	    HOFFSET(prob_domain_type_3d, hi_i), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "hi_j",  
	    HOFFSET(prob_domain_type_3d, hi_j), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "hi_k",  
	    HOFFSET(prob_domain_type_3d, hi_k), H5T_NATIVE_INT);

  } else {

  prob_domain_type_id = H5Tcreate (H5T_COMPOUND, sizeof(prob_domain_type_2d));
  
  H5Tinsert(prob_domain_type_id, "lo_i",  
	    HOFFSET(prob_domain_type_2d, lo_i), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "lo_j",  
	    HOFFSET(prob_domain_type_2d, lo_j), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "hi_i",  
	    HOFFSET(prob_domain_type_2d, hi_i), H5T_NATIVE_INT);
  H5Tinsert(prob_domain_type_id, "hi_j",  
	    HOFFSET(prob_domain_type_2d, hi_j), H5T_NATIVE_INT);

  }
  
  nx2 = *iu0 - *il0;
  ny2 = *ju0 - *jl0;
  kstart = *kl0;
  kend   = *ku0;
  jstart = *jl0;
  jend   = *ju0;
  istart = *il0;
  iend   = *iu0;
  k3d    = 1;
  if (*ndim < 3) {
    kstart = 0;
    kend   = 1;
    k3d    = 0;
  }
  if (*ndim < 2) {
    kstart = 0;
    kend   = 1;
    jstart = 0;
    jend   = 1;
    k3d    = 0;
  }
  
  /* compute no of blocks at each level */
  for (i = 0; i < *num_levels; i++) {
    n_at_level[i] = 0;
  }
  for (i = 0; i < *lnblocks; i++) {
    block_no_at_level[i] = 0;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  for (lb = 0; lb < *lnblocks; lb++) {
    block_no_at_level[lb] = n_to_left_level[lrefine[lb]-1] + 
                            n_at_level[lrefine[lb]-1];
    n_at_level[lrefine[lb]-1] = n_at_level[lrefine[lb]-1] + 1;
  }

  /* Find the minimum (left) sides of the domain */

  xl = 1.e30;
  yl = 1.e30;
  zl = 1.e30;
  for (lb = 0; lb < *lnblocks; lb++) {
    if (bnd_box[lb* *mdim*2 + 2*0 + 0] < xl) {
      xl = bnd_box[lb* *mdim*2 + 2*0 + 0];
    }
    if (bnd_box[lb* *mdim*2 + 2*1 + 0] < yl) {
      yl = bnd_box[lb* *mdim*2 + 2*1 + 0];
    }
    if (bnd_box[lb* *mdim*2 + 2*2 + 0] < zl) {
      zl = bnd_box[lb* *mdim*2 + 2*2 + 0];
    }
  }

  status = MPI_Allreduce(&xl, &xleft, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  status = MPI_Allreduce(&yl, &yleft, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  status = MPI_Allreduce(&zl, &zleft, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

  for (lb = 0; lb < *max_lnblocks; lb++) { 

    if (lb < *lnblocks) {
      lb2 = lb;
    } else {
      lb2 = *lnblocks - 1;
    }

    /* Open Group for this level */
    sprintf(name, "/level_%i",lrefine[lb2]-1);
    group_id = H5Gopen(*file_id, name);


    /* Open "boxes" dataset */
    dims[0] = 1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dopen(group_id, "boxes");
    
    offset[0] = block_no_at_level[lb2]; 
    count[0] = 1;
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			NULL, count, NULL);
    
    dx      = bnd_box[lb2* *mdim*2 + 2*0 + 1] - bnd_box[lb2* *mdim*2 + 2*0 + 0];
    dx      = dx/(float) *nxb;
    dy      = bnd_box[lb2* *mdim*2 + 2*1 + 1] - bnd_box[lb2* *mdim*2 + 2*1 + 0];
    dy      = dy/(float) *nyb;
    dz      = bnd_box[lb2* *mdim*2 + 2*2 + 1] - bnd_box[lb2* *mdim*2 + 2*2 + 0];
    dz      = dz/(float) *nzb;
    xmin    = bnd_box[lb2* *mdim*2 + 2*0 + 0] + dx/2.;
    xmax    = bnd_box[lb2* *mdim*2 + 2*0 + 1] - dx/2.;
    ymin    = bnd_box[lb2* *mdim*2 + 2*1 + 0] + dy/2.;
    ymax    = bnd_box[lb2* *mdim*2 + 2*1 + 1] - dy/2.;
    zmin    = bnd_box[lb2* *mdim*2 + 2*2 + 0] + dz/2.;
    zmax    = bnd_box[lb2* *mdim*2 + 2*2 + 1] - dz/2.;

    if (*ndim == 3) {

    boxes_3d[0].lo_i = ((int)((xmin - xleft)/dx));
    boxes_3d[0].lo_j = ((int)((ymin - yleft)/dy));
    boxes_3d[0].lo_k = ((int)((zmin - zleft)/dz));
    boxes_3d[0].hi_i = ((int)((xmax - xleft)/dx));
    boxes_3d[0].hi_j = ((int)((ymax - yleft)/dy));
    boxes_3d[0].hi_k = ((int)((zmax - zleft)/dz));

    status = H5Dwrite(dset_id, prob_domain_type_id, dataspace_id, filespace,
		      plist_id, boxes_3d);

    } else {

    boxes_2d[0].lo_i = ((int)((xmin - xleft)/dx));
    boxes_2d[0].lo_j = ((int)((ymin - yleft)/dy));
    boxes_2d[0].hi_i = ((int)((xmax - xleft)/dx));
    boxes_2d[0].hi_j = ((int)((ymax - yleft)/dy));

    status = H5Dwrite(dset_id, prob_domain_type_id, dataspace_id, filespace,
		      plist_id, boxes_2d);

    }

    /* close group data space and data */
    H5Sclose(dataspace_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);

    status = MPI_Barrier(MPI_COMM_WORLD);

    /*****************************************************************/

    /* Cell Centered Data */
    
    count[0]   = *nxb * *nyb * *nzb;
    offset[0]  = block_no_at_level[lb2] * count[0] * *num_components; 

    if (*nvar_chk_cc > 0) {
      
      ivar2 = 0;
      for (ivar=0; ivar < *nvar; ivar++) {

	if (checkp_on_cc[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0];
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		int udim[4]; int it[4];
		get_udim_chombo_r4 (*nvar, *nx, *ny, *nz, ivar, i, j, k, iorder, udim, it);

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
		      (udim[2] * udim[1] * udim[0] * it[3]) + 
		      (udim[1] * udim[0] * it[2]) + 
		      (udim[0] * it[1]) + 
		      it[0]];

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */


	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

    } /* end if (nvar_chk_cc */

    /* Face Centered Data */


    if (*nvar_chk_fc > 0) {
      
      /* facevarx */

      ivar2 = 0;
      for (ivar=0; ivar < *nbndvar; ivar++) {
	if (checkp_on_fc[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0]; 
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		int udim[4]; int it[4]; int it2[4];
		get_udim_chombo_r4 (*nbndvar, *nx+1, *ny, *nz, ivar, i, j, k, iorder, udim, it);
		get_udim_chombo_r4 (*nbndvar, *nx+1, *ny, *nz, ivar, i+1, j, k, iorder, udim, it2);

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  facevarx[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			   (udim[2] * udim[1] * udim[0] * it[3]) + 
			   (udim[1] * udim[0] * it[2]) + 
			   (udim[0] * it[1]) + 
			   it[0]] +
		  facevarx[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			   (udim[2] * udim[1] * udim[0] * it2[3]) + 
			   (udim[1] * udim[0] * it2[2]) + 
			   (udim[0] * it2[1]) + 
			   it2[0]];

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] / 2.;

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */

	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

      /* facevary */

      ivar2 = 0;
      for (ivar=0; ivar < *nbndvar; ivar++) {
	if (checkp_on_fc[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0]; 
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		int udim[4]; int it[4]; int it2[4];
		get_udim_chombo_r4 (*nbndvar, *nx, *ny+1, *nz, ivar, i, j, k, iorder, udim, it);
		get_udim_chombo_r4 (*nbndvar, *nx, *ny+1, *nz, ivar, i, j+1, k, iorder, udim, it2);

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  facevary[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			   (udim[2] * udim[1] * udim[0] * it[3]) + 
			   (udim[1] * udim[0] * it[2]) + 
			   (udim[0] * it[1]) + 
			   it[0]] +
		  facevary[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			   (udim[2] * udim[1] * udim[0] * it2[3]) + 
			   (udim[1] * udim[0] * it2[2]) + 
			   (udim[0] * it2[1]) + 
			   it2[0]];

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] / 2.;

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */

	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

      /* facevarz */

      ivar2 = 0;
      for (ivar=0; ivar < *nbndvar; ivar++) {
	if (checkp_on_fc[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0]; 
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		int udim[4]; int it[4]; int it2[4];
		get_udim_chombo_r4 (*nbndvar, *nx, *ny, *nz+k3d, ivar, i, j, k, iorder, udim, it);
		get_udim_chombo_r4 (*nbndvar, *nx, *ny, *nz+k3d, ivar, i, j, k+k3d, iorder, udim, it2);

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  facevarz[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			   (udim[2] * udim[1] * udim[0] * it[3]) + 
			   (udim[1] * udim[0] * it[2]) + 
			   (udim[0] * it[1]) + 
			   it[0]] +
		  facevarz[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			   (udim[2] * udim[1] * udim[0] * it2[3]) + 
			   (udim[1] * udim[0] * it2[2]) + 
			   (udim[0] * it2[1]) + 
			   it2[0]];

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] / 2.;

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */

	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

    } /* end if (nvar_chk_fc */



    /* Edge Centered Data */



    if (*nvar_chk_ec > 0) {
      
      /* unk_e_x */

      ivar2 = 0;
      for (ivar=0; ivar < *nbndvare; ivar++) {
	if (checkp_on_ec[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0]; 
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] = 0.;

		for (kk = 0; kk < 2; kk++) {
		for (jj = 0; jj < 2; jj++) {

		int udim[4]; int it[4];
		get_udim_chombo_r4 (*nbndvare, *nx, *ny+1, *nz+k3d, ivar, i, j+jj, k+kk*k3d, iorder, udim, it);

		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] =
		    unk_t[i2 + 
			  nx2 * j2 + 
			  nx2 * ny2 * k2] +
		    unk_e_x[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			    (udim[2] * udim[1] * udim[0] * it[3]) + 
			    (udim[1] * udim[0] * it[2]) + 
			    (udim[0] * it[1]) + 
			    it[0]];
		}
	        }

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] / 4.;

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */

	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

      /* unk_e_y */

      ivar2 = 0;
      for (ivar=0; ivar < *nbndvare; ivar++) {
	if (checkp_on_ec[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0]; 
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] = 0.;

		for (kk = 0; kk < 2; kk++) {
		for (ii = 0; ii < 2; ii++) {

		int udim[4]; int it[4];
		get_udim_chombo_r4 (*nbndvare, *nx+1, *ny, *nz+k3d, ivar, i+ii, j, k+kk*k3d, iorder, udim, it);

		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] =
		    unk_t[i2 + 
			  nx2 * j2 + 
			  nx2 * ny2 * k2] +
		    unk_e_y[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			    (udim[2] * udim[1] * udim[0] * it[3]) + 
			    (udim[1] * udim[0] * it[2]) + 
			    (udim[0] * it[1]) + 
			    it[0]];
		}
	        }

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] / 4.;

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */

	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

      /* unk_e_z */

      ivar2 = 0;
      for (ivar=0; ivar < *nbndvare; ivar++) {
	if (checkp_on_ec[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0]; 
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] = 0.;

		for (jj = 0; jj < 2; jj++) {
		for (ii = 0; ii < 2; ii++) {

		int udim[4]; int it[4];
		get_udim_chombo_r4 (*nbndvare, *nx+1, *ny+1, *nz, ivar, i+ii, j+jj, k, iorder, udim, it);

		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] =
		    unk_t[i2 + 
			  nx2 * j2 + 
			  nx2 * ny2 * k2] +
		    unk_e_z[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			    (udim[2] * udim[1] * udim[0] * it[3]) + 
			    (udim[1] * udim[0] * it[2]) + 
			    (udim[0] * it[1]) + 
			    it[0]];
		}
	        }

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] / 4.;

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */

	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

    } /* end if (nvar_chk_ec */


    /* Node Centered Data */


    if (*nvar_chk_nc > 0) {

      /* unk_n */

      ivar2 = 0;
      for (ivar=0; ivar < *nbndvarc; ivar++) {
	if (checkp_on_nc[ivar] == 1) {

	  /* Open "data:datatype=0" dataset */
	  dims[0] = *nxb * *nyb * *nzb;
	  dataspace_id = H5Screate_simple(1, dims, NULL);
	  dset_id = H5Dopen(group_id, "data:datatype=0");
	  
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
			      NULL, count, NULL);

	  offset[0] = offset[0] + count[0]; 
    
	  k2 = 0;
	  for (k = kstart; k < kend; k++) {
	    j2 = 0;
	    for (j = jstart; j < jend; j++) {
	      i2 = 0;
	      for (i = istart; i < iend; i++) {
		
		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] = 0.;

		for (kk = 0; kk < 2; kk++) {
		for (jj = 0; jj < 2; jj++) {
		for (ii = 0; ii < 2; ii++) {

		int udim[4]; int it[4];
		get_udim_chombo_r4 (*nbndvarc, *nx+1, *ny+1, *nz+k3d, ivar, i+ii, j+jj, k+kk*k3d, iorder, udim, it);

		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] =
		    unk_t[i2 + 
			  nx2 * j2 + 
			  nx2 * ny2 * k2] +
		    unk_n[(udim[3] * udim[2] * udim[1] * udim[0] * lb2) + 
			  (udim[2] * udim[1] * udim[0] * it[3]) + 
			  (udim[1] * udim[0] * it[2]) + 
			  (udim[0] * it[1]) + 
			  it[0]];
		}
	        }
		}

		unk_t[i2 + 
		      nx2 * j2 + 
		      nx2 * ny2 * k2] =
		  unk_t[i2 + 
			nx2 * j2 + 
			nx2 * ny2 * k2] / 8.;

		i2 = i2 + 1;
	      } /* end i loop */
	      j2 = j2 + 1;
	    } /* end j loop */
	    k2 = k2 + 1;
	  } /* end k loop */

	  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, filespace, 
			    plist_id, unk_t);

	  H5Sclose(dataspace_id);
	  H5Dclose(dset_id);
	  H5Sclose(filespace);

	  ivar2 = ivar2 + 1;
	} /* end if (checkp */
      } /* end ivar loop */

    } /* end if (nvar_chk_nc */

    /* close group data space and data */
    H5Gclose(group_id);

  } /* end loop over blocks */

  H5Pclose(plist_id);
  H5Tclose(prob_domain_type_id);

}

#ifdef UNDERSCORE
void write_blocks_chombo_r4_ (char* file_name_in, 
#endif
#ifdef DOUBLE_UNDERSCORE
void write_blocks_chombo_r4__ (char* file_name_in, 
#endif
#ifndef UNDERSCORE 
#ifndef DOUBLE_UNDERSCORE
void write_blocks_chombo_r4 (char* file_name_in, 
#endif
#endif
			   int*  num_components,
			   int*  num_levels,
			   char* compNames,
                           int*  ndim,
                           int*  nxb,
			   int*  nyb,
			   int*  nzb,
			   int   no_at_level[],
  			   int*  tot_blocks,
			   int*  tot_blocks_wr,
			   int*  max_lnblocks,
			   int*  lnblocks,
			   int*  n_to_left,
			   int   n_to_left_level[],
			   int*  mdim,
			   int*  mflags,
			   int   lrefine[], 
			   int   nodetype[], 
			   int   which_child[], 
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
			   int   iorder[]) 
{
  hid_t  file_id, dset_block;  /* file and dataset ids */
  hid_t  block_type_id;

  CREATE_File4(file_name_in, &file_id);

  WRITE_Header4(&file_id, 
		num_components, 
		num_levels,
		compNames,
		ndim,
		mdim,
		nxb,
		nyb,
		nzb,
                lnblocks,
                bnd_box,
		no_at_level);
  
  WRITE_Blocks4(&file_id, 
		num_levels,
		num_components,
		nxb,
		nyb,
		nzb,
		max_lnblocks, 
		n_to_left, 
		n_to_left_level,
		lnblocks, 
		mdim, 
		ndim, 
		mflags, 
		lrefine, 
		nodetype, 
		which_child, 
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
		nx, ny, nz,
		il0, iu0, jl0, ju0, kl0, ku0,
		iorder);

  /* Close resources */
  H5Fclose(file_id);

}



     
