/* header file for PARAMESH c_interface which gives acces to data in
   the fortran header module tree.F */

#include "underscore.h"

#ifdef REAL8
#define float double
#endif

extern void* c_amr_get_pointer(char*);

int two_tree;
two_tree = 2;
int three_tree;
three_tree = 3;

int* maxblocks_tr;
maxblocks_tr = (int *) c_amr_get_pointer("maxblocks_tr");
*maxblocks_tr = *maxblocks_tr;

int* nchild;
nchild = (int *) c_amr_get_pointer("nchild");
*nchild = *nchild;

int* mchild;
mchild = (int *) c_amr_get_pointer("mchild");
*mchild = *mchild;

int* nfaces;
nfaces = (int *) c_amr_get_pointer("nfaces");
*nfaces = *nfaces;

int* mfaces;
mfaces = (int *) c_amr_get_pointer("mfaces");
*mfaces = *mfaces;

int* mdim;
mdim = (int *) c_amr_get_pointer("mdim");
*mdim = *mdim;

int* mflags;
mflags = (int *) c_amr_get_pointer("mflags");
*mflags = *mflags;

int* nboundaries;
nboundaries = (int *) c_amr_get_pointer("nboundaries");
*nboundaries = *nboundaries;

int* mboundaries;
mboundaries = (int *) c_amr_get_pointer("mboundaries");
*mboundaries = *mboundaries;

int (*neigh)[*maxblocks_tr][*mfaces][two_tree];
neigh = (int (*)[][][]) c_amr_get_pointer("neigh");
(*neigh)[0][0][0] = (*neigh)[0][0][0];

int (*child)[*maxblocks_tr][*mchild][two_tree];
child = (int (*)[][][]) c_amr_get_pointer("child");
(*child)[0][0][0] = (*child)[0][0][0];

int (*which_child)[*maxblocks_tr];
which_child = (int (*)[]) c_amr_get_pointer ("which_child");
(*which_child)[0] = (*which_child)[0];

int (*parent)[*maxblocks_tr][two_tree];
parent = (int (*)[][]) c_amr_get_pointer("parent");
(*parent)[0][0] = (*parent)[0][0];

int (*lrefine)[*maxblocks_tr];
lrefine = (int (*)[]) c_amr_get_pointer("lrefine");
(*lrefine)[0] = (*lrefine)[0];

int *lnblocks;
lnblocks = (int *) c_amr_get_pointer("lnblocks");
*lnblocks = *lnblocks;

int *new_lnblocks;
new_lnblocks = (int *) c_amr_get_pointer("new_lnblocks");
*new_lnblocks = *new_lnblocks;

int (*nodetype)[*maxblocks_tr];
nodetype = (int (*)[]) c_amr_get_pointer("nodetype");
(*nodetype)[0] = (*nodetype)[0];

int (*empty)[*maxblocks_tr];
empty = (int (*)[]) c_amr_get_pointer("empty");
(*empty)[0] = (*empty)[0];

int (*bflags)[*maxblocks_tr][*mflags];
bflags = (int (*)[][]) c_amr_get_pointer("bflags");
(*bflags)[0][0] = (*bflags)[0][0];

int (*newchild)[*maxblocks_tr];
newchild = (int (*)[]) c_amr_get_pointer("newchild");
(*newchild)[0] = (*newchild)[0];

int (*derefine)[*maxblocks_tr];
derefine = (int (*)[]) c_amr_get_pointer("derefine");
(*derefine)[0] = (*derefine)[0];

int (*refine)[*maxblocks_tr];
refine = (int (*)[]) c_amr_get_pointer("refine");
(*refine)[0] = (*refine)[0];

int (*stay)[*maxblocks_tr];
stay = (int (*)[]) c_amr_get_pointer("stay");
(*stay)[0] = (*stay)[0];

float (*work_block)[*maxblocks_tr];
work_block = (float (*)[]) c_amr_get_pointer("work_block");
(*work_block)[0] = (*work_block)[0];

float (*coord)[*maxblocks_tr][*mdim];
coord = (float (*)[][]) c_amr_get_pointer("coord");
(*coord)[0][0] = (*coord)[0][0];

float (*bsize)[*maxblocks_tr][*mdim];
bsize = (float (*)[][]) c_amr_get_pointer("bsize");
(*bsize)[0][0] = (*bsize)[0][0];

float (*bnd_box)[*maxblocks_tr][*mdim][two_tree];
bnd_box = (float (*)[][][]) c_amr_get_pointer("bnd_box");
(*bnd_box)[0][0][0] = (*bnd_box)[0][0][0];

int *grid_xmin;
grid_xmin = (int *) c_amr_get_pointer("grid_xmin");
*grid_xmin = *grid_xmin;

int *grid_ymin;
grid_ymin = (int *) c_amr_get_pointer("grid_ymin");
*grid_ymin = *grid_ymin;

int *grid_zmin;
grid_zmin = (int *) c_amr_get_pointer("grid_zmin");
*grid_zmin = *grid_zmin;

int *grid_xmax;
grid_xmax = (int *) c_amr_get_pointer("grid_xmax");
*grid_xmax = *grid_xmax;

int *grid_ymax;
grid_ymax = (int *) c_amr_get_pointer("grid_ymax");
*grid_ymax = *grid_ymax;

int *grid_zmax;
grid_zmax = (int *) c_amr_get_pointer("grid_zmax");
*grid_zmax = *grid_zmax;

int *lrefine_max;
lrefine_max = (int *) c_amr_get_pointer("lrefine_max");
*lrefine_max = *lrefine_max;

int *lrefine_min;
lrefine_min = (int *) c_amr_get_pointer("lrefine_min");
*lrefine_min = *lrefine_min;

float (*level_cell_sizes)[*maxlevels][*mdim];
level_cell_sizes = (float (*)[][]) c_amr_get_pointer("level_cell_sizes");
(*level_cell_sizes)[0][0] = (*level_cell_sizes)[0][0];

int *grid_changed;
grid_changed = (int *) c_amr_get_pointer("grid_changed");
*grid_changed = *grid_changed;

int *grid_analysed_mpi;
grid_analysed_mpi = (int *) c_amr_get_pointer("grid_analysed_mpi");
*grid_analysed_mpi = *grid_analysed_mpi;

float (*boundary_box)[*mboundaries][*mdim][two_tree];
boundary_box = (float (*)[][][]) c_amr_get_pointer("boundary_box");
(*boundary_box)[0][0][0] = (*boundary_box)[0][0][0];

int (*boundary_index)[*mboundaries];
boundary_index = (int (*)[]) c_amr_get_pointer("boundary_index");
(*boundary_index)[0] = (*boundary_index)[0];

int (*surr_blks)[*maxblocks_alloc][3+2**k2d][1+2**k3d][three_tree];
surr_blks = (int (*)[][][][]) c_amr_get_pointer("surr_blks");
(*surr_blks)[0][0][0][0] = (*surr_blks)[0][0][0][0];

