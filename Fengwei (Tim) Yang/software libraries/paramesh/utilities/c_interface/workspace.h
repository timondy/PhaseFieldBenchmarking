/* header file for PARAMESH c_interface which gives acces to data in
   the fortran header module workspace.F */

#include "underscore.h"

#ifdef REAL8
#define float double
#endif

extern void* c_amr_get_pointer(char*);

float (*work)[*kuw][*juw][*iuw][*maxblocksw][*nvar_work];
work = (float (*) [][][][][]) c_amr_get_pointer("work");
(*work)[0][0][0][0][0] = (*work)[0][0][0][0][0];

int (*interp_mask_work)[*nvar_work];
interp_mask_work = (int (*)[]) c_amr_get_pointer("interp_mask_work");
(*interp_mask_work)[0] = (*interp_mask_work)[0];

int (*interp_mask_work_res)[*nvar_work];
interp_mask_work_res = (int (*)[]) c_amr_get_pointer("interp_mask_work_res");
(*interp_mask_work_res)[0] = (*interp_mask_work_res)[0];

float (*work1)[*npblks][*kuw1][*juw1][*iuw1];
work1 = (float (*)[][][][]) c_amr_get_pointer("work1");
(*work1)[0][0][0][0] = (*work1)[0][0][0][0];

if (*curvilinear == 1) {
float (*cell_vol_w) [*kuw1][*juw1][*iuw1];
cell_vol_w = (float (*)[][][]) c_amr_get_pointer("cell_vol_w");
(*cell_vol_w)[0][0][0] = (*cell_vol_w)[0][0][0];
}

