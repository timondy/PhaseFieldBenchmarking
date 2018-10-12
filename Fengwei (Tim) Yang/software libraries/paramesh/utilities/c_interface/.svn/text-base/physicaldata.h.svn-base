/* header file for PARAMESH c_interface which gives acces to data in
   the fortran header module physicaldata.F */

#include "underscore.h"

#ifdef REAL8
#define float double
#endif

extern void* c_amr_get_pointer(char*);

int two_physicaldata;
two_physicaldata = 2;
int three_physicaldata;
three_physicaldata = 3;

float (*unk)[*maxblocks][*ku_bnd][*ju_bnd][*iu_bnd][*nvar];
unk = (float (*)[][][][][]) c_amr_get_pointer("unk");
(*unk)[0][0][0][0][0] = (*unk)[0][0][0][0][0];

int (*interp_mask_unk)[*nvar];
interp_mask_unk = (int (*)[]) c_amr_get_pointer("interp_mask_unk");
(*interp_mask_unk)[0] = (*interp_mask_unk)[0];

int (*interp_mask_unk_res)[*nvar];
interp_mask_unk_res = (int (*)[]) c_amr_get_pointer("interp_mask_unk_res");
(*interp_mask_unk_res)[0] = (*interp_mask_unk_res)[0];

int (*gcell_on_cc)[*nvar];
gcell_on_cc = (int (*)[]) c_amr_get_pointer("gcell_on_cc");
(*gcell_on_cc)[0] = (*gcell_on_cc)[0];

int (*checkp_on_cc)[*nvar];
checkp_on_cc = (int (*)[]) c_amr_get_pointer("checkp_on_cc");
(*checkp_on_cc)[0] = (*checkp_on_cc)[0];

float (*facevarx)[*maxblocksf][*ku_bnd][*ju_bnd][*iu_bnd+1][*nbndvar];
facevarx = (float (*)[][][][][]) c_amr_get_pointer("facevarx");
(*facevarx)[0][0][0][0][0] = (*facevarx)[0][0][0][0][0];

float (*facevary)[*maxblocksf][*ku_bnd][*ju_bnd+*k2d][*iu_bnd][*nbndvar];
facevary = (float (*)[][][][][]) c_amr_get_pointer("facevary");
(*facevary)[0][0][0][0][0] = (*facevary)[0][0][0][0][0];

float (*facevarz)[*maxblocksf][*ku_bnd+*k3d][*ju_bnd][*iu_bnd][*nbndvar];
facevarz = (float (*)[][][][][]) c_amr_get_pointer("facevarz");
(*facevarz)[0][0][0][0][0] = (*facevarz)[0][0][0][0][0];

int (*interp_mask_facex)[*nbndvar];
interp_mask_facex = (int (*)[]) c_amr_get_pointer("interp_mask_facex");
(*interp_mask_facex)[0] = (*interp_mask_facex)[0];

int (*interp_mask_facey)[*nbndvar];
interp_mask_facey = (int (*)[]) c_amr_get_pointer("interp_mask_facey");
(*interp_mask_facey)[0] = (*interp_mask_facey)[0];

int (*interp_mask_facez)[*nbndvar];
interp_mask_facez = (int (*)[]) c_amr_get_pointer("interp_mask_facez");
(*interp_mask_facez)[0] = (*interp_mask_facez)[0];

int (*interp_mask_facex_res)[*nbndvar];
interp_mask_facex_res = (int (*)[]) c_amr_get_pointer("interp_mask_facex_res");
(*interp_mask_facex_res)[0] = (*interp_mask_facex_res)[0];

int (*interp_mask_facey_res)[*nbndvar];
interp_mask_facey_res = (int (*)[]) c_amr_get_pointer("interp_mask_facey_res");
(*interp_mask_facey_res)[0] = (*interp_mask_facey_res)[0];

int (*interp_mask_facez_res)[*nbndvar];
interp_mask_facez_res = (int (*)[]) c_amr_get_pointer("interp_mask_facez_res");
(*interp_mask_facez_res)[0] = (*interp_mask_facez_res)[0];

int (*gcell_on_fc)[*nbndvar][three_physicaldata];
gcell_on_fc = (int (*)[][]) c_amr_get_pointer("gcell_on_fc");
(*gcell_on_fc)[0][0] = (*gcell_on_fc)[0][0];

int (*checkp_on_fc)[*nbndvar][three_physicaldata];
checkp_on_fc = (int (*)[][]) c_amr_get_pointer("checkp_on_fc");
(*checkp_on_fc)[0][0] = (*checkp_on_fc)[0][0];

float (*unk_e_x)[*maxblocksue][*ku_bnd+*k3d][*ju_bnd+*k2d][*iu_bnd][*nbndvare];
unk_e_x = (float (*)[][][][][]) c_amr_get_pointer("unk_e_x");
(*unk_e_x)[0][0][0][0][0] = (*unk_e_x)[0][0][0][0][0];

float (*unk_e_y)[*maxblocksue][*ku_bnd+*k3d][*ju_bnd][*iu_bnd+1][*nbndvare];
unk_e_y = (float (*)[][][][][]) c_amr_get_pointer("unk_e_y");
(*unk_e_y)[0][0][0][0][0] = (*unk_e_y)[0][0][0][0][0];

float (*unk_e_z)[*maxblocksue][*ku_bnd][*ju_bnd+*k2d][*iu_bnd+1][*nbndvare];
unk_e_z = (float (*)[][][][][]) c_amr_get_pointer("unk_e_z");
(*unk_e_z)[0][0][0][0][0] = (*unk_e_z)[0][0][0][0][0];

int (*interp_mask_ec)[*nbndvare];
interp_mask_ec = (int (*)[]) c_amr_get_pointer("interp_mask_ec");
(*interp_mask_ec)[0] = (*interp_mask_ec)[0];

int (*interp_mask_ec_res)[*nbndvare];
interp_mask_ec_res = (int (*)[]) c_amr_get_pointer("interp_mask_ec_res");
(*interp_mask_ec_res)[0] = (*interp_mask_ec_res)[0];

int (*gcell_on_ec)[*nbndvare][three_physicaldata];
gcell_on_ec = (int (*)[][]) c_amr_get_pointer("gcell_on_ec");
(*gcell_on_ec)[0][0] = (*gcell_on_ec)[0][0];

int (*checkp_on_ec)[*nbndvare][three_physicaldata];
checkp_on_ec = (int (*)[][]) c_amr_get_pointer("checkp_on_ec");
(*checkp_on_ec)[0][0] = (*checkp_on_ec)[0][0];

float (*unk_n)[*maxblocksn][*ku_bnd+*k3d][*ju_bnd+*k2d][*iu_bnd+1][*nbndvarc];
unk_n = (float (*)[][][][][]) c_amr_get_pointer("unk_n");
(*unk_n)[0][0][0][0][0] = (*unk_n)[0][0][0][0][0];

int (*interp_mask_nc)[*nbndvarc];
interp_mask_nc = (int (*)[]) c_amr_get_pointer("interp_mask_nc");
(*interp_mask_nc)[0] = (*interp_mask_nc)[0];

int (*interp_mask_nc_res)[*nbndvarc];
interp_mask_nc_res = (int (*)[]) c_amr_get_pointer("interp_mask_nc_res");
(*interp_mask_nc_res)[0] = (*interp_mask_nc_res)[0];

int (*gcell_on_nc)[*nbndvarc];
gcell_on_nc = (int (*)[]) c_amr_get_pointer("gcell_on_nc");
(*gcell_on_nc)[0] = (*gcell_on_nc)[0];

int (*checkp_on_nc)[*nbndvarc];
checkp_on_nc = (int (*)[]) c_amr_get_pointer("checkp_on_nc");
(*checkp_on_nc)[0] = (*checkp_on_nc)[0];

float (*time_loc)[*maxblocks_alloc];
time_loc = (float (*)[]) c_amr_get_pointer("time_loc");
(*time_loc)[0] = (*time_loc)[0];

float (*dtlevel)[*maxlevels];
dtlevel = (float (*)[]) c_amr_get_pointer("dtlevel");
(*dtlevel)[0] = (*dtlevel)[0];

float (*phase_dt)[*maxlevels];
phase_dt = (float (*)[]) c_amr_get_pointer("phase_dt");
(*phase_dt)[0] = (*phase_dt)[0];

int (*loc_cycle)[*maxlevels];
loc_cycle = (int (*)[]) c_amr_get_pointer("loc_cycle");
(*loc_cycle)[0] = (*loc_cycle)[0];

int (*ncyc_local)[*maxlevels];
ncyc_local = (int (*)[]) c_amr_get_pointer("ncyc_local");
(*ncyc_local)[0] = (*ncyc_local)[0];


int (*ldtcomplete)[*maxblocks_alloc];
ldtcomplete = (int (*)[]) c_amr_get_pointer("ldtcomplete");
(*ldtcomplete)[0] = (*ldtcomplete)[0];

int* var_dt;
var_dt = (int *) c_amr_get_pointer("var_dt");
*var_dt = *var_dt;

int* pred_corr;
pred_corr = (int *) c_amr_get_pointer("pred_corr");
*pred_corr = *pred_corr;

if (*var_dt == 1 || *pred_corr == 1) {

float (*t_unk)[*maxblocks][*ku_bnd][*ju_bnd][*iu_bnd][*nvar];
t_unk = (float (*) [][][][][]) c_amr_get_pointer("t_unk");
(*t_unk)[0][0][0][0][0] = (*t_unk)[0][0][0][0][0];

float (*tfacevarx)[*maxblocksf][*ku_bnd][*ju_bnd][*iu_bnd+1][*nbndvar];
tfacevarx = (float (*)[][][][][]) c_amr_get_pointer("tfacevarx");
(*tfacevarx)[0][0][0][0][0] = (*tfacevarx)[0][0][0][0][0];

float (*tfacevary)[*maxblocksf][*ku_bnd][*ju_bnd+*k2d][*iu_bnd][*nbndvar];
tfacevary = (float (*)[][][][][]) c_amr_get_pointer("tfacevary");
(*tfacevary)[0][0][0][0][0] = (*tfacevary)[0][0][0][0][0];

float (*tfacevarz)[*maxblocksf][*ku_bnd+*k3d][*ju_bnd][*iu_bnd][*nbndvar];
tfacevarz = (float (*)[][][][][]) c_amr_get_pointer("tfacevarz");
(*tfacevarz)[0][0][0][0][0] = (*tfacevarz)[0][0][0][0][0];

float (*t_unk_e_x)[*maxblocksue][*ku_bnd+*k3d][*ju_bnd+*k2d][*iu_bnd][*nbndvare];
t_unk_e_x = (float (*)[][][][][]) c_amr_get_pointer("t_unk_e_x");
(*t_unk_e_x)[0][0][0][0][0] = (*t_unk_e_x)[0][0][0][0][0];

float (*t_unk_e_y)[*maxblocksue][*ku_bnd+*k3d][*ju_bnd][*iu_bnd+1][*nbndvare];
t_unk_e_y = (float (*)[][][][][]) c_amr_get_pointer("t_unk_e_y");
(*t_unk_e_y)[0][0][0][0][0] = (*t_unk_e_y)[0][0][0][0][0];

float (*t_unk_e_z)[*maxblocksue][*ku_bnd][*ju_bnd+*k2d][*iu_bnd+1][*nbndvare];
t_unk_e_z = (float (*)[][][][][]) c_amr_get_pointer("t_unk_e_z");
(*t_unk_e_z)[0][0][0][0][0] = (*t_unk_e_z)[0][0][0][0][0];

float (*t_unk_n)[*maxblocksn][*ku_bnd+*k3d][*ju_bnd+*k2d][*iu_bnd+1][*nbndvarc];
t_unk_n = (float (*) [][][][][]) c_amr_get_pointer("t_unk_n");
(*t_unk_n)[0][0][0][0][0] = (*t_unk_n)[0][0][0][0][0];

}

float (*unk1)[*npblks][*ku_bnd1][*ju_bnd1][*iu_bnd1][*nvar];
unk1 = (float (*) [][][][][]) c_amr_get_pointer("unk1");
(*unk1)[0][0][0][0][0] = (*unk1)[0][0][0][0][0];

float (*facevarx1)[*npblks][*ku_bnd1][*ju_bnd1][*iu_bnd1+1][*nbndvar];
facevarx1 = (float (*)[][][][][]) c_amr_get_pointer("facevarx1");
(*facevarx1)[0][0][0][0][0] = (*facevarx1)[0][0][0][0][0];

float (*facevary1)[*npblks][*ku_bnd1][*ju_bnd1+*k2d][*iu_bnd1][*nbndvar];
facevary1 = (float (*)[][][][][]) c_amr_get_pointer("facevary1");
(*facevary1)[0][0][0][0][0] = (*facevary1)[0][0][0][0][0];

float (*facevarz1)[*npblks][*ku_bnd1+*k3d][*ju_bnd1][*iu_bnd1][*nbndvar];
facevarz1 = (float (*)[][][][][]) c_amr_get_pointer("facevarz1");
(*facevarz1)[0][0][0][0][0] = (*facevarz1)[0][0][0][0][0];

float (*unk_e_x1)[*npblks][*ku_bnd1+*k3d][*ju_bnd1+*k2d][*iu_bnd1][*nbndvare];
unk_e_x1 = (float (*)[][][][][]) c_amr_get_pointer("unk_e_x1");
(*unk_e_x1)[0][0][0][0][0] = (*unk_e_x1)[0][0][0][0][0];

float (*unk_e_y1)[*npblks][*ku_bnd1+*k3d][*ju_bnd1][*iu_bnd1+1][*nbndvare];
unk_e_y1 = (float (*)[][][][][]) c_amr_get_pointer("unk_e_y1");
(*unk_e_y1)[0][0][0][0][0] = (*unk_e_y1)[0][0][0][0][0];

float (*unk_e_z1)[*npblks][*ku_bnd1][*ju_bnd1+*k2d][*iu_bnd1+1][*nbndvare];
unk_e_z1 = (float (*)[][][][][]) c_amr_get_pointer("unk_e_z1");
(*unk_e_z1)[0][0][0][0][0] = (*unk_e_z1)[0][0][0][0][0];

float (*unk_n1)[*npblks][*ku_bnd1+*k3d][*ju_bnd1+*k2d][*iu_bnd1+1][*nbndvarc];
unk_n1 = (float (*) [][][][][]) c_amr_get_pointer("unk_n1");
(*unk_n1)[0][0][0][0][0] = (*unk_n1)[0][0][0][0][0];

int* nfluxvar;
nfluxvar = (int *) c_amr_get_pointer("nfluxvar");
*nfluxvar = *nfluxvar;

int* nfluxes;
nfluxes = (int *) c_amr_get_pointer("nfluxes");
*nfluxes = *nfluxes;

float (*flux_x)[*maxblocksfl][*ku_bndi-*kl_bndi+1][*ju_bndi-*jl_bndi+1][two_physicaldata][*nfluxes];
flux_x = (float (*)[][][][][]) c_amr_get_pointer("flux_x");
(*flux_x)[0][0][0][0][0] = (*flux_x)[0][0][0][0][0];

float (*flux_y)[*maxblocksfl][*ku_bndi-*kl_bndi+1][two_physicaldata][*iu_bndi-*il_bndi+1][*nfluxes];
flux_y = (float (*)[][][][][]) c_amr_get_pointer("flux_y");
(*flux_y)[0][0][0][0][0] = (*flux_y)[0][0][0][0][0];

float (*flux_z)[*maxblocksfl][two_physicaldata][*ju_bndi-*jl_bndi+1][*iu_bndi-*il_bndi+1][*nfluxes];
flux_z = (float (*)[][][][][]) c_amr_get_pointer("flux_z");
(*flux_z)[0][0][0][0][0] = (*flux_z)[0][0][0][0][0];

float (*tflux_x)[*maxblocksfl][*ku_bndi-*kl_bndi+1][*ju_bndi-*jl_bndi+1][two_physicaldata][*nfluxes];
tflux_x = (float (*)[][][][][]) c_amr_get_pointer("tflux_x");
(*tflux_x)[0][0][0][0][0] = (*tflux_x)[0][0][0][0][0];

float (*tflux_y)[*maxblocksfl][*ku_bndi-*kl_bndi+1][two_physicaldata][*iu_bndi-*il_bndi+1][*nfluxes];
tflux_y = (float (*)[][][][][]) c_amr_get_pointer("tflux_y");
(*tflux_y)[0][0][0][0][0] = (*tflux_y)[0][0][0][0][0];

float (*tflux_z)[*maxblocksfl][two_physicaldata][*ju_bndi-*jl_bndi+1][*iu_bndi-*il_bndi+1][*nfluxes];
tflux_z = (float (*)[][][][][]) c_amr_get_pointer("tflux_z");
(*tflux_z)[0][0][0][0][0] = (*tflux_z)[0][0][0][0][0];

int* nedgevar1;
nedgevar1 = (int *) c_amr_get_pointer("nedgevar1");
*nedgevar1 = *nedgevar1;

int* nedgevar;
nedgevar = (int *) c_amr_get_pointer("nedgevar");
*nedgevar = *nedgevar;

int* nedges;
nedges = (int *) c_amr_get_pointer("nedges");
*nedges = *nedges;

float (*bedge_facex_y)[*maxblockse][*ku_bnd+1][*ju_bnd+1][two_physicaldata][*nedges];
bedge_facex_y = (float (*)[][][][][]) c_amr_get_pointer("bedge_facex_y");
(*bedge_facex_y)[0][0][0][0][0] = (*bedge_facex_y)[0][0][0][0][0];

float (*bedge_facex_z)[*maxblockse][*ku_bnd+1][*ju_bnd+1][two_physicaldata][*nedges];
bedge_facex_z = (float (*)[][][][][]) c_amr_get_pointer("bedge_facex_z");
(*bedge_facex_z)[0][0][0][0][0] = (*bedge_facex_z)[0][0][0][0][0];

float (*bedge_facey_x)[*maxblockse][*ku_bnd+1][two_physicaldata][*iu_bnd+1][*nedges];
bedge_facey_x = (float (*)[][][][][]) c_amr_get_pointer("bedge_facey_x");
(*bedge_facey_x)[0][0][0][0][0] = (*bedge_facey_x)[0][0][0][0][0];

float (*bedge_facey_z)[*maxblockse][*ku_bnd+1][two_physicaldata][*iu_bnd+1][*nedges];
bedge_facey_z = (float (*)[][][][][]) c_amr_get_pointer("bedge_facey_z");
(*bedge_facey_z)[0][0][0][0][0] = (*bedge_facey_z)[0][0][0][0][0];

float (*bedge_facez_x)[*maxblockse][two_physicaldata][*ju_bnd+1][*iu_bnd+1][*nedges];
bedge_facez_x = (float (*)[][][][][]) c_amr_get_pointer("bedge_facez_x");
(*bedge_facez_x)[0][0][0][0][0] = (*bedge_facez_x)[0][0][0][0][0];

float (*bedge_facez_y)[*maxblockse][two_physicaldata][*ju_bnd+1][*iu_bnd+1][*nedges];
bedge_facez_y = (float (*)[][][][][]) c_amr_get_pointer("bedge_facez_y");
(*bedge_facez_y)[0][0][0][0][0] = (*bedge_facez_y)[0][0][0][0][0];

float (*tbedge_facex_y)[*maxblockse][*ku_bnd+1][*ju_bnd+1][two_physicaldata][*nedges];
tbedge_facex_y = (float (*)[][][][][]) c_amr_get_pointer("tbedge_facex_y");
(*tbedge_facex_y)[0][0][0][0][0] = (*tbedge_facex_y)[0][0][0][0][0];

float (*tbedge_facex_z)[*maxblockse][*ku_bnd+1][*ju_bnd+1][two_physicaldata][*nedges];
tbedge_facex_z = (float (*)[][][][][]) c_amr_get_pointer("tbedge_facex_z");
(*tbedge_facex_z)[0][0][0][0][0] = (*tbedge_facex_z)[0][0][0][0][0];

float (*tbedge_facey_x)[*maxblockse][*ku_bnd+1][two_physicaldata][*iu_bnd+1][*nedges];
tbedge_facey_x = (float (*)[][][][][]) c_amr_get_pointer("tbedge_facey_x");
(*tbedge_facey_x)[0][0][0][0][0] = (*tbedge_facey_x)[0][0][0][0][0];

float (*tbedge_facey_z)[*maxblockse][*ku_bnd+1][two_physicaldata][*iu_bnd+1][*nedges];
tbedge_facey_z = (float (*)[][][][][]) c_amr_get_pointer("tbedge_facey_z");
(*tbedge_facey_z)[0][0][0][0][0] = (*tbedge_facey_z)[0][0][0][0][0];

float (*tbedge_facez_x)[*maxblockse][two_physicaldata][*ju_bnd+1][*iu_bnd+1][*nedges];
tbedge_facez_x = (float (*)[][][][][]) c_amr_get_pointer("tbedge_facez_x");
(*tbedge_facez_x)[0][0][0][0][0] = (*tbedge_facez_x)[0][0][0][0][0];

float (*tbedge_facez_y)[*maxblockse][two_physicaldata][*ju_bnd+1][*iu_bnd+1][*nedges];
tbedge_facez_y = (float (*)[][][][][]) c_amr_get_pointer("tbedge_facez_y");
(*tbedge_facez_y)[0][0][0][0][0] = (*tbedge_facez_y)[0][0][0][0][0];

int* curvilinear;
curvilinear = (int *) c_amr_get_pointer("curvilinear");
*curvilinear = *curvilinear;

if (*curvilinear == 1) {

float (*cell_vol)[*ku_bnd1][*ju_bnd1][*iu_bnd1];
cell_vol = (float (*)[][][]) c_amr_get_pointer("cell_vol");
(*cell_vol)[0][0][0] = (*cell_vol)[0][0][0];

float (*cell_area1)[*ku_bnd1][*ju_bnd1][*iu_bnd1+1];
cell_area1 = (float (*)[][][]) c_amr_get_pointer("cell_area1");
(*cell_area1)[0][0][0] = (*cell_area1)[0][0][0];

float (*cell_area2)[*ku_bnd1][*ju_bnd1+*k2d][*iu_bnd1];
cell_area2 = (float (*)[][][]) c_amr_get_pointer("cell_area2");
(*cell_area2)[0][0][0] = (*cell_area2)[0][0][0];

float (*cell_area3)[*ku_bnd1+*k3d][*ju_bnd1][*iu_bnd1];
cell_area3 = (float (*)[][][]) c_amr_get_pointer("cell_area3");
(*cell_area3)[0][0][0] = (*cell_area3)[0][0][0];

float (*cell_leng1)[*ku_bnd1+*k3d][*ju_bnd1+*k2d][*iu_bnd1];
cell_leng1 = (float (*)[][][]) c_amr_get_pointer("cell_leng1");
(*cell_leng1)[0][0][0] = (*cell_leng1)[0][0][0];

float (*cell_leng2)[*ku_bnd1+*k3d][*ju_bnd1][*iu_bnd1+1];
cell_leng2 = (float (*)[][][]) c_amr_get_pointer("cell_leng2");
(*cell_leng2)[0][0][0] = (*cell_leng2)[0][0][0];

float (*cell_leng3)[*ku_bnd1][*ju_bnd1+*k2d][*iu_bnd1+1];
cell_leng3 = (float (*)[][][]) c_amr_get_pointer("cell_leng3");
(*cell_leng3)[0][0][0] = (*cell_leng3)[0][0][0];

}

if (*nfield_divf > 0) { 
int (*i_divf_fc_vars)[*nfield_divf][three_physicaldata];
i_divf_fc_vars = (int (*)[][]) c_amr_get_pointer("i_divf_fc_vars");
(*i_divf_fc_vars)[0][0] = (*i_divf_fc_vars)[0][0];
}

/*****************/

int* diagonals;
diagonals = (int *) c_amr_get_pointer("diagonals");
*diagonals = *diagonals;

int* amr_error_checking;
amr_error_checking = (int *) c_amr_get_pointer("amr_error_checking");
*amr_error_checking = *amr_error_checking;

int* no_permanent_guardcells;
no_permanent_guardcells = (int *) c_amr_get_pointer("no_permanent_guardcells");
*no_permanent_guardcells = *no_permanent_guardcells;

int* advance_all_levels;
advance_all_levels = (int *) c_amr_get_pointer("advance_all_levels");
*advance_all_levels = *advance_all_levels;

int* force_consistency;
force_consistency = (int *) c_amr_get_pointer("force_consistency");
*force_consistency = *force_consistency;

int* consv_fluxes;
consv_fluxes = (int *) c_amr_get_pointer("consv_fluxes");
*consv_fluxes = *consv_fluxes;

int* consv_flux_densities;
consv_flux_densities = (int *) c_amr_get_pointer("consv_flux_densities");
*consv_flux_densities = *consv_flux_densities;

int* edge_value;
edge_value = (int *) c_amr_get_pointer("edge_value");
*edge_value = *edge_value;

int* edge_value_integ;
edge_value_integ = (int *) c_amr_get_pointer("edge_value_integ");
*edge_value_integ = *edge_value_integ;

int* empty_cells;
empty_cells = (int *) c_amr_get_pointer("empty_cells");
*empty_cells = *empty_cells;

int* conserve;
conserve = (int *) c_amr_get_pointer("conserve");
*conserve = *conserve;

int* clean_divb;
clean_divb = (int *) c_amr_get_pointer("clean_divb");
*clean_divb = *clean_divb;

int* divergence_free;
divergence_free = (int *) c_amr_get_pointer("divergence_free");
*divergence_free = *divergence_free;

int* cartesian_pm;
cartesian_pm = (int *) c_amr_get_pointer("cartesian_pm");
*cartesian_pm = *cartesian_pm;

int* cylindrical_pm;
cylindrical_pm = (int *) c_amr_get_pointer("cylindrical_pm");
*cylindrical_pm = *cylindrical_pm;

int* spherical_pm;
spherical_pm = (int *) c_amr_get_pointer("spherical_pm");
*spherical_pm = *spherical_pm;

int* polar_pm;
polar_pm = (int *) c_amr_get_pointer("polar_pm");
*polar_pm = *polar_pm;
