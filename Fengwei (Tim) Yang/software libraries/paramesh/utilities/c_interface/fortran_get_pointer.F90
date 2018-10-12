!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine fortran_get_pointer(ptr,name_in,string_len)

      use physicaldata
      use tree
      use paramesh_dimensions
      use workspace
      use io

      implicit none

      integer(kind=selected_int_kind(14)) :: ptr
      integer :: string_len

      ! These variables must be declared w. save (static) so that they
      ! occupy identifiable memory addresses by C
      integer, save :: mmdim, mmchild, mmfaces
      integer, save :: nnfluxes, nnfluxvar
      integer, save :: nnedgevar1, nnedgevar, nnedges
      integer, save :: nnchild, nnfaces, nnboundaries, mmboundaries
      integer, save :: mmflags
      integer, save :: nndim,ll2p5d,nnbedges,kk3d,kk2d,kk1d
      integer, save :: nnxb, nnyb, nnzb
      integer, save :: mmaxdim,ggc_off_x,ggc_off_y,ggc_off_z
      integer, save :: mmaxblocks_tr,nnguard,nnpgs
      integer, save :: nnvar,nnfacevar,nnvaredge,nnvarcorn
      integer, save :: iiface_off
      integer, save :: iil_bnd,jjl_bnd,kkl_bnd
      integer, save :: iiu_bnd,jju_bnd,kku_bnd
      integer, save :: iil_bndi,jjl_bndi,kkl_bndi
      integer, save :: iiu_bndi,jju_bndi,kku_bndi
      integer, save :: nnbndvar,nnbndvare,nnbndvarc,nnpblks
      integer, save :: iil_bnd1,jjl_bnd1,kkl_bnd1
      integer, save :: iiu_bnd1,jju_bnd1,kku_bnd1
      integer, save :: nnvar_work, nnguard_work
      integer, save :: iilw,jjlw,kklw,iiuw,jjuw,kkuw
      integer, save :: iilw1,jjlw1,kklw1,iiuw1,jjuw1,kkuw1
      integer, save :: mmaxlevels,nnfield_divf
      integer, save :: mmaxblocks, mmaxblocks_alloc
      integer, save :: mmaxblocksf, mmaxblocksue, mmaxblocksn
      integer, save :: mmaxblocksfl, mmaxblockse, mmaxblocksw
      logical, save :: aamr_error_checking, nno_permanent_guardcells
      logical, save :: aadvance_all_levels, fforce_consistency
      logical, save :: cconsv_fluxes, cconsv_flux_densities
      logical, save :: eedge_value, eedge_value_integ
      logical, save :: vvar_dt, ppred_corr, eempty_cells
      logical, save :: cconserve, ddivergence_free
      logical, save :: ccurvilinear, ccartesian_pm
      logical, save :: ccylindrical_pm
      logical, save :: sspherical_pm, ppolar_pm
      logical, save :: ddiagonals

      character (len=string_len) name_in,name

      name = trim(name_in)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA from module PHYSICALDATA !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (name == 'unk') then
         call c_get_pointer(unk,ptr)

      elseif (name == 'interp_mask_unk') then
         call c_get_pointer(interp_mask_unk,ptr)
         
      elseif (name == 'interp_mask_unk_res') then
         call c_get_pointer(interp_mask_unk_res,ptr) 
         
      elseif (name == 'gcell_on_cc') then
         call c_get_pointer(gcell_on_cc,ptr)

      elseif (name == 'checkp_on_cc') then
         call c_get_pointer(checkp_on_cc,ptr) 

      elseif (name == 'facevarx') then
         call c_get_pointer(facevarx,ptr) 
      elseif (name == 'facevary') then
         call c_get_pointer(facevary,ptr) 
      elseif (name == 'facevarz') then
         call c_get_pointer(facevarz,ptr)

      elseif (name == 'interp_mask_facex') then
         call c_get_pointer(interp_mask_facex,ptr)
      elseif (name == 'interp_mask_facey') then
         call c_get_pointer(interp_mask_facey,ptr)
      elseif (name == 'interp_mask_facez') then
         call c_get_pointer(interp_mask_facez,ptr)

      elseif (name == 'interp_mask_facex_res') then
         call c_get_pointer(interp_mask_facex_res,ptr) 
      elseif (name == 'interp_mask_facey_res') then
         call c_get_pointer(interp_mask_facey_res,ptr) 
      elseif (name == 'interp_mask_facez_res') then
         call c_get_pointer(interp_mask_facez_res,ptr) 

      elseif (name == 'gcell_on_fc') then
         call c_get_pointer(gcell_on_fc,ptr) 

      elseif (name == 'checkp_on_fc') then
         call c_get_pointer(checkp_on_fc,ptr) 

      elseif (name == 'unk_e_x') then
         call c_get_pointer(unk_e_x,ptr)
      elseif (name == 'unk_e_y') then
         call c_get_pointer(unk_e_y,ptr)
      elseif (name == 'unk_e_z') then
         call c_get_pointer(unk_e_z,ptr)

      elseif (name == 'interp_mask_ec') then
         call c_get_pointer(interp_mask_ec,ptr)

      elseif (name == 'interp_mask_ec_res') then
         call c_get_pointer(interp_mask_ec_res,ptr)

      elseif (name == 'gcell_on_ec') then
         call c_get_pointer(gcell_on_ec,ptr)

      elseif (name == 'checkp_on_ec') then
         call c_get_pointer(checkp_on_ec,ptr)

      elseif (name == 'unk_n') then
         call c_get_pointer(unk_n,ptr)

      elseif (name == 'interp_mask_nc') then
         call c_get_pointer(interp_mask_nc,ptr) 

      elseif (name == 'interp_mask_nc_res') then
         call c_get_pointer(interp_mask_nc_res,ptr)

      elseif (name == 'gcell_on_nc') then
         call c_get_pointer(gcell_on_nc,ptr)

      elseif (name == 'checkp_on_nc') then
         call c_get_pointer(checkp_on_nc,ptr)

      elseif (name == 'time_loc') then
         call c_get_pointer(time_loc,ptr)

      elseif (name == 'dtlevel') then
         call c_get_pointer(dtlevel,ptr)

      elseif (name == 'phase_dt') then
         call c_get_pointer(phase_dt,ptr)

      elseif (name == 'loc_cycle') then
         call c_get_pointer(loc_cycle,ptr)

      elseif (name == 'ncyc_local') then
         call c_get_pointer(ncyc_local,ptr)

      elseif (name == 'ldtcomplete') then
         call c_get_pointer(ldtcomplete,ptr)

      elseif (name == 't_unk') then
         call c_get_pointer(t_unk,ptr)

      elseif (name == 'tfacevarx') then
         call c_get_pointer(tfacevarx,ptr)
      elseif (name == 'tfacevary') then
         call c_get_pointer(tfacevary,ptr)
      elseif (name == 'tfacevarz') then
         call c_get_pointer(tfacevarz,ptr)

      elseif (name == 't_unk_e_x') then
         call c_get_pointer(t_unk_e_x,ptr)
      elseif (name == 't_unk_e_y') then
         call c_get_pointer(t_unk_e_y,ptr)
      elseif (name == 't_unk_e_z') then
         call c_get_pointer(t_unk_e_z,ptr)

      elseif (name == 't_unk_n') then
         call c_get_pointer(t_unk_n,ptr)

      elseif (name == 'unk1') then
         call c_get_pointer(unk1,ptr)

      elseif (name == 'facevarx1') then
         call c_get_pointer(facevarx1,ptr)
      elseif (name == 'facevary1') then
         call c_get_pointer(facevary1,ptr)
      elseif (name == 'facevarz1') then
         call c_get_pointer(facevarz1,ptr)

      elseif (name == 'unk_e_x1') then
         call c_get_pointer(unk_e_x1,ptr)
      elseif (name == 'unk_e_y1') then
         call c_get_pointer(unk_e_y1,ptr)
      elseif (name == 'unk_e_z1') then
         call c_get_pointer(unk_e_z1,ptr)

      elseif (name == 'unk_n1') then
         call c_get_pointer(unk_n1,ptr)

      elseif (name == 'nfluxvar') then
         nnfluxvar = nfluxvar
         call c_get_pointer(nnfluxvar,ptr)
         
      elseif (name == 'nfluxes') then
         nnfluxes = nfluxes
         call c_get_pointer(nnfluxes,ptr)

      elseif (name == 'flux_x') then
         call c_get_pointer(flux_x,ptr)
      elseif (name == 'flux_y') then
         call c_get_pointer(flux_y,ptr)
      elseif (name == 'flux_z') then
         call c_get_pointer(flux_z,ptr)

      elseif (name == 'tflux_x') then
         call c_get_pointer(tflux_x,ptr)
      elseif (name == 'tflux_y') then
         call c_get_pointer(tflux_y,ptr)
      elseif (name == 'tflux_z') then
         call c_get_pointer(tflux_z,ptr)

      elseif (name == 'nedgevar1') then
         nnedgevar1 = nedgevar1
         call c_get_pointer(nnedgevar1,ptr)

      elseif (name == 'nedgevar') then
         nnedgevar = nedgevar
         call c_get_pointer(nnedgevar,ptr)

      elseif (name == 'nedges') then
         nnedges = nedges
         call c_get_pointer(nnedges,ptr)

      elseif (name == 'bedge_facex_y') then
         call c_get_pointer(bedge_facex_y,ptr)
      elseif (name == 'bedge_facex_z') then
         call c_get_pointer(bedge_facex_z,ptr)
      elseif (name == 'bedge_facey_x') then
         call c_get_pointer(bedge_facey_x,ptr)
      elseif (name == 'bedge_facey_z') then
         call c_get_pointer(bedge_facey_z,ptr)
      elseif (name == 'bedge_facez_x') then
         call c_get_pointer(bedge_facez_x,ptr)
      elseif (name == 'bedge_facez_y') then
         call c_get_pointer(bedge_facez_y,ptr)

      elseif (name == 'tbedge_facex_y') then
         call c_get_pointer(tbedge_facex_y,ptr)
      elseif (name == 'tbedge_facex_z') then
         call c_get_pointer(tbedge_facex_z,ptr)
      elseif (name == 'tbedge_facey_x') then
         call c_get_pointer(tbedge_facey_x,ptr)
      elseif (name == 'tbedge_facey_z') then
         call c_get_pointer(tbedge_facey_z,ptr)
      elseif (name == 'tbedge_facez_x') then
         call c_get_pointer(tbedge_facez_x,ptr)
      elseif (name == 'tbedge_facez_y') then
         call c_get_pointer(tbedge_facez_y,ptr)

      elseif (name == 'cell_vol') then
         call c_get_pointer(cell_vol,ptr)
         
      elseif (name == 'cell_area1') then
         call c_get_pointer(cell_area1,ptr)
      elseif (name == 'cell_area2') then
         call c_get_pointer(cell_area2,ptr)
      elseif (name == 'cell_area3') then
         call c_get_pointer(cell_area3,ptr)

      elseif (name == 'cell_leng1') then
         call c_get_pointer(cell_leng1,ptr)
      elseif (name == 'cell_leng2') then
         call c_get_pointer(cell_leng2,ptr)
      elseif (name == 'cell_leng3') then
         call c_get_pointer(cell_leng3,ptr)
         
      elseif (name == 'i_divf_fc_vars') then
         call c_get_pointer(i_divf_fc_vars,ptr)

      elseif (name == 'diagonals') then
         ddiagonals = diagonals
         call c_get_pointer(ddiagonals,ptr)

      elseif (name == 'amr_error_checking') then
         aamr_error_checking = amr_error_checking
         call c_get_pointer(aamr_error_checking,ptr);

      elseif (name == 'no_permanent_guardcells') then
         nno_permanent_guardcells = no_permanent_guardcells
         call c_get_pointer(nno_permanent_guardcells,ptr)

      elseif (name == 'advance_all_levels') then
         aadvance_all_levels = advance_all_levels
         call c_get_pointer(aadvance_all_levels,ptr)

      elseif (name == 'force_consistency') then
         fforce_consistency = force_consistency
         call c_get_pointer(fforce_consistency,ptr)

      elseif (name == 'consv_fluxes') then
         cconsv_fluxes = consv_fluxes
         call c_get_pointer(cconsv_fluxes,ptr)

      elseif (name == 'consv_flux_densities') then
         cconsv_flux_densities = consv_flux_densities
         call c_get_pointer(cconsv_flux_densities,ptr)

      elseif (name == 'edge_value') then
         eedge_value = edge_value
         call c_get_pointer(eedge_value,ptr)

      elseif (name == 'edge_value_integ') then
         eedge_value_integ = edge_value_integ
         call c_get_pointer(eedge_value_integ,ptr)

      elseif (name == 'var_dt') then
         vvar_dt = var_dt
         call c_get_pointer(vvar_dt,ptr)

      elseif (name == 'pred_corr') then
         ppred_corr = pred_corr
         call c_get_pointer(ppred_corr,ptr)

      elseif (name == 'empty_cells') then
         eempty_cells = empty_cells
         call c_get_pointer(eempty_cells,ptr)

      elseif (name == 'conserve') then
         cconserve = conserve
         call c_get_pointer(cconserve,ptr)

      elseif (name == 'divergence_free') then
         ddivergence_free = divergence_free
         call c_get_pointer(ddivergence_free,ptr)

      elseif (name == 'curvilinear') then
         ccurvilinear = curvilinear
         call c_get_pointer(ccurvilinear,ptr)

      elseif (name == 'cartesian_pm') then
         ccartesian_pm = cartesian_pm
         call c_get_pointer(ccartesian_pm,ptr)

      elseif (name == 'cylindrical_pm') then
         ccylindrical_pm = cylindrical_pm
         call c_get_pointer(ccylindrical_pm,ptr)

      elseif (name == 'spherical_pm') then
         sspherical_pm = spherical_pm
         call c_get_pointer(sspherical_pm,ptr)

      elseif (name == 'polar_pm') then
         ppolar_pm = polar_pm
         call c_get_pointer(ppolar_pm,ptr)
         

! There are many more variables in physicaldata, but I'm pretty sure
! they will not and SHOULD not be accessed by users except perhaps by
! reading those parameters in in amr_runtime_parameters
         
!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA from module TREE !
!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif (name == 'maxblocks_tr') then
         mmaxblocks_tr = maxblocks_tr
         call c_get_pointer(mmaxblocks_tr,ptr)

      elseif (name == 'nchild') then
         nnchild = nchild
         call c_get_pointer(nnchild,ptr)

      elseif (name == 'mchild') then
         mmchild = mchild
         call c_get_pointer(mmchild,ptr)

      elseif (name == 'nfaces') then
         nnfaces = nfaces
         call c_get_pointer(nnfaces,ptr)

      elseif (name == 'mfaces') then
         mmfaces = mfaces
         call c_get_pointer(mmfaces,ptr)

      elseif (name == 'mdim') then
         mmdim = mdim
         call c_get_pointer(mmdim,ptr)

      elseif (name == 'mflags') then
         mmflags = mflags
         call c_get_pointer(mmflags,ptr)

      elseif (name == 'nboundaries') then
         nnboundaries = nboundaries
         call c_get_pointer(nnboundaries,ptr)

      elseif (name == 'mboundaries') then
         mmboundaries = mboundaries
         call c_get_pointer(mmboundaries,ptr)

      elseif (name == 'neigh') then
         call c_get_pointer(neigh,ptr)

      elseif (name == 'child') then
         call c_get_pointer(child,ptr)

      elseif (name == 'which_child') then
         call c_get_pointer(which_child,ptr)

      elseif (name == 'parent') then
         call c_get_pointer(parent,ptr)

      elseif (name == 'lrefine') then
         call c_get_pointer(lrefine,ptr)

      elseif (name == 'lnblocks') then
         call c_get_pointer(lnblocks,ptr)

      elseif (name == 'new_lnblocks') then
         call c_get_pointer(new_lnblocks,ptr)

      elseif (name == 'nodetype') then
         call c_get_pointer(nodetype,ptr)

      elseif (name == 'empty') then
         call c_get_pointer(empty,ptr)

      elseif (name == 'bflags') then
         call c_get_pointer(bflags,ptr)

      elseif (name == 'newchild') then
         call c_get_pointer(newchild,ptr)

      elseif (name == 'derefine') then
         call c_get_pointer(derefine,ptr)

      elseif (name == 'refine') then
         call c_get_pointer(refine,ptr)

      elseif (name == 'stay') then
         call c_get_pointer(stay,ptr)

      elseif (name == 'work_block') then
         call c_get_pointer(work_block,ptr)

      elseif (name == 'coord') then
         call c_get_pointer(coord,ptr)

      elseif (name == 'bsize') then
         call c_get_pointer(bsize,ptr)

      elseif (name == 'bnd_box') then
         call c_get_pointer(bnd_box,ptr)

      elseif (name == 'grid_xmin') then
         call c_get_pointer(grid_xmin,ptr)
      elseif (name == 'grid_ymin') then
         call c_get_pointer(grid_ymin,ptr)
      elseif (name == 'grid_zmin') then
         call c_get_pointer(grid_zmin,ptr)
      elseif (name == 'grid_xmax') then
         call c_get_pointer(grid_xmax,ptr)
      elseif (name == 'grid_ymax') then
         call c_get_pointer(grid_ymax,ptr)
      elseif (name == 'grid_zmax') then
         call c_get_pointer(grid_zmax,ptr)

      elseif (name == 'lrefine_max') then
         call c_get_pointer(lrefine_max,ptr)

      elseif (name == 'lrefine_min') then
         call c_get_pointer(lrefine_min,ptr)

      elseif (name == 'level_cell_sizes') then
         call c_get_pointer(level_cell_sizes,ptr)

      elseif (name == 'grid_changed') then
         call c_get_pointer(grid_changed,ptr)

      elseif (name == 'grid_analysed_mpi') then
         call c_get_pointer(grid_analysed_mpi,ptr)

      elseif (name == 'boundary_box') then
         call c_get_pointer(boundary_box,ptr)

      elseif (name == 'boundary_index') then
         call c_get_pointer(boundary_index,ptr)

      elseif (name == 'surr_blks') then
         call c_get_pointer(surr_blks,ptr)

#ifdef SAVE_MORTS
      elseif (name == 'surr_morts') then
         call c_get_pointer(surr_morts,ptr)
#endif

      elseif (name == 'bc_block_neighs') then
         call c_get_pointer(bc_block_neighs,ptr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA from module PARAMESH_DIMENSIONS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif (name == 'maxblocks') then
         mmaxblocks = maxblocks
         call c_get_pointer(mmaxblocks,ptr)

      elseif (name == 'maxblocksf') then
         mmaxblocksf = maxblocksf
         call c_get_pointer(mmaxblocksf,ptr)

      elseif (name == 'maxblocksue') then
         mmaxblocksue = maxblocksue
         call c_get_pointer(mmaxblocksue,ptr)

      elseif (name == 'maxblocksn') then
         mmaxblocksn = maxblocksn
         call c_get_pointer(mmaxblocksn,ptr)

      elseif (name == 'maxblocksfl') then
         mmaxblocksfl = maxblocksfl
         call c_get_pointer(mmaxblocksfl,ptr)

      elseif (name == 'maxblockse') then
         mmaxblockse = maxblockse
         call c_get_pointer(mmaxblockse,ptr)

      elseif (name == 'maxblocksw') then
         mmaxblocksw = maxblocksw
         call c_get_pointer(mmaxblocksw,ptr)

      elseif (name == 'maxblocks_alloc') then
         mmaxblocks_alloc = maxblocks_alloc
         call c_get_pointer(mmaxblocks_alloc,ptr)

      elseif (name == 'ndim') then
         nndim = ndim
         call c_get_pointer(nndim,ptr)

      elseif (name == 'l2p5d') then
         ll2p5d = l2p5d
         call c_get_pointer(ll2p5d,ptr)

      elseif (name == 'nbedges') then
         nnbedges = nbedges
         call c_get_pointer(nnbedges,ptr)

      elseif (name == 'k3d') then
         kk3d = k3d
         call c_get_pointer(kk3d,ptr)
      elseif (name == 'k2d') then
         kk2d = k2d
         call c_get_pointer(kk2d,ptr)
      elseif (name == 'k1d') then
         kk1d = k1d
         call c_get_pointer(kk1d,ptr)

      elseif (name == 'nxb') then
         nnxb = nxb
         call c_get_pointer(nnxb,ptr)
      elseif (name == 'nyb') then
         nnyb = nyb
         call c_get_pointer(nnyb,ptr) 
      elseif (name == 'nzb') then
         nnzb = nzb
         call c_get_pointer(nnzb,ptr)

      elseif (name == 'maxdim') then
         mmaxdim = maxdim
         call c_get_pointer(mmaxdim,ptr)

      elseif (name == 'gc_off_x') then 
         ggc_off_x = gc_off_x
         call c_get_pointer(ggc_off_x,ptr)
      elseif (name == 'gc_off_y') then
         ggc_off_y = gc_off_y
         call c_get_pointer(ggc_off_y,ptr)
      elseif (name == 'gc_off_z') then 
         ggc_off_z = gc_off_z
         call c_get_pointer(ggc_off_z,ptr)

      elseif (name == 'nguard') then
         nnguard = nguard
         call c_get_pointer(nnguard,ptr)

      elseif (name == 'npgs') then
         nnpgs = npgs
         call c_get_pointer(nnpgs,ptr)

      elseif (name == 'nvar') then
         nnvar = nvar
         call c_get_pointer(nnvar,ptr)

      elseif (name == 'nfacevar') then
         nnfacevar = nfacevar
         call c_get_pointer(nnfacevar,ptr)

      elseif (name == 'nvaredge') then
         nnvaredge = nvaredge
         call c_get_pointer(nnvaredge,ptr)

      elseif (name == 'nvarcorn') then
         nnvarcorn = nvarcorn
         call c_get_pointer(nnvarcorn,ptr)

      elseif (name == 'iface_off') then
         iiface_off = iface_off
         call c_get_pointer(iiface_off,ptr)

      elseif (name == 'il_bnd') then
         iil_bnd = il_bnd
         call c_get_pointer(iil_bnd,ptr)
      elseif (name == 'jl_bnd') then
         jjl_bnd = jl_bnd
         call c_get_pointer(jjl_bnd,ptr)
      elseif (name == 'kl_bnd') then
         kkl_bnd = kl_bnd
         call c_get_pointer(kkl_bnd,ptr)
      elseif (name == 'iu_bnd') then
         iiu_bnd = iu_bnd
         call c_get_pointer(iiu_bnd,ptr)
      elseif (name == 'ju_bnd') then
         jju_bnd = ju_bnd
         call c_get_pointer(jju_bnd,ptr)
      elseif (name == 'ku_bnd') then
         kku_bnd = ku_bnd
         call c_get_pointer(kku_bnd,ptr)

      elseif (name == 'il_bndi') then
         iil_bndi = il_bndi
         call c_get_pointer(iil_bndi,ptr)
      elseif (name == 'jl_bndi') then
         jjl_bndi = jl_bndi
         call c_get_pointer(jjl_bndi,ptr)
      elseif (name == 'kl_bndi') then
         kkl_bndi = kl_bndi
         call c_get_pointer(kkl_bndi,ptr)
      elseif (name == 'iu_bndi') then
         iiu_bndi = iu_bndi
         call c_get_pointer(iiu_bndi,ptr)
      elseif (name == 'ju_bndi') then
         jju_bndi = ju_bndi
         call c_get_pointer(jju_bndi,ptr)
      elseif (name == 'ku_bndi') then
         kku_bndi = ku_bndi
         call c_get_pointer(kku_bndi,ptr)

      elseif (name == 'nbndvar') then
         nnbndvar = nbndvar
         call c_get_pointer(nnbndvar,ptr)
         
      elseif (name == 'nbndvare') then
         nnbndvare = nbndvare
         call c_get_pointer(nnbndvare,ptr)

      elseif (name == 'nbndvarc') then
         nnbndvarc = nbndvarc
         call c_get_pointer(nnbndvarc,ptr)

      elseif (name == 'npblks') then
         nnpblks = npblks
         call c_get_pointer(nnpblks,ptr)

      elseif (name == 'il_bnd1') then
         iil_bnd1 = il_bnd1
         call c_get_pointer(iil_bnd1,ptr)
      elseif (name == 'jl_bnd1') then
         jjl_bnd1 = jl_bnd1
         call c_get_pointer(jjl_bnd1,ptr)
      elseif (name == 'kl_bnd1') then
         kkl_bnd1 = kl_bnd1
         call c_get_pointer(kkl_bnd1,ptr)
      elseif (name == 'iu_bnd1') then
         iiu_bnd1 = iu_bnd1
         call c_get_pointer(iiu_bnd1,ptr)
      elseif (name == 'ju_bnd1') then
         jju_bnd1 = ju_bnd1
         call c_get_pointer(jju_bnd1,ptr)
      elseif (name == 'ku_bnd1') then
         kku_bnd1 = ku_bnd1
         call c_get_pointer(kku_bnd1,ptr)

      elseif (name == 'nguard_work') then
         nnguard_work = nguard_work
         call c_get_pointer(nnguard_work,ptr)

      elseif (name == 'nvar_work') then
         nnvar_work = nvar_work
         call c_get_pointer(nnvar_work,ptr)

      elseif (name == 'ilw') then
         iilw = ilw
         call c_get_pointer(iilw,ptr)
      elseif (name == 'jlw') then
         jjlw = jlw
         call c_get_pointer(jjlw,ptr)
      elseif (name == 'klw') then
         kklw = klw
         call c_get_pointer(kklw,ptr)
      elseif (name == 'iuw') then
         iiuw = iuw
         call c_get_pointer(iiuw,ptr)
      elseif (name == 'juw') then
         jjuw = juw
         call c_get_pointer(jjuw,ptr)
      elseif (name == 'kuw') then
         kkuw = kuw
         call c_get_pointer(kkuw,ptr)

      elseif (name == 'ilw1') then
         iilw1 = ilw1
         call c_get_pointer(iilw1,ptr)
      elseif (name == 'jlw1') then
         jjlw1 = jlw1
         call c_get_pointer(jjlw1,ptr)
      elseif (name == 'klw1') then
         kklw1 = klw1
         call c_get_pointer(kklw1,ptr)
      elseif (name == 'iuw1') then
         iiuw1 = iuw1
         call c_get_pointer(iiuw1,ptr)
      elseif (name == 'juw1') then
         jjuw1 = juw1
         call c_get_pointer(jjuw1,ptr)
      elseif (name == 'kuw1') then
         kkuw1 = kuw1
         call c_get_pointer(kkuw1,ptr)

      elseif (name == 'maxlevels') then
         mmaxlevels = maxlevels
         call c_get_pointer(mmaxlevels,ptr)

      elseif (name == 'nfield_divf') then
         nnfield_divf = nfield_divf
         call c_get_pointer(nnfield_divf,ptr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DATA FROM MODULE WORKSPACE !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      elseif (name == 'work') then
         call c_get_pointer(work,ptr)
         
      elseif (name == 'interp_mask_work') then
         call c_get_pointer(interp_mask_work,ptr)

      elseif (name == 'interp_mask_work_res') then
         call c_get_pointer(interp_mask_work_res,ptr)

      elseif (name == 'work1') then
         call c_get_pointer(work1,ptr)

      elseif (name == 'cell_vol_w') then
         call c_get_pointer(cell_vol_w,ptr)

      elseif (name == 'output_dir') then
         call c_get_pointer(output_dir,ptr)

      else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IF THE NAME ISN'T FOUND, PRINT AN ERROR MESSAGE !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         print *,' ERROR in fortran_get_pointer: variable ',name, & 
     &           ' not found'
      end if

      return
      end subroutine fortran_get_pointer
