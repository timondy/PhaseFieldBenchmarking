!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_close
!! NAME
!!
!!   amr_close
!!
!! SYNOPSIS
!!
!!   Call amr_close()
!!
!! AGUMENTS
!!
!!   No arguments
!!
!! INCLUDES
!!
!!   paramesh_preprocessorh.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   mpi_morton
!!   tree
!!   prolong_arrays
!!   timings
!!
!! CALLS
!! 
!!   comm_finish
!!
!! RETURNS
!!
!!   Does not return anything.
!!
!! DESCRIPTION
!!
!!   This subroutine closes the PARAMESH amr package and deallocates all arrays which 
!!   were first allocated be a call to 'amr_initialize'.
!!
!! AUTHORS
!!
!!   Peter MacNeice and Kevin Olson.
!!
!!***

# include "paramesh_preprocessor.fh"
      Subroutine amr_close

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use mpi_morton
      Use tree
      Use prolong_arrays
      Use timings

      Use paramesh_interfaces, Only : comm_finish

      Implicit None

!----------------------------------------------------------------
      integer tmpset(8)

!-----Deallocate cell-centered variables

      Deallocate(unk)
      Deallocate(interp_mask_unk)
      Deallocate(interp_mask_unk_res)
      Deallocate(gcell_on_cc_pointer)
      Deallocate(gcell_on_cc)
      Deallocate(int_gcell_on_cc)
      Deallocate(checkp_on_cc)
      If (nvar > 0) Then
       Deallocate(unk1)
       Deallocate(gt_unk)
       If (var_dt .or. pred_corr) Then
          Deallocate(t_unk)
       End If  ! End of If (var_dt .or. pred_corr)
      End If  ! End If (nvar > 0)

!-----Deallocate face variables

      Deallocate(facevarx)
      Deallocate(facevary)
      Deallocate(facevarz)

      If (nfacevar > 0) Then
       Deallocate(facevarx1)
       Deallocate(facevary1)
       Deallocate(facevarz1)
       Deallocate(gt_facevarx)
       Deallocate(gt_facevary)
       Deallocate(gt_facevarz)
       If (var_dt .or. pred_corr) Then
        Deallocate(tfacevarx)
        Deallocate(tfacevary)
        Deallocate(tfacevarz)
       End if  ! End if (var_dt .or. pred_corr)
      End If  ! End If (nfacevar > 0)

      Deallocate(gcell_on_fc_pointer)
      Deallocate(gcell_on_fc)
      Deallocate(int_gcell_on_fc)
      Deallocate(interp_mask_facex)
      Deallocate(interp_mask_facey)
      Deallocate(interp_mask_facez)
      Deallocate(interp_mask_facex_res)
      Deallocate(interp_mask_facey_res)
      Deallocate(interp_mask_facez_res)
      Deallocate(checkp_on_fc)

!-----Deallocate edge variables

      Deallocate(unk_e_x)
      Deallocate(unk_e_y)
      Deallocate(unk_e_z)

      If (nvaredge > 0) then
       Deallocate(unk_e_x1)
       Deallocate(unk_e_y1)
       Deallocate(unk_e_z1)
       Deallocate(gt_unk_e_x)
       Deallocate(gt_unk_e_y)
       Deallocate(gt_unk_e_z)
       If (var_dt .or. pred_corr) Then
         Deallocate(t_unk_e_x)
         Deallocate(t_unk_e_y)
         Deallocate(t_unk_e_z)
       End If  ! End If (var_dt .or. pre_corr)
      End If  ! End If (nvaredge > 0)

      Deallocate(interp_mask_ec)
      Deallocate(interp_mask_ec_res)
      Deallocate(gcell_on_ec)
      Deallocate(gcell_on_ec_pointer)
      Deallocate(int_gcell_on_ec)
      Deallocate(checkp_on_ec)

!-----Deallocate corner variables

      Deallocate(unk_n)

      If (nvarcorn > 0) then
       Deallocate(unk_n1)
       Deallocate(gt_unk_n)
       If (var_dt .or. pred_corr) Then
          Deallocate(t_unk_n)
       End If  ! End of If (var_dt .or. pred_corr)
      End If  ! End If (nvarcorn > 0) 

      Deallocate(interp_mask_nc)
      Deallocate(interp_mask_nc_res)
      Deallocate(gcell_on_nc)
      Deallocate(gcell_on_nc_pointer)
      Deallocate(int_gcell_on_nc)
      Deallocate(checkp_on_nc)

!-----Deallocate variable timestep support arrays

      Deallocate(time_loc)
      Deallocate(ldtcomplete)

!-----Deallocate flux fix-up arrays

      Deallocate(flux_x)
      Deallocate(flux_y)
      Deallocate(flux_z)
      Deallocate(tflux_x)
      Deallocate(tflux_y)
      Deallocate(tflux_z)
      If (var_dt) Then
         Deallocate(ttflux_x)
         Deallocate(ttflux_y)
         Deallocate(ttflux_z)
      End If  ! End of If (var_dt)

!-----Deallocate edge fix-up arrays

      Deallocate(bedge_facex_y)
      Deallocate(bedge_facex_z)
      Deallocate(bedge_facey_x)
      Deallocate(bedge_facey_z)
      Deallocate(bedge_facez_x)
      Deallocate(bedge_facez_y)
      Deallocate(recvarx1e)
      Deallocate(recvary1e)
      Deallocate(recvarz1e)
      Deallocate(recvarx2e)
      Deallocate(recvary2e)
      Deallocate(recvarz2e)
      Deallocate(tbedge_facex_y)
      Deallocate(tbedge_facex_z)
      Deallocate(tbedge_facey_x)
      Deallocate(tbedge_facey_z)
      Deallocate(tbedge_facez_x)
      Deallocate(tbedge_facez_y)
      If (var_dt) Then
         Deallocate(ttbedge_facex_y)
         Deallocate(ttbedge_facex_z)
         Deallocate(ttbedge_facey_x)
         Deallocate(ttbedge_facey_z)
         Deallocate(ttbedge_facez_x)
         Deallocate(ttbedge_facez_y)
      End If  ! End of If (var_dt)

      If (curvilinear) Then
         Deallocate(cell_vol)
         Deallocate(cell_area1)
         Deallocate(cell_area2)
         Deallocate(cell_area3)
         Deallocate(cell_leng1)
         Deallocate(cell_leng2)
         Deallocate(cell_leng3)
      End If  ! End of If (curvilinear)

      Deallocate(recvarxf)
      Deallocate(recvaryf)
      Deallocate(recvarzf)
      Deallocate(bndtempx1)
      Deallocate(bndtempy1)
      Deallocate(bndtempz1)

!-----Deallocate tree data

      Deallocate(neigh)
      Deallocate(child)
      Deallocate(which_child)
      Deallocate(parent)
      Deallocate(lrefine)
      Deallocate(nodetype)
      Deallocate(empty)
      Deallocate(bflags)
      Deallocate(newchild)
      Deallocate(derefine)
      Deallocate(refine)
      Deallocate(stay)
      Deallocate(work_block)
      Deallocate(coord)
      Deallocate(bsize)
      Deallocate(bnd_box)
      Deallocate(level_cell_sizes)
      Deallocate(laddress)
      Deallocate(surr_blks)
#ifdef SAVE_MORTS
      Deallocate(surr_morts)
#endif
      Deallocate(boundary_box)
      Deallocate(boundary_index)

! workspace data

      Deallocate(work)
      Deallocate(interp_mask_work)
!next added by CEG as forgotten
      deallocate(interp_mask_work_res)
      If (nvar_work > 0) then
       Deallocate(recvw)
       Deallocate(sendw)
       Deallocate(tempw)
       Deallocate(work1)
       Deallocate(recvw1)
       Deallocate(tempw1)
      End If  ! End If (nvar_work > 0)
      If (curvilinear) Then
         Deallocate(cell_vol_w)
      End If  ! End of If (curvilinear)

! morton data

      Deallocate(mortonbnd)
      Deallocate(laddress_guard)
      Deallocate(laddress_prol)
      Deallocate(laddress_flux)
      Deallocate(laddress_restrict)

! prolong_arrays data

      Deallocate(prol_dx)
      Deallocate(prol_dy)
      Deallocate(prol_dz)
      Deallocate(prol_indexx)
      Deallocate(prol_indexy)
      Deallocate(prol_indexz)
      Deallocate(prol_f_dx)
      Deallocate(prol_f_dy)
      Deallocate(prol_f_dz)
      Deallocate(prol_f_indexx)
      Deallocate(prol_f_indexy)
      Deallocate(prol_f_indexz)
      Deallocate(prolw_dx)
      Deallocate(prolw_dy)
      Deallocate(prolw_dz)
      Deallocate(prolw_indexx)
      Deallocate(prolw_indexy)
      Deallocate(prolw_indexz)


      Deallocate(ladd_strt)
      Deallocate(ladd_end)

      Deallocate(timer_amr_1blk_to_perm)


      Deallocate(i_divf_fc_vars)


! CEG extra clean up calls
      call mpi_amr_comm_free()

      ! allocs from mpi_source/mpi_amr_store_comm_info.F90   
      deallocate(to_be_received)
      deallocate(to_be_sent)
      deallocate(pe_source_guard)
      deallocate(commatrix_guard)
      deallocate(pe_source_prol)
      deallocate(commatrix_prol)
      deallocate(pe_source_flux)
      deallocate(commatrix_flux)
      deallocate(pe_source_restrict)
      deallocate(commatrix_restrict)
      deallocate(bc_block_neighs)
      deallocate(bc_block_neighs_send)

      ! allocs from mpi_source/mpi_pack_blocks.F90   
      deallocate(mess_segment_loc)

      deallocate(block_starts)

      ! deallocate MPI types
      call amr_redist_blk_mpitypeset(tmpset)

! Call the machine/software environment specific closure routine.

      Call comm_finish()

      Return
      End Subroutine amr_close



