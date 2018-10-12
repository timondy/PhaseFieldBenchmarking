!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_flux_conserve
!! NAME
!!
!!   amr_flux_conserve
!!
!! SYNOPSIS
!!
!!   Call amr_flux_conserve (mype, nsub)
!!   Call amr_flux_conserve (mype, nsub, flux_dir)
!!
!!   Call amr_flux_conserve (integer, integer, optional integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The Calling processor.
!!
!!   integer, intent(in) :: nsub          
!!     The current time subcycle. If this is 1 then this info is used to 
!!     reset the temporary boundary flux arrays to 0. This argument only has
!!     an effect if variable time steps are being used.
!!
!!   optional, integer, intent(in) :: flux_dir
!!     Option integer which selects which coordinate direction to apply
!!     the flux conservation operation to:
!!     If flux_dir = 1 -> x direction
!!        flux_dir = 2 -> y direction
!!        flux_dir = 3 -> z direction
!!     If this argument is not specified, then the default behaviour is
!!     to operate on all directions.  Using this argument can be useful
!!     for Strang-split schemes to improve performance.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!! 
!!   amr_flux_conserve_udt
!!   amr_flux_conserve_vdt
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit the fluxes stored in the arrays
!!   flux_x, flux_y, or flux_z are corrected at jumps in refinement.  This
!!   is either an averaging proceedure or a sum as selected by the user
!!   by adjusting the logical variables which control this behaviour
!!   that are read at runtime from the file 'amr_runtime_parameters'.
!!
!! DESCRIPTION
!!
!!   This is a wrapper routine which makes the appropriate Call to the
!!   routines which manage flux conservation at the boundaries between
!!   grid blocks of different refinement level.
!! 
!!   These routines get block boundary data from neighbors who are
!!   parents of leaf blocks. This is required in flux conserving schemes
!!   where the coarser block needs to use the same fluxes and mean pressures
!!   as will be used on the finer blocks across their shared boundary.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson for 
!!   directional guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_flux_conserve(mype,nsub,flux_dir)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use paramesh_interfaces, only : amr_flux_conserve_udt,           & 
                                      amr_flux_conserve_vdt

      Implicit None

!-----Input/Output arguments.
      Integer, intent(in)  ::  mype,nsub
      Integer, optional, intent(in) :: flux_dir

!-----Local variables
      Integer :: lb

!-----Begin executable code.
      If (lnblocks > 0) Then
      Do lb = 1,lnblocks

      If (nodetype(lb) == 1 .or. advance_all_levels) Then

!------Store fluxes in temporary storage
       tflux_x(:,:,:,:,lb) = flux_x(:,:,:,:,lb)
       tflux_y(:,:,:,:,lb) = flux_y(:,:,:,:,lb)
       tflux_z(:,:,:,:,lb) = flux_z(:,:,:,:,lb)

      End If

      End Do
      End If

      If (var_dt) Then
       Call amr_flux_conserve_vdt(mype,nsub) ! Called if variable dt
      Else
       Call amr_flux_conserve_udt(mype,flux_dir) ! Called if uniform dt
      End If

      Return
      End Subroutine amr_flux_conserve
