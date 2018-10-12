!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_bcset_init
!! NAME
!!
!!   amr_bcset_init
!!
!! SYNOPSIS
!!
!!   Call amr_bcset_init()
!!
!! ARGUMENTS
!!
!!   No arguments.
!!
!! INCLUDES
!! 
!!   Noi includes.
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! CALLS
!!
!!   No calls made to other Paramesh routines.
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine sets index ranges which are later used when the
!!   boundary conditions are applied using amr_1blk_bcset.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 2001
!!
!!***

      Subroutine amr_bcset_init

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----base index ranges for unk, facevar^s, unk_e_^s and unk_n
 
      bc_index_i(1,1,:) = 1                        ! lo i on left boundary
      bc_index_i(2,1,:) = nguard                   ! hi i on left boundary
      bc_index_i(1,2,:) = nguard+1                 ! lo i of middle range
      bc_index_i(2,2,:) = nguard+nxb               ! hi i of middle range
      bc_index_i(1,3,:) = nguard+nxb+1             ! lo i on right boundary
      bc_index_i(2,3,:) = nxb+2*nguard             ! hi i on right boundary
 
      bc_index_j(1,1,:) = 1                        ! lo j on left boundary
      bc_index_j(2,1,:) = 1+(nguard-1)*k2d         ! hi j on left boundary
      bc_index_j(1,2,:) = nguard*k2d+1             ! lo j of middle range
      bc_index_j(2,2,:) = nguard*k2d+nyb           ! hi j of middle range
      bc_index_j(1,3,:) = (nguard+nyb)*k2d+1       ! lo j on right boundary
      bc_index_j(2,3,:) = nyb+2*nguard*k2d         ! hi j on right boundary
 
      bc_index_k(1,1,:) = 1                        ! lo k on left boundary
      bc_index_k(2,1,:) = 1+(nguard-1)*k3d         ! hi k on left boundary
      bc_index_k(1,2,:) = nguard*k3d+1             ! lo k of middle range
      bc_index_k(2,2,:) = nguard*k3d+nzb           ! hi k of middle range
      bc_index_k(1,3,:) = (nguard+nzb)*k3d+1       ! lo k on right boundary
      bc_index_k(2,3,:) = nzb+2*nguard*k3d         ! hi k on right boundary

!-----base index ranges for work
      
      bc_index_i(1,1,5) = 1                        ! lo i on left boundary
      bc_index_i(2,1,5) = nguard_work              ! hi i on left boundary
      bc_index_i(1,2,5) = nguard_work+1            ! lo i of middle range
      bc_index_i(2,2,5) = nguard_work+nxb          ! hi i of middle range
      bc_index_i(1,3,5) = nguard_work+nxb+1        ! lo i on right boundary
      bc_index_i(2,3,5) = nxb+2*nguard_work        ! hi i on right boundary

      bc_index_j(1,1,5) = 1                        ! lo j on left boundary
      bc_index_j(2,1,5) = 1+(nguard_work-1)*k2d    ! hi j on left boundary
      bc_index_j(1,2,5) = nguard_work*k2d+1        ! lo j of middle range
      bc_index_j(2,2,5) = nguard_work*k2d+nyb      ! hi j of middle range
      bc_index_j(1,3,5) = (nguard_work+nyb)*k2d+1  ! lo j on right boundary
      bc_index_j(2,3,5) = nyb+2*nguard_work*k2d    ! hi j on right boundary

      bc_index_k(1,1,5) = 1                        ! lo k on left boundary
      bc_index_k(2,1,5) = 1+(nguard_work-1)*k3d    ! hi k on left boundary
      bc_index_k(1,2,5) = nguard_work*k3d+1        ! lo k of middle range
      bc_index_k(2,2,5) = nguard_work*k3d+nzb      ! hi k of middle range
      bc_index_k(1,3,5) = (nguard_work+nzb)*k3d+1  ! lo k on right boundary
      bc_index_k(2,3,5) = nzb+2*nguard_work*k3d    ! hi k on right boundary

      Return
      End Subroutine amr_bcset_init
