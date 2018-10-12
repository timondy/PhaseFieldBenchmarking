!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_guardcell_reset
!!
!! NAME
!!
!!   amr_1blk_guardcell_reset
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_guardcell_reset()
!!
!! ARGUMENTS
!!
!!  None
!!
!! INCLUDES
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! CALLS
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!! This routine resets some variables which manage the guardcell
!! filling operation when operating in 1-block mode. It should be called
!! from amr_initialize, and at a synchronization point which separates
!! guardcell filling of data with a given time stamp, and guardcell
!! filling at the next time stamp.
!!
!! AUTHORS
!!
!!   Peter MacNeice, February 1999.
!!
!!***

      Subroutine amr_1blk_guardcell_reset

      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Begin Executable Code

! reset id recording parent block which is currently in cache
          lnew_parent = .True.
          pcache_blk_u = -1
          pcache_pe_u  =  -1
          pcache_blk_w = -1
          pcache_pe_w  =  -1


      Return
      End Subroutine amr_1blk_guardcell_reset
