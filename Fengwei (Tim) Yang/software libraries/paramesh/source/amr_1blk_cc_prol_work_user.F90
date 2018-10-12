!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_work_user
!! NAME
!!
!!   amr_1blk_cc_prol_work_user
!!
!! SYNOPSIS
!! 
!!   Call amr_1blk_cc_prol_work_user(???)
!!
!! ARGUMENTS
!!
!!   User defined.
!!
!! INCLUDES
!! 
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   prolong_arrays
!!
!! CALLS
!!
!!   No calls to other PARAMESH routines are made.
!!
!! DESCRIPTION
!!
!!   This is a stub routine and is meant to be a place holder to allow
!!   a user to write their own prolocation routine for cell centered 
!!   data.  
!!
!!   NOTE: To use this feature you must define interp_mask_unk to be >= 20.
!!
!!   NOTE2: Use one of the other routines which are provided for doing this
!!   operation as an example.
!!
!! AUTHORS
!!
!!   Peter MacNeice, Kevin Olson.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_cc_prol_work_user()

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use prolong_arrays

      Implicit None

      Include 'mpif.h'

!-----Begin Executable code

      Return
      End Subroutine amr_1blk_cc_prol_work_user
