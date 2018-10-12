!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


!!****f* source/amr_1blk_cc_prol_user
!! NAME
!!
!!   amr_1blk_cc_prol_user
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_user(????)
!!
!! ARGUMENTS
!!
!!   Arguments to this routine are user defined.
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
!! RETURNS
!!
!!   Returns whatever the user defines.
!!
!! DESCRIPTION
!!
!!   This is a stub routine and is meant to be a place holder to allow
!!   a user to write their own prologation routine for cell centered 
!!   data.  You can use on of the *prol* routines as a guide to writing 
!!   your own, application specific prolongation routine.
!!
!!   NOTE: To use this feature you must define interp_mask_unk to be >= 20.
!!
!!   NOTE2: Use one of the other routines which are provided for doing this
!!   operation as an example.
!!
!!***

      Subroutine amr_1blk_ec_prol_user()

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use prolong_arrays

      Implicit None

!-----Include Statements
      Include 'mpif.h'

      Return
      End Subroutine amr_1blk_ec_prol_user
