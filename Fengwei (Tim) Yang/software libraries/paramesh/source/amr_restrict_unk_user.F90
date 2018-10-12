!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_unk_user
!! NAME
!!
!!   amr_restrict_unk_user
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_user()
!!
!! ARGUMENTS
!!
!!   None or user defined.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!! CALLS
!!
!! DESCRIPTION
!!
!!   This is a stub routine to hold the place of a user written 
!!   interpolation routine to used during restriction of cell centered data 
!!   from fine to course meshes.
!! 
!! AUTHORS
!!
!!   YOU, the users of Paramesh.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_unk_user()

      Return
      End Subroutine amr_restrict_unk_user
