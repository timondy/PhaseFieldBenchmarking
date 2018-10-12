!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_unk_fun
!! NAME
!!
!!   amr_restrict_unk_fun
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_fun(datain, dataout)
!!   Call amr_restrict_unk_fun(real array, real array)
!!
!! ARGUMENTS
!!
!!   Real, Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real, Intent(inout) :: dataout(:,:,:,:) restricted data to return
!!   
!! INCLUDES
!! 
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_restrict_unk_genorder 
!!   amr_restrict_unk_user
!!
!! RETURNS
!!
!!   Restricted data returned in array 'dataout'.
!!
!! DESCRIPTION
!!   
!!   This routine performs restriction on the array datain and
!!   returns the result in dataout. Note that this does not update
!!   guard cell elements of dataout.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 1997
!!   Modified by Kevin Olson for high order restriction, 2004.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_unk_fun(datain,dataout)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use paramesh_interfaces, only : amr_restrict_unk_genorder,       & 
                                      amr_restrict_unk_user

      Implicit None

!-----Input/Output arguments.
      Real, Intent(in)    :: datain(:,:,:,:)
      Real, Intent(inout) :: dataout(:,:,:,:)

!-----Local variables.
      Integer :: ivar, order

!-----Begin Executable code.

      Do ivar = 1, nvar

         If (int_gcell_on_cc(ivar)) Then

         If (interp_mask_unk_res(ivar) < 20) Then

!-----------Call the default interpolation routine for interpolation 
            order = interp_mask_unk_res(ivar)
            If (order <=0 .or. order > 5) order = 1
            Call amr_restrict_unk_genorder(datain,dataout,order,ivar)

         ElseIf (interp_mask_unk_res(ivar) >= 20) Then

!-----------Call a user defined routine for restriction
            Call amr_restrict_unk_user()

         End If

         End If

      End Do  ! End Do ivar = 1, nvar

      Return
      End Subroutine amr_restrict_unk_fun




