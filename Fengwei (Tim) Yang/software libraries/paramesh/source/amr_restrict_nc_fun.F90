!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_nc_fun
!! NAME
!!
!!   amr_restrict_nc_fun
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_nc_fun(datain,dataout)
!!   Call amr_restrict_nc_fun(real array,real array)
!!
!! ARGUMENTS
!!
!!   Real, Intent(in) :: datain(:,:,:,:)
!!   Real, Intent(inout) :: dataout(:,:,:,:)
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
!!   amr_restrict_nc_genorder
!!   amr_restrict_nc_user
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
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_restrict_nc_fun(datain,dataout)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use paramesh_interfaces, Only : amr_restrict_nc_genorder,        & 
                                      amr_restrict_nc_user

      Implicit None

!-----Input/Output arguements.
      Real, Intent(in)    :: datain(:,:,:,:)
      Real, Intent(inout) :: dataout(:,:,:,:)

!-----Local variables
      integer :: ivar

!-----Begin executable code.

      If (nvarcorn > 0) Then

      Do ivar = 1, nvarcorn

       If (interp_mask_nc_res(ivar) < 20) Then

          Call amr_restrict_nc_genorder(datain,dataout,ivar)

       Else

          Call amr_restrict_nc_user()

       End If

      End Do  ! End Do ivar = 1, nvarcorn

      End If  ! End If (nvarcorn > 0)

      Return
      End Subroutine amr_restrict_nc_fun
