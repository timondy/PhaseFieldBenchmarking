!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_work_fun
!! NAME
!!
!!   source/amr_restrict_work_fun
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_work_fun(datainw,dataout2,iopt)
!!   Call amr_restrict_work_fun(real array,real array,integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datainw(:,:,:)  parent data to restrict
!!   Real,    Intent(inout) :: dataoutw(:,:,:) restricted data which is returned
!!   Integer, Intent(in)    :: iopt  iopt-1 indicate which variable in work to
!!                                   restrict
!! INCLUDES
!! 
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   amr_restrict_work_genorder 
!!
!! CALLS
!!
!!   amr_restrict_work_genorder
!!   amr_restrict_work_user
!!
!! RETURNS
!!
!!   Restricted data returned in array 'dataoutw'
!!
!! DESCRIPTION
!!
!!   This routine performs restriction on the array datainw and
!!   returns the result in dataoutw. Note that this does not updata
!!   guard cell elements of dataoutw.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 1997
!!   Modified by Kevin Olson for high order restriction, 2004.
!!
!!***

!!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_work_fun(datainw,dataoutw,iopt)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use paramesh_interfaces, Only : amr_restrict_work_genorder,      & 
                                      amr_restrict_work_user

      Implicit None

!-----Input/Output arguments.
      Real,    Intent(in)    :: datainw(:,:,:)
      Real,    Intent(inout) :: dataoutw(:,:,:)
      integer, Intent(in)    :: iopt

!-----Local arguments.
      Integer :: order

!-----Begin executable code.
      If (interp_mask_work_res(iopt-1) < 20) Then
         
!--------Default interpolation routine for restriction of 'work' array
         order = interp_mask_work_res(iopt-1)
         If (order <= 0 .or. order > 5) order = 1
         Call amr_restrict_work_genorder (datainw,dataoutw,iopt,order)

      Else

!-------User defined routine for restriction of 'work' array 
        Call amr_restrict_work_User()

      End If

      Return
      End Subroutine amr_restrict_work_fun


