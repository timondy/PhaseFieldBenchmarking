!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source amr_restrict_fc_fun
!! NAME
!! 
!!   amr_restrict_fc_fun
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_fc_fun(recv,temp,icoord)
!!   Call amr_restrict_fc_fun(real array, real array, integer)
!!
!! ARGUMENTS
!!
!!   Real, Intent(in)    :: recv(:,:,:,:) face centered data to restrict
!!   Real, Intent(inout) :: temp(:,:,:,:) array holding restricted face-centered data
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
!!   amr_restrict_fc_genorder
!!   amr_restrict_fc_user
!!
!! RETURNS
!!
!!   Restricted data returned in arrays 'temp'.
!!
!! DESCRIPTION
!!
!!   This routine performs a default or user defined reduction operation on the 
!!   array recv and returns the result in temp.
!!
!!   Note that this does not update guard cell elements of temp.
!!
!!   Also note that we use stride 2 along each dimension when computing
!!   reduced data values on cell faces, so not all values of temp
!!   have been updated.
!!
!! AUTHORS
!!
!!    Written :     Peter MacNeice          July 1997
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_fc_fun(recv,temp,icoord)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata

      Use paramesh_interfaces, Only : amr_restrict_fc_genorder,        & 
                                      amr_restrict_fc_user

      Implicit None

!-----Inout/Output arguments
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      integer, intent(in)    :: icoord

!-----Local variables
      integer :: order, ivar

!-----Begin executable code.

      Do ivar = 1, nbndvar

         If (icoord == 1) Then
            order = interp_mask_facex_res(ivar)
         ElseIf (icoord == 2) Then
            order = interp_mask_facey_res(ivar)
         ElseIf (icoord == 3) Then
            order = interp_mask_facez_res(ivar)
         End if

         If (order < 20) Then

            If (order <= 0 .or. order > 5) order = 1
            Call amr_restrict_fc_genorder(recv,temp,icoord,order,ivar)

         Else

            Call amr_restrict_fc_user()

         End if

      End Do  ! End Do ivar = 1,nbndvar

      Return
      End Subroutine amr_restrict_fc_fun

