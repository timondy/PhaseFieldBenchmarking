!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_ec_fun
!! NAME
!!   
!!   amr_restrict_ec_fun
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_ec_fun(recv,temp,icoord)
!!   Call amr_restrict_ec_fun(real array,real array,integer)
!!   
!! ARGUMENTS
!!
!!   Real,    intent(in)    :: recv(:,:,:,:)  edge centered data to restrict
!!   Real,    intent(inout) :: temp(:,:,:,:)  array holding restricted edge centered
!!                                              data
!!   Integer, intent(in)    :: icoord         which coordinate to restrict
!!
!! INCLUDES
!!
!!   paramesh_preprocessorh.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_restrict_ec_genorder
!!   amr_restrict_ec_user
!!
!! RETURNS
!!
!!   Restricted data returned in array 'temp'.
!!
!! DESCRIPTION
!!
!!   This routine performs a restriction operation on the 
!!   array recv and returns the result in temp for edge centered data.
!!
!!   The restriction is done either by calling the default routine
!!   'amr_restrict_ec_genorder' or a user supplied routine
!!   'amr_restrict_ec_user'.  Which routine is used depends on the
!!   settings defined be the user in the array 'interp_mask_ec_res'.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          December 2000
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_ec_fun(recv,temp,icoord)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use paramesh_interfaces, only : amr_restrict_ec_genorder, & 
                                      amr_restrict_ec_user

      Implicit None

!-----Input/Output arguments.
      Real,    Intent(in)    :: recv(:,:,:,:)
      Real,    Intent(inout) :: temp(:,:,:,:)
      Integer, Intent(in)    :: icoord

!-----Local variables and arrays.
      Integer :: ivar, order

!-----Begin executable code.

      Do ivar = 1, nbndvare
         
         order = interp_mask_ec_res(ivar)
         
         If (order < 20) Then

            If (order <= 0 .or. order > 5) order = 1
            Call amr_restrict_ec_genorder(recv,temp,icoord,order,ivar)
            
         Else
            
            Call amr_restrict_ec_user()
            
         End If
         
      End Do

      Return
      End Subroutine amr_restrict_ec_fun
