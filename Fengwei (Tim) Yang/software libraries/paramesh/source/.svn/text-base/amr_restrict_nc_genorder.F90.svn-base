!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_nc_genorder
!! NAME
!!
!!   amr_restrict_nc_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_nc_genorder(datain, dataout, ivar)
!!   Call amr_restrict_nc_genorder(real array, real array, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)   data to restrict
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which has been restricted
!!   Integer, Intent(in)    :: ivar              variable to restrict
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!! 
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! RETURNS
!!
!!   Restricted node-centered data are returned in the array 'dataout'.
!!
!! DESCRIPTION
!!
!!   This routine performs restriction on the array datain and
!!   returns the result in dataout via direct injection.  This is adequate
!!   for most applications since the locations of data points on parent blocks
!!   are coincident with points on their child blocks.
!!   
!! AUTHORS
!!
!!   Written :     Kevin Olson          March 2004
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_nc_genorder(datain,dataout,ivar)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Input/Output Arguments.
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: ivar

!-----Local variables.
      integer :: i,j,k

!-----Begin executable code.

      Do k=1+nguard*k3d,nzb+nguard*k3d+k3d
      Do j=1+nguard*k2d,nyb+nguard*k2d+k2d
      Do i=1+nguard,nxb+nguard+1

         dataout(ivar,i,j,k) = datain(ivar,i,j,k)

      End Do
      End Do
      End Do

      Return
      End subroutine amr_restrict_nc_genorder
