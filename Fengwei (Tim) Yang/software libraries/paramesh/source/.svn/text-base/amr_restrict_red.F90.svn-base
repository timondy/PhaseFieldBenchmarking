!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_red
!! NAME
!! 
!!   amr_restrict_red
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_red(icoord)
!!   Call amr_restrict_red(integer)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in) :: icoord  Coordinate direction of face centered data
!!                                  to restrict
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
!! DESCRIPTION
!!
!!   This routine performs a reduction operation on the 
!!   array recvarx(y)(z) and returns the result in bndtempx(y)(z).
!!   These data arrays are defined on block boundaries only and is
!!   called at part of the overall flux-fixup function.
!!
!!   Note that this does not update guard cell elements of bndtempx(y)(z).
!!
!!   Also note that we use stride 2 along each dimension when computing
!!   reduced data values on block faces, so not all values of dataout
!!   have been updated.
!!
!!   This particular version is only appropriate for 2nd order schemes 
!!   using linear interpolation with even number of mesh points along 
!!   each block axis.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          February 1997
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f

      Subroutine amr_restrict_red(icoord)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Input/Output arguements.
      integer, intent(in)    :: icoord

!-----Local variables.
      integer :: i,j,k,ivar
      integer :: nguard0

!-----Begin executable codes.

      nguard0 = nguard*npgs

      If (icoord == 1) then                         ! x-comp of flux

         Do k=1+nguard0*k3d,nzb+nguard0*k3d,2
           Do j=1+nguard0*k2d,nyb+nguard0*k2d,2
             Do i=1,2
               Do ivar=1,nfluxes
                 bndtempx1(ivar,i,j,k) = (                             & 
                              recvarxf(ivar,i,j,k) +                   & 
                              recvarxf(ivar,i,j+k2d,k) +               & 
                              recvarxf(ivar,i,j,k+k3d) +               & 
                              recvarxf(ivar,i,j+k2d,k+k3d))            & 
                             *red_f
               End Do
             End Do
           End Do
         End Do

      ElseIf (icoord == 2) then                     ! y-comp of flux

        Do k=1+nguard0*k3d,nzb+nguard0*k3d,2
          Do j=1,2
            Do i=1+nguard0,nxb+nguard0,2
              Do ivar=1,nfluxes
                bndtempy1(ivar,i,j,k) = (                              & 
                             recvaryf(ivar,i,j,k) +                    & 
                             recvaryf(ivar,i+1,j,k) +                  & 
                             recvaryf(ivar,i,j,k+k3d) +                & 
                             recvaryf(ivar,i+1,j,k+k3d))               & 
                            *red_f
               End Do
             End Do
           End Do
         End Do

      ElseIf (icoord == 3) then                      ! z-comp of flux

        Do k=1,2
          Do j=1+nguard0,nyb+nguard0,2
            Do i=1+nguard0,nxb+nguard0,2
              Do ivar=1,nfluxes
                bndtempz1(ivar,i,j,k) = (                              & 
                             recvarzf(ivar,i,j,k) +                    & 
                             recvarzf(ivar,i+1,j,k) +                  & 
                             recvarzf(ivar,i,j+1,k) +                  & 
                             recvarzf(ivar,i+1,j+1,k))                 & 
                            *red_f
               End Do
             End Do
           End Do
         End Do

      End If  ! End If (icoord == 1)

      Return
      End Subroutine amr_restrict_red
