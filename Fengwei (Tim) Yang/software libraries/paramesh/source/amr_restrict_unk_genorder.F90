!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_unk_genorder
!! NAME
!!
!!   amr_restrict_unk_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_genorder (datain, dataout, order, ivar)
!!   Call amr_restrict_unk_genorder (real array, real array, integer, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned
!!   Integer, Intent(in)    :: order order of interpolating polynomial to use
!!   Integer, Intent(in)    :: ivar  variable number in unk to restrict
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
!! CALLS
!! 
!! RETURNS
!!
!!   Restricted data returned in array 'dataout'.  
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   cell centered data stored in 'unk'.  It uses a lagrange polynomial
!!   interpolation scheme.  Also, the interpolation stencil is automatically
!!   shifted to avoid using data in guardcells. CAUTION: you must realize that
!!   use of this routine with 'order' set to values higher than 1 MAY 
!!   cause asymmetric interpolation and your results may lose symmetry.
!!
!!   Data is passed in via the array 'datain' and returned in the array 
!!   'dataout'.  The order of the interpolating polynomial is also passed 
!!   in the variable 'order' and can take on value ranging from 1 to 5.  
!!   The last argument 'ivar' specifies which variable in 'unk' to apply 
!!   the interpolation to.
!!
!! AUTHORS
!!
!!   Written :     Kevin Olson          March 2004
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_restrict_unk_genorder(datain,dataout,order,ivar) 

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Input/Output arguments.
      Real,    Intent(in)    :: datain(:,:,:,:)
      Real,    Intent(inout) :: dataout(:,:,:,:)
      Integer, Intent(in)    :: order, ivar

!-----Local arrays and variables.
      Real       :: xi, xj, www
      Real, Save :: weight(5,3,-4:5)
      Integer :: i,j,k
      Integer :: i0, j0, k0, is, js, ks, iw, jw, kw
      Integer :: iii, jjj, kkk
      Integer, Save :: iparmin,iparmax
      Integer, Save :: jparmin,jparmax
      Integer, Save :: kparmin,kparmax
      Integer :: istart, jstart, kstart
      Integer :: iend, jend, kend
      Integer :: order2
      Logical, Save :: first = .True.

!-----Being executable code.

      If (first) then

        first = .False.

        Do order2 = 1, 5

!----left

          xi = 0.-.5
          Do i = 0,order2
            weight(order2,1,i) = 1.
            xj = 0.-.5
            Do j = 0,order2
              If (i .ne. j) Then
                weight(order2,1,i) =                                    & 
                    weight(order2,1,i)*(0.-xj)/(xi-xj)
              End If
              xj = xj + 1.
            End Do
            xi = xi + 1.
          End Do

!-----middle

          istart = -int(order2/2)
          iend = istart + order2
          xi = real(istart)-.5
          Do i = istart,iend
            weight(order2,2,i) = 1.
            xj = real(istart)-.5
            Do j = istart,iend
              If (i .ne. j) Then
                weight(order2,2,i) =                                    & 
                    weight(order2,2,i)*(0.-xj)/(xi-xj)
              End If
              xj = xj + 1.
            End Do
            xi = xi + 1.
          End Do

!-----right

          istart = -order2 + 1
          iend = istart + order2
          xi = real(istart)-.5
          Do i = istart,iend
            weight(order2,3,i) = 1.
            xj = real(istart)-.5
            Do j = istart,iend
              If (i .ne. j) Then
                weight(order2,3,i) =                                    & 
                    weight(order2,3,i)*(0.-xj)/(xi-xj)
              End If
              xj = xj + 1.
            End Do
            xi = xi + 1.
          End Do

        End Do

        iparmin = 1+nguard
        iparmax = nxb+nguard
        jparmin = 1+nguard*k2d
        jparmax = nyb+nguard*k2d
        kparmin = 1+nguard*k3d
        kparmax = nzb+nguard*k3d

      End If  ! End If (first)

!-----Loops over parent cells which are restricted to
!CEG Tried zeroing all in one go - slower
!      dataout(ivar,:,:,:) = 0.

      Do k0 = kparmin,kparmax,2
        Do j0 = jparmin,jparmax,2
          Do i0 = iparmin,iparmax,2

! CEG Bizarrely seems faster to do indiviual elements rather than zeroing them all as tried above 
            dataout(ivar,i0,j0,k0) = 0.

            If (ndim == 3) Then
              If (k0 == kparmin) Then
                kstart = 0
                kw = 1
              Else if (k0 == kparmax-1) Then
                kstart = -order + 1
                kw = 3
              Else
                kstart = -int(order/2)
                kw = 2
              End If
              ks = k0+kstart
              kend = kstart + order
            Else
              ks     = 1
              kstart = 1
              kend   = 1
            End If  ! End If (ndim == 3)

            If (ndim >= 2) Then
              If (j0 == jparmin) Then
                jstart = 0
                jw = 1
              Else if (j0 == jparmax-1) Then
                jstart = -order + 1
                jw = 3
              Else
                jstart = -int(order/2)
                jw = 2
              End If
              js = j0+jstart
              jend = jstart + order
            Else
              js     = 1
              jstart = 1
              jend   = 1
            End If  ! End If (ndim >= 2)

            If (i0 == iparmin) Then
              istart = 0
              iw = 1
            Else if (i0 == iparmax-1) Then
              istart = -order + 1
              iw = 3
            Else
              istart = -int(order/2)
              iw = 2
            End If  ! End If (i0 == iparmin)
            is = i0+istart
            iend = istart + order

!-------Loops over child cells that are restricted from to parent point
!-------(i0, j0, k0)
            k = ks
            Do kkk = kstart,kend
              j = js
              Do jjj = jstart,jend
                i = is
                Do iii = istart,iend

                  If (ndim == 1) Then
                    www = weight(order,iw,iii)
                  Else if (ndim == 2) Then
                    www = weight(order,iw,iii)*                        & 
                          weight(order,jw,jjj)
                  Else if (ndim == 3) Then
                    www = weight(order,iw,iii)*                        & 
                          weight(order,jw,jjj)*                        & 
                          weight(order,kw,kkk)
                  End If

                  If (curvilinear_conserve) Then
                    dataout(ivar,i0,j0,k0) =                           & 
                      dataout(ivar,i0,j0,k0) +                         & 
                      (datain(ivar,i,j,k))
                  Else
                    dataout(ivar,i0,j0,k0) =                           & 
                      dataout(ivar,i0,j0,k0) +                         & 
                      (www*datain(ivar,i,j,k))
                  End If

                  i = i + 1
                End Do  ! End Do iii = istart,iend
                j = j + 1
              End Do  ! End Do jjj = jstart,jend
              k = k + 1
            End Do  ! End Do kkk = kstart,kend

          End Do  ! End Do i0 = iparmin,iparmax,2
        End Do  ! End Do j0 = jparmin,jparmax,2
      End Do  ! End Do k0 = kparmin,kparmax,2

      Return
      End Subroutine amr_restrict_unk_genorder
