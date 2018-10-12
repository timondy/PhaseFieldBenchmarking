!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_work_genorder
!! NAME
!!
!!   amr_restrict_work_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_work_genorder(datainw,dataoutw,order)
!!   Call amr_restrict_work_genorder(real array, real array, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict                         
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned   
!!   Integer, Intent(in)    :: order order of interpolating polynomial
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!
!! CALLS
!!
!! RETURNS
!!
!!   Restricted data returned in array 'dataoutw'.
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   cell centered data stored in 'work'.  It uses a lagrange polynomial
!!   interpolation scheme.  Also, the interpolation stencil is automatically
!!   shifted to avoid using data in guardcells. CAUTION: you must realize that
!!   use of this routine with 'order' set to values higher than 1 MAY
!!   cause asymmetric interpolation and your results may lose symmetry.
!!
!!   Data is passed in in the array 'datainw' and returned in the array
!!   'dataoutw'.  The order of the interpolating polynomial is also passed
!!   in the variable 'order' and can take on value ranging from 1 to 5.
!!   The last argument 'ivar' specifies which variable in 'work' to apply
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

      Subroutine amr_restrict_work_genorder(datainw,dataoutw,iopt,order)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use workspace

      Implicit None

!-----Input/Output arguments.
      Real,    Intent(in)    :: datainw(:,:,:)
      Real,    Intent(inout) :: dataoutw(:,:,:)
      Integer, Intent(in)    :: iopt, order

!-----Local arrays and variables.
      Real    :: xi, xj, www
      Real,    Save :: weight(5,3,-4:5)
      Integer :: i,j,k
      Integer :: i0, j0, k0, is, js, ks
      Integer :: iw, jw, kw, iii, jjj, kkk
      Integer, Save :: iparmin,iparmax
      Integer, Save :: jparmin,jparmax
      Integer, Save :: kparmin,kparmax
      Integer :: istart, jstart, kstart
      Integer :: iend, jend, kend
      Integer :: order2
      Logical, Save :: first = .True.

!-----Begin executable code.

      If (first) then

      first = .False.

      Do order2 = 1, 5

!-----left

      xi = 0.-.5
      Do i = 0,order2
         weight(order2,1,i) = 1.
         xj = 0.-.5
         Do j = 0,order2
            If (i .ne. j) then
               weight(order2,1,i) =                                    & 
                    weight(order2,1,i)*(0.-xj)/(xi-xj)
            end If
            xj = xj + 1.
         End Do
         xi = xi + 1.
      End Do

!-----middle

      istart = -int(order2/2)
      iend = istart + order2
      xi = Real(istart)-.5
      Do i = istart,iend
         weight(order2,2,i) = 1.
         xj = Real(istart)-.5
         Do j = istart,iend
            If (i .ne. j) then
               weight(order2,2,i) =                                    & 
                   weight(order2,2,i)*(0.-xj)/(xi-xj)
            end If
            xj = xj + 1.
         End Do
         xi = xi + 1.
      End Do

!-----right

      istart = -order2 + 1
      iend = istart + order2
      xi = Real(istart)-.5
      Do i = istart,iend
         weight(order2,3,i) = 1.
         xj = Real(istart)-.5
         Do j = istart,iend
            If (i .ne. j) then
               weight(order2,3,i) =                                    & 
                    weight(order2,3,i)*(0.-xj)/(xi-xj)
            end If
            xj = xj + 1.
         End Do
         xi = xi + 1.
      End Do

      End Do  ! End Do order2 = 1, 5

      iparmin = 1+nguard_work
      iparmax = nxb+nguard_work
      jparmin = 1+nguard_work*k2d
      jparmax = nyb+nguard_work*k2d
      kparmin = 1+nguard_work*k3d
      kparmax = nzb+nguard_work*k3d

      End If  ! End If (first)

      Do k0 = kparmin,kparmax,2
         Do j0 = jparmin,jparmax,2
            Do i0 = iparmin,iparmax,2

               dataoutw(i0,j0,k0) = 0.

               If (ndim == 3) then
                  If (k0 == kparmin) then
                     kstart = 0
                     kw = 1
                  Else If (k0 == kparmax-1) then
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
               End If

               If (ndim >= 2) then
                  If (j0 == jparmin) then
                     jstart = 0
                     jw = 1
                  Else If (j0 == jparmax-1) then
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
               End If

               If (i0 == iparmin) then
                  istart = 0
                  iw = 1
               Else If (i0 == iparmax-1) then
                  istart = -order + 1
                  iw = 3
               Else
                  istart = -int(order/2)
                  iw = 2
               End If
               is = i0+istart
               iend = istart + order

               k = ks
               Do kkk = kstart,kend
                  j = js
                  Do jjj = jstart,jend
                     i = is
                     Do iii = istart,iend

                        If (ndim == 1) then
                           www = weight(order,iw,iii)
                        Else If (ndim == 2) then
                           www = weight(order,iw,iii)*                        & 
                                 weight(order,jw,jjj)
                        Else If (ndim == 3) then
                           www = weight(order,iw,iii)*                        & 
                                 weight(order,jw,jjj)*                        & 
                                 weight(order,kw,kkk)
                        End If

                        If (curvilinear_conserve) then
                           dataoutw(i0,j0,k0) =                               & 
                             dataoutw(i0,j0,k0) +                             & 
                             (datainw(i,j,k))
                        Else
                           dataoutw(i0,j0,k0) =                               & 
                             dataoutw(i0,j0,k0) +                             & 
                             (www*datainw(i,j,k))
                        End If

                        i = i + 1
                     End Do  ! End Do iii = istart,iend
                     j = j + 1
                  End Do  ! End Do jjj = jstart,jend
                  k = k + 1
               End Do  ! End Do kkk = kstart,kend

            End Do  ! End Do i0 = iparmin,iparmax,2
         End Do  ! End Do j0 = jparmin.jparmax,2
      End Do  ! End Do k0 = kparmin,kparmax,2

      Return
      End Subroutine amr_restrict_work_genorder
