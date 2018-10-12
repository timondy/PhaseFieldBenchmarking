!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_ec_genorder
!! NAME
!!
!!   amr_restrict_ec_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_ec_genorder(recv,temp,icoord,order,ivar)
!!   Call amr_restrict_ec_genorder(real array,real array,integer,integer,integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: recv(:,:,:,:)  Edge-centered data to restrict
!!   Real,    Intent(inout) :: temp(:,:,:,:)  Restricted edge-centered data
!!   Integer, Intent(in)    :: icoord         selects which edge to operate on
!!                                            iccord = 1, selects x-edge 
!!                                            iccord = 2, selects y-edge 
!!                                            iccord = 3, selects z-edge 
!!   Integer, Intent(in)    :: order          order of Lagrange polynomial to use
!!   Integer, Intent(in)    :: ivar           which variable in unk_e_x(y,z) to
!!                                            operate on.
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
!!   Restrict edge-centered data are returned in the array 'temp'.
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   edge centered data stored in 'unk_e_?'.  It uses a lagrange polynomial
!!   interpolation scheme.  Also, the interpolation stencil is automatically
!!   shifted to avoid using data in guardcells. CAUTION: you must realize that
!!   use of this routine with 'order' set to values higher than 1 MAY
!!   cause asymmetric interpolation and your results may loose symmetry.
!!
!!   Data is passed in in the array 'recv' and returned in the array
!!   'temp'.  The order of the interpolating polynomial is also passed
!!   in the variable 'order' and can take on value ranging from 1 to 5.
!!   The last argument 'ivar' specifies which variable in 'unk_e_?' to apply
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

      Subroutine amr_restrict_ec_genorder(recv,temp,icoord,order,ivar)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Input/Output arguments.
      Real,    Intent(in)    :: recv(:,:,:,:)
      Real,    Intent(inout) :: temp(:,:,:,:)
      Integer, Intent(in)    :: icoord, order, ivar

!-----Logcal variables and arrays.
      Real    :: xi, xj, www
      Real, Save :: weight(5,3,-4:5)
      Integer :: i,j,k
      Integer :: iparmin,jparmin,kparmin
      Integer :: iparmax,jparmax,kparmax
      Integer :: iii,jjj,kkk
      Integer :: is,js,ks,iw,jw,kw
      Integer :: i0, j0, k0
      Integer :: istart, jstart, kstart
      Integer :: iend, jend, kend
      Integer :: order2
      Logical, Save :: first = .true.

!-----Begin executable code.

      If (first) Then

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
      xi = Real(istart)-.5
      Do i = istart,iend
         weight(order2,2,i) = 1.
         xj = Real(istart)-.5
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
      xi = Real(istart)-.5
      Do i = istart,iend
         weight(order2,3,i) = 1.
         xj = Real(istart)-.5
         Do j = istart,iend
            If (i .ne. j) Then
               weight(order2,3,i) =                                    & 
                    weight(order2,3,i)*(0.-xj)/(xi-xj)
            End If
            xj = xj + 1.
         End Do
         xi = xi + 1.
      End Do

      End Do  ! End Do order2 = 1, 5

      End If  ! End If (first)

      If (icoord == 1) Then                         ! x-edge variables

        iparmin = 1+nguard
        iparmax = nxb+nguard
        jparmin = 1+nguard*k2d
        jparmax = nyb+(nguard+1)*k2d
        kparmin = 1+nguard*k3d
        kparmax = nzb+(nguard+1)*k3d


        Do k0 = kparmin,kparmax,2
        Do j0 = jparmin,jparmax,2
        Do i0 = iparmin,iparmax,2


           If (i0 == iparmin) Then
              istart = 0
              iw = 1
           ElseIf (i0 == iparmax-1) Then
              istart = -order + 1
              iw = 3
           Else
              istart = -int(order/2)
              iw = 2
           End If
           is = i0+istart
           iend = istart + order
           k = k0
           j = j0

           temp(ivar,i0,j0,k0) = 0.

           i = is
           Do iii = istart,iend

           www = weight(order,iw,iii)

           If (curvilinear_conserve) Then
           temp(ivar,i0,j0,k0) =                                       & 
              temp(ivar,i0,j0,k0) + recv(ivar,i,j,k)
           Else
           temp(ivar,i0,j0,k0) =                                       & 
              temp(ivar,i0,j0,k0) + (www*recv(ivar,i,j,k))
           End If

           i = i + 1
           End Do  ! End DO iii = istart,iend

        End Do  ! End Do i0 = iparmin,iparmax,2
        End Do  ! End Do j0 = jparmin,jparmax,2
        End Do  ! End Do k0 = kparmin,kparmax,2

      ElseIf (icoord == 2) Then                     ! y-edge variables

        iparmin = 1+nguard
        iparmax = nxb+nguard+1
        jparmin = 1+nguard*k2d
        jparmax = nyb+nguard*k2d
        kparmin = 1+nguard*k3d
        kparmax = nzb+(nguard+1)*k3d

        Do k0 = kparmin,kparmax,2
        Do j0 = jparmin,jparmax,2
        Do i0 = iparmin,iparmax,2

           If (ndim >= 2) Then
              If (j0 == jparmin) Then
                 jstart = 0
                 jw = 1
              ElseIf (j0 == jparmax-1) Then
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
           
           k = k0
           i = i0

           temp(ivar,i0,j0,k0) = 0.

           j = js
           Do jjj = jstart,jend

           www = weight(order,jw,jjj)

           If (curvilinear_conserve) Then
           temp(ivar,i0,j0,k0) =                                       &
              temp(ivar,i0,j0,k0) + recv(ivar,i,j,k)
           Else
           temp(ivar,i0,j0,k0) =                                       & 
              temp(ivar,i0,j0,k0) + (www*recv(ivar,i,j,k))

           End If

           j = j + 1
           End Do  ! End Do jjj = jstart,jend

        End Do  ! End Do i0 = iparmin,iparmax,2
        End Do  ! End Do j0 = jparmin,jparmax,2
        End Do  ! End Do k0 = kparmin,kparmax,2

      ElseIf (icoord == 3) Then                     ! z-edge variables

        iparmin = 1+nguard
        iparmax = nxb+nguard+1
        jparmin = 1+nguard*k2d
        jparmax = nyb+(nguard+1)*k2d
        kparmin = 1+nguard*k3d
        kparmax = nzb+nguard*k3d

        Do k0 = kparmin,kparmax,2
        Do j0 = jparmin,jparmax,2
        Do i0 = iparmin,iparmax,2

           If (ndim == 3) Then
              If (k0 == kparmin) Then
                 kstart = 0
                 kw = 1
              ElseIf (k0 == kparmax-1) Then
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
           
           i = i0
           j = j0

           temp(ivar,i0,j0,k0) = 0.

           k = ks
           Do kkk = kstart,kend

           www = weight(order,kw,kkk)

           If (curvilinear_conserve) Then
           temp(ivar,i0,j0,k0) =                                       & 
              temp(ivar,i0,j0,k0) + recv(ivar,i,j,k)
           Else
           temp(ivar,i0,j0,k0) =                                       & 
              temp(ivar,i0,j0,k0) + (www*recv(ivar,i,j,k))
           End If

           k = k + 1
           End Do  ! End Do kkk = kstart,kend

        End Do  ! End Do i0 = iparmin,iparmax,2
        End Do  ! End Do j0 = jparmin,jparmax,2
        End Do  ! End Do k0 = kparmin,kparmax,2

      End If

      Return
      End Subroutine amr_restrict_ec_genorder
