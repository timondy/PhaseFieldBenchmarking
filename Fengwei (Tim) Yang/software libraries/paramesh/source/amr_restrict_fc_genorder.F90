!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_fc_genorder.F90
!! NAME
!!   amr_restrict_fc_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_fc_genorder(recv,temp,icoord,order,ivar)
!!   Call amr_restrict_fc_genorder(real array,real array,integer,integer,integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: recv(:,:,:,:)  Face-centered data to restrict             
!!   Real,    Intent(inout) :: temp(:,:,:,:)  Restricted face-centered data              
!!   Integer, Intent(in)    :: icoord         selects which edge to operate on           
!!                                            iccord = 1, selects x-face                 
!!                                            iccord = 2, selects y-face                 
!!                                            iccord = 3, selects z-face                 
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
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   edge centered data stored in 'facevar?'.  It uses a lagrange polynomial
!!   interpolation scheme.  Also, the interpolation stencil is automatically
!!   shifted to avoid using data in guardcells. CAUTION: you must realize that
!!   use of this routine with 'order' set to values higher than 1 MAY
!!   cause asymmetric interpolation and your results may loose symmetry.
!!
!!   Data is passed in in the array 'recv' and returned in the array
!!   'temp'.  The order of the interpolating polynomial is also passed
!!   in the variable 'order' and can take on value ranging from 1 to 5.
!!   The last argument 'ivar' specifies which variable in 'facevar?' to apply
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

      Subroutine amr_restrict_fc_genorder(recv,temp,icoord,order,ivar)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Input/Output arguments
      Real,    Intent(in)    :: recv(:,:,:,:)
      Real,    Intent(inout) :: temp(:,:,:,:)
      Integer, Intent(in)    :: icoord, order, ivar

!-----Local variables and arrays
      Real,Save :: weight(5,3,-4:5)
      Real      :: xi, xj, www
      Integer :: i,j,k
      Integer :: iparmin,jparmin,kparmin
      Integer :: iparmax,jparmax,kparmax
      Integer :: iii,jjj,kkk
      Integer :: is,js,ks,iw,jw,kw
      Integer :: i0, j0, k0
      Integer :: istart, jstart, kstart
      Integer :: iend, jend, kend
      Integer :: order2
      Logical,Save :: first = .True.

!-----Begin executable code

      If (first) Then

      first = .False.

      Do order2 = 1, 5

!-----left

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

      End Do  ! End Do order2 = 1, 5

      End If  ! End If (first)

      If (icoord == 1) Then                         ! x-face variables

        iparmin = 1+nguard
        iparmax = nxb+nguard+1
        jparmin = 1+nguard*k2d
        jparmax = nyb+nguard*k2d
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

           i = i0
        
           temp(ivar,i0,j0,k0) = 0.
           
           k = ks
           Do kkk = kstart,kend
              j = js
              Do jjj = jstart,jend
                 
                 If (ndim == 1) Then
                    www = 1.
                 ElseIf (ndim == 2) Then
                    www = weight(order,jw,jjj)
                 ElseIf (ndim == 3) Then
                    www = weight(order,jw,jjj) * weight(order,kw,kkk)
                 End If
                 
                 If (curvilinear_conserve) Then
                    temp(ivar,i0,j0,k0) =                              & 
                         temp(ivar,i0,j0,k0) + (recv(ivar,i,j,k))
                 Else
                    temp(ivar,i0,j0,k0) =                              &  
                         temp(ivar,i0,j0,k0) + (www*recv(ivar,i,j,k))
                 End If
                 
                 j = j + 1
              End Do  ! End Do jjj = jstart,jend
              k = k + 1
           End Do  ! End Do kkk = kstart,kend
           
        End Do  ! End Do i0 = iparmin,iparmax,2
        End Do  ! End Do j0 = jparmin,jparmax,2
        End Do  ! End Do k0 = kparmin,kparmax,2

      ElseIf (icoord == 2) Then                     ! y-face variables
           
        iparmin = 1+nguard
        iparmax = nxb+nguard
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

           If (i0 == iparmin) Then
              istart = 0
              iw = 1
           ElseIf (i0 == iparmax-1) Then
              istart = -order + 1
              iw = 3
           Else
              istart = -int(order/2)
              iw = 2
           End If  ! End If i0 = iparmin,iparmax,2

           is = i0+istart
           iend = istart + order
           
           j = j0

           temp(ivar,i0,j0,k0) = 0.

           k = ks
           Do kkk = kstart,kend
              i = is
              Do iii = istart,iend
                 
                 If (ndim == 1) Then
                    www = 1.
                 ElseIf (ndim == 2) Then
                    www = weight(order,iw,iii)
                 ElseIf (ndim == 3) Then
                    www = weight(order,iw,iii) * weight(order,kw,kkk)
                 End If
              
                 If (curvilinear_conserve) Then
                    temp(ivar,i0,j0,k0) =                              & 
                         temp(ivar,i0,j0,k0) + (recv(ivar,i,j,k))
                 Else
                    temp(ivar,i0,j0,k0) =                              & 
                         temp(ivar,i0,j0,k0) + (www*recv(ivar,i,j,k))
                 End If  
              
                 i = i + 1
              End Do  ! End Do iii = istart,iend
              k = k + 1
           End Do  ! End Do kkk = kstart,kend

        End Do  ! End Do i0 = iparmin,iparmax,2
        End Do  ! End Do j0 = jparmin,jparmax,2
        End Do  ! End Do k0 = kparmin,kparmax,2

      ElseIf (icoord == 3) Then                     ! z-face variables

        iparmin = 1+nguard
        iparmax = nxb+nguard
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
           
           If (i0 == iparmin) Then
              istart = 0
              iw = 1
           ElseIf (i0 == iparmax-1) Then
              istart = -order + 1
              iw = 3
           Else
              istart = -int(order/2)
              iw = 2
           End If  ! End If (i0 == iparmin)

           is = i0+istart
           iend = istart + order

           k = k0

           temp(ivar,i0,j0,k0) = 0.
           
           j = js
           Do jjj = jstart,jend
              i = is
              Do iii = istart,iend
                 
                 www = weight(order,iw,iii) * weight(order,jw,jjj)
                 
                 If (curvilinear_conserve) Then
                    temp(ivar,i0,j0,k0) =                              & 
                         temp(ivar,i0,j0,k0) + (recv(ivar,i,j,k))
                 Else
                    temp(ivar,i0,j0,k0) =                              & 
                         temp(ivar,i0,j0,k0) + (www*recv(ivar,i,j,k))
                 End If
              
                 i = i + 1
              End Do  ! End Do iii = istart,iend
              j = j + 1
           End Do  ! End Do jjj = jstart,jend
        
        End Do  ! End Do i0 = iparmin,iparmax,2
        End Do  ! End Do j0 = jparmin,jparmax,2
        End Do  ! End Do k0 = kparmin,kparmax,2

      End If  ! End If (icoord == 1)

      End Subroutine amr_restrict_fc_genorder
