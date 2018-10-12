!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_prolong_cc_fun_init
!! NAME
!!
!!   amr_prolong_cc_fun_init
!!
!! SYNOPSIS
!!
!!   Call amr_prolong_cc_fun_init()
!!   
!! ARGUMENTS
!!
!!   No arguments.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   workspace
!!   prolong_arrays
!!
!! CALLS
!!
!!   No other Paramesh routines called.
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine computes the values of dx,dy and dz used during the
!!   interpolation process. These are used inside the prolongation routines
!!   saving needless repetitive computation at the cost of minimal storage
!!   space.
!!
!!   This particular prolongation is simple linear interpolation. It can
!!   be used for blocks with an even or odd number of grid cells.
!!   If CONSERVE is defined then the new mesh points immediately adjacent
!!   to the old block boundaries are treated specially, in a way which
!!   guarantees conservation.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          June 1997
!!
!!**

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_prolong_cc_fun_init

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use prolong_arrays

      Implicit None

!-----Local variables and arrays
      Real :: xc0,yc0,zc0,xc,yc,zc
      Integer :: i,j,k,ii,jj,kk,loop,jchild,ioff,joff,koff
      Integer :: i1,j1,k1,i1p,j1p,k1p

!-----Begin executable code.

!-----Conditional constants. Vary depending on whether a block has an even or
!-----odd number of cells along a given axis.
      xc0 = .75
      yc0 = .75
      zc0 = .75
      If (Mod(nxb,2) == 1) xc0=.5
      If (Mod(nyb,2) == 1) yc0=.5
      If (Mod(nzb,2) == 1) zc0=.5

!-----Initialize the values of dx,dy,dz needed in prolong_unk_fun.
      Do k=kl_bnd1,ku_bnd1
        kk = k-nguard+2*nzb
        zc = zc0+real(kk)*.5
        prol_dz(k) = Mod(zc,1.)
      End Do
      Do j=jl_bnd1,ju_bnd1
        jj = j-nguard+2*nyb
        yc = yc0+real(jj)*.5
        prol_dy(j) = Mod(yc,1.)
      End Do
      Do i=il_bnd1,iu_bnd1
        ii = i-nguard+2*nxb
        xc = xc0+real(ii)*.5
        prol_dx(i) = Mod(xc,1.)
      End Do

!-----Compute the indeces used in the interpolation
!-----This includes the special conservative treatment near the block boundary.
!-----The outer loop selects which offset value is being used for indexing (ie
!-----which side of the parent block the child will be on.)

      Do loop = 1,2

! compute the offset in the parent block appropriate for the different children
      If (loop == 1) jchild=1
      If (loop == 2) jchild=nchild
      ioff = Mod(jchild-1,2)*nxb/2
      joff = Mod((jchild-1)/2,2)*nyb/2
      koff = Mod((jchild-1)/4,2)*nzb/2

                                                       ! note the 2*nxb and
      Do i=il_bnd1,iu_bnd1                             ! nxb components in the
        ii = i-nguard                                  ! expression for i1 are
        i1 = (ii+nxb*2)/2-nxb+ioff+nguard              ! included so i1 will
        i1p= i1+1                                      ! be correct for -ve
        prol_indexx(1,i,loop) = i1                     ! values of i also.
        prol_indexx(2,i,loop) = i1p                    ! (true also for j1,k1)
        If (conserve) Then
        If (Mod(nxb,2) == 0) Then
          If (ioff == 0) Then
            If (i == nguard) prol_indexx(2,i,loop) = i1
            If (i == nguard+1) prol_indexx(1,i,loop) = i1p
          Else
            If (i == iu_bnd1-nguard) prol_indexx(2,i,loop) = i1
            If (i == iu_bnd1-nguard+1) prol_indexx(1,i,loop) = i1p
          End If
        End If
        End If  ! End If (conserve)
      End Do  ! End Do i=il_bnd1,iu_bnd1

      prol_indexy(:,:,loop)=jl_bnd1
      If (ndim >= 2) Then
      Do j=jl_bnd1,ju_bnd1
        jj = j-nguard
        j1 = (jj+nyb*2)/2-nyb+joff+nguard
        j1p= j1+1
        prol_indexy(1,j,loop) = j1
        prol_indexy(2,j,loop) = j1p
        If (conserve) Then
        If (Mod(nyb,2) == 0) Then
          If (joff == 0) Then
            If (j == nguard) prol_indexy(2,j,loop) = j1
            If (j == nguard+1) prol_indexy(1,j,loop) = j1p
          Else
            If (j == ju_bnd1-nguard) prol_indexy(2,j,loop) = j1
            If (j == ju_bnd1-nguard+1) prol_indexy(1,j,loop) = j1p
          End If
        End If
        End If  ! End If (conserve)
      End Do  ! End Do j=jl_bnd1,ju_bnd1
      End If  ! End If (ndim >= 2) 

      prol_indexz(:,:,loop)=kl_bnd1
      If (ndim == 3) Then
      Do k=kl_bnd1,ku_bnd1
        kk = k-nguard
        k1 = (kk+nzb*2)/2-nzb+koff+nguard
        k1p= k1+1
        prol_indexz(1,k,loop) = k1
        prol_indexz(2,k,loop) = k1p
        If (conserve) Then
        If (Mod(nzb,2) == 0) Then
          If (koff == 0) Then
            If (k == nguard) prol_indexz(2,k,loop) = k1
            If (k == nguard+1) prol_indexz(1,k,loop) = k1p
          Else
            If (k == ku_bnd1-nguard) prol_indexz(2,k,loop) = k1
            If (k == ku_bnd1-nguard+1) prol_indexz(1,k,loop) = k1p
          End If
        End If
        End If  ! End If (conserve)
      End Do  ! End Do k=kl_bnd1,ku_bnd1
      End If  ! End If (ndim == 3)

      End Do  ! End Do loop = 1,2

!-----set flag to pass error check at the start of prolong_unk_fun
      prol_init = 100

!-------Initialize the values of dx,dy,dz needed in prolong_work_fun
      Do k=klw1,kuw1
        kk = k-nguard_work+2*nzb
        zc = zc0+real(kk)*.5
        prolw_dz(k) = Mod(zc,1.)
        Do j=jlw1,juw1
           jj = j-nguard_work+2*nyb
           yc = yc0+real(jj)*.5
           prolw_dy(j) = Mod(yc,1.)
           Do i=ilw1,iuw1
               ii = i-nguard_work+2*nxb
               xc = xc0+real(ii)*.5
               prolw_dx(i) = Mod(xc,1.)
           End Do
        End Do
      End Do

!-------Compute the indeces used in the interpolation
!-------This includes the special conservative treatment near the block boundary.
!-------The outer loop selects which offset value is being used for indexing (ie
!-------which side of the parent block the child will be on.)

      Do loop = 1,2

!-------compute the offset in the parent block appropriate for the different children
        If (loop == 1) jchild=1
        If (loop == 2) jchild=nchild
        ioff = Mod(jchild-1,2)*nxb/2
        joff = Mod((jchild-1)/2,2)*nyb/2
        koff = Mod((jchild-1)/4,2)*nzb/2

        Do i=ilw1,iuw1
          ii = i-nguard_work
          i1 = (ii+nxb*2)/2-nxb+ioff+nguard_work
          i1p= i1+1
          prolw_indexx(1,i,loop) = i1
          prolw_indexx(2,i,loop) = i1p
          If (conserve) Then
          If (Mod(nxb,2) == 0) Then
            If (ioff == 0) Then
              If (i == nguard_work) prolw_indexx(2,i,loop) = i1
              If (i == nguard_work+1) prolw_indexx(1,i,loop) = i1p
            Else
              If (i == iu_bnd1-nguard_work) prolw_indexx(2,i,loop) = i1
              If (i == iu_bnd1-nguard_work+1)                          &
                       prolw_indexx(1,i,loop) = i1p
            End If
          End If
          End If  ! End If (conserve)
        End Do  ! End Do i=ilw1,iuw1

        prolw_indexy(:,:,loop)=jlw1
        If (ndim >= 2) Then
        Do j=jlw1,juw1
          jj = j-nguard_work
          j1 = (jj+nyb*2)/2-nyb+joff+nguard_work
          j1p= j1+1
          prolw_indexy(1,j,loop) = j1
          prolw_indexy(2,j,loop) = j1p
          If (conserve) Then
          If (Mod(nyb,2) == 0) Then
            If (joff == 0) Then
              If (j == nguard_work) prolw_indexy(2,j,loop) = j1
              If (j == nguard_work+1) prolw_indexy(1,j,loop) = j1p
            Else
              If (j == ju_bnd1-nguard_work) prolw_indexy(2,j,loop) = j1
              If (j == ju_bnd1-nguard_work+1)                          &
                       prolw_indexy(1,j,loop) = j1p
            End If
          End If
          End If  ! End If (conserve)
        End Do  ! End Do j=jlw1,juw1
        End If  ! End If (ndim >= 2)

        prolw_indexz(:,:,loop)=klw1
        If (ndim == 3) Then
        Do k=klw1,kuw1
          kk = k-nguard_work
          k1 = (kk+nzb*2)/2-nzb+koff+nguard_work
          k1p= k1+1
          prolw_indexz(1,k,loop) = k1
          prolw_indexz(2,k,loop) = k1p
          If (conserve) Then
          If (Mod(nzb,2) == 0) Then
            If (koff == 0) Then
              If (k == nguard_work) prolw_indexz(2,k,loop) = k1
              If (k == nguard_work+1) prolw_indexz(1,k,loop) = k1p
            Else
              If (k == ku_bnd1-nguard_work) prolw_indexz(2,k,loop) = k1
              If (k == ku_bnd1-nguard_work+1)                          &
                       prolw_indexz(1,k,loop) = k1p
            End If
          End If
          End If  ! End If (conserve)
        End Do  ! End Do k=klw1,kuw1
        End If  ! End If (ndim == 3)

      End Do  ! End Do loop = 1,2

!-----set flag to pass error check at the start of prolong_work_fun
      prolw_init = 100

      Return
      End Subroutine amr_prolong_cc_fun_init
