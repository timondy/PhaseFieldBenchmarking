!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_t_to_perm
!! NAME
!!
!!   amr_1blk_t_to_perm
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_t_to_perm (lcc,lfc,lec,lnc,lb,idest)
!!   Call amr_1blk_t_to_perm (logical,logical,logical,logical,
!!                             integer,integer)
!!
!! ARGUMENTS
!!
!!   Logical, Intent(in) :: lcc  copies cell centered data if true
!!   Logical, Intent(in) :: lfc  copies cell face-centered data if true
!!   Logical, Intent(in) :: lec  copies cell edge-centered data if true
!!   Logical, Intent(in) :: lnc  copies cell corner data if true
!!   Integer, Intent(in) :: lb   block into which data is to be copied
!!   Integer, Intent(in( :: idest sets value for last dimension index
!!                                 in the 1-blk data arrays
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
!!   This routine copies data from the 1-block working arrays with guardcells
!!   to the permanent data arrays, which may or may not have permanent
!!   guardcells, depending on whether NO_PERMANENT_GUARDCELLS is set to true. 
!!
!! AUTHORS
!!
!! Written :     Peter MacNeice          February 1999
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_t_to_perm (lcc,lfc,lec,lnc,lb,idest)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace

      Implicit None

!-----Input/Output Arguments
      integer, intent(in) :: lb,idest
      logical, intent(in) :: lcc,lfc,lec,lnc

!-----Local arrays and variables
      integer :: nguard0

!-----Begin Executable Code

      If (var_dt .or. pred_corr) Then

      nguard0 = nguard*(1-npgs)

!------cell-centered data
       If (lcc) Then

           If (no_permanent_guardcells) Then
           t_unk(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)       & 
             = unk1(:,il_bnd+nguard0:iu_bnd+nguard0,                   & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
           Else
           t_unk(:,:,:,:,lb) = unk1(:,:,:,:,idest)
           End If

       End If  ! End If (lcc)

!------cell face-centered data
       If (lfc) Then
!--------x-face
         If (no_permanent_guardcells) Then
         tfacevarx(1:nfacevar,il_bnd:iu_bnd+1,                         & 
                             jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)           & 
             = facevarx1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0+1,   & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         Else
         tfacevarx(1:nfacevar,:,:,:,lb) =                              & 
                              facevarx1(1:nfacevar,:,:,:,idest)
         End If

         If (ndim > 1) Then
!--------y-face
         If (no_permanent_guardcells) Then
         tfacevary(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,         & 
                               kl_bnd:ku_bnd,lb)                       & 
             = facevary1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0,     & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d+k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         Else
         tfacevary(1:nfacevar,:,:,:,lb) =                              & 
                              facevary1(1:nfacevar,:,:,:,idest)
         End If
         End If  ! End If (ndim > 1)

         If (ndim == 3) Then
!--------z-face
         If (no_permanent_guardcells) Then
         tfacevarz(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd,             & 
                               kl_bnd:ku_bnd+k3d,lb)                   & 
             = facevarz1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0,     & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d+k3d,idest)
         Else
         tfacevarz(1:nfacevar,:,:,:,lb) =                              & 
                              facevarz1(1:nfacevar,:,:,:,idest)
         End If
         End If  ! End If (ndim == 3)

        End If  ! End If (lfc)


!------cell edge-centered data
       If (lec) Then
         If (ndim > 1) Then
!--------x-edge
         If (no_permanent_guardcells) Then
         t_unk_e_x(1:nvaredge,il_bnd:iu_bnd,                           & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)           & 
             = unk_e_x1(1:nvaredge,il_bnd+nguard0:iu_bnd+nguard0,      & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         Else
         t_unk_e_x(1:nvaredge,:,:,:,lb) =                              & 
                              unk_e_x1(1:nvaredge,:,:,:,idest)
         End If
!--------y-edge
         If (no_permanent_guardcells) Then
         t_unk_e_y(1:nvaredge,il_bnd:iu_bnd+1,                         & 
                     jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)               & 
             = unk_e_y1(1:nvaredge,il_bnd+nguard0:iu_bnd+1+nguard0,    & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         Else
         t_unk_e_y(1:nvaredge,:,:,:,lb) =                              & 
                              unk_e_y1(1:nvaredge,:,:,:,idest)
         End If

         If (ndim == 3) Then
!--------z-edge
         If (no_permanent_guardcells) Then
         t_unk_e_z(1:nvaredge,il_bnd:iu_bnd+1,                         & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb)               & 
             = unk_e_z1(1:nvaredge,il_bnd+nguard0:iu_bnd+1+nguard0,    & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         Else
         t_unk_e_z(1:nvaredge,:,:,:,lb) =                              & 
                              unk_e_z1(1:nvaredge,:,:,:,idest)
         End If
         End If  ! End If (ndim == 3)

         End If  ! End If (ndim > 1)
        End If  ! End If (lec)

!------cell corner data
       If (lnc) Then
         If (no_permanent_guardcells) Then
         t_unk_n(1:nvarcorn,il_bnd:iu_bnd+1,                           & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)           & 
             = unk_n1(1:nvarcorn,il_bnd+nguard0:iu_bnd+1+nguard0,      & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         Else
         t_unk_n(1:nvarcorn,:,:,:,lb) =                                & 
                              unk_n1(1:nvarcorn,:,:,:,idest)
        End If
       End If  ! End If (lnc)

      End If  ! End If (var_dt .or. pred_corr)

      Return
      End Subroutine amr_1blk_t_to_perm
