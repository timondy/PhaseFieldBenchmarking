!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
!-----------------------------------------------------------------
! prolong_arrays Module



      Module prolong_arrays

      Use paramesh_dimensions

      Private


      Public :: prol_dx,prol_dy,prol_dz
      Public :: prol_indexx,prol_indexy,prol_indexz,prol_init
      Public :: prol_f_dx,prol_f_dy,prol_f_dz
      Public :: prol_f_indexx,prol_f_indexy,prol_f_indexz
      Public :: prol_f_init
      Public :: prolw_dx,prolw_dy,prolw_dz
      Public :: prolw_indexx,prolw_indexy,prolw_indexz,prolw_init
      Real, Save, Allocatable :: prol_dx(:)
      Real, Save, Allocatable :: prol_dy(:) 
      Real, Save, Allocatable :: prol_dz(:)
      Integer, Save, Allocatable :: prol_indexx(:,:,:)
      Integer, Save, Allocatable :: prol_indexy(:,:,:)
      Integer, Save, Allocatable :: prol_indexz(:,:,:)
      Real, Save, Allocatable :: prol_f_dx(:)
      Real, Save, Allocatable :: prol_f_dy(:)
      Real, Save, Allocatable :: prol_f_dz(:)
      Integer, Save, Allocatable :: prol_f_indexx(:,:,:)
      Integer, Save, Allocatable :: prol_f_indexy(:,:,:)
      Integer, Save, Allocatable :: prol_f_indexz(:,:,:)
      Real, Save, Allocatable :: prolw_dx(:)
      Real, Save, Allocatable :: prolw_dy(:)
      Real, Save, Allocatable :: prolw_dz(:)
      Integer, Save, Allocatable :: prolw_indexx(:,:,:)
      Integer, Save, Allocatable :: prolw_indexy(:,:,:)
      Integer, Save, Allocatable :: prolw_indexz(:,:,:)
      Integer, Save :: prol_init
      Integer, Save :: prol_f_init
      Integer, Save :: prolw_init

      Public :: prol_fc_dbz, prol_fc_dbz_ivar, prol_fc_dbz_n
      Logical, Save :: prol_fc_dbz = .false.
      Integer, Allocatable, Save :: prol_fc_dbz_ivar(:,:)
      Integer, Save :: prol_fc_dbz_n = 0

      Public :: prol_fc_clean_divb, prol_fc_clean_divb_ivar,  & 
     &     prol_fc_clean_divb_n
      Logical, Save :: prol_fc_clean_divb = .false.
      Integer, Allocatable, Save :: prol_fc_clean_divb_ivar(:,:)
      Integer, Save :: prol_fc_clean_divb_n = 0

      end Module prolong_arrays
