!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! workspace module


      Module workspace

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Use paramesh_dimensions

      Private

! workspace arrays
      Public :: work, recvw, sendw, tempw, interp_mask_work, & 
     &          interp_mask_work_res
      Real, Allocatable, Save :: work(:,:,:,:,:)
      Real, Allocatable, Save :: recvw(:,:,:)
      Real, Allocatable, Save :: sendw(:,:,:)
      Real, Allocatable, Save :: tempw(:,:,:)
      Integer, Allocatable, Save :: interp_mask_work(:)
      Integer, Allocatable, Save :: interp_mask_work_res(:)

! common block storing the solution for cell-centered quantities.
      Public :: work1, recvw1, tempw1, tempw2
      Real, Allocatable, Save :: work1(:,:,:,:)
      Real, Allocatable, Save :: recvw1(:,:,:,:)
      Real, Allocatable, Save :: tempw1(:,:,:)
      Real, Allocatable, Save :: tempw2(:,:,:)

! arrays used to store geometry information for the working block
      Public :: cell_vol_w
      Real, Allocatable :: cell_vol_w(:,:,:)

! Index arrays used to record destination data values for fine layer
! neighbor guardcells
      Integer,Public :: f2c_ind_work(2,3,27)

      End Module workspace


