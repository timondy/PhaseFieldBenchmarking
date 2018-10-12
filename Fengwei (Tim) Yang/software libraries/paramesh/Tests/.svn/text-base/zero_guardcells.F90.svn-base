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


      subroutine zero_guardcells(ioptw)



! This routine sets all guard cell data values to zero.
!
! Written:      Peter MacNeice August 1997
!
!


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace


      integer ioptw

      integer ka,kb,kc,kd
      integer kbw,kcw,kdw
      integer nguard0, nguard_work0


      if (.not.no_permanent_guardcells) then

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      ka=1
      kb=1+(nguard0-1)*k3d
      kc=1+(nzb+nguard0)*k3d
      kd=nzb+2*nguard0*k3d
      kbw=1+(nguard_work0-1)*k3d
      kcw=1+(nzb+nguard_work0)*k3d
      kdw=nzb+2*nguard_work0*k3d

      if (nvar > 0) then

      unk(:,1:nguard0,:,:,:) = 0.
      unk(:,1+nxb+nguard0:nxb+2*nguard0,:,:,:) = 0.
      unk(:,:,1:nguard0,:,:) = 0.
      unk(:,:,1+nyb+nguard0:nyb+2*nguard0,:,:) = 0.
      if(ndim.eq.3) then
      unk(:,:,:,ka:kb,:) = 0.
      unk(:,:,:,kc:kd,:) = 0.
      endif

      end if

      if (nvarcorn > 0) then

      unk_n(:,1:nguard0,:,:,:) = 0.
      unk_n(:,2+nxb+nguard0:nxb+2*nguard0+1,:,:,:) = 0.
      unk_n(:,:,1:nguard0,:,:) = 0.
      unk_n(:,:,2+nyb+nguard0:nyb+2*nguard0+1,:,:) = 0.
      if(ndim.eq.3) then
      unk_n(:,:,:,ka:kb,:) = 0.
      unk_n(:,:,:,kc+k3d:kd+k3d,:) = 0.
      endif

      end if

      if (nvaredge > 0) then

      unk_e_x(:,1:nguard0,:,:,:) = 0.
      unk_e_x(:,1+nxb+nguard0:nxb+2*nguard0,:,:,:) = 0.
      unk_e_x(:,:,1:nguard0,:,:) = 0.
      unk_e_x(:,:,2+nyb+nguard0:nyb+2*nguard0+1,:,:) = 0.
      if(ndim.eq.3) then
      unk_e_x(:,:,:,ka:kb,:) = 0.
      unk_e_x(:,:,:,kc+k3d:kd+k3d,:) = 0.
      endif

      unk_e_y(:,1:nguard0,:,:,:) = 0.
      unk_e_y(:,2+nxb+nguard0:nxb+2*nguard0+1,:,:,:) = 0.
      unk_e_y(:,:,1:nguard0,:,:) = 0.
      unk_e_y(:,:,1+nyb+nguard0:nyb+2*nguard0,:,:) = 0.
      if(ndim.eq.3) then
      unk_e_y(:,:,:,ka:kb,:) = 0.
      unk_e_y(:,:,:,kc+k3d:kd+k3d,:) = 0.
      endif

      unk_e_z(:,1:nguard0,:,:,:) = 0.
      unk_e_z(:,2+nxb+nguard0:nxb+2*nguard0+1,:,:,:) = 0.
      unk_e_z(:,:,1:nguard0,:,:) = 0.
      unk_e_z(:,:,2+nyb+nguard0:nyb+2*nguard0+1,:,:) = 0.
      if(ndim.eq.3) then
      unk_e_z(:,:,:,ka:kb,:) = 0.
      unk_e_z(:,:,:,kc:kd,:) = 0.
      endif

      end if

      if (nvar_work > 0) then

      work(1:nguard_work0,:,:,:,ioptw-1) = 0.
      work(1+nxb+nguard_work0:nxb+2*nguard_work0, & 
     &                                     :,:,:,ioptw-1) = 0.
      work(:,1:nguard_work0,:,:,ioptw-1) = 0.
      work(:,1+nyb+nguard_work0:nyb+2*nguard_work0, & 
     &                                       :,:,ioptw-1) = 0.
      if(ndim.eq.3) then
      work(:,:,1:kbw,:,ioptw-1) = 0.
      work(:,:,kcw:kdw,:,ioptw-1) = 0.
      endif

      end if

      if (nfacevar > 0) then

      facevarx(:,1:nguard0,:,:,:) = 0.
      facevarx(:,2+nxb+nguard0:1+nxb+2*nguard0,:,:,:) = 0.
      facevarx(:,:,1:nguard0,:,:) = 0.
      facevarx(:,:,1+nyb+nguard0:nyb+2*nguard0,:,:) = 0.
      if(ndim.eq.3) then
      facevarx(:,:,:,ka:kb,:) = 0.
      facevarx(:,:,:,kc:kd,:) = 0.
      endif

      facevary(:,1:nguard0,:,:,:) = 0.
      facevary(:,1+nxb+nguard0:nxb+2*nguard0,:,:,:) = 0.
      facevary(:,:,1:nguard0,:,:) = 0.
      facevary(:,:,2+nyb+nguard0:1+nyb+2*nguard0,:,:) = 0.
      if(ndim.eq.3) then
      facevary(:,:,:,ka:kb,:) = 0.
      facevary(:,:,:,kc:kd,:) = 0.
      endif

      facevarz(:,1:nguard0,:,:,:) = 0.
      facevarz(:,1+nxb+nguard0:nxb+2*nguard0,:,:,:) = 0.
      facevarz(:,:,1:nguard0,:,:) = 0.
      facevarz(:,:,1+nyb+nguard0:nyb+2*nguard0,:,:) = 0.
      if(ndim.eq.3) then
      facevarz(:,:,:,ka:kb,:) = 0.
      facevarz(:,:,:,1+(1+nzb+nguard0)*k3d:1+(nzb+2*nguard0)*k3d,:)=0.
      endif

      end if

      endif

      return
      end
