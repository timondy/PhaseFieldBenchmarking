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


      subroutine amr_1blk_fc_prol_user( & 
     &       recv,ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff, & 
     &       mype,lb,pe_p,lb_p,iface,ivar)


!
! Arguments :
!    interp_mask      integer mask vector     this marks the cell-face-centered
!                                             variables which are to be
!                                             prolonged conservatively
!

!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted 
! from one of the arrays facevarx(y)(z), and performs a prolongation
! operation on it. The data in recv is from a parent block and the
! result of the prolongation operation is written into layer idest of
! the working block arrays facevarx(y)(z)1.
! The position of the child block, isg, within the 
! parent block is specified by the ioff, joff and koff arguments.
!
! iface controls which array is updated, ie facevarx1 if iface=1,
! facevary1 if iface=2, and facevarz1 if iface=3.
!
! This particular prolongation implements a Muscl style monotonization
! which guarantees a conservative second order interpolation. It can
! only be used for blocks with an even number of grid cells.

! This routine implements linear interpolation in the direction perpendicular
! to the face in question, ie for facevarx1 this is the x direction, and
! a Muscl style monotonic interpolant which guarantees a conservative second
! order interpolation in the other directions.
!


! It will only work if nguard > 1.


! The 1D definition of the Muscl algorithm is as follows:
!
! if child grid cell ic1 and ic2 are the two child cells corresponding
! to cell ip of the parent, then
!
!                 U(ic1) = U(ip) - gradm * .5 * dxc
!                 U(ic2) = U(ip) + gradm * .5 * dxc
!
! where dxc is the cell size of the children, and
!
!        gradm = s * max( 0. , min( abs(gradc) , 2*s*gradl, 2*s*gradr ) )
!
! with
!        gradc = (U(ip+1) - U(ip-1) )/ (2.*dxp)
!        gradl = (U(ip)   - U(ip-1) )/ dxp
!        gradr = (U(ip+1) - U(ip)   )/ dxp
!
! where dxp (=2*dxc) is the cell size of the children.
! The multidimensional implementation applies this 1D operation successively
! in each direction as required.
!
!
! This particular prolongation can only be used for blocks with an even 
! number of grid cells.
!
! Note: before using this routine in your program, make sure that the
! routine amr_prolong_face_fun_init has been called. This will normally
! be called from amr_initialize.
!
!
! Written :     Peter MacNeice          September 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use prolong_arrays

      implicit none

      include 'mpif.h'

!------------------------------------
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: lb,lb_p,pe_p
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)

!------------------------------------
! local arrays and variables

      real :: recv1(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &          kl_bnd1:ku_bnd1+k3d)
      real :: temp(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &          kl_bnd1:ku_bnd1+k3d)
      real :: gradl,gradr,gradc,gradm,ss
      real :: dxpr, dypr, dzpr, dxx, dyy, dzz
      real :: cxx, cyy, czz, sdx, sdy, sdz
      real :: factor

      integer, parameter :: largei=100
      integer, parameter :: i_muscl = 3
      integer :: ii, jj, kk, parent_blk, parent_pe
      integer :: icl, icu, jcl, jcu, kcl, kcu
      integer :: i_ind, j_ind, k_ind
      integer :: ilow, ihi, jlow, jhi, klow, khi
      integer :: i, j, k, i1, j1, k1
      integer :: i1p, j1p, k1p, is
      integer :: ierr, ierrorcode

!------------------------------------

      if(prol_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : amr_1blk_fc_prol_muscl. ', & 
     &       'You must call amr_prolong_face_fun_init ', & 
     &       'before you can use this routine!'
       call amr_abort
      endif

        if(nguard.le.1) then
          write(*,*) ' Error - the muscl interpolation version of', & 
     &               ' the routine amr_1blk_fc_prol_gen_unk_fun ', & 
     &               ' requires nguard > 1.'
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif

!------------------------------------


      parent_blk= lb_p
      parent_pe = mype
      if (curvilinear) then
         call amr_block_geometry(parent_blk,parent_pe)
      endif


! Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb


        i_ind = 1
        j_ind = 1
        k_ind = 1
        if(ioff.gt.0) i_ind = 2
        if(joff.gt.0) j_ind = 2
        if(koff.gt.0) k_ind = 2


! reciprocal of parent cell size
        dxpr = .5
        dypr = .5
        dzpr = .5


        ilow = (icl-nguard-1+largei)/2+nguard+ioff-largei/2
        ihi  = (icu-nguard-1+largei)/2+nguard+ioff-largei/2+2
        jlow = ((jcl-nguard-1+largei)/2+nguard+joff-largei/2-1)*k2d+1
        jhi  = ((jcu-nguard-1+largei)/2+nguard+joff-largei/2+1)*k2d+1
        klow = ((kcl-nguard-1+largei)/2+nguard+koff-largei/2-1)*k3d+1
        khi  = ((kcu-nguard-1+largei)/2+nguard+koff-largei/2+1)*k3d+1


! Interpolation loop.
!
! Note that the range of indeces used in the facevar plane differs
! depending on the value of iface_off. This assumes that the face values
! corresponding to index 0 (ie nguard faces to the left of the block
! boundary) are never needed, when iface_off=-1.

        if(iface.eq.1) then

      if (curvilinear) then
! Multiply parent solution by the area of the parent cells y face
      do k=kl_bnd1,ku_bnd1
      do j=jl_bnd1,ju_bnd1
      do i=il_bnd1,iu_bnd1+1
        recv(ivar,i,j,k) = recv(ivar,i,j,k)*cell_area1(i,j,k)
      enddo
      enddo
      enddo
      endif
 
!
! Perform sweep in x direction
        do k=klow,khi
        do j=jlow,jhi
           do i=icl,icu+iface_off
              i1  = prol_f_indexx(1,i,i_ind)
              i1p = prol_f_indexx(2,i,i_ind)
              dxx = prol_f_dx(i)
              cxx = 1.-dxx

! compute interpolated values at location (i,j,k)
                 temp(ivar,i,j,k) = & 
     &                          dxx*recv(ivar,i1,j,k) + & 
     &                          cxx*recv(ivar,i1p,j,k)
           enddo
        enddo
        enddo

!
! Perform sweep in y direction
        if( ndim.ge.2) then

        do k=klow,khi
        do j=jlow,jhi
        do i=icl,icu+iface_off
        recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo

        do k=klow,khi
        do i=icl,icu+iface_off
        do j=jcl,jcu
           jj = (j-nguard-1+largei)/2+nguard+1+joff-largei/2
              gradl = (recv1(ivar,i,jj,k)-recv1(ivar,i,jj-k2d,k)) & 
     &                                         *dypr
              gradr = (recv1(ivar,i,jj+k2d,k)-recv1(ivar,i,jj,k)) & 
     &                                         *dypr
              gradc = (recv1(ivar,i,jj+k2d,k)-recv1(ivar,i,jj-k2d,k)) & 
     &                                         *dypr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss, & 
     &                                       2.*gradl*ss))
              is = mod(j-nguard-1+largei,2)
              sdy = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,i,jj,k)+gradm*sdy
        enddo
        enddo
        enddo

        endif

!
! Perform sweep in z direction
        if(ndim.eq.3) then

        do k=klow,khi
        do j=jcl,jcu
        do i=icl,icu+iface_off
        recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo

        do j=jcl,jcu
        do i=icl,icu+iface_off
        do k=kcl,kcu
           kk = (k-nguard-1+largei)/2+nguard+1+koff-largei/2
              gradl = (recv1(ivar,i,j,kk)-recv1(ivar,i,j,kk-k3d)) & 
     &                                         *dzpr
              gradr = (recv1(ivar,i,j,kk+k3d)-recv1(ivar,i,j,kk)) & 
     &                                         *dzpr
              gradc = (recv1(ivar,i,j,kk+k3d)-recv1(ivar,i,j,kk-k3d)) & 
     &                                         *dzpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss, & 
     &                                       2.*gradl*ss))
              is = mod(k-nguard-1+largei,2)
              sdz = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,i,j,kk)+gradm*sdz
        enddo
        enddo
        enddo

        endif


! save new interpolated solution

      if (curvilinear) then
! Generate cell volume data for the current child block
      call amr_block_geometry(lb,mype)

! Divide interpolated solution by the 2**ndim*volume of the child cells
      factor = .5**ndim
      endif

        do k=kcl,kcu
        do j=jcl,jcu
        do i=icl,icu+iface_off
        if (curvilinear) then
           facevarx1(ivar,i,j,k,idest)=temp(ivar,i,j,k) & 
     &                                *factor/cell_area1(i,j,k)
        else
           facevarx1(ivar,i,j,k,idest)=temp(ivar,i,j,k)
        endif

        enddo
        enddo
        enddo


        endif                              ! end of iface=1 



        if(iface.eq.2) then

      if (curvilinear) then
! Multiply parent solution by the area of the parent cells y face
      do k=kl_bnd1,ku_bnd1
      do j=jl_bnd1,ju_bnd1+k2d
      do i=il_bnd1,iu_bnd1
        recv(ivar,i,j,k) = recv(ivar,i,j,k)*cell_area2(i,j,k)
      enddo
      enddo
      enddo
      endif

!
! Perform sweep in y direction
        do k=klow,khi
        do i=ilow,ihi
           do j=jcl,jcu+iface_off
              j1  = prol_f_indexy(1,j,j_ind)
              j1p = prol_f_indexy(2,j,j_ind)
              dyy = prol_f_dy(j)
              cyy = 1.-dyy

! compute interpolated values at location (i,j,k)
                 temp(ivar,i,j,k) = & 
     &                          dyy*recv(ivar,i,j1,k) + & 
     &                          cyy*recv(ivar,i,j1p,k)
           enddo
        enddo
        enddo


        do k=klow,khi
        do j=jcl,jcu+iface_off
        do i=ilow,ihi
        recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo

!
! Perform sweep in x direction
        do k=klow,khi
        do j=jcl,jcu+iface_off
        do i=icl,icu
           ii = (i-nguard-1+largei)/2+nguard+1+ioff-largei/2
              gradl = (recv1(ivar,ii,j,k)-recv1(ivar,ii-1,j,k)) & 
     &                                         *dxpr
              gradr = (recv1(ivar,ii+1,j,k)-recv1(ivar,ii,j,k)) & 
     &                                         *dxpr
              gradc = (recv1(ivar,ii+1,j,k)-recv1(ivar,ii-1,j,k)) & 
     &                                         *dxpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss, & 
     &                                       2.*gradl*ss))
              is = mod(i-nguard-1+largei,2)
              sdx = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,ii,j,k)+gradm*sdx
        enddo
        enddo
        enddo


! Perform sweep in z direction
        if(ndim.eq.3) then

        do k=klow,khi
        do j=jcl,jcu+iface_off
        do i=icl,icu
        recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo

        do j=jcl,jcu+iface_off
        do i=icl,icu
        do k=kcl,kcu
           kk = (k-nguard-1+largei)/2+nguard+1+koff-largei/2
              gradl = (recv1(ivar,i,j,kk)-recv1(ivar,i,j,kk-k3d)) & 
     &                                         *dzpr
              gradr = (recv1(ivar,i,j,kk+k3d)-recv1(ivar,i,j,kk)) & 
     &                                         *dzpr
              gradc = (recv1(ivar,i,j,kk+k3d)-recv1(ivar,i,j,kk-k3d)) & 
     &                                         *dzpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss, & 
     &                                       2.*gradl*ss))
              is = mod(k-nguard-1+largei,2)
              sdz = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,i,j,kk)+gradm*sdz
        enddo
        enddo
        enddo

        endif

! save new interpolated solution

      if (curvilinear) then
! Generate cell volume data for the current child block
      call amr_block_geometry(lb,mype)

! Divide interpolated solution by the 2**ndim*volume of the child cells
      factor = .5**ndim
      endif

        do k=kcl,kcu
        do j=jcl,jcu+iface_off
        do i=icl,icu
        if (curvilinear) then
           facevary1(ivar,i,j,k,idest)=temp(ivar,i,j,k) & 
     &                                *factor/cell_area2(i,j,k)
        else
           facevary1(ivar,i,j,k,idest)=temp(ivar,i,j,k)
        endif
        enddo
        enddo
        enddo


        endif                              ! end of iface=2 



        if(iface.eq.3) then

      if (curvilinear) then
! Multiply parent solution by the area of the parent cells y face
      do k=kl_bnd1,ku_bnd1+k3d
      do j=jl_bnd1,ju_bnd1
      do i=il_bnd1,iu_bnd1
        recv(ivar,i,j,k) = recv(ivar,i,j,k)*cell_area3(i,j,k)
      enddo
      enddo
      enddo
      endif
!
! Perform sweep in z direction
        do j=jlow,jhi
        do i=ilow,ihi
           do k=kcl,kcu+iface_off
              k1  = prol_f_indexz(1,k,k_ind)
              k1p = prol_f_indexz(2,k,k_ind)
              dzz = prol_f_dz(k)
              czz = 1.-dzz

! compute interpolated values at location (i,j,k)
                 temp(ivar,i,j,k) = & 
     &                          dzz*recv(ivar,i,j,k1) + & 
     &                          czz*recv(ivar,i,j,k1p)
           enddo
        enddo
        enddo

!
! Perform sweep in x direction
        do k=kcl,kcu+iface_off
        do j=jlow,jhi
        do i=ilow,ihi
        recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo

        do k=kcl,kcu+iface_off
        do j=jlow,jhi
        do i=icl,icu
           ii = (i-nguard-1+largei)/2+nguard+1+ioff-largei/2
              gradl = (recv1(ivar,ii,j,k)-recv1(ivar,ii-1,j,k)) & 
     &                                         *dxpr
              gradr = (recv1(ivar,ii+1,j,k)-recv1(ivar,ii,j,k)) & 
     &                                         *dxpr
              gradc = (recv1(ivar,ii+1,j,k)-recv1(ivar,ii-1,j,k)) & 
     &                                         *dxpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss, & 
     &                                       2.*gradl*ss))
              is = mod(i-nguard-1+largei,2)
              sdx = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,ii,j,k)+gradm*sdx
        enddo
        enddo
        enddo

!
! Perform sweep in y direction
        do k=kcl,kcu+iface_off
        do j=jlow,jhi
        do i=icl,icu
        recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo

        do k=kcl,kcu+iface_off
        do i=icl,icu
        do j=jcl,jcu
           jj = (j-nguard-1+largei)/2+nguard+1+joff-largei/2
              gradl = (recv1(ivar,i,jj,k)-recv1(ivar,i,jj-k2d,k)) & 
     &                                         *dypr
              gradr = (recv1(ivar,i,jj+k2d,k)-recv1(ivar,i,jj,k)) & 
     &                                         *dypr
              gradc = (recv1(ivar,i,jj+k2d,k)-recv1(ivar,i,jj-k2d,k)) & 
     &                                         *dypr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss, & 
     &                                       2.*gradl*ss))
              is = mod(j-nguard-1+largei,2)
              sdy = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,i,jj,k)+gradm*sdy
        enddo
        enddo
        enddo

! save new interpolated solution

      if (curvilinear) then
! Generate cell volume data for the current child block
      call amr_block_geometry(lb,mype)

! Divide interpolated solution by the 2**ndim*volume of the child cells
      factor = .5**ndim

      endif


        do k=kcl,kcu+iface_off
        do j=jcl,jcu
        do i=icl,icu
        if (curvilinear) then
           facevarz1(ivar,i,j,k,idest)=temp(ivar,i,j,k) & 
     &                                *factor/cell_area3(i,j,k)
        else
           facevarz1(ivar,i,j,k,idest)=temp(ivar,i,j,k)
        endif   
        enddo
        enddo
        enddo


        endif                              ! end of iface=3 


      return
      end subroutine amr_1blk_fc_prol_user





