      subroutine amr_bc_block(jface,ibc,iopt,l,mype)

!------------------------------------------------------------------------
! This is a template. User editing is required before this routine 
! can be used.

! This routine sets the guard cell elements of the solution arrays 
! or workspace array on face jface of block l, using a boundary
! condition algorithm which the user must insert. Places where user
! editing is required are clearly indicated.
!
! Written :     Peter MacNeice          September 1997
!------------------------------------------------------------------------
!
! Arguments:
!       jface           number designating selected block face
!       ibc             integer specifying the particular boundary
!                       condition selected for this face.
!       iopt            a switch to control which data is updated
!                       iopt=1 will use 'unk' and 'facevarx(y)(z)'
!                       iopt=2 will use 'work'
!       l               number designating selected block 
!       mype            local processor number
!
!------------------------------------

      integer jface,iopt,l,mype


! include file to define physical qualities of the model and mesh
#include "physicaldata.fh"
      include 'workspace.fh'

! include file defining the tree
        include 'tree.fh'

! include file required for shmem library.
        include 'mpp/shmem.fh'
        include 'shmem_reduce.fh'

!------------------------------------



         ilbnd = 1
         iubnd = nxb+2*nguard
         jlbnd = 1
         jubnd = nyb+2*nguard*k2d
         klbnd = 1
         kubnd = nzb+2*nguard*k3d


! Limit the application of this routine to leaf blocks
! or the parents of leaf blocks.
      if(nodetype(l).le.2) then


! act on unk and facevar^s or on work?
      if(iopt.eq.1) then

        if(ibc.eq.-21) then


        if(jface.eq.1) then


! Loop over the cells of the selected face setting unk.
          do i=1,nguard
            unk(:,nguard+1-i,:,:,l)=unk(:,nguard+i,:,:,l)
          enddo

        elseif(jface.eq.2) then

          do i=1,nguard
            unk(:,nxb+nguard+i,:,:,l)=unk(:,nxb+nguard+1-i,:,:,l)
          enddo

        elseif(jface.eq.3) then

          do j=1,nguard
            unk(:,:,nguard+1-j,:,l)=unk(:,:,nguard+j,:,l)
          enddo

        elseif(jface.eq.4) then

          do j=1,nguard
            unk(:,:,nyb+nguard+j,:,l)=unk(:,:,nyb+nguard+1-j,:,l)
          enddo

        elseif(jface.eq.5) then

          do k=1,nguard
            unk(:,:,:,nguard+1-k,l)=unk(:,:,:,nguard+k,l)
          enddo

        elseif(jface.eq.6) then

          do k=1,nguard
            unk(:,:,:,nzb+nguard+k,l)=unk(:,:,:,nzb+nguard+1-k,l)
          enddo

        endif


! Now repeat for the facevar arrays if necessary.


       if(nfacevar.gt.0) then


       if(jface.eq.1) then
         do i=1,nguard
           facevarx(:,nguard+1-i,:,:,l)=
     .                               facevarx(:,nguard+1+i,:,:,l)
         enddo
         if(ndim.ge.2) then
         do i=1,nguard
           facevary(:,nguard+1-i,:,:,l)=facevary(:,nguard+i,:,:,l)
         enddo
         endif
         if(ndim.eq.3) then
         do i=1,nguard
           facevarz(:,nguard+1-i,:,:,l)=facevarz(:,nguard+i,:,:,l)
         enddo
         endif


       elseif(jface.eq.2) then
         do i=1,nguard
           facevarx(:,nxb+nguard+1+i,:,:,l) =
     .                             facevarx(:,nxb+nguard+1-i,:,:,l)
         enddo
         if(ndim.ge.2) then
         do i=1,nguard
           facevary(:,nxb+nguard+i,:,:,l) =
     .                             facevary(:,nxb+nguard+1-i,:,:,l)
         enddo
         endif
         if(ndim.eq.3) then
         do i=1,nguard
           facevarz(:,nxb+nguard+i,:,:,l) =
     .                            facevarz(:,nxb+nguard+1-i,:,:,l)
         enddo
         endif

       elseif(jface.eq.3) then
         do j=1,nguard
           facevarx(:,:,nguard+1-j,:,l)=facevarx(:,:,nguard+j,:,l)
         enddo
         do j=1,nguard
           facevary(:,:,nguard+1-j,:,l)=
     .                            facevary(:,:,nguard+1+j,:,l)
         enddo
         if(ndim.eq.3) then
         do j=1,nguard
           facevarz(:,:,nguard+1-j,:,l)=facevarz(:,:,nguard+j,:,l)
         enddo
         endif


       elseif(jface.eq.4) then
         do j=1,nguard
           facevarx(:,:,nyb+nguard+j,:,l) =
     .                             facevarx(:,:,nyb+nguard+1-j,:,l)
         enddo
         do j=1,nguard
           facevary(:,:,nyb+nguard+1+j,:,l) =
     .                             facevary(:,:,nyb+nguard+1-j,:,l)
         enddo
         if(ndim.eq.3) then
         do j=1,nguard
           facevarz(:,:,nyb+nguard+j,:,l) =
     .                            facevarz(:,:,nyb+nguard+1-j,:,l)
         enddo
         endif

       elseif(jface.eq.5) then
         do k=1,nguard
           facevarx(:,:,:,nguard+1-k,l)=facevarx(:,:,:,nguard+k,l)
         enddo
         do k=1,nguard
           facevary(:,:,:,nguard+1-k,l)=facevary(:,:,:,nguard+k,l)
         enddo
         if(ndim.eq.3) then
         do k=1,nguard
           facevarz(:,:,:,nguard+1-k,l)=
     .                    facevarz(:,:,:,nguard+1+k,l)
         enddo
         endif


       elseif(jface.eq.6) then
         do k=1,nguard
           facevarx(:,:,:,nzb+nguard+k,l) =
     .                             facevarx(:,:,:,nzb+nguard+1-k,l)
         enddo
         do k=1,nguard
           facevary(:,:,:,nzb+nguard+k,l) =
     .                             facevary(:,:,:,nzb+nguard+1-k,l)
         enddo
         if(ndim.eq.3) then
         do k=1,nguard
           facevarz(:,:,:,nzb+nguard+1+k,l) =
     .                            facevarz(:,:,:,nzb+nguard+1-k,l)
         enddo
         endif

       endif

       endif                            ! end of nfacevar if test
       endif                            ! end of ibc if test


      elseif(iopt.eq.2) then



        ilbnd = 1
        iubnd = nxb+2*nguard_work
        jlbnd = 1
        jubnd = nyb+2*nguard_work*k2d
        klbnd = 1
        kubnd = nzb+2*nguard_work*k3d



        if(ibc.eq.-21) then


        if(jface.eq.1) then

! Loop over the cells of the selected face setting work.
          do i=1,nguard_work
            work(nguard_work+1-i,:,:,l) =
     .                      work(nguard_work+i,:,:,l)
          enddo

        elseif(jface.eq.2) then

          do i=1,nguard_work
            work(nxb+nguard_work+i,:,:,l) =
     .                      work(nxb+nguard_work+1-i,:,:,l)
          enddo

        elseif(jface.eq.3) then

          do j=1,nguard_work
            work(:,nguard_work+1-j,:,l) = 
     .                      work(:,nguard_work+j,:,l)
          enddo

        elseif(jface.eq.4) then

          do j=1,nguard_work
            work(:,nyb+nguard_work+j,:,l) = 
     .                      work(:,nyb+nguard_work+1-j,:,l)
          enddo

        elseif(jface.eq.5) then

          do k=1,nguard_work
            work(:,:,nguard_work+1-k,l) = 
     .                      work(:,:,nguard_work+k,l)
          enddo

        elseif(jface.eq.6) then

          do k=1,nguard_work
            work(:,:,nzb+nguard_work+k,l) =
     .                      work(:,:,nzb+nguard_work+1-k,l)
          enddo

        endif                              ! end of jface if test
        endif                              ! end of ibc if test
       

       endif                               ! end of iopt if test

       endif                               ! end of nodetype if test


      return
      end
