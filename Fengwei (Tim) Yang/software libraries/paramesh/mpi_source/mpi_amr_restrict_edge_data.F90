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
!#define DEBUG

      subroutine amr_restrict_edge_data(mype)




!------------------------------------------------------------------------
!
! This routine does the data averaging on cell edges required when a 
! child block passes data back to its parent. The parent receives data 
! at the block boundary only.
!
! This routine provides a mechanism for passing data defined at block
! boundaries from leaf blocks back to their parents.
! The averaging rules used to combine interface values on the finer
! mesh to construct interface values on the coarser parent mesh are
! specified by the user who provides a function called amr_restrict_edge
! to do this.
!
! This routine is only relevant for schemes with even number of grid points.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      use paramesh_interfaces, only : amr_restrict_edge
      use paramesh_mpi_interfaces, only : mpi_put_edge_buffer_1blk

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: mype

!------------------------------------
! local arrays

      integer :: remote_pe,remote_block
      integer :: remote_pe2,remote_block2
      integer,save :: anodetype(1)
      integer :: cnodetype,cneigh(2,6)
      save cnodetype,cneigh
      integer :: dtype, vtype, index, index0
      integer :: i, j, k, n
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ia0, ib0, ja0, jb0, ka0, kb0
      integer :: nprocs
      integer :: nguard0
      integer :: ioff, joff, koff, iface, jface, kface, ivar, iblk
      integer :: ii, jj, kk
      integer :: lb, ierrorcode, ierr, ich, icoord, jchild


      logical :: lfound

!------------------------------------

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      nguard0 = nguard*npgs

      if(mpi_pattern_id.ne.40 .and. nprocs.gt.1) then
        write(*,*) 'Paramesh error : wrong pattern being', & 
     &' used for pre-communication for edge restrict : ', & 
     &' Fix - insert appropriate call to mpi_amr_comm_setup'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a parent block of at least one leaf node?
      if(nodetype(lb).eq.2) then

! If yes then cycle through its children.
       do ich=1,nchild

       jchild = ich
       remote_pe = child(2,ich,lb)
       remote_block  = child(1,ich,lb)
       remote_pe2 = child(2,ich,lb)
       remote_block2  = child(1,ich,lb)

! if (remote_block2,remote_pe2) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.
          if(remote_pe2.ne.mype) then
             lfound = .false.
             iblk = ladd_strt(remote_pe2)
             do while(.not.lfound.and.iblk.le.ladd_end(remote_pe2))
#ifdef DEBUG
             write(*,*) 'pe ',mype,' searching buffer for ', & 
     &            remote_block2,remote_pe2,' current buffer entry ', & 
     &          ' iblk ',iblk,' laddress ',laddress(:,iblk)
#endif /* DEBUG */
              if(remote_block2.eq.laddress(1,iblk).and. & 
     &             remote_pe2 .eq.laddress(2,iblk) ) then
                remote_block2 = iblk
                remote_pe2    = mype
                lfound = .true.
#ifdef DEBUG
             write(*,*) 'pe ',mype,' remote block ', & 
     &          remote_block2,remote_pe2,' located in buffer slot ', & 
     &          iblk
#endif /* DEBUG */
              else
                iblk = iblk+1
              endif
            enddo
          endif


! Is this child a leaf block(nodetype=1)? 
! If it is then fetch its data.
       cnodetype = nodetype(remote_block2)
       if(cnodetype.eq.1) then

! fetch child's neighbors. This info will be needed when we are
! advancing the solution on all levels.
         cneigh(:,:) = neigh(:,:,remote_block2)

         if (remote_pe == mype .and. remote_block <= lnblocks) then
         
         else

! find remote block in buffer

           call mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)

         end if

         do icoord=1,ndim

           if (remote_pe == mype .and. remote_block <= lnblocks) then

           if(icoord.eq.1) then
             recvarx1e(:,:,:,:) = bedge_facex_y(:,:,:,:,remote_block)
             if((ndim.eq.3).or.(l2p5d.eq.1)) then
               recvarx2e(:,:,:,:) = bedge_facex_z(:,:,:,:,remote_block)
             end if
           elseif(icoord.eq.2) then
             recvary1e(:,:,:,:) = bedge_facey_x(:,:,:,:,remote_block)
             if((ndim.eq.3).or.(l2p5d.eq.1)) then
                recvary2e(:,:,:,:) = bedge_facey_z(:,:,:,:,remote_block)
             end if
           elseif(icoord.eq.3) then
             recvarz1e(:,:,:,:) = bedge_facez_x(:,:,:,:,remote_block)
             recvarz2e(:,:,:,:) = bedge_facez_y(:,:,:,:,remote_block)
           endif

         end if                    ! end if (remote_pe


! compute the offset in the parent block appropriate for this child
           iface = mod(jchild-1,2)
           jface = mod((jchild-1)/2,2)
           kface = mod((jchild-1)/4,2)
           ioff = iface*nxb/2
           joff = jface*nyb/2
           koff = kface*nzb/2

! Compute restricted data from the data in the buffer and
! update only boundary values on the parent block
           if(icoord.eq.1) then

             i = iface+1
! apply, only if the appropriate child neighbor does not exist
             if(cneigh(1,i).gt.-20.and.cneigh(1,i).lt.0) then

             call amr_restrict_edge(icoord)
             do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
               kk = (k-nguard0*k3d)/2+nguard0*k3d+1
               do j=1+nguard0,nyb+nguard0,2
                 jj = (j-nguard0)/2+nguard0+1
                 bedge_facex_y(:,i,jj+joff,kk+koff, & 
     &                                    lb) = recvarx1e(:,i,j,k)
               enddo
             enddo

             if((ndim.eq.3).or.(l2p5d.eq.1)) then
               do k=1+nguard0*k3d,nzb+nguard0*k3d,2
                 kk = (k-nguard0*k3d)/2+nguard0*k3d+1
                 do j=1+nguard0,nyb+nguard0+1,2
                   jj = (j-nguard0)/2+nguard0+1
                   bedge_facex_z(:,i,jj+joff,kk+koff, & 
     &                                    lb) = recvarx2e(:,i,j,k)
                 enddo
               enddo
             endif

             endif


           elseif(icoord.eq.2) then

             j = jface+1
! apply, only if the appropriate child neighbor does not exist
             if(cneigh(1,j+2).gt.-20.and.cneigh(1,j+2).lt.0) then

             call amr_restrict_edge(icoord)
             do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
               kk = (k-nguard0*k3d)/2+nguard0*k3d+1
               do i=1+nguard0,nxb+nguard0,2
                 ii = (i-nguard0)/2+nguard0+1
                 bedge_facey_x(:,ii+ioff,j,kk+koff, & 
     &                                    lb) = recvary1e(:,i,j,k)
               enddo
             enddo

             if((ndim.eq.3).or.(l2p5d.eq.1)) then
               do k=1+nguard0*k3d,nzb+nguard0*k3d,2
                 kk = (k-nguard0*k3d)/2+nguard0*k3d+1
                 do i=1+nguard0,nxb+nguard0+1,2
                   ii = (i-nguard0)/2+nguard0+1
                   bedge_facey_z(:,ii+ioff,j,kk+koff, & 
     &                                    lb) = recvary2e(:,i,j,k)
                 enddo
               enddo
             endif

             endif

           elseif(icoord.eq.3) then

             k = kface+1
! apply, only if the appropriate child neighbor does not exist
             if(cneigh(1,k+4).gt.-20.and.cneigh(1,k+4).lt.0) then

             call amr_restrict_edge(icoord)
             do j=1+nguard0,nyb+nguard0+1,2
               jj = (j-nguard0)/2+nguard0+1
               do i=1+nguard0,nxb+nguard0,2
                 ii = (i-nguard0)/2+nguard0+1
                 bedge_facez_x(:,ii+ioff,jj+joff,k, & 
     &                                    lb) = recvarz1e(:,i,j,k)
               enddo
             enddo
             do j=1+nguard0,nyb+nguard0,2
               jj = (j-nguard0)/2+nguard0+1
               do i=1+nguard0,nxb+nguard0+1,2
                 ii = (i-nguard0)/2+nguard0+1
                 bedge_facez_y(:,ii+ioff,jj+joff,k, & 
     &                                    lb) = recvarz2e(:,i,j,k)
               enddo
             enddo

           endif

           endif

       enddo

       endif

       enddo

      endif

      enddo
      endif

      return
      end subroutine amr_restrict_edge_data
