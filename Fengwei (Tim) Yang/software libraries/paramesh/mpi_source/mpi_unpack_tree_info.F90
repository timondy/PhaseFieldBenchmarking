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

      subroutine mpi_unpack_tree_info(mype,nprocs,buf_dim_bytes, & 
     &                             buf_dim,R_buffer)

!------------------------------------------------------------------------
!
! This subroutine unpacks tree information which has been sent.
!
!
! Written :     Peter MacNeice, Maharaj Bhat & Michael Gehmeyr  June 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      nprocs         number of processors
!      buf_dim        dimension of buffer
!      buf_dim_bytes  size of message in buffer, in bytes
!      R_buffer       receive buffer 
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none

      include 'mpif.h'

      integer, intent(in) :: mype,nprocs,buf_dim,buf_dim_bytes
      real,    intent(inout) ::  R_buffer(buf_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: i,j,k,ierrorcode,ierr,position
      integer :: lblk, lnumb, lb, ii, i_pe, iii, jj, n
      integer :: len_surr_blks,kkkk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
#ifdef SAVE_MORTS & 
     &                +  6*3*(1+2*k2d)*(1+2*k3d)
#endif

      lnumb = sum(commatrix_recv(:))
      if(lnumb.gt.maxblocks_alloc) then
            call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

!debugging only
      do kkkk=0,nprocs-1
      print *,'CEG - MPI_BARRIER in mpi_unpack_tree_info.F90', kkkk
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.kkkk) then
!debugging only

      position = 1

      ii = 0
      lblk = 0 
      do i=1,nprocs
        if(commatrix_recv(i).gt.0) then
          ii = ii + 1
          i_pe = pe_source(ii)

          do j=1,commatrix_recv(i)

            lblk = lblk + 1
            lb = lblk + strt_buffer - 1

                                  ! unpack all arrays from buffer into lb
#ifdef DEBUG
            write(*,*) 'pe ',mype,' lblk ',lblk, & 
     &                 ' unpacking tree info ' & 
     &  ,' from pe ',i,' expecting ',commatrix_recv(i)
#endif /* DEBUG */

        do iii = 1,mdim
           coord(iii,lb) = R_buffer(position)
           position = position + 1
        end do
#ifdef DEBUG
            write(*,*) 'pe ',mype,' lblk ',lblk,' coord ' & 
     &         ,coord(:,lb)
#endif /* DEBUG */

        do iii = 1,mdim
           bsize(iii,lb) = R_buffer(position)
           position = position + 1
        end do

        do iii = 1,2
           do jj = 1,mdim
              bnd_box(iii,jj,lb) = R_buffer(position)
              position = position + 1
           end do
        end do
#ifdef DEBUG
            write(*,*) 'pe ',mype,' lblk ',lblk,' bnd_box ' & 
     &         ,bnd_box(:,:,lb)
#endif /* DEBUG */

        do iii = 1,2
           parent(iii,lb) = int(R_buffer(position))
           position = position + 1
        end do
#ifdef DEBUG
            write(*,*) 'pe ',mype,' lblk ',lblk,' parent ' & 
     &         ,parent(:,lb)
#endif /* DEBUG */
#ifdef DEBUG
        write(*,*) 'pe ',mype,' bef child position ' & 
     &        ,position,' shape R_buf ',shape(R_buffer)
#endif /* DEBUG */

        do iii = 1,2
           do jj = 1,mchild
             child(iii,jj,lb) = int(R_buffer(position))
             position = position + 1
           end do
        end do

#ifdef DEBUG
            write(*,*) 'pe ',mype,' lblk ',lblk,' child'
#endif /* DEBUG */
        which_child(lb) = int(R_buffer(position))
        position = position + 1

        if (int(R_buffer(position)) == 1) then
           newchild(lb) = .true.
        else
           newchild(lb) = .false.
        end if
        position = position + 1
#ifdef DEBUG
            write(*,*) 'pe ',mype,' lblk ',lblk,' new_child'
#endif /* DEBUG */

        do iii = 1,2
           do jj = 1,mfaces
              neigh(iii,jj,lb) = int(R_buffer(position))
              position = position + 1
           end do
        end do

        lrefine(lb) = int(R_buffer(position))
        position = position + 1

        nodetype(lb) = int(R_buffer(position))
        position = position + 1

        do n = 1,3
           do k = 1,1+2*k3d
              do jj = 1,1+2*k2d
                 do iii = 1,3
                    surr_blks(n,iii,jj,k,lb) = & 
     &                   int(R_buffer(position))
                    position = position + 1
                 end do
              end do
           end do
        end do

#ifdef DEBUG
         write(*,*) 'pe ',mype,' lblk ',lblk,' surr_blks' & 
     &          ,surr_blks(:,:,:,1,lb) & 
     &          ,' current position ',position
#endif /* DEBUG */

#ifdef SAVE_MORTS
      do n = 1,6
         do k = 1,1+2*k3d
            do jj = 1,1+2*k2d
               do iii = 1,3
                  surr_morts(n,iii,jj,k,lb) = & 
     &                 int(R_buffer(position))
                  position = position + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

#ifdef DEBUG
        write(*,*) 'pe ',mype,' bef empty position ' & 
     &        ,position,' shape R_buf ',shape(R_buffer)
#endif /* DEBUG */

        empty(lb) = int(R_buffer(position))
        position = position + 1
#ifdef DEBUG
            write(*,*) 'pe ',mype,' lblk ',lblk,' empty'
#endif /* DEBUG */

        if (int(R_buffer(position)) == 1) then
           ldtcomplete(lb) = .true.
        else
           ldtcomplete(lb) = .false.
        end if
        position = position + 1
 
#ifdef DEBUG
! record original address of the source of this tree info
!            laddress(:,lb) = to_be_received(:,j,i_pe)

            write(*,*) 'pe ',mype,' lblk ',lblk & 
     &     ,' tree info unpacked into ',lb,' coord ',coord(:,lb), & 
     &     ' size ',bsize(:,lb),' bnd_box ',bnd_box(:,:,lb), & 
     &     ' parent ',parent(:,lb), & 
     &     ' child ',child(:,:,lb), & 
     &     ' which_child ',which_child(lb) & 
     &     ,' neigh ',neigh(:,:,lb),' lrefine ',lrefine(lb) & 
     &     ,' nodetype ',nodetype(lb) & 
     &     ,' surr_blks ',surr_blks(:,:,:,:,lb)
#ifdef SAVE_MORTS & 
     &     ,' surr_morts ',surr_morts(:,:,:,:,lb)
#endif & 
     &     ,' empty ',empty(lb) & 
     &     ,' original address ',laddress(:,lb)
#endif /* DEBUG */

          enddo
        endif
      enddo

!debugging only
      endif ! (mype.eq.kkkk) then
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo ! kkkk=0,nprocs-1
!debugging only
      return
      end subroutine mpi_unpack_tree_info
