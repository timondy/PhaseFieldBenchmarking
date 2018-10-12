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

      subroutine mpi_pack_tree_info(mype,nprocs,buf_dim_bytes, & 
     &                           buf_dim,S_buffer)


!#define DEBUG
!
! DESIGN ISSUES :
!   Improve communication in this routine by defining an mpi datatype
!   for the combined real+integer message info.
!
!------------------------------------------------------------------------
!
! This subroutine packs all blocks which are to be sent from mype.
!
!
! Written :     Peter MacNeice, Maharaj Bhat & Michael Gehmeyr   June 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      nprocs         number of processors
!      buf_dim_bytes  the no of bytes which are sent in this message
!      buf_dim        dimension of buffer
!      S_buffer       send buffer 
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      implicit none

      include 'mpif.h'

      integer, intent(in)  :: mype,nprocs
      integer, intent(in)  :: buf_dim,buf_dim_bytes
      real,    intent(out) :: S_buffer(buf_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: array_size,position
      integer :: i,j,k,lb,irpe,isize,next_pe,iblk,ii
      integer :: index, len_surr_blks
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! define starting end ending indices for the send and recv buffers


      is_buf = 0
      ir_buf = 0

      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
#ifdef SAVE_MORTS & 
     &                +6*3*(1+2*k2d)*(1+2*k3d)
#endif

      array_size = 32 + 16 + len_surr_blks

      index = 0
      do irpe = 1,nprocs      ! define send buffer indices
        if (commatrix_send(irpe).gt.0) then

           isize = array_size*commatrix_send(irpe)
           is_buf(1,irpe) = index + 1
           is_buf(2,irpe) = index + isize
           
           index = index + isize
        endif
      enddo
           index = 0
      do irpe = 1,nprocs      ! define recv buffer indices
        if (commatrix_recv(irpe).gt.0) then

           isize = array_size*commatrix_recv(irpe)
           ir_buf(1,irpe) = index + 1
           ir_buf(2,irpe) = index + isize
           
           index = index + isize
        endif
      enddo

#ifdef DEBUG
      write(*,*) 'pack tree ',mype,' is and ir_buf set' & 
     &    ,' commatrix_send ', commatrix_send, & 
     &     ' to_be_sent ',to_be_sent(1,:,:)
#endif /* DEBUG */

! pack buffer succinctly with all arrays of all blocks that are 
! to be packed for irpe

!      S_buffer = 0
      position = 1

      next_pe = 0
      do irpe = 1,nprocs      ! define recv buffer indices
        if (commatrix_send(irpe).gt.0) then
          next_pe = next_pe+1
          do iblk = 1,commatrix_send(irpe)
            if(to_be_sent(1,iblk,next_pe).gt.0) then

              lb = to_be_sent(1,iblk,next_pe)

                                  ! pack all arrays for lb into buffer
#ifdef DEBUG
      write(*,*) 'pack tree ',mype,' packing blk ',lb
#endif /* DEBUG */

              do i = 1,mdim
                 S_buffer(position) = coord(i,lb)
                 position = position + 1
              end do

              do i = 1,mdim
                 S_buffer(position) = bsize(i,lb)
                 position = position + 1
              end do

              do i = 1,2
                 do j = 1,mdim
                    S_buffer(position) = bnd_box(i,j,lb)
                    position = position + 1
                 end do
              end do

              do i = 1,2
                 S_buffer(position) = real(parent(i,lb))
                 position = position + 1
              end do

              do i = 1,2
                 do j = 1,mchild
                    S_buffer(position) = child(i,j,lb)
                    position = position + 1
                 end do
              end do


              S_buffer(position) = real(which_child(lb))
              position = position + 1

              if (newchild(lb)) then
                 S_buffer(position) = real(1)
              else
                 S_buffer(position) = real(0)
              end if
              position = position + 1

              do i = 1,2
                 do j = 1,mfaces
                    S_buffer(position) = real(neigh(i,j,lb))
                    position = position + 1
                 end do
              end do

              S_buffer(position) = real(lrefine(lb))
              position = position + 1

              S_buffer(position) = real(nodetype(lb))
              position = position + 1

              do ii = 1,3
                 do k = 1,1+2*k3d
                    do j = 1,1+2*k2d
                       do i = 1,3
                          S_buffer(position) = & 
     &                         real(surr_blks(ii,i,j,k,lb))
                          position = position + 1
                       end do
                    end do
                 end do
              end do

#ifdef SAVE_MORTS
      do ii = 1,6
         do k = 1,1+2*k3d
            do j = 1,1+2*k2d
               do i = 1,3
                  S_buffer(index) = & 
     &                 real(surr_morts(ii,i,j,k,lb))
                  index = index + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

              S_buffer(position) = real(empty(lb))
              position = position + 1

              if (ldtcomplete(lb)) then
                 S_buffer(position) = real(1)
              else
                 S_buffer(position) = real(0)
              end if

              position = position + 1

            endif
          enddo
        endif
      enddo

#undef DEBUG
      return
      end subroutine mpi_pack_tree_info
