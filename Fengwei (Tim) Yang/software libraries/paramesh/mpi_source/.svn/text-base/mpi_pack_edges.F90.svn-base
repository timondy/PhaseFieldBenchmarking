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
!#define DEBUGX

      subroutine mpi_pack_edges(mype,nprocs,buf_dim,S_buffer, & 
     &                          offset)

!------------------------------------------------------------------------
!
! This subroutine packs all edge data which are to be sent from mype.
!
!
! Written :     Peter MacNeice, Maharaj Bhat & Michael Gehmeyr    July 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      nprocs         number of processors
!      buf_dim        dimension of buffer
!      S_buffer       send buffer 
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_get_edge_buffer

      implicit none

      include 'mpif.h'

      integer, intent(in)  ::  mype,nprocs
      integer, intent(in)  ::  buf_dim,offset
      real,    intent(out) ::  S_buffer(buf_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: loc_message_size(27),dtype
      integer :: tot_no_blocks_to_be_received,lindex
      integer :: nguard0 
      integer :: lb,irpe
      integer :: index, array_size, ksw
      integer :: edgepack_debug(1000)
      integer :: jrpe, isize, iblk, itype, iseg, next_pe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs

! define starting end ending indices for the send and recv buffers

      is_buf = 0
      ir_buf = 0

      array_size = nedges*2*2*( (nyb+k2d)*(nzb+k3d) & 
     &                         +(nxb+1)  *(nzb+k3d)*k2d & 
     &                         +(nxb+1)  *(nyb+k3d)*k3d )
      array_size = array_size + 2

#ifdef DEBUG
      write(*,*) 'pe ',mype,' array_size ',array_size
#endif /* DEBUG */


        loc_message_size(:) = 3
! size of messages for single faces
        loc_message_size(13) = nedges*( (nyb+k2d)*k2d*(nzb+k3d) & 
     &                                 +(nyb+1)*(nzb+k3d)*k3d   )+ 3
        loc_message_size(15) = nedges*( (nyb+k2d)*k2d*(nzb+k3d) & 
     &                                 +(nyb+1)*(nzb+k3d)*k3d   )+ 3
        loc_message_size(11) = nedges*( (nxb+1)*(nzb+k3d) & 
     &                                 +(nxb+1)*(nzb+k3d)*k3d   )+ 3
        loc_message_size(17) = nedges*( (nxb+1)*(nzb+k3d) & 
     &                                 +(nxb+1)*(nzb+k3d)*k3d   )+ 3
        loc_message_size( 5) = nedges*( (nxb+1)*(nyb+k2d) & 
     &                                 +(nxb+1)*(nyb+k2d)*k2d   )+ 3
        loc_message_size(23) = nedges*( (nxb+1)*(nyb+k2d) & 
     &                                 +(nxb+1)*(nyb+k2d)*k2d   )+ 3

        ksw = 0
        if(ndim.eq.3) ksw = 1
        if(l2p5d.eq.1) then
          loc_message_size(13) = nedges*2*(nyb+1)*k2d + 3
          loc_message_size(15) = nedges*2*(nyb+1)*k2d + 3
          loc_message_size(11) = nedges*2*(nxb+1) + 3
          loc_message_size(17) = nedges*2*(nxb+1) + 3
        endif
        if (l2p5d.eq.1 .or. k3d.eq.1) ksw = 1

! size of messages for single edges (note the second factor of 2 is
! included because each edge is included in 2 arrays, and it is easier
! to send data for both, than to ensure consistency later, since if
! we send to only 1 array we would not know which of the two 
! had the most up to date data.)
        loc_message_size(2) = nedges*2*(nxb+1) + 3
        loc_message_size(8) = nedges*2*(nxb+1) + 3
        loc_message_size(20) = nedges*2*(nxb+1) + 3
        loc_message_size(26) = nedges*2*(nxb+1) + 3
        loc_message_size(4) = nedges*2*(nyb+1)*k2d + 3
        loc_message_size(6) = nedges*2*(nyb+1)*k2d + 3
        loc_message_size(22) = nedges*2*(nyb+1)*k2d + 3
        loc_message_size(24) = nedges*2*(nyb+1)*k2d + 3
        loc_message_size(10) = nedges*2*(nzb+k3d)*ksw + 3
        loc_message_size(12) = nedges*2*(nzb+k3d)*ksw + 3
        loc_message_size(16) = nedges*2*(nzb+k3d)*ksw + 3
        loc_message_size(18) = nedges*2*(nzb+k3d)*ksw + 3
        loc_message_size(14) = nedges*2*(1+ksw)* & 
     &                           ( (nyb+k2d)*(nzb+k3d) & 
     &                           +(nxb+1)  *(nzb+k3d)*k2d & 
     &                           +(nxb+1)  *(nyb+k2d)*k3d ) + 3
        loc_message_size(:) = loc_message_size(:) + offset



      index = 0
      jrpe = 0
      do irpe = 1,nprocs      ! define send buffer indices
        if (commatrix_send(irpe).gt.0) then

           jrpe = jrpe + 1
           isize = 0
           do iblk = 1,commatrix_send(irpe)
            itype = to_be_sent(3,iblk,jrpe)
            isize = isize + loc_message_size(itype)
#ifdef DEBUGX
       write(*,*) 'edge: pe ',mype,' sizing send buf to pe ',irpe, & 
     &          ' adding message type ',itype,' size ', & 
     &        loc_message_size(itype), & 
     &        ' accumulated size ',isize
#endif /* DEBUGX */
           enddo
           is_buf(1,irpe) = index + 1
           is_buf(2,irpe) = index + isize
           
           index = index + isize
        endif
      enddo
           index = 0

! set up a pointer list to the start address in recv_buffer for each
! block of information in the received messages
      tot_no_blocks_to_be_received = sum(commatrix_recv(:))
      if(allocated(mess_segment_loc)) deallocate(mess_segment_loc)
      allocate(mess_segment_loc(tot_no_blocks_to_be_received))
      mess_segment_loc = 0
      iseg = 0
      lindex = 0
      jrpe =  0
      do irpe = 1,nprocs      ! define recv buffer indices
        if (commatrix_recv(irpe).gt.0) then

           jrpe = jrpe + 1
           isize = 0
           do iblk = 1,commatrix_recv(irpe)
             itype = to_be_received(3,iblk,jrpe)
             isize = isize + loc_message_size(itype)
             iseg = iseg+1
             mess_segment_loc(iseg) = lindex+1
             lindex = lindex+loc_message_size(itype)
#ifdef DEBUGX
       write(*,*) 'edge: pe ',mype,' sizing recv buf from pe ',irpe, & 
     &          ' adding message type ',itype,' size ', & 
     &        loc_message_size(itype), & 
     &        ' accumulated size ',isize,' iseg ',iseg, & 
     &        ' mess_segment_loc ',mess_segment_loc(iseg), & 
     &        ' lindex ',lindex
#endif /* DEBUGX */
           enddo
           ir_buf(1,irpe) = index + 1
           ir_buf(2,irpe) = index + isize
           
           index = index + isize
        endif
      enddo

! pack buffer succinctly with all arrays of all blocks that are 
! to be packed for irpe

      index = 1 


      S_buffer = 0
#ifdef DEBUG
        write(*,*) 'pe ',mype,' nprocs ',nprocs, & 
     &             ' start edgepacking'
#endif /* DEBUG */

      next_pe = 0
      do irpe = 1,nprocs      ! define recv buffer indices

        edgepack_debug(irpe) = index
#ifdef DEBUG
        write(*,*) 'pe ',mype,' irpe ',irpe,' commatrix ', & 
     &        commatrix(mype+1,irpe)
#endif /* DEBUG */
        if (commatrix_send(irpe).gt.0) then
          next_pe = next_pe+1
          do iblk = 1,commatrix_send(irpe)
            if(to_be_sent(1,iblk,next_pe).gt.0) then


              lb = to_be_sent(1,iblk,next_pe)
              dtype = to_be_sent(3,iblk,next_pe)

#ifdef DEBUGX
        write(*,*) 'pe ',mype,' :edgepack for rempe ',irpe, & 
     &     ' in buffer layer ', next_pe,' blk ', iblk, & 
     &        ' from local lb ',lb,' index ',index,' dtype ',dtype
#endif /* DEBUGX */
                                  ! pack all arrays for lb into buffer
              call mpi_get_edge_buffer( mype,lb,dtype,index, & 
     &                             buf_dim,S_buffer)

            endif
          enddo
        endif
      enddo


#ifdef DEBUGX
      do j=1,nprocs
      write(*,*) 'pe ',mype,' edgepack_debug(',j,') = ', & 
     &                   edgepack_debug(j)
      enddo
#endif /* DEBUGX */

#ifdef DEBUGX
      if (index .ne. buf_dim) then
      print *,' '
      print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
      print *,' Sbuffer is too big by ',size(S_buffer,dim=1)-index,mype
      print *,' index = ',index,' buf_dim = ',buf_dim,mype
      print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
      print *,' '
      end if
#endif

      return
      end subroutine mpi_pack_edges
