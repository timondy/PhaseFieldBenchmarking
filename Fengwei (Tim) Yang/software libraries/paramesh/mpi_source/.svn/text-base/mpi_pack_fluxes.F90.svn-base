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
      subroutine mpi_pack_fluxes(mype,nprocs, & 
     &                           buf_dim,S_buffer,offset,flux_dir)

!------------------------------------------------------------------------
!
! This subroutine packs all fluxes which are to be sent from mype.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
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

      use paramesh_mpi_interfaces, only : mpi_get_flux_buffer

      implicit none

      include 'mpif.h'

      integer, intent(in)  ::  buf_dim,offset
      real,    intent(out) ::  S_buffer(buf_dim)
      integer, intent(in)  ::  mype,nprocs
      integer, intent(in), optional :: flux_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lb,irpe
      integer :: index
      integer :: loc_message_size(27),dtype
      integer :: tot_no_blocks_to_be_received,lindex
      integer :: nguard0 
      integer :: jrpe, isize, iblk, itype, iseg, next_pe
      integer :: flux_dirt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs

!
! define starting end ending indices for the send and recv buffers

      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      end if

      is_buf = 0
      ir_buf = 0

      if (flux_dirt == 1) then

      loc_message_size(:) = 3
      loc_message_size(13) = nfluxes*nyb*nzb + 3
      loc_message_size(15) = nfluxes*nyb*nzb + 3
      loc_message_size(11) = 3
      loc_message_size(17) = 3
      loc_message_size( 5) = 3
      loc_message_size(23) = 3
      loc_message_size(14) = nfluxes*2* & 
     &                (nyb*nzb+0+0) + 3
      loc_message_size(:) = loc_message_size(:) + offset

      elseif (flux_dirt == 2) then

      loc_message_size(:) = 3
      loc_message_size(13) = 3
      loc_message_size(15) = 3
      loc_message_size(11) = nfluxes*nxb*nzb + 3
      loc_message_size(17) = nfluxes*nxb*nzb + 3
      loc_message_size( 5) = 3
      loc_message_size(23) = 3
      loc_message_size(14) = nfluxes*2* & 
     &                (0+nxb*nzb*k2d+0) + 3
      loc_message_size(:) = loc_message_size(:) + offset

      elseif (flux_dirt == 3) then

      loc_message_size(:) = 3
      loc_message_size(13) = 3
      loc_message_size(15) = 3
      loc_message_size(11) = 3
      loc_message_size(17) = 3
      loc_message_size( 5) = nfluxes*nxb*nyb + 3
      loc_message_size(23) = nfluxes*nxb*nyb + 3
      loc_message_size(14) = nfluxes*2* & 
     &                (0+0+nxb*nyb*k3d) + 3
      loc_message_size(:) = loc_message_size(:) + offset

      else

      loc_message_size(:) = 3
      loc_message_size(13) = nfluxes*nyb*nzb + 3
      loc_message_size(15) = nfluxes*nyb*nzb + 3
      loc_message_size(11) = nfluxes*nxb*nzb + 3
      loc_message_size(17) = nfluxes*nxb*nzb + 3
      loc_message_size( 5) = nfluxes*nxb*nyb + 3
      loc_message_size(23) = nfluxes*nxb*nyb + 3
      loc_message_size(14) = nfluxes*2* & 
     &                (nyb*nzb+nxb*nzb*k2d+nxb*nyb*k3d) + 3
      loc_message_size(:) = loc_message_size(:) + offset

      end if

      index = 0
      jrpe = 0
      do irpe = 1,nprocs      ! define send buffer indices
        if (commatrix_send(irpe).gt.0) then

           jrpe = jrpe + 1
           isize = 0
           do iblk = 1,commatrix_send(irpe)
            itype = to_be_sent(3,iblk,jrpe)
            isize = isize + loc_message_size(itype)
#ifdef DEBUG
       write(*,*) 'flux: pe ',mype,' sizing send buf to pe ',irpe, & 
     &          ' adding message type ',itype,' size ', & 
     &        loc_message_size(itype), & 
     &        ' accumulated size ',isize
#endif /* DEBUG */
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
       write(*,*) 'flux: pe ',mype,' sizing recv buf from pe ',irpe, & 
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
        write(*,*) 'pe ',mype,' nprocs ',nprocs,' start packing'
#endif /* DEBUG */

      next_pe = 0
      do irpe = 1,nprocs      ! define recv buffer indices
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

#ifdef DEBUG
        write(*,*) 'pe ',mype,' :pack for rempe ',irpe, & 
     &     ' in buffer layer ', next_pe,' blk ', iblk, & 
     &        ' from local lb ',lb,' index ',index
#endif /* DEBUG */
                                  ! pack all arrays for lb into buffer
              call mpi_get_flux_buffer( mype,lb,dtype,index, & 
     &                             buf_dim,S_buffer,flux_dir)

            endif
          enddo
        endif
      enddo


#ifdef DEBUG
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
      end subroutine mpi_pack_fluxes
