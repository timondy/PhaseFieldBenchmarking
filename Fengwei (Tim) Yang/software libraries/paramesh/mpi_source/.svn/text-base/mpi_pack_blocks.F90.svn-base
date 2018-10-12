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
      subroutine mpi_pack_blocks(mype,nprocs,iopt, & 
     &                           lcc,lfc,lec,lnc, & 
     &                           buf_dim,S_buffer,offset, & 
     &                           nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine packs all blocks which are to be sent from mype.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      nprocs         number of processors
!      iopt           option setting for work array
!      lcc            logical switch controlling whether unk data
!                     is packed
!      lfc            logical switch controlling whether facevar data
!                     is packed
!      lec            logical switch controlling whether unk_e_? data
!                     is packed
!      lnc            logical switch controlling whether unk_n data
!                     is packed
!      buf_dim        dimension of buffer
!      S_buffer       send buffer 
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_get_buffer

      implicit none

      include 'mpif.h'

      integer, intent(in)  ::  mype,nprocs,iopt
      logical, intent(in)  ::  lcc,lfc,lec,lnc
      integer, intent(in)  ::  buf_dim,offset
      real,    intent(out) ::  S_buffer(buf_dim)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: loc_message_size(2*27),dtype
      integer :: tot_no_blocks_to_be_received,lindex
      integer :: lb, irpe
      integer :: index 
      integer :: ierrorcode,ierr
      integer :: invar, ibndvar, ivaredge, ivarcorn
      integer :: jrpe, isize, iblk, itype, iseg, next_pe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(iopt.gt.1.and.(lfc.or.lec.or.lnc)) then
         write(*,*) 'Paramesh error : calling mpi_pack_blocks with ', & 
     &              'inconsistent argument list - iopt is > 1 while ', & 
     &              'one or more of lfc, lec and lnc are set to true.'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
!
! define starting end ending indices for the send and recv buffers

      is_buf = 0
      ir_buf = 0

      l_datapacked = .false.
      if(iopt.gt.1) then
        l_datapacked(1) = .true.
      elseif(iopt.eq.1) then
        if(lcc) l_datapacked(2) = .true.
        if(lfc) l_datapacked(3) = .true.
        if(lec) l_datapacked(4) = .true.
        if(lnc) l_datapacked(5) = .true.
      endif

      invar = 0
      if (lcc) invar = nvar
      if(lcc.and.lguard_in_progress) invar = ngcell_on_cc

      ibndvar = 0
      if (lfc) ibndvar = nbndvar
      if(lfc.and.lguard_in_progress) ibndvar =  & 
     &                                    maxval(ngcell_on_fc(1:ndim))

      ivaredge = 0
      if (lec) ivaredge = nvaredge
      if(lec.and.lguard_in_progress) ivaredge =  & 
     &                                    maxval(ngcell_on_ec(1:3))

      ivarcorn = 0
      if (lnc) ivarcorn = nvarcorn
      if(lnc.and.lguard_in_progress) ivarcorn = ngcell_on_nc


      if(iopt.eq.1) then
        loc_message_size(:) = invar*message_size_cc(:) & 
     &                     +ibndvar*message_size_fc(:) & 
     &                     +ivaredge*message_size_ec(:) & 
     &                     +ivarcorn*message_size_nc(:) + 3 & 
     &                     + offset
      else
        loc_message_size = message_size_wk + 3 + offset
      endif

#ifdef DEBUG
      write(*,*) 'pack_blocks : pe ',mype, & 
     &       ' lcc lfc lec lnc ',lcc,lfc,lec,lnc, & 
     &       ' lguard_in_progress ',lguard_in_progress,' iopt ',iopt & 
     &        ,' ngcell_on_cc ',ngcell_on_cc
      write(*,*) 'pack_blocks : pe ',mype,' loc_message_size(14) ', & 
     &                               loc_message_size(14)
      write(*,*) 'pack_blocks : pe ',mype,' loc_message_size(17) ', & 
     &                               loc_message_size(17)
#endif /* DEBUG */

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
            write(*,*) 'pe ',mype,' sizing send buf to pe ',irpe, & 
     &          ' adding message type ',itype,' size ', & 
     &        loc_message_size(itype), & 
     &        ' accumulated size ',isize, & 
     &' invar ',invar,' message_size_cc ',message_size_cc(itype) & 
     &,' ibndvar ',ibndvar,' message_size_fc ',message_size_fc(itype) & 
     &,' ivaredge ',ivaredge,' message_size_ec ',message_size_ec(itype) & 
     &,' ivarcorn ',ivarcorn,' message_size_nc ',message_size_nc(itype) & 
     &,' offset ',offset
#endif /* DEBUG */
           enddo
           is_buf(1,irpe) = index + 1
           is_buf(2,irpe) = index + isize
           
           index = index + isize
        endif
      enddo

! set up a pointer list to the start address in recv_buffer for each
! block of information in the received messages

      index = 0
      tot_no_blocks_to_be_received = sum(commatrix_recv(:))

      if(allocated(mess_segment_loc)) deallocate(mess_segment_loc)
      allocate(mess_segment_loc(tot_no_blocks_to_be_received))
      mess_segment_loc = 0
      iseg = 0
      lindex = 0
      jrpe =  0

#ifdef DEBUG
      write(*,*) 'pe ',mype,' tot_no_blocks_to_be_received ', & 
     &   tot_no_blocks_to_be_received
#endif /* DEBUG */
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

#ifdef DEBUG
              write(*,*) 'pe ',mype,' sizing recv buf from pe ',irpe, & 
     &          ' adding message type ',itype,' size ', & 
     &          loc_message_size(itype), & 
     &          ' accumulated size ',isize,' iseg ',iseg, & 
     &          ' mess_segment_loc ',mess_segment_loc(iseg), & 
     &          ' lindex ',lindex
              call amr_flush(6)
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
!      S_buffer = 0
#ifdef DEBUG
        write(*,*) 'pe ',mype,' nprocs ',nprocs,' start packing'
#endif /* DEBUG */

      next_pe = 0
      do irpe = 1,nprocs      ! define recv buffer indices
#ifdef DEBUG
        write(*,*) 'pe ',mype,' irpe ',irpe,' commatrix_send ', & 
     &        commatrix_send(irpe)
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
     &     ' from local lb ',lb,' dtype ',dtype,' index ',index & 
     &     ,' buf_dim ',buf_dim
#endif /* DEBUG */
                                  ! pack all arrays for lb into buffer
              call mpi_get_buffer( mype,lb,dtype,iopt,index, & 
     &                             lcc,lfc,lec,lnc, & 
     &                             buf_dim,S_buffer, & 
     &                             nlayersx,nlayersy,nlayersz)

            endif
          enddo
        endif
      enddo

!      write(*,*) 'pe ',mype,' tot_no_blocks_to_be_received ', tot_no_blocks_to_be_received

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
      end subroutine mpi_pack_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine pf_mpi_pack_blocks(mype,nprocs,iopt, & 
     &                           lcc,lfc,lec,lnc, & 
     &                           buf_dim,S_buffer,offset, &
     &                           nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine packs all blocks which are to be sent from mype.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      nprocs         number of processors
!      iopt           option setting for work array
!      lcc            logical switch controlling whether unk data
!                     is packed
!      lfc            logical switch controlling whether facevar data
!                     is packed
!      lec            logical switch controlling whether unk_e_? data
!                     is packed
!      lnc            logical switch controlling whether unk_n data
!                     is packed
!      buf_dim        dimension of buffer
!      S_buffer       send buffer 
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_get_buffer

      implicit none

      include 'mpif.h'

      integer, intent(in)  ::  mype,nprocs,iopt
      logical, intent(in)  ::  lcc,lfc,lec,lnc
      integer, intent(in)  ::  buf_dim,offset
      real,    intent(out) ::  S_buffer(buf_dim)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: loc_message_size(2*27),dtype
      integer :: tot_no_blocks_to_be_received,lindex
      integer :: lb, irpe
      integer :: index 
      integer :: ierrorcode,ierr
      integer :: invar, ibndvar, ivaredge, ivarcorn
      integer :: jrpe, isize, iblk, itype, iseg, next_pe
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! define starting end ending indices for the send and recv buffers

      is_buf = 0
      ir_buf = 0

      l_datapacked = .false.
      if(iopt.gt.1) then
        l_datapacked(1) = .true.
      elseif(iopt.eq.1) then
        if(lcc) l_datapacked(2) = .true.
      endif

      invar = 0
      if (lcc) invar = nvar
      if(lcc.and.lguard_in_progress) invar = ngcell_on_cc

      ibndvar = 0
      ivaredge = 0
      ivarcorn = 0

      if(iopt.eq.1) then
        loc_message_size(:) = invar*message_size_cc(:) & 
     &                     +ibndvar*message_size_fc(:) & 
     &                     +ivaredge*message_size_ec(:) & 
     &                     +ivarcorn*message_size_nc(:) + 3 & 
     &                     + offset
      else
        loc_message_size = message_size_wk + 3 + offset
      endif

      index = 0
      jrpe = 0
      do irpe = 1,nprocs      ! define send buffer indices
        if (commatrix_send(irpe).gt.0) then
          jrpe = jrpe + 1
          isize = 0
          do iblk = 1,commatrix_send(irpe)
            itype = to_be_sent(3,iblk,jrpe)
            isize = isize + loc_message_size(itype)
          enddo
          is_buf(1,irpe) = index + 1
          is_buf(2,irpe) = index + isize
          index = index + isize
        endif
      enddo

! set up a pointer list to the start address in recv_buffer for each
! block of information in the received messages

      index = 0
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
          enddo
          ir_buf(1,irpe) = index + 1
          ir_buf(2,irpe) = index + isize
          index = index + isize
        endif
      enddo

! pack buffer succinctly with all arrays of all blocks that are 
! to be packed for irpe

      index = 1 
      next_pe = 0
      do irpe = 1,nprocs      ! define recv buffer indices
        if (commatrix_send(irpe).gt.0) then
          next_pe = next_pe+1
          do iblk = 1,commatrix_send(irpe)
            if(to_be_sent(1,iblk,next_pe).gt.0) then
              lb = to_be_sent(1,iblk,next_pe)
              dtype = to_be_sent(3,iblk,next_pe)
                                  ! pack all arrays for lb into buffer
              call pf_mpi_get_buffer( mype,lb,dtype,iopt,index, & 
     &                             lcc,lfc,lec,lnc, & 
     &                             buf_dim,S_buffer, & 
     &                             nlayersx,nlayersy,nlayersz)
            endif
          enddo
        endif
      enddo

!      write(*,*) 'pe ',mype,' tot_no_blocks_to_be_received ', tot_no_blocks_to_be_received

      return
      end subroutine pf_mpi_pack_blocks



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mpi_Sbuffer_size(mype,nprocs,iopt, & 
     &                            lcc,lfc,lec,lnc, & 
     &                            buf_dim,offset, & 
     &                            block_sections, fluxes, edges, &
     &                            flux_dir, &
     &                            nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine computes the exact size of the S_buffer
!
!
! Written :     Kevin Olson          January 2007
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      nprocs         number of processors
!      iopt           option setting for work array
!      lcc            logical switch controlling whether unk data
!                     is packed
!      lfc            logical switch controlling whether facevar data
!                     is packed
!      lec            logical switch controlling whether unk_e_? data
!                     is packed
!      lnc            logical switch controlling whether unk_n data
!                     is packed
!      buf_dim        dimension of buffer
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_get_Sbuffer_size, &
     &                                    mpi_get_Sbuffer_size_fluxes, &
     &                                    mpi_get_Sbuffer_size_edges

      implicit none

      include 'mpif.h'

      integer, intent(in)  ::  mype,nprocs,iopt
      logical, intent(in)  ::  lcc,lfc,lec,lnc
      integer, intent(out) ::  buf_dim
      integer, intent(in)  ::  offset
      logical, intent(in)  ::  block_sections, fluxes, edges
      integer, intent(in), optional :: flux_dir
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: loc_message_size(2*27),dtype
      integer :: tot_no_blocks_to_be_received,lindex
      integer :: lb, irpe
      integer :: index 
      integer :: ierrorcode,ierr
      integer :: invar, ibndvar, ivaredge, ivarcorn
      integer :: jrpe, isize, iblk, itype, iseg, next_pe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(iopt.gt.1.and.(lfc.or.lec.or.lnc)) then
         write(*,*) 'Paramesh error : calling mpi_Sbuffer_size with ', & 
     &              'inconsistent argument list - iopt is > 1 while ', & 
     &              'one or more of lfc, lec and lnc are set to true.'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

! count up the number of entries needed in S_buffer

      index = 1 
#ifdef DEBUG
        write(*,*) 'pe ',mype,' nprocs ',nprocs,' start packing'
#endif /* DEBUG */

      next_pe = 0
      do irpe = 1,nprocs      ! define recv buffer indices
#ifdef DEBUG
        write(*,*) 'pe ',mype,' irpe ',irpe,' commatrix_send ', & 
     &        commatrix_send(irpe)
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
     &     ' from local lb ',lb,' dtype ',dtype,' index ',index & 
     &     ,' buf_dim ',buf_dim
#endif /* DEBUG */

              if (block_sections) then

                call mpi_get_Sbuffer_size( mype,lb,dtype,iopt,index, & 
     &                                   lcc,lfc,lec,lnc, & 
     &                                   nlayersx,nlayersy,nlayersz)

              elseif (fluxes) then

                call mpi_get_Sbuffer_size_fluxes( mype,lb,dtype,index, &
     &                                          flux_dir)

              elseif (edges) then

                call mpi_get_Sbuffer_size_edges( mype,lb,dtype,index)

              end if

            endif
          enddo
        endif
      enddo

      if (index > 0) then
         buf_dim = index 
      else
         buf_dim = 1
      end if

      return
      end subroutine mpi_Sbuffer_size




