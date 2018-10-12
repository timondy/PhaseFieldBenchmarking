!!!#define DEBUG

!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_mpi_find_blk_in_buffer
!! NAME
!!
!!   amr_mpi_find_blk_in_buffer
!!
!! SYNOPSIS
!!
!!    Call amr_mpi_find_blk_in_buffer(mype,remote_block,remote_pe,idest,
!!                                    dtype,index0,lfound)
!!    Call amr_mpi_find_blk_in_buffer(integer,integer,integer,integer,
!!                                    integer,integer,logical)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in)  :: mype           local processor
!!   Integer, Intent(in)  :: remote_block   remote block 
!!   Integer, Intent(in)  :: remote_pe      remote processor
!!   Integer, Intent(in)  :: idest
!!   Integer, Intent(out) :: dtype          message type - an integer between 1 and 27
!!                                          indicating the section of a blck contained in
!!                                          the message segment from (remote_block,remote_pe).
!!   Integer, Intent(out) :: index0         the address index0+1 is where you should start
!!                                           reading the data from this message.
!!   Logical, Intent(out) :: lfound         logical indicating if the block is found or not
!!    
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   mpi_morton
!!
!! CALLS
!!
!!   amr_abort
!!
!! RETURNS
!!
!!   Returns message type in dtype.
!!   Returns starting address in buffer where to find the requested block 
!!    in the variable index0+1.
!!   Returns if the block is found or not in lfound.
!!
!! DESCRIPTION
!!
!!   This routine finds where data for a remote block with address
!!   (remote_block,remote_pe) is in the recv buffer. It returns
!!   the message type, dtype, and the address in recv of the data word
!!   preceeding the first word of real data.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          May 2001
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_mpi_find_blk_in_buffer(                           & 
              mype,remote_block,remote_pe,idest,dtype,index0,lfound)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      implicit none

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in)  :: mype,remote_pe,remote_block,idest
      Integer, Intent(out) :: dtype,index0
      Logical, Intent(out) :: lfound

!-----Local arrays and variables.
      integer :: jseg,seg_no,iaddress,no_of_comms,jpe,jpe0
      integer :: ierrorcode,ierr,no_of_segments
      integer :: rem_pe,rem_blk,seg_offset
      integer :: iseg_no,jj
      logical :: llfound

!-----Begin executable code.

      If (remote_pe.ne.mype) Then

        rem_blk = remote_block
        rem_pe  = remote_pe

        llfound = .False.
        jj = ladd_strt(rem_pe)
        iseg_no = 0
        Do While(.Not.llfound.and.jj <= ladd_end(rem_pe))
          If (rem_blk == laddress(1,jj).and.                           & 
             rem_pe == laddress(2,jj)) Then
            llfound = .True.
            iseg_no = jj - strt_buffer + 1
          Else
            jj = jj+1
          End If
        End Do

      ElseIf (remote_pe == mype.and.remote_block > lnblocks) Then

        rem_blk = laddress(1,remote_block)
        rem_pe  = laddress(2,remote_block)
        iseg_no = remote_block - strt_buffer + 1
        llfound = .True.

      End If  ! End If (remote_pe.ne.mype)

      If (rem_pe.ne.mype) Then

#ifdef DEBUG
!-------locate rem_pe in the list of sending processors
        no_of_comms = size(pe_source)
        lfound = .False.
        jpe  = 0
        jpe0 = 0
        Do While((.Not.lfound).and.(jpe0 < no_of_comms))
          jpe0 = jpe0+1
          If (rem_pe == pe_source(jpe0)-1) Then
            lfound = .True.
            jpe = jpe0
          End If
        End Do
!-------If rem_pe is not located stop with error message
        If (jpe == 0) Then
          If (idest == 2) return
          Write(*,*) 'Paramesh error : pe ',mype,                      &  
           ' pe address of required data is not in the list of ',      & 
           'communicating pes. ',                                      & 
           ' remote_block ',remote_block,                              & 
           ' remote_pe ',remote_pe,                                    & 
           ' rem_pe ',rem_pe,                                          & 
           ' laddress ',laddress
          Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        End If
        no_of_segments = size(to_be_received,2)
        lfound = .False.
        jseg = 0
        seg_no = jseg
        Do While((.Not.lfound).and.(jseg < no_of_segments))
          jseg = jseg+1
          If (to_be_received(1,jseg,jpe) == rem_blk .and.              & 
             to_be_received(2,jseg,jpe) == rem_pe+1 ) Then
            lfound = .True.
            seg_no = jseg
          End If
        End Do
!-------Now compute where the list of segments from proc rem_pe actually begins
!-------in the complete list of message segments received on this processor
        seg_offset = 0
        If (rem_pe > 0) seg_offset = sum(commatrix_recv(1:rem_pe))
        seg_no = seg_no + seg_offset
#ifdef DEBUGX
        If (seg_no.ne.iseg_no) Then
          Write(*,*) 'seg_no and iseg_no are different ',              & 
                      seg_no,iseg_no
          Call amr_abort()
        End If
#endif /* DEBUGX */
#endif /* DEBUG */

        seg_no = iseg_no
        lfound = llfound

!-------If the requested segment is not located stop with error message
        If (seg_no == 0) Then
          If (idest == 2) return
#ifdef DEBUG
          Write(*,*) 'Paramesh error : ',                              & 
           'message segment required is not in the list of ',          & 
           'segments received.: proc ',mype,' to_be_received ',        & 
           to_be_received(:,:,jpe),' mpi_pattern_id ',mpi_pattern_id
#else
          Write(*,*) 'Paramesh error : ',                              & 
           'message segment required is not in the list of ',          & 
           'segments received.'
#endif
          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
        End If  ! End If (seg_no == 0)

!-------set start address for this segment in R_buffer
        iaddress = mess_segment_loc(seg_no)

!-------Read out message into appropriate part of unk1 or work1
        dtype = anint(temprecv_buf(iaddress+2))

!-------We must Write the message into recv first in case it is larger
!-------than the range we actually need.

        index0 = iaddress+2

!-------NOTE: start address of data is now index0+1

      End If  ! If (rem_pe.ne.mype)

      Return
      End Subroutine amr_mpi_find_blk_in_buffer
