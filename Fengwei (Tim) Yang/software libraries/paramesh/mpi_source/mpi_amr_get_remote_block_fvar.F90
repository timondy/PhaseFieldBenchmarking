!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_get_remote_block_fvar
!! NAME
!!
!!   mpi_amr_get_remote_block_fvar
!!
!! SYNOPSIS
!!
!!   call mpi_amr_get_remote_block_fvar (mype, remote_pe, remote_block, icoord,
!!                                       recvx, recvy, recvz, idest)
!!
!!   call mpi_amr_get_remote_block_fvar (integer, integer, integer, integer,
!!                                       real, real, real, integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype             
!!     The local processor
!!
!!   integer, intent(in) :: remote_pe        
!!     The remote processor.
!!
!!   integer, intent(in) :: remote_block     
!!     The local block id of the block to be copied from
!!     the remote processor.
!!    
!!   integer, intent(in) :: icoord           
!!     Coordinate to fetch, ie facevarx, facveary or facevarz.
!!
!!   real, intent(out) :: recvx            
!!     Output array if icoord = 1.
!!  
!!   real, intent(out) :: recvy
!!     Output array if icoord = 2.
!!
!!   real, intent(out) :: recvz
!!     Output array if icoord = 3.
!!
!!   integer, intent(in) :: idest            
!!     Selects the storage space in the 1blk data structures which is to
!!     be used in this call. If the leaf node is having its
!!     guardcells filled then set this to 1, if its parent
!!     is being filled set it to 2.
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
!!   workspace
!!   mpi_morton
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_mpi_find_blk_in_buffer
!!   mpi_set_message_limits
!!
!! RETURNS
!!
!!   Upon return the data from 'remote_block' on 'remote_pe' is placed locally on the
!!   calling processor in the arrays 'recvx', 'recvy', or 'recvz'.
!!
!! DESCRIPTION
!! 
!!  This routine copies guard cell information to face iface in layer
!!  idest of the working block, from the appropriate face of the neighboring 
!!  block, assuming that the neighboring block is on a different processor.
!! 
!! AUTHORS
!!
!! Written by Peter MacNeice (August 2001).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_amr_get_remote_block_fvar(mype,                   & 
                           remote_pe,remote_block,icoord,              & 
                           recvx,recvy,recvz,idest)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use mpi_morton
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in) :: mype,remote_pe,remote_block
      Integer, Intent(in) :: icoord,idest
      Real, Intent(out)   :: recvx(:,:,:,:)
      Real, Intent(out)   :: recvy(:,:,:,:)
      Real, Intent(out)   :: recvz(:,:,:,:)

!-----Local arrays and variables.
      Integer :: nguard0
      Integer :: nguard_work0
      Integer :: ierrorcode,ierr
      Integer :: dtype
      Integer :: vtype
      Integer :: index,index0
      Integer :: i, j, k, ia, ib, ja, jb, ka, kb
      Logical :: lfound,lerror

!-----Begin executable code.
      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      If (remote_block <= lnblocks.And.remote_pe == mype) Then

!-------Copy complete remote block into a buffer block Called recv.
        If (no_permanent_guardcells) Then
        If (icoord == 1) Then
          recvx(:,:,:,:) = gt_facevarx(:,:,:,:,remote_block)
        ElseIf (icoord == 2) Then
          If (ndim >= 2) Then
             recvy(:,:,:,:) = gt_facevary(:,:,:,:,remote_block)
          End If
        ElseIf (icoord == 3) Then
          If (ndim == 3) Then
             recvz(:,:,:,:) = gt_facevarz(:,:,:,:,remote_block)
          End If
        End If
        Else ! no_permanent_guardcells
        If (icoord == 1) Then
          recvx(:,:,:,:) = facevarx(:,:,:,:,remote_block)
        ElseIf (icoord == 2) Then
          If (ndim >= 2) Then
             recvy(:,:,:,:) = facevary(:,:,:,:,remote_block)
          End If
        ElseIf (icoord == 3) Then
          If (ndim == 3) Then
             recvz(:,:,:,:) = facevarz(:,:,:,:,remote_block)
          End If
        End If  

        End If  ! End If (no_permanent_guardcells)

      Else

        Call amr_mpi_find_blk_in_buffer(mype,remote_block,             & 
                                        remote_pe,idest,dtype,index0,  &
                                        lfound)

        lerror=.False.
        If (.Not.lfound) lerror = .True.
        If ( icoord == 1 .And.                                         & 
            (dtype.ne.13.And.dtype.ne.14.And.dtype.ne.15.And.          & 
             dtype.ne.13+27.And.dtype.ne.14+27.And.dtype.ne.15+27)     & 
            ) lerror = .True.
        If ( icoord == 2 .And.                                         & 
            (dtype.ne.11.And.dtype.ne.14.And.dtype.ne.17.And.          & 
             dtype.ne.11+27.And.dtype.ne.14+27.And.dtype.ne.17+27)     & 
            ) lerror = .True.
        If ( icoord == 3 .And.                                         & 
            (dtype.ne.5.And.dtype.ne.14.And.dtype.ne.23.And.           & 
             dtype.ne.5+27.And.dtype.ne.14+27.And.dtype.ne.23+27)      & 
            ) lerror = .True.

        If (lerror) Then
          write(*,*) 'Paramesh error : pe ',mype,                      & 
            ' needed remote blk facevar data ',                        & 
            remote_block,remote_pe,' but could not find it or only ',  & 
            ' found part of it in the message buffer.',                & 
            '  Contact developers for help.',                          & 
            ' lfound ',lfound,' dtype ',dtype,' icoord ',icoord
          Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        End If

!-------starting index if cell-centered data is also included in recv_buf
        index = index0
        If (l_datapacked(2)) index =                                   & 
                             index + nvar*message_size_cc(dtype)

        vtype = 2
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype)

        If (icoord == 1) Then
          Do k = ka,kb
          Do j = ja,jb
          Do i = ia,ib
            recvx(1:nbndvar,i,j,k) =                                   & 
                    temprecv_buf(index+1:index+nbndvar)
            index = index+nbndvar
          End Do
          End Do
          End Do
        Else
          index  = index + nbndvar*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
        End If

        If (ndim >= 2) Then
          vtype = 3
          Call mpi_set_message_limits(                                 & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype)

          If (icoord == 2) Then
            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              recvy(1:nbndvar,i,j,k) =                                 & 
                    temprecv_buf(index+1:index+nbndvar)
              index = index+nbndvar
            End Do
            End Do
            End Do
          Else
            index  = index + nbndvar*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
          End If

        End If

        If (ndim == 3) Then
         vtype = 4
         Call mpi_set_message_limits(                                  & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype)

         If (icoord == 3) Then
           Do k = ka,kb
           Do j = ja,jb
           Do i = ia,ib
             recvz(1:nbndvar,i,j,k) =                                  & 
                    temprecv_buf(index+1:index+nbndvar)
             index = index+nbndvar
           End Do
           End Do
           End Do
         End If
        Else
         index  = index + nbndvar*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
        End If

      End If  ! End If (remote_block <= lnblocks.And.remote_pe == mype)

      Return
      End Subroutine mpi_amr_get_remote_block_fvar
