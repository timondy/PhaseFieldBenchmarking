!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/amr_sort_morton_reorder_grid
!! NAME
!!
!!   amr_sort_morton_reorder_grid
!!
!! SYNOPSIS
!!
!!   call amr_sort_morton_reorder_grid(mort_no, new_loc, nprocs)
!!
!!   call amr_sort_morton_reorder_grid(integer, integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(out) ::  mort_no(:,:)
!!     This is a list of the morton numbers computed and returned by this subroutine.  Note
!!     that this is a 2 dimensional array.  For high levels of refinement, more bits are required
!!     to store the morton numbers, we spread these bit patterns over several 32 (or 64) bit
!!     integers.  The first dimension represents the number of these integers which are used and
!!     the 2nd dimension is the length of the block list.
!!
!!   integer, intent(inout) :: new_loc(:,:)
!!     The new locations (processor and local id) that the block is to migrate to during
!!     during reordering.
!!
!!   integer, intent(in)    :: nprocs
!!     The number of processors.
!!
!! INCLUDES
!!
!!   mpif.h
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   io
!!   paramesh_mpi_interfaces
!!
!! CALLS
!! 
!!   morton_sort
!!    
!! RETURNS
!!
!!   Returns the computed morton number for each block.
!!
!! DESCRIPTION
!!
!!   Given a list of morton numbers for blocks this subroutine sorts them into 
!!   morton order.  This particular version of this routine is used when reorder of the
!!   grid is requested and it is assumed that the blocks are not in any particular 
!!   order.  It does this by collecting the lists of blocks from other processors
!!   to processor 0 (zero) and doing a sort there.  The results are Then transfered
!!   back to the other processors.
!!
!!   This routine does not move the block data, but only computes where each
!!   block should move to to acheive a morton order list.
!!
!! AUTHORS
!!
!!   Kevin Olson (2005).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_sort_morton_reorder_grid (mort_no,new_loc,nprocs)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_interfaces, only : morton_sort

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(inout) :: mort_no(:,:)
      Integer, Intent(inout) :: new_loc(:,:)
      Integer, Intent(in)    :: nprocs

!-----Local variables and arrays.
      Integer :: statr(MPI_STATUS_SIZE),reqr
      Integer :: lnblocks_per_proc(0:nprocs-1)
      Integer :: itemp
      Integer :: lnblocks2,tot_blocks,no_per_proc,idi,idp
      Integer :: i,j
      Integer :: excess,nprocs_y,nprocs_x,irnkg_s
      Integer :: ierr
      Integer :: lnblocks_left
!      Integer :: lreflevel(2*maxblocks_tr)
      Integer,allocatable :: lreflevel(:)
      Integer :: ireduce_datain(1),ireduce_dataout(1)
      Integer :: idatain(1),idataout(1)
      Integer :: iproc, mype, nsend, nrecv
      Integer, Allocatable :: mort_no_across(:,:)
      Integer, Allocatable :: lreflevel_across(:)
      Integer, Allocatable :: lreflevel_old(:)
      Integer, Allocatable :: ix_across(:)
      Integer, Allocatable :: irnkg(:)
      Logical :: lswap
      Logical, save :: first = .True.

!-----Begin executable code.
! CEG allocate memory
  allocate(lreflevel(2*maxblocks_tr))

      Call MPI_COMM_RANK (MPI_COMM_WORLD, mype, ierr)

!-----compute total no. of blocks across all processors
      ireduce_datain(1) = lnblocks
      Call MPI_ALLREDUCE (ireduce_datain,ireduce_dataout,              & 
                          1,MPI_INTEGER,                               & 
                          MPI_SUM,MPI_COMM_WORLD,ierr)
      tot_blocks = ireduce_dataout(1)

!-----allocate an array to hold all processor numbers
      Allocate(mort_no_across(size(mort_no,1),tot_blocks))
      Allocate(lreflevel_across(tot_blocks))
      Allocate(lreflevel_old(tot_blocks))
      Allocate(ix_across(tot_blocks))
      Allocate(irnkg(tot_blocks))

!-----sort morton number array mort_no and the local list of refinement levels
!-----and also return the index associated with this permutation
      lreflevel(:) = 0
      lreflevel_old(:) = 0

      Do i = 1,lnblocks
         lreflevel(i) = lrefine(i)
      End Do
      Do i = 1, tot_blocks
         ix_across(i) = i
      End Do

!----Stuff all the morton numbers into one array on processor 0
!----First collect the number of blocks on each process to processor 0.

      Call MPI_GATHER(lnblocks,          1, MPI_INTEGER,               & 
                      lnblocks_per_proc(0), 1, MPI_INTEGER,            & 
                      0, MPI_COMM_WORLD, ierr)

!-----Proc. 0 recvs from all other processes.
      nrecv = 1
      Do iproc = 0, nprocs - 1
         If (mype == 0) Then

            If (lnblocks_per_proc(iproc) > 0) Then
            If (iproc > 0) Then

               Call MPI_IRECV(                                         & 
                    mort_no_across(1,nrecv),                           & 
                    size(mort_no,1)*lnblocks_per_proc(iproc),          & 
                    MPI_INTEGER,                                       & 
                    iproc,1,                                           & 
                    MPI_COMM_WORLD,                                    & 
                    reqr,ierr)
               Call MPI_WAIT(reqr,statr,ierr)

               Call MPI_IRECV(                                         & 
                    lreflevel_across(nrecv),                           & 
                    lnblocks_per_proc(iproc),                          & 
                    MPI_INTEGER,                                       & 
                    iproc,2,                                           & 
                    MPI_COMM_WORLD,                                    & 
                    reqr,ierr)
               Call MPI_WAIT(reqr,statr,ierr)

               nrecv = nrecv + lnblocks_per_proc(iproc)

            Else

               mort_no_across(:,1:lnblocks) = mort_no(:,1:lnblocks)
               lreflevel_across(1:lnblocks) = lreflevel(1:lnblocks)

               nrecv = nrecv + lnblocks_per_proc(iproc)

            End If ! End If (iproc > 0)
            End If ! End If (lnblocks_per_proc(iproc) > 0)

         Else

            If (lnblocks > 0 .and. iproc == mype) Then

               Call MPI_SSEND(mort_no(:,1:lnblocks),                   & 
                    size(mort_no,1)*lnblocks,                          & 
                    MPI_INTEGER,                                       & 
                    0,1,                                               & 
                    MPI_COMM_WORLD,ierr)
               
               Call MPI_SSEND(lreflevel(1:lnblocks),lnblocks,          & 
                    MPI_INTEGER,                                       & 
                    0,2,                                               & 
                    MPI_COMM_WORLD,ierr)
            End If

         End If  ! End If (mype == 0)

      End Do  ! End Do iproc = 0, nprocs - 1

!-----sort top level
      If (mype == 0) Then
         If (tot_blocks > 0) Then

            Call morton_sort(mort_no_across(:,1:tot_blocks),               & 
                         ix_across(1:tot_blocks),                      & 
                         tot_blocks)

            lreflevel_old(:) = lreflevel_across(:)
            Do i = 1,tot_blocks
               lreflevel_across(i) = lreflevel_old(ix_across(i))
            End Do

         End If  ! End If (tot_blocks > 0)
      End If  ! End If (mype == 0)

!-----order segments with same morton number in order of increasing
!-----refinement level
      If (mype == 0) Then
      lswap = .True.
      Do While (lswap)
        lswap = .False.
        Do i = 1,tot_blocks-1
          If (mort_no_across(1,i) == mort_no_across(1,i+1).and.        & 
              mort_no_across(2,i) == mort_no_across(2,i+1).and.        & 
              mort_no_across(3,i) == mort_no_across(3,i+1).and.        & 
              mort_no_across(4,i) == mort_no_across(4,i+1).and.        & 
              mort_no_across(5,i) == mort_no_across(5,i+1).and.        & 
              mort_no_across(6,i) == mort_no_across(6,i+1).and.        & 
              lreflevel_across(i) > lreflevel_across(i+1) ) Then
            lswap = .True.
            itemp = ix_across(i)
            ix_across(i) = ix_across(i+1)
            ix_across(i+1) = itemp
            itemp = lreflevel_across(i)
            lreflevel_across(i) = lreflevel_across(i+1)
            lreflevel_across(i+1) = itemp
          End If
        End Do
      End Do  

      Do i = 1,tot_blocks
         irnkg(ix_across(i)) = i
      End Do

      End If  ! End If (mype == 0)

!-----Now we need to send the array irnkg back to their original processor
!-----Proc. 0 recvs from all other processes.
      nsend = 1
      Do iproc = 0, nprocs - 1

         If (mype > 0) Then

            If (lnblocks > 0 .and. iproc == mype) Then
               Call MPI_IRECV(                                         & 
                    irnkg(1),                                          & 
                    lnblocks,                                          & 
                    MPI_INTEGER,                                       & 
                    0,3,                                               & 
                    MPI_COMM_WORLD,                                    & 
                    reqr,ierr)
               Call MPI_WAIT(reqr,statr,ierr)
            End If 

         Else

            If (lnblocks_per_proc(iproc) > 0) Then
               If (iproc > 0) Then
                  Call MPI_int_SSEND(                                  &  
                       irnkg(nsend),                                   & 
                       lnblocks_per_proc(iproc),                       & 
                       MPI_INTEGER,                                    & 
                       iproc,3,                                        & 
                       MPI_COMM_WORLD,ierr)
               Else
                  irnkg(1:lnblocks) = irnkg(1:lnblocks)
               End If
               nsend = nsend + lnblocks_per_proc(iproc) 
               
            End If

         End If  ! End If (mype > 0)

      End Do  ! End Do iproc = 0, nprocs - 1

!-----Compute total list length.
!-----I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save 
!-----statement.
      lnblocks2 = lnblocks 
      ireduce_datain(1) = lnblocks2
      Call MPI_ALLREDUCE (ireduce_datain,ireduce_dataout,              & 
                          1,MPI_INTEGER,                               & 
                          MPI_SUM,MPI_COMM_WORLD,ierr)
      tot_blocks = ireduce_dataout(1)

      no_per_proc = tot_blocks/nprocs

      excess = tot_blocks - no_per_proc*nprocs
      nprocs_y = (no_per_proc+1)*nprocs - tot_blocks
!-----no. of processors which will get no_per_proc + 1 blocks
      nprocs_x = nprocs - nprocs_y
!-----rank in list which divides those which go on processor with one number
!-----of blocks from those which go on another set of blocks w. a different
!-----no. of blocks
      irnkg_s = nprocs_x*(no_per_proc+1)

!-----Compute new_locs from rankings (irnkg) returned by amr_bi_sort.
!-----The following divides blocks evenly among processors without regard to
!-----work.

      Do i = 1,lnblocks

         idp = (irnkg(i)-1)/(no_per_proc+1) ! processor to send to
         If (irnkg(i) <= irnkg_s) Then
            idi = mod((irnkg(i)-1),no_per_proc+1) + 1 ! rank inside 
                                                      ! local array
                                                      ! to write to
         Else
            idp = (irnkg(i)-irnkg_s-1)/(no_per_proc) ! processor to send to
            idp = idp + nprocs_x
            idi = mod((irnkg(i)-irnkg_s-1),no_per_proc) + 1 ! rank inside 
                                                            ! local array
                                                            ! to write to
         End If

         new_loc(1,i) = idi
         new_loc(2,i) = idp
         
      End Do

      Deallocate(mort_no_across)
      Deallocate(lreflevel_across)
      Deallocate(lreflevel_old)
      Deallocate(ix_across)
      Deallocate(irnkg)

!CEG free memory
  deallocate(lreflevel)


      Return
      End Subroutine amr_sort_morton_reorder_grid

