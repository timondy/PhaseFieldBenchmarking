!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_sort_morton
!! NAME
!!
!!   amr_sort_morton
!!
!! SYNOPSIS
!!
!!   call amr_sort_morton (mort_no, new_loc, nprocs)
!!
!!   call amr_sort_morton (integer, integer, integer)
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
!!   morton order.  This routine assumes that new blocks have been created 
!!   or removed by the refine/derefine process on each local processor AND that the
!!   original list was itself morton ordered.  Hence, blocks only need to be 
!!   sorted locally to create a list which is morton sorted list across all processors.
!!
!!   This routine does not move the block data, but only computes where each
!!   block should move to to acheive a morton order list.
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2004).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_sort_morton (mort_no,new_loc,nprocs)

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
      Integer, intent(inout) :: mort_no(:,:)
      Integer, intent(inout) :: new_loc(:,:)
      Integer, intent(in)    :: nprocs

!-----Local variables and arrays.
      Integer :: itemp
      Integer :: lnblocks2,tot_blocks,no_per_proc,idi,idp
      Integer :: i,j, pid
      Integer :: excess,nprocs_y,nprocs_x,irnkg_s
      Integer :: ierr
      Integer :: lnblocks_left
!      Integer :: ix(2*maxblocks_tr)
!      Integer :: irnkg(2*maxblocks_tr)
!      Integer :: lreflevel(2*maxblocks_tr)
!      Integer :: lreflevel_old(2*maxblocks_tr)
      Integer,allocatable :: ix(:)
      Integer,allocatable :: irnkg(:)
      Integer,allocatable :: lreflevel(:)
      Integer,allocatable :: lreflevel_old(:)
      Integer :: ireduce_datain(1),ireduce_dataout(1)
      Integer :: idatain(1),idataout(1)
      Logical :: lswap
      Logical, save :: first = .True.

!-----Begin executable code.
! CEG allocate memory
  allocate(ix(2*maxblocks_tr))
  allocate(irnkg(2*maxblocks_tr))
  allocate(lreflevel(2*maxblocks_tr))
  allocate(lreflevel_old(2*maxblocks_tr))

!-----sort morton number array mort_no and the local list of refinement levels
!-----and also return the index associated with this permutation
      If (first) Then
         lreflevel(:) = 0
         lreflevel_old(:) = 0
         first = .False.
      End If
      Do i = 1,lnblocks
         ix(i) = i
         lreflevel(i) = lrefine(i)
      End Do

!-----sort top level
      If (lnblocks > 0) Then
        Call morton_sort(mort_no(:,1:lnblocks),ix(1:lnblocks),lnblocks)
        lreflevel_old(:) = lreflevel(:)
        Do i = 1,lnblocks
           lreflevel(i) = lreflevel_old(ix(i))
        End Do
      End If

!-----order segments with same morton number in order of increasing
!-----refinement level
      lswap = .True.
      Do while (lswap)
        lswap = .False.
        Do i = 1,lnblocks-1
          if(mort_no(1,i) == mort_no(1,i+1).And.                       & 
             mort_no(2,i) == mort_no(2,i+1).And.                       & 
             mort_no(3,i) == mort_no(3,i+1).And.                       & 
             mort_no(4,i) == mort_no(4,i+1).And.                       & 
             mort_no(5,i) == mort_no(5,i+1).And.                       & 
             mort_no(6,i) == mort_no(6,i+1).And.                       & 
             lreflevel(i) > lreflevel(i+1) ) Then
            lswap = .True.
            itemp = ix(i)
            ix(i) = ix(i+1)
            ix(i+1) = itemp
            itemp = lreflevel(i)
            lreflevel(i) = lreflevel(i+1)
            lreflevel(i+1) = itemp
          End If
        End Do
      End Do                      

      Do i = 1,lnblocks
         irnkg(ix(i)) = i
      End Do

      lnblocks_left = 0
      idatain(1) = lnblocks
      Call MPI_SCAN (idatain(1),idataout(1),1,                         & 
                     MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,               & 
                     ierr)
      lnblocks_left = idataout(1)
      lnblocks_left = lnblocks_left - lnblocks
      Do j = 1,lnblocks
        irnkg(j) = irnkg(j) + lnblocks_left
      End Do


!-----Compute total list length.
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

!----Compute new_locs from rankings (irnkg) returned by amr_bi_sort.
!----The following divides blocks evenly among processors without regard to
!----work.

      call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr);

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

!         if (pid.eq.1) print *,i,idi,idp
         
      End Do

!      call ceg_block_balancing(pid, nprocs, new_loc)

! CEG free memory
  deallocate(ix)
  deallocate(irnkg)
  deallocate(lreflevel)
  deallocate(lreflevel_old)


      Return
      End Subroutine amr_sort_morton
