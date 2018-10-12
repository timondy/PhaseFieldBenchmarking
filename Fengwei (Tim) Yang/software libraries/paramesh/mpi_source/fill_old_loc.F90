!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/fill_old_loc
!! NAME
!!
!!   fill_old_loc
!!
!! SYNOPSIS
!!
!!   call fill_old_loc(new_loc,old_loc,nprocs,mype)
!!   call fill_old_loc(integer array, integer array, integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(inout) :: new_loc(:,:)
!!     new locations (processor and location in morton list) 
!!   integer, intent(out) :: old_loc(:,:)
!!     old locations (processor and location in morton list) 
!!   integer, intent(in) :: nprocs, mype
!!     number of procs. and proc. id.
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
!!   io
!!
!! CALLS
!! 
!! RETURNS
!!
!!   Returns a set of locations in the array new_loc(:,:) which are linked
!!   with the 'old' locations in old_loc.  
!!
!! DESCRIPTION
!!
!!   This routine sorts morton numbers stored locally on a single processor.
!!   Here, the blocks are moving from their old locations (stored a proc., 
!!   list id pair) to new locations.
!!   A block can be sent from old_loc (by posting am MPI_SEND) to new_loc (by
!!   posting an MPI_RECV).
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2001).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine fill_old_loc(new_loc,old_loc,nprocs,mype)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      integer, intent(inout) :: new_loc(:,:)
      integer, intent(out)   :: old_loc(:,:)
      integer, intent(in)    :: nprocs,mype

!-----Local variables and arrays.
      Integer :: nrecv,nsend
      Integer :: nsend_to_proc(0:nprocs), nrecv_pack(0:nprocs)
      Integer :: ierr
      Integer :: i
!      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
!      Integer :: reqr(maxblocks_tr)
!      Integer :: kk(maxblocks_tr)
      Integer,allocatable :: statr(:,:)
      Integer,allocatable :: reqr(:)
      Integer,allocatable :: kk(:)
      Integer :: idummy_array(1)

!-----Begin executable code.
! CEG allocate memory
  allocate(statr(MPI_STATUS_SIZE,maxblocks_tr))
  allocate(reqr(maxblocks_tr))
  allocate(kk(maxblocks_tr))

!-----fill `old_loc' (pointer from new block location back to
!-----its old, unsorted location)

!-----count no. of receives to post
!-----count no. of sends on each proc. to all other procs
      nsend_to_proc(:) = 0
      Do i = 1,maxblocks_tr
        If (new_loc(1,i) > 0) Then
        If (new_loc(2,i).ne.mype) Then ! its a send
          nsend_to_proc(new_loc(2,i)) =                                & 
           nsend_to_proc(new_loc(2,i)) + 1
        End If
        End If
      End Do

!-----collect data for `this' proc from other procs so that
!-----the total no. of receives to post can be computed
!-----(Changed by M. Zingale and J. Dursi)
      nrecv = 0
      Call MPI_ALLREDUCE(nsend_to_proc, nrecv_pack, nprocs,            & 
                         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      nrecv = nrecv_pack(mype)

      old_loc(:,:) = -1
      Do i = 1,nrecv
        Call MPI_IRECV(kk(i),                                          & 
                       1,                                              & 
                       MPI_INTEGER,                                    & 
                       MPI_ANY_SOURCE,                                 & 
                       MPI_ANY_TAG,                                    & 
                       MPI_COMM_WORLD,                                 & 
                       reqr(i),                                        & 
                       ierr)
      End Do

      nsend = 0
      Do i = 1,maxblocks_tr
        If (new_loc(1,i) > 0) Then
        If (new_loc(2,i).ne.mype) Then
          nsend = nsend + 1
          idummy_array(1) = i
          Call MPI_SSEND(idummy_array(1),                              & 
                         1,                                            & 
                         MPI_INTEGER,                                  & 
                         new_loc(2,i),                                 & 
                         new_loc(1,i),                                 & 
                         MPI_COMM_WORLD,                               & 
                         ierr)
        Else
          old_loc(1,new_loc(1,i)) = i
          old_loc(2,new_loc(1,i)) = mype
        End If
        End If
      End Do

      If (nrecv > 0) Then
         Call MPI_WAITALL (nrecv, reqr, statr, ierr)
         Do i = 1,nrecv
           old_loc(1,statr(MPI_TAG,i)) = kk(i)
           old_loc(2,statr(MPI_TAG,i)) = statr(MPI_SOURCE,i)
         End Do
      End If

      Call MPI_BARRIER (MPI_COMM_WORLD,ierr)

! CEG deallocate memory
  deallocate(statr)
  deallocate(reqr)
  deallocate(kk)

      Return
      End Subroutine fill_old_loc


