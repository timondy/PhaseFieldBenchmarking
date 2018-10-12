!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_morton_bnd_restrict
!! NAME
!!
!!   mpi_morton_bnd_restrict
!!
!! SYNOPSIS
!!
!!   Call mpi_morton_bnd_restrict(mype, nprocs, tag_offset)
!!   Call mpi_morton_bnd_restrict(integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in)    :: mype       Local processor id.
!!   Integer, Intent(in)    :: nprocs     Number of processors.
!!   Integer, Intent(inout) :: tag_offset A unique id used in marking messages.
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
!!   timings
!!   mpi_morton
!!   constants
!!
!! CALLS
!!
!!    mpi_amr_write_restrict_comm
!!    process_fetch_list
!!
!! RETURNS
!!
!!    Does not return anything.
!!
!! DESCRIPTION
!!
!!   This routine does a communications analysis for data restriction and
!!   constructs and stores lists of off-processor blocks which are to be 
!!   communicated during restriction.  Also stored are which sections of 
!!   the blocks to be fetched.
!!
!! AUTHORS
!!
!!    Written by Peter MacNeice  and Michael Gehmeyr, February 2000.
!!    Major simplification and rewrite by Kevin Olson, August 2007.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_morton_bnd_restrict (mype,                        &
                                          nprocs,                      &
                                          tag_offset)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use mpi_morton
      Use constants

      Use paramesh_mpi_interfaces, only : mpi_amr_write_restrict_comm, & 
                                          process_fetch_list

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Variables
      Integer, Intent(in)    ::  mype,nprocs
      Integer, Intent(inout) ::  tag_offset

!-----Local variables
      Integer :: lb,i,j,k,j00
      Integer :: ierror
      Integer :: istack
      Integer :: iproc
      Integer :: npts_neigh1,npts_neigh2, sendmsg
      Integer,Dimension (:),  Allocatable :: n_to_left
      Integer,Dimension (:,:),Allocatable :: fetch_list
      Integer,Dimension (:,:),Allocatable :: tfetch_list

!-----Begin executable code.
      npts_neigh1 = npts_neigh
      npts_neigh2 = npts_neigh+100
      Allocate(fetch_list(3,npts_neigh2))
      Allocate(n_to_left(0:nprocs-1))

!----COMPUTE the number of blocks to the 'left' (ie. stored on processors with
!----smaller process ids) of every other processor
      n_to_left(mype) = lnblocks
      sendmsg = n_to_left(mype)
      Call MPI_ALLGATHER(sendmsg,                                      &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         n_to_left,                                    &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         MPI_COMM_WORLD,                               &
                         ierror)
                        
      Do iproc = nprocs-1, 1, -1
         n_to_left(iproc) = n_to_left(iproc-1)
      End Do
      n_to_left(iproc) = 0
      Do iproc = 2, nprocs-1
         n_to_left(iproc) = n_to_left(iproc) + n_to_left(iproc-1)
      End Do

!-----Initializations
      commatrix_send = 0
      commatrix_recv = 0

!-----Construct a list of potential neighbors of all blocks on this
!-----processor, and potential neighbors of their parents.
!-----Exclude any which are on this processor.

      istack = 0
      Do lb = 1, lnblocks

      If (nodetype(lb) == 2 .or.                                       &
          (advance_all_levels .and. nodetype(lb) == 3)) Then

!------ADD OFF PROCESSOR CHILDREN OF BLOCK 'lb' TO FETCH LIST
        Do i = 1,nchild
        If (child(1,i,lb) > 0 .and.                            & 
            child(2,i,lb) .ne. mype) Then

            istack = istack + 1
            If (istack > npts_neigh1) Call expand_fetch_list
            fetch_list(1,istack) = child(1,i,lb)
            fetch_list(2,istack) = child(2,i,lb)
!-----------Fetch entire block
            fetch_list(3,istack) = 14

         End If  ! End If child
         End Do  ! End Do i = 1,nchild

      End If  ! End If (nodetype(lb) <= 2 .or. advance_all_levels)

      End Do  ! End Do lb = 1, lnblocks

      Call process_fetch_list(fetch_list,                              &
                              istack,                                  &
                              mype,                                    &
                              nprocs,                                  &
                              n_to_left,                               &
                              tag_offset)

!------Mark morton data up to date
       morton_limits_set = .True.

!------Store communication info for future use
       Call mpi_amr_write_restrict_comm(nprocs)

!------Deallocate any memory which was dynamically allocated for local 
!------use in this routine.
       If (Allocated(fetch_list)) deallocate(fetch_list)
       If (Allocated(n_to_left)) deallocate(n_to_left)

      Return

      Contains
        Subroutine expand_fetch_list

               If (Allocated(tfetch_list)) deallocate(tfetch_list)
               Allocate(tfetch_list(3,npts_neigh2))
               tfetch_list(:,:istack-1) = fetch_list(:,:istack-1)
               npts_neigh1 = npts_neigh1 + 3000
               npts_neigh2 = npts_neigh2 + 3000
               deallocate(fetch_list)
               Allocate(fetch_list(3,npts_neigh2))
               fetch_list(:,:istack-1) = tfetch_list(:,:istack-1)
               Deallocate(tfetch_list)

        End Subroutine expand_fetch_list
      End Subroutine mpi_morton_bnd_restrict


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

      Subroutine pf_morton_bnd_restrict (mype,                        &
                                          nprocs,                      &
                                          tag_offset)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use mpi_morton
      Use constants

      Use paramesh_mpi_interfaces, only : mpi_amr_write_restrict_comm, & 
                                          process_fetch_list

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Variables
      Integer, Intent(in)    ::  mype,nprocs
      Integer, Intent(inout) ::  tag_offset

!-----Local variables
      Integer :: lb,i,j,k,j00
      Integer :: ierror
      Integer :: istack
      Integer :: iproc
      Integer :: npts_neigh1,npts_neigh2, sendmsg
      Integer,Dimension (:),  Allocatable :: n_to_left
      Integer,Dimension (:,:),Allocatable :: fetch_list
      Integer,Dimension (:,:),Allocatable :: tfetch_list

!-----Begin executable code.
      npts_neigh1 = npts_neigh
      npts_neigh2 = npts_neigh+100
      Allocate(fetch_list(3,npts_neigh2))
      Allocate(n_to_left(0:nprocs-1))

!----COMPUTE the number of blocks to the 'left' (ie. stored on processors with
!----smaller process ids) of every other processor
      n_to_left(mype) = lnblocks
      sendmsg = n_to_left(mype)
      Call MPI_ALLGATHER(sendmsg,                                      &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         n_to_left,                                    &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         MPI_COMM_WORLD,                               &
                         ierror)
                        
      Do iproc = nprocs-1, 1, -1
         n_to_left(iproc) = n_to_left(iproc-1)
      End Do
      n_to_left(iproc) = 0
      Do iproc = 2, nprocs-1
         n_to_left(iproc) = n_to_left(iproc) + n_to_left(iproc-1)
      End Do

!-----Initializations
      commatrix_send = 0
      commatrix_recv = 0

!-----Construct a list of potential neighbors of all blocks on this
!-----processor, and potential neighbors of their parents.
!-----Exclude any which are on this processor.

      istack = 0
      Do lb = 1, lnblocks

      If (nodetype(lb) == 2 .or.                                       &
          (advance_all_levels .and. nodetype(lb) == 3)) Then

!------ADD OFF PROCESSOR CHILDREN OF BLOCK 'lb' TO FETCH LIST
        Do i = 1,nchild
        If (child(1,i,lb) > 0 .and.                            & 
            child(2,i,lb) .ne. mype) Then

            istack = istack + 1
            If (istack > npts_neigh1) Call expand_fetch_list
            fetch_list(1,istack) = child(1,i,lb)
            fetch_list(2,istack) = child(2,i,lb)
!-----------Fetch entire block
            fetch_list(3,istack) = 14

         End If  ! End If child
         End Do  ! End Do i = 1,nchild

      End If  ! End If (nodetype(lb) <= 2 .or. advance_all_levels)

      End Do  ! End Do lb = 1, lnblocks

      Call process_fetch_list(fetch_list,                              &
                              istack,                                  &
                              mype,                                    &
                              nprocs,                                  &
                              n_to_left,                               &
                              tag_offset)

!------Mark morton data up to date
       morton_limits_set = .True.

!------Store communication info for future use
       Call mpi_amr_write_restrict_comm(nprocs)

!------Deallocate any memory which was dynamically allocated for local 
!------use in this routine.
       If (Allocated(fetch_list)) deallocate(fetch_list)
       If (Allocated(n_to_left)) deallocate(n_to_left)

      Return

      Contains
        Subroutine expand_fetch_list

               If (Allocated(tfetch_list)) deallocate(tfetch_list)
               Allocate(tfetch_list(3,npts_neigh2))
               tfetch_list(:,:istack-1) = fetch_list(:,:istack-1)
               npts_neigh1 = npts_neigh1 + 3000
               npts_neigh2 = npts_neigh2 + 3000
               deallocate(fetch_list)
               Allocate(fetch_list(3,npts_neigh2))
               fetch_list(:,:istack-1) = tfetch_list(:,:istack-1)
               Deallocate(tfetch_list)

        End Subroutine expand_fetch_list
      End Subroutine pf_morton_bnd_restrict

