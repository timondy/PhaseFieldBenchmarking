!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_morton_bnd_fluxcon
!! NAME
!!
!!   mpi_morton_bnd_fluxcon
!!
!! SYNOPSIS
!!
!!   Call mpi_morton_bnd_fluxon(mype, nprocs, tag_offset)
!!   Call mpi_morton_bnd_fluxcon(integer, integer, integer)
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
!!    mpi_amr_write_guard_comm
!!    process_fetch_list
!!
!! RETURNS
!!
!!    Does not return anything.
!!
!! DESCRIPTION
!!
!!   This routine does a communications analysis for the flux and edge fixups
!!   and constructs and stores lists of off-processor blocks which are to be 
!!   communicated during these steps.  Also stored are which sections of 
!!   the blocks to be fetched.
!!
!! AUTHORS
!!
!!    Written by Peter MacNeice  and Michael Gehmeyr, February 2000.
!!    Major simplification and rewrite by Kevin Olson August 2007.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use mpi_morton
      Use constants

      Use paramesh_mpi_interfaces, only : mpi_amr_write_guard_comm,    & 
                                          process_fetch_list

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      integer, intent(in)    ::  mype,nprocs
      integer, intent(inout) ::  tag_offset

!-----Local Variables
      real    :: eps,accuracy

      Integer :: lb,i,j,k,j00
      Integer :: ierror
      Integer :: max_no_of_blocks
      Integer :: istack, ioff, joff, koff
      Integer :: isize, isrc, idest, itag, nrecv, itag_count
      Integer :: ii, jj, kk
      Integer :: nguarda 
      Integer :: iproc, gid
      Integer :: ie, ie_max, if1, jf1, kf1, if2, jf2, kf2
      Integer :: north, south
      Integer :: npts_neigh1,npts_neigh2, sendmsg

      Integer,Dimension(:),  Allocatable :: n_to_left
      Integer,Dimension(:),  Allocatable :: recvrequest
      Integer,Dimension(:,:),Allocatable :: recvstatus
      Integer,Dimension(:,:),Allocatable :: fetch_list
      Integer,Dimension(:,:),Allocatable :: tfetch_list

!-----Begin executable code.

      accuracy = 100./10.**precision(accuracy)
      eps = accuracy                                                                                                     
      nguarda = max(nguard,nguard_work)

      npts_neigh1 = npts_neigh
      npts_neigh2 = npts_neigh+100
      Allocate(fetch_list(3,npts_neigh2))
      Allocate(n_to_left(0:nprocs-1))

!-----store the max no of blocks on any one processor
      Call MPI_ALLREDUCE(lnblocks,                                     & 
                         max_no_of_blocks,                             & 
                         1,                                            & 
                         MPI_INTEGER,                                  & 
                         MPI_MAX,                                      & 
                         MPI_COMM_WORLD,                               & 
                         ierror)

!-----COMPUTE the number of blocks to the 'left' (ie. stored on processors with
!-----smaller process ids) of every other processor
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

!-----Construct a list of potential neighbors of all blocks on this
!-----processor.
!-----Exclude any which are on this processor.

      istack = 0

      Do lb = 1, lnblocks

      If (nodetype(lb) == 1) Then

!------ADD OFF PROCESSOR NEIGHBORS OF BLOCK 'lb' TO FETCH LIST
       Do k = 1,1+2*k3d
        Do j = 1,1+2*k2d
        Do i = 1,3
        If (surr_blks(1,i,j,k,lb) > 0 .and.                            & 
            surr_blks(2,i,j,k,lb) .ne. mype .and.                      &
            surr_blks(3,i,j,k,lb) == 2) Then

            istack = istack + 1
            If (istack > npts_neigh1) Call expand_fetch_list
            fetch_list(1,istack) = surr_blks(1,i,j,k,lb)
            fetch_list(2,istack) = surr_blks(2,i,j,k,lb)
            j00 = j

!-----------if this block is a polar block then change the way j is applied
!-----------in the formula
            
            j00 = j00 + 1 - k2d
            If(spherical_pm) Then
               If(lsingular_line) Then
                  If(abs(bnd_box(1,2,lb)) < eps.and.j+1-k2d == 1) Then
                     j00 = 3
                  Else If(abs(bnd_box(2,2,lb)-pi) < eps .and.             &
                          j+1-k2d == 3) Then
                     j00 = 1
                  End If
               End If
            End If

!-----------compute message type - note this index is computed to reflect the part
!-----------of the remote block to be acquired, not the part of the local blocks
!-----------guardcells which will be filled.
            kk = k + 1 - k3d
            fetch_list(3,istack) = (4-i)+((4-j00)-1)*3+((4-kk)-1)*9
            if(nguarda.gt.nmax_lays) fetch_list(3,istack) = 14

         End If  ! End If surr_blks(1,i,j,k,lb) > 0 .and. ...

         End Do  ! End Do i = 1,3
         End Do  ! End Do j = 1,1+2*k2d
        End Do  ! End Do k = 1,1+2*k3d

      End If  ! End If (nodetype(lb) == 1)

      End Do  ! End Do lb = 1, lnblocks

      If (nedgevar1 > 0) Then
!-----Mark Edges-------!
!-----If edge averaging is required then we need to identify leaf blocks
!-----which have refinement off their diagonal edges but not in the neighbors
!-----bounding that diagonal element. In this case we will have to make
!-----sure that the edge integral on the local block match those of the
!-----refined diagonal neighbor, if we wish to preserve a div B constraint.

      no_of_diagonal_edges = 0
      Do lb = 1, lnblocks

         If (nodetype(lb) == 1) Then

!-------Cycle through the 12 edges of this block
            ie_max = 2
            if (ndim == 2) ie_max = 4
            if (ndim == 3) ie_max = 12
            Do ie = 1,ie_max
               i = 2
               j = 1+k2d
               k = 1+k3d
               if1 = 2
               jf1 = 1+k2d
               kf1 = 1+k3d
               if2 = 2
               jf2 = 1+k2d
               kf2 = 1+k3d
               If (ie == 1) Then
                  i = 1
                  j = 1
                  if1 = 1
                  jf2 = 1
               Elseif (ie == 2) Then
                  i = 1
                  j = 1+2*k2d
                  if1 = 1
                  jf2 = 1+2*k2d
               Elseif (ie == 3) Then
                  i = 3
                  j = 1
                  if1 = 3
                  jf2 = 1
               Elseif (ie == 4) Then
                  i = 3
                  j = 3
                  if1 = 3
                  jf2 = 3
               Elseif (ie == 5) Then
                  j = 1
                  k = 1
                  jf1 = 1
                  kf2 = 1
               Elseif (ie == 6) Then
                  j = 3
                  k = 1
                  jf1 = 1+2*k2d
                  kf2 = 1
               Elseif (ie == 7) Then
                  j = 1
                  k = 3
                  jf1 = 1
                  kf2 = 3
               Elseif (ie == 8) Then
                  j = 3
                  k = 3
                  jf1 = 1+2*k2d
                  kf2 = 3
               Elseif (ie == 9) Then
                  i = 1
                  k = 1
                  if1 = 1
                  kf2 = 1
               Elseif (ie == 10) Then
                  i = 1
                  k = 3
                  if1 = 1
                  kf2 = 3
               Elseif (ie == 11) Then
                  i = 3
                  k = 1
                  if1 = 3
                  kf2 = 1
               Elseif (ie == 12) Then
                  i = 3
                  k = 3
                  if1 = 3
                  kf2 = 3
               End If
!--------------If corner block is a nodetype 2 and the blocks on the 
!--------------faces adjacent to this block are not refined (nodetype = 1) 
!--------------then mark the edge
               If (surr_blks(3,i,j,k,lb) == 2 .and.                   &
                   surr_blks(3,if1,jf1,kf1,lb) == 1 .and.             &
                   surr_blks(3,if2,jf2,kf2,lb) == 1) Then
                  no_of_diagonal_edges = no_of_diagonal_edges + 1
                  edge_mark(6,1,no_of_diagonal_edges) = ie
                  edge_mark(6,2,no_of_diagonal_edges) = lb
                  If (surr_blks(1,i,j,k,lb) < -1) Then
                     edge_mark(:,3:4,no_of_diagonal_edges) = -1
                  Else
                     edge_mark(6,3,no_of_diagonal_edges) =             &
                       surr_blks(1,i,j,k,lb)
                     edge_mark(6,4,no_of_diagonal_edges) =             &
                       surr_blks(2,i,j,k,lb)
                  End if
               End If

            End Do  ! End Do ie = 1,12

         End If  ! End If (nodetype(lb) == 1)

      End Do  ! End Do lb = 1, lnblocks
      End If  ! End If (nedgevar1 > 0)

      Call process_fetch_list(fetch_list,                              &
                              istack,                                  &
                              mype,                                    &
                              nprocs,                                  &
                              n_to_left,                               &
                              tag_offset)

! Mark morton data up to date
       morton_limits_set = .True.

! Store communication info for future use
       Call mpi_amr_write_flux_comm(nprocs)

! Deallocate any memory which was dynamically allocated for local use in this
! routine.
      If (Allocated(fetch_list))  Deallocate(fetch_list)
      If (Allocated(n_to_left))   Deallocate(n_to_left)
      If (Allocated(recvrequest)) Deallocate( recvrequest )
      If (Allocated(recvstatus))  Deallocate( recvstatus )

      Return

      Contains
        Subroutine expand_fetch_list

               if (Allocated(tfetch_list)) Deallocate(tfetch_list)
               Allocate(tfetch_list(3,npts_neigh2))
               tfetch_list(:,:istack-1) = fetch_list(:,:istack-1)
               npts_neigh1 = npts_neigh1 + 3000
               npts_neigh2 = npts_neigh2 + 3000
               Deallocate(fetch_list)
               Allocate(fetch_list(3,npts_neigh2))
               fetch_list(:,:istack-1) = tfetch_list(:,:istack-1)
               Deallocate(tfetch_list)

        End Subroutine expand_fetch_list
      End Subroutine mpi_morton_bnd_fluxcon



