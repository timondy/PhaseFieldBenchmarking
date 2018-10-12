!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_boundary_block_info
!! NAME
!!
!!   mpi_amr_boundary_block_info
!!
!! SYNOPSIS
!!
!!   Call mpi_amr_boundary_block_info (mype, nprocs)
!!
!!   Call mpi_amr_boundary_block_info (integer, integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype           
!!        The local processor number.
!!
!!   integer, intent(in) :: nprocs
!!        The number of MPI processes.
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
!!   Does not Call any other Paramesh routines.
!!
!! RETURNS
!!
!!   Does not return anything.  
!!
!! DESCRIPTION
!!
!!   This routine constructs a list of block faces which are external
!!   boundaries.
!!
!!
!! AUTHORS
!!
!!   Peter MacNeice (November 2002)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_amr_boundary_block_info(mype,nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer,Intent(in) :: mype,nprocs

!-----Local arrays and variables
      Integer,Allocatable :: ib_global(:)
      Integer,Allocatable :: ib_count_send(:),ib_count_recv(:)
      Integer,Allocatable :: ib_global_recv(:)
      Integer,Allocatable :: ib_recv_index(:)
      Integer,Allocatable :: ib_send_index(:)

      Integer :: ib,ib0,ib1,lb,ib_sum
      Integer, Allocatable :: ib_list(:,:)
      Integer :: i,j,k,i1,i2,iproc,ierror


!-----Begin executable code.

!-----compute size for ib_list and then allocate it.
      ib0 = 0
      Do lb = 1,lnblocks
        ib = minval(surr_blks(1,:,:,:,lb))
        If (ib <= -20) Then
          Do k = 1,1+2*k3d
          Do j = 1,1+2*k2d
          Do i = 1,3
            ib1 = surr_blks(1,i,j,k,lb)
            If (ib1 <= -20) Then
               ib0 = ib0 + 1
            End If
          End Do
          End Do
          End Do
        End If
      End Do
      If (Allocated(ib_list)) deallocate(ib_list)
      allocate(ib_list(6,ib0))

!-----construct on pe list of blocks next to boundary
      ib0 = 0
      Do lb = 1,lnblocks
        ib = minval(surr_blks(1,:,:,:,lb))
        If (ib <= -20) Then
          Do k = 1,1+2*k3d
          Do j = 1,1+2*k2d
          Do i = 1,3
            ib1 = surr_blks(1,i,j,k,lb)
            If (ib1 <= -20) Then
            ib0 = ib0 + 1
            ib_list(1,ib0) = lb
            ib_list(2,ib0) = mype
            ib_list(3,ib0) = i
            ib_list(4,ib0) = j
            ib_list(5,ib0) = k
            ib_list(6,ib0) = ib
            End If
          End Do
          End Do
          End Do
        End If
      End Do

      If (.not.allocated(ib_global))                                   & 
                               Allocate(ib_global(0:nprocs-1))
      If (.not.allocated(ib_global_recv))                              & 
                               Allocate(ib_global_recv(0:nprocs-1))
      If (.not.allocated(ib_count_send))                               & 
                               Allocate(ib_count_send(0:nprocs-1))
      If (.not.allocated(ib_count_recv))                               & 
                               Allocate(ib_count_recv(0:nprocs-1))
      If (.not.allocated(ib_recv_index))                               & 
                        Allocate(ib_recv_index(0:nprocs-1))
      If (.not.allocated(ib_send_index))                               & 
                        Allocate(ib_send_index(0:nprocs-1))
      ib_global_recv = 0
      ib_global = ib0

      Call MPI_ALLTOALL (ib_global     ,1,MPI_INTEGER,                 & 
                         ib_global_recv,1,MPI_INTEGER,                 & 
                         MPI_COMM_WORLD,ierror)
      ib_count_recv = ib_global_recv*6
      Do iproc = 0,nprocs-1
        ib_count_send(iproc) = ib_count_recv(mype) 
      End Do

!-----Compute displacements in all to all message buffer, in bytes.
      ib_recv_index(0) = 0
      If (nprocs > 1) Then
      Do iproc = 1,nprocs-1
        ib_recv_index(iproc) = ib_recv_index(iproc-1)                  & 
                              + ib_global_recv(iproc-1)*6
      End Do
      End If
      Do iproc = 0,nprocs-1
        ib_send_index(iproc) = ib_recv_index(mype) 
      End Do

!-----exchange no of boundary blocks on each processor
      Call comm_int_sum_to_all(ib_sum,ib0)
      bc_block_neighs_length = ib_sum

!-----allocate storage for global list
      If (Allocated(bc_block_neighs)) Deallocate(bc_block_neighs)
      Allocate(bc_block_neighs(6,ib_sum))
      If (Allocated(bc_block_neighs_send))                             & 
                                     deallocate(bc_block_neighs_send)
      Allocate(bc_block_neighs_send(6,ib_sum))

!-----exchange info between procs
!-----Put local data into is correct place on the list
      bc_block_neighs_send = -1
      bc_block_neighs = -5
      i1 = ib_recv_index(mype)/6 + 1
      i2 = i1 + ib0 - 1
      
      If (ib0 > 0) bc_block_neighs_send(:,i1:i2) = ib_list(:,1:ib0)

!-----Exchange lists between all procs
      Call MPI_ALLTOALLV (bc_block_neighs_send,ib_count_send,          & 
                                          ib_send_index,MPI_INTEGER,   & 
                          bc_block_neighs,ib_count_recv,               & 
                                          ib_recv_index,MPI_INTEGER,   & 
                          MPI_COMM_WORLD,ierror)

      Deallocate(ib_list)
      Deallocate(ib_global)
      Deallocate(ib_global_recv)
      Deallocate(ib_count_send)
      Deallocate(ib_count_recv)
      Deallocate(ib_recv_index)
      Deallocate(ib_send_index)

!-----Set status flag
      bc_block_neighs_status = 100

      Return
      End Subroutine mpi_amr_boundary_block_info
