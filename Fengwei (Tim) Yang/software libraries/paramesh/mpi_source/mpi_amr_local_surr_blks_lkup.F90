!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_local_surr_blks_lkup
!! NAME
!!
!!   mpi_amr_local_surr_blks_lkup
!!
!! SYNOPSIS
!!
!!   call mpi_amr_local_surr_blks_lkup(mype, lb, 
!!                                     surrblks, l_parent, psurrblks)
!!
!!   call mpi_amr_local_surr_blks_lkup(integer, integer, 
!!                                     integer, logical, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   integer, intent(in) :: lb
!!     Block for which neighboring blocks are to be found.
!!
!!   integer, intent(out) :: surrblks(:,:,:,:)
!!     List of the surrounding blocks which is returned.
!!
!!   logical, intent(in) :: l_parent
!!     Logical flag which indicates of parents of surrounding blocks are
!!     also to be found and returned.
!!   
!!   integer, intent(in) :: psurrblks
!!     The list of the surrounding blocks of the parent of block 'lb'.
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
!!   No other Paramesh routines are called
!!    
!! RETURNS
!!
!!   A list of the surrounding blocks of a block 'lb' on the local processor in the
!!   array 'surrblks' and also a list of the surrounding blocks of its parent block in
!!   the arrays 'psurrblks'.
!!
!! DESCRIPTION
!!
!!   This routine finds the addresses of surrounding blocks of 
!!   the block lb on the local processor, and of its parent from the list
!!   surrounding blocks computed by the subroutine mpi_amr_local_surr_blks.
!!   This routine finds these block addresses by searching the list list of 
!!   off-processor blocks which have previously been communicated.
!!
!! AUTHORS
!!
!!   Peter MacNeice (2000).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_amr_local_surr_blks_lkup(mype,lb,                 & 
                                              surrblks,l_parent,       &
                                              psurrblks)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      integer, intent(in)    ::  mype,lb
      integer, intent(out)   ::  surrblks(:,:,:,:)
      integer, intent(out)   ::  psurrblks(:,:,:,:)
      logical, intent(in)    ::  l_parent

!-----Local arrays and variables.
      Integer :: iblk
      Integer :: remote_block,remote_pe
      Integer :: ierrorcode,ierr
      Logical :: lfound

!-----This routine assumes that the grid blocks are ordered by morton
!-----number and that any blocks with different refinement levels but
!-----the same morton number are ordered from coarse to fine.

!-----Begin executable code.
      surrblks = -1
      psurrblks = -1

      surrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =                          & 
             surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb)

      If ( l_parent .And. (parent(1,lb) > 0) ) Then

          lfound = .False.
          If (parent(2,lb).ne.mype) Then

            iblk = ladd_strt(parent(2,lb))
            Do While(.Not.lfound.And.                                  & 
                     iblk <= ladd_end(parent(2,lb)))
              If (parent(1,lb) == laddress(1,iblk).And.                & 
                  parent(2,lb) == laddress(2,iblk) ) Then
                remote_block = iblk
                remote_pe    = mype
                lfound = .True.
              Else
                iblk = iblk+1
              End If
            End Do

          Else

            remote_block = parent(1,lb)
            remote_pe    = mype
            If (remote_block <= lnblocks) lfound = .True.

          End If

          If (.Not.lfound) Then
          Write(*,*) 'Error in mpi_amr_local_surr_blks_lkup : ',       & 
                     'remote block ',parent(:,lb),                     & 
                     ' not located on pe ',mype,                       & 
                     ' while processing blk ',lb,mype
          Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
          End If

          If (lfound) Then
          psurrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =                     & 
               surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,remote_block)
          End If

      End If  ! End If ( l_parent .And. (parent(1,lb) > 0) )

      Return
      End Subroutine mpi_amr_local_surr_blks_lkup
