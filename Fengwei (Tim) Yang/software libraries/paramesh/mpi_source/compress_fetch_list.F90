!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/compress_fetch_list
!! NAME
!!
!!   compress_fetch_list
!!
!! SYNOPSIS
!!
!!   Call compress_fetch_list (fetch_list,                     
!!                             istack,                          
!!                             no_of_remote_neighs,             
!!                             mype,                            
!!                             nprocs,                          
!!                             n_to_left)
!!   Call compress_fetch_list (Integer array,                     
!!                             Integer,                          
!!                             Integer,             
!!                             Integer,                            
!!                             Integer,                          
!!                             Integer)
!!
!! ARGUMENTS
!!
!!   fetch_list -> A list of block ids (for a particular blocks) 
!!                 to be compressed
!!   istack     -> Size of the list initially.
!!   no_of_remote_neighs -> Size of list after compression.
!!   mype       -> processor id of calling processor.
!!   nprocs     -> no. of processors.
!!   n_to_left  -> of blocks in the morton ordered list with morton numbers
!!                 less than that of the block calling this routine.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   mpi_morton
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_q_sort
!!   rationalize_fetch_list
!!
!! RETURNS
!!
!!   Returns a compressed list in the array 'fetch_list'
!!
!! DESCRIPTION
!!
!!   This routine returns a compressed list of off processor blocks which
!!   has been constucted in the mpi_morton_bnd_XXX routines for a particular
!!   block.  It sorts and eliminates any redundent entries in the list.
!!   It does this by sorting the list in order of increasing morton number, 
!!   and then sorting each sub-list with the same morton number in order of
!!   increasing refinement level. Then it removes any identical entries.
!!   This routine takes the place of the old routine 'compress_list'.
!!
!! AUTHORS
!!
!!    Based on the routine 'compress_list' originally written by Peter MacNeice
!!    and
!!    rewritten by Kevin Olson, 2007,2008.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine compress_fetch_list (fetch_list,                      &
                                      istack,                          &
                                      no_of_remote_neighs,             & 
                                      mype,                            &
                                      nprocs,                          &
                                      n_to_left)

!-----Use Statements
      Use mpi_morton
      Use paramesh_interfaces, only : amr_q_sort

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(inout), Dimension(:,:) :: fetch_list
      Integer, intent(out)   :: no_of_remote_neighs
      Integer, intent(in)    :: istack, mype, nprocs
      Integer, intent(in)    :: n_to_left(0:nprocs-1)

!-----Local Variables
      Integer, Dimension(istack) :: t1fetch_list, t2fetch_list
      Integer, Dimension(istack) :: indx
      Integer :: i, itemp, istart, iend, level, j_pe
      Integer :: jstack, j
      Integer :: ierrorcode, ierr
      Integer :: itest

!-----Begin Executable Code.

!-----set indx so fetch_list(3,:) can be permuted.
      Do i=1,istack
         indx(i) = i
         !-----Compute a unique global Ide of the blocks in the list.
         ! CEG moved from separate loop
         t1fetch_list(i) = n_to_left(fetch_list(2,i)) + fetch_list(1,i)
      End Do

!-----Sort Global Ids
      Call amr_q_sort(t1fetch_list(:), istack, ia=indx(:))

!-----Now reorder the morton number part of fetch_list.
      fetch_list(1,1:istack) = t1fetch_list(1:istack)

      t2fetch_list(1:istack) = fetch_list(2,1:istack)
      Do i=1,istack
         fetch_list(2,i) = t2fetch_list(indx(i))
      End Do

      t2fetch_list(1:istack) = fetch_list(3,1:istack)
      Do i=1,istack
         fetch_list(3,i) = t2fetch_list(indx(i))
      End Do

!-----If any entries have the same block id
!-----but multiple data request types then mark them to fetch the complete
!-----blocks.
      istart = 1
      iend = 1
      i = 2
      Do While(i <= istack)
! CEG rewritten to stop incorrect memory access
        Do While(i.le.istack)
           if (fetch_list(1,i).eq.fetch_list(1,i-1)) then
              i = i+1
           else
              exit
           endif
        End Do
        iend = i-1
        If (istart < iend) Then

          call rationalize_fetch_list(fetch_list,istart,iend,          & 
                                      size(fetch_list,2))

        End If
        istart = iend+1
        i = istart+1

      End Do

!-----Finally remove any repetition of elements in this list.
      jstack = 0
      itest = -1
      Do i = 1,istack
         If (fetch_list(1,i).ne.itest) Then
           jstack = jstack + 1
           fetch_list(:,jstack) = fetch_list(:,i)
           itest = fetch_list(1,i)
         End If
      End Do

      no_of_remote_neighs = jstack

!-----Convert from global id, back to local id, proc id
      Do i = 1, no_of_remote_neighs
         fetch_list(1,i) = fetch_list(1,i) -                          &
             n_to_left(fetch_list(2,i))
      End do

      Return
      End Subroutine compress_fetch_list

