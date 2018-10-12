!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/find_surrblks
!! NAME
!!
!!   find_surrblks
!!
!! SYNOPSIS
!!
!!   Call find_surrblks()
!!
!! ARGUMENTS
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   local_tree_common
!!
!! CALLS
!!
!!   free_local_tree
!!   local_tree_build
!!   tree_search_for_surrblks
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a local tree has been created on
!!   each processor and the 'surr_blks' has been filled.
!!
!! DESCRIPTION
!!
!!   This routine finds the surrounding blocks at the same refinement level
!!   of all blocks.  It does this by first constucting a local tree on each 
!!   processor which represents the entire domain.  This local tree is then 
!!   searched by the list of local blocks to construct the surrounding block
!!   list.
!!
!! AUTHORS
!!
!!   Kevin Olson
!!
!!***

!-----A local data module used by the routine find_surrblks and other routines
#include "paramesh_preprocessor.fh"

      Module local_tree_common

      Use local_tree_module
      Type(node), Pointer, Save :: local_tree 

      End Module local_tree_common
      

      Subroutine find_surrblks()

!-----Use Statements
      Use local_tree_common

!-----Inlcude Statements
      Include 'mpif.h'

!-----Local variables
      Integer :: ierr
      Logical, Save :: first_call = .True.

      If (first_call) Nullify(local_tree)

      If (associated(local_tree)) Then
         Call free_local_tree(local_tree)
      End If

!-----Build the local tree
      Call local_tree_build()

!-----Search the local tree
      Call tree_search_for_surrblks()

      Call free_local_tree(local_tree)

      first_call = .false.

      End Subroutine find_surrblks

