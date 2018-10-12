!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_refine_derefine
!! NAME
!!
!!   amr_refine_derefine
!!
!! SYNOPSIS
!!
!!   call amr_refine_derefine()
!!
!!
!! ARGUMENTS
!!
!!   NO ARGUMENTS
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
!!   timings
!!   io
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_check_refine
!!   amr_check_derefine
!!   amr_refine_blocks
!!   amr_derefine_blocks
!!   amr_morton_order
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit blocks marked for refinement have 
!!   been created, blocks marked for derefinement have been eliminated, and 
!!   the blocks have been reordered to acheive load balance.
!!
!! DESCRIPTION
!! 
!!  Subroutine to refine or derefine blocks.  This routine is called by the 
!!  user who sets the refine or derefine flags to be true or false.  These 
!!  flags are logical variables called 'refine' and 'derefine' and are stored 
!!  for each block.  If a block is marked for refinement, amr_refine_derefine 
!!  will create its new child blocks.  Also, tests will be executed to see if 
!!  any other blocks need to also refine (by calling amr_check_refine) to 
!!  ensure that the a jump in refinement of more than one level is not created.
!!
!!  If a block is marked for derefinement, amr_refine_derefine first checks 
!!  to make that the block can derefine and not create a jump in refinement of
!!  more than one level.  A check is also run to check that all the siblings of 
!!  the derefining block's siblings are also marked for derefinement.  If these 
!!  tests succeed, the block is removed from the list of blocks.
!!
!!  Once these operations are completed, the routine 'amr_morton_order' is 
!!  called and the tree data structure is reorganized to acheive load balance 
!!  using a morton space filling curve.  After this routine is called, the 
!!  routine 'amr_redist_blk' is called, which actually moves the block data 
!!  into the correct positions in the morton order list of blocks.
!!
!!  Finally, the routine 'amr_morton_process' is called.  This routine 
!!  computes the communications patterns needed for guardcell filling, 
!!  restriction, and prologation and stores them for later use. 
!!
!! AUTHORS
!!
!!  Kevin Olson (1997)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_refine_derefine

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_interfaces, Only : amr_check_refine,                & 
                                      amr_check_derefine,              & 
                                      amr_refine_blocks,               & 
                                      amr_derefine_blocks,             & 
                                      amr_morton_order
      Use paramesh_mpi_interfaces, Only : mpi_amr_singular_line

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Local variables and arrays.
      Integer :: lnblocks2,tot_blocks,tot_blocksa,icontinue
      Integer :: icontinue_ref,icontinue_deref
      Integer :: icontinue2,max_blocks
      Integer :: min_blocks
      Integer :: nprocs,mype
      Integer :: i,l
      Integer :: lnblocks_old=0
      Integer :: istrategy
      Integer :: ierrorcode, ierr
      Logical :: l_move_solution
!      Logical :: refinet(maxblocks_tr)
      Logical,allocatable :: refinet(:)
      Logical,save :: first_Call = .True.

!-----Begin executable code.
! CEG allocate memory
  allocate(refinet(maxblocks_tr))

      Call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      Call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

!-----error trap for lrefine_max and lrefine_min
      If (lrefine_max < 1.Or.lrefine_max > 100) Then
        Write(*,*) 'PARAMESH error : lrefine_max has a bad value'      & 
                   ,lrefine_max
        Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      If (lrefine_min < 1.Or.lrefine_min > 100.Or.                     & 
         lrefine_min > lrefine_max) Then
        Write(*,*) 'PARAMESH error : lrefine_min or lrefine_max ',     & 
         'has a bad value : lrefine_min= ',lrefine_min, & 
         ' lrefine_max= ',lrefine_max
        Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      End If

!-----enforce refinement level limits
!-----first upper limit
      Where (lrefine==lrefine_max) refine = .False.
!-----Then lower limit
      Where (lrefine==lrefine_min) derefine = .False.
!-----finally force grid to refine toward base level if too coarse
      Where ( (lrefine<lrefine_min) .And.                              & 
              (nodetype==1) ) refine = .True.

      refinet(1:lnblocks) = refine(1:lnblocks)
      newchild(:) = .FALSE.

!-------CHECK derefinements and refinements

!-------test to see if any refinements have been requested.
      icontinue=0
      icontinue_deref = 0
      If (lnblocks > 0) Then
         Do l = 1,lnblocks
            If (nodetype(l) == 1.And.refine(l)) Then
               icontinue=1
               Exit
            End If
         End Do
      End If

! CEG added for growing domains to also include the case of new blocks being added

      if (lnblocks_old.ne.lnblocks .and. .not.first_Call .and.  icontinue.eq.0) then
!         print*, 'Grid has grown spotted in amr_refine_derefine by processor ', mype, lnblocks_old, lnblocks
         icontinue = 1
      endif
         

!   print *,'CEG AllReduce in mpi_amr_refine_derefine'
!     called once per processor at end of timestep
      Call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,          & 
                          MPI_MAX,MPI_COMM_WORLD,ierr)

      icontinue_ref = icontinue2

      If (spherical_pm) Then
         istrategy = 0
         If (ndim == 3.And.lsingular_line) istrategy = 1
         Call mpi_amr_singular_line(istrategy,nprocs)
      End If

      Call amr_check_refine (nprocs,mype,icontinue_ref)
!      print *,'test 1'

!-----test to see if any derefinements have been requested
      icontinue=0
      icontinue_deref = 0
      If (lnblocks > 0) Then
         Do l=1,lnblocks
            If (nodetype(l) == 1.And.derefine(l)) Then
               icontinue=1
               Exit
            End If
         End Do
      End If
      Call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,          & 
                          MPI_MAX,MPI_COMM_WORLD,ierr)
      icontinue_deref = icontinue2

      If (icontinue_deref > 0) Call amr_check_derefine (mype)
!      print *,'test 2'

!-----test to see if any derefinements have been requested which passed 
!-----the tests in amr_check_derefine
      If (icontinue_deref > 0) Then

      icontinue=0
      icontinue_deref = 0
      If (lnblocks > 0) Then
         Do l=1,lnblocks
            If (nodetype(l) == 1.And.derefine(l)) Then
               icontinue=1
               Exit
            End If
         End Do
      End If
      Call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,          & 
                          MPI_MAX,MPI_COMM_WORLD,ierr)
      icontinue_deref = icontinue2

      End If

!      print *,'test 2b', lnblocks_old, lnblocks
      If (icontinue_ref   == 0 .And.                                   & 
          icontinue_deref == 0 .And.                                   & 
          .Not.first_Call) then
!CEG added deallocate in here
           deallocate(refinet)
           return
      endif
      first_Call = .False.
        
!      print *,'test 3', lnblocks_old, lnblocks

!-----NOW Actually refine and derefine the mesh

      lnblocks_old = lnblocks
      If (icontinue_ref > 0) Call amr_refine_blocks (nprocs,mype)
!      print *,'test 4'

      If (icontinue_deref > 0.Or.icontinue_ref > 0)                    & 
                     Call amr_derefine_blocks(lnblocks_old,mype)

      Call MPI_ALLREDUCE (lnblocks,tot_blocks,1,MPI_INTEGER,           & 
                          MPI_SUM,MPI_COMM_WORLD,ierr)
!      print *,'test 5', lnblocks_old, lnblocks


!-----set work values
      work_block(:) = 0.
      Do i = 1,lnblocks
         if (nodetype(i) == 1) work_block(i) = 2.        !<<< USER EDIT
         if (nodetype(i) >= 2) work_block(i) = 1.        !<<< USER EDIT
      end do
      l_move_solution = .True.
!      print *,'test 6'

      Call amr_morton_order (lnblocks_old,nprocs,mype,                 & 
                             l_move_solution)
!      print *,'test 7', lnblocks_old, lnblocks
!-----I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      Call MPI_ALLREDUCE (lnblocks2,tot_blocksa,1,MPI_INTEGER,         & 
                          MPI_SUM,MPI_COMM_WORLD,ierr)
      Call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,          & 
                          MPI_MAX,MPI_COMM_WORLD,ierr)
      Call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,          & 
                          MPI_MIN,MPI_COMM_WORLD,ierr)

      Call amr_morton_process()
!      print *,'test 8', lnblocks_old, lnblocks

!-----Set up an array of cell sizes for each grid refinement level.
!-----These can be Used to minimize variation due to roundoff, but
!-----should ONLY be used with a uniformly spaced grid.
      level_cell_sizes = 0.
      level_cell_sizes(1,1) = (grid_xmax-grid_xmin)/real(nxb)
      If (ndim > 1)                                                    & 
        level_cell_sizes(2,1) = (grid_ymax-grid_ymin)/real(nyb)
      If (ndim == 3)                                                   & 
        level_cell_sizes(3,1) = (grid_zmax-grid_zmin)/real(nzb)
      Do i=2,lrefine_max
        level_cell_sizes(1:ndim,i) = .5*level_cell_sizes(1:ndim,i-1)
      End Do

!-----set grid modification flag
      grid_changed = 1
      grid_analysed_mpi = 1

!      print *,'test 9', lnblocks_old, lnblocks
      Call mpi_amr_boundary_block_info(mype,nprocs)

! CEG deallocate memory
      deallocate(refinet)

      lnblocks_old = lnblocks

      Return
      End Subroutine amr_refine_derefine

