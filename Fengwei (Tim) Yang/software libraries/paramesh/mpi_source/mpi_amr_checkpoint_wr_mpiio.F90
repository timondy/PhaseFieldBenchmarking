!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_checkpoint_wr_mpiio
!! NAME
!!
!!   amr_checkpoint_wr_mpiio
!!
!! SYNOPSIS
!!
!!   call amr_checkpoint_wr_mpiio(file_num)
!!   call amr_checkpoint_wr_mpiio(file_num, l_with_guardcells)
!!
!!   call amr_checkpoint_wr_mpiio(integer, optional logical)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: file_num
!!     An integer number which will be appended to the end of the file name.
!!
!!   optional, logical, intent(in) :: l_with_guardcells
!!     If true, then guardcells are included in the checkpoint file.  Otherwise 
!!     (the default) they are not included.
!!   
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!! CALLS
!!
!! RETURNS
!!
!!   Does not return anything. 
!!
!! DESCRIPTION
!! 
!!  This is a dummy, placeholer routine for mpiio checkpointing.  This routine
!!  will return an error message if the mpiio checkpointing capability has
!!  not been installed.  The routine which actually writes an mpiio checkpoint
!!  file can be found in utilities/io/checkpointing/mpiio.
!!
!! AUTHORS
!!
!!   Kevin Olson (2004)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_checkpoint_wr_mpiio (file_num, l_with_guardcells)

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output statements
      Integer, intent(in) :: file_num
      Logical, optional, intent(in)  :: l_with_guardcells

!-----Local arrays and variables
      Integer :: mype, ierr

!-----Begin executable code.
      Call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

      If (mype == 0) Then
         Print *,' WARNING: you are calling amr_checkpoint_wr_mpiio '
         Print *,'          but your version of paramesh is not '
         Print *,'          yet configured to do this.          '
         Print *,'          Go to utilities/io/checkpoint/mpiio '
         Print *,'          in the main paramesh directory, run  '
         Print *,'          the INSTALL script, and recompile !!!'
      End If

      Return
      End Subroutine amr_checkpoint_wr_mpiio

