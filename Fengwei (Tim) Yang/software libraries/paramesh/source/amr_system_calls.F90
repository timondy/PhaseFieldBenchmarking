!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


!!****f* source/amr_abort
!! NAME
!!   amr_abort
!!
!! SYNOPSIS
!!
!!   Call amr_abort()
!!
!! ARGUMENTS
!! 
!!   None
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
!! DESCRIPTION
!!
!!   A wrapper routine for the MPI_ABORT command.
!!
!! AUTHORS
!!
!!   Peter MacNeice
!!
!!***

      Subroutine amr_abort()

      Implicit None

!-----Include statments.
      Include 'mpif.h'

!-----Local variables.
      Integer :: ierrorcode,ierr

!-----Begin executable code.
      Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)

      Return
      End Subroutine amr_abort


!!****f* source/amr_abort
!! NAME
!!   amr_abort
!!
!! SYNOPSIS
!!
!!   Call amr_abort()
!!
!! ARGUMENTS
!! 
!!   Integer, Intent(in) :: iunit  The unit number to flush the io buffer.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!! CALLS
!!
!! RETURNS
!!
!! DESCRIPTION
!!
!!   A stub routine for whatever flush call you have on your specific system.
!!
!! AUTHORS
!!
!!   Peter MacNeice
!!
!!***

      Subroutine amr_flush(iunit)

      Implicit None

!-----Local variables
      Integer, Intent(in) :: iunit

      Return
      End Subroutine amr_flush

