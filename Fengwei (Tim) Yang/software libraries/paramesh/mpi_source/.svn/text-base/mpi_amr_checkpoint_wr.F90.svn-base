!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_checkpoint_wr
!! NAME
!!
!!   amr_checkpoint_wr
!!
!! SYNOPSIS
!!
!!   Call amr_checkpoint_wr(file_num)
!!   Call amr_checkpoint_wr(file_num, l_with_guardcells, check_format, 
!!                          user_attr1, user_attr2, user_attr3, 
!!                          user_attr4, user_attr5)
!!
!!   Call amr_checkpoint_wr(integer, optional logical, optional char*80,
!!                          optional real, optional real, optional real, 
!!                          optional real, optional real)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: file_num
!!     An integer number which be appended to the end of the file name.
!!
!!   optional, logical, intent(in) :: l_with_guardcells
!!     If true, then guardcells are included in the checkpoint file.  Otherwise 
!!     (the default) they are not included.
!!   
!!   optional, character(len=80), intent(in) :: check_format
!!     Argument describing what type of output to use.  Currently supports
!!     fortran binary(the default), parallel hdf5 output, or parallel mpiio
!!     in native binary.
!!     To produce parallel hdf5 output:
!!       character(len=80) :: checkf
!!       checkf = 'hdf5'
!!       Call amr_checkpoint_wr(...., check_format=checkf, ....) 
!!
!!   optional, real, intent(in) :: user_attr1(2,3,4,5)
!!     Arguments which allow the user to add some extra information to the file.  
!!     Currently only 5 real numbers can be added.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_checkpoint_wr_default
!!   amr_checkpoint_wr_hdf5
!!   amr_checkpoint_wr_mpiio
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a checkpoint file has been written.
!!
!! DESCRIPTION
!! 
!!  Subroutine to checkpoint runs using PARAMESH.
!!  It writes out the tree data structure and data stored in PARAMESH blocks.
!!  Optionally, a user may add a small amout of attribute data to the files written.
!!  If the default bevaviour is selected, writes are done serially by processor 0. 
!!  I.e. data is collected from other processors and then written out.  
!!  The default behaviour USES UNFORMATTED DIRECT I/O.
!!  Optionally, the user can specify if the file will be written using the HDF5
!!  portable data format.  If supported on the system, selecting hdf5 output will 
!!  also result in the file will be written in parallel.
!!
!!  The files produced will have names of the form 'paramesh_chk_######' or 
!!  'paramesh_chk_######.hdf5'. where '######' is the file_num argument passed into
!!  this routine.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson for parallel hdf5 and 
!!   mpiio (2004-2005).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_checkpoint_wr(file_num,                           & 
                                   l_with_guardcells,                  & 
                                   check_format,                       & 
                                   user_attr_1,                        & 
                                   user_attr_2,                        & 
                                   user_attr_3,                        & 
                                   user_attr_4,                        & 
                                   user_attr_5)

!-----Use Statements
      Use paramesh_interfaces, Only : amr_checkpoint_wr_default,       & 
                                      amr_checkpoint_wr_hdf5,          & 
                                      amr_checkpoint_wr_mpiio

      implicit none

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguemetns
      Integer, Intent(in)                      :: file_num
      Logical, Optional, Intent(in)            :: l_with_guardcells
      Character (len=80), Optional, Intent(in) :: check_format
      Real, Optional, Intent(in)               :: user_attr_1,         & 
                                                  user_attr_2,         & 
                                                  user_attr_3,         & 
                                                  user_attr_4,         & 
                                                  user_attr_5
!-----Local variables and arrays.
      Character (len=80) :: check_format_in
      Integer            :: mype, ierr
      
!-----Begin executable code.
      Call MPI_COMM_RANK (MPI_COMM_WORLD, mype, ierr)

      check_format_in = 'default'
      If (Present(check_format)) check_format_in = check_format
      check_format_in = trim(check_format_in)

      If (check_format_in(1:7) == 'default') Then
         Call amr_checkpoint_wr_default (file_num, l_with_guardcells,  & 
           user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      ElseIf (check_format_in(1:4) == 'hdf5') Then
         Call amr_checkpoint_wr_hdf5 (file_num, l_with_guardcells,     & 
           user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      ElseIf (check_format_in(1:5) == 'mpiio') Then
         Call amr_checkpoint_wr_mpiio (file_num, l_with_guardcells,    & 
           user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      Else
         If (mype == 0) Then
            Print *,' UNRECOGNIZED I/O FORMAT, CHECKPOINT NOT WRITTEN '
         End If
      End If

      End Subroutine amr_checkpoint_wr


