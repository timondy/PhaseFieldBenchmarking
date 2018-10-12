!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_morton_process
!! NAME
!!
!!   amr_morton_process
!!
!! SYNOPSIS
!!
!!   Call amr_morton_process()
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
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   io
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   mpi_morton_bnd
!!   mpi_morton_bnd_prolong
!!   mpi_morton_bnd_fluxcon
!!   mpi_morton_bnd_restrict
!!   mpi_setup
!!
!! RETURNS
!!
!!   Does not return anything.
!!
!! DESCRIPTION
!!
!!   This routine manages the construction of list of blocks which are
!!   communicated for guardcell filling, prolongation, restriction, and
!!   flux and edge fixups.  It does this by simply calling the routines
!!   mpi_morton_bnd, mpi_morton_bnd_prolong, mpi_morton_bnd_restrict and
!!   mpi_morton_bnd_fluxcon.  
!!
!! AUTHORS
!!
!!    Kevin Olson
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_morton_process()

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use io

      Use paramesh_mpi_interfaces, Only : mpi_morton_bnd,              & 
                                          mpi_morton_bnd_prolong,      & 
                                          mpi_morton_bnd_fluxcon,      & 
                                          mpi_morton_bnd_restrict,     & 
                                          mpi_setup


      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Local Variables
      Integer :: nprocs,mype,tag_offset,ierr

!-----Begin Executable code.

!      print *,'morton mpi'
      Call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      Call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

!-----call setup routines in preparation for calling 
!-----all the mpi_morton_bnd_XXX routines.

!-----Find the coordinate ranges
!      print *,'morton 0'
      Call mpi_amr_global_domain_limits

!      print *,'morton 1'
      Call mpi_setup(mype,nprocs)
!      print *,'morton 2'

!-----Set up surrounding blocks of all local blocks (must not precede
!-----setting of grid_xmin,... etc)

      Call find_surrblks()
!      print *,'morton 3'

!-----Create guardcell filling communications information
      tag_offset = 100
      Call mpi_morton_bnd(mype,nprocs,tag_offset)
!      print *,'morton 4'
!-----Create prolongation communications information
      tag_offset = 100
      Call mpi_morton_bnd_prolong(mype,nprocs,tag_offset)
!      print *,'morton 5'
!-----Create flux fix and edge fix communications information
      tag_offset = 100
      Call mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)
!      print *,'morton 6'
!-----Create restriction communications information
      tag_offset = 100
      Call mpi_morton_bnd_restrict(mype,nprocs,tag_offset)
!      print *,'morton 7'

      Return
      End Subroutine amr_morton_process

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine pf_morton_process()

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use io

      Use paramesh_mpi_interfaces, Only : mpi_morton_bnd,              & 
                                          mpi_morton_bnd_prolong,      & 
                                          mpi_morton_bnd_fluxcon,      & 
                                          mpi_morton_bnd_restrict,     & 
                                          mpi_setup

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Local Variables
      Integer :: nprocs,mype,tag_offset,ierr

!-----Begin Executable code.

!      print *,'morton mpi'
      Call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      Call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

!-----call setup routines in preparation for calling 
!-----all the mpi_morton_bnd_XXX routines.

!-----Find the coordinate ranges
!      print *,'morton 0'
!      Call mpi_amr_global_domain_limits

!      print *,'morton 1', mype
!      Call mpi_setup(mype,nprocs)
      Call pf_mpi_setup(mype,nprocs)
!      print *,'morton 2', mype

!-----Set up surrounding blocks of all local blocks (must not precede
!-----setting of grid_xmin,... etc)

      Call find_surrblks()
!      print *,'morton 3', mype

!-----Create guardcell filling communications information
      tag_offset = 100
!      Call mpi_morton_bnd(mype,nprocs,tag_offset)
      Call pf_morton_bnd(mype,nprocs,tag_offset)
!      print *,'morton 4', mype
!-----Create prolongation communications information
      tag_offset = 100
!      Call mpi_morton_bnd_prolong(mype,nprocs,tag_offset)
      Call pf_morton_bnd_prolong(mype,nprocs,tag_offset)
!      print *,'morton 5', mype
!-----Create flux fix and edge fix communications information
      tag_offset = 100
      Call mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)
!      print *,'morton 6', mype
!-----Create restriction communications information
      tag_offset = 100
!      Call mpi_morton_bnd_restrict(mype,nprocs,tag_offset)
      Call pf_morton_bnd_restrict(mype,nprocs,tag_offset)
!      print *,'morton 7', mype

      Return
      End Subroutine pf_morton_process
