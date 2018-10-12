!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_set_runtime_parameters
!! NAME
!!
!!   amr_set_runtime_parameters
!!
!! SYNOPSIS
!!
!!   Call amr_set_runtime_parameters()
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
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   io
!!
!! CALLS
!!
!! RETURNS
!!
!!   Nothing returned
!!
!! DESCRIPTION
!!
!!   This routine reads in the runtime parameters which are used by
!!   PARAMESH to set array sizes.  The runtime parameters are stored
!!   in the file 'amr_runtime_parameters' which you must create.
!!   A copy of the file 'amr_runtime_parameters' must be available for 
!!   each processor to open and then to read from.
!!
!! AUTHORS
!!
!!   Kevin Olson
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

        Subroutine amr_set_runtime_parameters()

!-------Use statements.
        Use paramesh_dimensions
        Use physicaldata
        Use tree
        Use timings
        Use io

        Implicit None

!-------Include statments
        Include 'mpif.h'

!-------Local variables
        Integer :: ierr
        Integer :: mype, npes, iproc

!-------Begin executable code.

        Call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
        Call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierr)

        Do iproc = 0, npes-1
        If (iproc == mype) Then

        Open (unit=35,                                                 & 
              file='amr_runtime_parameters',                           & 
              status='old',                                            & 
              action='READ',                                           & 
              form='formatted')
!-------Integers
        Read (35,*) maxblocks
        Read (35,*) ndim
        Read (35,*) l2p5d
        Read (35,*) nxb
        Read (35,*) nyb
        Read (35,*) nzb
        Read (35,*) nvar
        Read (35,*) nfacevar
        Read (35,*) nvaredge
        Read (35,*) nvarcorn
        Read (35,*) nvar_work
        Read (35,*) nguard
        Read (35,*) nguard_work
        Read (35,*) nfluxvar
        Read (35,*) nedgevar1
        Read (35,*) iface_off
        Read (35,*) mflags
        Read (35,*) nfield_divf
        Read (35,*) nboundaries
!-------Logicals
        Read (35,*) diagonals
        Read (35,*) amr_error_checking
        Read (35,*) no_permanent_guardcells
        Read (35,*) advance_all_levels
        Read (35,*) force_consistency
        Read (35,*) consv_fluxes
        Read (35,*) consv_flux_densities
        Read (35,*) edge_value
        Read (35,*) edge_value_integ
        Read (35,*) var_dt
        Read (35,*) pred_corr
        Read (35,*) empty_cells
        Read (35,*) conserve
        Read (35,*) divergence_free
        Read (35,*) curvilinear
        Read (35,*) curvilinear_conserve
        Read (35,*) cartesian_pm
        Read (35,*) cylindrical_pm
        Read (35,*) spherical_pm
        Read (35,*) polar_pm
        Read (35,*) lsingular_line
        Read (35,*) timing_mpi
        Read (35,*) timing_mpix
!-------characters
        Read (35,*) output_dir

        Close(35)

        End If  ! End If (iproc == mype)

        Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        End Do  ! End Do iproc =

        amr_log_file = trim(output_dir) // 'amr.log'

        Return
        End Subroutine amr_set_runtime_parameters
