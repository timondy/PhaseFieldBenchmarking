!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! timings Module
!------------------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Module timings

        Use paramesh_dimensions

        Double Precision,Public :: timer_amr_initialize
        Double Precision,Public :: start_time, end_time


        Double Precision,Public :: timer_amr_refine_derefine
        Double Precision,Public :: timer_amr_check_derefine(0:4)
        Double Precision,Public :: timer_amr_check_refine
        Double Precision,Public :: timer_amr_refine_blocks
        Double Precision,Public :: timer_amr_derefine_blocks
        Double Precision,Public :: timer_amr_morton_order
        Double Precision,Public :: timer_amr_morton_process
        Double Precision,Public :: timer_amr_boundary_block_info
        Double Precision,Public :: timer_amr_global_domain_limits
        Double Precision,Public :: timer_mort_comm_for_surrblks
        Double Precision,Public :: timer_mpi_setup
        Double Precision,Public :: timer_amr_morton_limits
        Double Precision,Public :: timer_mpi_morton_bnd(1:4,0:20)
        Double Precision,Public :: timer_mpi_morton_bnd3(1:4,1:7)
        Double Precision,Public :: timer_mpi_morton_bnd_prolong1
        Double Precision,Public :: timer_mpi_morton_bnd_fluxcon
        Double Precision,Public :: timer_mpi_morton_bnd_restrict

        Double Precision,Public :: timer_amr_guardcell
        Double Precision,Public :: timer_amr_1blk_guardcell(0:3)
        Double Precision,Public :: timer_amr_1blk_guardcell_c_to_f
        Double Precision,Public :: timer_amr_1blk_guardcell_srl

        Double Precision,Public :: timer_amr_1blk_cc_cp_remote(0:3)

        Double Precision,Public :: timer_amr_1blk_copy_soln
        Double Precision,Allocatable,Public :: timer_amr_1blk_to_perm(:)
        Double Precision,Public :: timer_amr_comm_setup(0:9)
        Double Precision,Public :: timer_amr_1blk_cc_prol_gen_unk
        Double Precision,Public :: timer_amr_1blk_cc_prol_gen_work
        Double Precision,Public :: timer_amr_1blk_fc_prol_gen
        Double Precision,Public :: timer_amr_1blk_ec_prol_gen
        Double Precision,Public :: timer_amr_1blk_nc_prol_gen

        Double Precision,Public :: timer_amr_prolong
        Double Precision,Public :: timer_amr_restrict
        Double Precision,Public :: timer_amr_1blk_restrict
        Double Precision,Public :: timer_amr_test_refinement(0:4)
        Double Precision,Public :: timer_advance_soln(0:6)
        Double Precision,Public :: no_of_flops_advance

        Integer,Public,Parameter  :: addflops = 1
        Integer,Public,Parameter  :: mulflops = 1
        Integer,Public,Parameter  :: divflops = 1

        Integer,Public          :: no_of_calls_check_derefine
        Integer,Public          :: mess_counter_chk_deref


! timing_mpi flag
      Public :: timing_mpi
      Logical, Save :: timing_mpi

! timing_mpix flag
      Public :: timing_mpix
      Logical, Save :: timing_mpix

      End Module timings
!-----------------------------------------------------------------
