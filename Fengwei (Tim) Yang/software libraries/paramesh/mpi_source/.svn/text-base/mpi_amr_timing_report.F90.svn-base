!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2002
! By: Peter J. Macneice, Drexell University.
!     Kevin M. Olson, Univ. of MD Baltimore Campus.
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


        subroutine amr_timing_init

!----------------------------------------------------------------
! include file to define physical qualities of the model and mesh
        use paramesh_dimensions
        use physicaldata

! include file defining the tree
        use tree

        use timings

! identify interface blocks for routines called in this program
        use paramesh_interfaces, only : comm_start, & 
     &                                  amr_initialize, & 
     &                                  amr_refine_derefine, & 
     &                                  amr_guardcell, & 
     &                                  amr_prolong, & 
     &                                  amr_restrict, & 
     &                                  amr_test_refinement, & 
     &                                  amr_close

        implicit none

        include 'mpif.h'

! local amr variables
        integer           :: nprocs, mype
        integer           :: minstp, maxstp, ngcell_on_cc_l

        double precision  :: gno_of_flops_advance
        double precision  :: time1

!---------------------------------------------------------------

        if (timing_mpi) then
        timer_amr_initialize = 0.
        timer_amr_refine_derefine = 0.
        timer_amr_check_derefine = 0.
        timer_amr_check_refine = 0.
        timer_amr_refine_blocks = 0.
        timer_amr_derefine_blocks = 0.
        timer_amr_morton_order = 0.
        timer_amr_morton_process = 0.
          timer_amr_global_domain_limits = 0.
          timer_mort_comm_for_surrblks = 0.
          timer_mpi_setup = 0.
          timer_amr_morton_limits = 0.
          timer_mort_comm_for_surrblks = 0.

        timer_mpi_morton_bnd = 0.
        timer_mpi_morton_bnd3= 0.
        timer_mpi_morton_bnd_prolong1 = 0.
        timer_mpi_morton_bnd_fluxcon = 0.
        timer_mpi_morton_bnd_restrict = 0.

        timer_amr_boundary_block_info = 0.

        timer_amr_1blk_restrict = 0.

        timer_amr_guardcell = 0.
        timer_amr_1blk_guardcell = 0.
        timer_amr_1blk_guardcell_srl = 0.
        timer_amr_1blk_guardcell_c_to_f = 0.

        timer_amr_1blk_cc_cp_remote = 0.

        timer_amr_1blk_copy_soln = 0.
        timer_amr_1blk_to_perm = 0.

        timer_amr_comm_setup = 0.
        timer_amr_1blk_cc_prol_gen_unk = 0.
        timer_amr_1blk_cc_prol_gen_work = 0.
        timer_amr_1blk_fc_prol_gen = 0.
        timer_amr_1blk_ec_prol_gen = 0.
        timer_amr_1blk_nc_prol_gen = 0.

        timer_amr_prolong = 0.

        timer_amr_restrict = 0.
        timer_amr_test_refinement = 0.
        timer_advance_soln = 0.

        no_of_flops_advance = 0.

        no_of_calls_check_derefine = 0
        mess_counter_chk_deref = 0

        start_time = mpi_wtime()

        endif

        return
        end subroutine amr_timing_init




      subroutine amr_timing_report(minstp, maxstp)




!---------------------------------------------------------------
!
! Ouput timings
!
!---------------------------------------------------------------
! include file to define physical qualities of the model and mesh
        use paramesh_dimensions
        use physicaldata

! include file defining the tree
        use tree

        use timings

! identify interface blocks for routines called in this program
        use paramesh_interfaces, only : comm_start, & 
     &                                  amr_initialize, & 
     &                                  amr_refine_derefine, & 
     &                                  amr_guardcell, & 
     &                                  amr_prolong, & 
     &                                  amr_restrict, & 
     &                                  amr_test_refinement, & 
     &                                  amr_close

        implicit none

        integer, intent(in) :: minstp, maxstp

        include 'mpif.h'

! local amr variables
        integer           :: nprocs, mype, tot_blocks, ierr
        integer           :: i,j,k,l
        integer           :: io, ivar, ngcell_on_cc_l
        double precision  :: gno_of_flops_advance
        double precision  :: time1

!---------------------------------------------------------------
      if (timing_mpi) then

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)


! Sum flops
      call comm_int_sum_to_all(tot_blocks,lnblocks)
      write(*,*)'no of flops in advance = ',no_of_flops_advance
      call comm_dble_sum_to_all(gno_of_flops_advance, & 
     &                            no_of_flops_advance)

! Record no. of cell centered variables n guardcell filling
      ngcell_on_cc_l = 0
      do ivar = 1,nvar
        if(gcell_on_cc(ivar)) ngcell_on_cc_l = ngcell_on_cc_l + 1
      enddo


      write(*,*) 'pe ',mype,' lnblocks ',lnblocks

      call amr_flush(6)
!      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      write(*,*) 'CEG in mpi_amr_timing_report'

      if(mype.eq.0) then

      timer_amr_1blk_to_perm(0) = 0.
      do i=1,1+nvar_work
        timer_amr_1blk_to_perm(0) = timer_amr_1blk_to_perm(0) + & 
     &                              timer_amr_1blk_to_perm(i)
      enddo


      io = 11
      open(unit=io,status='unknown',position='append', & 
     &                              file='profile.dat')

      end_time = mpi_wtime()
      write(io,*)' '
      write(io,*)' '
      write(io,*) '===================================================='
      write(io,*) ' Performance measurements'
      write(io,*) '===================================================='
      write(io,*) 'nprocs = ',nprocs
      write(io,*) 'ndim   = ',ndim
      write(io,*) 'iteration range  = ',minstp,maxstp
      write(io,*) 'no_of_calls_check_derefine = ', & 
     &                          no_of_calls_check_derefine
      if (no_permanent_guardcells) then
      write(io,*) 'no_permanent_guardcell = true'
      else
      write(io,*) 'no_permanent_guardcell = false'
      endif 
      if (advance_all_levels) then
         write(io,*) 'advance_all_levels     = true'
      else
         write(io,*) 'advance_all_levels     = false'
      end if
      write(io,*)'nvar       = ',nvar
      write(io,*)'no of cc variables advanced = ',ngcell_on_cc_l
      write(io,*)'nfacevar   = ',nfacevar
      write(io,*)'nvaredge   = ',nvaredge
      write(io,*)'nvarcorn   = ',nvarcorn
      write(io,*)'total blocks = ',tot_blocks
      write(io,*)'blocks per pe = ',real(tot_blocks)/real(nprocs)
      write(io,*)'no of flops in advance = ',gno_of_flops_advance
      write(io,*) '===================================================='
      write(io,*) 'Total time = ',end_time-start_time
      write(io,*)'timer_amr_initialize      = ',timer_amr_initialize
      write(io,*)' '
      write(io,*)'timer_amr_refine_derefine = ', & 
     &                                     timer_amr_refine_derefine
      write(io,*)'  timer_amr_check_derefine(0)        = ', & 
     &                                  timer_amr_check_derefine(0)
      do i = 1,4
      write(io,*)'    timer_amr_check_derefine(',i,')  = ', & 
     &                                  timer_amr_check_derefine(i)
      enddo
      write(io,*)'    mess_counter_chk_deref = ', & 
     &                                  mess_counter_chk_deref
      write(io,*)'  timer_amr_check_refine    = ', & 
     &                                     timer_amr_check_refine
      write(io,*)'  timer_amr_refine_blocks = ', & 
     &                                     timer_amr_refine_blocks
      write(io,*)'  timer_amr_derefine_blocks = ', & 
     &                                     timer_amr_derefine_blocks
      write(io,*)'  timer_amr_morton_order    = ', & 
     &                                     timer_amr_morton_order
      write(io,*)'  timer_amr_morton_process  = ', & 
     &                                     timer_amr_morton_process
      write(io,*)'    timer_amr_global_domain_limits  = ', & 
     &                                   timer_amr_global_domain_limits
      write(io,*)'    timer_mort_comm_for_surrblks  = ', & 
     &                                     timer_mort_comm_for_surrblks
      write(io,*)'    timer_mpi_setup               = ', & 
     &                                     timer_mpi_setup
      write(io,*)'    timer_amr_morton_limits       = ', & 
     &                                     timer_amr_morton_limits

! morton bnd routine timings
      write(io,*)'timer_mpi_morton_bnd (total) = ', & 
     &                               timer_mpi_morton_bnd(1,0)
      do i = 1,20
        write(io,*)'  timer_mpi_morton_bnd(1,',i,') = ', & 
     &                               timer_mpi_morton_bnd(1,i)
        if(i.eq.3) then
          do j = 1,7
            write(io,*)'    timer_mpi_morton_bnd3(1,',j,') = ', & 
     &                               timer_mpi_morton_bnd3(1,j)
          enddo
        endif
      enddo
      write(io,*)'timer_mpi_morton_bnd_prolong1 (total) = ', & 
     &                               timer_mpi_morton_bnd(2,0)
      do i = 1,20
        write(io,*)'  timer_mpi_morton_bnd(2,',i,') = ', & 
     &                               timer_mpi_morton_bnd(2,i)
        if(i.eq.3) then
          do j = 1,7
            write(io,*)'    timer_mpi_morton_bnd3(2,',j,') = ', & 
     &                               timer_mpi_morton_bnd3(2,j)
          enddo
        endif
      enddo
      write(io,*)'timer_mpi_morton_bnd_restrict (total) = ', & 
     &                               timer_mpi_morton_bnd(3,0)
      do i = 1,20
        write(io,*)'  timer_mpi_morton_bnd(3,',i,') = ', & 
     &                               timer_mpi_morton_bnd(3,i)
        if(i.eq.3) then
          do j = 1,7
            write(io,*)'    timer_mpi_morton_bnd3(3,',j,') = ', & 
     &                               timer_mpi_morton_bnd3(3,j)
          enddo
        endif
      enddo
      write(io,*)'timer_mpi_morton_bnd_fluxcon (total) = ', & 
     &                               timer_mpi_morton_bnd(4,0)
      do i = 1,20
        write(io,*)'  timer_mpi_morton_bnd(4,',i,') = ', & 
     &                               timer_mpi_morton_bnd(4,i)
        if(i.eq.3) then
          do j = 1,7
            write(io,*)'    timer_mpi_morton_bnd3(4,',j,') = ', & 
     &                               timer_mpi_morton_bnd3(4,j)
          enddo
        endif
      enddo


      write(io,*)' '
      write(io,*)'timer_amr_guardcell       = ',timer_amr_guardcell
      write(io,*)' '
      write(io,*)'timer_amr_1blk_guardcell(0) = ', & 
     &                                timer_amr_1blk_guardcell(0)
      write(io,*)'  timer_amr_1blk_guardcell(setup) = ', & 
     &                                timer_amr_1blk_guardcell(1)
      write(io,*)'  timer_amr_1blk_guardcell(coarse neigh) = ', & 
     &                                timer_amr_1blk_guardcell(2)
      write(io,*)'  timer_amr_1blk_guardcell(local srl call) = ', & 
     &                                timer_amr_1blk_guardcell(3)
      write(io,*)'  timer_amr_1blk_guardcell_srl= ', & 
     &                               timer_amr_1blk_guardcell_srl
      write(io,*)'    timer_amr_1blk_cc_cp_remote(0)= ', & 
     &                               timer_amr_1blk_cc_cp_remote(0)
      write(io,*)'      timer_amr_1blk_cc_cp_remote(1)= ', & 
     &                               timer_amr_1blk_cc_cp_remote(1)
      write(io,*)'      timer_amr_1blk_cc_cp_remote(2)= ', & 
     &                               timer_amr_1blk_cc_cp_remote(2)
      write(io,*)'      timer_amr_1blk_cc_cp_remote(3)= ', & 
     &                               timer_amr_1blk_cc_cp_remote(3)
      write(io,*)'  timer_amr_1blk_guardcell_c_to_f= ', & 
     &                               timer_amr_1blk_guardcell_c_to_f
      write(io,*)' '
      write(io,*)'timer_amr_1blk_copy_soln  = ',timer_amr_1blk_copy_soln
      write(io,*)'timer_amr_1blk_to_perm(0) = ', & 
     &                                         timer_amr_1blk_to_perm(0)
      do i=1,1+nvar_work
        write(io,*)'  timer_amr_1blk_to_perm(',i,') = ', & 
     &                                         timer_amr_1blk_to_perm(i)
      enddo
      write(io,*)'timer_amr_comm_setup(0)   = ',timer_amr_comm_setup(0)
      write(io,*)'  timer_amr_comm_setup(init        )   = ', & 
     &                                          timer_amr_comm_setup(1)
      write(io,*)'  timer_amr_comm_setup(read cntrl  )   = ', & 
     &                                          timer_amr_comm_setup(2)
      write(io,*)'  timer_amr_comm_setup(set buf size)   = ', & 
     &                                          timer_amr_comm_setup(3)
      write(io,*)'  timer_amr_comm_setup(pack        )   = ', & 
     &                                          timer_amr_comm_setup(4)
      write(io,*)'  timer_amr_comm_setup(xchange     )   = ', & 
     &                                          timer_amr_comm_setup(5)
      write(io,*)'  timer_amr_comm_setup(unpack      )   = ', & 
     &                                          timer_amr_comm_setup(6)
      write(io,*)'  timer_amr_comm_setup(set tr buf size)   = ', & 
     &                                          timer_amr_comm_setup(7)
      write(io,*)'  timer_amr_comm_setup(tree pack  )   = ', & 
     &                                          timer_amr_comm_setup(8)
      write(io,*)'  timer_amr_comm_setup(tree unpack)   = ', & 
     &                                          timer_amr_comm_setup(9)
      write(io,*)'timer_amr_1blk_cc_prol_gen_unk = ', & 
     &                               timer_amr_1blk_cc_prol_gen_unk
      write(io,*)'timer_amr_1blk_cc_prol_gen_work = ', & 
     &                               timer_amr_1blk_cc_prol_gen_work
      write(io,*)'timer_amr_1blk_fc_prol_gen = ', & 
     &                               timer_amr_1blk_fc_prol_gen
      write(io,*)'timer_amr_1blk_ec_prol_gen = ', & 
     &                               timer_amr_1blk_ec_prol_gen
      write(io,*)'timer_amr_1blk_nc_prol_gen = ', & 
     &                               timer_amr_1blk_nc_prol_gen
      write(io,*)' '
      write(io,*)'timer_amr_prolong         = ',timer_amr_prolong
      write(io,*)'timer_amr_restrict        = ',timer_amr_restrict
      write(io,*)'  timer_amr_1blk_restrict   = ', & 
     &                                           timer_amr_1blk_restrict
      write(io,*)'timer_amr_test_refinement(0) = ', & 
     &                                      timer_amr_test_refinement(0)
      write(io,*)'  timer_amr_test_refinement(copy to work  ) = ', & 
     &                                      timer_amr_test_refinement(1)
      write(io,*)'  timer_amr_test_refinement(comm setup    ) = ', & 
     &                                      timer_amr_test_refinement(2)
      write(io,*)'  timer_amr_test_refinement(1blk_guardcell) = ', & 
     &                                      timer_amr_test_refinement(3)
      write(io,*)'  timer_amr_test_refinement(error eval    ) = ', & 
     &                                      timer_amr_test_refinement(4)
      write(io,*)'timer_advance_soln(0)         = ', & 
     &                                             timer_advance_soln(0)
      write(io,*)'  timer_advance_soln(timestep      )= ', & 
     &                                             timer_advance_soln(1)
      write(io,*)'  timer_advance_soln(copy to 1blk  )= ', & 
     &                                             timer_advance_soln(2)
      write(io,*)'  timer_advance_soln(comm setup    )= ', & 
     &                                             timer_advance_soln(3)
      write(io,*)'  timer_advance_soln(1blk guardcell)= ', & 
     &                                             timer_advance_soln(4)
      write(io,*)'  timer_advance_soln(advance       )= ', & 
     &                                             timer_advance_soln(5)
      write(io,*)'  timer_advance_soln(copy to perm  )= ', & 
     &                                             timer_advance_soln(6)
      write(io,*) '===================================================='

      close(unit=io)

      endif              ! pe==0
!---------------------------------------------------------------

      endif


      return
      end subroutine amr_timing_report
