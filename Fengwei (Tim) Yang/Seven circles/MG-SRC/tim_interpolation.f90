subroutine tim_interpolation(pid, noprocs)
  use paramesh_interfaces
  use paramesh_dimensions
  use tree
  use time_dep_parameters
  use multigrid_parameters
  use generic_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer ::  lb, loop_count, iopt, nlayers

  iopt = 1
  nlayers = nguard
  
  ! Initialise everything
  refine(1:lnblocks) = .false.
  call amr_refine_derefine

  lrefine_max = 7
  do loop_count=desired_lvl, 7
    if(pid.eq.0.and.verbose.ge.2) print *,'Pre refinement', loop_count, lnblocks
    call pf_guardcell(pid,iopt,loop_count)
!!    call amr_initial_soln
!!    derefine(1:lnblocks) = .false.
    write(*,*) "lnblocks: ", lnblocks, loop_count

    do lb = 1, lnblocks
      if(nodetype(lb).eq.1) refine(lb) = .true.
    enddo
 
    ! refine grid
    call amr_refine_derefine
    call amr_prolong(pid,iopt,nlayers)
    call pf_guardcell(pid,iopt,loop_count+1)


!!    call amr_mg_get_max_refinement(pid,noprocs)
!!    call amr_multigrid_child_set()
!!    call reset_nodetype(noprocs,pid,loop_count)

!!    call pf_guardcell(pid,iopt,loop_count+1)

  enddo
!!    call amr_mg_init()


end subroutine tim_interpolation

