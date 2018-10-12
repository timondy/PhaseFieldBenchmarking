!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! amr_initial_soln
!!!! * Stub that calls user specific function at a blockwise level
!!!!      -> app_initial_soln_blk
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_initial_soln
  
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  
  integer :: lb

!--------------------------------------------------------------

! loop over leaf grid blocks
  if(lnblocks.gt.0) then
!!!$OMP PARALLEL PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,100) 
     do lb=1,lnblocks
        if(nodetype(lb).eq.1 .or. advance_all_levels) then
              call app_initial_soln_blk(lb)
        endif
     enddo ! end loop over grid blocks
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
  endif

  call app_initial_soln_end

  return
end subroutine amr_initial_soln

