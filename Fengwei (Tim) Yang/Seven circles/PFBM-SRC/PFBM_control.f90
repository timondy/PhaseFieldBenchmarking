!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_read_checkpoint                                          REQUIRED
!!!!  * Calls amr_checkpoint_re with application specific extra parameters
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_pretimeloop_output                                       REQUIRED
!!!!  * Any processing or output needed before first timestep
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_close_files                                              REQUIRED
!!!!  * Closes any files opened in application section
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_read_checkpoint(pid, noprocs)
  use paramesh_interfaces
  use time_dep_parameters
  use checkpoint_parameters

  use paramesh_dimensions
  use tree
  use multigrid_parameters
  use physicaldata
  implicit none
  integer,intent(in) :: pid, noprocs
  
  if(pid.eq.0) write(6,*) 'Checkpoint on : ', chk_chkpt
  call amr_checkpoint_re(chk_chkpt, user_attr_1=chk_dt, user_attr_2=chk_t, user_attr_3=chk_dtold)
  if(pid.eq.0) write(6,*) 'Checkpoint read with ', chk_t, chk_dt, chk_dtold

  !!--tim in most of cases of re-run failing program, time-step should have
  !reached 20 or more, we reserve file names from 1 - 19 for
  !initial condition files.
!!  chk_dt = 0.2/50.0

!!  mode_1 = 2
  chk_dtold = chk_dt
  if(chk_chkpt.lt.20)then
    chk_chkpt = 0
    mode_1 = 1
  endif
end subroutine app_read_checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_pretimeloop_output(pid, noprocs)
!!--tim  use multigrid_parameters  ! for min_dx, total_vars
!!--tim  use solution_parameters   ! for most of the parameters
!!--tim  use tree                  ! for grid_xmax
!!--tim  use paramesh_dimensions   ! for nxb
!!--tim  use checkpoint_parameters ! for chk_chkpt
!!--tim  implicit none
!!--tim  integer,intent(in) :: pid, noprocs

  !CEG Add README output
!!--tim  if(pid.eq.0.and.chk_chkpt.eq.1) then
!!--tim     open (unit=61,file="README")
!--tim     write(61,*) 'Parameter list:'
 !--tim    write(61,*) 'Mc-Inf      ',  mcinf
 !--tim    write(61,*) 'kappa_E     ',  ke
 !--tim    write(61,*) 'epsilon     ',  epsilon
 !--tim    write(61,*) 'lambda      ',  lambda
 !--tim    write(61,*) 'D_solute    ', D_solute
 !--tim    write(61,*) 'D_therm     ', D_therm
 !--tim    write(61,*) 'Lewis no.   ',  le
 !--tim    write(61,*) 'Delta       ',  delta
 !--tim    write(61,*) 'nucleus rad ', nuc_radius
!     write(61,*) 'End time    ',  end_time
!!--tim     write(61,*) 'Domain size ', grid_xmax
!!--tim     write(61,*) 'Finest dx   ', min_dx
!!--tim     write(61,*) 'Block size  ', nxb
!!--tim     write(61,*) 'total_vars  ', total_vars
 !--tim    write(61,*) 'Solute on?  ', solute
!!--tim     write(61,*) ' '
!!--tim     flush(61)
!!--tim  endif

!!--tim  if(pid.eq.0) then
 !--tim    print *,"Lambda: ",lambda,"Dtherm",D_therm,"Dsolute",D_solute
!!--tim     open (unit=125,file="tip_rad.txt")
!!--tim     open (unit=126,file="tip_loc.txt")
!!--tim     open (unit=133,file="real_time.txt")
     !open (unit=301,file="Vcycnt.txt")
!!--tim  endif


end subroutine app_pretimeloop_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_close_files(pid, noprocs)
 implicit none
 integer,intent(in) :: pid, noprocs

!!--tim  if(pid.eq.0) then
!!--tim     close(61)
!!--tim     close(125)
!!--tim     close(126)
!!--tim     close(133)
!!--tim  endif

end subroutine app_close_files  

