!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  amr_mg_control
!!!!  * Shell for calling main time-stepping code
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  amr_mg_time_dep_control
!!!!  * Main time-stepping loop
!!!!  * Calls various application specific functions:
!!!!       -> app_pretimeloop_output
!!!!       -> app_output_perstep
!!!!       -> app_output_occasional
!!!!       -> app_close_files
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  get_input_parameters
!!!!  * Processes command line inputs to read checkpoint data
!!!!  * Calls application specific function
!!!!          app_read_checkpoint
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  timestep_adjust
!!!!  * Changes timestep size
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!!  nodebound
!!!!  * [unused] constraint on very small values in unk() arrays
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine amr_mg_control(pid,noprocs)
! Wrapper for controlling whole solver process
  integer, intent(in) :: pid, noprocs

  call amr_mg_time_dep_control(pid,noprocs)

end subroutine amr_mg_control

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_time_dep_control(pid,noprocs)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_interfaces
  use checkpoint_parameters
  use generic_parameters
  use paramesh_mpi_interfaces
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer :: back_step_ct,time_step,i_step,i_step_total,ierr,level,iopt, lb, npts
  double precision defect,end_time,old_defect_1, rmsres(total_vars)
  integer failer, finegridblocks(2)
  double precision :: timeprev, timenew
  double precision local_r
  double precision :: wtime, wtimeprev
  double precision :: defect_1
  integer :: i,j,k,g,nguard0, lvl, input_lvl, nlayers
  double precision :: defect_defect, convrate,total_convrate
  double precision :: dx, xi, yi
 
  failer=0
  rmsres(:) = 0.0
  npts = 0

  end_time=time+simulation_time
  incset=0
  if(pid.eq.0.and.verbose.ge.2) write(6,*) 'end time is ',end_time

!!  call tim_interpolation_convergence(pid, noprocs,lrefine_max, 7)
!!  stop

! Pre timestep output
  
  call app_pretimeloop_output(pid, noprocs)

  if (verbose.ge.2.and.chk_chkpt.gt.1) then
     if (pid.eq.0) write(6,*) '----------Loaded solution--------------'
     call app_output_perstep(pid, noprocs, chk_chkpt-1)
     if (pid.eq.0) write(6,*) '---------------------------------------'
  endif

  if(pid.eq.0.and.verbose.ge.2)print *,"Entered multigrid control"
  iopt=1
  back_step_ct=0
  allocate(has_children(1:maxblocks),stat=time_step)!Using time_step as dud var
  allocate(nodetype_copy(1:maxblocks),stat=time_step)!Using time_step as dud var
!  call amr_guardcell(pid,iopt,nguard)
  !has_children(:) = .false.
  
!  call pf_guardcell(pid,iopt,mg_max_level)
  call amr_mg_init()
  call amr_multigrid_block_types(pid,noprocs,finegridblocks)
  
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!  call amr_guardcell(pid,iopt,nguard)
!  call pf_guardcell(pid,iopt,mg_max_level)

  if(pid.eq.0.and.verbose.ge.1)print *,"Number of 'real' variables:",total_vars
  
  i_step_total=0
  if(pid.eq.0.and.verbose.ge.2)print *,"Starting timestepping"
!  call ceg_block_balancing(pid, noprocs)
!  call amr_multigrid_block_types(pid,noprocs,finegridblocks)
!  stop
  wtimeprev = MPI_Wtime()


!!  call tim_visual_paraview(chk_chkpt)
!!  call app_output_occasional(pid, noprocs, lrefine_max)
!!  stop


!!  call tim_para_csv(pid, noprocs, 0)
!!  call tim_visual_paraview(chk_chkpt-1)

  do time_step=chk_chkpt,max_time_iter+back_step_ct


    if(time_step.eq.1) mode_1 = 1
    ! Little bit of code to ensure the last time step ends on simulation_time
    !!--modified by tim
    !!--may cause tiny dt to be taken on the very last time step.
    if(time.ge.end_time)exit
    if((time+dt).gt.end_time)then
      dt=end_time-time
      if(dt.lt.1.e-16)then
        exit
      endif
    end if
    time=time+dt

    if(pid.eq.0.and.verbose.ge.2)print *,"max_time_iter, back_step_ct:: ", max_time_iter,back_step_ct
    nodetype_copy(1:lnblocks)=nodetype(1:lnblocks)
    call amr_multigrid_child_set()
    if(pid.eq.0.and.verbose.ge.1)print *,""!Blank line between timesteps
!!new    if(pid.eq.0.and.verbose.ge.0)write(6,'(A,I6,A,E12.4,A,E12.4)'),"Starting time step no.",time_step-back_step_ct,".  Time:",time-dt," to ",time
    if(pid.eq.0.and.verbose.ge.0)write(6,*),"No.",time_step-back_step_ct,".  Time:",time-dt," to ",time
    if(pid.eq.0.and.verbose.ge.1)write(*,'(A,f12.4,A,f12.4)')"Current dt:",dt,"  smallest h:",bsize(1,lnblocks)/real(nxb)
!!    if(pid.eq.0.and.verbose.ge.1)write(*,'(A,f12.4,A,f12.4,A,f12.4)')"interface:",ep,'  points-in-region:',ep/(bsize(1,lnblocks)/real(nxb))
    if(pid.eq.0.and.verbose.ge.1)print *,"mode_1:",mode_1
 
    call amr_mg_get_max_refinement(pid,noprocs)
    call reset_nodetype(noprocs,pid,mg_max_level)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Adaptivity tests
    if(local_adaptation.ge.1.and.failer.eq.0.and.chk_restart.ne.1)then
    ! DO NOT DO THIS CALL IF grow=1
    !call amr_guardcell(pid,iopt,nguard)

      ! Just do main variable guardcells
      gcell_on_cc(:)=.false.
      do current_var=1,total_vars
        current_unk=1+(current_var-1)*nunkvbles
        gcell_on_cc(current_unk)=.true.
        current_unk=4+(current_var-1)*nunkvbles
        gcell_on_cc(current_unk)=.true.
        current_unk=5+(current_var-1)*nunkvbles
        gcell_on_cc(current_unk)=.true.
      end do

      ! Loop over grids to match guardcells
      do g = lrefine_max, 2, -1
        call reset_nodetype(noprocs,pid,mg_max_level)
!       call reset_nodetype(noprocs,pid,g)
!       call pf_guardcell(pid,1,g)

        call pf_restrict(pid, 1, 0, .true., g)

        call reset_nodetype(noprocs, pid, g-1)
        call pf_guardcell(pid,iopt,g-1)
      end do
      call reset_nodetype(noprocs, pid, mg_max_level)

      if(pid.eq.0.and.verbose.ge.1)print *,pid,"Testing refinement"
      if(local_adaptation.eq.1) then
           ! Default adaptation
        call amr_test_refinement(pid, lrefine_min, lrefine_max)
      else
           ! User supplied adaptation
        call app_test_refinement(pid, lrefine_min, lrefine_max)
      endif
      if(pid.eq.0.and.verbose.ge.2)print *,"Refining"
      if(pid.eq.0.and.verbose.ge.2) call system('date +"%k:%M:%S:%N"')
      if(pid.eq.0) wtimeprev = MPI_Wtime()
      if(pid.eq.0) then
        call cpu_time(timenew)
        timeprev=timenew
      endif
        
      call amr_refine_derefine
                        
      if(pid.eq.0) then
         wtime = MPI_Wtime()
         call cpu_time(timenew)
         if (verbose.ge.2) print*,'Sys/Real adapt time was ', timenew-timeprev, wtime-wtimeprev
         timeprev=timenew
         wtimeprev=wtime
      endif
      if(pid.eq.0.and.verbose.ge.2) call system('date +"%k:%M:%S:%N"')

!     call amr_multigrid_block_types(pid,noprocs,finegridblocks)
      call amr_prolong(pid,iopt,nguard)
     
      ! CEG try domain growing
      if (grow.eq.1) call amr_grow_domain(pid, noprocs, time_step-back_step_ct)


      !call amr_multigrid_block_types(pid,noprocs)
      if(pid.eq.0.and.verbose.ge.2)print *,"Reinitialising multigrid"
      call amr_mg_init() !NOTE FOR CHRIS, THIS IS NEW THING WHICH STOPPED SEG FAULT
      nodetype_copy(1:lnblocks)=nodetype(1:lnblocks)
      call amr_mg_get_max_refinement(pid,noprocs)
      call amr_multigrid_child_set()
      call reset_nodetype(noprocs,pid,mg_max_level)
      !call amr_mg_init() !NOTE FOR CHRIS, THIS IS OLD THING WHICH NEEDED TO BE MOVED UP
      call amr_multigrid_block_types(pid,noprocs,finegridblocks)
    else if (verbose.ge.3) then
      call amr_multigrid_block_types(pid,noprocs,finegridblocks)
    end if  ! end if test ref

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Copy existing soln(s) to rhs
!    call amr_guardcell(pid,iopt,nguard)
    call pf_guardcell(pid,iopt,mg_max_level)
    call amr_mg_runtime_rhs()
    call pf_guardcell(pid,iopt,mg_max_level)

    if(pid.eq.0.and.verbose.ge.2)print *,"Refinement/guardcell update/RHS setting done, start v-cycles"

    defect=0.0
    old_defect_1=0.0

    if(time_step.eq.1)then
      max_v_cycle_count = 1000
    else
      max_v_cycle_count = 10
    endif

    do i_step=1,max_v_cycle_count
      i_step_total=i_step_total+1
      ! Back up copies of old defects
      old_defect_1=defect
      defect=0.0

      call amr_multigrid_v_cycle(pid,noprocs,mg_max_level,defect, rmsres, npts)
      
      call amr_vcycle_summary(pid, noprocs, i_step, defect, old_defect_1, npts, rmsres, convrate)

      !!--tim
      !!--record the defect right after the very first V-cycle,
      !!--this is then used for stopping criterion.
      if(i_step.eq.1) defect_1 = defect

      !!--tim
      !!--relative stopping criterion.
      if(i_step.gt.2.and.defect.le.defect_1*1.e-8) exit
       
      !!-tim
      !!--absolute stopping criterion.
      if(defect.le.1.e-11) exit
        
    end do

    local_r = -1.0

    ! Timer end code
    if(pid.eq.0) then
      wtime = MPI_Wtime()
      call cpu_time(timenew)
      if (verbose.ge.2) write(6,*)'Sys/Real solve time was ', timenew-timeprev, wtime-wtimeprev
        timeprev=timenew
        wtimeprev=wtime
        if (verbose.ge.4) print *,"Vcycle: ",i_step,"Defect=",defect
!!new        if (verbose.ge.1) write(6,'(A,I7,3X,A,F4.1)')"Total V-cycles used:",i_step_total,"Avg V-cycles per timestep: ",(i_step_total*1.0)/&
!!new              (1.0*(time_step-back_step_ct))
        if (verbose.ge.1) write(6,*)''
!        write (unit=301,fmt=*)time,i_step,(i_step_total*1.0)/(1.0*(time_step-back_step_ct))
    endif
    failer=0

!!--tim     call timestep_adjust(defect, local_r, incset, failer, i_step, pid)

    if (failer.eq.1) then
      i_step_total=i_step_total-i_step
      back_step_ct=back_step_ct+1
      if (dt.lt.min_dt) then
        if(pid.eq.0) print*,'Tiny dt - aborting'
        call MPI_Abort(MPI_COMM_WORLD, failer, ierr)
      endif
    else ! i.e. if (failer.eq.0) then
      ! Output on every step
      call app_output_perstep(pid, noprocs, time_step-back_step_ct)

      if(mod(time_step-back_step_ct,output_rate).eq.0) then
        !Output every output_rate steps
        call app_output_occasional(pid, noprocs, time_step-back_step_ct)

        call cpu_time(timenew)
        if (pid.eq.0) then
          wtime = MPI_Wtime()
          if (verbose.ge.2) print*,'Output took ', timenew-timeprev, wtime-wtimeprev
          timeprev=timenew
        endif
      endif
      ! Set to BDF2 now we have passed the first timestep
      if(time_step.eq.1) mode_1 = 2
!!      if(time_step.ge.2) mode_1 = 3
!!      if(time_step.ge.3) mode_1 = 4
!        call nodebound
    endif
    chk_restart = 0
!!    call tim_true_solution(pid, noprocs, time, 1.5, 0.0)
!!    stop
! #WARNING DEBUGGING get outs on
! DEBUGGING get outs
!     if (mode_1.eq.1.and.dt.gt.0.02) dt = 0.02
!     if (dt.gt.0.1) dt = 0.1
!     if (mode_1.eq.2) return
!     if (mode_1.eq.2) exit
!     if (time_step-back_step_ct.eq.6) return!585
!     if (time_step-back_step_ct.eq.1003) exit!585
!     if (time_step-back_step_ct.eq.10000) return!585


    call tim_evaluation(pid, noprocs, time, time_step)
    
  end do

!!  call tim_visual_paraview(time_step)

!!  call tim_true_solution(pid, noprocs, simulation_time, 1.5, 0.0)
!!  call app_output_occasional(pid, noprocs, time_step-1)
!!  call tim_para_csv(pid, noprocs, time_step-1)
!!  stop

! Close all output text files
!!--tim
!! checkpoint output for final result

!!  call app_output_occasional(pid, noprocs, 10)
!!--tim  write(*,*) "Average convrate =", total_convrate/(i_step_total-time_step-back_step_ct-1)
!!--tim  if(pid.eq.0) write(*,*) file_max_1, file_min_1
!!--tim  call amr_infinity_norm(time, lrefine_max, pid)

  call app_close_files(pid, noprocs)

  !!
!!--tim  if(pid.eq.0)close(101250)
!!--tim  if(pid.eq.0)close(101251)
!!--tim  if(pid.eq.0)close(101252)
!!--tim  if(pid.eq.0)close(101253)

!  if(pid.eq.0)print *,"Mode: ",mode_1
!  if(pid.eq.0)print *,"Min level: ",mg_min_lvl
!  call amr_multigrid_block_types(pid,noprocs,finegridblocks)

  deallocate(has_children)
  deallocate(nodetype_copy)

! print out the final grid state
!  call amr_plotfile_chombo(time_step-back_step_ct)
!  call amr_checkpoint_wr(time_step-back_step_ct, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold)
!  write(61,*) 'Step, time', time_step-back_step_ct, time

!101 format(e19.13)

end subroutine amr_mg_time_dep_control

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_input_parameters(pid, noprocs)
  use multigrid_parameters
  use time_dep_parameters  
  use generic_parameters
  use tree                 
  use paramesh_dimensions  
  use paramesh_interfaces
  use checkpoint_parameters
  implicit none
  include 'mpif.h'
  integer,intent(in) :: pid, noprocs
  integer :: ierr
  !BUFFER holds the command line arguments
  CHARACTER *100 BUFFER  
  integer, intrinsic :: iargc
  logical :: there

  if (IArgC().gt.0) then            ! Do we have command line arguments?
     if (IArgC().eq.2) then         ! Are we expecting to specify which checkpoint to load
        call getarg(1,BUFFER)
        read(BUFFER,*) chk_restart
        if (chk_restart.eq.1) then  ! Have we done it right?
           call getarg(2,BUFFER)
           read(BUFFER,*) chk_chkpt
        else
           if (pid.eq.0.and.verbose.ge.0) then
              write(6,*) 'ERROR: Command line usage is:'
              write(6,*) './executable          for new run'
              write(6,*) './executable 1 <n>    to load checkpoint <n>'
              write(6,*) './executable 2        to load latest checkpoint, saved in CHK.out'
           endif
           stop
        endif
     else if (IArgC().eq.1) then
        if(pid.eq.0) then
           inquire(file="CHK.out",exist=there)
           if (there) then
!!              open (unit=99,file="CHK.out")
!!              read(99,*) chk_chkpt
!!              close(99)
           else
              if (verbose.ge.1) write(6,*) 'File CHK.out not found.  Starting new simulation.'
              chk_chkpt=0
           endif
        endif
        call MPI_Bcast(chk_chkpt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     else
        if (pid.eq.0.and.verbose.ge.0) then
           write(6,*) 'ERROR: Command line usage is:'
           write(6,*) './executable          for new run'
           write(6,*) './executable 1 <n>    to load checkpoint <n>'
           write(6,*) './executable 2        to load latest checkpoint, saved in CHK.out'
        endif
        stop
     endif

     if (chk_chkpt.gt.0) then
        chk_checkf = "hdf5"

        call app_read_checkpoint(pid, noprocs)
!        call amr_guardcell(pid,1,1)

!        call amr_refine_derefine
!        if(pid.eq.0) write(6,*) 'Checkpoint load reordered again'
        time = chk_t
        dt = chk_dt
        dtold = chk_dtold
        if (verbose.ge.3) call graphviz_output(pid, noprocs, chk_chkpt)

        if (pid.eq.0.and.verbose.ge.2 ) then  
           write(6,'(A)') 'Global domain limits'
           write(6,'(F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X)') grid_xmin,grid_xmax,grid_ymin,grid_ymax,grid_zmin,grid_zmax
        endif

!        call amr_multigrid_block_types(pid, noprocs, finegridblocks)
!        call amr_multigrid_block_timer(pid,noprocs)
        chk_chkpt=chk_chkpt+1

!!--tim
!!        mode_1 = 2
!!        mode_1 = 1
     endif
  endif

  call mpi_amr_global_domain_limits() ! This routine sets grid_xmax, grid_xmin, ...
  call set_domain_limits()            ! This is a local routine to turn output from mpi_amr_global_domain_limits into boundary_box array

  ! Initialise everything
  refine(1:lnblocks) = .false.
  call amr_refine_derefine

!!  if(interpolating.eq.1)then
!!    call tim_interpolation(pid, noprocs)
!!    call tim_two_checkpoints(pid, noprocs)
!!  endif


  if (pid.eq.0.and.verbose.ge.2) then  
     write(6,'(A)') 'Global domain limits'
     write(6,'(F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X)') grid_xmin,grid_xmax,grid_ymin,grid_ymax,grid_zmin,grid_zmax
  endif

  ! For new jobs do initial adaptation
  if (chk_chkpt.eq.0) then
     if(pid.eq.0.and.verbose.ge.1) write(6,*) ' starting new job'
     chk_chkpt=1
     call amr_initial_meshing(pid, noprocs)
  endif

end subroutine get_input_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine timestep_adjust(defect, local_r, incset, failer, i_step, pid)
  use generic_parameters
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_interfaces!, only :amr_guardcell, & 
       !&                         amr_prolong, & 
       !&                         amr_restrict
  implicit none
  integer, intent(in) :: i_step, pid
  integer, intent(inout) :: incset, failer
  double precision, intent(inout) :: defect, local_r

! Not done Local error test
  if (local_r.lt.0.0) then
     ! James's method:
     !   In all failure cases dt->.75*dt and timestep is restarted (failure is lack of convergence or NANs)
     !   If converged in small (<4) cycles dt->10/9*dt then next cycle
     !   If converged in large (>8) cycles dt->9/10*dt then next cycle (note, not reset since it converged)
     if(defect.eq.0.0)then
        ! If it's = 0, then in all likelyhood either everything's melted, or
        ! the grid's full of NANs
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[37;1;48m"
        if(pid.eq.0.and.verbose.ge.1)print *,"Zero defect - let's ignore for now - treat as tiny"
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        if (incset.eq.0) then
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[38;1;47m"
           if(pid.eq.0.and.verbose.ge.1)print *,"Low v-cycle count, increasing dt"
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
           dtold=dt
           dt=10.0/9.0*dt
           if(dt.gt.max_dt)then
              dt=max_dt
           end if
!           incset=1
           incset=1
        else
           incset = incset+1
        endif

!        if(pid.eq.0)print *,"Zero defect, asumming failure"
!        call amr_mg_time_step_reset()
!        back_step_ct=back_step_ct+1
!        time=time-dt
!        dt=dt*.75
!        incset = 0
!        failer=1
     else if(defect.gt.defect_too_big)then
        ! Pretty obvious really
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[39;1;46m"
        if(pid.eq.0.and.verbose.ge.1)print *,"Large defect, asumming failure"
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.75
        incset = 0
        failer=1
     else if(i_step.le.low_vcycle_count.and.(defect.lt.defect_tol_min.or.shampine_test.eq.1))then
        if (MOD(incset,stepsize_increase_frequency).eq.0) then
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[40;1;45m"
           if(pid.eq.0.and.verbose.ge.1)print *,"Low v-cycle count, increasing dt"
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
           dtold=dt
           dt=10./9.*dt
           if(dt.gt.max_dt)then
              dt=max_dt
           end if
           incset = 1
        else
           incset = incset+1
        endif
     else if((i_step.ge.max_v_cycle_count).or.(defect.gt.defect_tol_min.and.shampine_test.eq.0))then
!     else if((i_step.ge.10).or.(defect.gt.defect_tol))then
!         if(pid.eq.0)print *,"Timestep decreased removed"
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[41;1;44m"
        if(pid.eq.0.and.verbose.ge.1)print *,"Excessive v-cycle count/lack of convergence, decreasing dt and restarting time step"
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.75
        incset = 1
        failer=1
     else if(i_step.gt.high_vcycle_count)then
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[42;1;43m"
        if(pid.eq.0.and.verbose.ge.1)print *,"High v-cycle count, decreasing dt"   
if (pid.eq.0.and.verbose.ge.1)write(6,'(A,A)') 27,"[0m"
        dtold=dt
        dt=0.9*dt
        incset = 1
     end if
  else     ! Local error estimation techniques
     if(defect.gt.defect_too_big.or.(i_step.ge.max_v_cycle_count).or.(defect.gt.defect_too_big.and.shampine_test.eq.0))then
        ! Pretty obvious really
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: Large defect, asumming failure", defect
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.5
        incset = 0
        failer=1
     else if (local_r.lt.lee_tol_failing) then
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: Tiny local_r, asumming failure", local_r
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.5
        incset = 0
        failer=1
     else if (local_r.lt.lee_tol_halving) then
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: Small local_r, halving stepsize", local_r
        dtold=dt
        dt=dt*0.5
        incset = 0
     else if (local_r.lt.lee_tol_decrease) then
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: local_r requests decreased timestep", local_r
        dtold=dt
        dt=dt*local_r
        incset = 0
     else if (local_r.gt.lee_tol_increase.and.dt.lt.max_dt) then
        if(pid.eq.0.and.verbose.ge.1)print *,"LEE: local_r requests increased timestep", local_r
        dtold=dt
        dt=dt*2.0
        if (dt.gt.max_dt) dt = max_dt
        incset = 1
     endif
  endif

end subroutine timestep_adjust


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nodebound
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none

  integer i, j, k, lb, v


!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,v,k,j,i)
!!!$OMP DO SCHEDULE(DYNAMIC,8) 
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        do k=kl_bnd+1,ku_bnd-1
           do j=jl_bnd+1,ju_bnd-1
              do i=il_bnd+1,iu_bnd-1
                 do v = 1, nvar
                    if (ABS(unk(v,i,j,k,lb)).lt.1.0e-15) then
!                       print *,v,i,j,k,lb,unk(v,i,j,k,lb), ABS(unk(v,i,j,k,lb)
                       unk(v,i,j,k,lb)=0.0
                    endif
                 end do
              end do
           end do
        end do
     end if
   end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

end subroutine nodebound

