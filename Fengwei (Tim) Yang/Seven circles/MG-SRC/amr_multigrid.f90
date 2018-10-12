!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_runtime_rhs
!!!!  * Sets up right hand side variable per unknown using previous timesteps
!!!!  * Also predicts initial guess at solution for new timestep
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_time_step_reset
!!!!  * Rolls back to start of previous timestep if last attempt failed
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_multigrid_v_cycle
!!!!  * Recursive 2-grid V-cycle call
!!!!  * Calls user supplied function
!!!!             -> app_mg_smooth
!!!! 
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_get_defect
!!!!  * Calculates defect contribution and the maximum defect and RMS residual
!!!!  * Calls user supplied function
!!!!             -> app_mg_get_defect_var
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_multigrid_child_set
!!!!  * Finds out which blocks have children and sets has_children() array
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  reset_nodetype
!!!!  * This resets nodetypes so that MLAT can work alongside guardcell call
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_mg_get_max_refinement
!!!!  * Finds maximum refinement level on a processor
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_runtime_rhs()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  integer :: lb,i,j,k,v,u
  double precision :: rfactor, rf1, rf2, rf3
 
  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
!  rf1=dt
!  rf2=1.0
!  rf3=0.0

  ! Note this could be done more simply, however this is left as is for additional functionality
  if(lnblocks.gt.0) then
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k,u,v)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
    do lb=1,lnblocks
!if(nodetype(lb).eq.1) then
!if(nvar.gt.0) then
      do k=kl_bnd,ku_bnd
        do j=jl_bnd,ju_bnd
          do i=il_bnd,iu_bnd
            do v=1,total_vars
! set the start point for the real variable and its workspace data
              u=1+(v-1)*nunkvbles
              if(mode_1.eq.1)then
                if(v.eq.1)then
                  unk(u+1,i,j,k,lb)=unk(u,i,j,k,lb)
                  ! back up previous timesteps
                  unk(u+4,i,j,k,lb)=unk(u+3,i,j,k,lb)
                  unk(u+3,i,j,k,lb)=unk(u,i,j,k,lb)
                endif
                if(v.eq.2)then
                  unk(u+1,i,j,k,lb)=0.0
                  ! back up previous timesteps
                  unk(u+4,i,j,k,lb)=unk(u+3,i,j,k,lb)
                  unk(u+3,i,j,k,lb)=unk(u,i,j,k,lb)
                endif

! set initial guess as first order prediction
!!                             unk(u,i,j,k,lb)=unk(u+3,i,j,k,lb) + &
!!                                rfactor*(unk(u+3,i,j,k,lb)-unk(u+4,i,j,k,lb))
              else if(mode_1.eq.2) then
! back up previous timesteps
                if(v.eq.1)then
                  unk(u+4,i,j,k,lb)=unk(u+3,i,j,k,lb)
                  unk(u+3,i,j,k,lb)=unk(u,i,j,k,lb)
! set up RHS based on previous solutions
                  unk(u+1,i,j,k,lb)= rf2*unk(u+3,i,j,k,lb) - rf3*unk(u+4,i,j,k,lb)
                endif
                if(v.eq.2.or.v.eq.3)then
                  unk(u+1,i,j,k,lb) = 0.0
                  ! back up previous timesteps
                  unk(u+4,i,j,k,lb)=unk(u+3,i,j,k,lb)
                  unk(u+3,i,j,k,lb)=unk(u,i,j,k,lb)
                endif


              end if
            end do
          end do
        end do
      end do
!end if
!end if
    end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
  end if

!!--fwy  if (shampine_test.eq.1) call shampine_store()
end subroutine amr_mg_runtime_rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_time_step_reset()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use tree
!  use solution_parameters
  implicit none
  integer :: lb,i,j,k,u,v
  ! Note this could be done simpler, however this is left as is for additional functionality
  if(lnblocks.gt.0) then
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k,u,v)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
     do lb=1,lnblocks
!if(nvar.gt.0) then
       do k=kl_bnd,ku_bnd
         do j=jl_bnd,ju_bnd
           do i=il_bnd,iu_bnd
             do v=1,total_vars
               ! set the start point for the real variable
               u=1+(v-1)*nunkvbles
               if(mode_1.eq.1)then
!!               unk(u,i,j,k,lb)=unk(u+1,i,j,k,lb)
                 unk(u,i,j,k,lb) = unk(u+3,i,j,k,lb)
               else if(mode_1.eq.2)then
                 unk(u,i,j,k,lb)=unk(u+3,i,j,k,lb)
                 unk(u+3,i,j,k,lb)=unk(u+4,i,j,k,lb)
               end if
             end do
           end do
         end do
       end do
!end if
     end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
  end if
end subroutine amr_mg_time_step_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine amr_multigrid_v_cycle(pid,noprocs,level,global_max_defect, rmsres, npts)
!
  use paramesh_interfaces
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use tree
  use workspace
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  implicit none  
  include 'mpif.h'
  integer, intent(in) :: pid,noprocs, level
  integer, intent(inout) :: npts
  double precision, intent(inout) :: global_max_defect, rmsres(total_vars)
  integer :: smooth_count,lb,nlayers,iopt,ierr
  double precision :: local_max_defect, global_total
  integer :: i,j,k,new_level
  integer :: nguard0, core
  ! Plan:
  ! Need to cycle over all the leaf (node) blocks smoothing
  ! Next job is to find the defect
  ! Now do the FAS bit...
  ! Restrict if leaf nodes are not the coarsest level
  ! Possibly have to swap lots of data about here, 
  !    probably put this in it's own sub for my sanity
  ! Call solver again (recursive)
  ! More swapping? If so own sub for sanity
  ! Refine back
  ! Reverse FAS bit
  ! Smooth again
  ! Find defect/residual
  ! Setup a couple of parameters
   
  iopt = 1
  nlayers = nguard

  ! Calculate initial residual if wanted for output purposes only
  if (verbose.ge.2.and.level.eq.lrefine_max) then
    do current_var=1,total_vars
      current_unk=1+(current_var-1)*nunkvbles
      current_work=1
      local_max_defect=0.0
      call amr_mg_get_defect(pid, level, local_max_defect, rmsres, npts)
    end do
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Multigrid bit
  ! Note to self, this isn't as simple as suggested on the website
  ! this restriction (and prolongation) only operates on 'work'
  if(level.gt.mg_min_lvl)then

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Pre-smooth
 
    do smooth_count=1,smooth_count_max
   
      call app_mg_smooth(pid, level)
      call MPI_Allreduce(mg_update_total, global_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
      if (global_total.le.0.0) then
!!new        if (pid.eq.0.and.verbose.ge.4) print '(A,I,A)', 'Pre Update get-out after ',smooth_count, ' iterations'
        exit
      endif
!        call amr_guardcell(pid,iopt,level)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! CEG - i.e. coarsening

    ! CEG In order to fix the MG coarsening and defect calculations we need to coarsen 
    !     the solutions first of ALL variables, then calculate the residuals on both 
    !     fine and coarse grids without the additional RHS contribution, before adding 
    !     the defect into unk(2) etc
     
    ! * Coarsen all solutions
    gcell_on_cc(:)=.false.
    do current_var=1,total_vars
      current_unk=1+(current_var-1)*nunkvbles
      gcell_on_cc(current_unk)=.true.
      gcell_on_cc(current_unk+1)=.true.
      gcell_on_cc(current_unk+3)=.true.
      gcell_on_cc(current_unk+4)=.true.
      gcell_on_cc(current_unk+5)=.true.
      gcell_on_cc(current_unk+6)=.true.
    end do
    call pf_restrict(pid, 1, 0, .false., level)

    call reset_nodetype(noprocs, pid, level-1)
    call pf_guardcell(pid,iopt,level-1)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Start looping code block 1
    ! Note to self, within the "looping code" only make guardcell calls on work array
    ! put the full guardcell call at the end
    do current_var=1,total_vars
     
      ! set the start point for the real variable and it's workspace data
      current_unk=1+(current_var-1)*nunkvbles
      current_work=1

      ! Reset to current ref level, this only really needs to be done for current var>1
      ! Defect calculation
      if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)print *, 'Residual to coarsening'
      call zero_level_work_array(level)
      local_max_defect=0.0
      call amr_mg_get_defect(pid, level, local_max_defect, rmsres, npts)

      ! Need guardcells of defect setting up, defect is in work(2) so call guardcell with iopt=3
      call reset_nodetype(noprocs,pid,level)
      call pf_guardcell(pid,current_work+2,level)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
50    format(1x,i3,12(2x,f6.3))
!!      if(current_var.eq.3)then
!!        do j = 1, nyb+2*nguard
!!          write(*,50) j, (unk(12,i,j,1,1),i=1,nxb+2*nguard)
!!        enddo !do j=1, nyb+2*nuard
!!      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! CEG Restrict work(2)
      call pf_mg_restrict(noprocs,pid,level)
      call reset_nodetype(noprocs,pid,level-1)

      ! CEG Now calculate RHS defect
      if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)print *, 'Residual of coarsened solution'
      local_max_defect=0.0
      call amr_mg_get_defect(pid, level-1,local_max_defect, rmsres, npts)

      ! CEG  Finally add FAS RHS to solution

      !!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k)
      !!!$OMP DO SCHEDULE(DYNAMIC,1)
      do lb=block_starts(level-1), block_starts(level)-1
        do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
          do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
            do i=il_bnd+nguard,iu_bnd-nguard
              if(has_children(lb)) then
                ! RHS <-- previous timesteps - defect of coarsened solution
                unk(current_unk+1,i,j,k,lb)=unk(current_unk+1,i,j,k,lb)-work(i,j,k,lb,2)
                unk(current_unk+2,i,j,k,lb)=unk(current_unk,i,j,k,lb)   ! Store coarsened value as no longer done via work
              else
                unk(current_unk+2,i,j,k,lb)=unk(current_unk,i,j,k,lb)
              end if
            end do
          end do
        end do
      end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
        gcell_on_cc(current_unk:current_unk+4)=.false.  ! Preparatory for setting at end of loop for each variable
        gcell_on_cc(current_unk)=.true.
        gcell_on_cc(current_unk+1)=.true.
        gcell_on_cc(current_unk+3)=.true.
        gcell_on_cc(current_unk+4)=.true.
        gcell_on_cc(current_unk+5)=.true.
        gcell_on_cc(current_unk+6)=.true.
     end do
     call pf_guardcell(pid,iopt,level-1)
     gcell_on_cc(:)=.true.

    ! End looping code block 1
            
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Ok Think everything is now assembled let's try that recursive call
    new_level=level-1
    ! CEG - i.e. coarsened
    call reset_nodetype(noprocs, pid, level-1)
     
    call amr_multigrid_v_cycle(pid, noprocs, new_level, global_max_defect, rmsres, npts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ! Start looping code block 2
    ! Store coarse grid correction work(...1)
    ! Now we need to subtract the restricted version of v from the current (unk(3))

    do current_var=1,total_vars,1
      call reset_nodetype(noprocs,pid,level-1)
      ! set the start point for the real variable and its workspace data
      current_unk=1+(current_var-1)*nunkvbles
      current_work=1

      !!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k)
      !!!$OMP DO SCHEDULE(DYNAMIC,1)
      do lb=block_starts(level-1), block_starts(level)-1
        ! Put result straight into work ready for prolongation
        do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
          do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
            do i=il_bnd+nguard,iu_bnd-nguard
              if(has_children(lb))then
              ! work(...,1) is change in solution over the solve
                work(i,j,k,lb,current_work) = unk(current_unk,i,j,k,lb) - unk(current_unk+2,i,j,k,lb)
              else
                work(i,j,k,lb,current_work) = 0.0
              endif
            end do
          end do
        end do
      end do

      !!!$OMP END DO NOWAIT
      !!!$OMP END PARALLEL 
      call pf_guardcell(pid,1+current_work,level)
      ! Reset to current ref level, this only really needs to be done for current var>1
      call amr_mg_prolong(noprocs,pid,level)  ! Prolongs work array

      call reset_nodetype(noprocs, pid, level)
      call pf_guardcell(pid,1+current_work,level+1)

      call reset_nodetype(noprocs,pid,level)
      !!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k)
      !!!$OMP DO SCHEDULE(DYNAMIC,1)

      do lb=block_starts(level), block_starts(level+1)-1
        do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
          do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
            do i=il_bnd+nguard,iu_bnd-nguard
              ! Fine grid solution <-- pre-coarsened solution + coarse grid correction
              unk(current_unk,i,j,k,lb) = unk(current_unk,i,j,k,lb) + work(i,j,k,lb,current_work)
            end do
          end do
        end do
      end do
      !!!$OMP END DO NOWAIT
      !!!$OMP END PARALLEL 
      gcell_on_cc(:)=.false.
      gcell_on_cc(current_unk)=.true.
      call pf_guardcell(pid,iopt,level)
      gcell_on_cc(:)=.true.
    end do
    ! End looping code block 2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Post Smooth
    do smooth_count=1,smooth_count_max
      call app_mg_smooth(pid,level)
      call MPI_Allreduce(mg_update_total, global_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         MPI_COMM_WORLD, ierr)
      if (global_total.le.0.0) then
!!new        if (pid.eq.0.and.verbose.ge.4) print '(A,I,A)', 'Post Update get out after ',smooth_count, ' iterations'
        exit
      endif
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else
    ! CEG Coarse grid solve
    global_max_defect = 1.0
    do smooth_count=1,solve_count_max
      if(global_max_defect.lt.defect_tol_min*defect_tol_min)then
        exit
      end if
           
      call app_mg_smooth(pid,level)
      call MPI_Allreduce(mg_update_total, global_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         MPI_COMM_WORLD, ierr)
      if (global_total.le.0.0) then
!!new        if (pid.eq.0.and.verbose.ge.4) print '(A,I,A)', 'CG Update get out after ',smooth_count, ' iterations'
        exit
      endif
    end do

    local_max_defect=0.0
    do current_var=1,total_vars
      ! set the start point for the real variable and its workspace data
      current_unk=1+(current_var-1)*nunkvbles
      current_work=1
      local_max_defect=0.0
      if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)print *, 'Coarse grid solve defect', level
      call amr_mg_get_defect(pid, level,local_max_defect, rmsres, npts)
    end do
    call MPI_AllReduce(local_max_defect,global_max_defect,&
                      1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)
  end if

  ! Defect calculation
  local_max_defect=0.0
  do current_var=1,total_vars
    current_unk=1+(current_var-1)*nunkvbles
    current_work=1
    local_max_defect=0.0
    if (pid.eq.0.and.verbose.ge.4.and.current_var.eq.1)print *, 'Post smooths defect'
    call amr_mg_get_defect(pid, level,local_max_defect, rmsres, npts)
  end do

  ! NOTE: max of defects only needed at the end of v-cycle
  call MPI_AllReduce(local_max_defect,global_max_defect,&
                    1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)


!!!  call app_output_occasional(pid, noprocs, 6)
end subroutine amr_multigrid_v_cycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_get_defect(pid, level,local_max_defect, rms, npts)
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, level
  double precision, intent(inout) :: local_max_defect
  integer :: lb, npts, Gnpts, ierr, i
  double precision, intent(inout) ::  rms(total_vars)
  double precision :: global_max_defect, outres, Grms(total_vars), locrms

  rms(current_var) = 0.0
  npts = 0
  locrms=0.0
  global_max_defect = 0.0
!$OMP PARALLEL SHARED(level, local_max_defect, npts, locrms) PRIVATE(lb)
! REDUCTION(MAX:locrms)
!$OMP DO SCHEDULE(DYNAMIC,1) REDUCTION(+:locrms) REDUCTION(+:npts) REDUCTION(MAX:local_max_defect)
  do lb=block_starts(level), block_starts(level+1)-1
     call app_mg_get_defect_var(lb,local_max_defect, locrms, npts)
  end do
!$OMP END DO NOWAIT
  rms(current_var) = locrms
!$OMP END PARALLEL 
!  print *, locrms, rms(current_var), npts
  ! Output if verbose >= 3
  if ((verbose.eq.3.and.current_var.eq.total_vars).or.verbose.ge.4) then

     call MPI_Reduce(local_max_defect,global_max_defect,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_Reduce(rms, Grms, total_vars, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     call MPI_Reduce(npts, Gnpts, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     outres = 0.0
     if (pid.eq.0) then
        do i = 1, total_vars
           Grms(i) = Grms(i)/real(Gnpts)
           outres = outres + Grms(i)
           Grms(i) = sqrt(Grms(i))
        end do
     endif

     if (pid.eq.0) print '(A,I2,X,A,X,I2,X,A,E10.3,4X,A,E10.3)', 'Grid ', &
                level, 'Vble', current_var, 'Max defect: ', global_max_defect, 'RMS:', Grms(current_var)
!!fwy     write(*,*) Grms(current_var)
  end if

  return
end subroutine amr_mg_get_defect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zero_level_work_array(level)
  use multigrid_parameters
  use generic_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  include 'mpif.h'
  integer, intent(in) :: level
  integer :: lb

!!!$OMP PARALLEL SHARED(level) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
  do lb=block_starts(level), block_starts(level+1)-1
     work(:,:,:,lb,:) = 0.0
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

  return

end subroutine zero_level_work_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_multigrid_child_set()
  ! Sweeps through all (local) blocks and checks if they're parents or not
  ! This needs to be called once at start and once after all adaptivity events
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer lb
  ! Only operating on blocks not on finest possible as these should never have chldren set true
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
  do lb=1,lnblocks
!  do lb=1,block_starts(lrefine_max)-1
     if(child(1,1,lb).eq.-1)then
        has_children(lb)=.false.
     else if(child(1,1,lb).gt.0)then
        has_children(lb)=.true.
     else
        print *,"Odd value in child array for block ", lb, "Value: ", child(1,lb,1)
     end if
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 
end subroutine amr_multigrid_child_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reset_nodetype(noprocs,pid,level)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  integer, intent(in) :: noprocs,pid,level
  integer :: lb
  ! This resets nodetypes so that mlat can work alongside guardcell call
  ! Basically it's a fudge and only really sorts which blocks are nodetype 1
  ! It does have one obvious bonus, that is nodetype becomes redundant outside of GC
  ! calls and this can reset it if needed
  ! Conditions for nodetypes
  ! if a block is more refined than the current "level" do nothing
  ! if a block is same as current level set to 1
  ! if a block is less refined than current level AND has no children, set to 1
  ! if a block is less refined than current level AND has children, leave alone
  ! <snip>
  ! James: Guess what, I have a NEW way of doing this now
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
  do lb=1,lnblocks
     nodetype(lb)=nodetype_copy(lb)
     if(lrefine(lb).gt.level)then
        nodetype(lb) = -1
     else if(lrefine(lb).eq.level)then
        nodetype(lb) = 1
     else if((lrefine(lb).eq.level-1).and.(nodetype(lb).ne.1))then
        nodetype(lb) = 2
     end if
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

  call amr_get_new_nodetypes (noprocs, pid, level)

end subroutine reset_nodetype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_mg_get_max_refinement(pid,noprocs)
  ! This subroutine sweeps through all blocks and finds the max ref level on the current processor
  ! Then it passes the max to all procs
  ! Why is this needed? Well multigrid needs to know what level to start at
  ! Should be called after all refinement events
  use tree
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid,noprocs
  integer :: lb,local_max
  local_max=0
!!!$OMP PARALLEL SHARED(lnblocks,local_max) PRIVATE(lb)
!!!$OMP DO SCHEDULE(DYNAMIC,1) REDUCTION (MAX:local_max)
  do lb=1,lnblocks
     if(lrefine(lb).gt.local_max)then
        local_max=lrefine(lb)
     end if
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

  call MPI_AllReduce(local_max,mg_max_level,&
       1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,lb)
end subroutine amr_mg_get_max_refinement

