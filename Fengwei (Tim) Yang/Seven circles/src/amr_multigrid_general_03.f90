!#include "paramesh_preprocessor.fh"

subroutine amr_mg_runtime_rhs()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  use solution_parameters
  implicit none
  integer :: lb,i,j,k
  double precision :: rfactor, rf1, rf2, rf3

  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
!  rf1=dt
!  rf2=1.0
!  rf3=0.0


  ! Note this could be done simpler, however this is left as is for additional functionality
  if(lnblocks.gt.0) then
     do lb=1,lnblocks
        !if(nodetype(lb).eq.1) then
!           if(nvar.gt.0) then
              do k=kl_bnd,ku_bnd
                 do j=jl_bnd,ju_bnd
                    do i=il_bnd,iu_bnd
                       do current_var=1,total_vars
                          ! set the start point for the real variable and its workspace data
                          current_unk=1+(current_var-1)*nunkvbles
                          if(mode_1.eq.1)then
                             unk(current_unk+1,i,j,k,lb)=unk(current_unk,i,j,k,lb)
! back up previous timesteps
                             unk(current_unk+4,i,j,k,lb)=unk(current_unk+3,i,j,k,lb)
                             unk(current_unk+3,i,j,k,lb)=unk(current_unk,i,j,k,lb)
! set initial guess as first order prediction
                             unk(current_unk,i,j,k,lb)=unk(current_unk+3,i,j,k,lb) + &
                                rfactor*(unk(current_unk+3,i,j,k,lb)-unk(current_unk+4,i,j,k,lb))
                          else if(mode_1.eq.2) then
! back up previous timesteps
                             unk(current_unk+4,i,j,k,lb)=unk(current_unk+3,i,j,k,lb)
                             unk(current_unk+3,i,j,k,lb)=unk(current_unk,i,j,k,lb)
! set up RHS based on previous solutions
!                             unk(current_unk+1,i,j,k,lb)=unk(current_unk,i,j,k,lb)
                             unk(current_unk+1,i,j,k,lb)= rf2*unk(current_unk+3,i,j,k,lb) - rf3*unk(current_unk+4,i,j,k,lb)
! set initial guess as first order prediction
                             unk(current_unk,i,j,k,lb)=unk(current_unk+3,i,j,k,lb) + &
                                rfactor*(unk(current_unk+3,i,j,k,lb)-unk(current_unk+4,i,j,k,lb))
!                          else
!                             unk(current_unk+1,i,j,k,lb)=unk(current_unk,i,j,k,lb)/dt
                          end if
                       end do
                    end do
                 end do
              end do
           !end if
        !end if
     end do
  end if

  if (shampine_test.eq.1) call shampine_store()


end subroutine amr_mg_runtime_rhs

subroutine amr_mg_time_step_reset()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use tree
  use solution_parameters
  implicit none
  integer :: lb,i,j,k
  ! Note this could be done simpler, however this is left as is for additional functionality
  if(lnblocks.gt.0) then
     do lb=1,lnblocks
        !if(nvar.gt.0) then
           do k=kl_bnd,ku_bnd
              do j=jl_bnd,ju_bnd
                 do i=il_bnd,iu_bnd
                    do current_var=1,total_vars
                       ! set the start point for the real variable
                       current_unk=1+(current_var-1)*nunkvbles
                       if(mode_1.eq.1)then
                          unk(current_unk,i,j,k,lb)=unk(current_unk+1,i,j,k,lb)
                       else if(mode_1.eq.2)then
                          unk(current_unk,i,j,k,lb)=unk(current_unk+3,i,j,k,lb)
                          unk(current_unk+3,i,j,k,lb)=unk(current_unk+4,i,j,k,lb)
                       end if
                    end do
                 end do
              end do
           end do
        !end if
     end do
  end if

end subroutine amr_mg_time_step_reset


recursive subroutine amr_multigrid_v_cycle(mype,nprocs,level,global_max_defect)
!
!    
!
  use paramesh_interfaces
  use multigrid_parameters
  use paramesh_dimensions
  use solution_parameters
  use time_dep_parameters
  use physicaldata
  use tree
  use workspace
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  implicit none  
  include 'mpif.h'
  integer, intent(in) :: mype,nprocs
  integer, intent(in) :: level
  real, intent(inout) :: global_max_defect
  integer :: smooth_count,lb,nlayers,iopt,ierr
  real :: local_max_defect
  integer :: i,j,k,new_level
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
  ! Pre-smooth
  do smooth_count=1,smooth_count_max
     call amr_mg_smooth(mype,level)
!     call amr_guardcell(mype,iopt,level)
  end do
  
  ! Multigrid bit
  ! Note to self, this isn't as simple as suggested on the website
  ! this restriction (and prolongation) only operates on 'work'
  if(level.gt.mg_min_lvl)then
  ! CEG - i.e. coarsening
     
     ! Start looping code block 1
     ! Note to self, within the "looping code" only make guardcell calls on work array
     ! put the full guardcell call at the end
     do current_var=1,total_vars
        ! set the start point for the real variable and it's workspace data
        current_unk=1+(current_var-1)*nunkvbles
        current_work=1
        ! Reset to current ref level, this only really needs to be done for current var>1
        ! Defect calculation
        local_max_defect=0.0
        call amr_mg_get_defect(level,local_max_defect)
        ! Need guardcells of defect setting up, defect is in work(2) so call guardcell with iopt=3
        call reset_nodetype(nprocs,mype,level)
        call pf_guardcell(mype,current_work+2,level)
!        call amr_guardcell(mype,current_work+2,nlayers)
        
        ! at this stage guardcells are correct in unk and work(2)
!        do lb=1,lnblocks
        do lb=block_starts(level), block_starts(level+1)-1
!           if(lrefine(lb).eq.level)then
              do k=kl_bnd,ku_bnd
                 do j=jl_bnd,ju_bnd
                    do i=il_bnd,iu_bnd
                       work(i,j,k,lb,1)=unk(current_unk,i,j,k,lb)
                       work(i,j,k,lb,2)=unk(current_unk+1,i,j,k,lb)
                       if (nunkvbles.gt.3) then
                          work(i,j,k,lb,3)=unk(current_unk+3,i,j,k,lb)
                          work(i,j,k,lb,4)=unk(current_unk+4,i,j,k,lb)
                       endif
                    end do
                 end do
              end do
!           end if
        end do
        ! guardcells correct in both work elements
        call amr_mg_restrict(nprocs,mype,level)
!        call pf_mg_restrict(nprocs,mype,level)
        call reset_nodetype(nprocs,mype,level-1)
        ! ok work on level-1 now contains v (work(1)) and d (work(2))
        ! Now copy these into unk(1) and unk(2)
!        do lb=1,lnblocks
        do lb=block_starts(level-1), block_starts(level)-1
           if(has_children(lb))then
              do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
                 do j=jl_bnd+nguard,ju_bnd-nguard
                    do i=il_bnd+nguard,iu_bnd-nguard
                       unk(current_unk,i,j,k,lb)=work(i,j,k,lb,1)
                       unk(current_unk+1,i,j,k,lb)=work(i,j,k,lb,2)
                       unk(current_unk+2,i,j,k,lb)=work(i,j,k,lb,1)
!                       if (nunkvbles.gt.3) then
!                          unk(current_unk+3,i,j,k,lb)=work(i,j,k,lb,3)
!                          unk(current_unk+4,i,j,k,lb)=work(i,j,k,lb,4)
!                       endif
                    end do
                 end do
              end do
           end if
        end do
        ! Selective guardcell filling
        gcell_on_cc(:)=.false.
        gcell_on_cc(current_unk)=.true.
!        call amr_guardcell(mype,iopt,nguard)
        call pf_guardcell(mype,iopt,level-1)
        ! Now reset to fill all vars guardcells
        gcell_on_cc(:)=.true.
     end do
     !call amr_guardcell(mype,iopt,nlayers)
     ! End looping code block 1
          
     ! Start looping code block 2
     do current_var=1,total_vars
        ! set the start point for the real variable and it's workspace data
        current_unk=1+(current_var-1)*nunkvbles
!        current_work=1
        ! Now we need to add N(v) to the rhs (unk(2))
        ! get defect will find f-N(v) and put it in work(2)
        call amr_mg_get_defect(level-1,local_max_defect)
!        do lb=1,lnblocks
        do lb=block_starts(level-1), block_starts(level)-1
           if(has_children(lb))then
              do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
                 do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                    do i=il_bnd+nguard,iu_bnd-nguard
                       unk(current_unk+1,i,j,k,lb)=2.0*unk(current_unk+1,i,j,k,lb)-work(i,j,k,lb,2)
                       if (nunkvbles.gt.3) then
                          unk(current_unk+3,i,j,k,lb)=2.0*unk(current_unk+3,i,j,k,lb)-work(i,j,k,lb,2)
                          unk(current_unk+4,i,j,k,lb)=2.0*unk(current_unk+4,i,j,k,lb)-work(i,j,k,lb,2)
                       endif
                    end do
                 end do
              end do
           end if
        end do
        gcell_on_cc(:)=.false.
        gcell_on_cc(current_unk+1)=.true.
        if (nunkvbles.gt.3) then
           gcell_on_cc(current_unk+3)=.true.
           gcell_on_cc(current_unk+4)=.true.
        endif
!        call amr_guardcell(mype,iopt,nlayers)
        call pf_guardcell(mype,iopt,level-1)
        gcell_on_cc(:)=.true.
     end do
     ! End looping code block 2
             
     ! Ok Think everything is now assembled (sort of) let's try that recursive call
     new_level=level-1
     ! CEG - i.e. coarsened
     call reset_nodetype(nprocs,mype,level-1)
     call amr_multigrid_v_cycle(mype,nprocs,new_level,global_max_defect)
     ! Ok we got this far, what now?
     
     ! Start looping code block 3
     ! Now we need to subtract the restricted version of v from the current (unk(3))
     do current_var=1,total_vars,1
        call reset_nodetype(nprocs,mype,level-1)
        ! set the start point for the real variable and it's workspace data
        current_unk=1+(current_var-1)*nunkvbles
        current_work=1
!        do lb=1,lnblocks
        do lb=block_starts(level-1), block_starts(level)-1
           if(has_children(lb))then
              ! Put result straight into work ready for prolongation
              do k=kl_bnd,ku_bnd
                 do j=jl_bnd,ju_bnd
                    do i=il_bnd,iu_bnd
                       work(i,j,k,lb,current_work)=unk(current_unk,i,j,k,lb)-unk(current_unk+2,i,j,k,lb)
                       !if(current_var+1.le.total_vars)then
                       !   work(i,j,k,lb,2)=unk(current_unk+3,i,j,k,lb)-unk(current_unk+5,i,j,k,lb)
                       !end if
                    end do
                 end do
              end do
           !else
           !  do k=kl_bnd,ku_bnd
           !     do j=jl_bnd,ju_bnd
           !        do i=il_bnd,iu_bnd
           !           work(i,j,k,lb,2)=0.0
           !        end do
           !     end do
           !  end do
          end if
        end do
        ! Reset to current ref level, this only really needs to be done for current var>1
        call amr_mg_prolong(nprocs,mype,level)
!        call pf_mg_prolong(nprocs,mype,level)
        call reset_nodetype(nprocs,mype,level)
!        do lb=1,lnblocks
        do lb=block_starts(level), block_starts(level+1)-1
!           if(lrefine(lb).eq.level)then
              !do k=kl_bnd,ku_bnd
              !   do j=jl_bnd,ju_bnd
              !      do i=il_bnd,iu_bnd
              do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
                 do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                    do i=il_bnd+nguard,iu_bnd-nguard
                       unk(current_unk,i,j,k,lb)=unk(current_unk,i,j,k,lb)+work(i,j,k,lb,current_work)
                       !if(current_var+1.le.total_vars)then! Test to see whether second var is present
                       !   unk(current_unk+3,i,j,k,lb)=unk(current_unk+3,i,j,k,lb)+work(i,j,k,lb,2)
                       !end if
                    end do
                 end do
              end do
!           end if
        end do
        gcell_on_cc(:)=.false.
        gcell_on_cc(current_unk)=.true.
!        call amr_guardcell(mype,iopt,nlayers)
        call pf_guardcell(mype,iopt,level)
        gcell_on_cc(:)=.true.
     end do
     !call amr_guardcell(mype,iopt,nlayers)
     ! End looping code block 3
  else
     ! CEG Coarse grid solve
     ! Should be a solve call here.....
     if(full_solve)then
        global_max_defect = 1.0
        do smooth_count=1,solve_count_max
           if(global_max_defect.lt.defect_tol**2)then
              exit
           end if
           call amr_mg_smooth(mype,level)
           !call amr_guardcell(mype,iopt,nlayers)
           local_max_defect=0.0
           do current_var=1,total_vars
              ! set the start point for the real variable and it's workspace data
              current_unk=1+(current_var-1)*nunkvbles
              current_work=1
              call amr_mg_get_defect(level,local_max_defect)
           end do
           call MPI_AllReduce(local_max_defect,global_max_defect,&
                1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)
        end do
     end if
  end if
  ! Post Smooth
  do smooth_count=1,smooth_count_max
     call amr_mg_smooth(mype,level)
     !call amr_guardcell(mype,iopt,nlayers)
  end do
  ! Defect calculation
  local_max_defect=0.0
  do current_var=1,total_vars
     current_unk=1+(current_var-1)*nunkvbles
     current_work=1
     call amr_mg_get_defect(level,local_max_defect)
  end do
  ! NOTE: max of defects only needed at the end of v-cycle
  call MPI_AllReduce(local_max_defect,global_max_defect,&
       1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)
end subroutine amr_multigrid_v_cycle

subroutine amr_mg_smooth(mype,level)
  use paramesh_dimensions
  use tree
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use paramesh_interfaces
  implicit none
  include 'mpif.h'
  integer, intent(in) :: mype,level
  integer lb, iopt
  iopt=1
  ! Cycle over all blocks
  do current_var=1,total_vars
     ! set the start point for the real variable and it's workspace data
     current_unk=1+(current_var-1)*nunkvbles
     current_work=1
!     do lb=1,lnblocks
      do lb=block_starts(level), block_starts(level+1)-1
         call amr_mg_smooth_var(lb)
     end do
     ! Update gcells for "current_unk" only
     gcell_on_cc(:)=.false.
     gcell_on_cc(current_unk)=.true.
!     call amr_guardcell(mype,iopt,nguard)
     call pf_guardcell(mype,iopt,level)
     gcell_on_cc(:)=.true.
  end do
end subroutine amr_mg_smooth

subroutine amr_mg_get_defect(level,local_max_defect)
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer, intent(in) :: level
  double precision, intent(inout) :: local_max_defect
  integer :: lb
!  do lb=1,lnblocks
  do lb=block_starts(level), block_starts(level+1)-1
!     if(lrefine(lb).eq.level)then
!        if(nvar.gt.0)then
     call amr_mg_get_defect_var(lb,local_max_defect)
!        end if
!        if(nvarcorn.gt.0)then
!           print *,"No Code for corner variables yet"
!        end if
!        if(nfacevar.gt.0) then
!           print *,"No Code for face variables yet"
!        end if
!        if(nvaredge.gt.0) then
!           print *,"No Code for edge variables yet"
!        end if
!     end if
  end do
end subroutine amr_mg_get_defect

subroutine amr_mg_smooth_var(lb)
  ! Smoother for node centered variables
  ! lb is the block number
  ! n is the variable number being smoothed
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer,intent(in) :: lb
  integer :: i,j,k
  real :: new_soln(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
  real :: dx,Fi,dFi_dvi
  real :: LapMult,tau_val
  double precision :: rfactor, rf1, rf2, rf3, dx2inv

  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)

!  rf1=dt
!  rf2=1.0
!  rf3=0.0

!  new_soln(:,:,:)=unk(current_unk,:,:,:,lb)
  dx = bsize(1,lb)/real(nxb)
  dx2inv = 1.0/(dx*dx)
  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
     do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
           if (ndim.eq.3) then
              call get_stencil(i,j,k,dx,lb,1,Fi,LapMult)
              call get_stencil_diff(i,j,k,dx,lb,1,dFi_dvi,LapMult)
           else
              call get_stencil_2d(i,j,dx,lb,1,Fi,dx2inv,LapMult,tau_val)
              call get_stencil_diff_2d(i,j,dx,lb,1,dFi_dvi,LapMult,tau_val)
           endif
! CEG : these 4 lines removed by James
           !call amr_multigrid_get_g(lb,i,j,k,temp)
!           Fi=Fi+temp
           !call amr_multigrid_get_c(lb,i,j,k,temp)
!           dFi_dvi=dFi_dvi+temp
           if(mode_1.eq.1)then
!           if(mode_1.ge.1)then
              Fi=unk(current_unk,i,j,k,lb)-dt*Fi-unk(current_unk+1,i,j,k,lb)
              dFi_dvi=1.0-dt*dFi_dvi
           else if(mode_1.eq.2)then
              Fi=unk(current_unk,i,j,k,lb)-rf1*Fi-unk(current_unk+1,i,j,k,lb)
!              Fi=unk(current_unk,i,j,k,lb)-rf1*Fi-(rf2*unk(current_unk+3,i,j,k,lb) - rf3*unk(current_unk+4,i,j,k,lb)) 
              dFi_dvi=1.0-rf1*dFi_dvi
!           else
!              Fi=unk(current_unk,i,j,k,lb)/dt-Fi-unk(current_unk+1,i,j,k,lb)
!              dFi_dvi=1.0/dt-dFi_dvi
           end if
           new_soln(i,j,k)=unk(current_unk,i,j,k,lb)-weight*Fi/dFi_dvi
!           if (gauss_seidel) unk(current_unk,i,j,k,lb) = new_soln(i,j,k)
        end do
     end do
  end do
!  if (.not.gauss_seidel) unk(current_unk,:,:,:,lb)=new_soln(:,:,:)
  unk(current_unk,:,:,:,lb)=new_soln(:,:,:)
end subroutine amr_mg_smooth_var

subroutine amr_mg_get_defect_var(lb,max_defect)
  ! Defect for node centered variables
  ! lb is the block number
  ! n is the variable number whose defect is being calculated
  ! defect will be put into work(i,j,k,lb,n+1)
  use multigrid_parameters
  use paramesh_dimensions
  use time_dep_parameters
  use solution_parameters
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: lb
  real, intent(inout) :: max_defect
  integer :: i,j,k
  real :: dx,Fi, dx2inv
  double precision :: rfactor, rf1, rf2, rf3
  real :: LapMult,tau_val

  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
!  rf1=dt
!  rf2=1.0
!  rf3=0.0

  dx = bsize(1,lb)/real(nxb)
  dx2inv = 1.0/(dx*dx)

  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
     do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
           if (ndim.eq.3) then
              call get_stencil(i,j,k,dx,lb,1,Fi,LapMult)
           else
              call get_stencil_2d(i,j,dx,lb,1,Fi,dx2inv,LapMult,tau_val)
           endif
! CEG : Call to amr_multigrid_get_g removed by James, hence no need for temp addition
           !call amr_multigrid_get_g(lb,i,j,k,temp)
!           Fi=Fi+temp
           if(mode_1.eq.1)then
!           if(mode_1.ge.1)then
              Fi=unk(current_unk,i,j,k,lb)-dt*Fi-unk(current_unk+1,i,j,k,lb)
           else if(mode_1.eq.2)then
              Fi=unk(current_unk,i,j,k,lb)-rf1*Fi-unk(current_unk+1,i,j,k,lb)
!              Fi=unk(current_unk,i,j,k,lb)-rf1*Fi-(rf2*unk(current_unk+3,i,j,k,lb) - rf3*unk(current_unk+4,i,j,k,lb))
!           else
!              Fi=unk(current_unk,i,j,k,lb)/dt-Fi-unk(current_unk+1,i,j,k,lb)
           end if
           work(i,j,k,lb,current_work+1)=-Fi
           if(Fi.lt.0.0)then
              Fi=-Fi
           end if
           if(Fi.gt.max_defect)then
              max_defect=Fi
           end if
        end do
     end do
  end do
end subroutine amr_mg_get_defect_var

subroutine amr_mg_get_rmsres(pid, noprocs)
  use paramesh_dimensions
  use physicaldata
  use tree
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer :: lb, i, j, k, npts, Gnpts, ierr
  integer, intent(inout) :: pid, noprocs
  double precision :: rms(3), Fi, dFi, res, Grms(3), dx, dx2inv, rf1, rf2, rf3, rfactor, LapMult, tau_val
  npts=0
  rms(:) = 0.0
  Gnpts=0
  Grms(:) = 0.0
  rfactor = dt/dtold
  rf1 = (rfactor+1.0)*dt/(2.0*rfactor+1.0)
!  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
!  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
  do current_var = 1, total_vars
     current_unk=1+(current_var-1)*nunkvbles
     do lb=1,lnblocks
        if(.not.has_children(lb)) then
           dx = bsize(1,lb)/real(nxb)
           dx2inv = 1.0/(dx*dx)
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if (ndim.eq.3) then
                       call get_stencil(i, j, k, dx, lb, 1, Fi,LapMult)
                    else
                       call get_stencil_2d(i, j, dx, lb, 1, Fi, dx2inv, LapMult, tau_val)
                    endif
                    Fi = unk(current_unk,i,j,k,lb)-rf1*Fi-unk(current_unk+1,i,j,k,lb)
                    res = Fi

                    rms(current_var) = rms(current_var) + res*res
                    npts = npts + 1
                 end do
              end do
           end do
        endif
     end do
  end do
  call MPI_Reduce(rms, Grms, total_vars, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(npts, Gnpts, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (pid.eq.0) then
     do i = 1, total_vars
        Grms(i) = Grms(i)/real(Gnpts)
        rms(i) = sqrt(Grms(i))
     end do
     print *, 'RMS = ', Grms(1:total_vars)
  endif
end subroutine amr_mg_get_rmsres


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
!  do lb=1,lnblocks
  do lb=1,block_starts(lrefine_max)-1
     if(child(1,1,lb).eq.-1)then
        has_children(lb)=.false.
     else if(child(1,1,lb).gt.0)then
        has_children(lb)=.true.
     else
        print *,"Odd value in child array for block ",lb,"Value: ",child(1,lb,1)
     end if
  end do
end subroutine amr_multigrid_child_set

subroutine reset_nodetype(nprocs,mype,level)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  integer, intent(in) :: nprocs,mype,level
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

!  if (mype==0) print *, level
  call amr_get_new_nodetypes (nprocs, mype, level)

end subroutine reset_nodetype


subroutine amr_mg_get_max_refinement(mype,nprocs)
  ! This subroutine sweeps through all blocks and finds the max ref level on the current processor
  ! Then it passes the max to all procs
  ! Why is this needed? Well multigrid needs to know what level to start at
  ! Should be called after all refinement events
  use tree
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: mype,nprocs
  integer :: lb,local_max
  local_max=0
  do lb=1,lnblocks
     if(lrefine(lb).gt.local_max)then
        local_max=lrefine(lb)
     end if
  end do
  call MPI_AllReduce(local_max,mg_max_level,&
       1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,lb)
end subroutine amr_mg_get_max_refinement

