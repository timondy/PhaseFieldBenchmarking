!#include "paramesh_preprocessor.fh"


subroutine amr_mg_control(mype,nprocs)
! Wrapper for controlling whole solver process
  integer, intent(in) :: mype, nprocs

  call amr_mg_time_dep_control(mype,nprocs)
end subroutine amr_mg_control

subroutine amr_mg_time_dep_control(mype,nprocs)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use solution_parameters
  use paramesh_interfaces!, only :amr_guardcell, & 
       !&                         amr_prolong, & 
       !&                         amr_restrict
  implicit none
  include 'mpif.h'
  integer, intent(in) :: mype, nprocs
  integer :: back_step_ct,time_step,i_step,i_step_total,ierr,level,iopt, lb
  double precision defect,end_time,old_defect_1,old_defect_2
  integer output_pe, failer, finegridblocks(2)
  double precision :: timeprev, timenew
  logical chk_gd
  character(len=80) :: chk_checkf
  double precision local_r, convrate

!CEG additions for continuation
  INTEGER::chk_restart=0, chk_chkpt, incset
  DOUBLE PRECISION:: chk_t, chk_dt, chk_dtold, wtime, wtimeprev

!DEFINE BUFFER HOLDS THE COMMAND LINE ARGUMENT
  CHARACTER *100 BUFFER  

  integer, external     :: iargc

  output_pe=0
  failer=0

  end_time=t0+simulation_time
  incset=0
  if(mype.eq.output_pe) write(6,*) 'end time is ',end_time

!GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT
  if (IArgC().gt.0) then
     if (IArgC().eq.2) then
        call getarg(1,BUFFER)
        read(BUFFER,*) chk_restart
        if (chk_restart.eq.1) then
           call getarg(2,BUFFER)
           read(BUFFER,*) chk_chkpt
        endif
     else 
        if(mype.eq.0) then
           open (unit=99,file="CHK.out")
           read(99,*) chk_chkpt
           close(99)
        endif
        call MPI_Bcast(chk_chkpt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     endif

     if (chk_chkpt.gt.0) then
           chk_checkf = "hdf5"

        if(mype.eq.output_pe) write(6,*) 'Checkpoint on : ', chk_chkpt
        call amr_checkpoint_re(chk_chkpt, user_attr_1=chk_dt, user_attr_2=chk_t, user_attr_3=chk_dtold)
        if(mype.eq.output_pe) write(6,*) 'Checkpoint read with ', chk_t, chk_dt, chk_dtold
!        call amr_refine_derefine
!        if(mype.eq.output_pe) write(6,*) 'Checkpoint load reordered again'
        t0 = chk_t
        dt = chk_dt
        dtold = chk_dtold
        call graphviz_output(mype, nprocs, chk_chkpt)

        if (mype.eq.0) then  
           write(6,'(A)') 'Global domain limits'
           write(6,'(F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X)') grid_xmin,grid_xmax,grid_ymin,grid_ymax,grid_zmin,grid_zmax
        endif

        if(mype.eq.output_pe) write(6,*) '---------------------------------------'
        call amr_multigrid_block_types(mype,nprocs,finegridblocks)
!        call amr_multigrid_block_timer(mype,nprocs)
        call phase_field_check(1)
        if(mype.eq.output_pe) write(6,*) '---------------------------------------'
        chk_chkpt=chk_chkpt+1

        mode_1 = 2
     else
        if(mype.eq.output_pe) write(6,*) ' starting new job'
        chk_chkpt=1
     endif
  else
     if(mype.eq.output_pe) write(6,*) IArgC(), ' command line arguments not found - not 2'
     chk_chkpt=1
  endif

  time=t0
  D_solute=a_2*lambda
  D_therm=Le*D_solute
!  if(total_vars.lt.3)then
!     D_therm=a_2*lambda
!  end if

!CEG Add README output
  if(mype.eq.output_pe.and.chk_chkpt.eq.1) then
     open (unit=61,file="README")
     write(61,*) 'Parameter list:'
     write(61,*) 'Mc-Inf      ',  mcinf
     write(61,*) 'kappa_E     ',  ke
     write(61,*) 'epsilon     ',  epsilon
     write(61,*) 'lambda      ',  lambda
     write(61,*) 'D_solute    ', D_solute
     write(61,*) 'D_therm     ', D_therm
     write(61,*) 'Lewis no.   ',  le
     write(61,*) 'Delta       ',  delta
     write(61,*) 'nucleus rad ', nuc_radius
     write(61,*) 'End time    ',  end_time
     write(61,*) 'Domain size ', grid_xmax
     write(61,*) 'Finest dx   ', min_dx
     write(61,*) 'Block size  ', nxb
     write(61,*) 'total_vars  ', total_vars
     write(61,*) 'Solute on?  ', solute
     write(61,*) ' '
     flush(61)
  endif

  if(mype.eq.output_pe)print *,"Entered multigrid control"
  iopt=1
  back_step_ct=0
  allocate(has_children(1:maxblocks),stat=time_step)!Using time_step as dud var
  allocate(nodetype_copy(1:maxblocks),stat=time_step)!Using time_step as dud var
!  call amr_guardcell(mype,iopt,nguard)
  has_children(:) = .false.
  
  call pf_guardcell(mype,iopt,mg_max_level)
  call amr_mg_init()
  call amr_multigrid_block_types(mype,nprocs,finegridblocks)
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!  call amr_guardcell(mype,iopt,nguard)
  call pf_guardcell(mype,iopt,mg_max_level)
  if(mype.eq.output_pe)print *,"Complete mg_init"
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  if(mype.eq.output_pe)print *,"Number of 'real' variables:",total_vars
  if(mype.eq.output_pe)print *,"Lambda: ",lambda,"Dtherm",D_therm,"Dsolute",D_solute
  if(mype.eq.0)open (unit=125,file="tip_rad.txt")
  if(mype.eq.0)open (unit=126,file="tip_loc.txt")
  if(mype.eq.0)open (unit=133,file="real_time.txt")

!  if(mype.eq.0)open (unit=301,file="Vcycnt.txt")
  
  i_step_total=0
  if(mype.eq.output_pe)print *,"Starting timestepping"
!  call ceg_block_balancing(mype, nprocs)
  call amr_multigrid_block_types(mype,nprocs,finegridblocks)
!  stop
  wtimeprev = MPI_Wtime()
  do time_step=chk_chkpt,max_time_iter+back_step_ct
     ! Little bit of code to ensure the last time step ends on simulation_time
!     if((time+dt).gt.end_time)then
!        dt=end_time-time
!     end if
     time=time+dt
     if(time.ge.end_time)exit

     if(mype.eq.output_pe)print *,"max_time_iter, back_step_ct:: ", max_time_iter,back_step_ct
     nodetype_copy(1:lnblocks)=nodetype(1:lnblocks)
     call amr_multigrid_child_set()
     if(mype.eq.output_pe)print *,""!Blank line between timesteps
     if(mype.eq.output_pe)print *,"Starting time step no.",time_step-back_step_ct,"times:",time-dt,"to",time
     if(mype.eq.output_pe)print *,"Current dt:",dt
     call amr_mg_get_max_refinement(mype,nprocs)
     call reset_nodetype(nprocs,mype,mg_max_level)
     if(failer.eq.0.and.chk_restart.ne.1)then
        if(mype.eq.output_pe)print *,mype,"Testing refinement"
        call amr_test_refinement(mype,lrefine_min,lrefine_max)
        if(mype.eq.output_pe)print *,"Refining"
        if(mype.eq.0) call system('date +"%k:%M:%S:%N"')
        if(mype.eq.0) wtimeprev = MPI_Wtime()
        if(mype.eq.0) then
           call cpu_time(timenew)
           timeprev=timenew
        endif
        call amr_refine_derefine
        if(mype.eq.0) then
           wtime = MPI_Wtime()
           call cpu_time(timenew)
           print*,'Sys/Real adapt time was ', timenew-timeprev, wtime-wtimeprev
           timeprev=timenew
           wtimeprev=wtime
        endif
        if(mype.eq.0)call system('date +"%k:%M:%S:%N"')

        call amr_multigrid_block_types(mype,nprocs,finegridblocks)
        call amr_prolong(mype,iopt,nguard)

        ! CEG try domain growing
        if (grow.eq.1) call amr_grow_domain(mype, nprocs, time_step-back_step_ct)

        if(mype.eq.output_pe)print *,"Prolonging results"
        call amr_prolong(mype,iopt,nguard)
!        call pf_prolong(mype,iopt,nguard,mg_max_level)

        !call amr_multigrid_block_types(mype,nprocs)
        if(mype.eq.output_pe)print *,"Reinitialising multigrid"
        call amr_mg_init() !NOTE FOR CHRIS, THIS IS NEW THING WHICH STOPPED SEG FAULT
        nodetype_copy(1:lnblocks)=nodetype(1:lnblocks)
        call amr_mg_get_max_refinement(mype,nprocs)
        call amr_multigrid_child_set()
        call reset_nodetype(nprocs,mype,mg_max_level)
        !call amr_mg_init() !NOTE FOR CHRIS, THIS IS OLD THING WHICH NEEDED TO BE MOVED UP
        call amr_multigrid_block_types(mype,nprocs,finegridblocks)

     end if  ! end if test ref
     ! Copy existing soln(s) to rhs
!     call amr_guardcell(mype,iopt,nguard)
     call pf_guardcell(mype,iopt,mg_max_level)
     call amr_mg_runtime_rhs()
     if(mype.eq.output_pe)print *,"Refinement/guardcell update/RHS setting done, start v-cycles"

     ! Removed FMG bit
!     if (time_step-back_step_ct.eq.2) then
!        call output_2d(mype, time_step-back_step_ct)
!        stop
!     endif  


     defect=0.0
     old_defect_1=0.0
     old_defect_2=0.0
     do i_step=1,max_v_cycle_count
        i_step_total=i_step_total+1
        ! Back up copies of old defects
        old_defect_2=old_defect_1
        old_defect_1=defect
        defect=0.0
        call amr_multigrid_v_cycle(mype,nprocs,mg_max_level,defect)

        if(mype.eq.output_pe)then
           ! Second arguement to mod controls output frequency
           if(i_step.eq.1)then
 2001         format ("   Vcycle: ",I8,"  Defect=",E14.7) 
              if(mype.eq.output_pe)write(*,2001) i_step, defect
           else
 2000         format ("   Vcycle: ",I8,"  Defect=",E14.7,"  rate=",F13.7) 
              convrate = old_defect_1/defect
              if(mype.eq.output_pe)write(*,2000) i_step, defect, convrate
           end if
        end if
        call amr_mg_get_rmsres(mype, nprocs)

        if (shampine_test.eq.1) then 
           call shampine_convergence(i_step, failer, mype)
           if(defect.lt.1.0e-12)exit
           if (failer.eq.1) exit           
        else
           ! Final check, if converged, stop cycles
           if(defect.lt.defect_tol)exit
        end if

        ! Start testing defect to see whether it's growing/too big
        if(defect.gt.10.0)then
           if(mype.eq.output_pe)print *,"Defect too large, aborting"
           exit
        end if
        !if(mype.eq.output_pe)write (unit=12,fmt=101) defect
!CEG changed 2 to 20
        if(i_step.gt.20)then
           ! If defects grow, abort
           if(defect.gt.old_defect_1)then
              if(mype.eq.output_pe)print *,"Defects growing, aborting v-cycles"
              exit
           end if
        end if
        if (i_step.ge.10)then
           ! Don't waste time on cycles going nowhere
           exit
        endif
     end do

!     if (time_step-back_step_ct.gt.2) then
!        call ceg_local_timestep(mype, nprocs, local_r)
!     else
        local_r = -1.0
!     endif

     ! Timer end code, and run stat
     if(mype.eq.0) then
        wtime = MPI_Wtime()
        call cpu_time(timenew)
        write(6,*)'Sys/Real solve time was ', timenew-timeprev, wtime-wtimeprev
!        write(6,*)wtime,wtimeprev
        write (133,'(E12.5,X,E12.5,X,6I,X,6I,X,3I)')time, wtime-wtimeprev, finegridblocks(1), &
                      finegridblocks(2), i_step
        timeprev=timenew
        wtimeprev=wtime
        call flush(133)
        print *,"Vcycle: ",i_step,"Defect=",defect
        print *,"Total V-cycles used:",i_step_total,"Avg v-cycles per timestep",(i_step_total*1.0)/&
              (1.0*(time_step-back_step_ct))
!        write (unit=301,fmt=*)time,i_step,(i_step_total*1.0)/(1.0*(time_step-back_step_ct))
     endif
!     call amr_multigrid_block_timer(mype,nprocs)
     ! Start testing final defect
     ! In all failure cases dt->.75*dt and timestep is restarted (failure is lack of convergence or NANs)
     ! If converged in small (<4) cycles dt->10/9*dt then next cycle
     ! If converged in large (>8) cycles dt->9/10*dt then next cycle (note, not reset since it converged)

     failer=0

     call timestep_adjust(defect, local_r, incset, failer, i_step, mype)
     if (failer.eq.1) then
        i_step_total=i_step_total-i_step
        back_step_ct=back_step_ct+1
     endif

     if (failer.eq.0) then
! Some generic terminal output
        if(mod(time_step,1).eq.0)then
           ! Guardcell call should not be needed hence commented
!           call amr_guardcell(mype,iopt,nguard)
           call phase_field_check(mype,time_step-back_step_ct)
        end if
! Chombo output
        if(mod(time_step-back_step_ct,100).eq.-1)then
           if(mype.eq.output_pe)then
              print *,"Outputting current state"
           end if
           ! Write file to disk
           call amr_plotfile_chombo(time_step-back_step_ct)
           write(61,*) 'Step, time', time_step-back_step_ct, time
           flush(61) 
        end if
!CEG adding checkpointing output
        if(mod(time_step-back_step_ct,100).eq.0)then
           if (ndim.eq.2) call output_2d(mype, time_step-back_step_ct)
           if(mype.eq.output_pe)then
              print *,"Checkpoint!"
           end if
           chk_gd = .FALSE.
           chk_checkf = 'hdf5'
           call amr_checkpoint_wr(time_step-back_step_ct, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold)
           if(mype.eq.0) then
              open (unit=99,file="CHK.out")
              write(99,*) time_step-back_step_ct
              close(99)
              call cpu_time(timenew)
              wtime = MPI_Wtime()
              print*,'Checkpointing took ', timenew-timeprev, wtime-wtimeprev
              timeprev=timenew
           endif
        endif

!        if (mode_1.eq.2) return
        mode_1 = 2
!         print*,'mode 1!!!!'
        call nodebound

     else
        if (dt.lt.1.0d-10) then
           if(mype.eq.output_pe) print*,'Tiny dt - aborting'
           call MPI_Abort(MPI_COMM_WORLD, failer, ierr)
        endif
     endif

     chk_restart = 0

!     if (mode_1.eq.1.and.dt.gt.0.02) dt = 0.02
!     if (dt.gt.0.1) dt = 0.1

!     if (mode_1.eq.2) return
!     if (mode_1.eq.2) exit
!     if (time_step-back_step_ct.eq.90) return!585
!     if (time_step-back_step_ct.eq.1000) return!585
!     if (time_step-back_step_ct.eq.10000) return!585

  end do

  if(mype.eq.output_pe)print *,"Mode: ",mode_1
  if(mype.eq.output_pe)print *,"Min level: ",mg_min_lvl
  call amr_multigrid_block_types(mype,nprocs,finegridblocks)
  if(mype.eq.0)close(61)
  if(mype.eq.0)close(125)
  if(mype.eq.0)close(126)
  if(mype.eq.0)close(127)
  if(mype.eq.0)close(128)
  if(mype.eq.0)close(129)
  if(mype.eq.0)close(130)
  if(mype.eq.0)close(131)
  if(mype.eq.0)close(132)

  deallocate(has_children)
  deallocate(nodetype_copy)

  ! print out the final grid state
!  call amr_plotfile_chombo(time_step-back_step_ct)
!  call amr_checkpoint_wr(time_step-back_step_ct, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold)
!  write(61,*) 'Step, time', time_step-back_step_ct, time

!101 format(e19.13)
end subroutine amr_mg_time_dep_control


subroutine timestep_adjust(defect, local_r, incset, failer, i_step, pid)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use solution_parameters
  use paramesh_interfaces!, only :amr_guardcell, & 
       !&                         amr_prolong, & 
       !&                         amr_restrict
  implicit none
  integer, intent(inout) :: incset, failer, i_step, pid
  double precision, intent(inout) :: defect, local_r

! Not done Local error test
  if (local_r.lt.0.0) then
     ! James's methods
     if(defect.eq.0.0)then
        ! If it's = 0, then in all likelyhood either everything's melted, or
        ! the grid's full of NANs
        if(pid.eq.0)print *,"Zero defect - let's ignore for now - treat as tiny"
        if (incset.eq.0) then
           if(pid.eq.0)print *,"Low v-cycle count, increasing dt"
           dtold=dt
           dt=10.0/9.0*dt
           if(dt.gt.max_dt)then
              dt=max_dt
           end if
!           incset=1
           incset=0
        else
           incset = 0
        endif

!        if(pid.eq.0)print *,"Zero defect, asumming failure"
!        call amr_mg_time_step_reset()
!        back_step_ct=back_step_ct+1
!        time=time-dt
!        dt=dt*.75
!        incset = 0
!        failer=1
     else if(defect.gt.1.0E-7)then
        ! Pretty obvious really
        if(pid.eq.0)print *,"Large defect, asumming failure"
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*.75
        incset = 0
        failer=1
     else if((i_step.lt.10.and.defect.lt.defect_tol).or.(shampine_test.eq.1.and.i_step.lt.5))then
        if (incset.eq.0) then
           if(pid.eq.0)print *,"Low v-cycle count, increasing dt"
           dtold=dt
           dt=10.0/9.0*dt
           if(dt.gt.max_dt)then
              dt=max_dt
           end if
!           incset=1
           incset=0
        else
           incset = 0
        endif
     else if((i_step.ge.max_v_cycle_count).or.(defect.gt.defect_tol.and.shampine_test.eq.0))then
!     else if((i_step.ge.10).or.(defect.gt.defect_tol))then
!         if(pid.eq.0)print *,"Timestep decreased removed"
        if(pid.eq.0)print *,"Excessive v-cycle count/lack of convergence, decreasing dt and restarting time step"
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*.75
        incset = 0
        failer=1
     else if(i_step.gt.10)then
        if(pid.eq.0)print *,"High v-cycle count, decreasing dt"   
        dtold=dt
        dt=.9*dt
        incset = 0
     end if
  else     ! Local error estimation techniques
     if(defect.gt.1.0E-7.or.(i_step.ge.max_v_cycle_count).or.(defect.gt.defect_tol.and.shampine_test.eq.0))then
        ! Pretty obvious really
        if(pid.eq.0)print *,"LEE: Large defect, asumming failure", defect
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.5
        incset = 0
        failer=1
     else if (local_r.lt.0.1) then
        if(pid.eq.0)print *,"LEE: Tiny local_r, asumming failure", local_r
        call amr_mg_time_step_reset()
        time=time-dt
        dt=dt*0.5
        incset = 0
        failer=1
     else if (local_r.lt.0.5) then
        if(pid.eq.0)print *,"LEE: Small local_r, halving stepsize", local_r
        dtold=dt
        dt=dt*0.5
        incset = 0
     else if (local_r.lt.0.9) then
        if(pid.eq.0)print *,"LEE: local_r requests decreased timestep", local_r
        dtold=dt
        dt=dt*local_r
        incset = 0
     else if (local_r.gt.1.25.and.dt.lt.max_dt) then
        if(pid.eq.0)print *,"LEE: local_r requests increased timestep", local_r
        dtold=dt
        dt=dt*2.0
        if (dt.gt.max_dt) dt = max_dt
        incset = 1
     endif
  endif

end subroutine timestep_adjust





subroutine nodebound
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none

  integer i, j, k, lb, v


  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        do k=kl_bnd+1,ku_bnd-1
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
             do v = 1, 15
                if (ABS(unk(v,i,j,k,lb)).lt.1.0e-15) then
!                   print *,v,i,j,k,lb,unk(v,i,j,k,lb), ABS(unk(v,i,j,k,lb)
                   unk(v,i,j,k,lb)=0.0
                endif
              end do
            end do
         end do
       end do
     end if
   end do

end subroutine nodebound

subroutine output_2d(pid, fnum)
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use io

  implicit none

  integer pid
  integer i, j, k, v, lb, fnum
  character*80 fname
  character (len=2)  :: proc_string
  character (len=2)  :: grid_string
  character (len=5)  :: fnum_string
  double precision dx
  real x0,y0,p0,n0
  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
  real, dimension(2) :: dH
  integer, dimension(2,8) :: H

  H = RESHAPE(Hdum,shape=(/2,8/))
     
!  write(fname,'(A,I0.2)'),"/tmp/2dout_",lrefine_max,"_",fnum 
   Write (proc_string, '(i2.2)') pid
   Write (grid_string, '(i2.2)') lrefine_max
   Write (fnum_string, '(i5.5)') fnum
!  fname = trim(output_dir) // '2dout_' // grid_string // "_" // fnum_string
!  open (unit=140,file=fname)
 
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)

        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
              do v = 1, 15
                 if (ABS(unk(v,i,j,1,lb)).lt.1.0e-15) unk(v,i,j,1,lb)=0.0
              end do
!              write(140,'(E12.5,X,E12.5,X,E12.5,X,E12.5,X,E12.5,X)') &
!                 bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dx, unk(1,i,j,1,lb), unk(6,i,j,1,lb), unk(11,i,j,1,lb)
           end do
        end do
     end if
  end do
  !  close(140)
  
  fname = trim(output_dir)// '/2dcontour_' // grid_string // "_" // proc_string // "_" // fnum_string
!  fname = 'output/'// '2dcontour_' // grid_string // "_" // !proc_string // "_" // fnum_string
  
  !  write(fname,'(A,I0.2)'),"/tmp/2dcontour_",lrefine_max,"_",fnum
  open (unit=140,file=fname)
  
  
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        
        dx = bsize(1,lb)/real(nxb)
        
        do j=jl_bnd+1,ju_bnd-1
           do i=il_bnd+1,iu_bnd-1
                 if (unk(1,i,j,1,lb).ge.0.0.and.&
             (   unk(1,i+1,j  ,1,lb).lt.0.0.or.&
                 unk(1,i+1,j+1,1,lb).lt.0.0.or.&
                 unk(1,i  ,j+1,1,lb).lt.0.0.or.&
                 unk(1,i-1,j+1,1,lb).lt.0.0.or.&
                 unk(1,i-1,j  ,1,lb).lt.0.0.or.&
                 unk(1,i-1,j-1,1,lb).lt.0.0.or.&
                 unk(1,i  ,j-1,1,lb).lt.0.0.or.&
                 unk(1,i+1,j-1,1,lb).lt.0.0)&
                    )then 
                x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
                y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
                n0=0.0
                p0= unk(1,i,j,1,lb)
                do k=1,8
                   if( unk(1,i+H(1,k),j+H(2,k),1,lb)<n0)then
                       n0=unk(1,i+H(1,k),j+H(2,k)  ,1,lb)
                       dH=H(:,k)
                   end if
                end do
                  
                write(140,'(5(E12.5,X))') &
                (-n0*x0+p0*(x0+dH(1)*dx))/(p0-n0),&
                (-n0*y0+p0*(y0+dH(2)*dx))/(p0-n0),&
                unk(1,i  ,j  ,1,lb),&    
                unk(6,i,j,1,lb), unk(11,i,j,1,lb)
                 
              endif
           end do
        end do
     end if
  end do
  
end subroutine output_2d



                         

