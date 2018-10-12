!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_output_perstep                          REQUIRED
!!!!  * Post V-cycles anything to be done every step
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_output_occasional                       REQUIRED
!!!!  * Post V-cycles anything to be done every 'output_frequency' steps
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  phase_field_check
!!!!  * Main routine for obtaining tip location and radius
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  phase_loc_rad
!!!!  * Computes tip location and radius using quadratic interpolation
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  get_quadratic_interp
!!!!  * Does quadratic interpolation of supplied data
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  pf_output_2d
!!!!  * Prints points showing 2-d phase field boundary
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Generic routines for output

subroutine app_output_perstep(pid, noprocs, inumber)
  implicit none
  integer, intent(in) :: inumber,pid, noprocs

!--tim  call phase_field_check(pid, inumber)

  return
end subroutine app_output_perstep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_output_occasional(pid, noprocs, inumber)
  use paramesh_dimensions
  use time_dep_parameters
  use paramesh_interfaces
  use checkpoint_parameters
  implicit none
  integer, intent(in) :: inumber, pid, noprocs
  logical :: chombo = .false.
  logical chk_gd

  ! 2-d points around front
  if (ndim.eq.2) call pf_output_2d(pid, inumber)

  ! Chombo output
  if (chombo)then
     if(pid.eq.0)then
        print *,"Output Chombo file"
     end if
     ! Write file to disk
     call amr_plotfile_chombo(inumber)
     write(61,*) 'Step, time', inumber, time
     flush(61) 
  end if

  ! Checkpointing output
  if(pid.eq.0) print *,"Checkpoint!"
  chk_gd = .FALSE.
  chk_checkf = 'hdf5'
  call amr_checkpoint_wr(inumber, user_attr_1=dt, user_attr_2=time, user_attr_3=dtold)
  if (pid.eq.0) then
     open (unit=99,file="CHK.out")
     write(99,*) inumber
     close(99)
  endif

  return
end subroutine app_output_occasional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!subroutine phase_field_check(pid, inumber)
  ! Not really a user edit but since it's specific to the problem
  ! being simulated it's sat here.
!!!  use paramesh_dimensions
!!!  use physicaldata
!!!  use tree
!!!  use paramesh_comm_data ! Needed for definition of amr_mpi_real
!!!  use paramesh_interfaces
!!!  use solution_parameters
!!!  use time_dep_parameters
!!!  use multigrid_parameters
!!!  implicit none
!!!  include 'mpif.h'
!!!  integer, intent(in) :: inumber,pid
!!!  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, found
!!!  double precision :: phase,umin,umax,phimin,phimax,u,phi,temp,gradphi,gradphi2,theta1,theta2
!!!  double precision :: dx,dist,phi_val_1,phi_val_2,phi_val_3,phi_val_4
!!!  double precision :: new_loc, new_rad, temp2
!!!  double precision :: xl, xm, xr, phil, phim, phir, thetal, thetam, thetar, sigmastartemp
!!!  logical :: on_axisx=.false.
!  character (len=19) :: filename
  !print *,"Guardcell axis node count:",axis_node_count
  !axis_node_count=0
!!!  phase=0.0
!!!  umin=1e30
!!!  umax=-1e30
!!!  phimin=1e30
!!!  phimax=-1e30
!!!  mid_x=il_bnd+nxb/2+nguard
!!!  mid_y=jl_bnd+nyb/2+nguard
!!!  if (ndim.eq.3) mid_z=kl_bnd+nzb/2+nguard
!!!  found=0
!!!  if(pid.eq.0)print *,"Phase field summary for time step: ",inumber
! CEG removed profiles
!  if(pid.gt.0)print *,"Warning, solute profiler not yet coded for parallelism"
!  if(mod(inumber,100).eq.0)then
!     write (filename,1001) "solute_profile_",inumber
!     if(pid.eq.0)open (unit=136,file=filename)
!     write (filename,1001) "phase_profile_",inumber
!     if(pid.eq.0)open (unit=137,file=filename)
!     write (filename,1001) "thermal_profile_",inumber
!     if(pid.eq.0)open (unit=138,file=filename)
!  end if
  

!!!  dist = 0.0
!!!  sigmastartemp = 0.0
!!!  if (nvar.eq.2.and.solute) found=1
!!!  do lb=1,lnblocks
!!!     if(nodetype(lb).eq.1) then
!!!        dx = bsize(1,lb)/real(nxb)
        ! First find out if we're "on axis"
        ! how do we know what the "axis" is?
        ! Simple the nucleate_(x,y,z) parameter holds the co-ords
        ! To be on an axis 2 out of 3 of the block faces must match these co-ords
        ! Only testing lower co-ords so for corner nucleation the nucleus MUST BE at min co-ords for each dimension
!!!        on_axisx=.false.
        
!!!        if(bnd_box(1,2,lb).eq.nucleate_y)then
           ! lower y bnd matches
!!!           if(ndim.eq.2.or.bnd_box(1,3,lb).eq.nucleate_z)then
              ! lower z bnd matches
              ! y and z co-ords match so we're on the x axis
!!!              on_axisx=.true.
              !axis_node_count=axis_node_count+1
!!!           end if
!!!        end if
!!!        j=jl_bnd+nguard
!!!        k=kl_bnd+nguard*k3d

!!!        if(on_axisx)then
           ! On x axis now test whether the interface is within this block
           ! Note testing from il_bnd to iu_bnd
           ! This ensure that if interface is between 2 blocks we still find it
           ! Basically we're checking guardcells as well as computational domain
!!!           if(unk(1,il_bnd+nguard,j,k,lb).gt.0.0.and.&
!!!              unk(1,iu_bnd-nguard+1,j,k,lb).le.0.0)then
!              print *, lb, unk(1,il_bnd+nguard,j,k,lb),unk(1,iu_bnd-nguard+1,j,k,lb)
!!!              call phase_loc_rad(lb, new_loc, new_rad)
!!!              found=found+1

!              print *,'fnd', unk(1,iu_bnd,j,k,lb)

!!!           end if

!!!           if (unk(1,il_bnd+nguard,j,k,lb).gt.-0.9.and.&
!!!               unk(1,iu_bnd-nguard+1,j,k,lb).le.-0.9)then
!!!              if (nvar.eq.3.or.(.not.solute)) then
!!!                 do i = il_bnd+nguard, iu_bnd-nguard
!!!                    if (unk(1,i+1,j,k,lb).le.-0.9) exit
!!!                 end do
                 ! Set up local variables to do this double variable interpolation, first in X for xm then in theta for thetam
!!!                 xl = bnd_box(1,1,lb)+(i-1.5)*dx
!!!                 xr = bnd_box(1,1,lb)+(i+1-1.5)*dx
!!!                 phil = unk(1,i,j,k,lb)
!!!                 phim = -0.9
!!!                 phir = unk(1,i+1,j,k,lb)
!!!                 thetal = unk(1+nunkvbles,i,j,k,lb)
!!!                 thetar = unk(1+nunkvbles,i+1,j,k,lb)

!!!                 xm = xl + (phim-phil)*(xr-xl)/(phir-phil)
!!!                 thetam = thetal + (xm-xl)*(thetar-thetal)/(xr-xl)
!!!                 found=found+1
!!!              endif
!!!           end if
!!!        end if
!!!     end if
!     if (found.eq.2) exit
!!!  end do

!  print *,"Phase_field_check axis node count:",axis_node_count
!!!  call MPI_REDUCE(new_rad,temp2,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!!!  if(pid.eq.0)write (6,'(A,A)') 27, "[36;1;47m"
!!!  if(pid.eq.0)print *,"Tip Radius",temp2
!!!  if(pid.eq.0)write(unit=125,fmt=*) time, temp2, dt

!!!  call MPI_REDUCE(dist,temp,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!!!  call MPI_REDUCE(new_loc,temp2,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!!!  if(pid.eq.0)print *,"Tip position",temp2
!!!  if(pid.eq.0)write(unit=126,fmt=*) time, temp2
!  call MPI_REDUCE(phase,temp,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!  if(pid.eq.0)print *,"Amount of solid",temp

!!!  call MPI_REDUCE(thetam,temp,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
!!!  if ((nvar.eq.3.or.(.not.solute)).and.pid.eq.0)print *,"Sigma star temperature", temp

!!!  if(pid.eq.0)write (6,'(A,A)') 27, "[0m"
!!!  if(pid.eq.0)call flush(125)
!!!  if(pid.eq.0)call flush(126)
! CEG removed profiles and p2, p3, p4, out1, out2, out3
!!!1001 format(A15,I4) 
!!!end subroutine phase_field_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CEG routine for calculating tip location and radius using a series of quadratics

subroutine  phase_loc_rad(lb, cur_loc, cur_rad)
!!--tim  use paramesh_dimensions
!!--tim  use physicaldata
!!--tim  use tree
!!--tim  implicit none
!!--tim  integer i, j, k, lb
!!--tim  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
!!--tim  double precision x0, x1, x2, xm1

!!--tim  dx = bsize(1,lb)/real(nxb)

!!--tim  j = jl_bnd+nguard
!!--tim  k = kl_bnd+nguard*k3d
!!--tim  do i=il_bnd+nguard,iu_bnd-nguard
!!--tim     if (unk(1,i,j,k,lb).gt.0.0 .and. unk(1,i+1,j,k,lb).le.0.0) then
! Linearly interpolated points
        ! Get x1a - x-position of interface on line y=h/2,z=h/2 using x-y plane
!        x1a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k,lb)/(unk(1,i+1,j,k,lb)-unk(1,i,j,k,lb)))
        ! Get x1b - x-position of interface on line y=3h/2,z=h/2 using x-y plane
!        x1b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k,lb)/(unk(1,i+1,j+1,k,lb)-unk(1,i,j+1,k,lb)))
        ! Get x2a - x-position of interface on line y=h/2,z=3h/2 using x-y plane
!        x2a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k+1,lb)/(unk(1,i+1,j,k+1,lb)-unk(1,i,j,k+1,lb)))
        ! Get x2b - x-position of interface on line y=3h/2,z=3h/2 using x-y plane
!        x2b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k+1,lb)/(unk(1,i+1,j+1,k+1,lb)-unk(1,i,j+1,k+1,lb)))

!!--tim        x1 = bnd_box(1,1,lb)+(i-1.5)*dx
!!--tim        x0 = x1-dx
!!--tim        xm1 = x0-dx
!!--tim        x2 = x1+dx

!!--tim        call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k,lb), unk(1,i,j,k,lb), unk(1,i+1,j,k,lb), x1a)
!!--tim        if (unk(1,i,j+1,k,lb).gt.0.0.or.i.le.2) then
!!--tim           call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), unk(1,i+1,j+1,k,lb), x1b)
!!--tim        else
!!--tim           call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k,lb), unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), x1b)
!!--tim        endif
!!--tim        if (ndim.eq.3) then
!!--tim           if (unk(1,i,j,k+1,lb).gt.0.0.or.i.le.2) then
!!--tim              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), unk(1,i+1,j,k+1,lb), x2a)
!!--tim           else
!!--tim              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j,k+1,lb), unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), x2a)
!!--tim           endif
!!--tim           if (unk(1,i,j+1,k+1,lb).gt.0.0.or.i.le.2) then
!!--tim              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), unk(1,i+1,j+1,k+1,lb), x2b)
!!--tim           else
!!--tim              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k+1,lb), unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), x2b)
!!--tim           endif
!!--tim        else
!!--tim           x2a = x1a
!!--tim           x2b = x1b
!!--tim        endif
!!--tim        exit
!!--tim     endif
!!--tim  enddo

!!--tim  b1 = 0.125*(9.0*x1a-x2a)
!!--tim  b3 = 0.125*(9.0*x1b-x2b)
!!--tim  b5 = 0.125*(9.0*b1-b3)

!!--tim  a5 = 0.5*(b3-b1)/(dx*dx)

!!--tim  cur_loc=b5
! Wolfram-Alpha gives radius of curvature of y=f(x) as R=( (1+(y')^2)^1.5 )/|y''|
!!--tim  cur_rad = 1.0/abs(2.0*a5)

end subroutine phase_loc_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_quadratic_interp(x0, x1, x2, y0, y1, y2, xstar)
!!--tim  double precision x0, x1, x2, y0, y1, y2, xstar, a, b, c

!!--tim  b = ((y0-y2)*(x1*x1-x2*x2)-(y1-y2)*(x0*x0-x2*x2))/((x0-x2)*(x1*x1-x2*x2)-(x1-x2)*(x0*x0-x2*x2))
!!--tim  a = (y0-b*(x0-x2)-y2) / (x0*x0-x2*x2)
!!--tim  c = y2-a*x2*x2-b*x2

! Choose -b MINUS (...)/2a as want left sided root
!!--tim  xstar = (-b-sqrt(b*b-4.0*a*c))/(2.0*a)

end subroutine get_quadratic_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rewritten by Peter Bollada

subroutine pf_output_2d(pid, fnum)
!!--tim  use time_dep_parameters
!!--tim  use solution_parameters
!!--tim  use multigrid_parameters
!!--tim  use paramesh_dimensions
!!--tim  use physicaldata
!!--tim  use tree
!!--tim  use io

!!--tim  implicit none

!!--tim  integer pid
!!--tim  integer i, j, k, v, lb, fnum
!!--tim  character*80 fname
!!--tim  character (len=2)  :: proc_string
!!--tim  character (len=2)  :: grid_string
!!--tim  character (len=5)  :: fnum_string
!!--tim  double precision dx
!!--tim  double precision x0,y0,p0,n0
!!--tim  integer, dimension(16) :: Hdum= (/ 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1,  0, -1 , 1, -1 /) 
!!--tim  double precision, dimension(2) :: dH
!!--tim  integer, dimension(2,8) :: H

!!--tim  H = RESHAPE(Hdum,shape=(/2,8/))
     
!  write(fname,'(A,I0.2)'),"/tmp/2dout_",lrefine_max,"_",fnum 
!!--tim   Write (proc_string, '(i2.2)') pid
!!--tim   Write (grid_string, '(i2.2)') lrefine_max
!!--tim   Write (fnum_string, '(i5.5)') fnum
!  fname = trim(output_dir) // '2dout_' // grid_string // "_" // fnum_string
!  open (unit=140,file=fname)
 
!!--tim  do lb=1,lnblocks
!!--tim     if(nodetype(lb).eq.1) then
!!--tim        dx = bsize(1,lb)/real(nxb)

!!--tim        do j=jl_bnd+1,ju_bnd-1
!!--tim           do i=il_bnd+1,iu_bnd-1
!!--tim              do v = 1, 15
!!--tim                 if (ABS(unk(v,i,j,1,lb)).lt.1.0e-15) unk(v,i,j,1,lb)=0.0
!!--tim              end do
!              write(140,'(E12.5,X,E12.5,X,E12.5,X,E12.5,X,E12.5,X)') &
!                 bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dx, unk(1,i,j,1,lb), unk(6,i,j,1,lb), unk(11,i,j,1,lb)
!!--tim           end do
!!--tim        end do
!!--tim     end if
!!--tim  end do
  !  close(140)
  
!!--tim  fname = trim(output_dir)// '/2dcontour_' // grid_string // "_" // proc_string // "_" // fnum_string
!  fname = 'output/'// '2dcontour_' // grid_string // "_" // !proc_string // "_" // fnum_string
  
  !  write(fname,'(A,I0.2)'),"/tmp/2dcontour_",lrefine_max,"_",fnum
!!--tim  open (unit=140,file=fname)
  
  
!!--tim  do lb=1,lnblocks
!!--tim     if(nodetype(lb).eq.1) then
        
!!--tim        dx = bsize(1,lb)/real(nxb)
        
!!--tim        do j=jl_bnd+1,ju_bnd-1
!!--tim           do i=il_bnd+1,iu_bnd-1
!!--tim                 if (unk(1,i,j,1,lb).ge.0.0.and.&
!!--tim             (   unk(1,i+1,j  ,1,lb).lt.0.0.or.&
!!--tim                 unk(1,i+1,j+1,1,lb).lt.0.0.or.&
!!--tim                 unk(1,i  ,j+1,1,lb).lt.0.0.or.&
!!--tim                 unk(1,i-1,j+1,1,lb).lt.0.0.or.&
!!--tim                 unk(1,i-1,j  ,1,lb).lt.0.0.or.&
!!--tim                 unk(1,i-1,j-1,1,lb).lt.0.0.or.&
!!--tim                 unk(1,i  ,j-1,1,lb).lt.0.0.or.&
!!--tim                 unk(1,i+1,j-1,1,lb).lt.0.0)&
!!--tim                    )then 
!!--tim                x0 =  bnd_box(1,1,lb)+(i-1.5)*dx
!!--tim                y0 =  bnd_box(1,2,lb)+(j-1.5)*dx
!!--tim                n0=0.0
!!--tim                p0= unk(1,i,j,1,lb)
!!--tim                do k=1,8
!!--tim                   if( unk(1,i+H(1,k),j+H(2,k),1,lb)<n0)then
!!--tim                       n0=unk(1,i+H(1,k),j+H(2,k)  ,1,lb)
!!--tim                       dH=H(:,k)
!!--tim                   end if
!!--tim                end do
                  
!!--tim                write(140,'(5(E12.5,X))') &
!!--tim                (-n0*x0+p0*(x0+dH(1)*dx))/(p0-n0),&
!!--tim                (-n0*y0+p0*(y0+dH(2)*dx))/(p0-n0),&
!!--tim                unk(1,i  ,j  ,1,lb),&    
!!--tim                unk(6,i,j,1,lb), unk(11,i,j,1,lb)
                 
!!--tim              endif
!!--tim           end do
!!--tim        end do
!!--tim     end if
!!--tim  end do
  
end subroutine pf_output_2d


