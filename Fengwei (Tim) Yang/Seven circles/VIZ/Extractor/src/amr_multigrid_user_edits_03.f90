!#include "paramesh_preprocessor.fh"
subroutine phase_field_check(mype,inumber)
  ! Not really a user edit but since it's specific to the problem
  ! being simulated it's sat here.
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use paramesh_interfaces
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use io

  implicit none
  include 'mpif.h'
  integer, intent(in) :: inumber,mype
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, found, cnt, quadtot, noprocs
  integer, allocatable :: quadprocs(:), quadoff(:)
  double precision, allocatable :: quadlist(:),  quadlistloc(:)
  double precision  lhs(3), mat(3,3)
  double precision :: phase,umin,umax,phimin,phimax,u,phi,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dy,dz,phi_val_1,phi_val_2,phi_val_3,phi_val_4
  double precision :: new_loc, new_rad, temp2, sigmastartemp, locmin, locmax
  double precision :: xl, xm, xr, phil, phim, phir, thetal, thetam, thetar
  double precision :: Divisor, quad_a, quad_b, quad_c
  logical :: on_axisx, on_plane_xy
  character*80 fname
  character (len=2)  :: proc_string

  call MPI_COMM_SIZE(MPI_COMM_WORLD, noprocs, ierr)
  phase=0.0
  umin=1e30
  umax=-1e30
  phimin=1e30
  phimax=-1e30
  mid_x=il_bnd+nxb/2+nguard
  mid_y=jl_bnd+nyb/2+nguard
  if (ndim.eq.3) mid_z=kl_bnd+nzb/2+nguard
  found=0
  if(mype.eq.0)print *,"Phase field summary for time step: ",inumber
! CEG removed profiles
!  if(mype.gt.0)print *,"Warning, solute profiler not yet coded for parallelism"
!  if(mod(inumber,100).eq.0)then
!     write (filename,1001) "solute_profile_",inumber
!     if(mype.eq.0)open (unit=136,file=filename)
!     write (filename,1001) "phase_profile_",inumber
!     if(mype.eq.0)open (unit=137,file=filename)
!     write (filename,1001) "thermal_profile_",inumber
!     if(mype.eq.0)open (unit=138,file=filename)
!  end if
  
  new_rad=0.0
  new_loc=0.0

  sigmastartemp = 0.0
  if (nvar.eq.2.and.solute) found=1
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        if (ndim.eq.3) dz = bsize(3,lb)/real(nzb)
        ! First find out if we're "on axis"
        ! how do we know what the "axis" is?
        ! Simple the nucleate_(x,y,z) parameter holds the co-ords
        ! To be on an axis 2 out of 3 of the block faces must match these co-ords
        ! Only testing lower co-ords so for corner nucleation the nucleus MUST BE at min co-ords for each dimension
        on_axisx=.false.
        
        if(bnd_box(1,2,lb).eq.nucleate_y)then
           ! lower y bnd matches
           if(ndim.eq.2.or.bnd_box(1,3,lb).eq.nucleate_z)then
              ! lower z bnd matches
              ! y and z co-ords match so we're on the x axis
              on_axisx=.true.
              !axis_node_count=axis_node_count+1
           end if
        end if
        j=jl_bnd+nguard
        k=kl_bnd+nguard*k3d

        if(on_axisx)then
           ! On x axis now test whether the interface is within this block
           ! Note testing from il_bnd to iu_bnd
           ! This ensure that if interface is between 2 blocks we still find it
           ! Basically we're checking guardcells as well as computational domain
           if(unk(1,il_bnd+nguard,j,k,lb).gt.0.0.and.&
              unk(1,iu_bnd-nguard+1,j,k,lb).le.0.0)then
              print *, lb, unk(1,il_bnd+nguard,j,k,lb),unk(1,iu_bnd-nguard+1,j,k,lb)
              call phase_loc_rad(lb, new_loc, new_rad)
              found=found+1

              print *,'fnd', unk(1,iu_bnd,j,k,lb), nodetype(lb)

           end if

           if (unk(1,il_bnd+nguard,j,k,lb).gt.-0.9.and.&
               unk(1,iu_bnd-nguard+1,j,k,lb).le.-0.9)then
              if (nvar.eq.3.or.(.not.solute)) then
                 do i = il_bnd+nguard, iu_bnd-nguard
                    if (unk(1,i+1,j,k,lb).le.-0.9) exit
                 end do
                 ! Set up local variables to do this double variable interpolation, first in X for xm then in theta for thetam
                 xl = bnd_box(1,1,lb)+(i-1.5)*dx
                 xr = bnd_box(1,1,lb)+(i+1-1.5)*dx
                 phil = unk(1,i,j,k,lb)
                 phim = -0.9
                 phir = unk(1,i+1,j,k,lb)
                 thetal = unk(1+nunkvbles,i,j,k,lb)
                 thetar = unk(1+nunkvbles,i+1,j,k,lb)

                 xm = xl + (phim-phil)*(xr-xl)/(phir-phil)
                 thetam = thetal + (xm-xl)*(thetar-thetal)/(xr-xl)
                 found=found+1
              endif
           end if
        end if
     end if
!     if (found.eq.2) exit
  end do

!  print *,"Phase_field_check axis node count:",axis_node_count
  call MPI_AllReduce(new_rad,temp2,1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(mype.eq.0)write (6,'(A,A)'), 27, "[36;1;47m"
  if(mype.eq.0)print *,"Tip Radius", temp2
  new_rad = temp2

  call MPI_AllReduce(new_loc,temp2,1,amr_mpi_real,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(mype.eq.0)print *,"Tip position",temp2
  new_loc = temp2

  if(mype.eq.0)write (6,'(A,A)'), 27, "[0m"


  write (proc_string, '(i2.2)') mype
  fname = trim(output_dir) // 'parabolicpts_' // proc_string 

  open(unit=54, file=fname, status='unknown')

! Now for parabolic bit

  locmin = new_loc - 10.0*new_rad
  locmax = new_loc - 2.0*new_rad

!  if(mype.eq.0)write (6,*) 'Parbolic bit : range', locmin, locmax
  if (mype.eq.0) write (6,*) 'Parbolic bit : range', locmin, locmax
  if (locmax.lt.nucleate_x) return
  cnt=0
  do lb=1,lnblocks
!     if(has_children(lb).eqv..false.) then
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        if (ndim.eq.3) dz = bsize(3,lb)/real(nzb)

        ! First find out if we're "on plane"
        on_plane_xy = .false.        
        if (ndim.eq.2) on_plane_xy = .true.
        if(ndim.eq.3.and.abs(bnd_box(1,3,lb)).lt.nucleate_z+1.0e-8) on_plane_xy = .true.

        if(on_plane_xy)then
           j=jl_bnd+nguard
           k=kl_bnd+nguard*k3d
!           if(mype.eq.0)write (6,'(A)') 'on plane'
           if((bnd_box(2,1,lb).gt.locmin .and. bnd_box(1,1,lb).lt.locmax) ) then
! .and.   	&
!               if (mype.eq.0)write (6,*) 'in range', mype, unk(1,il_bnd+nguard,j,k,lb), unk(1,iu_bnd-nguard+1,j,k,lb)
!           if (mype.eq.0)write (6,*) 'test', mype, bnd_box(1,1,lb), bnd_box(2,1,lb)
!           if(mype.eq.0)write (6,'(A)') 'in range'
             if(unk(1,il_bnd+nguard,j,k,lb).gt.0.0.and.unk(1,iu_bnd-nguard+1,ju_bnd-nguard,k,lb).le.0.0)then
!               write (6,*) 'phi diff'
              do j=jl_bnd+nguard,ju_bnd-nguard
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if (unk(1,i,j,k,lb).gt.0.0.and.unk(1,i+1,j,k,lb).le.0.0) then
                       cnt = cnt+1
                       write (54,*) bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dy, i, j
                    endif
                 end do
              end do
             end if
           end if

        end if
     end if
  end do
  cnt = cnt*2

  close(54)

  if (mype.eq.0) then
     allocate(quadprocs(noprocs))  
     allocate(quadoff(noprocs))  
  endif
  call MPI_Gather(cnt, 1, MPI_INTEGER, quadprocs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!  if (mype.eq.0)print *,'Nos per proc'
!  if (mype.eq.0)print *,quadprocs
  if (mype.eq.0) then
     quadoff(1) = 0
     do i = 2, noprocs
        quadtot = quadtot+quadprocs(i-1)
        quadoff(i) = quadtot
     end do
     quadtot = quadtot+quadprocs(noprocs)
!     write(6,*),quadtot
     allocate(quadlist(quadtot))
  endif
  allocate(quadlistloc(cnt))

  cnt=0
  do lb=1,lnblocks
!     if(has_children(lb).eqv..false.) then
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        if (ndim.eq.3) dz = bsize(3,lb)/real(nzb)

        ! First find out if we're "on plane"
        on_plane_xy = .false.        
        if (ndim.eq.2) on_plane_xy = .true.
        if(ndim.eq.3.and.abs(bnd_box(1,3,lb)).lt.nucleate_z+1.0e-8) on_plane_xy = .true.

        if(on_plane_xy)then
           j=jl_bnd+nguard
           k=kl_bnd+nguard*k3d
           if((bnd_box(2,1,lb).gt.locmin .and. bnd_box(1,1,lb).lt.locmax)) then
             if(unk(1,il_bnd+nguard,j,k,lb).gt.0.0.and.unk(1,iu_bnd-nguard+1,ju_bnd-nguard,k,lb).le.0.0)then
              do j=jl_bnd+nguard,ju_bnd-nguard
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if (unk(1,i,j,k,lb).gt.0.0.and.unk(1,i+1,j,k,lb).le.0.0) then
                       ! Note swapped x and y over so parabola is oriented around the right axis
                       quadlistloc(cnt*2+2) = bnd_box(1,1,lb)+(i-1.5)*dx
                       quadlistloc(cnt*2+1) = bnd_box(1,2,lb)+(j-1.5)*dy
                       cnt = cnt+1
                    endif
                 end do
              end do
             end if
           end if
        end if
     end if
  end do

  cnt = cnt*2

  print*,quadprocs
  print*,quadoff
  call MPI_Gatherv(quadlistloc, cnt, MPI_DOUBLE_PRECISION, quadlist, quadprocs, quadoff, &
                   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!  do i = 1, quadtot
!     quadlist(i) = quadlistloc(i)
!  end do

  if (mype.eq.0) then
     lhs(:) = 0.0
     mat(:,:) = 0.0
     mat(1,1) = quadtot/2

     do i = 1, quadtot/2
!        print*,quadlist(2*i-1),quadlist(2*i)
        lhs(1) = lhs(1) + quadlist(2*i)                                    ! y_i
        lhs(2) = lhs(2) + quadlist(2*i)*quadlist(2*i-1)                    ! y_i x_i
        lhs(3) = lhs(3) + quadlist(2*i)*quadlist(2*i-1)*quadlist(2*i-1)    ! y_i x_i^2
        mat(1,2) = mat(1,2) + quadlist(2*i-1)
        mat(1,3) = mat(1,3) + quadlist(2*i-1)*quadlist(2*i-1)
        mat(2,3) = mat(2,3) + quadlist(2*i-1)*quadlist(2*i-1)*quadlist(2*i-1)
        mat(3,3) = mat(3,3) + quadlist(2*i-1)*quadlist(2*i-1)*quadlist(2*i-1)*quadlist(2*i-1)
     end do
     mat(2,1) = mat(1,2)
     mat(3,1) = mat(1,3)
     mat(2,2) = mat(1,3)
     mat(3,2) = mat(2,3)

!     print *, mat(:,:), lhs(:)

! from http://www.microchip.com/forums/tm.aspx?m=461694&mpage=&print=true
     Divisor = mat(1,1)*mat(2,2)*mat(3,3) - mat(1,1)*mat(2,3)*mat(3,2) - mat(1,2)*mat(2,1)*mat(3,3) + &
               mat(1,2)*mat(3,1)*mat(2,3) + mat(2,1)*mat(1,3)*mat(3,2) - mat(1,3)*mat(2,2)*mat(3,1)

     quad_a = (-lhs(1)*mat(2,3)*mat(3,2) + lhs(1)*mat(2,2)*mat(3,3) - lhs(2)*mat(1,2)*mat(3,3) +   &
               lhs(2)*mat(1,3)*mat(3,2) + lhs(3)*mat(1,2)*mat(2,3) - lhs(3)*mat(1,3)*mat(2,2))/Divisor
     quad_b = (lhs(1)*mat(3,1)*mat(2,3) - lhs(1)*mat(2,1)*mat(3,3) + lhs(2)*mat(1,1)*mat(3,3) -    &
               lhs(2)*mat(1,3)*mat(3,1) - lhs(3)*mat(1,1)*mat(2,3) + lhs(3)*mat(2,1)*mat(1,3))/Divisor
     quad_c = (lhs(1)*mat(2,1)*mat(3,2) - lhs(1)*mat(2,2)*mat(3,1) - lhs(2)*mat(1,1)*mat(3,2) +    &
               lhs(2)*mat(1,2)*mat(3,1) + lhs(3)*mat(1,1)*mat(2,2) - lhs(3)*mat(1,2)*mat(2,1))/Divisor

     write(6,'(A,E14.7,A,E13.6,A,E13.6)') 'y=',quad_a, '+x*',quad_b, '+x*x*', quad_c
     ! As below R=( (1+(y')^2)^1.5 )/|y''| and y'=0 hence only need c
     print*,'Parabolic curvature is thus ', 1.0/abs(2*quad_c)
     print*
     print*,'Appropriately scaled quantities are now:'
     print*,'Tip radius actual      ', lambda/a_1*new_rad
     print*,'Tip radius parabolic   ', lambda/a_1/abs(2*quad_c)
     print*,'Tip location actual    ', lambda/a_1*new_loc
     print*,'Tip location parabolic ', lambda/a_1*quad_a
     deallocate(quadlist)
     deallocate(quadprocs)
     deallocate(quadoff)

  end if
  deallocate(quadlistloc)

1001 format(A15,I4) 
end subroutine phase_field_check

! CEG routine for calculating tip location and radius using a series of quadratics

subroutine  phase_loc_rad(lb, cur_loc, cur_rad)
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer i, j, k, lb
  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
  double precision x0, x1, x2, xm1

  dx = bsize(1,lb)/real(nxb)

  j = jl_bnd+nguard
  k = kl_bnd+nguard*k3d
  do i=il_bnd+nguard,iu_bnd-nguard
     if (unk(1,i,j,k,lb).gt.0.0 .and. unk(1,i+1,j,k,lb).le.0.0) then
! Linearly interpolated points
        ! Get x1a - x-position of interface on line y=h/2,z=h/2 using x-y plane
!        x1a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k,lb)/(unk(1,i+1,j,k,lb)-unk(1,i,j,k,lb)))
        ! Get x1b - x-position of interface on line y=3h/2,z=h/2 using x-y plane
!        x1b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k,lb)/(unk(1,i+1,j+1,k,lb)-unk(1,i,j+1,k,lb)))
        ! Get x2a - x-position of interface on line y=h/2,z=3h/2 using x-y plane
!        x2a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k+1,lb)/(unk(1,i+1,j,k+1,lb)-unk(1,i,j,k+1,lb)))
        ! Get x2b - x-position of interface on line y=3h/2,z=3h/2 using x-y plane
!        x2b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k+1,lb)/(unk(1,i+1,j+1,k+1,lb)-unk(1,i,j+1,k+1,lb)))

        x1 = bnd_box(1,1,lb)+(i-1.5)*dx
        x0 = x1-dx
        xm1 = x0-dx
        x2 = x1+dx

        call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k,lb), unk(1,i,j,k,lb), unk(1,i+1,j,k,lb), x1a)
        if (unk(1,i,j+1,k,lb).gt.0.0.or.i.le.2) then
           call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), unk(1,i+1,j+1,k,lb), x1b)
        else
           call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k,lb), unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), x1b)
        endif
        if (ndim.eq.3) then
           if (unk(1,i,j,k+1,lb).gt.0.0.or.i.le.2) then
              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), unk(1,i+1,j,k+1,lb), x2a)
           else
              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j,k+1,lb), unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), x2a)
           endif
           if (unk(1,i,j+1,k+1,lb).gt.0.0.or.i.le.2) then
              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), unk(1,i+1,j+1,k+1,lb), x2b)
           else
              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k+1,lb), unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), x2b)
           endif
        else
           x2a = x1a
           x2b = x1b
        endif
        exit
     endif
  enddo

  b1 = 0.125*(9.0*x1a-x2a)
  b3 = 0.125*(9.0*x1b-x2b)
  b5 = 0.125*(9.0*b1-b3)

  a5 = 0.5*(b3-b1)/(dx*dx)

  cur_loc=b5
! Wolfram-Alpha gives radius of curvature of y=f(x) as R=( (1+(y')^2)^1.5 )/|y''|
  cur_rad = 1.0/abs(2.0*a5)

end subroutine phase_loc_rad

subroutine get_quadratic_interp(x0, x1, x2, y0, y1, y2, xstar)
  double precision x0, x1, x2, y0, y1, y2, xstar, a, b, c

  b = ((y0-y2)*(x1*x1-x2*x2)-(y1-y2)*(x0*x0-x2*x2))/((x0-x2)*(x1*x1-x2*x2)-(x1-x2)*(x0*x0-x2*x2))
  a = (y0-b*(x0-x2)-y2) / (x0*x0-x2*x2)
  c = y2-a*x2*x2-b*x2

! Choose -b MINUS (...)/2a as want left sided root
  xstar = (-b-sqrt(b*b-4.0*a*c))/(2.0*a)

end subroutine get_quadratic_interp

