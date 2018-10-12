!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! ceg_local_timestep
!!!! * Uses method from my thesis
!!!! * Not default at the minute
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This all needs fixing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local error estimation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Operations to do for local error control:
!  1. Back up existing solutions
!  2. Calculate new RHS
!  3. Perform 2/3 V-cycles
!  4. Compute norms to get step-change factor
!  5. Put old solution back

subroutine ceg_local_timestep(pid, noprocs, r)
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  integer, intent(inout) :: pid, noprocs
  double precision, intent(inout) :: r
  integer i_step, npts
  double precision defect, rmsres(total_vars)

  print *,'Into ceg_local_timestep'

!  1. Back up existing solutions
  call ceg_local_backup

!  2. Calculate new RHS
  call ceg_local_rhs

!  3. Perform 2/3 V-cycles
  do i_step=1,3
     call amr_multigrid_v_cycle(pid, noprocs, mg_max_level, defect, rmsres, npts)
     if(pid.eq.0) write(*,*) i_step, defect
     if (defect.lt.1.0e-15) exit
  end do

!  4. Compute norms to get step-change factor
  call ceg_local_factor(pid, noprocs, r)

!  5. Put old solution back
  call ceg_local_restore

  print *,'Out of ceg_local_timestep'

end subroutine ceg_local_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to store the existing solutions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ceg_local_backup()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none

  integer lb, v, k, j, i, u

!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,v,u,k,j,i)
!!!$OMP DO SCHEDULE(DYNAMIC,8) 
  do lb = 1, lnblocks
     do v = 1, total_vars
        u = 1+(v-1)*nunkvbles
        do k=nguard*k3d+1, nzb+nguard*k3d
           do j=nguard+1, nyb+nguard
              do i=nguard+1, nxb+nguard
                 previt(i, j, k, lb, v) = unk(u, i, j, k, lb)
              end do
           end do
        end do
     end do
  end do
!!!$OMP END DO
!!!$OMP END PARALLEL 
  

end subroutine ceg_local_backup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to evaluate the RHS for the local error calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ceg_local_rhs()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none

  integer lb, v, k, j, i
  double precision :: rfactor, rf2, rf3, dti
  rfactor = dt/dtold
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)
  dti=1.0/dt
  
!0.5*dti
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,v,k,j,i)
!!!$OMP DO SCHEDULE(DYNAMIC,8) 
  do lb = 1, lnblocks
     do v = 2,2
        do k=nguard*k3d+1, nzb+nguard*k3d
           do j=nguard+1, nyb+nguard
              do i=nguard+1, nxb+nguard
                 unk(v, i, j, k, lb) =&
                 unk(v, i, j, k, lb) +&
                 (unk(v-1, i, j, k, lb)&
                 - (rf2*unk(v+2,i,j,k,lb)&
                 - rf3*unk(v+3,i,j,k,lb)))
              end do
           end do
        end do
     end do
  end do
!!!$OMP END DO
!!!$OMP END PARALLEL 
  
end subroutine ceg_local_rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to examine the local error calculation and see what step change factor is returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ceg_local_factor(pid, noprocs, r)
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  include 'mpif.h'

  integer, intent(inout) :: pid, noprocs
  integer lb, v, k, j, i, nblocks, npts, nptsG, ierr
  double precision, intent(inout) :: r
  double precision :: omega, atol, rtol, normPhi, normPhiG
  double precision :: rfactor, rf2, rf3, dti
  rfactor = dt/dtold
  rf2 = (rfactor+1.0)*(rfactor+1.0)/(2.0*rfactor+1.0)
  rf3 = rfactor*rfactor/(2.0*rfactor+1.0)

! 1. Need to calculate ATOL based on intergrid differences

  call ceg_local_atol(pid, noprocs, atol)

  rtol = atol

! 2. Loop over fine grid points to calculate norm-Phi
   
  nblocks = 0
  normPhi = 0.0

  do lb = 1, lnblocks
     if (nodetype(lb).eq.1) then
        do v = 1, 1
           do k=nguard*k3d+1, nzb+nguard*k3d
              do j=nguard+1, nyb+nguard
                 do i=nguard+1, nxb+nguard
                     omega = atol + rtol*(rf2*unk(v+3,i,j,k,lb) - rf3*unk(v+4,i,j,k,lb))
                     if (abs(omega).gt.1.0e-16) then
                        normPhi = normPhi +&
                        (unk(v,i,j,k,lb)-previt(i, j, k, lb, v))&
                        *(unk(v,i,j,k,lb)-previt(i, j, k, lb, v))/(omega*omega)
!                        print *, unk(v,i,j,k,lb),previt(i, j, k, lb, v),omega
                     end if
!                    unk(v, i, j, k, lb) = unk(v, i, j, k, lb) + unk(v-1, i, j, k, lb) - (rf2*unk(v+2,i,j,k,lb) - rf3*unk(v+3,i,j,k,lb))
                 end do
              end do
           end do
        end do
        nblocks = nblocks+1
     endif
  end do
  
! 3. Synchronise norm-Phi and N between processors and square root
  if (ndim.eq.2) then
     npts = nblocks * (nxb*nyb) / 4
  else
     npts = nblocks * (nxb*nyb*nzb) / 8
  endif

  call MPI_AllReduce(npts, nptsG, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  call MPI_AllReduce(normPhi, normPhiG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

! 4. Calculate stepchange factor r

  normPhiG = sqrt(normPhiG/nptsG)

! note exponent is (-1)/(k+1) where k is order of method here BDF2 => 2
  r = (2.0*normPhiG)**(-0.333333333333333)

  print *, npts, normPhiG, r


end subroutine ceg_local_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to put saved solutions back into the main unk arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ceg_local_restore()
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none

  integer lb, v, k, j, i

  do lb = 1, lnblocks
     do v = 1, total_vars
        do k=nguard*k3d+1, nzb+nguard*k3d
           do j=nguard+1, nyb+nguard
              do i=nguard+1, nxb+nguard
                 unk(1+(v-1)*nunkvbles, i, j, k, lb) = previt(i, j, k, lb, v)
              end do
           end do
        end do
     end do
  end do
  
end subroutine ceg_local_restore

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation of ATOL based on difference between grids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ceg_local_atol(pid, noprocs, atol)
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  include 'mpif.h'

  integer, intent(inout) :: pid, noprocs
  double precision, intent(inout) :: atol
  integer lb, v, k, j, i, nblocks, N_c, GN_c, ierr, is, ie, js, je, ks, kend, k_f, k_c, i_f, i_c, j_f, j_c
  double precision phidiff, phidiffG, interp

  nblocks=0
  phidiff = 0.0
  interp=0.0
  
  do lb = 1, lnblocks
     if (nodetype(lb).eq.1) then
! Now have fine grid block - need parent block and which child we are
!        parent(1,lb)-- id number      parent(2,lb) -- processor
!        print *, parent(1,lb), which_child(lb)
        if (parent(2,lb).ne.pid) then
! add this one to the list
           print *,'Not done interprocessor bits yet'
        else
           call ceg_get_parent_match_indicies(which_child(lb), is, ie, js, je, ks, kend)
           do k=1, nzb, 2
              k_f = nguard*k3d+k
              k_c = ks+k3d*(k/2)
              do j=nguard+1, nyb+nguard, 2
                 j_f = nguard+j
                 j_c = js+k2d*j/2
                 do i=1, nxb, 2
                    i_f = nguard+i
                    i_c = is+i/2
                    if (ndim.eq.3) then
                       interp = 0.125*(unk(1, i, j, k, lb)+unk(1, i+1, j, k, lb)+unk(1, i, j+1, k, lb)+unk(1, i+1, j+1, k, lb) + &
                                  unk(1, i, j, k+1, lb)+unk(1, i+1, j, k+1, lb)+unk(1, i, j+1, k+1, lb)+unk(1, i+1, j+1, k+1, lb))
                    else
                       interp = 0.25*(unk(1, i, j, k, lb)+unk(1, i+1, j, k, lb)+unk(1, i, j+1, k, lb)+unk(1, i+1, j+1, k, lb))
                    endif
                    phidiff = phidiff + (interp - unk(1, i, j, k, parent(1,lb)))*(interp - unk(1, i, j, k, parent(1,lb)))
                 end do
              end do
           end do

           nblocks = nblocks+1
        endif

!        nblocks = nblocks+1
     endif
  end do

  if (ndim.eq.2) then
     N_c = nblocks * (nxb*nyb) / 4
  else
     N_c = nblocks * (nxb*nyb*nzb) / 8
  endif

  call MPI_AllReduce(N_c, GN_c, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  call MPI_AllReduce(phidiff, phidiffG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  atol = 0.1*sqrt(GN_c*phidiffG)
  print *, N_c, phidiffG, atol


end subroutine ceg_local_atol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ceg_get_parent_match_indicies(switch, is, ie, js, je, ks, kend)
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  integer, intent(inout) :: switch, is, ie, js, je, ks, kend

  select case(switch)
     case (1)                      !  xmin ymin  zmin
        is = 1+nguard
        ie = nguard+nxb/2
        js = 1+nguard
        je = nguard+nyb/2
        ks = 1+nguard*k3d
        kend = (nguard+nzb/2)*k3d
     case (2)                      !  xmax ymin  zmin
        is = 1+nguard+nxb/2
        ie = nguard+nxb
        js = 1+nguard
        je = nguard+nyb/2
        ks = 1+nguard*k3d
        kend = (nguard+nzb/2)*k3d
     case (3)                      !  xmin ymax  zmin
        is = 1+nguard
        ie = nguard+nxb/2
        js = 1+nguard+nyb/2
        je = nguard+nyb
        ks = 1+nguard*k3d
        kend = (nguard+nzb/2)*k3d
     case (4)                      !  xmax ymax  zmin
        is = 1+nguard+nxb/2
        ie = nguard+nxb
        js = 1+nguard+nyb/2
        je = nguard+nyb
        ks = 1+nguard*k3d
        kend = (nguard+nzb/2)*k3d
     case (5)                      !  xmin ymin  zmax
        is = 1+nguard
        ie = nguard+nxb/2
        js = 1+nguard
        je = nguard+nyb/2
        ks = 1+(nguard+nzb/2)*k3d
        kend = (nguard+nzb)*k3d
     case (6)                      !  xmax ymin  zmax
        is = 1+nguard+nxb/2
        ie = nguard+nxb
        js = 1+nguard
        je = nguard+nyb/2
        ks = 1+(nguard+nzb/2)*k3d
        kend = (nguard+nzb)*k3d
     case (7)                      !  xmin ymax  zmax
        is = 1+nguard
        ie = nguard+nxb/2
        js = 1+nguard+nyb/2
        je = nguard+nyb
        ks = 1+(nguard+nzb/2)*k3d
        kend = (nguard+nzb)*k3d
     case (8)                      !  xmax ymax  zmax
        is = 1+nguard+nxb/2
        ie = nguard+nxb
        js = 1+nguard+nyb/2
        je = nguard+nyb
        ks = 1+(nguard+nzb/2)*k3d
        kend = (nguard+nzb)*k3d
  end select

end subroutine
