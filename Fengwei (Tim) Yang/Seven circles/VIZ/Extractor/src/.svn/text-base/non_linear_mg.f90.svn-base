!
!
!  Control file for paramesh time dependant non-linear multigrid solver
!  Written by J. Green
!  Based upon: http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr_users_guide.html#main_p
!
!
program phasefield
!#include "paramesh_preprocessor.fh"

  ! include file to define physical qualities of the model and mesh
  use paramesh_dimensions
  use physicaldata
  ! include workspace so interp_mask_work can be set
  use workspace
  ! include file defining the tree
  use tree
  use solution_parameters
  use multigrid_parameters
  ! IMPORTANT: Expose the PARAMESH interface blocks to your code
  use paramesh_interfaces
  use paramesh_mpi_interfaces
  
  ! Set implicit to none, not needed but is sensible
  implicit none

  ! include file required for mpi library.
  include 'mpif.h'
  
  integer pid, noprocs, ierr
  
! sdtart off by setting up Paramesh
  call amr_initialize
!CEG moved from control as needed for initial refinement

#ifdef MPE
  call MPE_Init_log()
#endif
  ! Set values of nprocs and mype with appropriate mpi calls
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, noprocs, ierr)

  call phase_parameter_set()
  call domain_setup(pid, noprocs)

  ! Call main MG loop routine
  call amr_mg_control(pid, noprocs)

  ! FIND A WAY AROUND THIS PROBLEM
  ! If prolong doesn't get called then for some reason the linker fails on paramesh code
  !../libs/libparamesh.a(amr_mg_prolong.o): In function `amr_mg_prolong_':
  !amr_mg_prolong.F90:(.text+0x206): undefined reference to `amr_prolong_'
  !collect2: ld returned 1 exit status
  !problem is probably caused because mg_prolong doesn't have the right 'use' statement
!  do loop_count=1,1
!     if (pid.eq.-1) call amr_prolong(pid, ierr, noprocs) !note arguments are dummy-rubbish as call is never made
!  end do

!  if (shampine_test.eq.1) call shampine_memory_dealloc
  call amr_mg_end()
!  call amr_close()
  
  stop
end program phasefield

subroutine phase_parameter_set

  ! include file to define physical qualities of the model and mesh
  use paramesh_dimensions
  use physicaldata
  ! include workspace so interp_mask_work can be set
  use workspace
  ! include file defining the tree
  use tree
  use solution_parameters
  use multigrid_parameters
  ! IMPORTANT: Expose the PARAMESH interface blocks to your code
  use paramesh_interfaces
  use paramesh_mpi_interfaces
  use multigrid_parameters

  ! Set implicit to none, not needed but is sensible
  implicit none

  ! Run case - variables declared in amr_multigrid_user_modules.f90
  tau0 = 1.0
  w0 = 1.0
  lambda = 2.0               !  Ramirez1 test
!  lambda = 2.0/0.626663849  !  Jan solute test
  delta = -0.55
!  epsilon = 0.03
  epsilon = 0.02
!  epsilon = 0.0171
  ke = 0.15
!  le = 1.0e50
  le = 40.0  !  Ramirez1 test
!  le = 1.0  !  Jan solute test
  mcinf = 0.07

  total_vars=nvar/nunkvbles

  if (total_vars.lt.3) then
     solute=.true.
     if (solute) then
        le = 1.0e50
        mcinf = 1.0-(1.0-ke)*abs(delta)
!        print *, 'Sett mcinf INCORRECTLY for verification test'
!        mcinf = 1.0-(1.0-ke)*(delta)
     else
        le = 1.0
        mcinf = 0.0
     endif
  endif

end subroutine phase_parameter_set

subroutine domain_setup(mype, nprocs)
  ! include file to define physical qualities of the model and mesh
  use paramesh_dimensions
  use physicaldata
  ! include workspace so interp_mask_work can be set
  use workspace
  ! include file defining the tree
  use tree
  use solution_parameters
  use multigrid_parameters
  use time_dep_parameters

  ! IMPORTANT: Expose the PARAMESH interface blocks to your code
  use paramesh_interfaces
  use paramesh_mpi_interfaces
  
  ! Set implicit to none, not needed but is sensible
  implicit none
  
  ! include file required for mpi library.
  include 'mpif.h'
  
  integer, intent(in) :: mype, nprocs
  
  ! VARIABLES NEEDED LOCALLY GO HERE
  integer :: numglobalref=1, numlocalref, tmpI(2)
  integer :: ierr,iopt,nlayers,five,six,loop_count,bc_cond_1
!  real :: grid_min=-3200.0,grid_max=3200.0      ! grid 11
  real :: grid_min=-1600.0,grid_max=1600.0      ! grid 10
!  real :: grid_min=-1200.0,grid_max=1200.0      ! grid 9
!  real :: grid_min=-800.0,grid_max=800.0      ! grid 9
!  real :: grid_min=-400.0,grid_max=400.0      ! grid 8
!  real :: grid_min=-200.0,grid_max=200.0     ! grid 7
!  real :: grid_min=-100.0,grid_max=100.0     ! grid 6
!  real :: grid_min=-50.0,grid_max=50.0       ! grid 5
!  real :: zbdry=6.25
!  real :: zbdry=12.5
!  real :: zbdry=25.0
!  real :: zbdry=800.0
  real :: zbdry=1200.0
!  real :: zbdry=200.0
  ! initialize package 
  ! amr package initialization - initializes all model independent stuff
  integer nBx, nBy, nBz, i, j, k
  real dBx, dBy, dBz

  eighthdmn = 1
  if (eighthdmn.eq.1) grid_min=0.0
  mg_min_lvl = 4

  nBx = 1
  nBy = nBx
  nBz = nBx

  dBx = (grid_max-grid_min)/nBx
  dBy = (grid_max-grid_min)/nBy
  dBz = (grid_max-grid_min)/nBz

  interp_mask_unk(:)=1
  interp_mask_work(:)=1

!  if (shampine_test.eq.1) call shampine_memory_alloc

  ! set a limit on the refinement level
  lrefine_min = 1
!  lrefine_max = 12
!  lrefine_max = 11
  lrefine_max = 10
!  lrefine_max = 9
!  lrefine_max = 8
!  lrefine_max = 7
!  lrefine_max = 6
!  lrefine_max = 5
!  lrefine_max = 4
!  lrefine_max = 3
!  lrefine_max = 2
  
  if (lrefine_max-1.lt.numglobalref) then
     numglobalref = lrefine_max - 1
     numlocalref = 1
  else
     numlocalref = lrefine_max - numglobalref - 1
     if (numlocalref.eq.0) numlocalref = 1
  endif

  ! set bound condition
  bc_cond_1=-21

  ! set up a single block covering the whole cube domain
  lnblocks = 0
  if(mype.eq.0) then
     do k = 1, nBz
        do j = 1, nBy
           do i = 1, nBx
              lnblocks = lnblocks+1
              coord(1,lnblocks) = grid_min+(i-0.5)*dBx
              coord(2,lnblocks) = grid_min+(j-0.5)*dBy

              bnd_box(1,1,lnblocks) = grid_min+(i-1.0)*dBx
              bnd_box(2,1,lnblocks) = grid_min+i*dBx

              bnd_box(1,2,lnblocks) = grid_min+(j-1.0)*dBy
              bnd_box(2,2,lnblocks) = grid_min+j*dBy

              if (ndim.eq.3) then
                 coord(3,lnblocks) = grid_min+(k-0.5)*dBz        ! Fix z coord at middle of layer
!                 bnd_box(1,3,lnblocks) = grid_min+0.5*zbdry
                 bnd_box(1,3,lnblocks) = grid_min+(k-1.0)*dBy
                 bnd_box(2,3,lnblocks) = grid_min+k*dBy
              endif

              bsize(:,lnblocks) = bnd_box(2,:,lnblocks) - bnd_box(1,:,lnblocks)
              nodetype(lnblocks) = 1
              lrefine(lnblocks) = 1
     
              ! all boundary conditions the same
              neigh(2,:,lnblocks) = 0
              if (i.eq.1) then
                 neigh(1,1,lnblocks) = bc_cond_1
              else
                 neigh(1,1,lnblocks) = lnblocks-1
              endif
              if (i.eq.nBx) then
                 neigh(1,2,lnblocks) = bc_cond_1
              else
                 neigh(1,2,lnblocks) = lnblocks+1
              endif
              if (j.eq.1) then
                 neigh(1,3,lnblocks) = bc_cond_1
              else
                 neigh(1,3,lnblocks) = lnblocks-nBx
              endif
              if (j.eq.nBy) then
                 neigh(1,4,lnblocks) = bc_cond_1
              else
                 neigh(1,4,lnblocks) = lnblocks+nBx
              endif

              if (ndim.eq.3) then
                 if (k.eq.1) then
                    neigh(1,5,lnblocks) = bc_cond_1
                 else
                    neigh(1,5,lnblocks) = lnblocks-nBx*nBy
                 endif
                 if (k.eq.nBz) then
                    neigh(1,6,lnblocks) = bc_cond_1
                 else
                    neigh(1,6,lnblocks) = lnblocks+nBx*nBy
                 endif
              endif
              refine(lnblocks)=.true.
           end do
        end do
     end do
     
     min_dx = (bnd_box(2,1,1)-bnd_box(1,1,1))/(2**(lrefine_max-1)*nxb)
     print *,"Min available dx=", min_dx
     max_dt = min_dx
     max_dt = 0.1
     print *,"Min available dy=",(bnd_box(2,2,1)-bnd_box(1,2,1))/(2**(lrefine_max-1)*nyb)
     if (ndim.eq.3) print *,"Min available dz=",(bnd_box(2,3,1)-bnd_box(1,3,1))/(2**(lrefine_max-1)*nzb)

!     do i = 1, lnblocks
!        print*,i,neigh(1,:,i),bnd_box(:,:,i)
!     end do

  endif

!ensure all processors have min_dx and max_dt set
  call MPI_Bcast(min_dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  nucleate_x=grid_min*0
  nucleate_y=grid_min*0
  nucleate_z=grid_min*0
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)          ! parallel synchronization

  ! Set boundary_box and boundary_index information
  ! (NOTE: if you are using periodic boundary conditions, then the arrays boundary_box
  !        and boundary_index can be set to zero).
  call mpi_amr_global_domain_limits() ! This routine sets grid_xmax, grid_xmin, ...
  
  write(6,'(A)') 'Global domain limits'
  write(6,'(F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X,F9.2,X)') grid_xmin,grid_xmax,grid_ymin,grid_ymax,grid_zmin,grid_zmax

  ! x boundaries (boundary boxes 1 and 2)
  boundary_box(1,2:3,1:2) = -1.e30 ! effectively this is -infinity
  boundary_box(2,2:3,1:2) =  1.e30
  boundary_box(1,1,1)     = -1.e30
  boundary_box(2,1,1)     =  grid_xmin ! min dimension of the comp. domain
  boundary_box(1,1,2)     =  grid_xmax ! max     "     "   "   "      "
  boundary_box(2,1,2)     =  1.e30
  
  boundary_index(1)       = bc_cond_1
  boundary_index(2)       = bc_cond_1
  
  ! y boundaries (boundary boxes 3 and 4)
  if (ndim .ge. 2) then
     boundary_box(1,1,3:4) = -1.e30
     boundary_box(2,1,3:4) =  1.e30
     boundary_box(1,3,3:4) = -1.e30
     boundary_box(2,3,3:4) =  1.e30
     boundary_box(1,2,3)   = -1.e30
     boundary_box(2,2,3)   =  grid_ymin
     boundary_box(1,2,4)   =  grid_ymax
     boundary_box(2,2,4)   =  1.e30
     
     boundary_index(3)     = bc_cond_1
     boundary_index(4)     = bc_cond_1
  end if
  
  ! z boundaries
  if (ndim .eq. 3) then
     ! five and six used since in 1d and 2d codes boundary_box only goes to 4
     five = 1+(4*k3d)
     six  = five + k3d
     boundary_box(1,1:2,five:six) = -1.e30
     boundary_box(2,1:2,five:six) =  1.e30
     boundary_box(1,3,five)       = -1.e30
     boundary_box(2,3,five)       = grid_zmin
     boundary_box(1,3,six)        = grid_zmax
     boundary_box(2,3,six)        = 1.e30

     boundary_index(five)         = bc_cond_1
     boundary_index(six)          = bc_cond_1
  end if


!  print *,'do refinement'
  
  do loop_count=1, numglobalref
     if(mype.eq.0) print *,'Pre refinement', loop_count
!     call amr_initial_soln
!     refine(1:lnblocks) = .true.
     
     ! refine grid and apply morton reordering to grid blocks if necessary
!     call amr_refine_derefine

!     call amr_multigrid_block_types(mype,nprocs, tmpI)
  enddo
!  call amr_multigrid_block_types(mype,nprocs, tmpI)

  ! File opening moved to here from amr_mg_time_dep_control
!  if(mype.eq.0)open (unit=134,file="total_data_count.txt")
!  if(mype.eq.0)open (unit=135,file="type_1_data_count.txt")
!  call amr_multigrid_block_timer(mype,nprocs)
  epsilon_tilda=4.0*epsilon/(1.0-3.0*epsilon)
  W0=W0*(1.0-3.0*epsilon)
  !set up initial solution on the grid blocks 
  !call amr_initial_soln()
  !call phase_field_interface_generator()
  ! exchange guardcell information - the call to guardcell also causes the
  ! guard cells at external boundaries to be filled using the user defined
  ! boundary conditions which the user must code into the routine amr_1blk_bcset.
  iopt = 1
  nlayers = nguard
  
  if (numlocalref.gt.0) then
     do loop_count=1,numlocalref
        if(mype.eq.0) print *,'Initial refinement', lrefine_min,lrefine_max, loop_count
!        call amr_initial_soln
     
!        call amr_test_refinement(mype,lrefine_min,lrefine_max, 1)
!        call amr_refine_derefine
!        call amr_prolong(mype,iopt,nlayers)
!  
!        call amr_guardcell(mype,iopt,nlayers)
!
!        call amr_multigrid_block_types(mype,nprocs, tmpI)

     enddo
  endif
    
!  call amr_multigrid_block_types(mype,nprocs, tmpI)
  
!  call phase_field_check(1)
  
end subroutine domain_setup
