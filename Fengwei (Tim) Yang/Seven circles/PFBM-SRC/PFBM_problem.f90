!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_parameter_set                                            REQUIRED
!!!!  * Called at start of program execution.
!!!!  * Paramesh setup needed first to get nvar set
!!!!  * Sets case specific parameters
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!  
!!!!  app_domain_setup                                             REQUIRED
!!!!  * Sets up the bounding boxes and boundary conditions for the initial 
!!!!    coarsest block distribution
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!
!  Control file for paramesh time dependant non-linear multigrid solver
!  Written by J. Green
!  Based upon: http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr_users_guide.html#main_p
!

subroutine app_parameter_set
  ! include modules for governing solution parameters
  use multigrid_parameters
  use generic_parameters
  ! include file to define physical quantities of the model
  use solution_parameters
  ! include paramesh data for nvar
  use paramesh_dimensions

  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set parameters that affect global generic parameters
  verbose=1            ! How much teerminal output
  local_adaptation = 0 ! Are we using local_adaptation?  0=no 1=yes 2=user supplied
  numglobalref=1       ! How many levels of global refinement are needed initially
  !-- growing domain will not be used here.
  grow = 0             ! Are we letting the domain grow as adaptation spreads?  1=yes
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  ! Do we really need all these values?
  total_vars=nvar/nunkvbles


end subroutine app_parameter_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_domain_setup(pid, noprocs)
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
  use generic_parameters

  ! IMPORTANT: Expose the PARAMESH interface blocks to your code
  use paramesh_interfaces
  use paramesh_mpi_interfaces
  
  ! Set implicit to none
  implicit none
  
  ! include file required for mpi library.
  include 'mpif.h'
  
  integer, intent(in) :: pid, noprocs
  ! VARIABLES NEEDED LOCALLY GO HERE
  integer :: ierr
  integer nBx, nBy, nBz, i, j, k
  double precision :: dBx, dBy, dBz, pi, pikk=1.0
!  double precision :: grid_min=-3200.0,grid_max=3200.0      ! grid 11
!  double precision :: grid_min=-1600.0,grid_max=1600.0      ! grid 10
!  double precision :: grid_min=-1200.0,grid_max=1200.0      ! grid 9
!  double precision :: grid_min=-800.0,grid_max=800.0      ! grid 9
!  double precision :: grid_min=-400.0,grid_max=400.0      ! grid 8
!  double precision :: grid_min=-200.0,grid_max=200.0     ! grid 7
!  double precision :: grid_min=-100.0,grid_max=100.0     ! grid 6   !Recommended coarse block resolution
!  double precision :: grid_min=-50.0,grid_max=50.0       ! grid 5
!  double precision :: grid_min=-25.0,grid_max=25.0       ! grid 4
!  double precision :: grid_min=-12.50,grid_max=12.50       ! grid 3
  pi=4.0*datan(pikk)

  grid_min=0.0
  grid_max=2.0*pi

  eighthdmn = 1

  ep = 0.025
  weight = 1.00
  dt = 0.01/8.0
  dtold = dt                    ! initializing previous dt
  solve_count_max = 50           ! coarsest grid solver
  smooth_count_max = 5          ! pre/post smoothers
  defect_too_big = 1.0e-2       ! tolerant for adaptive time-stepping
  min_dt = 1.0e-14              ! tolerant for adaptive time-stepping
  simulation_time = 10000.0     ! ending time
  low_vcycle_count = 6          ! tolerant for controlling V-cycle
  max_v_cycle_count = 10        ! max V-cycle per one time step
  output_rate = 3200000         ! one output on each many time-step
  
  
  if (eighthdmn.eq.1) grid_min = 0.0
  mg_min_lvl = 1                ! coarsest grid for multigrid

  nBx = 1
  nBy = nBx
  nBz = nBx

  dBx = (grid_max-grid_min)/nBx
  dBy = (grid_max-grid_min)/nBy
  dBz = (grid_max-grid_min)/nBz

  if (ndim.eq.2) nBz = 1


  ! set a limit on the refinement level
  lrefine_min = 1

!  lrefine_max = 12
!  lrefine_max = 11
!  lrefine_max = 10
!  lrefine_max = 9
!  lrefine_max = 8
!  lrefine_max = 7
!  lrefine_max = 6                       ! For G6 this is dx=0.39
!  lrefine_max = 5                       ! For G6 this is dx=0.78
!  lrefine_max = 4
!  lrefine_max = 3
  lrefine_max = 3

   numglobalref=2

  if(interpolating.eq.1) lrefine_max = 7


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!          Leave untouched from here down for phase field code                                       !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


! CEG Initialise my startlist for the blocklist for that grid level
  allocate(block_starts(lrefine_max+1))
  
  ! set up a single block covering the whole cube domain
  lnblocks = 0
  block_starts(1) = 1
  if (pid.eq.0) then
     do k = 1, nBz
        do j = 1, nBy
           do i = 1, nBx
              lnblocks = lnblocks+1
              coord(1,lnblocks) = nBx*grid_min + (i-0.5)*dBx
              coord(2,lnblocks) = nBy*grid_min + (j-0.5)*dBy

              bnd_box(1,1,lnblocks) = nBx*grid_min + (i-1.0)*dBx
              bnd_box(2,1,lnblocks) = nBx*grid_min + i*dBx

              bnd_box(1,2,lnblocks) = nBy*grid_min + (j-1.0)*dBy
              bnd_box(2,2,lnblocks) = nBy*grid_min + j*dBy

              if (ndim.eq.3) then
                 coord(3,lnblocks) = nBz*grid_min+(k-0.5)*dBz        ! Fix z coord at middle of layer
                 bnd_box(1,3,lnblocks) = nBz*grid_min+(k-1.0)*dBy
                 bnd_box(2,3,lnblocks) = nBz*grid_min+k*dBy
              endif

              bsize(:,lnblocks) = bnd_box(2,:,lnblocks) - bnd_box(1,:,lnblocks)
              nodetype(lnblocks) = 1
              lrefine(lnblocks) = 1
     
              ! boundary conditions need to represent whether far field or symmetry
!              neigh(2,:,lnblocks) = 0
              if (i.eq.1.and.eighthdmn.eq.1) then
                 neigh(:,1,lnblocks) = 1.0 !bc_cond_far
              else if (i.eq.1) then
                 neigh(:,1,lnblocks) = 1.0 !bc_cond_far
              else
                 neigh(1,1,lnblocks) = lnblocks-1
                 neigh(2,2,lnblocks) = 0
              endif
              if (i.eq.nBx) then
                 neigh(:,2,lnblocks) = 1.0 !bc_cond_far
              else
                 neigh(1,2,lnblocks) = lnblocks+1
                 neigh(2,2,lnblocks) = 0
              endif

              if (j.eq.1.and.eighthdmn.eq.1) then
                 neigh(:,3,lnblocks) = 1.0 !bc_cond_sym
              else if (j.eq.1) then
                 neigh(:,3,lnblocks) = 1.0 !bc_cond_sym
              else
                 neigh(1,3,lnblocks) = lnblocks-nBx
                 neigh(2,3,lnblocks) = 0
              endif
              if (j.eq.nBy) then
                 neigh(:,4,lnblocks) = 1.0 !bc_cond_sym
              else
                 neigh(1,4,lnblocks) = lnblocks+nBx
                 neigh(2,4,lnblocks) = 0
              endif

              if (ndim.eq.3) then
                 if (k.eq.1.and.eighthdmn.eq.1) then
                    neigh(:,5,lnblocks) = bc_cond_sym
                 else if (k.eq.1) then
                    neigh(:,5,lnblocks) = bc_cond_sym
                 else
                    neigh(1,5,lnblocks) = lnblocks-nBx*nBy
                    neigh(2,5,lnblocks) = 0
                 endif
                 if (k.eq.nBz) then
                    neigh(:,6,lnblocks) = bc_cond_sym
                 else
                    neigh(1,6,lnblocks) = lnblocks+nBx*nBy
                    neigh(2,6,lnblocks) = 0
                 endif
              endif
              neigh(2,:,lnblocks) = 0
              refine(lnblocks)=.true.
           end do
        end do
     end do
     
     min_dx = (bnd_box(2,1,1)-bnd_box(1,1,1))/(2**(lrefine_max-1)*nxb)
     if (pid.eq.0.and.verbose.ge.2) print *,"Min available dx=", min_dx
     max_dt = min_dx
     max_dt = 0.1
  endif

  ! Set boundary conditions
  if (eighthdmn.eq.1) then
     boundary_index(1)       = 1.0 !bc_cond_far
  else  
     boundary_index(1)       = 1.0 !bc_cond_far
  endif
  boundary_index(2)       = 1.0 !bc_cond_far
  
  ! y boundaries (boundary boxes 3 and 4)
  if (ndim .ge. 2) then
     if (eighthdmn.eq.1) then
        boundary_index(3)       = 1.0 !bc_cond_sym
     else  
        boundary_index(3)       = 1.0 !bc_cond_sym
     endif
     boundary_index(4)     = 1.0 !bc_cond_sym
  end if
  
  ! z boundaries
  if (ndim .eq. 3) then
     if (eighthdmn.eq.1) then
        boundary_index(5)       = bc_cond_sym
     else  
        boundary_index(5)       = bc_cond_sym
     endif
     boundary_index(6)          = bc_cond_sym
  end if

end subroutine app_domain_setup

