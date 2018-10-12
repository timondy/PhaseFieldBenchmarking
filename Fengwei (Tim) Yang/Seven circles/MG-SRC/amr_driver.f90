!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                          C A M P F I R E 
!!!! 
!!!!             Campfire: Adaptive, Multilevel, Parallel, 
!!!!                Fully Implicit Research Environment.
!!!! 
!!!! This forms a two part system that complements the Paramesh adaptive
!!!! meshing software.  The two parts are the locally adaptive FAS multigrid 
!!!! code contained in the libMG library, and application specific routines
!!!! that are called from here to perform local smoothing of the solution.
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!! Top level routine
!!!! * Calls set up functions for Paramesh
!!!! * Calls application specific functions to govern case
!!!!       -> app_parameter_set
!!!!       -> app_domain_setup
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!! set_domain_limits
!!!! * Assigns maximum and minimum extents of domain
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 
!!!! input_error_checking
!!!! * Checks input parameter file (amr_runtime_parameters) for consistency
!!!! * Terminates if badness is found
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! History: 
!  Control file for paramesh time dependent non-linear multigrid solver
!  Written by J. Green
!  Based upon: http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr_users_guide.html#main_p
!
!
program campfire

  ! include file to define physical qualities of the model and mesh
  use paramesh_dimensions
  use physicaldata
  ! include workspace so interp_mask_work can be set
  use workspace
  ! include file defining the tree
  use tree
  use multigrid_parameters
  use generic_parameters
  ! IMPORTANT: Expose the PARAMESH interface blocks to your code
  use paramesh_interfaces
  use paramesh_mpi_interfaces
  use time_dep_parameters
  
  ! Set implicit to none, not needed but is sensible
  implicit none

  ! include file required for mpi library.
  include 'mpif.h'
  
  integer pid, noprocs, ierr
  integer :: tid, nthreads, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

!  print *,'Into code. Calling amr_initialize' 
! start off by setting up Paramesh
  call amr_initialize

  total_vars=nvar/nunkvbles

! Call application specific inital options 
  call app_parameter_set()

  dtold=dt

!#ifdef MPE
!  call MPE_Init_log()
!#endif
  ! Set values of noprocs and pid with appropriate mpi calls
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, noprocs, ierr)

  if (pid.eq.0) print *,'Number of MPI processes = ', noprocs
!$OMP PARALLEL PRIVATE(tid,nthreads)
      tid = OMP_GET_THREAD_NUM
!     Only master thread does this
      IF (tid.eq.0.and.pid.eq.0) THEN
        nthreads = OMP_GET_NUM_THREADS
        print *, 'Number of OpenMP threads = ', nthreads
      END IF
!$OMP END PARALLEL
  
  ! Check data specified via amr_runtime_parameters is sensible
  call input_error_checking(pid, noprocs)

  if (pid.eq.0.and.verbose.ge.3) print *,'Domain set-up'
  call app_domain_setup(pid, noprocs)

  ! Get the parameters from the command line arguments, then either load checkpoints or do initial refinements
  call get_input_parameters(pid, noprocs)

  if (shampine_test.eq.1) call shampine_memory_alloc
  interp_mask_unk(:) = 0
  interp_mask_work(:) = 0

  ! Set up variables for which work arrays to iterate over in restriction and prolongation
  wvar_rest_loop_begin = 1
  wvar_rest_loop_end = 2
  wvar_prol_loop_begin = 1
  wvar_prol_loop_end = 1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Call main MG loop routine
  if (pid.eq.0.and.verbose.ge.3) print *,'Off to the transient code...'

  call amr_mg_control(pid, noprocs)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Clean up memory
!!-fwy  if (shampine_test.eq.1) call shampine_memory_dealloc
  call amr_mg_end()
!!-fwy  call amr_close()

  stop
end program campfire

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Extra rountine to process output from mpi_amr_global_domain_limits function 
! into boundary_box array

subroutine set_domain_limits
  ! include file to define physical qualities of the model and mesh
  use paramesh_dimensions
  use physicaldata
  ! include workspace so interp_mask_work can be set
  use workspace
  ! include file defining the tree
  use tree
  use multigrid_parameters

  ! Set implicit to none
  implicit none
  
  ! x boundaries (boundary boxes 1 and 2)
  boundary_box(1,2:3,1:2) = -1.e30 ! effectively this is -infinity
  boundary_box(2,2:3,1:2) =  1.e30
  boundary_box(1,1,1)     = -1.e30
  boundary_box(2,1,1)     =  grid_xmin ! min dimension of the comp. domain
  boundary_box(1,1,2)     =  grid_xmax ! max     "     "   "   "      "
  boundary_box(2,1,2)     =  1.e30
  
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
  end if
     
  ! z boundaries
  if (ndim .eq. 3) then
     ! five and six used since in 1d and 2d codes boundary_box only goes to 4
     boundary_box(1,1:2,5:6) = -1.e30
     boundary_box(2,1:2,5:6) =  1.e30
     boundary_box(1,3,5)       = -1.e30
     boundary_box(2,3,5)       = grid_zmin
     boundary_box(1,3,6)        = grid_zmax
     boundary_box(2,3,6)        = 1.e30
  end if

end subroutine set_domain_limits

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Consistency checking of amr_runtime_parameters

subroutine input_error_checking(pid, noprocs)
  ! include file to define physical qualities of the model and mesh
  use paramesh_dimensions
  use physicaldata
  ! include workspace so interp_mask_work can be set
  use workspace
  ! include file defining the tree
  use tree
  use multigrid_parameters

  ! Set implicit to none
  implicit none
  integer, intent(in) :: pid, noprocs
  
  if (ndim.eq.2) then
     if (nboundaries.ne.4) then
        print *, 'ERROR : In amr_runtime_parameters'
        print *, 'Running with ndim=2 so nboundaries should be 4'
        stop
     endif

!!     if (nxb.ne.nyb) then
!!        print *, 'ERROR : In amr_runtime_parameters'
!!        print *, 'nxb should be the same as nyb'
!!        stop
!!     endif

     if (MOD(nxb,2).eq.1) then
        print *, 'ERROR : In amr_runtime_parameters'
        print *, 'nxb should be even'
        stop
     endif

  else if (ndim.eq.3) then

     if (nboundaries.ne.6) then
        print *, 'ERROR : In amr_runtime_parameters'
        print *, 'Running with ndim=3 so nboundaries should be 6'
        stop
     endif

     if (nxb.ne.nyb) then
        print *, 'ERROR : In amr_runtime_parameters'
        print *, 'nxb should be the same as nyb'
        stop
     endif

     if (nxb.ne.nzb) then
        print *, 'ERROR : In amr_runtime_parameters'
        print *, 'nxb should be the same as nzb'
        stop
     endif

     if (MOD(nxb,2).eq.1) then
        print *, 'ERROR : In amr_runtime_parameters'
        print *, 'nxb should be even'
        stop
     endif

  else
!!     print *, 'ERROR : In amr_runtime_parameters'
!!     print *, 'ndim should be 2 or 3'
!!     stop
  endif

  if (MOD(nvar,nunkvbles).ne.0) then
     print *, 'ERROR : In amr_runtime_parameters'
     write (6,'(A,I2,A,I1,A)') 'nvar (',nvar,') should be a product of nunkvbles (',nunkvbles,')'
     stop
  endif

!  if (nunkvbles-1.gt.nvar_work) then
!     print *, 'ERROR : In amr_runtime_parameters'
!     write (6,'(A,I2,A,I1,A)') 'nvar_work (',nvar_work,') should be at least nunkvbles-1 (',nunkvbles-1,')'
!     stop
!  endif



end subroutine input_error_checking
