!#include "paramesh_preprocessor.fh"

module multigrid_parameters
  ! module that should control everything to do with the multigrid
  ! note this module was the first and will contain anything which hasn't 
  ! subsequently found a more appropriate module
  ! Weighting, probably optimal at 2/3 for linear elliptic, .7-.9 for most others
  real, parameter :: weight=0.885!2.0/3.0
  ! Selector whether Jacobi or Gauss-Seidel smooth is used, note not yet implemented
!  logical, parameter :: gauss_seidel=.false.
!  logical, parameter :: gauss_seidel=.true.
  ! How many pre/post smooths
  integer, parameter :: smooth_count_max=4
  ! Lowest (coarsest) refinement level which v-cycles go to
!  integer, parameter :: mg_min_lvl=4
  integer :: mg_min_lvl
  ! Highest (finest) ref level, set by amr_mg_get_max_refinement
  integer :: mg_max_level
  ! Whether to do a full solve on mg_min_lvl
  logical, parameter :: full_solve=.true.
  ! How many smooths constitute a full solve
  integer, parameter :: solve_count_max=8
  ! How many v-cycles to do before giving up and doing something else
  integer, parameter :: max_v_cycle_count=40
  ! Number of cycles leading up to the full multigrid, keep small (<4)
  integer, parameter :: pre_cycle_count=1
  ! tolerance value which the max abs of the defect will be tested against
!  real, parameter :: defect_tol=1e-10
  real, parameter :: defect_tol=1e-9
  ! Not used yet, might be needed for mlat scheme, not sure yet.
  integer, allocatable, dimension (:) :: nodetype_copy                    ! CEG : seems to be used
  logical, allocatable, dimension (:) :: has_children                     ! CEG : used to limit operations
  ! couple of variable which are used by the code, lots
  integer :: current_var, total_vars, current_unk, current_work
  ! Whether or not to do a full multigrid solve or just lots of v cycles (not really working)
  logical, parameter :: full_multigrid_cycle=.false.
  ! delta x on finest grid
  real :: min_dx
  ! How many variables per variable
  integer, parameter :: nunkvbles=5
  !Boundary conditions
  integer, parameter :: bc_cond_far = -21    ! Far field
  integer, parameter :: bc_cond_sym = -22    ! Symmetry
end module multigrid_parameters

module time_dep_parameters
  ! Module for time depenant stuff
  ! Mode controls the time advance code, leave at 1 for now
  integer :: mode_1=1
  ! current value for dt, not a parameter for variable dt later
  double precision :: dt=1.e-6
  ! old value for dt
  double precision :: dtold=1.e-6
  ! Max possible value of dt, not implemented yet, mostly for variable dt
  double precision :: max_dt
  ! total simulation time (end time-start time, not just end time)
  ! Note start time is set in solution parameters
!  double precision, parameter :: simulation_time=0.40
!  double precision, parameter :: simulation_time=0.5
!  double precision, parameter :: simulation_time=2.30
!  double precision, parameter :: simulation_time=10.0!200
!  double precision, parameter :: simulation_time=20000.0
  double precision, parameter :: simulation_time=1000.0
  ! current value of time
  double precision :: time
  ! max number of dt steps to take, leave this as a biiiiiiiiig number
  integer, parameter :: max_time_iter=300000000
end module time_dep_parameters

module solution_parameters
  ! Module for storing all the parameters for the problem being solved.
  ! Note this is full of old/dud parameters which are no longer used.
  ! End user should keep track of which are and aren't used
  ! The only critical ones are 
  ! t0: start time. For most cases this is 0 but if you need to use a non-zero 
  ! start time then look no further. This is here in case start time depends on
  ! other solution parameters.
  ! nucleate_x/y/z are the co-ordinates of the nucleation point
  ! nuc_radius is the nucleation radius
  ! anti_trapping_mod is a multiplier for the antitrapping terms, normally 1.0
  ! solute sets whether code is phase/thermal (.false.) or phase/solute (.true.) IFF only 2 real vars
  ! For thermal/phase le still NOT used to set D_therm
  double precision, parameter :: a_2=0.626663849
  double precision, parameter :: a_1=5.0*sqrt(2.0)/8.0
  double precision :: t0=0.0
!  double precision :: D_therm,A_0,lambda=5.0,delta=-0.2,epsilon=0.02,epsilon_tilde
!CEG run 24/2/10 modified epsilon up by a lot, delta from 0.2 to 0.6 and Lewis number up from 100 to 500 for 8017
!  double precision :: D_therm,A_0,lambda=5.0,delta=-0.6,epsilon=0.02,epsilon_tilde
!CEG run 25/2/10 to reproduce James's case 8075
!  double precision :: D_therm,A_0,lambda=5.0,delta=-0.6,epsilon=0.02,epsilon_tilde
!CEG run 25/2/10 to reproduce James's case 8075
!  double precision :: D_therm,A_0,lambda=3.2,delta=-0.65,epsilon=0.05,epsilon_tilde
!  double precision :: D_solute,ke=0.3,le=100,mcinf=0.05
!  double precision :: D_solute,ke=0.3,le=100,mcinf=0.05

  integer :: eighthdmn

!  double precision :: D_therm,A_0,lambda=2.0,delta=-.6,epsilon=0.02,epsilon_tilde
!  double precision :: D_solute,ke=.15,le=40,mcinf=0.1


! New test case run
  double precision :: D_therm, A_0, lambda, delta, epsilon, epsilon_tilde
  double precision :: D_solute, ke, le, mcinf

  double precision :: nucleate_x, nucleate_y, nucleate_z, nuc_radius
  double precision :: anti_trapping_mod=1.000
!  double precision :: anti_trapping_mod=0.000
!  logical, parameter :: solute=.false.
  logical :: solute

  integer :: grow

  integer, parameter :: shampine_test=1
  double precision :: shampine_beta, shampine_tol

contains
!  subroutine get_h_psi(psi,h)
!    implicit none
!    double precision, intent(in) :: psi
!    double precision, intent(out) :: h
!    h=psi
!  end subroutine get_h_psi
  subroutine get_dh_dpsi(psi,dh)
    implicit none
    double precision, intent(in) :: psi
    double precision, intent(out) :: dh
    dh=1.0
  end subroutine get_dh_dpsi
!  subroutine get_d2h_dpsi2(psi,d2h)
!    implicit none
!    double precision, intent(in) :: psi
!    double precision, intent(out) :: d2h
!    d2h=0.0
!  end subroutine get_d2h_dpsi2
end module solution_parameters

