!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! shampine_memory_alloc
!!!! * Allocates previt array in Paramesh structure
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! shampine_memory_dealloc
!!!! * Deallocates previt array in Paramesh structure
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! shampine_convergence
!!!! * Assesses convergence on this step
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!! shampine_store
!!!! * Stores the solution into previt
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shampine_memory_alloc
  use paramesh_dimensions
  use time_dep_parameters
  use generic_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none

  if (verbose.ge.3) print*,'previt allocated in shampine_memory_alloc'

  allocate(previt(il_bnd:iu_bnd, jl_bnd:ju_bnd, kl_bnd:ku_bnd, maxblocks, total_vars))

  ! Now set in MG-SRC/amr_modules.f90
  !shampine_tol = 3.0E-10

end subroutine shampine_memory_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shampine_memory_dealloc
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none

  deallocate(previt)

end subroutine shampine_memory_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine shampine_convergence(m, pass, mype)
  use paramesh_dimensions
  use time_dep_parameters
  use generic_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(inout) :: m, pass, mype
  double precision :: val, sigma, calc, gval, gcalc
  integer :: i, j, k, lb, v, gblks, ierr

  v = 1

  pass = 0
  val = 0.0

!!!$OMP PARALLEL SHARED(val) PRIVATE(lb, i, j, k)
!!!$OMP DO SCHEDULE(DYNAMIC,10) REDUCTION(+:val)
  do lb=1,lnblocks
    if(nodetype(lb).eq.1) then
       do k=nguard*k3d+kl_bnd, ku_bnd-nguard*k3d
           do j=nguard+1, nyb+nguard
              do i=nguard+1, nxb+nguard
                 val = val + (unk(v, i, j,k, lb) - previt(i, j, k, lb, 1))*(unk(v, i, j,k, lb) - previt(i, j, k, lb, 1))
              end do
           end do
        end do
     endif
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

  call MPI_AllReduce(val, gval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_AllReduce(lnblocks, gblks, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  val = sqrt(gval/gblks)

  if (m.eq.1) then

     shampine_beta = val

  else

     sigma = (val/shampine_beta)**(1.0/m)

     calc = sigma/(1.0-sigma) * val

     call MPI_AllReduce(calc, gcalc, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

     if (calc.lt.0.0) then
        if (mype.eq.0.and.verbose.ge.2) print*,'Shampine: Uurgh', sigma, val, shampine_beta

        calc = -calc
     endif

     if (mype.eq.0.and.verbose.ge.2) print*,'Shampine test:', calc, 0.33*shampine_tol
     if (calc.lt.0.33*shampine_tol) pass = 1

  endif

! Finally store current solution over previous iterations 
  call shampine_store()

  

end subroutine shampine_convergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shampine_store
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  implicit none
  integer :: i, j, k, lb, v

  v = 1

!!!$OMP PARALLEL SHARED(lnblocks,v) PRIVATE(lb,k,j,i)
!!!$OMP DO SCHEDULE(DYNAMIC,8) 
  do lb=1,lnblocks
!     print *,lb
     do k=nguard*k3d+kl_bnd, ku_bnd-nguard*k3d
        do j=nguard+1, nyb+nguard
           do i=nguard+1, nxb+nguard
              previt(i, j, k, lb, v) = unk(v, i, j, k, lb)              
           end do
        end do
     end do
  end do
!!!$OMP END DO NOWAIT
!!!$OMP END PARALLEL 

end subroutine shampine_store
