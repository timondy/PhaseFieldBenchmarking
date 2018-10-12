!#include "paramesh_preprocessor.fh"

subroutine shampine_memory_alloc
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  use solution_parameters
  implicit none

  allocate(previt(il_bnd:iu_bnd, jl_bnd:ju_bnd, kl_bnd:ku_bnd, maxblocks, total_vars))

  shampine_tol = 3.0E-10

end subroutine shampine_memory_alloc

subroutine shampine_memory_dealloc
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use tree
  use solution_parameters
  implicit none

  deallocate(previt)

end subroutine shampine_memory_dealloc


subroutine shampine_convergence(m, pass, mype)
  use paramesh_dimensions
  use time_dep_parameters
  use physicaldata
  use multigrid_parameters
  use paramesh_mpi_interfaces
  use solution_parameters
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(inout) :: m, pass, mype
  double precision :: val, sigma, calc, gval, gcalc
  integer :: i, j, k, lb, v, gblks, ierr

  v = 1

  pass = 0
  val = 0.0
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

  call MPI_AllReduce(val, gval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_AllReduce(lnblocks, gblks, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  val = sqrt(gval/gblks)

  if (m.eq.1) then

     shampine_beta = val

  else

     sigma = (val/shampine_beta)**(1.0/m)

     calc = sigma/(1.0-sigma) * val

     call MPI_AllReduce(calc, gcalc, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

     if (mype.eq.0) print*,'Shampine test:', calc, 0.33*shampine_tol
     if (calc.lt.0.33*shampine_tol) pass = 1

  endif

! Finally store current solution over previous iterations 
  call shampine_store()

  

end subroutine shampine_convergence

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

  do lb=1,lnblocks
     do k=nguard*k3d+kl_bnd, ku_bnd-nguard*k3d
        do j=nguard+1, nyb+nguard
           do i=nguard+1, nxb+nguard
              previt(i, j, k, lb, v) = unk(v, i, j,k, lb)              
           end do
        end do
     end do
  end do

end subroutine shampine_store
