
! Tim Yang
! 22/09/2015
! This subrountine uses amr_initial_soln as a template,
! to compute the norm(s) of solution, when the true solution is
! avaliable.
#include "paramesh_preprocessor.fh"
  subroutine tim_true_solution(pid,noprocs,f_time, max_domain, min_domain)
!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
  use paramesh_interfaces
  use paramesh_dimensions
  use tree
  use time_dep_parameters
  use multigrid_parameters
  use generic_parameters
  use workspace
  use physicaldata
  use checkpoint_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  double precision, intent(in) :: f_time, max_domain, min_domain
  integer :: nguard0, lb,k,j,i
  integer :: iopt, nlayers, ierr
  double precision :: dx, xi, yi, pi, pikk=1.0
  double precision :: positive_dist, max_positive_dist
  double precision :: negative_dist, min_negative_dist
  double precision :: true_radius
  double precision :: G_max_positive_dist, G_min_negative_dist
  double precision :: negative_dist2, min_negative_dist2, G_min_negative_dist2
  double precision :: energy
  double precision :: G_energy
  character(len=80) :: file_write, file_write_energy

!--------------------------------------------------------------
  nguard0 = nguard*npgs
  iopt = 1
  pi=4.0*datan(pikk)

  positive_dist = 0.0
  max_positive_dist = 0.0
  negative_dist = 0.0
  min_negative_dist = 0.0
  negative_dist2 = 0.0
  min_negative_dist2 = 100000.0
  do lb=1, lnblocks
    dx = bsize(1,lb)/real(nxb)
    if(lrefine(lb).eq.lrefine_max)then
      do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
          do i=il_bnd+nguard,iu_bnd-nguard
            xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
            yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)

            if(unk(1,i,j,k,lb).ge.0.0) then
              positive_dist = sqrt( (xi)*(xi) + (yi)*(yi) )
              if(positive_dist.ge.max_positive_dist) then
                max_positive_dist = positive_dist
              endif
            endif

            if(unk(1,i,j,k,lb).ge.0.99) then
              negative_dist= sqrt( (xi)*(xi) + (yi)*(yi) )
              if(negative_dist.ge.min_negative_dist) then
                min_negative_dist = negative_dist
              endif
            endif

            if(unk(1,i,j,k,lb).le.-0.99) then
              negative_dist2= sqrt( (xi)*(xi) + (yi)*(yi) )
              if(negative_dist2.le.min_negative_dist2) then
                min_negative_dist2 = negative_dist2
              endif
            endif

          enddo
        enddo
      enddo
    endif
  enddo
  call MPI_Reduce(max_positive_dist, G_max_positive_dist, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(min_negative_dist, G_min_negative_dist, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(min_negative_dist2, G_min_negative_dist2, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

  if(pid.eq.0) then
21668 format(A20)
21655 format(f21.16,',',f21.16)
    write(file_write,21668) 'width_error_ep_1.txt'
    open(unit=12315, file=file_write, position="append", action="write")
    true_radius = sqrt( 2.0*( (1.0/2.0) - f_time ) )
    write(*,*) "dist", abs(G_min_negative_dist-G_min_negative_dist2)
    write(*,*) "error", abs(true_radius-G_max_positive_dist)
    write(12315,21655) abs(G_min_negative_dist-G_min_negative_dist2), abs(true_radius-G_max_positive_dist)
    close(12315)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  energy = 0.0
  G_energy = 0.0
  do lb=1, lnblocks
    dx = bsize(1,lb)/real(nxb)
    if(lrefine(lb).eq.lrefine_max)then
      do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
          do i=il_bnd+nguard,iu_bnd-nguard
            energy = energy + dx*dx*(ep*&
                     ( ((unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb))/(2.0*dx))*((unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb))/(2.0*dx))&
                      +((unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb))/(2.0*dx))*((unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb))/(2.0*dx)))&
                      +((unk(1,i,j,k,lb)*unk(1,i,j,k,lb)-1.0)*(unk(1,i,j,k,lb)*unk(1,i,j,k,lb)-1.0))/(4.0*ep))
          enddo
        enddo
      enddo
    endif
  enddo
  call MPI_Reduce(energy, G_energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if(pid.eq.0) then
31668 format(A25)
31655 format(f30.16)
    write(file_write_energy,31668) 'interface_energy_ep_1.txt'
    open(unit=32315, file=file_write_energy, position="append", action="write")
    write(32315,31655) G_energy
    close(32315)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine tim_true_solution
