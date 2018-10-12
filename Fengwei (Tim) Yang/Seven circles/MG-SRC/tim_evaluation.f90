
! Tim Yang
! 22/09/2015
! This subrountine uses amr_initial_soln as a template,
! to compute the norm(s) of solution, when the true solution is
! avaliable.
#include "paramesh_preprocessor.fh"
  subroutine tim_evaluation(pid,noprocs,f_time, time_step)
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
  double precision, intent(in) :: f_time, time_step
  integer :: nguard0, lb,k,j,i
  integer :: iopt, nlayers, ierr
  double precision :: dx, xi, yi, pi, results, pikk=1.0
  double precision :: left, right, up, down
  double precision :: Gleft, Gright, Gup, Gdown
  integer :: counting = 0

!--------------------------------------------------------------
  nguard0 = nguard*npgs
  iopt = 1
  pi=4.0*datan(pikk)

  counting = 0
  do lb=1, lnblocks
    dx = bsize(1,lb)/real(nxb)
    if(lrefine(lb).eq.lrefine_max)then
      do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
          do i=il_bnd+nguard,iu_bnd-nguard
            xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
            yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)

            if(abs(xi-(2.0*pi/4.0)).lt.dx.and.abs(yi-(2.0*pi/4.0)).lt.dx)then
              if(xi.lt.(2.0*pi/4.0).and.yi.lt.(2.0*pi/4.0))then
                left = unk(1,i,j,k,lb)
              elseif((2.0*pi/4.0).lt.xi.and.(2.0*pi/4.0).lt.yi)then
                right = unk(1,i,j,k,lb)
              elseif(xi.lt.(2.0*pi/4.0).and.(2.0*pi/4.0).lt.yi)then
                up = unk(1,i,j,k,lb)
              elseif((2.0*pi/4.0).lt.xi.and.yi.lt.(2.0*pi/4.0))then
                down = unk(1,i,j,k,lb)
              endif
              counting = counting + 1
            endif

          enddo
        enddo
      enddo
    endif
  enddo
  call MPI_Reduce(left, Gleft, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(right, Gright, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(up, Gup, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(down, Gdown, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  results = ((Gleft + Gright + Gup + Gdown) / 4.0)
  if(pid.eq.0) write(*,*) '-------at 2PI/4-------', results
  if(results.lt.0.0.and.switch1.eq.0)then
    if(pid.eq.0) write(*,*) 'phase 1 time:', f_time
    call tim_visual_paraview(time_step)
    switch1 = 1
  endif

  counting = 0
  do lb=1, lnblocks
    dx = bsize(1,lb)/real(nxb)
    if(lrefine(lb).eq.lrefine_max)then
      do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
          do i=il_bnd+nguard,iu_bnd-nguard
            xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
            yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)

            if(abs(xi-((6.0*pi)/4.0)).lt.dx.and.abs(yi-(6.0*pi/4.0)).lt.dx)then
              if(xi.lt.((6.0*pi)/4.0))then
                left = unk(1,i,j,k,lb)
              elseif(((6.0*pi)/4.0).lt.xi)then
                right = unk(1,i,j,k,lb)
              elseif(xi.lt.(6.0*pi/4.0).and.(6.0*pi/4.0).lt.yi)then
                up = unk(1,i,j,k,lb)
              elseif((6.0*pi/4.0).lt.xi.and.yi.lt.(6.0*pi/4.0))then
                down = unk(1,i,j,k,lb)
              endif
              counting = counting + 1
            endif

          enddo
        enddo
      enddo
    endif
  enddo
  call MPI_Reduce(left, Gleft, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(right, Gright, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(up, Gup, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(down, Gdown, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  results = ((Gleft + Gright) / 2.0)
  if(pid.eq.0) write(*,*) '-------at 3PI/4-------', results
  if(results.lt.0.0)then
    if(pid.eq.0) write(*,*) 'phase 2 time:', f_time
    call tim_visual_paraview(time_step)
    stop
  endif

  end subroutine tim_evaluation
