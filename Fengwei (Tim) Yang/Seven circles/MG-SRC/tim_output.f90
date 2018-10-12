
! Tim Yang
! 01/02/2012
! This subrountine uses amr_initial_soln as a template,
! to compute the norm(s) of solution, when the true solution is
! avaliable.
#include "paramesh_preprocessor.fh"
      subroutine tim_output(lvl)
!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr
      use paramesh_dimensions
      use physicaldata
      use tree
      implicit none
      include 'mpif.h'
      integer :: nguard0, lb,k,j,i
!--------------------------------------------------------------
      double precision :: xi, yi, dx
      integer :: lvl, counting, base_case
      character(len=80) :: filename_one
      character(len=80) :: filename_two
      character(len=80) :: filename_mu1
      character(len=80) :: filename_mu2
      character(len=80) :: filename_p1
      character(len=80) :: filename_p2
      nguard0 = nguard*npgs
      counting = 0

      base_case = 3

      if(lrefine_max.gt.base_case)write(filename_one,*) lvl*10+1,"phi"
      if(lrefine_max.gt.base_case)write(filename_mu1,*) lvl*10+1,"mu"
      if(lrefine_max.gt.base_case)write(filename_p1,*) lvl*10+1,"p"

      write(filename_two,*) lvl*10+2,"phi"
      write(filename_mu2,*) lvl*10+2,"mu"
      write(filename_p2,*) lvl*10+2,"p"

      if(lrefine_max.gt.base_case)open(unit = 65776, file=filename_one)
      if(lrefine_max.gt.base_case)open(unit = 86645, file=filename_mu1)
      if(lrefine_max.gt.base_case)open(unit = 78887, file=filename_p1)

      open(unit = 65777, file=filename_two)
      open(unit = 86646, file=filename_mu2)
      open(unit = 78888, file=filename_p2)

666   format(3(2x,f25.16))

      do lb=1, lnblocks
        if(lrefine(lb).eq.lrefine_max)then
          do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
            do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
              do i=il_bnd+nguard,iu_bnd-nguard
                counting = counting + 1
              enddo
            enddo
          enddo
        endif
      enddo
      write(65777,*) counting
      write(86646,*) counting
      write(78888,*) counting
      do lb=1, lnblocks
        dx = bsize(1,lb)/real(nxb)
        if(lrefine(lb).eq.lrefine_max)then
          do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
            do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
              do i=il_bnd+nguard,iu_bnd-nguard
                xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
                yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
                write(65777,666) xi, yi, unk(1,i,j,k,lb)
                write(86646,666) xi, yi, unk(6,i,j,k,lb)
                write(78888,666) xi, yi, unk(11,i,j,k,lb)
              enddo
            enddo
          enddo
        endif
      enddo

      if(lrefine_max.gt.base_case)then
        do lb=1,lnblocks
          refine(lb) = .false.
          derefine(lb) = .true.
        enddo
        call amr_refine_derefine

        counting = 0
        do lb=1, lnblocks
          if(lrefine(lb).eq.lrefine_max-1)then
            do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                do i=il_bnd+nguard,iu_bnd-nguard
                  counting = counting + 1
                enddo
              enddo
            enddo
          endif
        enddo
        write(65776,*) counting
        write(86645,*) counting
        write(78887,*) counting
        do lb=1, lnblocks
          dx = bsize(1,lb)/real(nxb)
          if(lrefine(lb).eq.lrefine_max-1)then
            do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                do i=il_bnd+nguard,iu_bnd-nguard
                  xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
                  yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
                  write(65776,666) xi, yi, unk(1,i,j,k,lb)
                  write(86645,666) xi, yi, unk(6,i,j,k,lb)
                  write(78887,666) xi, yi, unk(11,i,j,k,lb)
                enddo
              enddo
            enddo
          endif
        enddo
      endif

      if(lrefine_max.gt.base_case)close(65776)
      if(lrefine_max.gt.base_case)close(78887)
      if(lrefine_max.gt.base_case)close(86645)
      close(65777)
      close(86646)
      close(78888)
      write(*,*) "-------------------------------------------------"
      write(*,*)  "outputs for convergence test done"
      write(*,*) "-------------------------------------------------"
      stop
      end subroutine tim_output
