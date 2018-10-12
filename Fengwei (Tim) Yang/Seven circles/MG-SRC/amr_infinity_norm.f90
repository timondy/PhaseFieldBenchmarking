
! Tim Yang
! 01/02/2012
! This subrountine uses amr_initial_soln as a template,
! to compute the norm(s) of solution, when the true solution is
! avaliable.
#include "paramesh_preprocessor.fh"
      subroutine amr_infinity_norm(t,lvl,mype)
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
      double precision :: Pi, pikk = 1.0
      REAL :: xi, yi, zi, maximum, diff, t, maxdiff, tmp
      real :: restricted(1:lnblocks*nxb*nyb)
      real :: array(1:lnblocks*nxb*nyb)
      integer :: lvl, count_node_1, count_node_2, mype, count_node_3
      integer :: count_node_4
      character(len=80) :: file_write, file_read
      nguard0 = nguard*npgs
      Pi=4.0*datan(pikk)
      maximum = 0.0
      diff = 0.0
      count_node_1 = 0
      count_node_2 = 1
      count_node_3 = 0
      count_node_4 = 1

      if(mype.eq.0)write(file_write,*) lrefine_max
      if(mype.eq.0)write(file_read,*)  lrefine_max-1
      if(mype.eq.0) open(unit = 10001, file=file_write)
      if(mype.eq.0) open(unit = 10002, file=file_read)

      if(lnblocks.gt.0) then
        do lb=1,lnblocks
          !! nodetype = 1 when on finest grid.
          if(nodetype(lb).eq.1) then
            do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                do i=il_bnd+nguard,iu_bnd-nguard
                  if(mype.eq.0)write(10001,*) unk(1,i,j,k,lb)
                enddo ! do i
              enddo ! do j
            enddo ! do k
          endif ! if(nodetype(lb).eq.1 .or. advance_all_levels) then
        enddo ! do lb=1, lnblocks
      endif ! if(lnblocks.gt.0) then

      if(lrefine_max.gt.2.and.mype.eq.0)then
        do lb=1,lnblocks
          refine(lb) = .false.
          derefine(lb) = .true.
        enddo
        call amr_refine_derefine
!!        do lb=1, lnblocks
!!          if(nodetype(lb).eq.1) then
!!            tmp = 0.0
!!            count_node_3 = 0
!!            do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
!!              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
!!                do i=il_bnd+nguard,iu_bnd-nguard
!!                  if(count_node_3.lt.16)then
!!                    tmp = tmp + unk(1,i,j,k,lb)
!!                    count_node_3 = count_node_3 + 1
!!                  else if (count_node_3.eq.16)then
!!                    array(count_node_4) = tmp / 16.0
!!                    count_node_4 = count_node_4 + 1
!!                  endif
!!                enddo
!!              enddo
!!            enddo
!!          endif
!!        enddo

        do lb=1, lnblocks
          if(nodetype(lb).eq.1)then
            count_node_1 = count_node_1 + nxb*nyb 
          endif
        enddo
  
        do lb=1, count_node_1
          read(10002,*) restricted(lb)
        enddo
!!        do lb=1, count_node_4
!!          diff = abs(array(lb) - restricted(lb))
!!          if(diff.gt.maxdiff)then
!!            maxdiff=diff
!!          endif
!!        enddo
        do lb=1, lnblocks
          if(nodetype(lb).eq.1)then
            do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                do i=il_bnd+nguard,iu_bnd-nguard
                  diff = abs(unk(1,i,j,k,lb) - restricted(count_node_2))
                  count_node_2 = count_node_2 + 1
                  if(diff.gt.maxdiff)then
                    maxdiff = diff
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
        if(mype.eq.0)write(*,*) "infinity norm: ", maxdiff
      else
        if(mype.eq.0)write(*,*) "Coarsest level, no output."
      endif

      if(mype.eq.0) close(10001)
      if(mype.eq.0) close(10002)
      return
      end subroutine amr_infinity_norm
