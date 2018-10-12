!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine conserve_check_1blk(mype,istep,ivar)


!-----------------------------------------------------------------------
! This routine sums the solution to test whether fluxes are
! properly applied to ensure conservation.
!
! Written:      Peter MacNeice          July 1999
!
!
      use physicaldata
      use tree

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"
        include 'shmem_reduce.fh'


!------------------------------------------

! local variables
      integer remote_blk,remote_pe,ich
      real rho_total_lb(maxblocks),sum,total
      real rho_total,rho_total_local
      save rho_total,rho_total_local,rho_total_lb,sum

!-----------------------------------------------------------------------


      rho_total = 0.
      rho_total_local = 0.

! loop over leaf blocks on this processor
      if(lnblocks.gt.0) then
      do lb=1,lnblocks

        dx = size(1,lb)/(real(nxb/2)*2)
        dy = 1.
        dz = 1.
        if(ndim.ge.2) dy = size(2,lb)/(real(nyb/2)*2)
        if(ndim.eq.3) dz = size(3,lb)/(real(nzb/2)*2)

        rho_total_lb(lb) = 0.
        do k = 1+nguard*npgs*k3d,nzb+nguard*npgs*k3d
        do j = 1+nguard*npgs*k2d,nyb+nguard*npgs*k2d
        do i = 1+nguard*npgs,nxb+nguard*npgs
        rho_total_lb(lb) =  & 
     &            rho_total_lb(lb)+unk(ivar,i,j,k,lb)*dx*dy*dz
        enddo
        enddo
        enddo

      enddo
      endif


      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then
        rho_total_local = rho_total_local + rho_total_lb(lb)

      endif
      enddo
      endif

      call shmem_barrier_all()



!------------------------------------------

        lb = 7
      if(nodetype(lb).eq.2) then

        total = 0.
        do ich = 1,nchild
        remote_blk = child(1,ich,lb)
        remote_pe  = child(2,ich,lb)
        call shmem_real_get( sum,rho_total_lb(remote_blk), & 
     &                   1,remote_pe)
        total = total + sum
        enddo
        write(*,*) 'Comparison of parent and sum of children ', & 
     &       lb,rho_total_lb(lb),total,' ivar ',ivar

      endif

      call shmem_barrier_all()

      return
      end
