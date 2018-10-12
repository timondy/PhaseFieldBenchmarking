!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_timestep(dt,dtmin,dtmax,mype)


!-----------------------------------------------------------------------
! This routine computes the hydro timestep, and distributes
! it to all the processors.
!
! Written:      Peter MacNeice February 1997
!
!
      use physicaldata
      use tree

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"
      include 'shmem_reduce.fh'


      real  kappa
      parameter( courant=.1, kappa=1.0 )

!------------------------------------------

! local variables
      real tpwrk(max(maxlevels/2+1,SHMEM_REDUCE_MIN_WRKDATA_SIZE))

      real dtl,dtmaxl
      save dtl,tpwrk,dtmaxl

      integer npes,shmem_n_pes

!-----------------------------------------------------------------------


      call shmem_barrier_all()

      npes = shmem_n_pes()

      dtmax = 0.
      dtmin = 1.e30
      dtl = 1.e30
      dtlevel(:) = 1.e30

!------------------------------------------


! loop over leaf blocks on this processor
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif

        dx = size(1,lb)/(real(nxb/2)*2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! users timestep calaculation

        dtl = courant*dx*dx/kappa

! end of users timestep calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        dtlevel(lrefine(lb)) = min(dtlevel(lrefine(lb)),dtl)



#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo

      endif
      call shmem_barrier_all()



!------------------------------------------


! find smallest timesteps for each refinement level  across all processors
      do i=1,maxlevels
        call comm_real_min_to_all(dtlevel(i),dtlevel(i))
      enddo

      call shmem_barrier_all()

! ensure that the timestep does not increase as refinement level increases.
      do i=2,maxlevels
        dtlevel(i) = min(dtlevel(i),dtlevel(i-1))
      enddo

! ensure that each timestep is either equal to or a factor 2 larger
! than the timestep at the next higher refinement level.
      dtmin = dtlevel(maxlevels)
      do i = maxlevels-1,1,-1
        ratio  = dtlevel(i)/dtmin
        iratio = max((int(ratio)/2)*2,1)
        iratio = min(iratio,2)
        dtlevel(i) = real(iratio)*dtmin
      enddo

      dtmaxl = 0.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
        dtmaxl = max(dtmaxl,dtlevel(lrefine(lb)))
#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo
      endif
      call shmem_barrier_all()

! find largest timestep for any leaf node across all processors
      call comm_real_max_to_all(dtmaxl,dtmaxl)
      call shmem_barrier_all()

      dtmax = dtmaxl


#ifdef VAR_DT
      dt = dtmax
#else
      dt = dtmin
      dtmax = dtmin
#endif

      if(mype.eq.0) write(*,*) 'proc ',mype,' dt ',dt
      call shmem_barrier_all()

      return
      end
