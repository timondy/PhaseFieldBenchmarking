!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f

       subroutine amr_timestep(dt,dtmin,dtmax,mype)


!-----------------------------------------------------------------------
! This routine computes the hydro timestep, and distributes
! it to all the processors.
!
!
! This template shows an example of how pointers may be used in PARAMESH.
!
!
! Written:      Peter MacNeice February 1997
!
!
      use paramesh_dimensions
      use physicaldata
      use tree

      include 'mpif.h'

      real, parameter :: courant = .4

!------------------------------------------

! local variables

      real :: speed2(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real :: press(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real :: maxspeed

      real :: dtl,dtmaxl,dtl0
      save :: dtl,tpwrk,dtmaxl,dtl0

      integer :: npes, ierr

      include 'pointers.fh'

      integer large_integer,intgr
      parameter(large_integer=huge(intgr)/100)

!-----------------------------------------------------------------------

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)

      dtmax = 0.
      dtmin = 1.e30
      dtl = 1.e30
      dtlevel(:) = 1.e30

!------------------------------------------


! loop over leaf blocks on this processor
      if(lnblocks.gt.0) then
      do l=1,lnblocks

      if(nodetype(l).eq.1 .or. advance_all_levels) then


        dx = bsize(1,l)/(real(nxb/2)*2)

        rho => unk(1,:,:,:,l)
        vx => unk(2,:,:,:,l)
        vy => unk(3,:,:,:,l)
        vz => unk(4,:,:,:,l)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! users timestep calaculation

        speed2(:,:,:) = (vx(:,:,:)*vx(:,:,:)+ & 
     &       vy(:,:,:)*vy(:,:,:)+ & 
     &       vz(:,:,:)*vz(:,:,:))/ & 
     &       (rho(:,:,:)*rho(:,:,:)+1.e-10)


        maxspeed=0.
        do k=1+nguard*k3d,nzb+nguard*k3d
        do j=1+nguard*k2d,nyb+nguard*k2d
        do i=1+nguard,nxb+nguard
          maxspeed=max(maxspeed,speed2(i,j,k))
        enddo
        enddo
        enddo


        dtl = courant*dx/sqrt(maxspeed+1.e-10)

! end of users timestep calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        dtlevel(lrefine(l)) = min(dtlevel(lrefine(l)),dtl)

      endif

      enddo
      endif

!------------------------------------------


! find smallest timesteps for each refinement level  across all processors
      do i=1,maxlevels
        dtl0 = dtlevel(i)
        call comm_real_min_to_all(dtlevel(i),dtl0)
      enddo

! ensure that the timestep does not increase as refinement level increases.
      do i=2,maxlevels
        dtlevel(i) = min(dtlevel(i),dtlevel(i-1))
      enddo

! ensure that each timestep is either equal to or a factor 2 larger
! than the timestep at the next higher refinement level.
      dtmin = dtlevel(maxlevels)
      dtmaxl = large_integer*dtmin
      do i=1,maxlevels-1
        dtlevel(i) = min(dtlevel(i),dtmaxl)
      enddo

      do i = maxlevels-1,1,-1
        ratio  = dtlevel(i)/dtlevel(i+1)
        iratio = max((int(ratio)/2)*2,1)
        iratio = min(iratio,2)
        dtlevel(i) = real(iratio)*dtlevel(i+1)
      enddo

      dtmaxl = 0.
      if(lnblocks.gt.0) then
      do l=1,lnblocks
      if(nodetype(l).eq.1) then
        dtmaxl = max(dtmaxl,dtlevel(lrefine(l)))
      endif
      enddo
      endif

! find largest timestep for any leaf node across all processors
      dtl0 = dtmaxl
      call comm_real_max_to_all(dtmaxl,dtl0)

      dtmax = dtmaxl

      if (var_dt) then
         dt = dtmax
      else
         dt = dtmin
         dtmax = dtmin
      end if

      if(mype.eq.0) write(*,*) 'proc ',mype,' dt ',dt

      return
      end subroutine amr_timestep
