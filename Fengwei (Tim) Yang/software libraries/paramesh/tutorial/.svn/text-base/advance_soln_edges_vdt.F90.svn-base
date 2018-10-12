!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!#define NOTNOW

      subroutine advance_soln(mype,time,dt,istep)

!------------------------------------------------------------------------

! This is a template for a variable timestep integration
! routine, with a multi-step integration algorithm, which
! advances the solution on all refinement levels, and applies
! conservation constraints.
!
! Please read the instructions here, and in the Users Manual before
! using this routine.
!
! The user is required to modify this in the following ways:
!   1. set NPHASE, the number of steps in a timestep.
!   2. specify the association between block boundary fluxes and
!      the physical variables.
!   3. modify a loop within a routine called SOURCE which computes source
!      terms at any required times.
!   4. modify the loop which computes fluxes in a routine called ADVANCE,
!      which advances the solution through any stage of a timestep.
! The places where these mods are required are clearly identified.


!------------------------------------------------------------------------

!
! To use different timesteps at different levels, you must 
! advance the solution on all refinement levels, not just leaf
! blocks.

! It is also assumed that the timesteps on blocks at a given refinement
! level are either the same or a factor of 2 larger than on blocks which 
! are 1 level more refined.
!

! This code example illustrates how to advance the solution using 
! different timesteps for each refinement level, for an algorithm
! which itself has multiple time levels, ie predictor-corrector,
! Runge-Kutta, etc.
! In the comments, we use P to designate the phase of the timestep.
! For example, predictor-corrector timestep has 2 phases, the predictor
! step which is phase 1, and the corrector step which is phase 2.
! We need to keep multiple copies of the solution for each of these
! phases. Copy P is the solution at the end of phase P. 
!
!      eg. predictor-corrector
!
!            dt
!       <--------->
!
!       |----|----|--->   t 
!
!       0    1    2       solution copy
!
!
! If our algorithm has N phases, then we need N+1 copies of the solution.
!
! The N+1 copies are stored by setting NVAR to be (N+1)*NVARP, where
! NVARP is the actual number of physical variables which constitute
! the solution at any given time.
!
! You must ensure that all N+1 copies have the same initial solution
! at the beginning of the first timestep, so remember to initialize
! all copies in your initialization routine.
!
!------------------------------------------------------------------------


!
! include files for amr

      use physicaldata
      use tree

      implicit none

      integer mype,istep
      real    time,dt

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"
        include         'shmem_reduce.fh'


      integer nguard0
      parameter(nguard0 = nguard*npgs)



! Store the assignment between block boundary fluxes and the physical
! variables.
      common/flux_assign/ iflux_target(nfluxvar)
      integer iflux_target


! local variables

      integer nphase,nvarp,msub,msub1,nsub,nsubs_per_phase
      integer ivar,ivar1,ivar2,ivar3,ivar4
      integer iphase,lcycle
      integer lb,lref,phase0,iopt,nlayers,idest
      integer icoord,ilev,iflx,ng0,np

      logical lcc,lfc,l_srl_only,ldiag

      real    dt,dt0,dtmin,dtmax
      real    time0,time1
      real    rnvarp
      real    area_xy,area_yz,area_zx
      real    dx,dy,dz,dxi,dyi,dzi,rvol
      real    xtest,ytest,ztest
      integer i,j,k,l,iii

      real    src(1:nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                   kl_bnd1:ku_bnd1)

      integer lref_max_local,lref_max
      save    lref_max_local,lref_max

!------------------------------------------------------------------------


#ifdef ADVANCE_ALL_LEVELS


!-----------------------------------------------------------------
!---------------------------------------------<<<< USER EDIT -----
!-----------------------------------------------------------------

! Set no. of stages(or phases) in a timestep
!     nphase = 1                                      ! single stage
      nphase = 2                                      ! predictor-corrector
!     nphase = 4                                      ! 4th order runge-kutta

!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------


! set the correspondence between fluxes (FLUX_X, FLUX_Y, FLUX_Z)
! stored at block boundaries and the solution variables to which
! they are to be applied.

!-----------------------------------------------------------------
!---------------------------------------------<<<< USER EDIT -----
!-----------------------------------------------------------------

      do iflx = 1,nfluxvar

      iflux_target(iflx) = iflx        ! In this example the i-th flux
                                       ! is associated with the i-th physical
                                       ! variable.
      enddo

!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------

! Compute the maximum refinement level
      lref_max_local = maxval(lrefine)
      call shmem_barrier_all()
      call comm_int_max_to_all(lref_max,lref_max_local)
      call shmem_barrier_all()


! Compute the no. of physical variables
      nvarp  = nfacevar/(nphase+1)

! error checking
      if(mype.eq.0) then
      rnvarp = real(nfacevar)/real(nphase+1)
      rnvarp = abs(rnvarp - real(nvarp))
      if(rnvarp.gt.0) then
        write(*,*) 'ERROR : nfacevar may not be set correctly! '
        write(*,*) 'ERROR : Are nfacevar and nphase consistent?'
        call abort()
      endif
      endif

!      write(*,*) 'nvarp ',nvarp

! Assume at this point all blocks are time synchronized


! Step 1
! Compute a new timestep
! The global timestep is dt = dtmax, and the finest timestep is dtmin
! Rewrite amr_timestep to ensure that timesteps are powers of 2 multiples of
! the minimum timestep. 
      call amr_timestep(dt,dtmin,dtmax,mype)

! Step 2
! Compute the ratio of the coarsest to finest timesteps
      msub = int( (dtmax+.1*dtmin)/dtmin )

! compute no of time increments for a global timestep
      msub1 = msub*nphase
      
! Step 3
! set the time to which each block has been advanced
      time_loc(:) = time

      call shmem_barrier_all()

! Step 4
! Now advance the solution msub1 times in timestep increments of dtmin/nphase
      do nsub = 1,msub1

      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'Starting nsub = ',nsub
      write(*,*) ' '
      write(*,*) ' '

!--------------------
! Step 4.1
! Set times and timestep phases

! Time at beginning of this sub-timestep
        time0 = time + (dtmin/real(nphase))*(real(nsub-1)+.01)

! Time at end of this sub-timestep
        time1 = time + (dtmin/real(nphase))*real(nsub)

! Determine timestep phases for all refinement levels
        do ilev = 1,maxlevels
          ncyc_local(ilev) = int((dtlevel(ilev)+.01*dtmin)/dtmin) & 
     &                       *nphase
          loc_cycle(ilev) = mod(nsub-1,ncyc_local(ilev))+1
          nsubs_per_phase = ncyc_local(ilev)/nphase
          iphase = (nsub-1)/nsubs_per_phase
          phase_dt(ilev) = mod(iphase,nphase) + 1
        enddo


!--------------------
! Step 4.2
! Loop over blocks, reseting the solutions wherever a block is
! about to start a new timestep.

        do lb = 1,lnblocks


! Solution copy 0 is where we store the solution at the 
! beginning of the current timestep, until after all finer blocks have 
! reached the end of the current timestep. If phase=1 then we are 
! starting a new timestep on this level and we can safely advance 
! copy 0 to the current time.
          lref = lrefine(lb)
          phase0 = phase_dt(lref)
          lcycle = loc_cycle(lref)
          if( lcycle.eq.1 ) then


! Set solution copies  P = 0 : nphase-1 to be the same as copy P = nphase.

          do np = 1,nphase

            ivar1 = (np-1)*nvarp + 1
            ivar2 = (np-1)*nvarp + nvarp
            ivar3 = nphase*nvarp + 1
            ivar4 = nphase*nvarp + nvarp
            if(nvar.gt.0) & 
     &        unk(ivar1:ivar2,:,:,:,lb) = unk(ivar3:ivar4,:,:,:,lb)
            if(nfacevar.gt.0) then
              facevarx(ivar1:ivar2,:,:,:,lb) = & 
     &          facevarx(ivar3:ivar4,:,:,:,lb)
              facevary(ivar1:ivar2,:,:,:,lb) = & 
     &          facevary(ivar3:ivar4,:,:,:,lb)
              facevarz(ivar1:ivar2,:,:,:,lb) = & 
     &          facevarz(ivar3:ivar4,:,:,:,lb)
            endif


          enddo

          endif                         ! end of lcycle if test

#ifdef VAR_DT
          if( lcycle.eq.1 ) then
!          write(*,*) 'reset ttbedge : blk ',lb,
!     .       '  phase0 time0 time_loc ', phase0,time0,time_loc(lb)
! Zero out fluxes accumulated at fine boundaries
           if(nvar.gt.0.and.nfluxes.gt.0) then
            ttflux_x(:,:,:,:,lb) = 0.
            ttflux_y(:,:,:,:,lb) = 0.
            ttflux_z(:,:,:,:,lb) = 0.
           endif
           if(nfacevar.gt.0.and.nedgevar.gt.0) then
            ttbedge_facex_y(:,:,:,:,lb) = 0.
            ttbedge_facex_z(:,:,:,:,lb) = 0.
            ttbedge_facey_x(:,:,:,:,lb) = 0.
            ttbedge_facey_z(:,:,:,:,lb) = 0.
            ttbedge_facez_x(:,:,:,:,lb) = 0.
            ttbedge_facez_y(:,:,:,:,lb) = 0.
           endif

          endif                         ! end of lcycle if test

#endif

        enddo


        call shmem_barrier_all()


!--------------------
! Step 4.3
! Perform a global guardcell filling, if permanent guardcell storage
! is allocated - otherwise make the necessary copy of the solution
! which the 1blk_guardcell routines use.
 
#ifdef NO_PERMANENT_GUARDCELLS
! Store a copy of the current solution in GT_UNK, GT_FACEVAR*
        call amr_1blk_copy_soln
#else
        iopt = 1
        nlayers = nguard
        call amr_guardcell(mype,iopt,nlayers)
#endif


!--------------------
! Step 4.4
! The main loop over blocks which computes the source terms and
! actually advances the solution through this phase of the timestep.


        do lb = 1,lnblocks

! Does the current block need to be advanced?
          if(time0.ge.time_loc(lb)) then

            idest = 1
            lcc = .true.
            lfc = .true.
            iopt = 1
#ifdef NO_PERMANENT_GUARDCELLS
! Copy data from current block into working block and fill its guardcells
            nlayers = nguard
            l_srl_only = .false.
            icoord = 0
            ldiag = .true.
            call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype, & 
     &                          lcc,lfc,l_srl_only,icoord,ldiag)
#else

            call amr_perm_to_1blk( lcc,lfc,lb,mype,iopt,idest)
#endif

! compute sub-timestep to use for this block
            lref = lrefine(lb)
            phase0 = phase_dt(lref)
            dt0 = dtlevel(lref)*real(phase0)/real(nphase)

!            write(*,*) 'dt0 nsub lb phase0 time_loc time0 ',
!     .       dt0,nsub,lb,phase0,time_loc(lb),time0

! The source terms can be computed using solutions at phases 0 : PHASE0 - 1.
! Guardcell data is available at all these time levels. Note it is
! likely that you will need to use some average of these timelevels for
! guardcells which overlie a coarser region.
! For example, with predictor-corrector,
!   (1) compute source terms inside the block using the solution level 
!       P = PHASE0 - 1.
!   (2) for guardcell info at a coarse boundary,
!       if PHASE0=1 on this block and PHASE0=1 on coarser neighbor, 
!          use copy 0
!       if PHASE0=2 on this block and PHASE0=1 on coarser neighbor, 
!          use (copy 0 + copy 1)/2
!       if PHASE0=1 on this block and PHASE0=2 on coarser neighbor, 
!          use copy 1
!       if PHASE0=2 on this block and PHASE0=2 on coarser neighbor, 
!          use (copy 1 + copy 2)/2

            call source( phase0, nphase, src, lb )


! Advance local solution on block lb through timestep DT0,
! Advance from copy 0.
! Save block boundary fluxes into FLUX_X, FLUX_Y, FLUX_Z for later use.
            call advance( phase0, nphase, dt0, src, lb,istep ) 

! Store solution in copy T_UNK, T_FACEVARX, etc.
            call amr_1blk_t_to_perm( lcc,lfc,lb,idest)

          endif                         ! end of TIME0 if test


        enddo


      call shmem_barrier_all()


!--------------------
! Step 4.5

! Loop over blocks
        do lb = 1,lnblocks

!pmn          if(time0.ge.time_loc(lb)) then

! Capture solution for current phase from temporary storage.
           lref = lrefine(lb)
           phase0 = phase_dt(lref)
           ivar1 = phase0*nvarp+1
           ivar2 = (phase0+1)*nvarp
           if(nvar.gt.0) & 
     &      unk(ivar1:ivar2,:,:,:,lb) = t_unk(ivar1:ivar2,:,:,:,lb)
           if(nfacevar.gt.0) then
           facevarx(ivar1:ivar2,:,:,:,lb)= & 
     &                              tfacevarx(ivar1:ivar2,:,:,:,lb)
           facevary(ivar1:ivar2,:,:,:,lb)= & 
     &                              tfacevary(ivar1:ivar2,:,:,:,lb)
           facevarz(ivar1:ivar2,:,:,:,lb)= & 
     &                              tfacevarz(ivar1:ivar2,:,:,:,lb)
           endif




!pmn          endif                         ! end of time0 if test
        enddo

! Record whether the timestep is now complete, for each block. 
! amr_flux_conserve needs this info to determine whether to 
! overwrite or accumulate fluxes at refinement jumps.
        do lb=1,lnblocks
           ilev = lrefine(lb)
           lcycle = loc_cycle(ilev)
           ldtcomplete(lb) = .false.
           if(lcycle.eq.ncyc_local(ilev)) ldtcomplete(lb) = .true.
           write(*,*) 'ldtcomplete lb ',lb , ldtcomplete(lb)
        enddo


      call shmem_barrier_all()

#ifdef NOTNOW
        do lb=1,lnblocks
           xtest = (bnd_box(1,1,lb)+2.05)*(bnd_box(2,1,lb)+2.05)
           ytest = (bnd_box(1,2,lb)+2.05)*(bnd_box(2,2,lb)+2.05)
           ztest = -1.
           if(ndim.eq.3) ztest = (bnd_box(1,3,lb)-.05) & 
     &                          *(bnd_box(2,3,lb)-.05)
           if(nodetype(lb).eq.1.and. & 
     &        xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
          write(*,*) 'block ',lb,nodetype(lb),bnd_box(1:2,1:ndim,lb)
      write(*,*) 'bef flux nsub ',nsub
             k = 1+nguard0*k3d
             do j=1,nyb+2*nguard0
!             write(*,50) j,(unk(2,i,j,k,lb),i=1,nxb+2*nguard0)
!50           format(1x,i3,6(2x,f7.4))
!50           format(1x,i3,7(2x,f7.4))
             enddo
           endif
        enddo
#endif

!--------------------
! Step 4.5
! Modify block boundary fluxes to ensure conservation
      call amr_flux_conserve(mype,nsub)

      call shmem_barrier_all()

! Modify block boundary edge data to ensure conservation
      call amr_edge_average(mype,nsub)

      call shmem_barrier_all()


!--------------------
! Step 4.6
! Apply any changes to the block boundary fluxes which may have 
! been made by AMR_FLUX_CONSERVE to enforce conservation.

! Loop over blocks
      do lb = 1,lnblocks

! Has the next finer level finished its timestep, requiring this block
! to apply corrected fluxes ?

        lref = lrefine(lb)
        phase0 = phase_dt(lref)
        ivar1 = phase0*nvarp+1
        ivar2 = (phase0+1)*nvarp


! Compute the time at the end of this blocks local timestep
        if(time0.ge.time_loc(lb)) then

          dt0 = dtlevel(lref)/real(nphase)
          time_loc(lb) = time_loc(lb) + dt0

        endif


! only apply flux correction if this block will begin a new timestep 
! for the next value of nsub. ??
        if(ldtcomplete(lb) & 
     &            .and.lrefine(lb).lt.lref_max) then

!          write(*,*) 'apply flux correct block ',lb,' nsub ',nsub,
!     .               ' phase ',phase0

! Adjust the solution for the current phase at the block boundaries using 
! the corrected conservative fluxes.
          dx = size(1,lb)/real(nxb-gc_off_x)
          dy = size(2,lb)/real(nyb-gc_off_y*k2d)
          dz = size(3,lb)/real(nzb-gc_off_z*k3d)
          if(ndim.lt.2) dy = 1.
          if(ndim.lt.3) dz = 1.
          dxi = 1./dx
          dyi = 1./dy
          dzi = 1./dz
          rvol = dxi*dyi*dzi
          area_xy = dx*dy
          area_yz = dy*dz
          area_zx = dz*dx

          ng0 = nguard*npgs


! Loop over the stored fluxes, applying them to correct the appropriate
! solution variables.
          do iflx = 1,nfluxvar

! modify flux/variable association for the current timestep phase
            ivar = phase0*nvarp + iflux_target(iflx)

!          write(*,*) 'flux correction on blk  : ',lb,
!     .        ' ivar ',ivar


! face 1
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarx(ivar,1+ng0,j,k,lb) = & 
     &                facevarx(ivar,1+ng0,j,k,lb) - dzi*dyi*( & 
     &         ( ( bedge_facex_z(iflx,1,j,k,lb) & 
     &            -bedge_facex_z(iflx,1,j+k2d,k,lb) )*dz & 
     &          -( bedge_facex_y(iflx,1,j,k,lb) & 
     &            -bedge_facex_y(iflx,1,j,k+k3d,lb) )*dy ) & 
     &       - ( ( tbedge_facex_z(iflx,1,j,k,lb) & 
     &            -tbedge_facex_z(iflx,1,j+k2d,k,lb) )*dz & 
     &          -( tbedge_facex_y(iflx,1,j,k,lb) & 
     &            -tbedge_facex_y(iflx,1,j,k+k3d,lb) )*dy ) & 
     &                                                       )
            enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              facevary(ivar,1+ng0,j,k,lb) = facevary(ivar,1+ng0,j,k,lb) & 
     &                - dzi*dxi*( & 
     &          -  bedge_facex_z(iflx,1,j,k,lb) & 
     &          + tbedge_facex_z(iflx,1,j,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarz(ivar,1+ng0,j,k,lb) = facevarz(ivar,1+ng0,j,k,lb) & 
     &                - dyi*dxi*( & 
     &          +  bedge_facex_y(iflx,1,j,k,lb) & 
     &          - tbedge_facex_y(iflx,1,j,k,lb) )*dy
            enddo
          enddo

! face 2
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarx(ivar,nxb+1+ng0,j,k,lb) = & 
     &            facevarx(ivar,nxb+1+ng0,j,k,lb) - dzi*dyi*( & 
     &         ( ( bedge_facex_z(iflx,2,j,k,lb) & 
     &            -bedge_facex_z(iflx,2,j+k2d,k,lb) )*dz & 
     &          -( bedge_facex_y(iflx,2,j,k,lb) & 
     &            -bedge_facex_y(iflx,2,j,k+k3d,lb) )*dy ) & 
     &       - ( ( tbedge_facex_z(iflx,2,j,k,lb) & 
     &            -tbedge_facex_z(iflx,2,j+k2d,k,lb) )*dz & 
     &          -( tbedge_facex_y(iflx,2,j,k,lb) & 
     &            -tbedge_facex_y(iflx,2,j,k+k3d,lb) )*dy ) & 
     &                                                       )
            enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              facevary(ivar,nxb+ng0,j,k,lb) = & 
     &          facevary(ivar,nxb+ng0,j,k,lb) & 
     &                - dzi*dxi*( & 
     &          +  bedge_facex_z(iflx,2,j,k,lb) & 
     &          - tbedge_facex_z(iflx,2,j,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarz(ivar,nxb+ng0,j,k,lb) = & 
     &              facevarz(ivar,nxb+ng0,j,k,lb) & 
     &                - dyi*dxi*( & 
     &          -  bedge_facex_y(iflx,2,j,k,lb) & 
     &          + tbedge_facex_y(iflx,2,j,k,lb) )*dy
            enddo
          enddo


! face 3
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(ivar,i,1+ng0,k,lb) = & 
     &                facevary(ivar,i,1+ng0,k,lb) - dzi*dxi*( & 
     &        ( ( bedge_facey_x(iflx,i,1,k,lb) & 
     &           -bedge_facey_x(iflx,i,1,k+k3d,lb) )*dx & 
     &         -( bedge_facey_z(iflx,i,1,k,lb) & 
     &           -bedge_facey_z(iflx,i+1,1,k,lb) )*dz ) & 
     &      - ( ( tbedge_facey_x(iflx,i,1,k,lb) & 
     &           -tbedge_facey_x(iflx,i,1,k+k3d,lb) )*dx & 
     &         -( tbedge_facey_z(iflx,i,1,k,lb) & 
     &           -tbedge_facey_z(iflx,i+1,1,k,lb) )*dz ) & 
     &                                                       )
              enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(ivar,i,1+ng0,k,lb) = & 
     &                facevarx(ivar,i,1+ng0,k,lb) - dzi*dxi*( & 
     &         +  bedge_facey_z(iflx,i,1,k,lb) & 
     &         - tbedge_facey_z(iflx,i,1,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(ivar,i,1+ng0,k,lb) = & 
     &                facevarz(ivar,i,1+ng0,k,lb) - dzi*dxi*( & 
     &         -  bedge_facey_x(iflx,i,1,k,lb) & 
     &         + tbedge_facey_x(iflx,i,1,k,lb) )*dx
            enddo
          enddo
! face 4
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(ivar,i,nyb+1+ng0,k,lb) = & 
     &            facevary(ivar,i,nyb+1+ng0,k,lb) - dzi*dxi*( & 
     &       ( ( bedge_facey_x(iflx,i,2,k,lb) & 
     &           -bedge_facey_x(iflx,i,2,k+k3d,lb) )*dx & 
     &         -( bedge_facey_z(iflx,i,2,k,lb) & 
     &           -bedge_facey_z(iflx,i+1,2,k,lb) )*dz ) & 
     &      - ( ( tbedge_facey_x(iflx,i,2,k,lb) & 
     &           -tbedge_facey_x(iflx,i,2,k+k3d,lb) )*dx & 
     &         -( tbedge_facey_z(iflx,i,2,k,lb) & 
     &           -tbedge_facey_z(iflx,i+1,2,k,lb) )*dz ) & 
     &                                                       )
              enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(ivar,i,nyb+ng0,k,lb) = & 
     &                facevarx(ivar,i,nyb+ng0,k,lb) - dzi*dxi*( & 
     &         -  bedge_facey_z(iflx,i,2,k,lb) & 
     &         + tbedge_facey_z(iflx,i,2,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(ivar,i,nyb+ng0,k,lb) = & 
     &                facevarz(ivar,i,nyb+ng0,k,lb) - dzi*dxi*( & 
     &         +  bedge_facey_x(iflx,i,2,k,lb) & 
     &         - tbedge_facey_x(iflx,i,2,k,lb) )*dx
            enddo
          enddo

          if(ndim.eq.3) then

! face 5
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(ivar,i,j,1+ng0,lb) = & 
     &                facevarz(ivar,i,j,1+ng0,lb) - dyi*dxi*( & 
     &        ( ( bedge_facez_y(iflx,i,j,1,lb) & 
     &           -bedge_facez_y(iflx,i+1,j,1,lb) )*dy & 
     &         -( bedge_facez_x(iflx,i,j,1,lb) & 
     &           -bedge_facez_x(iflx,i,j+k2d,1,lb) )*dx ) & 
     &      - ( ( tbedge_facez_y(iflx,i,j,1,lb) & 
     &           -tbedge_facez_y(iflx,i+1,j,1,lb) )*dy & 
     &         -( tbedge_facez_x(iflx,i,j,1,lb) & 
     &           -tbedge_facez_x(iflx,i,j+k2d,1,lb) )*dx ) & 
     &                                                       )
              enddo
            enddo
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(ivar,i,j,1+ng0,lb) = & 
     &                facevarx(ivar,i,j,1+ng0,lb) - dyi*dxi*( & 
     &         -  bedge_facez_y(iflx,i,j,1,lb) & 
     &         + tbedge_facez_y(iflx,i,j,1,lb) )*dy
              enddo
            enddo
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(ivar,i,j,1+ng0,lb) = & 
     &                facevary(ivar,i,j,1+ng0,lb) - dzi*dxi*( & 
     &         +  bedge_facez_x(iflx,i,j,1,lb) & 
     &         - tbedge_facez_x(iflx,i,j,1,lb) )*dz
              enddo
            enddo

! face 6
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(ivar,i,j,nzb+1+ng0,lb) = & 
     &            facevarz(ivar,i,j,nzb+1+ng0,lb) - dyi*dxi*( & 
     &        ( ( bedge_facez_y(iflx,i,j,2,lb) & 
     &           -bedge_facez_y(iflx,i+1,j,2,lb) )*dy & 
     &         -( bedge_facez_x(iflx,i,j,2,lb) & 
     &           -bedge_facez_x(iflx,i,j+k2d,2,lb) )*dx ) & 
     &      - ( ( tbedge_facez_y(iflx,i,j,2,lb) & 
     &           -tbedge_facez_y(iflx,i+1,j,2,lb) )*dy & 
     &         -( tbedge_facez_x(iflx,i,j,2,lb) & 
     &           -tbedge_facez_x(iflx,i,j+k2d,2,lb) )*dx ) & 
     &                                                       )
              enddo
            enddo
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(ivar,i,j,nzb+ng0,lb) = & 
     &                facevarx(ivar,i,j,nzb+ng0,lb) - dyi*dxi*( & 
     &         +  bedge_facez_y(iflx,i,j,2,lb) & 
     &         - tbedge_facez_y(iflx,i,j,2,lb) )*dy
              enddo
            enddo
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(ivar,i,j,nzb+ng0,lb) = & 
     &                facevary(ivar,i,j,nzb+ng0,lb) - dzi*dxi*( & 
     &         -  bedge_facez_x(iflx,i,j,2,lb) & 
     &         + tbedge_facez_x(iflx,i,j,2,lb) )*dz
              enddo
            enddo
          endif




          enddo                          


          endif


      enddo
!--------------------


#ifdef NOTNOW
        do lb=1,lnblocks
           xtest = (bnd_box(1,1,lb)+2.05)*(bnd_box(2,1,lb)+2.05)
           ytest = (bnd_box(1,2,lb)+2.05)*(bnd_box(2,2,lb)+2.05)
           ztest = -1.
           if(ndim.eq.3) ztest = (bnd_box(1,3,lb)-.05) & 
     &                          *(bnd_box(2,3,lb)-.05)
           if(nodetype(lb).eq.1.and. & 
     &        xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
          write(*,*) 'block ',lb,nodetype(lb),bnd_box(1:2,1:ndim,lb)
      write(*,*) 'nsub ',nsub
             k = 1+nguard0*k3d
             do j=1,nyb+2*nguard0
!             write(*,50) j,(unk(2,i,j,k,lb),i=1,nxb+2*nguard0)
             enddo
           endif
        enddo
#endif

      ivar = nphase*nvarp + 1
!      if(lcycle.eq.ncyc_local(ilev)) 
!     .            call conserve_check_1blk(mype,istep,ivar)


      enddo                   ! end of NSUB do loop

      call shmem_barrier_all()

!--------------------
! Step 5 
! increment global time
      time = time + dt
      write(*,*) 'Time = ',time

        do lb = 1,lnblocks


! Capture solution at the end of the global timestep into solution copies
! P = 0 : nphase-1 from copy P = nphase.
          do np = 1,nphase
            ivar1 = (np-1)*nvarp + 1
            ivar2 = (np-1)*nvarp + nvarp
            ivar3 = nphase*nvarp + 1
            ivar4 = nphase*nvarp + nvarp

            if(nvar.gt.0) & 
     &        unk(ivar1:ivar2,:,:,:,lb) = unk(ivar3:ivar4,:,:,:,lb)
            if(nfacevar.gt.0) then
            facevarx(ivar1:ivar2,:,:,:,lb)= & 
     &                               facevarx(ivar3:ivar4,:,:,:,lb)
            facevary(ivar1:ivar2,:,:,:,lb)= & 
     &                               facevary(ivar3:ivar4,:,:,:,lb)
            facevarz(ivar1:ivar2,:,:,:,lb)= & 
     &                               facevarz(ivar3:ivar4,:,:,:,lb)
            endif
          enddo

        enddo

#else
! error trapping
      if(mype.eq.0) then
      write(*,*)  'Error: this time advance only works if all' & 
     &           ,' refinement levels are being advanced! '
      endif
      call shmem_barrier_all()
      call abort()

#endif

! Test div B
      call gtest_neigh_data(mype,istep)


      return
      end subroutine advance_soln



      subroutine source(phase0, nphase, src, iblock)

!-----------------------------------------------------------
!
! This routine computes source terms for this block for any specified
! stage of the integration timestep.
! The stage is specified by PHASE0. The computed source term for this
! block is returned in the array SRC.
! SRC are source terms per unit volume per unit time.
!
!-----------------------------------------------------------
!
! Arguments:
!
!  phase0             integer          current phase of timestep
!  nphase             integer          no. of phases in each timestep
!  src                real array       source terms for equations
!  iblock             integer          the current block number
!
!-----------------------------------------------------------


! AMR include files

      use physicaldata
      use tree

!-------------------------------
      implicit none 

      integer phase0,iblock,nphase
      real    src(1:nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                   kl_bnd1:ku_bnd1)

!-------------------------------
! Local variables

      integer nvarp

!-------------------------------


! Compute number of physical variables
        nvarp = nvar/(nphase+1)

!-----------------------------------------------------------------
!---------------------------------------------<<<< USER EDIT -----
!-----------------------------------------------------------------
! Section to be modified by the user

! In this diffusion equation example there are no source terms
      src(1:nvar,:,:,:) = 0.

!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------


      return
      end subroutine source







      subroutine advance( phase0, nphase, dt0, src ,iblock,istep) 

!
! This routine advances the solution on the current grid block, iblock,
! through any specified stage of the integration timestep.
! The stage is specified by PHASE0. The computed source term for this
! block is returned in the array SRC.

!-----------------------------------------------------------
!
! Arguments:
!
!  phase0             integer          current phase of timestep
!  nphase             integer          no. of phases in each timestep
!  dt0                real             time increment associated with the
!                                        current phase of the timestep
!  src                real array       source terms for equations
!  iblock             integer          the current block number
!
!-----------------------------------------------------------

! AMR include files

      use physicaldata
      use tree

      implicit none 

      integer nguard0
      parameter(nguard0 = nguard*npgs)

!-------------------------------

! Store the assignment between block boundary fluxes and the physical
! variables.
      common/flux_assign/ iflux_target(nfluxvar)
      integer iflux_target

      integer phase0,iblock,nphase,istep
      real    dt0
      real    src(1:nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                   kl_bnd1:ku_bnd1)

!-------------------------------
!
! Local variables

      real    dx,dy,dz,dxi,dyi,dzi

      real    ex_y(nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                                  kl_bnd1:ku_bnd1+k3d)
      real    ex_z(nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                                  kl_bnd1:ku_bnd1)
      real    ey_x(nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                                  kl_bnd1:ku_bnd1+k3d)
      real    ey_z(nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                                  kl_bnd1:ku_bnd1)
      real    ez_x(nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                                  kl_bnd1:ku_bnd1+k3d)
      real    ez_y(nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                                  kl_bnd1:ku_bnd1+k3d)

      integer i,j,k,idest,nvarp,iv,iflx,ivar,ndel
      real pi,xi,yi,zi


      real xtest,ytest,ztest,rvol
      integer lb
!-------------------------------

        pi = 3.1415926535897932384

        dx = size(1,iblock)/real(nxb-gc_off_x)
        dy = size(2,iblock)/real(nyb-gc_off_y*k2d)
        dz = size(3,iblock)/real(nzb-gc_off_z*k3d)
        dxi = 1./dx
        dyi = 0.
        if(ndim.ge.2) dyi = 1./dy
        dzi = 0.
        if(ndim.eq.3) dzi = 1./dz
        if(ndim.eq.2) rvol = 1./(dx*dy) 
        if(ndim.eq.3) rvol = 1./(dx*dy*dz) 

        idest = 1

! Compute number of physical variables
        nvarp = nvar/(nphase+1)


        ex_y(:,:,:,:) = 0.
        ex_z(:,:,:,:) = 0.
        ey_x(:,:,:,:) = 0.
        ey_z(:,:,:,:) = 0.
        ez_x(:,:,:,:) = 0.
        ez_y(:,:,:,:) = 0.

! Loop over physical variables
        do iv=1,nvarp

          ivar = (phase0-1)*nvarp + iv

          
! compute edge data on current block, for timestep stage PHASE0
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d+k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard+1
                yi = bnd_box(1,2,iblock) + dy*(j-nguard-1)
                ex_y(:,i,j,k) = sin(2.*pi*yi/8.)
!                ex_y(:,i,j,k) = 0.
                if(ndim.eq.3) ez_y(:,i,j,k) = ex_y(:,i,j,k)
              enddo
            enddo
          enddo

          if(ndim.eq.3) then
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d+k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard+1
                zi = bnd_box(1,3,iblock) + dz*(k-nguard-1)
                ey_z(:,i,j,k) = sin(2.*pi*zi/8.)
!                ey_z(:,i,j,k) = 0.
                ex_z(:,i,j,k) = ey_z(:,i,j,k)
              enddo
            enddo
          enddo
          endif

          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d+k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d+k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard
                xi = bnd_box(1,1,iblock) + dx*(i-nguard-1)
                ey_x(:,i,j,k) = sin(2.*pi*xi/8.)
!                ey_x(:,i,j,k) = 0.
                if(ndim.eq.3) ez_x(:,i,j,k) = ey_x(:,i,j,k)
              enddo
            enddo
          enddo


#ifdef NOTNOW
           lb = iblock
           xtest = (bnd_box(1,1,lb)+2.05)*(bnd_box(2,1,lb)+2.05)
           ytest = (bnd_box(1,2,lb)+2.05)*(bnd_box(2,2,lb)+2.05)
           ztest = -1.
           if(ndim.eq.3) ztest = (bnd_box(1,3,lb)-.05) & 
     &                          *(bnd_box(2,3,lb)-.05)
           if(nodetype(lb).eq.1.and. & 
     &        xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
          write(*,*) 'block ',lb,nodetype(lb),bnd_box(1:2,1:ndim,lb)
          write(*,*) 'during advance before update'
             k = 1+nguard*k3d
             do j=1,nyb+2*nguard
!             write(*,50) j,(unk1(1,i,j,k,1),i=1,nxb+2*nguard)
50           format(1x,i3,6(2x,f7.4))
!50           format(1x,i3,7(2x,f7.4))
             enddo
           endif
#endif

           ndel = nguard*(1-npgs) 
! set values for the solution at the end of timestep stage PHASE0
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard+1
              facevarx1(ivar+nvarp,i,j,k,idest) =  & 
     &                facevarx1(iv,i,j,k,idest) - dzi*dyi*(  & 
     &          (ex_z(ivar,i,j,k) - ex_z(ivar,i,j+k2d,k))*dz & 
     &        - (ex_y(ivar,i,j,k) - ex_y(ivar,i,j,k+k3d))*dy )
              enddo
            enddo
          enddo
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d+k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard
              facevary1(ivar+nvarp,i,j,k,idest) =  & 
     &                facevary1(iv,i,j,k,idest) - dzi*dxi*(  & 
     &          (ey_x(ivar,i,j,k) - ey_x(ivar,i,j,k+k3d))*dx & 
     &        - (ey_z(ivar,i,j,k) - ey_z(ivar,i+1,j,k))*dz )

              enddo
            enddo
          enddo
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d+k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard
              facevarz1(ivar+nvarp,i,j,k,idest) =  & 
     &                facevarz1(iv,i,j,k,idest) - dyi*dxi*(  & 
     &          (ez_y(ivar,i,j,k) - ez_y(ivar,i+1,j,k))*dy & 
     &        - (ez_x(ivar,i,j,k) - ez_x(ivar,i,j+k2d,k))*dx )
              enddo
            enddo

          enddo


!#ifdef NOTNOW
           lb = iblock
           xtest = (bnd_box(1,1,lb)+2.05)*(bnd_box(2,1,lb)+2.05)
           ytest = (bnd_box(1,2,lb)+2.05)*(bnd_box(2,2,lb)+2.05)
           ztest = -1.
           if(ndim.eq.3) ztest = (bnd_box(1,3,lb)-.05) & 
     &                          *(bnd_box(2,3,lb)-.05)
           if(nodetype(lb).eq.1.and. & 
     &        xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
          write(*,*) 'block ',lb,nodetype(lb),bnd_box(1:2,1:ndim,lb)
          write(*,*) 'during advance after update'
             k = 1+nguard*k3d
             do j=1,nyb+2*nguard
!             write(*,50) j,(unk1(ivar+nvarp,i,j,k,1),
!     .                                             i=1,nxb+2*nguard)
50           format(1x,i3,7(2x,f7.4))
             enddo
           endif
!#endif

!        write(*,*) 'advanced ',iv,' to ',ivar+nvarp,' block ',iblock
        enddo                    ! end of loop over physical variables


! Capture fluxes at the block boundaries
        do iflx = 1,nfluxvar

        ivar = (phase0-1)*nvarp + iflux_target(iflx)

        bedge_facex_y(iflx,1,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ex_y(ivar,  1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d) 
        bedge_facex_y(iflx,2,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ex_y(ivar,  nxb+1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d) 

        bedge_facex_z(iflx,1,1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ex_z(ivar,  1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 
        bedge_facex_z(iflx,2,1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ex_z(ivar,  nxb+1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 

        bedge_facey_x(iflx,1+nguard0:nxb+nguard0,1, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ey_x(ivar,  1+nguard:nxb+nguard, & 
     &                1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d) 
        bedge_facey_x(iflx,1+nguard0:nxb+nguard0,2, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ey_x(ivar,  1+nguard:nxb+nguard, & 
     &                nyb+1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d) 

        bedge_facey_z(iflx,1+nguard0:nxb+nguard0+1,1, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ey_z(ivar,  1+nguard:nxb+nguard+1, & 
     &                1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 
        bedge_facey_z(iflx,1+nguard0:nxb+nguard0+1,2, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ey_z(ivar,  1+nguard:nxb+nguard+1, & 
     &                nyb+1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 

        bedge_facez_x(iflx,1+nguard0:nxb+nguard0, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d+k2d,1,iblock)= & 
     &    ez_x(ivar,  1+nguard:nxb+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                1+nguard*k3d)
        bedge_facez_x(iflx,1+nguard0:nxb+nguard0, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d+k2d,2,iblock)= & 
     &    ez_x(ivar,  1+nguard:nxb+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                nzb+(1+nguard)*k3d)

        bedge_facez_y(iflx,1+nguard0:nxb+nguard0+1, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d,1,iblock)= & 
     &    ez_y(ivar,  1+nguard:nxb+nguard+1, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d)
        bedge_facez_y(iflx,1+nguard0:nxb+nguard0+1, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d,2,iblock)= & 
     &    ez_y(ivar,  1+nguard:nxb+nguard+1, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                nzb+(1+nguard)*k3d)


        enddo


      return
      end subroutine advance

