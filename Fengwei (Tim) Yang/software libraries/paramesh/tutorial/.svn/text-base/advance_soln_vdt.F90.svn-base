!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!#define DEBUG

      subroutine advance_soln(mype,time,dt,istep)

!------------------------------------------------------------------------

! This is a template for a variable timestep integration
! routine, with a multi-step integration algorithm, which
! advances the solution on all refinement levels, and applies
! conservation constraints.
!
! BEWARE - This is difficult!  Be prepared to invest some time
! understanding the flow of this example.
 
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
! phases. Copy P+1 is the solution at the end of phase P. 
!
!      eg. predictor-corrector
!
!            dt
!       <--------->
!
!       |----|----|--->   t 
!
!       1    2    3       solution copy
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
! This example will work whether NO_PERMANENT_GUARDCELLS is defined
! or not.

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
      integer lb,lref,phase0,phase1,iopt,nlayers,idest
      integer icoord,ilev,iflx,ng0,np

      logical lcc,lfc,l_srl_only,ldiag

      real    dt,dt0,dtmin,dtmax
      real    time0,time1
      real    rnvarp
      real    area_xy,area_yz,area_zx
      real    dx,dy,dz,dxi,dyi,dzi,rvol
      real    xtest,ytest,ztest
      integer i,j,k,l

      real    src(1:nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                   kl_bnd1:ku_bnd1)

      integer lref_max_local,lref_max
      save    lref_max_local,lref_max

!------------------------------------------------------------------------


! You must advance the solution at all levels, if you want to use this
! example.
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
! they are to be applied. Since it is possible that not all solution
! variables will have a conservation constraint, we allow for the
! possibility that you will only allocate memory to store fluxes
! for those variables which are conserved. In this case, iflux_target
! defines the variable to which a given flux is to be applied.


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
      nvarp  = nvar/(nphase+1)

! error checking
      if(mype.eq.0) then
      rnvarp = real(nvar)/real(nphase+1)
      rnvarp = abs(rnvarp - real(nvarp))
      if(rnvarp.gt.0) then
        write(*,*) 'ERROR : nvar may not be set correctly! '
        write(*,*) 'ERROR : Are nvar and nphase consistent?'
        call abort()
      endif
      endif

!      write(*,*) 'nvarp ',nvarp

! Assume at this point all blocks are time synchronized. We are now
! about to start a global timestep. The next time all blocks are
! guaranteed to be synchronized is at the end of this global timestep.


! Step 1
! Compute a new timestep
! The global timestep is dt = dtmax, and the finest timestep is dtmin.
! It is assumed that all blocks at a given refinement level will use
! the same local timestep. Also, timesteps used for different refinement
! levels must be a power of 2 multiple of the smallest timestep, and
! that timestep is a monotonic function of refinement level with
! the shortest time associated with the finest refinement level.
      call amr_timestep(dt,dtmin,dtmax,mype)

! Step 2
! Compute the ratio of the coarsest to finest timesteps
      msub = int( (dtmax+.1*dtmin)/dtmin )

! compute no of time increments for a global timestep
      msub1 = msub*nphase
      
! Step 3
! set the time to which each block has been advanced. At this
! point all blocks have been advanced to the current time.
      time_loc(:) = time

      call shmem_barrier_all()

! Step 4
! Now advance the solution msub1 times in timestep increments of dtmin
      do nsub = 1,msub1

#ifdef DEBUG
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'Starting nsub = ',nsub
      write(*,*) ' '
      write(*,*) ' '
#endif

!--------------------
! Step 4.1
! Set times and timestep phases

! Time at beginning of this sub-timestep
        time0 = time + (dtmin/real(nphase))*(real(nsub-1)+.01)

! Time at end of this sub-timestep
        time1 = time + (dtmin/real(nphase))*real(nsub)

! Determine timestep phases for all refinement levels.
! phase_dt(ilev) stores the current phase for any blocks at refinement 
! level ilev.
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


! Solution copy 1 is where we store the solution at the 
! beginning of the current timestep, until after all finer blocks have 
! reached the end of the current timestep. If phase=1 then we are 
! starting a new timestep on this level and we can safely advance 
! copy 1 to the current time.
          lref = lrefine(lb)
          phase0 = phase_dt(lref)
          lcycle = loc_cycle(lref)
          if( lcycle.eq.1 ) then


! Set solution copies  P = 1 : nphase to be the same as copy P = nphase+1.

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
! Zero out fluxes accumulated at fine boundaries
            ttflux_x(:,:,:,:,lb) = 0.
            ttflux_y(:,:,:,:,lb) = 0.
            ttflux_z(:,:,:,:,lb) = 0.
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
! Store a copy of the current solution in GT_UNK
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
            lfc = .false.
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

! compute phase for level coarser than the current blocks level.
! This info may be needed to time-average guardcell data inside SOURCE
! and ADVANCE
            phase1 = phase0
            if(lref.gt.1) phase1 = phase_dt(lref-1)

!            write(*,*) 'dt0 nsub lb phase0 time_loc time0 ',
!     .       dt0,nsub,lb,phase0,time_loc(lb),time0

! The source terms can be computed using solutions at phases 1 : PHASE0 .
! Guardcell data is available at all these time levels. Note it is
! likely that you will need to use some average of these timelevels for
! guardcells which overlie a coarser region.
! For example, with predictor-corrector,
!   (1) compute source terms inside the block using the solution level 
!       P = PHASE0 .
!   (2) for guardcell info at a coarse boundary,
!       if PHASE0=1 on this block and PHASE0=1 on coarser neighbor, 
!          use copy 1
!       if PHASE0=2 on this block and PHASE0=1 on coarser neighbor, 
!          use (copy 1 + copy 2)/2
!       if PHASE0=1 on this block and PHASE0=2 on coarser neighbor, 
!          use copy 2
!       if PHASE0=2 on this block and PHASE0=2 on coarser neighbor, 
!          use (copy 2 + copy 3)/2

            call source( phase0, phase1, nphase, src, lb )


! Advance local solution on block lb through timestep DT0,
! Remember, the solution is advanced from copy 1, regardless of the
! current phase.
! Inside advance we must save the block boundary fluxes into FLUX_X, 
! FLUX_Y, FLUX_Z for later use when enforcing any conservation constraints.
! For guardcell info, comment (2) above applies here also.
            call advance( phase0, phase1, nphase, dt0, src, lb,istep) 

! Store solution in copy T_UNK, T_FACEVARX, etc.
! The solution at the end of the current phase cannot be stored in
! the appropriate layer of UNK etc , until all blocks updating at this
! cycle have completed, otherwise the guardcell data would not be
! time coherent. Therefore we save the update in T_UNK, etc, and
! is the next section of code we copy from T_UNK to UNK, etc.
            call amr_1blk_t_to_perm( lcc,lfc,lb,idest)

          endif                         ! end of TIME0 if test


        enddo
      call shmem_barrier_all()


!--------------------
! Step 4.5

! Loop over blocks
        do lb = 1,lnblocks


! Capture solution for current phase from temporary storage.
           lref = lrefine(lb)
           phase0 = phase_dt(lref)
           ivar1 = phase0*nvarp+1
           ivar2 = (phase0+1)*nvarp
           if(nvar.gt.0) & 
     &      unk(ivar1:ivar2,:,:,:,lb) = t_unk(ivar1:ivar2,:,:,:,lb)
           if(nfacevar.gt.0) then
           facevarx(ivar1:ivar2,:,:,:,lb)= & 
     &                                tfacevarx(ivar1:ivar2,:,:,:,lb)
           facevary(ivar1:ivar2,:,:,:,lb)= & 
     &                                tfacevary(ivar1:ivar2,:,:,:,lb)
           facevarz(ivar1:ivar2,:,:,:,lb)= & 
     &                                tfacevarz(ivar1:ivar2,:,:,:,lb)
           endif


#ifdef DEBUG
           l=lb
           xtest = (bnd_box(1,1,l)+2.05)*(bnd_box(2,1,l)+2.05)
           ytest = (bnd_box(1,2,l)+2.05)*(bnd_box(2,2,l)+2.05)
           ztest = -1.
           if(ndim.eq.3) ztest = (bnd_box(1,3,l)-.05) & 
     &                          *(bnd_box(2,3,l)-.05)
           if(nodetype(l).eq.1.and. & 
     &        xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
             write(*,*) 'block ',l,nodetype(l),bnd_box(1:2,1:ndim,l)
             write(*,*) 'ivar range ',ivar1,ivar2
             k = 1+nguard0*k3d
             do j=1,nyb+2*nguard0
             write(*,50) j,(unk(ivar1,i,j,k,l),i=1,nxb+2*nguard0)
50           format(1x,i3,6(2x,f7.4))
             enddo
           endif
#endif

        enddo

! Record whether the timestep is now complete, for each block. 
! amr_flux_conserve needs this info to determine whether to 
! overwrite or accumulate fluxes at refinement jumps.
        do lb=1,lnblocks
           ilev = lrefine(lb)
           lcycle = loc_cycle(ilev)
           ldtcomplete(lb) = .false.
           if(lcycle.eq.ncyc_local(ilev)) ldtcomplete(lb) = .true.
#ifdef DEBUG
           write(*,*) 'ldtcomplete lb ',lb , ldtcomplete(lb)
#endif

        enddo


      call shmem_barrier_all()

#ifdef DEBUG
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
             write(*,50) j,(unk(2,i,j,k,lb),i=1,nxb+2*nguard0)
             enddo
           endif
        enddo
#endif

!--------------------
! Step 4.5
! Modify block boundary fluxes to ensure conservation. amr_flux_conserve
! collects fluxes from fine neighbors, and accumulates them in
! internal storage for blocks. For any blocks which are completing 
! the last phase of their current timestep these accumulated fluxes
! are transferred to the T_FLUX_*  arrays.
      call amr_flux_conserve(mype,nsub)

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
! for the next value of nsub. 
        if(ldtcomplete(lb) & 
     &            .and.lrefine(lb).lt.lref_max) then

#ifdef DEBUG
          write(*,*) 'apply flux correct block ',lb,' nsub ',nsub, & 
     &               ' phase ',phase0
#endif

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


            unk(ivar,1+ng0,:,:,lb) =  & 
     &                        unk(ivar,1+ng0,:,:,lb) & 
     &                    - ( tflux_x(iflx,1,:,:,lb) & 
     &                       - flux_x(iflx,1,:,:,lb) )*rvol*area_yz
            unk(ivar,nxb+ng0,:,:,lb) =  & 
     &                        unk(ivar,nxb+ng0,:,:,lb) & 
     &                    + ( tflux_x(iflx,2,:,:,lb) & 
     &                       - flux_x(iflx,2,:,:,lb) )*rvol*area_yz
            if(ndim.ge.2) then
            unk(ivar,:,1+ng0*k2d,:,lb) =  & 
     &                        unk(ivar,:,1+ng0*k2d,:,lb) & 
     &                    - ( tflux_y(iflx,:,1,:,lb) & 
     &                       - flux_y(iflx,:,1,:,lb) )*rvol*area_zx
            unk(ivar,:,nyb+ng0*k2d,:,lb) =  & 
     &                        unk(ivar,:,nyb+ng0*k2d,:,lb) & 
     &                    + ( tflux_y(iflx,:,2,:,lb) & 
     &                       - flux_y(iflx,:,2,:,lb) )*rvol*area_zx
            endif
            if(ndim.eq.3) then
            unk(ivar,:,:,1+ng0*k3d,lb) =  & 
     &                        unk(ivar,:,:,1+ng0*k3d,lb) & 
     &                    - ( tflux_z(iflx,:,:,1,lb) & 
     &                       - flux_z(iflx,:,:,1,lb) )*rvol*area_xy
            unk(ivar,:,:,nzb+ng0*k3d,lb) =  & 
     &                        unk(ivar,:,:,nzb+ng0*k3d,lb) & 
     &                    + ( tflux_z(iflx,:,:,2,lb) & 
     &                       - flux_z(iflx,:,:,2,lb) )*rvol*area_xy
            endif


          enddo                          

          endif


      enddo
!--------------------


#ifdef DEBUG
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
             write(*,50) j,(unk(2,i,j,k,lb),i=1,nxb+2*nguard0)
             enddo
           endif
        enddo
#endif


      enddo                   ! end of NSUB do loop

      call shmem_barrier_all()

!--------------------
! Step 5 
! increment global time
      time = time + dt
      write(*,*) 'Time = ',time

        do lb = 1,lnblocks

! Capture solution at the end of the global timestep into solution copies
! P = 1 : nphase from copy P = nphase+1.
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


      return
      end subroutine advance_soln



      subroutine source(phase0, phase1, nphase, src, lb)

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
!  phase1             integer          current phase of timestep for
!                                        any coarse neighbors of the
!                                        current block
!  nphase             integer          no. of phases in each timestep
!  src                real array       source terms for equations
!  lb                 integer          the current block number
!
!-----------------------------------------------------------


! AMR include files

      use physicaldata
      use tree

!-------------------------------
      implicit none 

      integer phase0,phase1,lb,nphase
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

! In this example there are no source terms
      src(1:nvar,:,:,:) = 0.

!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------


      return
      end subroutine source







      subroutine advance( phase0, phase1, nphase, dt0, src ,lb,istep) 

!
! This routine advances the solution on the current grid block, lb,
! through any specified stage of the integration timestep.
! The stage is specified by PHASE0. The computed source term for this
! block is returned in the array SRC.

!-----------------------------------------------------------
!
! Arguments:
!
!  phase0             integer          current phase of timestep
!  phase1             integer          current phase of timestep for
!                                        any coarse neighbors of the
!                                        current block
!  nphase             integer          no. of phases in each timestep
!  dt0                real             time increment associated with the
!                                        current phase of the timestep
!  src                real array       source terms for equations
!  lb                 integer          the current block number
!
!-----------------------------------------------------------

!
! The fluxes computed in this example are appropriate for a set
! of NVARP 3D scalar diffusion equations with constant diffusion
! coefficients.
!
! For multi-phase algorithms the guardcell values at the interface
! with a coarser neighbor must be set with the appropriate time
! averaging. As an example, the following values would apply for
! a 2 step predictor-corrector algorithm.
!       if PHASE0=1 on this block and PHASE0=1 on coarser neighbor,
!          use copy 1
!       if PHASE0=2 on this block and PHASE0=1 on coarser neighbor,
!          use (copy 1 + copy 2)/2
!       if PHASE0=1 on this block and PHASE0=2 on coarser neighbor,
!          use copy 2
!       if PHASE0=2 on this block and PHASE0=2 on coarser neighbor,
!          use (copy 2 + copy 3)/2

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

      integer phase0,phase1,lb,nphase,istep
      real    dt0
      real    src(1:nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                   kl_bnd1:ku_bnd1)

!-------------------------------
!
! Local variables

      real    dx,dy,dz,dxi,dyi,dzi
      real    flxx(nvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                                  kl_bnd1:ku_bnd1)
      real    flxy(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                                  kl_bnd1:ku_bnd1)
      real    flxz(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                  kl_bnd1:ku_bnd1+k3d)
      integer i,j,k,idest,nvarp,iv,iflx,ivar,iface

      real xtest,ytest,ztest
      integer lb
!-------------------------------


        dx = size(1,lb)/real(nxb-gc_off_x)
        dy = size(2,lb)/real(nyb-gc_off_y*k2d)
        dz = size(3,lb)/real(nzb-gc_off_z*k3d)
        dxi = 1./dx
        dyi = 0.
        if(ndim.ge.2) dyi = 1./dy
        dzi = 0.
        if(ndim.eq.3) dzi = 1./dz

! This routine operates on 1 block at a time. The solution for this 
! block was copied into the working arrays UNK1, etc, during the
! guardcell call which preceded the call to this routine.
! By default the guardcell routines put this data into layer 1
! of the working arrays. Therefore we compute the fluxes from
! layer 1.
        idest = 1

! Compute number of independent physical variables
        nvarp = nvar/(nphase+1)


! Loop over physical variables
        do iv=1,nvarp

! compute the storage index appropriate for the current physical variable 
! for the current timestep phase.
          ivar = (phase0-1)*nvarp + iv

!------------------
!
! Construct time-averages for guardcell data as appropriate.
! Cycle over block faces. If a neighbor does not exist on a face and
! the face is not an external boundary then it is a coarser neighbor.

! To illustrate how this is done we have provided the code for the
! single step and predictor-corrector cases.


        if(nphase.eq.1) then

! For nphase=1 no action is needed here. The values currently
! in the guardcells correspond to the beginning of the local
! timestep, as required.



        elseif(nphase.eq.2) then
! Predictor-corrector
!       if PHASE0=1 on this block and PHASE0=1 on coarser neighbor,
!          use copy 1
!       if PHASE0=2 on this block and PHASE0=1 on coarser neighbor,
!          use (copy 1 + copy 2)/2
!       if PHASE0=1 on this block and PHASE0=2 on coarser neighbor,
!          use copy 2
!       if PHASE0=2 on this block and PHASE0=2 on coarser neighbor,
!          use (copy 2 + copy 3)/2
!
! the variable index ranges for each copy are as follows:
!          Copy 1    1         to   nvarp
!          Copy 2    nvarp+1   to   2*nvarp
!          Copy 3    2*nvarp+1 to   3*nvarp


! Cycle over block faces
          do iface = 1,nfaces
          if( (neigh(1,iface,lb).lt.1) .and. & 
     &        (neigh(1,iface,lb).gt.-20) ) then

! face 1
           if(iface.eq.1) then
             do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
             do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
               i=il_bnd1+nguard-1
               if(phase0.eq.2.and.phase1.eq.1) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                           ( unk1(iv,i,j,k,idest) +  & 
     &                             unk1(nvarp+iv,i,j,k,idest) )
               elseif(phase0.eq.1.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = unk1(nvarp+iv,i,j,k,idest) 
               elseif(phase0.eq.2.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                          ( unk1(nvarp+iv,i,j,k,idest) + & 
     &                            unk1(2*nvarp+iv,i,j,k,idest) )
               endif
             enddo
             enddo

! face 2 
           elseif(iface.eq.2) then
             do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
             do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
               i=iu_bnd1-nguard+1
               if(phase0.eq.2.and.phase1.eq.1) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                           ( unk1(iv,i,j,k,idest) +  & 
     &                             unk1(nvarp+iv,i,j,k,idest) )
               elseif(phase0.eq.1.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = unk1(nvarp+iv,i,j,k,idest) 
               elseif(phase0.eq.2.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                          ( unk1(nvarp+iv,i,j,k,idest) + & 
     &                            unk1(2*nvarp+iv,i,j,k,idest) )
               endif
             enddo
             enddo

! face 3 
           elseif(iface.eq.3) then
             do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
             do i=il_bnd1+nguard,iu_bnd1-nguard
               j=jl_bnd1+nguard-1
               if(phase0.eq.2.and.phase1.eq.1) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                           ( unk1(iv,i,j,k,idest) +  & 
     &                             unk1(nvarp+iv,i,j,k,idest) )
               elseif(phase0.eq.1.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = unk1(nvarp+iv,i,j,k,idest) 
               elseif(phase0.eq.2.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                          ( unk1(nvarp+iv,i,j,k,idest) + & 
     &                            unk1(2*nvarp+iv,i,j,k,idest) )
               endif
             enddo
             enddo

! face 4 
           elseif(iface.eq.4) then
             do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
             do i=il_bnd1+nguard,iu_bnd1-nguard
               j=ju_bnd1-nguard+1
               if(phase0.eq.2.and.phase1.eq.1) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                           ( unk1(iv,i,j,k,idest) +  & 
     &                             unk1(nvarp+iv,i,j,k,idest) )
               elseif(phase0.eq.1.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = unk1(nvarp+iv,i,j,k,idest) 
               elseif(phase0.eq.2.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                          ( unk1(nvarp+iv,i,j,k,idest) + & 
     &                            unk1(2*nvarp+iv,i,j,k,idest) )
               endif
             enddo
             enddo

! face 5 
           elseif(iface.eq.5) then
             do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
             do i=il_bnd1+nguard,iu_bnd1-nguard
               k=kl_bnd1+nguard-1
               if(phase0.eq.2.and.phase1.eq.1) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                           ( unk1(iv,i,j,k,idest) +  & 
     &                             unk1(nvarp+iv,i,j,k,idest) )
               elseif(phase0.eq.1.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = unk1(nvarp+iv,i,j,k,idest) 
               elseif(phase0.eq.2.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                          ( unk1(nvarp+iv,i,j,k,idest) + & 
     &                            unk1(2*nvarp+iv,i,j,k,idest) )
               endif
             enddo
             enddo

! face 6 
           elseif(iface.eq.6) then
             do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
             do i=il_bnd1+nguard,iu_bnd1-nguard
               k=ku_bnd1-nguard+1
               if(phase0.eq.2.and.phase1.eq.1) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                           ( unk1(iv,i,j,k,idest) +  & 
     &                             unk1(nvarp+iv,i,j,k,idest) )
               elseif(phase0.eq.1.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = unk1(nvarp+iv,i,j,k,idest) 
               elseif(phase0.eq.2.and.phase1.eq.2) then
                 unk1(ivar,i,j,k,idest) = .5* & 
     &                          ( unk1(nvarp+iv,i,j,k,idest) + & 
     &                            unk1(2*nvarp+iv,i,j,k,idest) )
               endif
             enddo
             enddo

           endif                     ! end of iface if test

          endif
          enddo                      ! end of loop over block faces

        endif                        ! end of nphase if test
!------------------


          
! compute fluxes on current block, for timestep stage PHASE0
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard+1

!
!-----------------------------------------------------------------
!---------------------------------------------<<<< USER EDIT -----
!-----------------------------------------------------------------
! Section to be modified by the user
!
! x-fluxes for the diffusion equation example
!
              flxx(ivar,i,j,k) = - (unk1(ivar,i,j,k,idest) - & 
     &                   unk1(ivar,i-1,j,k,idest))*dxi*dt0 & 
     &                   *real(iv)

!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------
              enddo
            enddo
          enddo

          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d+k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard
!-----------------------------------------------------------------
!---------------------------------------------<<<< USER EDIT -----
!-----------------------------------------------------------------
! Section to be modified by the user
!
! y-fluxes for the diffusion equation example
!
              flxy(ivar,i,j,k) = - (unk1(ivar,i,j,k,idest) - & 
     &                   unk1(ivar,i,j-k2d,k,idest))*dyi*dt0 & 
     &                   *real(iv)
!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------
              enddo
            enddo
          enddo
        

          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d+k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard
!-----------------------------------------------------------------
!---------------------------------------------<<<< USER EDIT -----
!-----------------------------------------------------------------
! Section to be modified by the user
!
! z-fluxes for the diffusion equation example
!
              flxz(ivar,i,j,k) = - (unk1(ivar,i,j,k,idest) - & 
     &                   unk1(ivar,i,j,k-k3d,idest))*dzi*dt0 & 
     &                   *real(iv)
!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------


              enddo
            enddo
          enddo


#ifdef DEBUG
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
             write(*,50) j,(unk1(1,i,j,k,1),i=1,nxb+2*nguard)
50           format(1x,i3,6(2x,f7.4))
!50           format(1x,i3,7(2x,f7.4))
             enddo
           endif
#endif

! set values for the solution at the end of timestep stage PHASE0
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard

                unk1(ivar+nvarp,i,j,k,idest) = unk1(iv,i,j,k,idest)  & 
     &          - (  & 
     &             (flxx(ivar,i+1,j,k)   - flxx(ivar,i,j,k) )*dxi & 
     &          +  (flxy(ivar,i,j+k2d,k) - flxy(ivar,i,j,k) )*dyi & 
     &          +  (flxz(ivar,i,j,k+k3d) - flxz(ivar,i,j,k) )*dzi & 
     &          + src(ivar,i,j,k)*dt0)


              enddo
            enddo
          enddo

#ifdef DEBUG
        write(*,*) 'advanced ',iv,' to ',ivar+nvarp,' block ',lb
#endif
        enddo                    ! end of loop over physical variables


! Capture fluxes at the block boundaries. Note that the arrays
! FLUX_* may or may not have guardcell memory allocated, hence
! the use of NGUARD) rather than NGUARD on the left side of each
! assignment. The temporary flux arrays flx*, however, do have
! memory allocated for guardcell storage and so the indexing on the
! right uses NGUARD.
        do iflx = 1,nfluxvar

        ivar = (phase0-1)*nvarp + iflux_target(iflx)
        flux_x(iflx,1,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,lb)=  & 
     &    flxx(ivar,  1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 
        flux_x(iflx,2,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,lb)=  & 
     &    flxx(ivar,  nxb+1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 
        flux_y(iflx,  1+nguard0:nxb+nguard0, & 
     &                1, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,lb)= & 
     &    flxy(ivar,  1+nguard:nxb+nguard, & 
     &                1+nguard*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 
        flux_y(iflx,  1+nguard0:nxb+nguard0, & 
     &                2, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,lb)= & 
     &    flxy(ivar,  1+nguard:nxb+nguard, & 
     &                nyb+(1+nguard)*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d) 
        flux_z(iflx,  1+nguard0:nxb+nguard0, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1,lb) =  & 
     &    flxz(ivar,  1+nguard:nxb+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d)
        flux_z(iflx,  1+nguard0:nxb+nguard0, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                2,lb) =  & 
     &    flxz(ivar,  1+nguard:nxb+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                nzb+(1+nguard)*k3d)


        enddo


      return
      end subroutine advance

