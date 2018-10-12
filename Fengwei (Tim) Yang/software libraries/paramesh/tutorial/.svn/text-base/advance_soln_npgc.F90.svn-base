!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine advance_soln(mype,time,dt)
!
!
!--------------------------------------------------------------
! include files for amr

        use physicaldata
        use tree

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"
        include         'shmem_reduce.fh'

      integer nguard0
      parameter(nguard0 = nguard*npgs)

!--------------------------------------------------------------

! Store the assignment between block boundary fluxes and the physical
! variables.
      common/flux_assign/ iflux_target(nfluxvar)
      integer iflux_target


      logical lcc,lfc,ldiag,l_srl_only

      real    dx,dy,dz,dxi,dyi,dzi
      real    flxx(nvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                                  kl_bnd1:ku_bnd1)
      real    flxy(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                                  kl_bnd1:ku_bnd1)
      real    flxz(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                  kl_bnd1:ku_bnd1+k3d)
      real    src(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                  kl_bnd1:ku_bnd1)
      integer i,j,k,idest,nvarp,iv,iflx,ivar,ndel

      real xtest,ytest,ztest
      integer lb

!--------------------------------------------------------------

      call amr_timestep(dt,dtmin,dtmax,mype)


#ifdef NO_PERMANENT_GUARDCELLS
! Store a copy of the current solution in gt_unk
        call amr_1blk_copy_soln
#endif

      call shmem_barrier_all()


      do ivar=1,nfluxvar
        iflux_target(ivar) = ivar
      enddo
      src(:,:,:,:) = 0.


! loop over leaf grid blocks
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif


#ifdef NO_PERMANENT_GUARDCELLS
! Copy data from current block into working block and fill its guardcells
        idest = 1
        iopt = 1
        nlayers = nguard
        lcc = .true.
        lfc = .false.
        l_srl_only = .false.
        icoord = 0
        ldiag = .true.
        call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype, & 
     &                          lcc,lfc,l_srl_only,icoord,ldiag)
#endif


        dx = size(1,lb)/real(nxb-gc_off_x)
        dy = size(2,lb)/real(nyb-gc_off_y*k2d)
        dz = size(3,lb)/real(nzb-gc_off_z*k3d)
        dxi = 1./dx
        dyi = 0.
        if(ndim.ge.2) dyi = 1./dy
        dzi = 0.
        if(ndim.eq.3) dzi = 1./dz

        idest = 1

! Loop over physical variables
        do ivar=1,nvar

! compute fluxes on current block
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
     &                   unk1(ivar,i-1,j,k,idest))*dxi*dt & 
     &                   *real(ivar)

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
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------
! Section to be modified by the user
!
! y-fluxes for the diffusion equation example
!
              flxy(ivar,i,j,k) = - (unk1(ivar,i,j,k,idest) - & 
     &                   unk1(ivar,i,j-k2d,k,idest))*dyi*dt & 
     &                   *real(ivar)
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
     &                   unk1(ivar,i,j,k-k3d,idest))*dzi*dt & 
     &                   *real(ivar)
!-----------------------------------------------------------------
!---------------------------------------------^^^^ USER EDIT -----
!-----------------------------------------------------------------


              enddo
            enddo
          enddo


           ndel = nguard*(1-npgs) 
! set values for the solution at the end of timestep stage PHASE0
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard

                unk1(ivar,i,j,k,idest) = unk1(ivar,i,j,k,idest)  & 
     &          - (  & 
     &             (flxx(ivar,i+1,j,k)   - flxx(ivar,i,j,k) )*dxi & 
     &          +  (flxy(ivar,i,j+k2d,k) - flxy(ivar,i,j,k) )*dyi & 
     &          +  (flxz(ivar,i,j,k+k3d) - flxz(ivar,i,j,k) )*dzi & 
     &          + src(ivar,i,j,k)*dt)

              enddo
            enddo
          enddo


        enddo                    ! end of loop over physical variables


! Capture fluxes at the block boundaries
        do iflx = 1,nfluxvar

        ivar = iflux_target(iflx)
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


#ifdef NO_PERMANENT_GUARDCELLS
! Store new solution for this block in secondary copy
        call amr_1blk_to_perm( lcc,lfc,lb,iopt,idest)
#endif


#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo             ! end loop over grid blocks
      endif

      call shmem_barrier_all()


!--------------------
! Modify block boundary fluxes to ensure conservation
      nsub = 1
      call amr_flux_conserve(mype,nsub)

      call shmem_barrier_all()


!--------------------
! Apply any changes to the block boundary fluxes which may have 
! been made by AMR_FLUX_CONSERVE to enforce conservation.

! Loop over blocks
      do lb = 1,lnblocks


! Adjust the solution at the block boundaries using 
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

! modify flux/variable association for the current timestep 
            ivar = iflux_target(iflx)


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

      enddo
!--------------------

      call shmem_barrier_all()

! Increment time
      time = time + dt


      return
      end
