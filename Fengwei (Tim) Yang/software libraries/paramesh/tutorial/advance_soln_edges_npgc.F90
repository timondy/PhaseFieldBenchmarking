!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine advance_soln(mype,time,dt,istep)
!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
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
      parameter(nguard0=nguard*npgs)

!--------------------------------------------------------------

      logical lcc,lfc,ldiag,l_srl_only

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

      integer i,j,k,idest,iv
      real pi,xi,yi,zi


      real xtest,ytest,ztest,rvol
      integer lb

!--------------------------------------------------------------


      call amr_timestep(dt,dtmin,dtmax,mype)


#ifdef NO_PERMANENT_GUARDCELLS
! Store a copy of the current solution in gt_unk
        call amr_1blk_copy_soln
#endif


      call shmem_barrier_all()

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
        lfc = .true.
        l_srl_only = .false.
        icoord = 0
        ldiag = .true.
        call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype, & 
     &                          lcc,lfc,l_srl_only,icoord,ldiag)
#endif

        pi = 3.1415926535897932384

        dx = size(1,lb)/real(nxb-gc_off_x)
        dy = size(2,lb)/real(nyb-gc_off_y*k2d)
        dz = size(3,lb)/real(nzb-gc_off_z*k3d)

        dxi = 1./dx
        dyi = 0.
        if(ndim.ge.2) dyi = 1./dy
        dzi = 0.
        if(ndim.eq.3) dzi = 1./dz


        ex_y(:,:,:,:) = 0.
        ex_z(:,:,:,:) = 0.
        ey_x(:,:,:,:) = 0.
        ey_z(:,:,:,:) = 0.
        ez_x(:,:,:,:) = 0.
        ez_y(:,:,:,:) = 0.


! Loop over physical variables
        do iv=1,nfacevar

          iblock = lb

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

! set values for the solution at the end of timestep stage PHASE0
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard+1
              facevarx1(iv,i,j,k,idest) = & 
     &                facevarx1(iv,i,j,k,idest) - dzi*dyi*( & 
     &          (ex_z(iv,i,j,k) - ex_z(iv,i,j+k2d,k))*dz & 
     &        - (ex_y(iv,i,j,k) - ex_y(iv,i,j,k+k3d))*dy )
              enddo
            enddo
          enddo
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d+k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard
              facevary1(iv,i,j,k,idest) = & 
     &                facevary1(iv,i,j,k,idest) - dzi*dxi*( & 
     &          (ey_x(iv,i,j,k) - ey_x(iv,i,j,k+k3d))*dx & 
     &        - (ey_z(iv,i,j,k) - ey_z(iv,i+1,j,k))*dz )

              enddo
            enddo
          enddo
          do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d+k3d
            do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
              do i=il_bnd1+nguard,iu_bnd1-nguard
              facevarz1(iv,i,j,k,idest) = & 
     &                facevarz1(iv,i,j,k,idest) - dyi*dxi*( & 
     &          (ez_y(iv,i,j,k) - ez_y(iv,i+1,j,k))*dy & 
     &        - (ez_x(iv,i,j,k) - ez_x(iv,i,j+k2d,k))*dx )
              enddo
            enddo
          enddo


! Capture fluxes at the block boundaries
        bedge_facex_y(iv,1,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ex_y(iv,  1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d)

        bedge_facex_y(iv,2,1+nguard0*k2d:nyb+nguard0*k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ex_y(iv,  nxb+1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d)

        bedge_facex_z(iv,1,1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ex_z(iv,  1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d)
        bedge_facex_z(iv,2,1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ex_z(iv,  nxb+1+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                1+nguard*k3d:nzb+nguard*k3d)

        bedge_facey_x(iv,1+nguard0:nxb+nguard0,1, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ey_x(iv,  1+nguard:nxb+nguard, & 
     &                1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d)
        bedge_facey_x(iv,1+nguard0:nxb+nguard0,2, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d+k3d,iblock)= & 
     &    ey_x(iv,  1+nguard:nxb+nguard, & 
     &                nyb+1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d+k3d)
        bedge_facey_z(iv,1+nguard0:nxb+nguard0+1,1, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ey_z(iv,  1+nguard:nxb+nguard+1, & 
     &                1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d)
        bedge_facey_z(iv,1+nguard0:nxb+nguard0+1,2, & 
     &                1+nguard0*k3d:nzb+nguard0*k3d,iblock)= & 
     &    ey_z(iv,  1+nguard:nxb+nguard+1, & 
     &                nyb+1+nguard, & 
     &                1+nguard*k3d:nzb+nguard*k3d)

        bedge_facez_x(iv,1+nguard0:nxb+nguard0, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d+k2d,1,iblock)= & 
     &    ez_x(iv,  1+nguard:nxb+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                1+nguard*k3d)
        bedge_facez_x(iv,1+nguard0:nxb+nguard0, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d+k2d,2,iblock)= & 
     &    ez_x(iv,  1+nguard:nxb+nguard, & 
     &                1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &                nzb+(1+nguard)*k3d)

        bedge_facez_y(iv,1+nguard0:nxb+nguard0+1, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d,1,iblock)= & 
     &    ez_y(iv,  1+nguard:nxb+nguard+1, & 
     &                1+nguard*k2d:nyb+nguard*k2d, & 
     &                1+nguard*k3d)
        bedge_facez_y(iv,1+nguard0:nxb+nguard0+1, & 
     &                1+nguard0*k2d:nyb+nguard0*k2d,2,iblock)= & 
     &    ez_y(iv,  1+nguard:nxb+nguard+1, & 
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


! Modify block boundary edge data to ensure conservation
      nsub=1
      call amr_edge_average(mype,nsub)

      call shmem_barrier_all()

! Apply any changes to the block boundary fluxes which may have
! been made by AMR_FLUX_CONSERVE to enforce conservation.

! Loop over blocks
      do lb = 1,lnblocks

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
          do iv = 1,nfluxvar

! face 1
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarx(iv,1+ng0,j,k,lb) = & 
     &                facevarx(iv,1+ng0,j,k,lb) - dzi*dyi*( & 
     &         ( ( bedge_facex_z(iv,1,j,k,lb) & 
     &            -bedge_facex_z(iv,1,j+k2d,k,lb) )*dz & 
     &          -( bedge_facex_y(iv,1,j,k,lb) & 
     &            -bedge_facex_y(iv,1,j,k+k3d,lb) )*dy ) & 
     &       - ( ( tbedge_facex_z(iv,1,j,k,lb) & 
     &            -tbedge_facex_z(iv,1,j+k2d,k,lb) )*dz & 
     &          -( tbedge_facex_y(iv,1,j,k,lb) & 
     &            -tbedge_facex_y(iv,1,j,k+k3d,lb) )*dy ) & 
     &                                                       )
            enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              facevary(iv,1+ng0,j,k,lb) = facevary(iv,1+ng0,j,k,lb) & 
     &                - dzi*dxi*( & 
     &          -  bedge_facex_z(iv,1,j,k,lb) & 
     &          + tbedge_facex_z(iv,1,j,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarz(iv,1+ng0,j,k,lb) = facevarz(iv,1+ng0,j,k,lb) & 
     &                - dyi*dxi*( & 
     &          +  bedge_facex_y(iv,1,j,k,lb) & 
     &          - tbedge_facex_y(iv,1,j,k,lb) )*dy
            enddo
          enddo

! face 2
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarx(iv,nxb+1+ng0,j,k,lb) = & 
     &            facevarx(iv,nxb+1+ng0,j,k,lb) - dzi*dyi*( & 
     &         ( ( bedge_facex_z(iv,2,j,k,lb) & 
     &            -bedge_facex_z(iv,2,j+k2d,k,lb) )*dz & 
     &          -( bedge_facex_y(iv,2,j,k,lb) & 
     &            -bedge_facex_y(iv,2,j,k+k3d,lb) )*dy ) & 
     &       - ( ( tbedge_facex_z(iv,2,j,k,lb) & 
     &            -tbedge_facex_z(iv,2,j+k2d,k,lb) )*dz & 
     &          -( tbedge_facex_y(iv,2,j,k,lb) & 
     &            -tbedge_facex_y(iv,2,j,k+k3d,lb) )*dy ) & 
     &                                                       )
            enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              facevary(iv,nxb+ng0,j,k,lb) =  & 
     &          facevary(iv,nxb+ng0,j,k,lb) & 
     &                - dzi*dxi*( & 
     &          +  bedge_facex_z(iv,2,j,k,lb) & 
     &          - tbedge_facex_z(iv,2,j,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarz(iv,nxb+ng0,j,k,lb) =  & 
     &              facevarz(iv,nxb+ng0,j,k,lb) & 
     &                - dyi*dxi*( & 
     &          -  bedge_facex_y(iv,2,j,k,lb) & 
     &          + tbedge_facex_y(iv,2,j,k,lb) )*dy
            enddo
          enddo


! face 3 
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(iv,i,1+ng0,k,lb) = & 
     &                facevary(iv,i,1+ng0,k,lb) - dzi*dxi*( & 
     &        ( ( bedge_facey_x(iv,i,1,k,lb) & 
     &           -bedge_facey_x(iv,i,1,k+k3d,lb) )*dx & 
     &         -( bedge_facey_z(iv,i,1,k,lb) & 
     &           -bedge_facey_z(iv,i+1,1,k,lb) )*dz ) & 
     &      - ( ( tbedge_facey_x(iv,i,1,k,lb) & 
     &           -tbedge_facey_x(iv,i,1,k+k3d,lb) )*dx & 
     &         -( tbedge_facey_z(iv,i,1,k,lb) & 
     &           -tbedge_facey_z(iv,i+1,1,k,lb) )*dz ) & 
     &                                                       )
              enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(iv,i,1+ng0,k,lb) = & 
     &                facevarx(iv,i,1+ng0,k,lb) - dzi*dxi*( & 
     &         +  bedge_facey_z(iv,i,1,k,lb) & 
     &         - tbedge_facey_z(iv,i,1,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(iv,i,1+ng0,k,lb) = & 
     &                facevarz(iv,i,1+ng0,k,lb) - dzi*dxi*( & 
     &         -  bedge_facey_x(iv,i,1,k,lb) & 
     &         + tbedge_facey_x(iv,i,1,k,lb) )*dx
            enddo
          enddo
! face 4 
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(iv,i,nyb+1+ng0,k,lb) = & 
     &            facevary(iv,i,nyb+1+ng0,k,lb) - dzi*dxi*( & 
     &       ( ( bedge_facey_x(iv,i,2,k,lb) & 
     &           -bedge_facey_x(iv,i,2,k+k3d,lb) )*dx & 
     &         -( bedge_facey_z(iv,i,2,k,lb) & 
     &           -bedge_facey_z(iv,i+1,2,k,lb) )*dz ) & 
     &      - ( ( tbedge_facey_x(iv,i,2,k,lb) & 
     &           -tbedge_facey_x(iv,i,2,k+k3d,lb) )*dx & 
     &         -( tbedge_facey_z(iv,i,2,k,lb) & 
     &           -tbedge_facey_z(iv,i+1,2,k,lb) )*dz ) & 
     &                                                       )
              enddo
          enddo
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(iv,i,nyb+ng0,k,lb) = & 
     &                facevarx(iv,i,nyb+ng0,k,lb) - dzi*dxi*( & 
     &         -  bedge_facey_z(iv,i,2,k,lb) & 
     &         + tbedge_facey_z(iv,i,2,k,lb) )*dz
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(iv,i,nyb+ng0,k,lb) = & 
     &                facevarz(iv,i,nyb+ng0,k,lb) - dzi*dxi*( & 
     &         +  bedge_facey_x(iv,i,2,k,lb) & 
     &         - tbedge_facey_x(iv,i,2,k,lb) )*dx
            enddo
          enddo

          if(ndim.eq.3) then

! face 5 
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(iv,i,j,1+ng0,lb) = & 
     &                facevarz(iv,i,j,1+ng0,lb) - dyi*dxi*( & 
     &        ( ( bedge_facez_y(iv,i,j,1,lb) & 
     &           -bedge_facez_y(iv,i+1,j,1,lb) )*dy & 
     &         -( bedge_facez_x(iv,i,j,1,lb) & 
     &           -bedge_facez_x(iv,i,j+k2d,1,lb) )*dx ) & 
     &      - ( ( tbedge_facez_y(iv,i,j,1,lb) & 
     &           -tbedge_facez_y(iv,i+1,j,1,lb) )*dy & 
     &         -( tbedge_facez_x(iv,i,j,1,lb) & 
     &           -tbedge_facez_x(iv,i,j+k2d,1,lb) )*dx ) & 
     &                                                       )
              enddo
            enddo
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(iv,i,j,1+ng0,lb) = & 
     &                facevarx(iv,i,j,1+ng0,lb) - dyi*dxi*( & 
     &         -  bedge_facez_y(iv,i,j,1,lb) & 
     &         + tbedge_facez_y(iv,i,j,1,lb) )*dy
              enddo
            enddo
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(iv,i,j,1+ng0,lb) = & 
     &                facevary(iv,i,j,1+ng0,lb) - dzi*dxi*( & 
     &         +  bedge_facez_x(iv,i,j,1,lb) & 
     &         - tbedge_facez_x(iv,i,j,1,lb) )*dz
              enddo
            enddo

! face 6 
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevarz(iv,i,j,nzb+1+ng0,lb) = & 
     &            facevarz(iv,i,j,nzb+1+ng0,lb) - dyi*dxi*( & 
     &        ( ( bedge_facez_y(iv,i,j,2,lb) & 
     &           -bedge_facez_y(iv,i+1,j,2,lb) )*dy & 
     &         -( bedge_facez_x(iv,i,j,2,lb) & 
     &           -bedge_facez_x(iv,i,j+k2d,2,lb) )*dx ) & 
     &      - ( ( tbedge_facez_y(iv,i,j,2,lb) & 
     &           -tbedge_facez_y(iv,i+1,j,2,lb) )*dy & 
     &         -( tbedge_facez_x(iv,i,j,2,lb) & 
     &           -tbedge_facez_x(iv,i,j+k2d,2,lb) )*dx ) & 
     &                                                       )

              enddo
            enddo
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0+1,iu_bnd-ng0
              facevarx(iv,i,j,nzb+ng0,lb) = & 
     &                facevarx(iv,i,j,nzb+ng0,lb) - dyi*dxi*( & 
     &         +  bedge_facez_y(iv,i,j,2,lb) & 
     &         - tbedge_facez_y(iv,i,j,2,lb) )*dy
              enddo
            enddo
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(iv,i,j,nzb+ng0,lb) = & 
     &                facevary(iv,i,j,nzb+ng0,lb) - dzi*dxi*( & 
     &         -  bedge_facez_x(iv,i,j,2,lb) & 
     &         + tbedge_facez_x(iv,i,j,2,lb) )*dz
              enddo
            enddo
          endif

          enddo

      enddo


      call shmem_barrier_all()

! Increment time
      time = time + dt

! Test div B
      call gtest_neigh_data(mype,istep)


      return
      end
