!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!#define MULTIPLE_LEVEL1
!#define NOTNOW

!
! This is a template main program for the multidimensional PARAMESH amr package.
!

! In physicaldata.fh set
! N_DIM  2
! nxb = nyb = 4
! nvar = 1
! maxblocks = 100
! nguard = 1
! nfacevar = 0

! In workspace.fh set
! nguard_work = 2
! nvar_work = 1

!----------------------------------------------------------------
!
! BLOCK 1


! include file to define physical qualities of the model and mesh
        use physicaldata

! include file defining the tree
        use tree

! include file required for shmem library.
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"


! local amr variables
        integer nprocs,iunit
        integer shmem_my_pe,shmem_n_pes

        logical restart

!---------------------------------------------------------------
!
! BLOCK 2


! make initialization call for amr package
#ifdef SGI_SHMEM
        call comm_start(maxprocs,nprocs,mype)
#endif
        call amr_initialize


        mype = shmem_my_pe()
        nprocs = shmem_n_pes()

        if(mype.eq.0) write(*,*) 'Running on ',nprocs,' processors'


#ifndef NO_PERMANENT_GUARDCELLS
        write(*,*)  & 
     &  'Error: NO_PERMANENT_GUARDCELLS must be defined ', & 
     &  'for this test! '
        stop
#endif /* NO_PERMANENT_GUARDCELLS */


!---------------------------------------------------------------
!
! BLOCK 3


! set a limit on the refinement level
        lrefine_min = 1
        lrefine_max = 4

        restart = .false.
!       restart = .true.

        iunit = 60
        if(.not.restart) then

! set up initial grid state.

#ifndef MULTIPLE_LEVEL1
! set a single block covering the whole square domain
        lnblocks = 0
        if(mype.eq.0.) then
                lnblocks = 1
                coord(:,1) = 0.
		bnd_box(1,:,1) = -4.
		bnd_box(2,:,1) = 4.
                size(:,1) = bnd_box(2,:,1) - bnd_box(1,:,1)
                nodetype(1) = 1
                lrefine(1) = 1

                neigh(1,:,1) = 1                  ! initial block is its own
                neigh(2,:,1) = 0                  ! neighbor ie periodic bc's
                neigh(1,:,1) = -21
                neigh(2,:,1) = -21

                refine(1)=.true.
        endif
#endif

#ifdef MULTIPLE_LEVEL1
! sets up multiple blocks at the coarsest refinement level

        if(mype.eq.0) then

          if(ndim.eq.2) then

            lnblocks = 4

            size(1:2,1:lnblocks) = 4.

            coord(1:2,1) = -2.
            bnd_box(1,1:2,1) = -4.
            bnd_box(2,1:2,1) = 0.
            neigh(1,1,1) = -21
            neigh(1,2,1) = 2
            neigh(1,3,1) = -21
            neigh(1,4,1) = 3

            coord(1,2) = 2.
            coord(2,2) = -2.
            bnd_box(1,1,2) = 0.
            bnd_box(2,1,2) = 4.
            bnd_box(1,2,2) = -4.
            bnd_box(2,2,2) = 0.
            neigh(1,1,2) = 1
            neigh(1,2,2) = -21
            neigh(1,3,2) = -21
            neigh(1,4,2) = 4

            coord(1,3) = -2.
            coord(2,3) = 2.
            bnd_box(1,1,3) = -4.
            bnd_box(2,1,3) = 0.
            bnd_box(1,2,3) = 0.
            bnd_box(2,2,3) = 4.
            neigh(1,1,3) = -21
            neigh(1,2,3) = 4
            neigh(1,3,3) = 1
            neigh(1,4,3) = -21

            coord(1:2,4) = 2.
            bnd_box(1,1:2,4) = 0.
            bnd_box(2,1:2,4) = 4.
            neigh(1,1,4) = 3
            neigh(1,2,4) = -21
            neigh(1,3,4) = 2
            neigh(1,4,4) = -21

        
        elseif(ndim.eq.3) then

            lnblocks =8 

            size(1:3,1:lnblocks) = 4.

            coord(1:3,1) = -2.
            bnd_box(1,1:3,1) = -4.
            bnd_box(2,1:3,1) = 0.
            neigh(1,1,1) = -21
            neigh(1,2,1) = 2
            neigh(1,3,1) = -21
            neigh(1,4,1) = 3
            neigh(1,5,1) = -21
            neigh(1,6,1) = 5

            coord(1,2) = 2.
            coord(2,2) = -2.
            coord(3,2) = -2.
            bnd_box(1,1,2) = 0.
            bnd_box(2,1,2) = 4.
            bnd_box(1,2,2) = -4.
            bnd_box(2,2,2) = 0.
            bnd_box(1,3,2) = -4.
            bnd_box(2,3,2) = 0.
            neigh(1,1,2) = 1
            neigh(1,2,2) = -21
            neigh(1,3,2) = -21
            neigh(1,4,2) = 4
            neigh(1,5,2) = -21
            neigh(1,6,2) = 6

            coord(1,3) = -2.
            coord(2,3) = 2.
            coord(3,3) = -2.
            bnd_box(1,1,3) = -4.
            bnd_box(2,1,3) = 0.
            bnd_box(1,2,3) = 0.
            bnd_box(2,2,3) = 4.
            bnd_box(1,3,3) = -4.
            bnd_box(2,3,3) = 0.
            neigh(1,1,3) = -21
            neigh(1,2,3) = 4
            neigh(1,3,3) = 1
            neigh(1,4,3) = -21
            neigh(1,5,3) = -21
            neigh(1,6,3) = 7

            coord(1:2,4) = 2.
            coord(3,4) = -2.
            bnd_box(1,1:2,4) = 0.
            bnd_box(2,1:2,4) = 4.
            bnd_box(1,3,4) = -4.
            bnd_box(2,3,4) = 0.
            neigh(1,1,4) = 3
            neigh(1,2,4) = -21
            neigh(1,3,4) = 2
            neigh(1,4,4) = -21
            neigh(1,5,4) = -21
            neigh(1,6,4) = 8


            coord(1:2,5) = -2.
            coord(3,5) = 2.
            bnd_box(1,1:2,5) = -4.
            bnd_box(2,1:2,5) = 0.
            bnd_box(1,3,5) = 0.
            bnd_box(2,3,5) = 4.
            neigh(1,1,5) = -21
            neigh(1,2,5) = 6
            neigh(1,3,5) = -21
            neigh(1,4,5) = 7
            neigh(1,5,5) = 1
            neigh(1,6,5) = -21

            coord(1,6) = 2.
            coord(2,6) = -2.
            coord(3,6) = 2.
            bnd_box(1,1,6) = 0.
            bnd_box(2,1,6) = 4.
            bnd_box(1,2,6) = -4.
            bnd_box(2,2,6) = 0.
            bnd_box(1,3,6) = 0.
            bnd_box(2,3,6) = 4.
            neigh(1,1,6) = 5
            neigh(1,2,6) = -21
            neigh(1,3,6) = -21
            neigh(1,4,6) = 8
            neigh(1,5,6) = 2
            neigh(1,6,6) = -21

            coord(1,7) = -2.
            coord(2,7) = 2.
            coord(3,7) = 2.
            bnd_box(1,1,7) = -4.
            bnd_box(2,1,7) = 0.
            bnd_box(1,2,7) = 0.
            bnd_box(2,2,7) = 4.
            bnd_box(1,3,7) = 0.
            bnd_box(2,3,7) = 4.
            neigh(1,1,7) = -21
            neigh(1,2,7) = 8
            neigh(1,3,7) = 5 
            neigh(1,4,7) = -21
            neigh(1,5,7) = 1 
            neigh(1,6,7) = -21

            coord(1:3,8) = 2.
            bnd_box(1,1:3,8) = 0.
            bnd_box(2,1:3,8) = 4.
            neigh(1,1,8) = 7
            neigh(1,2,8) = -21
            neigh(1,3,8) = 6
            neigh(1,4,8) = -21
            neigh(1,5,8) = 4
            neigh(1,6,8) = -21


        endif

         do lb = 1,lnblocks
          nodetype(lb) = 1
          lrefine(lb) = 1
          neigh(2,:,lb) = 0
          refine(lb)=.true.
         enddo

        endif

#endif /*  MULTIPLE_LEVEL1 */


! Now cycle over blocks, refining all existing leaf blocks.

        do loop_count=1,2

                refine(1:lnblocks) = .true.
                call shmem_barrier_all()

! refine grid and apply morton reordering to grid blocks if necessary
                call amr_refine_derefine

        enddo
        call shmem_barrier_all()

        else                                ! restart
        call amr_checkpoint_re(iunit)
        endif                               ! restart

        call shmem_barrier_all()


	if(mype.eq.0) write(*,*) 'pe / blk / blk-coords / blk-sizes'
        call shmem_barrier_all()
        do l=1,lnblocks
        write(*,*) mype,l,(coord(i,l),i=1,ndim),(size(j,l),j=1,ndim) & 
     &            ,bnd_box(1:2,1,l)
        enddo

!---------------------------------------------------------------
!
! BLOCK 4


! Establish solution on the initial grid.
	if(.not.restart) call amr_initial_soln

        call amr_test_refinement(mype,lrefine_min,lrefine_max)
        call amr_refine_derefine
        iopt = 1
        nlayers = nguard
!        call amr_prolong(mype,iopt,nlayers)
        write(*,*) 'new lnblocks ',lnblocks
	if(.not.restart) call amr_initial_soln
        istep=0
        call conserve_check(mype,istep)
!        unk(1:nvar,:,:,:,1:lnblocks)=1.

	if(mype.eq.0) write(*,*) 'pe / blk / blk-coords / blk-sizes'
        call shmem_barrier_all()
        do l=1,lnblocks
        write(*,*) mype,l,(coord(i,l),i=1,ndim),(size(j,l),j=1,ndim) & 
     &            ,bnd_box(1:2,1,l)
        enddo


!---------------------------------------------------------------
!
! BLOCK 5


! exchange guardcell information
!        iopt = 1
!        nlayers = nguard
!        call amr_guardcell(mype,iopt,nlayers)


!----------------------------------------------------------
!
! BLOCK 6

         time = 0.

         minstp = 1
         maxstp = 250

! loop over time steps.
        do istep = minstp, maxstp

!        if(istep.eq.2) then
!               unk(1:nvar,:,:,:,1:lnblocks) = 1.
!               time = 0.
!        endif

! advance the solution using a User provided routine
        call advance_soln(mype,time,dt,istep)

        if(istep.eq.maxstp) & 
     &  call output2d_tecplot(lrefine_min,lrefine_max,time,dt)
 
#ifndef ADVANCE_ALL_LEVELS
        iempty = 0 
        iopt = 1
        call amr_restrict(mype,iopt,iempty)
#endif
!
! test to see if additional refinement or derefinement is necessary
! note - a call to guardcell must come before this call to ensure that
! the refinement test can be done on parents of leafs also. This avoids
! a potential refinement/derefinement flip-flop happenning on successive
! timesteps.

        if(istep.lt.2) then
!        call amr_test_refinement(mype,lrefine_min,lrefine_max)
        endif

#ifdef NOTNOW 
       if(mype.eq.0) write(*,*) ' pe  blk    refine  derefine', & 
     &                           '  curr.ref.level'
       call shmem_barrier_all()
       do l=1,lnblocks
         write(*,51) mype,l,refine(l),derefine(l),lrefine(l)
51      format(1x,i3,2x,i3,2x,l8,2x,l8,10x,i3)
       enddo
#endif


! refine grid and apply morton reordering to grid blocks if necessary
        call amr_refine_derefine

! prolong solution to any new leaf blocks if necessary
        iopt = 1
        nlayers = nguard
        call amr_prolong(mype,iopt,nlayers)

        write(*,*) 'after prolong'
        call conserve_check(mype,istep)

! exchange guardcell information
!        call amr_guardcell(mype,iopt,nlayers)


        if(mype.eq.0) write(*,*) 'iteration ',istep, & 
     &                           ' no of blocks = ',lnblocks

        if(mype.eq.0) write(*,*) 'pe / blk / blk-coords / blk-sizes'
        call shmem_barrier_all()
        do l=1,lnblocks
        write(*,*) mype,l,(coord(i,l),i=1,ndim),(size(j,l),j=1,ndim)
        enddo


        do l=1,lnblocks
           xtest = (bnd_box(1,1,l)-.05)*(bnd_box(2,1,l)-.05)
           ytest = (bnd_box(1,2,l)-.05)*(bnd_box(2,2,l)-.05)
           ztest = -1.
           if(ndim.eq.3) ztest = (bnd_box(1,3,l)-.05) & 
     &                          *(bnd_box(2,3,l)-.05)
           if(nodetype(l).eq.1.and. & 
     &        xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
             write(*,*) 'block ',l,nodetype(l),bnd_box(1:2,1:ndim,l)
             k = 1+nguard*npgs*k3d
             do j=1,nyb+2*nguard*npgs
             write(*,51) j,(unk(1,i,j,k,l),i=1,nxb+2*nguard*npgs)
50           format(1x,i3,4(2x,f7.4))
51           format(1x,i3,5(2x,f7.4))
             enddo
           endif
        enddo

        call conserve_check(mype,istep)

        end do
!---------------------------------------------------------------

!        call output2d_tecplot(lrefine_min,lrefine_max,time,dt)

!        call amr_checkpoint_wr(iunit)

        call amr_close

        stop
        end
