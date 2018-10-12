!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
!
! This is a test program to verify the indexing in the boundary condition
! routine amr_1blk_bcset.
!

! In physicaldata.fh set
! N_DIM  2
! nxb = nyb = 4
! nvar = 2
! nguard = 2
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



#ifdef NO_PERMANENT_GUARDCELLS
        write(*,*)  & 
     &  'Error: NO_PERMANENT_GUARDCELLS must be undefined ', & 
     &  'for this test! '
        stop
#endif /* NO_PERMANENT_GUARDCELLS */


!---------------------------------------------------------------
!
! BLOCK 3


! set a limit on the refinement level
        lrefine_min = 1
        lrefine_max = 4

        iunit = 60

! set up initial grid state.

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

                neigh(1,:,1) = -21
                neigh(2,:,1) = -21

                refine(1)=.true.
        endif

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
	call amr_initial_soln_bc_test(mype)

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
        iopt = 1
        nlayers = nguard
        call amr_guardcell(mype,iopt,nlayers)


!----------------------------------------------------------
!
! BLOCK 6

! output unk
        call output2d_bc_test_tecplot( & 
     &                 lrefine_min,lrefine_max,time,dt)

! output facevarx
        call output2d_fbc_test_tecplot( & 
     &                 lrefine_min,lrefine_max,time,dt,1)
! output facevary
        call output2d_fbc_test_tecplot( & 
     &                 lrefine_min,lrefine_max,time,dt,2)
 

        call amr_close

        stop
        end

        subroutine amr_initial_soln_bc_test(mype)


! include file to define physical qualities of the model and mesh
        use physicaldata

! include file defining the tree
        use tree

! include file required for shmem library.
#include "amr_shmem.fh"

        integer mype


        if(mype.eq.0) then

        do k = kl_bnd,ku_bnd
        do j = jl_bnd,ju_bnd
        do i = il_bnd,iu_bnd
          unk(:,i,j,k,1) = real(i+j)
        enddo
        enddo
        enddo

        do k = kl_bnd,ku_bnd
        do j = jl_bnd,ju_bnd
        do i = il_bnd,iu_bnd+1
          facevarx(:,i,j,k,1) = real(i+j)
        enddo
        enddo
        enddo

        do k = kl_bnd,ku_bnd
        do j = jl_bnd,ju_bnd+k2d
        do i = il_bnd,iu_bnd
          facevary(:,i,j,k,1) = real(i+j)
        enddo
        enddo
        enddo

        if(ndim.eq.3) then
        do k = kl_bnd,ku_bnd+k3d
        do j = jl_bnd,ju_bnd
        do i = il_bnd,iu_bnd
          facevarz(:,i,j,k,1) = real(i+j)
        enddo
        enddo
        enddo
        endif

        endif

        return
        end
