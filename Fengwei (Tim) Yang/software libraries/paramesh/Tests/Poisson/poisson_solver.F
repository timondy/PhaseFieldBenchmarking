! This program uses the oct-tree/block-based code to solve
! Poisson's equation in 3D. It was written to demo the amr
! package.
!
! The program uses a conjugate gradient scheme to compute 
! an approximate solution on each grid. The potential is
! stored as variable 1 in the solution array unk.
!
! A point source mass is located inside a cubic domain, and the
! program starts by using a single grid block on one processor
! to determine the gravitational potential resulting from it.
!
! The program does a double iteration. Inside the conjugate
! gradient routine the program iterates on each grid block
! keeping the block's guard cell values fixed until some accuracy 
! threshold is achieved. Then the blocks update their
! guard cell data and test the same accuracy criterion, and
! repaeat this double loop until the guard cell update no
! longer causes the accuracy threshold to be exceeded.
!
! Written: 	Peter MacNeice, Kevin Olson	January 1997
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	program poisson_3d

! include file to define physical qualities of the model and mesh
#include "physicaldata.fh"

! include file defining the tree 
	include 'tree.fh'


! include file required for shmem library.
	include 'amr_shmem.fh'
	include 'shmem_reduce.fh'



! local variables
        real	source(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                 kl_bnd1:ku_bnd1,maxblocks)

	integer nprocs,my_proc
	integer	loop,niterations
	integer shmem_my_pe,shmem_n_pes,num_procs
	integer	iopt,nlayers,icoord

! storage space for diagnostic curves
	parameter(maxcurves=100,nmax=(nxb+2)*maxcurves)
	real	output(2,0:nmax-1)
	real	coutput(2,0:nmax-1)
	integer	nc,ncurves

	logical	lstop,lperiodic
        logical lcc,lfc,l_srl_only,ldiag

	save my_proc,output,coutput,nc,ncurves

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call amr_initialize


	my_proc = shmem_my_pe()
	nprocs = shmem_n_pes()

	if(my_proc.eq.0) write(*,*) 'Running on ',nprocs,' processors'

! Arguments for call to amr_1blk_guardcell
        lcc = .true.
        lfc = .false.
        l_srl_only = .false.
        ldiag = .false.
        icoord = 0

! set up initial grid, one block on processor 0, no blocks anywhere else.
	lnblocks = 0
	lperiodic=.false.
	if(my_proc.eq.0.) then
		lnblocks = 1
       		size(:,1)=1.
                coord(:,1) = .5
                bnd_box(1,:,1) = 0.
                bnd_box(2,:,1) = 1.
		nodetype(1) = 1
		neigh(:,:,1) = -20
		if(lperiodic) then
			neigh(1,:,1) = 1
			neigh(2,:,1) = 0
		endif
	endif


! set arguments to choose unk as the dataset to be used in guardcell
! and prolong.
	iopt=1
	nlayers=1


! initialize solution to poisson's equation
	unk(:,:,:,:,:) = 0.


	if(my_proc.eq.0) open(unit=9,status='new',file='idl1.dat')

	call shmem_barrier_all()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin iteration to solution
	ncurves=0
	niterations = 10
	do loop=1,niterations

! compute source terms on leaf blocks only
	call msource(my_proc,source)
	call shmem_barrier_all()


! perform a single iteration on each leaf block
	call conjugate_grad(my_proc,loop,iopt,nlayers,source)
	call shmem_barrier_all()

	if(my_proc.eq.0.and.loop.eq.1) then
#ifdef NO_PERMANENT_GUARDCELLS
        call amr_1blk_guardcell(my_proc,iopt,nlayers,1,my_proc,lcc,lfc,
     .                              l_srl_only,icoord,ldiag)
#else
        unk1(:,:,:,:,1) = unk(:,:,:,:,1)
#endif
	ncurves=ncurves+1
	j=nyb/2+1+nguard
	k=nzb/2+1+nguard
        dx = size(1,1)/real((nxb/2)*2)
	do i=0+nguard,nxb+1+nguard
	ii = i-nguard
        xi = coord(1,1) - dx*(real(nxb-1)*.5-real(ii-1))
	index = ii + (ncurves-1)*(nxb+2)
	output(1,index)=xi
!	output(2,index)=unk(1,i,j,k,1)
	output(2,index)=unk1(1,i,j,k,1)
	enddo
	endif

	call shmem_barrier_all()


	if((loop.gt.1).and.(lnblocks.gt.0)) then
	do l=1,lnblocks
		do kk=1+nguard,nzb+nguard
		do jj=1+nguard,nyb+nguard
		do ii=1+nguard,nxb+nguard
			if(nodetype(l).eq.1.and.
     .					source(ii,jj,kk,l).gt.1.) then
				j=jj
				k=kk
				refine(l)=.true.
			endif
		enddo
		enddo
		enddo
	if(refine(l)) then
	write(*,*) 'write: refine : proc ',my_proc,' block ',l
#ifdef NO_PERMANENT_GUARDCELLS
        call amr_1blk_guardcell(my_proc,iopt,nlayers,l,my_proc,lcc,lfc,
     .                              l_srl_only,icoord,ldiag)
#else
        unk1(:,:,:,:,1) = unk(:,:,:,:,l)
#endif
		ncurves=ncurves+1
        	dx = size(1,l)/real((nxb/2)*2)
		do i=0+nguard,nxb+1+nguard
		ii = i - nguard
       		xi = coord(1,l) - dx*(real(nxb-1)*.5-real(ii-1))
		index = ii + (ncurves-1)*(nxb+2)
		output(1,index)=xi
!		output(2,index)=unk(1,i,j,k,l)
		output(2,index)=unk1(1,i,j,k,1)
		enddo
	endif
	enddo
	endif
50	format(f9.5,2x,f12.5)


	if(loop.eq.niterations) goto 500
! test to see if additional refinement is necessary
	if(my_proc.eq.0.and.loop.eq.1) refine(1)=.true.

	call shmem_barrier_all()

! refine grid if necessary
	call amr_refine_derefine
	call shmem_barrier_all()

! reorder distribution of blocks across processors to balance workload
!	call amr_morton_order (nprocs,lnblocks_old)
!	call shmem_barrier_all()

! prolong solution to any new leaf blocks if necessary
	call amr_prolong(my_proc,iopt,nlayers)
	call shmem_barrier_all()

#ifndef NO_PERMANENT_GUARDCELLS
	call amr_guardcell(my_proc,iopt,nlayers)
	call shmem_barrier_all()
#endif

500	continue

	enddo
! end iteration to solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! collect diagnostic output from all processors to proc 0
	call shmem_barrier_all()
	if(my_proc.eq.0) then
		do k=0,nprocs-1
		if(k.eq.0) then
			nc=ncurves
			coutput(:,:)=output(:,:)
		else
			call shmem_get(nc,ncurves,1,k)
			call shmem_get(coutput,output,nmax*2,k)
		endif
!		write(*,*) 'output: k= ',k,' nc = ',nc
		if(nc.gt.0) then
			index=0
			do j=1,nc
			do i=0,nxb+1
			write(9,50) coutput(1,index),coutput(2,index)
			index=index+1
			enddo
			enddo
		endif
		enddo
	endif
	call shmem_barrier_all()

	if(my_proc.eq.0) close(unit=9)
	if(my_proc.eq.0) open(unit=9,status='new',file='idl2.dat')

	call shmem_barrier_all()
        ncurves=0
	if(lnblocks.gt.0) then
        do l=1,lnblocks
#ifdef NO_PERMANENT_GUARDCELLS
        call amr_1blk_guardcell(my_proc,iopt,nlayers,l,my_proc,lcc,lfc,
     .                              l_srl_only,icoord,ldiag)
#else
        unk1(:,:,:,:,1) = unk(:,:,:,:,l)
#endif
                refine(l)=.false.
                if(nodetype(l).eq.1) then
                dx = size(1,l)/real((nxb/2)*2)
                dy = size(2,l)/real((nyb/2)*2)
                dz = size(3,l)/real((nzb/2)*2)
                do kk=1+nguard,nzb+nguard
			k0 = kk-nguard
                        zi = coord(3,l) - dz*
     .				(real(nzb-1)*.5-real(k0-1))
                do jj=1+nguard,nyb+nguard
			j0 = jj-nguard
                        yi = coord(2,l) - dy*
     .				(real(nyb-1)*.5-real(j0-1))
		        if((yi.eq..5).and.(zi.eq..5)) then
                                j=jj
                                k=kk
                                refine(l)=.true.
                        endif
                enddo
                enddo
                endif
        if(refine(l)) then
                ncurves=ncurves+1
                do i=0+nguard,nxb+1+nguard
		i0 = i-nguard
                xi = coord(1,l) - dx*(real(nxb-1)*.5-real(i0-1))
                index = i0 + (ncurves-1)*(nxb+2)
                output(1,index)=xi
!                output(2,index)=unk(1,i,j,k,l)
                output(2,index)=unk1(1,i,j,k,1)
                enddo
        endif
        enddo
        endif

! collect diagnostic output from all processors to proc 0
	call shmem_barrier_all()
        if(my_proc.eq.0) then
                do k=0,nprocs-1
                if(k.eq.0) then
                        nc=ncurves
                        coutput(:,:)=output(:,:)
                else
                        call shmem_get(nc,ncurves,1,k)
                        call shmem_get(coutput,output,nmax*2,k)
                endif
!               write(*,*) 'output: k= ',k,' nc = ',nc
                if(nc.gt.0) then
                        index=0
                        do j=1,nc
                        do i=0,nxb+1
                        write(9,50) coutput(1,index),coutput(2,index)
                        index=index+1
                        enddo
                        enddo
                endif
                enddo
        endif

	call shmem_barrier_all()

	if(my_proc.eq.0) close(unit=9)

        write(*,*) 'npgs ',npgs

	call shmem_barrier_all()

        call amr_close()

	stop
	end
