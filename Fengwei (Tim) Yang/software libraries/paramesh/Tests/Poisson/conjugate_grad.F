      subroutine conjugate_grad (mype,loop,iopt,nlayers,source)

! A conjugate-gradient solver for Poisson's equation, written
! for a test of the amr code on the T3D, based on a routine
! written by Kevin Olson.

#include "physicaldata.fh"
	include 'tree.fh'

        include         'amr_shmem.fh'
        include         'shmem_reduce.fh'


c Solves grad^2 u = q using conjugant gradient assuming u is known on the
c boundries of the domain and that the size of all cells are the same.
c `h' is the cell size.

!	parameter(nstepsmax=100)
	real	h,hinv
	real	res(0:nxb+1+nguard,0:nyb+1+nguard,0:nzb+1+nguard)
	real	s(0:nxb+1+nguard,0:nyb+1+nguard,0:nzb+1+nguard)
	real	beta
	real	omega(maxblocks),eps(maxblocks)
	real	q_sq
	real	error(maxblocks),error_o(maxblocks),errorl
	integer	i,j,k,l,nsteps,mype
	integer	iopt,nlayers
	logical	start
        logical lcc,lfc,l_srl_only,ldiag
        integer icoord

	common/reduct/q_sq_tot,q_sq_sum
	real	q_sq_tot,q_sq_sum
!	real*8	q_sq_tot,q_sq_sum

	real	source(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                 kl_bnd1:ku_bnd1,	maxblocks)


	integer	ierror,ierrorg
	save	ierror,ierrorg


!
        lcc = .true.
        lfc = .false.
        l_srl_only =  .false.
        ldiag = .false.
        icoord = 0


c compute stopping error eps

	q_sq_tot=0.
	if(lnblocks.gt.0) then
	do l=1,lnblocks
	if(nodetype(l).eq.1) then

	q_sq = 0.
	do k = 1+nguard,nzb+nguard
         do j = 1+nguard,nyb+nguard
            do i = 1+nguard,nxb+nguard

               q_sq = q_sq + (source(i,j,k,l)**2)

            end do
         end do
	end do
	if(mype.eq.0) write(*,*) 'l ',l,' q_sq ',q_sq,' size ',
     .						size(1,l)

	q_sq_tot=q_sq_tot+q_sq*size(1,l)*size(2,l)*size(3,l)

	eps(l) = 1.e-5*q_sq

	endif
	enddo
	endif
	if(mype.eq.0) write(*,*) 'q_sq ',q_sq,' q_sq_tot ',q_sq_tot
        call shmem_barrier_all()
        call comm_real_sum_to_all(q_sq_sum,q_sq_tot)
        call shmem_barrier_all()

	if(lnblocks.gt.0) then
	do l=1,lnblocks
	if(nodetype(l).eq.1) then
		epsnorm=1.e-5*q_sq_sum/(size(1,l)*size(2,l)*size(3,l))
		eps(l)=max(eps(l),epsnorm)
	endif
	enddo
	endif

	ierrorg=10
	nstepmaxg=100
	nstepsg=0

	nstepsmax=100

! loop to allow block solutions to adjust to changing guardcell values
	do while (ierrorg.gt.0.and.nstepsg.lt.nstepmaxg)
	nstepsg=nstepsg+1

! exchange guardcell data between blocks
#ifndef NO_PERMANENT_GUARDCELLS
	call amr_guardcell(mype,iopt,nlayers)
        call shmem_barrier_all()
#endif

#ifdef NO_PERMANENT_GUARDCELLS
        call amr_1blk_copy_soln
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over all active leaf blocks on this processor
	if(lnblocks.gt.0) then
	do l=1,lnblocks
	if(nodetype(l).eq.1) then

	h=size(1,l)/(2*(nxb/2))
	hinv = 1./h**2

	error(l) = max(1.e10,2.*eps(l))


	start = .TRUE.
	nsteps = 0

#ifdef NO_PERMANENT_GUARDCELLS
        idest = 1
        call amr_1blk_guardcell(mype,iopt,nlayers,l,mype,
     .                          lcc,lfc,l_srl_only,icoord,ldiag)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c start conj. grad. iterations on block l
	do while (error(l).gt.eps(l).and.nsteps.lt.nstepsmax)

	nsteps = nsteps + 1
!	if(mype.eq.0) write(*,*) 'proc ',mype,' block ',l,
!     .					' Conj-grad step ',nsteps


c compute residuals = (grad^2 u)  - q via 7 point finite diff. formula
         res(:,:,:) = 0.
         s(:,:,:) = 0.

         idest = 1
         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard
                  
#ifdef NO_PERMANENT_GUARDCELLS
                  res(i,j,k) = hinv*( unk1(1,i+1,j,k,idest) + 
     &             unk1(1,i-1,j,k,idest) + unk1(1,i,j+1,k,idest) +
     &             unk1(1,i,j-1,k,idest) + unk1(1,i,j,k+1,idest) +
     &             unk1(1,i,j,k-1,idest) - (6.*unk1(1,i,j,k,idest)) )
     &             - source(i,j,k,l)
#else
                  res(i,j,k) = hinv*( unk(1,i+1,j,k,l) + 
     &                 unk(1,i-1,j,k,l) + unk(1,i,j+1,k,l) +
     &                 unk(1,i,j-1,k,l) + unk(1,i,j,k+1,l) +
     &                 unk(1,i,j,k-1,l) - (6.*unk(1,i,j,k,l)) )
     &                 - source(i,j,k,l)
#endif
               end do
            end do
         end do

c compute error

         error_o(l) = error(l)
         error(l) = 0.
         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard

                  error(l) = error(l) + (res(i,j,k)**2)

               end do
            end do
         end do

	if(error(l).gt.eps(l)) then
         print *,'proc ',mype,' block ',l,' **error ',error(l),
     .' eps ',eps(l),nsteps
	else
         print *,'proc ',mype,' block ',l,'   error ',error(l),
     .' eps ',eps(l),nsteps
	endif
         beta = error(l)/error_o(l)

c search direction

         if (start) then

            do k = 0+nguard,nzb+1+nguard
               do j = 0+nguard,nyb+1+nguard
                  do i = 0+nguard,nxb+1+nguard

                     s(i,j,k) = -res(i,j,k)
                     
                  end do
               end do
            end do

         else

            do k = 0+nguard,nzb+1+nguard
               do j = 0+nguard,nyb+1+nguard
                  do i = 0+nguard,nxb+1+nguard

                     s(i,j,k) = -res(i,j,k) + 
     .					(beta*s(i,j,k))
                     
                  end do
               end do
            end do

         end if

	if(start) start=.false.

c distance to minimum, res is used as a temporary here and = grad^2 s

         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard
                  
                  res(i,j,k) = hinv*( s(i+1,j,k) + 
     &                 s(i-1,j,k) + s(i,j+1,k) +
     &                 s(i,j-1,k) + s(i,j,k+1) +
     &                 s(i,j,k-1) - (6.*s(i,j,k)) )

               end do
            end do
         end do

         omega(l) = 0.
         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard

                  omega(l) = omega(l) + (s(i,j,k)*res(i,j,k))

               end do
            end do
         end do

         omega(l) = error(l)/omega(l)

c next iterate

         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard

#ifdef NO_PERMANENT_GUARDCELLS
                  unk1(1,i,j,k,idest) = unk1(1,i,j,k,idest) + 
     .					(omega(l)*s(i,j,k))
#else
                  unk(1,i,j,k,l) = unk(1,i,j,k,l) + 
     .					(omega(l)*s(i,j,k))
#endif

               end do
            end do
         end do


! end of while loop
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef NO_PERMANENT_GUARDCELLS
         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard
                  unk(1,i-nguard,j-nguard,k-nguard,l) = 
     .                           unk1(1,i,j,k,idest) 
               end do
            end do
         end do
#endif


! end of loop over blocks on this processor
	endif
	enddo
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! swap guardcells between blocks
        call shmem_barrier_all()
#ifndef NO_PERMANENT_GUARDCELLS
	call amr_guardcell(mype,iopt,nlayers)
        call shmem_barrier_all()
#endif
#ifdef NO_PERMANENT_GUARDCELLS
        call amr_1blk_copy_soln
#endif


! recompute residuals and errors following guard cell swap
        ierror = 0
	if(lnblocks.gt.0) then
	do l=1,lnblocks
	if(nodetype(l).eq.1) then

#ifdef NO_PERMANENT_GUARDCELLS
        call amr_1blk_guardcell(mype,iopt,nlayers,l,mype,
     .                          lcc,lfc,l_srl_only,icoord,ldiag)
#endif

	h=size(1,l)/(2*(nxb/2))
	hinv = 1./h**2

         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard

#ifdef NO_PERMANENT_GUARDCELLS
                res(i,j,k) = hinv*( unk1(1,i+1,j,k,idest) +
     &            unk1(1,i-1,j,k,idest) + unk1(1,i,j+1,k,idest) +
     &            unk1(1,i,j-1,k,idest) + unk1(1,i,j,k+1,idest) +
     &            unk1(1,i,j,k-1,idest) - (6.*unk1(1,i,j,k,idest)) )
     &            - source(i,j,k,l)
#else
                res(i,j,k) = hinv*( unk(1,i+1,j,k,l) +
     &            unk(1,i-1,j,k,l) + unk(1,i,j+1,k,l) +
     &            unk(1,i,j-1,k,l) + unk(1,i,j,k+1,l) +
     &            unk(1,i,j,k-1,l) - (6.*unk(1,i,j,k,l)) )
     &            - source(i,j,k,l)
#endif
               end do
            end do
         end do

c compute error

         errorl = 0.
         do k = 1+nguard,nzb+nguard
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard

                  errorl = errorl + (res(i,j,k)**2)

               end do
            end do
         end do

	if(errorl.gt.eps(l)) then
		write(*,*) 'not converged: pe ',mype,' blk ',l,
     .		' errorl/eps = ',
     .		errorl/eps(l),' eps = ',eps(l)
		ierror=1
	endif

	endif
	enddo
	endif
        call shmem_barrier_all()
        call comm_int_sum_to_all(ierrorg,ierror)
	if(mype.eq.0) write(*,*) 'nstepsg= ',nstepsg,
     .				' ierrorg = ',ierrorg
        call shmem_barrier_all()

	enddo


        call shmem_barrier_all()

	return
	end 
