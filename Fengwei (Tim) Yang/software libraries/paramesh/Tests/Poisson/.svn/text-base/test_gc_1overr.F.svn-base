	subroutine test_gc_1overr(mype)
	implicit none
!
! This routine sets 1/r boundary values on external faces on the
! first block only.


#include "physicaldata.fh"
	include 'tree.fh'

	integer	mype,l,i,j,k,ii,jj,kk
	real	x,y,z,xi,yi,zi,xposn,yposn,zposn
	real	dist,potnl,dx,dy,dz

	call shmem_barrier_all()

	if(lnblocks.gt.0) then
	do l=1,lnblocks

! translate indexes into coordinates at each location
		dx = size(1,l)/real((nxb/2)*2)
		dy = size(2,l)/real((nyb/2)*2)
		dz = size(3,l)/real((nzb/2)*2)
		x = coord(1,l)
		y = coord(2,l)
		z = coord(3,l)

		xposn = .5
		yposn = .5
		zposn = .5

! face 1
		i=il_bnd
		ii = i-nguard
		xi = x - dx*(real(nxb-1)*.5-real(ii-1))
		do k=kl_bnd,ku_bnd
		kk = k-nguard
		zi = z - dz*(real(nzb-1)*.5-real(kk-1))
		do j=jl_bnd,ju_bnd
		jj = j-nguard
		yi = y - dy*(real(nyb-1)*.5-real(jj-1))
		dist = 1.e-10+sqrt((xi-xposn)*(xi-xposn)+
     .			(yi-yposn)*(yi-yposn)+(zi-zposn)*(zi-zposn))
		potnl = -1./dist
		unk(1,i,j,k,l)=potnl
		enddo
		enddo

! face 2
		i=iu_bnd
		ii = i-nguard
		xi = x - dx*(real(nxb-1)*.5-real(ii-1))
		do k=kl_bnd,ku_bnd
		kk = k-nguard
		zi = z - dz*(real(nzb-1)*.5-real(kk-1))
		do j=jl_bnd,ju_bnd
		jj = j-nguard
		yi = y - dy*(real(nyb-1)*.5-real(jj-1))
		dist = 1.e-10+sqrt((xi-xposn)*(xi-xposn)+
     .			(yi-yposn)*(yi-yposn)+(zi-zposn)*(zi-zposn))
		potnl = -1./dist
		unk(1,i,j,k,l)=potnl
		enddo
		enddo

! face 3 
		j=jl_bnd
		jj = j-nguard
		yi = y - dy*(real(nyb-1)*.5-real(jj-1))
		do k=kl_bnd,ku_bnd
		kk = k-nguard
		zi = z - dz*(real(nzb-1)*.5-real(kk-1))
		do i=il_bnd,iu_bnd
		ii = i-nguard
		xi = x - dx*(real(nxb-1)*.5-real(ii-1))
		dist = 1.e-10+sqrt((xi-xposn)*(xi-xposn)+
     .			(yi-yposn)*(yi-yposn)+(zi-zposn)*(zi-zposn))
		potnl = -1./dist
		unk(1,i,j,k,l)=potnl
		enddo
		enddo

! face 4 
		j=ju_bnd
		jj = j-nguard
		yi = y - dy*(real(nyb-1)*.5-real(jj-1))
		do k=kl_bnd,ku_bnd
		kk = k-nguard
		zi = z - dz*(real(nzb-1)*.5-real(kk-1))
		do i=il_bnd,iu_bnd
		ii = i-nguard
		xi = x - dx*(real(nxb-1)*.5-real(ii-1))
		dist = 1.e-10+sqrt((xi-xposn)*(xi-xposn)+
     .			(yi-yposn)*(yi-yposn)+(zi-zposn)*(zi-zposn))
		potnl = -1./dist
		unk(1,i,j,k,l)=potnl
		enddo
		enddo

! face 5 
		k=kl_bnd
		kk = k-nguard
		zi = z - dz*(real(nzb-1)*.5-real(kk-1))
		do j=jl_bnd,ju_bnd
		jj = j-nguard
		yi = y - dy*(real(nyb-1)*.5-real(jj-1))
		do i=il_bnd,iu_bnd
		ii = i-nguard
		xi = x - dx*(real(nxb-1)*.5-real(ii-1))
		dist = 1.e-10+sqrt((xi-xposn)*(xi-xposn)+
     .			(yi-yposn)*(yi-yposn)+(zi-zposn)*(zi-zposn))
		potnl = -1./dist
		unk(1,i,j,k,l)=potnl
		enddo
		enddo

! face 6 
		k=ku_bnd
		kk = k-nguard
		zi = z - dz*(real(nzb-1)*.5-real(kk-1))
		do j=jl_bnd,ju_bnd
		jj = j-nguard
		yi = y - dy*(real(nyb-1)*.5-real(jj-1))
		do i=il_bnd,iu_bnd
		ii = i-nguard
		xi = x - dx*(real(nxb-1)*.5-real(ii-1))
		dist = 1.e-10+sqrt((xi-xposn)*(xi-xposn)+
     .			(yi-yposn)*(yi-yposn)+(zi-zposn)*(zi-zposn))
		potnl = -1./dist
		unk(1,i,j,k,l)=potnl
		enddo
		enddo

	enddo
	endif

	call shmem_barrier_all()

	return
	end
