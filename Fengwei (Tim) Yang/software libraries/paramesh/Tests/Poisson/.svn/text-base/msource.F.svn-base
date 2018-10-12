	subroutine msource(mype,source)
!
! This routine sets the density that results when a point source
! of fixed mass is located on a grid.

! Variables:
!	totmass		the total mass of the point source
!	xposn		the x coordinate of the point mass
!	yposn		the y coordinate of the point mass
!	zposn		the z coordinate of the point mass

#include "physicaldata.fh"
	include 'tree.fh'

	real	source(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                 kl_bnd1:ku_bnd1,	maxblocks)



	totmass = 1.
	xposn = .5
	yposn = .5
	zposn = .5

	if(lnblocks.gt.0) then
	do l=1,lnblocks

	dx = size(1,l)/(2*(nxb/2))
	dy = size(2,l)/(2*(nyb/2))
	dz = size(3,l)/(2*(nzb/2))
	vol = dx*dy*dz
!	write(*,*) 'proc ',mype,' block ',l,' dxyz ',dx,dy,dz,vol

	if(nodetype(l).eq.1) then

		density = (totmass/vol)*4.*3.14159265359

		do k=1+nguard,nzb+nguard
			kk = k - nguard
			zl = coord(3,l)-dz*(real(nzb/2)+.5-real(kk-1))
			zu = zl+dz
			zinside = (zposn-zl)*(zposn-zu)
			zf = 0.
			if(zinside.lt.0.) zf = 1.
		do j=1+nguard,nyb+nguard
			jj = j - nguard
			yl = coord(2,l)-dy*(real(nyb/2)+.5-real(jj-1))
			yu = yl+dy
			yinside = (yposn-yl)*(yposn-yu)
			yf = 0.
			if(yinside.lt.0.) yf = 1.
		do i=1+nguard,nxb+nguard
			ii = i - nguard
			xl = coord(1,l)-dx*(real(nxb/2)+.5-real(ii-1))
			xu = xl+dx
			xinside = (xposn-xl)*(xposn-xu)
			xf = 0.
			if(xinside.lt.0.) xf = 1.
			source(i,j,k,l) = density*xf*yf*zf
!			source(i,j,k,l) = sin(3.14159265359*(xl+xu))
		enddo
		enddo
		enddo

	endif

	enddo
	endif

	call shmem_barrier_all()

	return
	end
