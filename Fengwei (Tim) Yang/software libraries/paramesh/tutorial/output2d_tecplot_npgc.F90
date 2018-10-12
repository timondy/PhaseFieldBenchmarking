!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

        subroutine output2d_tecplot(lrefine_min,lrefine_max,time,dt)

! 
! This routine captures the required output from all processors
! and outputs it.
!
!
! Written:       Peter MacNeice
! January 1998
!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use physicaldata

! include file defining the tree 
      use tree

! include file required for shmem library.
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"
      include 'shmem_reduce.fh'

      integer lrefine_min,lrefine_max
      real    time,dt

      integer nguard0
      parameter(nguard0 = nguard*npgs)


      integer nprocs,mype
      integer shmem_my_pe,shmem_n_pes

      real	cunk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real      cfacevarx(nfacevar, & 
     &          il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real	ccoord(3),csize(3)
      real	cbnd_box(2,3)
       

	integer	cnodetype,nblocks,level,iproc,crefine,lout

	save mype,nblocks,cunk,cnodetype,ccoord,csize,crefine
	save cbnd_box,cfacevarx

	character(len=3) filenumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call shmem_barrier_all()

        mype = shmem_my_pe()
        nprocs = shmem_n_pes()

        write(*,*) 'lrefine_min,lrefine_max ', & 
     &              lrefine_min,lrefine_max

	if(mype.eq.0) then

	do level=lrefine_min,lrefine_max

	write(filenumber,"(I3.3)") level
	iout = 10+level
        open(unit=iout,status='new',file='tecplot.'//filenumber)
	write(*,*) 'level ',level,' file ','tecplot.'//filenumber

	write(iout,*) 'TITLE="Tecplot Test"'
	write(iout,*) 'VARIABLES="X" "Y" "Density"'

	enddo


	do iproc=0,nprocs-1

	call shmem_integer_get(nblocks,lnblocks,1,iproc)
	write(*,*) 'proc ',iproc,' has ',nblocks,' blocks'
	if(nblocks.gt.0) then
        do lb=1,nblocks

		call shmem_integer_get(cnodetype,nodetype(lb),1, & 
     &						iproc)

       		if(cnodetype.eq.1) then
		call shmem_integer_get(crefine,lrefine(lb),1,iproc)
		call shmem_real_get(csize(1),size(1,lb),3,iproc)
		call shmem_real_get(ccoord(1),coord(1,lb),3,iproc)
		call shmem_real_get(cbnd_box(1,1),bnd_box(1,1,lb),6, & 
     &                                                      iproc)
!		call shmem_real_get(cfacevarx(1,1,1,1),
!     .                              facevarx(1,1,1,1,lb),
!     .		               	len_blockfx*nfacevar,iproc)
		call shmem_real_get(cunk(1,1,1,1),unk(1,1,1,1,lb), & 
     &		                 		len_block,iproc)
		lout = 10+crefine
		write(lout,*) 'ZONE F=POINT, I=',nxb, ', J=',nyb
                dy = csize(2)/real((nyb/2)*2)
                dx = csize(1)/real((nxb/2)*2)

                k=1
		do j=nguard0+1,nyb+2*nguard0
                yi = cbnd_box(1,2)+dy*(real(j-nguard0-1)+.5)
                do i=nguard0+1,nxb+2*nguard0
                xi = cbnd_box(1,1)+dx*(real(i-nguard0-1)+.5)
		 den =  cunk(1,i,j,k)
		write(lout,50) xi,yi,den
                enddo
                enddo

        	endif

        enddo
        endif

        enddo

50      format(3(1x,1pe13.6))

	do level=lrefine_min,lrefine_max
	iout = 10+level
        close(unit=iout)
	enddo

	endif
        call shmem_barrier_all()

      return
      end
