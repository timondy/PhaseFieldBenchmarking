!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine check_data(mype,pe,lblock,ia,ib,ja,jb,ka,kb)




! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata
      use workspace

! include file defining the tree
      use tree

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "test_defs.fh"
      include 'mpif.h'

      integer :: pe,lblock,ia,ib,ja,jb,ka,kb
      integer :: ierr

      ax = 1.
      ay = 10.
      az = 100.

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if(mype.eq.pe) then

      l=lblock

      ilbnd=ia
      iubnd=ib
      jlbnd=ja
      jubnd=jb
      klbnd=ka
      kubnd=kb

      dz = 0.
      z0 = 0.

      if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
      dy = bsize(2,l)/real(nyb)
      dx = bsize(1,l)/real(nxb)
      if(mod(nxb,2).eq.1) then
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       dy = bsize(2,l)/real(nyb-1)
       dx = bsize(1,l)/real(nxb-1)
      endif

      do k=klbnd,kubnd
      if(ndim.eq.3) then
      z0 = coord(3,l)-.5*(bsize(3,l)+dz)
      if(mod(nxb,2).eq.1) z0 = coord(3,l)-(.5*bsize(3,l)+dz)
      endif
      zk = z0 + dz*real(k-nguard)
       do j=jlbnd,jubnd
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)
       if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       yj = y0 + dy*real(j-nguard)
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)
       if(mod(nxb,2).eq.1) x0 =  & 
     &       coord(1,l)-(.5*bsize(1,l)+dx)
       xi = x0 + dx*real(i-nguard)
       value = ax*xi+ay*yj+az*zk
#ifdef TESTALLDIR
                        do ivar=1,nvar
       vv = value*real(ivar)
                        if(vv.ne.unk(ivar,i,j,k,l))  write(*,998) & 
     &                   mype,l,ivar,i,j,k,unk(ivar,i,j,k,l),vv
                        enddo
#endif
#ifdef TESTZDIR
                        if (ndim == 3) then
                        do ivar=1,nvar
       vv = zk*real(ivar)
                        if(vv.ne.unk(ivar,i,j,k,l)) write(*,998) & 
     &                   mype,l,ivar,i,j,k,unk(ivar,i,j,k,l),vv
                        enddo
                        end if
#endif
#ifdef TESTYDIR
                        do ivar=1,nvar
       vv = yj*real(ivar)
                        if(vv.ne.unk(ivar,i,j,k,l)) write(*,998) & 
     &                   mype,l,ivar,i,j,k,unk(ivar,i,j,k,l),vv
                        enddo
#endif
#ifdef TESTXDIR
                        do ivar=1,nvar
       vv = xi*real(ivar)
                        if(vv.ne.unk(ivar,i,j,k,l)) write(*,998) & 
     &                   mype,l,ivar,i,j,k,unk(ivar,i,j,k,l),vv
                        enddo
#endif
       enddo
       enddo
      enddo

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


998      format('u:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)


      return
      end
