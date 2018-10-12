!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_bc_block(jface,ibc,iopt,l,mype)





! This routine sets the guard cell elements of the solution arrays 
! to a linear function of the local coordinate on face jface of block l.

      integer jface,iopt,l,mype,ibc


! include file to define physical qualities of the model and mesh
      use physicaldata
      use workspace

! include file defining the tree
       use tree

! Only necessary for programs in ./Tests
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "test_defs.fh"

! default values used if ndim < 3.
      dz = 0.
      z0 = 0.


      if(iopt.eq.1) then


                ilbnd = 1
                iubnd = nxb+2*nguard
                jlbnd = 1
                jubnd = nyb+2*nguard*k2d
                klbnd = 1
                kubnd = nzb+2*nguard*k3d

                if(jface.eq.1) iubnd = nguard
                if(jface.eq.2) ilbnd = 1+nxb+nguard
                if(jface.eq.3) jubnd = nguard
                if(jface.eq.4) jlbnd = 1+nyb+nguard
                if(jface.eq.5) kubnd = nguard
                if(jface.eq.6) klbnd = 1+nzb+nguard



! set the solution array to be the grid points x,y or z coordinates
      if(nodetype(l).eq.1.or.nodetype(l).eq.2) then

      if(mod(nxb,2).eq.0) then
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
       dy = bsize(2,l)/real(nyb)
       dx = bsize(1,l)/real(nxb)
       do k=klbnd,kubnd
       do j=jlbnd,jubnd
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard)
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard)
       value = ax*x0 + ay*y0 + az*z0
       do ivar=1,nvar
       unk(ivar,i,j,k,l) = value*real(ivar)
       enddo
       enddo
       enddo
       enddo


      else
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       dy = bsize(2,l)/real(nyb-1)
       dx = bsize(1,l)/real(nxb-1)
       do k=klbnd,kubnd
       do j=jlbnd,jubnd
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*bsize(1,l)-dx+dx*real(i-nguard)
       y0 = coord(2,l)-.5*bsize(2,l)-dy+dy*real(j-nguard)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*bsize(3,l)-dz+dz*real(k-nguard)
       value = ax*x0 + ay*y0 + az*z0
       do ivar=1,nvar
       unk(ivar,i,j,k,l) = value*real(ivar)
       enddo
       enddo
       enddo
       enddo

      endif

      endif

! set up data in facevarx etc

        if(mod(nxb,2).eq.0) then

                ilbnd = 1
                iubnd = nxb+2*nguard
                jlbnd = 1
                jubnd = nyb+2*nguard
                klbnd = 1
                kubnd = nzb+2*nguard*k3d

                if(jface.eq.1) iubnd = nguard
                if(jface.eq.2) ilbnd = 1+nxb+nguard
                if(jface.eq.3) jubnd = nguard
                if(jface.eq.4) jlbnd = 1+nyb+nguard
                if(jface.eq.5) kubnd = nguard
                if(jface.eq.6) klbnd = 1+nzb+nguard



      if(nodetype(l).eq.1.or.nodetype(l).eq.2) then

              if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
              dy = bsize(2,l)/real(nyb)
              dx = bsize(1,l)/real(nxb)
              if(mod(nxb,2).eq.1) then
                      if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                      dy = bsize(2,l)/real(nyb-1)
                      dx = bsize(1,l)/real(nxb-1)
              endif

#ifdef FACEX
       ione = 1
       jone = 0
       kone = 0
       if(jface.eq.1) ione = 0
       do k=klbnd,kubnd
              if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
              zk = z0 + dz*real(k-nguard)
       do j=jlbnd,jubnd
                      y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                      yj = y0 + dy*real(j-nguard)
       do i=ilbnd,iubnd+ione
                       x0 = coord(1,l)-.5*bsize(1,l)-dx
                       xi = x0 + dx*real(i-nguard)
                       value = ax*xi + ay*yj + az*zk
                       do ivar=1,nbndvar
                              facevarx(ivar,i,j,k,l)=value*real(ivar)
                              enddo
                              enddo
                      enddo
              enddo
#endif

#ifdef FACEY
       ione = 0
       jone = 1
       kone = 0
       if(jface.eq.3) jone = 0
       do k=klbnd,kubnd
              if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
              zk = z0 + dz*real(k-nguard)
       do i=ilbnd,iubnd
                      x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                      xi = x0 + dx*real(i-nguard)
       do j=jlbnd,jubnd+jone
                       y0 = coord(2,l)-.5*bsize(2,l)-dy
                       yj = y0 + dy*real(j-nguard)
                       value = ax*xi + ay*yj + az*zk
                       do ivar=1,nbndvar
                              facevary(ivar,i,j,k,l)=value*real(ivar)
                              enddo
                              enddo
                      enddo
              enddo
#endif

#ifdef FACEZ
       if (ndim == 3) then
       ione = 0
       jone = 0
       kone = 1
       if(jface.eq.5) kone = 0
       do j=jlbnd,jubnd
              y0 = coord(2,l)-.5*(bsize(2,l)+dy)
              yj = y0 + dy*real(j-nguard)
       do i=ilbnd,iubnd
                      x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                      xi = x0 + dx*real(i-nguard)
       do k=klbnd,kubnd+kone
                       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*bsize(3,l)-dz
                       zk = z0 + dz*real(k-nguard)
                       value = ax*xi + ay*yj + az*zk
                       do ivar=1,nbndvar
                              facevarz(ivar,i,j,k,l)=value*real(ivar)
                              enddo
                              enddo
                      enddo
              enddo
        endif
#endif

        endif

      endif



      elseif(iopt.ge.2) then



                ilbnd = 1
                iubnd = nxb+2*nguard_work
                jlbnd = 1
                jubnd = nyb+2*nguard_work
                klbnd = 1
                kubnd = nzb+2*nguard_work*k3d

                if(jface.eq.1) iubnd = nguard_work
                if(jface.eq.2) ilbnd = 1+nxb+nguard_work
                if(jface.eq.3) jubnd = nguard_work
                if(jface.eq.4) jlbnd = 1+nyb+nguard_work
                if(jface.eq.5) kubnd = nguard_work
                if(jface.eq.6) klbnd = 1+nzb+nguard_work


! set the work array to be the grid points x,y or z coordinates
      if(nodetype(l).eq.1.or.nodetype(l).eq.2) then

      if(mod(nxb,2).eq.0) then
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
       dy = bsize(2,l)/real(nyb)
       dx = bsize(1,l)/real(nxb)
       do k=klbnd,kubnd
       do j=jlbnd,jubnd
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard_work)
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard_work)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard_work)
       value = ax*x0 + ay*y0 + az*z0
       work(i,j,k,l,iopt-1) = value
       enddo
       enddo
       enddo

      else
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       dy = bsize(2,l)/real(nyb-1)
       dx = bsize(1,l)/real(nxb-1)
       do k=klbnd,kubnd
       do j=jlbnd,jubnd
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*bsize(1,l)-dx+dx*real(i-nguard_work)
       y0 = coord(2,l)-.5*bsize(2,l)-dy+dy*real(j-nguard_work)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*bsize(3,l)-dz+dz*real(k-nguard_work)
       value = ax*x0 + ay*y0 + az*z0
       work(i,j,k,l,iopt-1) = value
       enddo
       enddo
       enddo

      endif

      endif


      endif


      return
      end
