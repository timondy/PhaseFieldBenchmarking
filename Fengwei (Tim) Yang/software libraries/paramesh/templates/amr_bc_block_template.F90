!!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_bc_block(jface,ibc,iopt,l,mype)





!------------------------------------------------------------------------
! This is a template. User editing is required before this routine 
! can be used.

! This routine sets the guard cell elements of the solution arrays 
! or workspace array on face jface of block l, using a boundary
! condition algorithm which the user must insert. Places where user
! editing is required are clearly indicated.
!
! Note, if your algorithm uses diagonal elements (see the users guide for
! a definition of diagonal elements) then you are responsible for setting
! the data values in these diagonal guard cells on the corners and edges
! of external non-periodic boundaries.
!
!
! Written :     Peter MacNeice          September 1997
!------------------------------------------------------------------------
!
! Arguments:
!       jface           number designating selected block face
!       ibc             the boundary condition flag. This should be
!                       set equal to neigh(1,jface,l) by the calling
!                       routine. Its value can be used to vary the
!                       boundary conditions applied on different boundaries.
!       iopt            a switch to control which data is updated
!                       iopt=1 will use 'unk' and 'facevarx(y)(z)'
!                       iopt=2 will use 'work'
!       l               number designating selected block 
!       mype            local processor number
!
!------------------------------------



! include file to define physical qualities of the model and mesh
      use physicaldata
      use workspace

! include file defining the tree
      use tree

      include 'mpif.h'

      integer jface,iopt,l,mype,ibc

!------------------------------------



      if(iopt.eq.1) then


! The loop bounds set here include the diagonal elements. If these are
! not required then change the bounds set in the next 6 lines to 
! ilbnd=1+nguard, iubnd=nxb+nguard, and similarily for y and z.
                ilbnd = 1
                iubnd = nxb+2*nguard
                jlbnd = 1
                jubnd = nyb+2*nguard*k2d
                klbnd = 1
                kubnd = nzb+2*nguard*k3d

                if(jface.eq.1) then
                   ilbnd = 1
                   iubnd = nguard
                elseif(jface.eq.2) then
                   ilbnd = 1+nxb+nguard
                   iubnd = nxb+2*nguard
                elseif(jface.eq.3) then
                   jlbnd = 1
                   jubnd = nguard
                elseif(jface.eq.4) then
                   jlbnd = 1+nyb+nguard
                   jubnd = nyb+2*nguard
                elseif(jface.eq.5) then
                   klbnd = 1
                   kubnd = nguard
                elseif(jface.eq.6) then
                   klbnd = 1+nzb+nguard
                   kubnd = nzb+2*nguard
                endif



! Limit the application of this routine to leaf blocks
! or the parents of leaf blocks.
      if(nodetype(l).eq.1.or.nodetype(l).eq.2) then


! Loop over the cells of the selected face setting unk.
       do k=klbnd,kubnd
         do j=jlbnd,jubnd
           do i=ilbnd,iubnd
             do ivar=1,nvar
               unk(ivar,i,j,k,l) = ??????                       !<<<< USER EDIT
             enddo
           enddo
         enddo
       enddo


! Now repeat for the facevar arrays if necessary.

       if(nfacevar.gt.0) then

       ione = 1
       jone = 0
       kone = 0
       if(jface.eq.1) ione = 0
       do k=klbnd,kubnd
         do j=jlbnd,jubnd
           do i=ilbnd,iubnd+ione
             do ivar=1,nbndvar
               facevarx(ivar,i,j,k,l)=???                       !<<<< USER EDIT
             enddo
           enddo
         enddo
       enddo

       ione = 0
       jone = 1
       kone = 0
       if(jface.eq.3) jone = 0
       do k=klbnd,kubnd
         do i=ilbnd,iubnd
           do j=jlbnd,jubnd+jone
             do ivar=1,nbndvar
               facevary(ivar,i,j,k,l)=???                       !<<<< USER EDIT
             enddo
           enddo
         enddo
       enddo

       if(ndim.eq.3) then
       ione = 0
       jone = 0
       kone = 1
       if(jface.eq.5) kone = 0
       do j=jlbnd,jubnd
         do i=ilbnd,iubnd
           do k=klbnd,kubnd+kone
             do ivar=1,nbndvar
               facevarz(ivar,i,j,k,l)=???                       !<<<< USER EDIT
             enddo
           enddo
         enddo
       enddo
       endif

      endif

      endif



      elseif(iopt.eq.2) then


! The loop bounds set here include the diagonal elements. If these are
! not required then change the bounds set in the next 6 lines to 
! ilbnd=1+nguard_work, iubnd=nxb+nguard_work, and similarily for y and z.
                ilbnd = 1
                iubnd = nxb+2*nguard_work
                jlbnd = 1
                jubnd = nyb+2*nguard_work*k2d
                klbnd = 1
                kubnd = nzb+2*nguard_work*k3d


                if(jface.eq.1) then
                   ilbnd = 1
                   iubnd = nguard_work
                elseif(jface.eq.2) then
                   ilbnd = 1+nxb+nguard_work
                   iubnd = nxb+2*nguard_work
                elseif(jface.eq.3) then
                   jlbnd = 1
                   jubnd = nguard_work
                elseif(jface.eq.4) then
                   jlbnd = 1+nyb+nguard_work
                   jubnd = nyb+2*nguard_work
                elseif(jface.eq.5) then
                   klbnd = 1
                   kubnd = nguard_work
                elseif(jface.eq.6) then
                   klbnd = 1+nzb+nguard_work
                   kubnd = nzb+2*nguard_work
                endif


! Limit the application of this routine to leaf blocks
! or the parents of leaf blocks.
       if(nodetype(l).eq.1.or.nodetype(l).eq.2) then

! Loop over the cells of the selected face setting work.
       do k=klbnd,kubnd
         do j=jlbnd,jubnd
           do i=ilbnd,iubnd
             work(i,j,k,l) = ??????                             !<<<< USER EDIT
           enddo
         enddo
       enddo

       endif


      endif


      return
      end
