      subroutine amr_1blk_bcset(mype,ibc,lb,pe,
     .    idest,iopt,id,jd,kd,ilays,jlays,klays,ip1,jp1,kp1)


! $RCSfile $
! $Revision $
! $Date $


      implicit real(a-h,o-z)


      integer mype,ibc,lb,pe
      integer idest,iopt,id,jd,kd,ilays,jlays,klays,ip1,jp1,kp1

!------------------------------------------------------------------------
!
! This routine sets guardcell values at external boundaries in the case
! where a single block is having its guardcells filled.
!
! It can be assumed in writing this routine, that all guardcells for this
! block which are not across an external boundary have already been
! properly filled.
!
!
! Arguments:
!      mype             local processor
!      ibc              the integer specifying the particular boundary
!                        condition to be imposed
!      lb               block number of selected block
!      pe               processor on which block lb is located
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      id               lower limit of index range of points in x direction
!      jd               lower limit of index range of points in y direction
!      kd               lower limit of index range of points in z direction
!      ilay             no. of mesh points in x direction to be set
!      jlay             no. of mesh points in y direction to be set
!      klay             no. of mesh points in z direction to be set
!      ip1              offset added to index range defined by (id,ilay)
!                        0 if guardcells are at lower end of i index
!                        1 if guardcells are at upper end of i index
!      jp1              offset added to index range defined by (jd,jlay)
!                        0 if guardcells are at lower end of j index
!                        1 if guardcells are at upper end of j index
!      kp1              offset added to index range defined by (kd,klay)
!                        0 if guardcells are at lower end of k index
!                        1 if guardcells are at upper end of k index
!
!
!
! Written :     Peter MacNeice          August 1998
!------------------------------------------------------------------------

#include "physicaldata.fh"
      include 'tree.fh'
      include 'workspace.fh'

      include 'amr_shmem.fh'



!-------------------------

!
! Adjust index ranges
      il = ilays-1
      jl = (jlays-1)*k2d
      kl = (klays-1)*k3d


! Which boundary condition has been specified?
      if(ibc.le.-20) then

!
! Do cell-face-centered data
        if(nfacevar.gt.0.and.iopt.eq.1) then


       do k=kd,kd+kl
         do j=jd,jd+jl
           do i=id,id+il+ip1
             do ivar=1,nbndvar
               facevarx1(ivar,i,j,k,:)=0.
             enddo
           enddo
         enddo
       enddo


       if(ndim.ge.2) then
       do k=kd,kd+kl
         do j=jd,jd+jl+jp1
           do i=id,id+il
             do ivar=1,nbndvar
               facevary1(ivar,i,j,k,:)=0.
             enddo
           enddo
         enddo
       enddo
       endif

       if(ndim.eq.3) then
       do k=kd,kd+kl+kp1
         do j=jd,jd+jl
           do i=id,id+il
             do ivar=1,nbndvar
               facevarz1(ivar,i,j,k,:)=0.
             enddo
           enddo
         enddo
       enddo
       endif

        endif                            ! end of nfacevar if test


!
! Now do cell centered data

        if(iopt.eq.1) then


! translate indexes into coordinates at each location
	  dx = size(1,lb)/real((nxb/2)*2)
	  dy = size(2,lb)/real((nyb/2)*2)
	  dz = size(3,lb)/real((nzb/2)*2)
	  x = coord(1,lb)
	  y = coord(2,lb)
	  z = coord(3,lb)

! location of mass source
	  xposn = .5
	  yposn = .5
	  zposn = .5


          l = idest
          do k = kd,kd+kl
	    kk = k-nguard
            zi = z - dz*(real(nzb-1)*.5-real(kk-1))
            do j = jd,jd+jl
              jj = j-nguard
              yi = y - dy*(real(nyb-1)*.5-real(jj-1))
              do i = id,id+il
		ii = i-nguard
		xi = x - dx*(real(nxb-1)*.5-real(ii-1))

		dist = 1.e-10+sqrt((xi-xposn)*(xi-xposn)+
     .			(yi-yposn)*(yi-yposn)+(zi-zposn)*(zi-zposn))
		potnl = -1./dist
		unk1(1,i,j,k,l)=potnl

              enddo
            enddo
          enddo

        elseif(iopt.ge.2) then


          do k = kd,kd+kl
          do j = jd,jd+jl
          do i = id,id+il

            work1(i,j,k,iopt-1) = 0.

          enddo
          enddo
          enddo


        endif
!-------------------------



       endif                            ! end of test of bc flag

      return
      end
