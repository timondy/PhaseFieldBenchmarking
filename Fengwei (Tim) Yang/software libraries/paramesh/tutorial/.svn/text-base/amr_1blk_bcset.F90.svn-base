!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_bcset(mype,ibc,lb,pe, & 
     &    idest,iopt,id,jd,kd,ilays,jlays,klays,ip1,jp1,kp1)


! $RCSfile $
! $Revision $
! $Date $

      use physicaldata
      use tree
      use workspace

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

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"


      real ccoord(3),csize(3)
      save    ccoord,csize



!---------------------------------------------------------------------------
! Do not modify this section 

!
! Adjust index ranges
      il = ilays-1
      jl = (jlays-1)*k2d
      kl = (klays-1)*k3d

!
! For purposes of setting boundary conditions it is often useful to know which
! block boundaries are selected by the index ranges in the argument list.
!
! If the selected index range is not on an x face of the block      iface = 0
! If the selected index range is on the lower x face of the block   iface = -1
! If the selected index range is on the upper x face of the block   iface = +1
!
! Set jface and kface in similar fashion for y and z.

      iface = 0
      jface = 0
      kface = 0

      if(ilays.gt.0.and.ilays.lt.nxb) then
        if(id.le.nguard) then
           iface = -1
        else
           iface = 1
        endif
      endif
      if(ndim.ge.2) then
      if(jlays.gt.0.and.jlays.lt.nyb) then
        if(jd.le.nguard) then
           jface = -1
        else
           jface = 1
        endif
      endif
      endif
      if(ndim.eq.3) then
      if(klays.gt.0.and.klays.lt.nzb) then
        if(kd.le.nguard) then
           kface = -1
        else
           kface = 1
        endif
      endif
      endif


!---------------------------------------------------------------------------
! Section to be modified by user


! Which boundary condition has been specified?


      if(ibc.eq.-21) then

! reflecting

          lwb = idest

!
! Do cell-face-centered data
        if(nfacevar.gt.0.and.iopt.eq.1) then

          id1 = id + ip1
          jd1 = jd + jp1*k2d
          kd1 = kd + kp1*k3d


!          facevarx1(:,id1:id1+il,jd1:jd1+jl,kd1:kd1+kl,idest) =
!     .            ????


          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)-gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)-gc_off_y)*k2d
               endif
              do i = id1,id1+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+2
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1)-gc_off_x)
                 endif

                facevarx1(:nfacevar,i,j,k,lwb) =  & 
     &                            facevarx1(:nfacevar,is,js,ks,lwb)
              enddo
            enddo
          enddo


!              facevary1(:,id1:id1+il,jd1:jd1+jl,kd1:kd1+kl,idest) =
!     .            ????

          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)-gc_off_z)*k3d
            endif
            do j = jd1,jd1+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y+1)*k2d+1
               elseif(jface.eq.1) then
                   js = (nyb+(nguard+1)*k2d) & 
     &                 -(j-(nyb+nguard+1)-gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)-gc_off_x)
                 endif

                facevary1(:nfacevar,i,j,k,lwb) = & 
     &                            facevary1(:nfacevar,is,js,ks,lwb)
              enddo
            enddo
          enddo

!              facevarz1(:,id1:id1+il,jd1:jd1+jl,kd1:kd1+kl,idest) =
!     .            ????

          do k = kd1,kd1+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z+1)*k3d+1
            elseif(kface.eq.1) then
               ks = (nzb+(nguard+1)*k3d) & 
     &             -(k-(nzb+nguard+1)-gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)-gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)-gc_off_x)
                 endif

                facevarz1(:nfacevar,i,j,k,lwb) = & 
     &                            facevarz1(:nfacevar,is,js,ks,lwb)
              enddo
            enddo
          enddo


        endif                            ! end of nfacevar if test


!
! Now do cell centered data

        if(iopt.eq.1) then

!          unk1(:,id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

          lbw = idest
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)-gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)-gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard-(i-(nxb+nguard+1)-gc_off_x)
                 endif

                unk1(:nvar,i,j,k,lbw) = unk1(:nvar,is,js,ks,lbw)
              enddo
            enddo
          enddo

        elseif(iopt.ge.2) then

!          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

          lbw = idest
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard_work-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard_work*k3d- & 
     &               (k-(nzb+nguard_work+1)-gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard_work-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard_work*k2d- & 
     &                   (j-(nyb+nguard_work+1)-gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard_work-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard_work- & 
     &                     (i-(nxb+nguard_work+1)-gc_off_x)
                 endif

                work1(i,j,k,lbw) = work1(is,js,ks,lbw)

              enddo
            enddo
          enddo


        endif
!-------------------------

!       elseif(ibc.eq.??) then


       endif                            ! end of test of bc flag


! End of Section to be modified by user
!---------------------------------------------------------------------------

      return
      end
