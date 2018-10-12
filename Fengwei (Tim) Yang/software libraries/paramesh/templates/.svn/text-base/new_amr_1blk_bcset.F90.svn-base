!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_bcset(mype2,ibc,lb,pe, & 
     &    idest,iopt,id,jd,kd,ilays,jlays,klays,ip1,jp1,kp1)




!------------------------------------------------------------------------
!
! This routine sets guardcell values at external boundaries in the case
! where a single block is having its guardcells filled.
!
! It can be assumed in writing this routine, that all guardcells for this
! block which are not across an external boundary have already been
! properly filled.
!
! This routine provides code to set boundary conditions for a selection
! of common boundary conditions.
!
!   ibc = -20                       reflecting
!   ibc = -21                       zero gradient
!
!------------------------------------------------------------------------
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

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      implicit none

      include 'mpif.h'

      integer mype2,ibc,lb,pe
      integer idest,iopt,id,jd,kd,ilays,jlays,klays,ip1,jp1,kp1


      real ccoord(3),csize(3)
      save    ccoord,csize

      integer il,jl,kl,iface,jface,kface,id1,jd1,kd1,i,j,k
      integer if_l,if_r,jf_l,jf_r,kf_l,kf_r
      integer ks,js,is,lbw,ivar

      real :: fact, kine

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
        if(id.le.nguard.and.neigh(1,1,lb).le.-20) then
           iface = -1
        elseif (id.gt.nguard.and.neigh(1,2,lb).le.-20) then
           iface = 1
        endif
      endif
      if(ndim.ge.2) then
      if(jlays.gt.0.and.jlays.lt.nyb) then
        if(jd.le.nguard.and.neigh(1,3,lb).le.-20) then
           jface = -1
        elseif (jd.gt.nguard.and.neigh(1,4,lb).le.-20) then
           jface = 1
        endif
      endif
      endif
      if(ndim.eq.3) then
      if(klays.gt.0.and.klays.lt.nzb) then
        if(kd.le.nguard.and.neigh(1,5,lb).le.-20) then
           kface = -1
        elseif (kd.gt.nguard.and.neigh(1,6,lb).le.-20) then
           kface = 1
        endif
      endif
      endif

! When using an odd sized grid, you need to set an additional face
! bounding the last interior cell, since this face is outside the
! physical boundary.
          if_l = 0
          if_r = 0
          if(iface.eq.-1) if_l = gc_off_x
          if(iface.eq. 1) if_r = gc_off_x
          jf_l = 0
          jf_r = 0
          if(jface.eq.-1) jf_l = gc_off_y*k2d
          if(jface.eq. 1) jf_r = gc_off_y*k2d
          kf_l = 0
          kf_r = 0
          if(kface.eq.-1) kf_l = gc_off_z*k3d
          if(kface.eq. 1) kf_r = gc_off_z*k3d


!---------------------------------------------------------------------------
! Section to be modified by user


! Which boundary condition has been specified?


      if(ibc.eq.-20) then
!
!-------------------------
! reflecting - all variables

          lbw = idest

!
! Do cell-face-centered data
        if(nfacevar.gt.0.and.iopt.eq.1) then

          id1 = id + ip1
          jd1 = jd + jp1*k2d
          kd1 = kd + kp1*k3d

!          facevarx1(:,id1:id1+il,jd:jd+jl,kd:kd+kl,idest) =
!     .            ????


          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id1-if_r ,id1+il+if_l
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+2
                     fact = -1.
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1)+gc_off_x)
                     fact = -1.
                 endif

                facevarx1(:nfacevar,i,j,k,lbw) =  & 
     &                    facevarx1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo


!              facevary1(:,id:id+il,jd1:jd1+jl,kd:kd+kl,idest) =
!     .            ????

          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd1-jf_r,jd1+jl+jf_l
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   fact = -1.
                   js = (2*nguard-j+gc_off_y+1)*k2d+1
               elseif(jface.eq.1) then
                   fact = -1.
                   js = (nyb+(nguard+1)*k2d) & 
     &                 -(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevary1(:nfacevar,i,j,k,lbw) = & 
     &                     facevary1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo

!              facevarz1(:,id:id+il,jd:jd+jl,kd1:kd1+kl,idest) =
!     .            ????

          fact = 1.
          do k = kd1-kf_r,kd1+kl+kf_l
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               fact = -1.
               ks = (2*nguard-k+gc_off_z+1)*k3d+1
            elseif(kface.eq.1) then
               fact = -1.
               ks = (nzb+(nguard+1)*k3d) & 
     &             -(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevarz1(:nfacevar,i,j,k,lbw) = & 
     &                     facevarz1(:nfacevar,is,js,ks,lbw)*fact
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
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard-(i-(nxb+nguard+1)+gc_off_x)
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                do ivar = 1,nvar
                     fact= 1.0
                  if(ivar.eq.2.and.abs(iface).ne.0)then
                     fact=-1.0
                  elseif(ivar.eq.3.and.abs(jface).ne.0)then
                     fact=-1.0
#if NDIM == 3
                  elseif(ivar.eq.4.and.abs(kface).ne.0)then
                     fact=-1.0
#endif
                  endif
                  unk1(ivar,i,j,k,lbw) = fact*unk1(ivar,is,js,ks,lbw)
                enddo
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
     &               (k-(nzb+nguard_work+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard_work-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard_work*k2d- & 
     &                   (j-(nyb+nguard_work+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard_work-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard_work- & 
     &                     (i-(nxb+nguard_work+1)+gc_off_x)
                 endif

                work1(i,j,k,lbw) = work1(is,js,ks,lbw)

              enddo
            enddo
          enddo


        endif


!-------------------------


      elseif(ibc.eq.-21) then
!
!-------------------------
! zero gradient - all variables

          lbw = idest

!
! Do cell-face-centered data
        if(nfacevar.gt.0.and.iopt.eq.1) then

          id1 = id + ip1
          jd1 = jd + jp1*k2d
          kd1 = kd + kp1*k3d


!          facevarx1(:,id1:id1+il,jd:jd+jl,kd:kd+kl,idest) =
!     .            ????


          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id1-if_r,id1+il+if_l
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+2
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevarx1(:nfacevar,i,j,k,lbw) =  & 
     &                            facevarx1(:nfacevar,is,js,ks,lbw)
              enddo
            enddo
          enddo


!              facevary1(:,id:id+il,jd1:jd1+jl,kd:kd+kl,idest) =
!     .            ????

          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd1-jf_r,jd1+jl+jf_l
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y+1)*k2d+1
               elseif(jface.eq.1) then
                   js = (nyb+(nguard+1)*k2d) & 
     &                 -(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevary1(:nfacevar,i,j,k,lbw) = & 
     &                            facevary1(:nfacevar,is,js,ks,lbw)
              enddo
            enddo
          enddo

!              facevarz1(:,id:id+il,jd:jd+jl,kd1:kd1+kl,idest) =
!     .            ????

          do k = kd1-kf_r,kd1+kl+kf_l
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z+1)*k3d+1
            elseif(kface.eq.1) then
               ks = (nzb+(nguard+1)*k3d) & 
     &             -(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevarz1(:nfacevar,i,j,k,lbw) = & 
     &                            facevarz1(:nfacevar,is,js,ks,lbw)
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
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard-(i-(nxb+nguard+1)+gc_off_x)
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
     &               (k-(nzb+nguard_work+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard_work-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard_work*k2d- & 
     &                   (j-(nyb+nguard_work+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard_work-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard_work- & 
     &                     (i-(nxb+nguard_work+1)+gc_off_x)
                 endif

                work1(i,j,k,lbw) = work1(is,js,ks,lbw)

              enddo
            enddo
          enddo


        endif
!-------------------------

      elseif (ibc.eq.-23) then


          lbw = idest
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                 unk1(igame,i,j,k,lbw) = gamma
                 unk1(igamc,i,j,k,lbw) = gamma
                 
                 unk1(idens,i,j,k,lbw) =  rho_ambient
                 unk1(ivelx,i,j,k,lbw) = wind_vel
                 unk1(ively,i,j,k,lbw) = 0.
                 unk1(ivelz,i,j,k,lbw) = 0.
                 unk1(ipres,i,j,k,lbw) = p_ambient
                 kine = 0.5 * (unk1(ivelx,i,j,k,lbw)**2)
                 unk1(iener,i,j,k,lbw) =  & 
     &                unk1(ipres,i,j,k,lbw) / & 
     &                (unk1(igame,i,j,k,lbw)-1.)
                 unk1(iener,i,j,k,lbw) =  & 
     &                unk1(iener,i,j,k,lbw) / & 
     &                unk1(idens,i,j,k,lbw)
                 unk1(iener,i,j,k,lbw) =  & 
     &                unk1(iener,i,j,k,lbw)  & 
     &                + kine
                 unk1(iener,i,j,k,lbw) =  & 
     &                max(unk1(iener,i,j,k,lbw),smallp)
              enddo
           enddo
        enddo

      endif                     ! end of test of bc flag


! End of Section to be modified by user
!---------------------------------------------------------------------------

      return
      end






