
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f

      subroutine amr_1blk_bcset(mype,ibc,lb,pe, & 
     &    idest,iopt,ibnd,jbnd,kbnd,surrblks)

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
!   ibc = -21                       zero gradient
!   ibc = -22                       reflecting
!
!
! It is important to point out that any vector components stored
! in either unk or unk_n will have different reflection properties
! across different boundary planes.
! Note, in this example the following case is assumed:
!           unk(1,....)  is a scalar
!           unk(2,....)  is the x component of a vector
!           unk(3,....)  is the y component of a vector
!           unk(4,....)  is the z component of a vector
!           unk(5,....)  is a scalar
!
!           unk_n        all components are scalar.
!
! This is significant because it means that when ibc=-22, unk(2,...)
! changes it sign across x boundary planes but not across y or z
! boundary planes. Similarly unk(3,...) changes it sign across y 
! boundary planes but not across x or z boundary planes.
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
!      ibnd             a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the i index range, eg
!                           ibnd  = -1      left end
!                                 =  0      middle
!                                 = +1      right. For example, if iface=-1,
!                        the i index applied when filling unk will run
!                        from 1:nguard, if iface=0 from 1+nguard:nxb+nguard,
!                        and if iface=+1 from nxb+nguard+1:nxb+2*nguard.
!      jbnd             a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the j index range.
!      kbnd             a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the k index range.
!
!
! Written :     Peter MacNeice          August 1998
! Modified:     Peter MacNeice          January 2001
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      include 'mpif.h'

      integer, intent(in) :: mype,ibc,lb,pe
      integer, intent(in) :: idest,iopt,ibnd,jbnd,kbnd
      integer, intent(in) :: surrblks(:,:,:,:)

      real    :: ccoord(3),csize(3)
      save    ccoord,csize

      integer :: il,jl,kl,id1,jd1,kd1,i,j,k
      integer :: ks,js,is,lbw,ivar
      integer :: iface,jface,kface
      integer :: iface0,jface0,kface0
      real    :: fact

!---------------------------------------------------------------------------
! Do not modify this section 

      iface = ibnd
      jface = jbnd
      kface = kbnd

      if(iopt.eq.1) then
        nguard0 = nguard
      elseif(iopt.gt.1) then
        nguard0 = nguard_work
      endif

      ilays = nguard0
      jlays = nguard0*k2d
      klays = nguard0*k3d
      if(iface.eq.0) ilays = nxb
      if(jface.eq.0) jlays = nyb
      if(kface.eq.0) klays = nzb

      ip1 = 0
      jp1 = 0
      kp1 = 0
      if(iface.eq.1) ip1 = 1
      if(jface.eq.1) jp1 = k2d
      if(kface.eq.1) kp1 = k3d

      
      id = 1
      if(iface.eq.0)  id = 1+nguard0
      if(iface.eq.+1) id = nxb+1+nguard0
      jd = 1
      if(jface.eq.0)  jd = 1+nguard0*k2d
      if(jface.eq.+1) jd = nyb+(1+nguard0)*k2d
      kd = 1
      if(kface.eq.0)  kd = 1+nguard0*k3d
      if(kface.eq.+1) kd = nzb+(1+nguard0)*k3d

!
! Adjust index ranges
      il = ilays-1
      jl = (jlays-1)*k2d
      kl = (klays-1)*k3d


!
! Now reset iface,jface,kface so that blocks next to a boundary
! which treat their diagonal and corner guardcells correctly.
      iface0 = iface
      jface0 = jface
      kface0 = kface
      if(iface.ne.0.and.surrblks(1,iface+2,2,2).gt.-20) iface=0
      if(jface.ne.0.and.surrblks(1,2,jface+2,2).gt.-20) jface=0
      if(kface.ne.0.and.surrblks(1,2,2,kface+2).gt.-20) kface=0


!---------------------------------------------------------------------------
! Section to be modified by user


! Which boundary condition has been specified?


      if(ibc.eq.-22) then
!
!-------------------------
! reflecting - all variables

        lbw = idest

        if(iopt.eq.1) then

!
! Do cell-face-centered data
        if(nfacevar.gt.0) then

          id1 = id + ip1
          jd1 = jd + jp1*k2d
          kd1 = kd + kp1*k3d

! set the following index range in facevarx1
!
!         facevarx1(:,id1:id1+il,jd:jd+jl,kd:kd+kl,idest) =
!     .            ????

          il0 = 0
          if(iface0.eq.0) il0 = 1

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
              do i = id1 ,id1+il+il0
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


! set the following index range in facevary1
!
!              facevary1(:,id:id+il,jd1:jd1+jl,kd:kd+kl,idest) =
!     .            ????

          jl0 = 0
          if(jface0.eq.0) jl0 = k2d

          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd1,jd1+jl+jl0
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

! set the following index range in facevarz1
!
!              facevarz1(:,id:id+il,jd:jd+jl,kd1:kd1+kl,idest) =
!     .            ????

          kl0 = 0
          if(kface0.eq.0) kl0 = k3d

          fact = 1.
          do k = kd1,kd1+kl+kl0
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


        if(nvar.gt.0) then

! set the following index range in unk1
!
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
                  elseif(ndim.eq.3.and.ivar.eq.4 & 
     &                                  .and.abs(kface).ne.0)then
                    fact=-1.0
                  endif
                  unk1(ivar,i,j,k,lbw) = fact*unk1(ivar,is,js,ks,lbw)
                enddo
              enddo
            enddo
          enddo

        endif                     ! end of nvar if test

!
! Now do cell corner data

        if(nvarcorn.gt.0) then

! set the appropriate index range in unk_n1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvarcorn
                   fact= 1.0
                   unk_n1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_n1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          endif                    ! end of nvarcorn if test

!
! Now do cell edge centered data

        if (ndim > 1) then

        if(nvaredge.gt.0) then
! First unk_e_x1
!
          il_extra = 0
          iu_extra = 0
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = nguard+(nguard+1-i)
                     fact = -1.
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard))
                     fact = -1.
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvaredge
                   unk_e_x1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_e_x1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

! Now unk_e_y1
!
          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = 1+(nguard-(j-nguard))*k2d
                  fact = -1.
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard))*k2d
                  fact = -1.
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (1+nguard)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (1+(nxb+nguard))-(i-(nxb+nguard+1))
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvaredge
                   unk_e_y1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_e_y1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          if (ndim == 3) then
! finally unk_e_z1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard-k)*k3d
                fact = -1.
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard))*k3d
                fact = -1.
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvaredge
                   unk_e_z1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_e_z1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          end if
          end if

         endif                    ! end of nvaredge if test


        elseif(iopt.ge.2) then

! set the following index range in work1
!
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


!-------------------------
! zero gradient - all variables

        lbw = idest

        if(iopt.eq.1) then

!
! Do cell-face-centered data
        if(nfacevar.gt.0) then

          id1 = id + ip1
          jd1 = jd + jp1*k2d
          kd1 = kd + kp1*k3d

! set the following index range in facevarx1
!
!         facevarx1(:,id1:id1+il,jd:jd+jl,kd:kd+kl,idest) =
!     .            ????

          il0 = 0
          if(iface0.eq.0) il0 = 1

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
              do i = id1 ,id1+il+il0
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+2
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1)+gc_off_x)
                 endif
                facevarx1(:nfacevar,i,j,k,lbw) =  & 
     &                    facevarx1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo


! set the following index range in facevary1
!
!              facevary1(:,id:id+il,jd1:jd1+jl,kd:kd+kl,idest) =
!     .            ????

          jl0 = 0
          if(jface0.eq.0) jl0 = k2d

          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd1,jd1+jl+jl0
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
     &                     facevary1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo

! set the following index range in facevarz1
!
!              facevarz1(:,id:id+il,jd:jd+jl,kd1:kd1+kl,idest) =
!     .            ????

          kl0 = 0
          if(kface0.eq.0) kl0 = k3d

          fact = 1.
          do k = kd1,kd1+kl+kl0
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
     &                     facevarz1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo


        endif                            ! end of nfacevar if test


!
! Now do cell centered data


        if(nvar.gt.0) then

! set the following index range in unk1
!
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

                do ivar = 1,nvar
                  fact= 1.0
                  unk1(ivar,i,j,k,lbw) = fact*unk1(ivar,is,js,ks,lbw)
                enddo
              enddo
            enddo
          enddo

        endif                     ! end of nvar if test

!
! Now do cell corner data

        if(nvarcorn.gt.0) then

! set the appropriate index range in unk_n1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

                 do ivar = 1,nvarcorn
                   fact= 1.0
                   unk_n1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_n1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          endif                    ! end of nvarcorn if test

!
! Now do cell edge centered data

        if (ndim > 1) then

        if(nvaredge.gt.0) then
! First unk_e_x1
!
          il_extra = 0
          iu_extra = 0
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = nguard+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard))
                 endif

                 do ivar = 1,nvaredge
                   unk_e_x1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_e_x1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

! Now unk_e_y1
!
          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = 1+(nguard-(j-nguard))*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (1+nguard)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (1+(nxb+nguard))-(i-(nxb+nguard+1))
                 endif

                 do ivar = 1,nvaredge
                   unk_e_y1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_e_y1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          if (ndim == 3) then
! finally unk_e_z1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

                 do ivar = 1,nvaredge
                   unk_e_z1(ivar,i,j,k,lbw) =  & 
     &                     fact*unk_e_z1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          end if
          end if

         endif                    ! end of nvaredge if test


        elseif(iopt.ge.2) then

! set the following index range in work1
!
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


!
!-------------------------


       endif                            ! end of test of bc flag


! End of Section to be modified by user
!---------------------------------------------------------------------------

      return
      end subroutine amr_1blk_bcset
