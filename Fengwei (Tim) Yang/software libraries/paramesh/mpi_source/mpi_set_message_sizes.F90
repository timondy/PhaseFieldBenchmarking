!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine mpi_set_message_sizes(iopt, & 
     &                                 nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine sets up the sizes of messages required when only a part of
! a block is to be fetched.
!
!
! Written :     Peter MacNeice          April 2001
!------------------------------------------------------------------------
!
! Arguments:
!      iopt           integer switch controlling whether work dat or other data
!                     is sized 
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      implicit none

      integer, intent(in)  ::  iopt
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

! local variables

      integer :: nguard0
      integer :: nlayers0x,nlayers0y,nlayers0z
      integer :: i, lnxb, lnyb, lnzb, i2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(iopt.eq.1) then
        nguard0 = nguard
      else
        nguard0 = nguard_work
      endif

      if(iopt.eq.1) then
         if (.not.present(nlayersx)) then
            nlayers0x = nguard
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard
         else
            nlayers0z = nlayersz
         end if
      else
         if (.not.present(nlayersx)) then
            nlayers0x = nguard_work
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard_work
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard_work
         else
            nlayers0z = nlayersz
         end if
      endif

! loop over possible message types
      do i=1,27

! set x index extent
        if(mod(i,3).eq.1) then
          lnxb = nlayers0x
        elseif(mod(i,3).eq.2) then
          lnxb = nxb
        elseif(mod(i,3).eq.0) then
          lnxb = nlayers0x
        endif
! set y index extent
        if(mod((i-1)/3,3).eq.0) then
          lnyb = nlayers0y*k2d
        elseif(mod((i-1)/3,3).eq.1) then
          lnyb = nyb
        elseif(mod((i-1)/3,3).eq.2) then
          lnyb = nlayers0y*k2d
        endif
! set z index extent
        if(i.le.9) then
          lnzb = nlayers0z*k3d
        elseif(i.ge.10.and.i.le.18) then
          lnzb = nzb
        elseif(i.ge.19) then
          lnzb = nlayers0z*k3d
        endif

        if(iopt.gt.1) then
          message_size_wk(i) = lnxb*lnyb*lnzb 
        else
          message_size_cc(i) = lnxb*lnyb*lnzb
          message_size_nc(i) = (lnxb+1)*(lnyb+k2d)*(lnzb+k3d)
          if(ndim.eq.1) then
            message_size_fc(i) = (lnxb+1)*lnyb*lnzb
            message_size_ec(i) = 0
          elseif(ndim.eq.2.and.l2p5d.eq.0) then
            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
     &                           lnxb*(lnyb+k2d)*lnzb
            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
     &                           (lnxb+1)*lnyb*(lnzb+k3d)  
          elseif(ndim.eq.3.or.l2p5d.eq.1) then
            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
     &           lnxb*(lnyb+k2d)*lnzb + lnxb*lnyb*(lnzb+k3d) 
            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
     &                           (lnxb+1)*lnyb*(lnzb+k3d) + & 
     &                           (lnxb+1)*(lnyb+k2d)*lnzb
          endif
        endif
#ifdef DEBUG
        write(*,*) 'message sizes ',i,' cc/nc/fc/ec ', & 
     &  message_size_cc(i),message_size_nc(i), & 
     &  message_size_fc(i),message_size_ec(i)
#endif /* DEBUG */

      enddo


! loop over possible message types

      if (iopt == 1) then
         nlayers0x = min(2*nlayers0x+1,2*nguard)
         nlayers0y = min(2*nlayers0y+1,2*nguard)
         nlayers0z = min(2*nlayers0z+1,2*nguard)
      else
         nlayers0x = min(2*nlayers0x+1,2*nguard_work)
         nlayers0y = min(2*nlayers0y+1,2*nguard_work)
         nlayers0z = min(2*nlayers0z+1,2*nguard_work)
      end if

      do i=27+1,2*27

        i2 = i - 27
! set x index extent
        if(mod(i2,3).eq.1) then
          lnxb = nlayers0x
        elseif(mod(i2,3).eq.2) then
          lnxb = nxb
        elseif(mod(i2,3).eq.0) then
          lnxb = nlayers0x
        endif
! set y index extent
        if(mod((i2-1)/3,3).eq.0) then
          lnyb = nlayers0y*k2d
        elseif(mod((i2-1)/3,3).eq.1) then
          lnyb = nyb
        elseif(mod((i2-1)/3,3).eq.2) then
          lnyb = nlayers0y*k2d
        endif
! set z index extent
        if(i2.le.9) then
          lnzb = nlayers0z*k3d
        elseif(i2.ge.10.and.i2.le.18) then
          lnzb = nzb
        elseif(i2.ge.19) then
          lnzb = nlayers0z*k3d
        endif

        if(iopt.gt.1) then
          message_size_wk(i) = lnxb*lnyb*lnzb 
        else
          message_size_cc(i) = lnxb*lnyb*lnzb
          message_size_nc(i) = (lnxb+1)*(lnyb+k2d)*(lnzb+k3d)
          if(ndim.eq.1) then
            message_size_fc(i) = (lnxb+1)*lnyb*lnzb
            message_size_ec(i) = 0
          elseif(ndim.eq.2.and.l2p5d.eq.0) then
            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
     &                           lnxb*(lnyb+k2d)*lnzb
            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
     &                           (lnxb+1)*lnyb*(lnzb+k3d)  
          elseif(ndim.eq.3.or.l2p5d.eq.1) then
            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
     &           lnxb*(lnyb+k2d)*lnzb + lnxb*lnyb*(lnzb+k3d) 
            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
     &                           (lnxb+1)*lnyb*(lnzb+k3d) + & 
     &                           (lnxb+1)*(lnyb+k2d)*lnzb
          endif
        endif

#ifdef DEBUG
        write(*,*) 'message sizes ',i,' cc/nc/fc/ec ', & 
     &  message_size_cc(i),message_size_nc(i), & 
     &  message_size_fc(i),message_size_ec(i)
#endif /* DEBUG */

      end do

! Commented out by CEG
! 2    continue
      continue

      return
      end subroutine mpi_set_message_sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CEG version

      subroutine pf_mpi_set_message_sizes(iopt)

!------------------------------------------------------------------------
!
! This subroutine sets up the sizes of messages required when only a part of
! a block is to be fetched.
!
!
! Written :     Peter MacNeice          April 2001
!------------------------------------------------------------------------
!
! Arguments:
!      iopt           integer switch controlling whether work dat or other data
!                     is sized 
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      implicit none

      integer, intent(in)  ::  iopt

! local variables

      integer :: nguard0
      integer :: nlayers0x,nlayers0y,nlayers0z
      integer :: i, lnxb, lnyb, lnzb, i2

! CEG first time switch
      logical :: ioptfirst1=.true.
      logical :: ioptfirst2=.true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((ioptfirst1.and.iopt.eq.1).or.(ioptfirst2.and.iopt.gt.1)) then

      if(iopt.eq.1) then
         nguard0 = nguard
         nlayers0x = nguard
         nlayers0y = nguard
         nlayers0z = nguard
         ioptfirst1 = .false.
      else
         nguard0 = nguard_work
         nlayers0x = nguard_work
         nlayers0y = nguard_work
         nlayers0z = nguard_work
         ioptfirst2 = .false.
      endif

! loop over possible message types
      do i=1,27

! set x index extent
        if(mod(i,3).eq.1) then
          lnxb = nlayers0x
        elseif(mod(i,3).eq.2) then
          lnxb = nxb
        elseif(mod(i,3).eq.0) then
          lnxb = nlayers0x
        endif
! set y index extent
        if(mod((i-1)/3,3).eq.0) then
          lnyb = nlayers0y*k2d
        elseif(mod((i-1)/3,3).eq.1) then
          lnyb = nyb
        elseif(mod((i-1)/3,3).eq.2) then
          lnyb = nlayers0y*k2d
        endif
! set z index extent
        if(i.le.9) then
          lnzb = nlayers0z*k3d
        elseif(i.ge.10.and.i.le.18) then
          lnzb = nzb
        elseif(i.ge.19) then
          lnzb = nlayers0z*k3d
        endif

        if(iopt.gt.1) then
          message_size_wk(i) = lnxb*lnyb*lnzb 
        else
          message_size_cc(i) = lnxb*lnyb*lnzb
!          message_size_nc(i) = (lnxb+1)*(lnyb+k2d)*(lnzb+k3d)
          message_size_nc(i) = 0
          if(ndim.eq.1) then
!            message_size_fc(i) = (lnxb+1)*lnyb*lnzb
            message_size_fc(i) = 0
            message_size_ec(i) = 0
          elseif(ndim.eq.2.and.l2p5d.eq.0) then
!            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
!     &                           lnxb*(lnyb+k2d)*lnzb
!            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
!     &                           (lnxb+1)*lnyb*(lnzb+k3d)  
            message_size_fc(i) = 0
            message_size_ec(i) = 0
          elseif(ndim.eq.3.or.l2p5d.eq.1) then
!            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
!     &           lnxb*(lnyb+k2d)*lnzb + lnxb*lnyb*(lnzb+k3d) 
!            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
!     &                           (lnxb+1)*lnyb*(lnzb+k3d) + & 
!     &                           (lnxb+1)*(lnyb+k2d)*lnzb
            message_size_fc(i) = 0
            message_size_ec(i) = 0
          endif
        endif
#ifdef DEBUG
        write(*,*) 'message sizes ',i,' cc/nc/fc/ec ', & 
     &  message_size_cc(i),message_size_nc(i), & 
     &  message_size_fc(i),message_size_ec(i)
#endif /* DEBUG */

      enddo

! loop over possible message types

      if (iopt == 1) then
         nlayers0x = min(2*nlayers0x+1,2*nguard)
         nlayers0y = min(2*nlayers0y+1,2*nguard)
         nlayers0z = min(2*nlayers0z+1,2*nguard)
      else
         nlayers0x = min(2*nlayers0x+1,2*nguard_work)
         nlayers0y = min(2*nlayers0y+1,2*nguard_work)
         nlayers0z = min(2*nlayers0z+1,2*nguard_work)
      end if

      do i=27+1,2*27

        i2 = i - 27
! set x index extent
        if(mod(i2,3).eq.1) then
          lnxb = nlayers0x
        elseif(mod(i2,3).eq.2) then
          lnxb = nxb
        elseif(mod(i2,3).eq.0) then
          lnxb = nlayers0x
        endif
! set y index extent
        if(mod((i2-1)/3,3).eq.0) then
          lnyb = nlayers0y*k2d
        elseif(mod((i2-1)/3,3).eq.1) then
          lnyb = nyb
        elseif(mod((i2-1)/3,3).eq.2) then
          lnyb = nlayers0y*k2d
        endif
! set z index extent
        if(i2.le.9) then
          lnzb = nlayers0z*k3d
        elseif(i2.ge.10.and.i2.le.18) then
          lnzb = nzb
        elseif(i2.ge.19) then
          lnzb = nlayers0z*k3d
        endif

        if(iopt.gt.1) then
          message_size_wk(i) = lnxb*lnyb*lnzb 
        else
          message_size_cc(i) = lnxb*lnyb*lnzb
!          message_size_nc(i) = (lnxb+1)*(lnyb+k2d)*(lnzb+k3d)
          message_size_nc(i) = 0
          if(ndim.eq.1) then
!            message_size_fc(i) = (lnxb+1)*lnyb*lnzb
!            message_size_ec(i) = 0
            message_size_fc(i) = 0
            message_size_ec(i) = 0
          elseif(ndim.eq.2.and.l2p5d.eq.0) then
!            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
!     &                           lnxb*(lnyb+k2d)*lnzb
!            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
!     &                           (lnxb+1)*lnyb*(lnzb+k3d)  
            message_size_fc(i) = 0
            message_size_ec(i) = 0
          elseif(ndim.eq.3.or.l2p5d.eq.1) then
!            message_size_fc(i) = (lnxb+1)*lnyb*lnzb +  & 
!     &           lnxb*(lnyb+k2d)*lnzb + lnxb*lnyb*(lnzb+k3d) 
!            message_size_ec(i) = lnxb*(lnyb+k2d)*(lnzb+k3d) + & 
!     &                           (lnxb+1)*lnyb*(lnzb+k3d) + & 
!     &                           (lnxb+1)*(lnyb+k2d)*lnzb
            message_size_fc(i) = 0
            message_size_ec(i) = 0
          endif
        endif

#ifdef DEBUG
        write(*,*) 'message sizes ',i,' cc/nc/fc/ec ', & 
     &  message_size_cc(i),message_size_nc(i), & 
     &  message_size_fc(i),message_size_ec(i)
#endif /* DEBUG */

      end do

! Commented out by CEG
! 2    continue
      continue

    endif

      return
      end subroutine pf_mpi_set_message_sizes
