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


      subroutine amr_1blk_cc_prol_gen_work_fun(recv, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)



!
!------------------------------------------------------------------------
!
! This routine is a wrapper routine which calls the functions
! which prolong data for WORK. The argument idest selects the layer
! within WORK on which prolongation is required.
! The argument interp selects the actual interpolation function to be
! called.
!
!
! Written :     Peter MacNeice          January 2002
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use prolong_arrays

      use paramesh_interfaces, only :  & 
     &                  amr_1blk_cc_prol_work_inject, & 
     &                  amr_1blk_cc_prol_work_linear, & 
     &                  amr_1blk_cc_prol_work_user, & 
     &                  amr_1blk_cc_prol_work_genorder


      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      integer, intent(in) :: interp
      real,    intent(inout) :: recv(:,:,:)

!------------------------------------

      if (interp < 20) then

! Simple Injection
      if(interp.eq.0) & 
     &  call amr_1blk_cc_prol_work_inject(recv, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype)

! Linear interpolation
      if(interp.eq.1) & 
     &   call amr_1blk_cc_prol_work_linear(recv, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype)

! General Order
      
      if (interp > 1)  & 
     &   call amr_1blk_cc_prol_work_genorder(recv, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,interp)

      elseif (interp == 20) then

      call amr_1blk_cc_prol_work_user(recv, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p)

      end if

      return
      end subroutine amr_1blk_cc_prol_gen_work_fun
