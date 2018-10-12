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


      subroutine amr_1blk_cc_prol_gen_unk_fun & 
     &  (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &   mype,lb,pe_p,lb_p)


!
!------------------------------------------------------------------------
!
! This routine is a wrapper routine which calls the functions
! which prolong data for UNK. The local logical array lmask can
! be used to control which routine is actually operating on
! each variable stored within UNK.
!
!
! Written :     Peter MacNeice          January 2002
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use prolong_arrays

      use paramesh_interfaces, only :  & 
     &                  amr_1blk_cc_prol_inject, & 
     &                  amr_1blk_cc_prol_linear, & 
     &                  amr_1blk_cc_prol_genorder, & 
     &                  amr_1blk_cc_prol_user

      implicit none

!------------------------------------

      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)


      include 'mpif.h'
      double precision :: time1


!------------------------------------

! local variables

      integer :: ivar

!------------------------------------

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

      do ivar = 1, nvar
      if (int_gcell_on_cc(ivar)) then

      if (interp_mask_unk(ivar) < 20) then

      if (interp_mask_unk(ivar) == 0) then
! Simple Injection 
      call amr_1blk_cc_prol_inject & 
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &     mype,ivar)
      
      elseif (interp_mask_unk(ivar) == 1) then
! Default multi-linear interpolation  
      call amr_1blk_cc_prol_linear & 
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &     mype,ivar)

      elseif (interp_mask_unk(ivar) > 1) then
! General order
      call amr_1blk_cc_prol_genorder & 
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &      mype,ivar,interp_mask_unk(ivar))

      end if

      elseif (interp_mask_unk(ivar) == 20) then

! User defined interpolation to be used for prolocation

         call amr_1blk_cc_prol_user & 
     &        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &        mype,lb,pe_p,lb_p,ivar)

      end if

      end if
      end do

      if (timing_mpi) then
              timer_amr_1blk_cc_prol_gen_unk =  & 
     &             timer_amr_1blk_cc_prol_gen_unk & 
     &                          + mpi_wtime() - time1
      endif

      return
      end subroutine amr_1blk_cc_prol_gen_unk_fun
