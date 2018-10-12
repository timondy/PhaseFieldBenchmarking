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


      subroutine amr_1blk_fc_prol_gen_fun( & 
     &       recv,ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff, & 
     &       mype,lb,pe_p,lb_p,iface)



!
!------------------------------------------------------------------------
!
! This routine is a wrapper routine which calls the functions
! which prolong data for FACEVAR. The local logical array lmask can
! be used to control which routine is actually operating on
! each variable stored within FACEVAR.
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
     &                  amr_1blk_fc_prol_inject, & 
     &                  amr_1blk_fc_prol_linear, & 
     &                  amr_1blk_fc_prol_genorder, & 
     &                  amr_1blk_fc_prol_user

      implicit none

!------------------------------------

      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)

!------------------------------------
! local variables

      integer :: interp_mask_face(nbndvar), ivar



      include 'mpif.h'
      double precision :: time1

!------------------------------------

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

      if(iface.eq.1) then
        interp_mask_face = interp_mask_facex
      elseif(iface.eq.2) then
        interp_mask_face = interp_mask_facey
      elseif(iface.eq.3) then
        interp_mask_face = interp_mask_facez
      endif

      do ivar = 1, nbndvar

      if (any(int_gcell_on_fc(1:ndim,ivar))) then
         
      if (interp_mask_face(ivar) < 20) then

      if (interp_mask_face(ivar) == 0) then
! Injection - zeroth order prolongation 
         call amr_1blk_fc_prol_inject( & 
     &        recv,ia,ib,ja,jb,ka,kb, & 
     &        idest,ioff,joff,koff,mype,iface,ivar)
         
      elseif (interp_mask_face(ivar) == 1) then
! Tri-linear interpolation 
         call amr_1blk_fc_prol_linear( & 
     &        recv,ia,ib,ja,jb,ka,kb, & 
     &        idest,ioff,joff,koff,mype,iface,ivar)
         
      elseif (interp_mask_face(ivar) > 1) then
! General order
         call amr_1blk_fc_prol_genorder( & 
     &        recv,ia,ib,ja,jb,ka,kb, & 
     &        idest,ioff,joff,koff,mype,iface,ivar, & 
     &        interp_mask_face(ivar))

      end if

      elseif (interp_mask_face(ivar) == 20) then

! User defined
         call amr_1blk_fc_prol_user( & 
     &        recv,ia,ib,ja,jb,ka,kb, & 
     &        idest,ioff,joff,koff, & 
     &        mype,lb,pe_p,lb_p,iface,ivar)

      end if

      end if

      end do


      if (timing_mpi) then
              timer_amr_1blk_fc_prol_gen =   & 
     &                          timer_amr_1blk_fc_prol_gen & 
     &                          + mpi_wtime() - time1
      endif

      return
      end subroutine amr_1blk_fc_prol_gen_fun
