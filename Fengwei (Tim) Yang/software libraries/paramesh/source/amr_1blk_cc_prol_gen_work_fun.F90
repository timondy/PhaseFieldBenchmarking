!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_gen_work_fun
!! NAME
!!
!!   amr_1blk_cc_prol_gen_work_fun
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_gen_work_fun (recvt, ia, ib, ja, jb, ka, kb,
!!                                      ioff, joff, koff, mype,
!!                                      lb, pe_p, lb_p, interp)
!!   Call amr_1blk_cc_prol_gen_work_fun (real, 
!!                                      integer, integer, integer, integer,
!!                                      integer, integer, integer, integer,
!!                                      integer, integer. integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(inout) :: recvt(:,:,:,:)
!!     Data to prolong.
!!
!!   Integer, Intent(in) :: ia,ib,ja,jb,ka,kb,idest
!!     Indeces in unk1 arrar to place prolonged data
!!
!!   Integer, Intent(in) :: ioff,joff,koff,mype
!!     Offsets and local pe id.
!!
!!   Integer, Intent(in) :: lb,lb_p,pe_p
!!     local block to prolong to block and pe of parent block
!!
!!   Integer, Intent(in) :: interp
!!     order of polynomial to use
!!
!! INCLUDES
!!   
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   workspace
!!   prolong_arrays
!!
!! CALLS
!!
!!   amr_1blk_cc_prol_work_inject
!!   amr_1blk_cc_prol_work_linear
!!   amr_1blk_cc_prol_work_user
!!   amr_1blk_cc_prol_work_genorder
!!
!! RETURNS
!!
!!   Upon return prolonged data is placed in work1.
!!
!! DESCRIPTION
!!
!!   This routine is a wrapper routine which calls the functions
!!   which prolong data for WORK. The local logical array lmask can
!!   be used to control which routine is actually operating on
!!   each variable stored within WORK.
!!
!! AUTHORS
!!
!! Written by Peter MacNeice January 2002.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      Subroutine amr_1blk_cc_prol_gen_work_fun(recvt,       & 
             ia,ib,ja,jb,ka,kb,                             & 
             idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)


!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use prolong_arrays

      Use paramesh_interfaces, Only :                  & 
                        amr_1blk_cc_prol_work_inject,  & 
                        amr_1blk_cc_prol_work_linear,  & 
                        amr_1blk_cc_prol_work_user,    & 
                        amr_1blk_cc_prol_work_genorder


      Implicit None

!-----Input/Output varaibles
      Real,    Intent(inout) :: recvt(:,:,:)
      Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      Integer, Intent(in)    :: ioff,joff,koff,mype
      Integer, Intent(in)    :: lb,lb_p,pe_p
      Integer, Intent(in)    :: interp

!-----Begin Executable Code

      If (interp < 20) Then

!-----Simple Injection
      If (interp.eq.0)                              & 
        Call amr_1blk_cc_prol_work_inject(recvt,    & 
             ia,ib,ja,jb,ka,kb,                     & 
             idest,ioff,joff,koff,mype)

!-----Linear interpolation
      If (interp.eq.1)                              & 
         Call amr_1blk_cc_prol_work_linear(recvt,   & 
             ia,ib,ja,jb,ka,kb,                     & 
             idest,ioff,joff,koff,mype)

!-----High order Lagrange ploynomail interpolation
      If (interp > 1)                               & 
         Call amr_1blk_cc_prol_work_genorder(recvt, & 
             ia,ib,ja,jb,ka,kb,                     & 
             idest,ioff,joff,koff,mype,interp)

      Elseif (interp == 20) Then

      Call amr_1blk_cc_prol_work_user()

      End If  ! End If (interp < 20)

      Return
      End Subroutine amr_1blk_cc_prol_gen_work_fun
