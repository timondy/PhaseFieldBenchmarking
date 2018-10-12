!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_gen_unk_fun
!! NAME
!!
!!   amr_1blk_cc_prol_gen_unk_fun
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_gen_unk_fun (recv, ia, ib, ja, jb, ka, kb,
!!                                      ioff, joff, koff, mype,
!!                                      lb, pe_p, lb_p)
!!   Call amr_1blk_cc_prol_gen_unk_fun (real, 
!!                                      integer, integer, integer, integer,
!!                                      integer, integer, integer, integer,
!!                                      integer, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(inout) :: recv(:,:,:,:)
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
!! INCLUDES
!!   
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   prolong_arrays
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_1blk_cc_prol_inject 
!!   amr_1blk_cc_prol_linear
!!   amr_1blk_cc_prol_genorder
!!   amr_1blk_cc_prol_user
!!
!! RETURNS
!!
!!   Upon return prolonged data is placed in unk1.
!!
!! DESCRIPTION
!!
!!   This routine is a wrapper routine which calls the functions
!!   which prolong data for UNK. The local logical array lmask can
!!   be used to control which routine is actually operating on
!!   each variable stored within UNK.
!!
!! AUTHORS
!!
!! Written by Peter MacNeice January 2002.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      Subroutine amr_1blk_cc_prol_gen_unk_fun          & 
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  & 
         mype,lb,pe_p,lb_p)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use prolong_arrays

      Use paramesh_interfaces, only :              & 
                       amr_1blk_cc_prol_inject,    & 
                       amr_1blk_cc_prol_linear,    & 
                       amr_1blk_cc_prol_genorder,  & 
                       amr_1blk_cc_prol_user

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output variables
      Real,    Intent(inout) :: recv(:,:,:,:)
      Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      Integer, Intent(in)    :: ioff,joff,koff,mype
      Integer, Intent(in)    :: lb,lb_p,pe_p

!-----Local variables
      Double Precision :: time1
      Integer :: ivar

!-----Begin Executable code

      If (timing_mpi) Then
         time1 = mpi_wtime()
      End If  ! End If (timing_mpi)

      Do ivar = 1, nvar
      If (int_gcell_on_cc(ivar)) Then

      If (interp_mask_unk(ivar) < 20) Then

      If (interp_mask_unk(ivar) == 0) then
!-----Simple Injection 
      Call amr_1blk_cc_prol_inject                       & 
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
            mype,ivar)
      
      Elseif (interp_mask_unk(ivar) == 1) Then
!-----Default multi-linear interpolation  
      Call amr_1blk_cc_prol_linear                       & 
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
           mype,ivar)

      Elseif (interp_mask_unk(ivar) > 1) Then
!-----High order Largrange ploynomial interpolation
      Call amr_1blk_cc_prol_genorder                     & 
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
            mype,ivar,interp_mask_unk(ivar))

      End If  ! End If (interp_mask_unk(ivar) == 0)

      Elseif (interp_mask_unk(ivar) == 20) Then

!--------User defined interpolation to be used for prolocation

         Call amr_1blk_cc_prol_user()

      End If  ! End If (interp_mask_unk(ivar) < 20

      End If  ! Enf If (int_gcell_on_cc(ivar))
      End Do  ! End Do ivar = 1, nvar

      If (timing_mpi) Then
              timer_amr_1blk_cc_prol_gen_unk =        & 
                   timer_amr_1blk_cc_prol_gen_unk     & 
                                + mpi_wtime() - time1
      End If  ! End If (timing_mpi)

      Return
      End Subroutine amr_1blk_cc_prol_gen_unk_fun
