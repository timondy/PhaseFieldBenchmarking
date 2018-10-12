!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_fc_prol_gen_fun
!! NAME
!!
!!   amr_1blk_fc_prol_gen_fun
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_fc_prol_gen_fun (recv, ia, ib, ja, jb, ka, kb,
!!                                  ioff, joff, koff, mype,
!!                                  lb, pe_p, lb_p, iface)
!!   Call amr_1blk_fc_prol_gen_fun (real, 
!!                                  integer, integer, integer, integer,
!!                                  integer, integer, integer, integer,
!!                                  integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(inout) :: recv(:,:,:,:)
!!     Data to prolong.
!!
!!   Integer, Intent(in) :: ia,ib,ja,jb,ka,kb,idest
!!     Indeces in facevarx(y,z)1 arrar to place prolonged data
!!
!!   Integer, Intent(in) :: ioff,joff,koff,mype
!!     Offsets and local pe id.
!!
!!   Integer, Intent(in) :: lb,lb_p,pe_p
!!     lb is local block to prolong to and lb_p, pe_p are block and 
!!     pe of parent block
!!
!!   Integer, Intent(in) :: iface 
!!     is the face variable to prolong, 1 -> facevarx, 2 -> facevary, 
!!     3 -> facevarz
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
!!   amr_1blk_fc_prol_inject 
!!   amr_1blk_fc_prol_linear
!!   amr_1blk_fc_prol_genorder
!!   amr_1blk_fc_prol_user
!!
!! RETURNS
!!
!!   Upon return prolonged data is placed in facevarx(y,z)1.
!!
!! DESCRIPTION
!!
!!   This routine is a wrapper routine which calls the functions
!!   which prolong data for FACEVAR. The local logical array lmask can
!!   be used to control which routine is actually operating on
!!   each variable stored within FACEVAR.
!
!! AUTHORS
!!
!! Written by Peter MacNeice January 2002.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_fc_prol_gen_fun(                             & 
                                          recv,ia,ib,ja,jb,ka,kb,      & 
                                          idest,ioff,joff,koff,        & 
                                          mype,lb,pe_p,lb_p,iface)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use prolong_arrays

      Use paramesh_interfaces, only :                                  & 
                        amr_1blk_fc_prol_inject,                       & 
                        amr_1blk_fc_prol_linear,                       & 
                        amr_1blk_fc_prol_genorder,                     & 
                        amr_1blk_fc_prol_user

      Implicit None

!-----Include Statements
      include 'mpif.h'

!-----Input/Output Variables
      real,    intent(inout) :: recv(:,:,:,:)
      integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in)    :: ioff,joff,koff,mype,iface
      integer, intent(in)    :: lb,lb_p,pe_p

!-----Local variables
      Double Precision :: time1

      Integer :: interp_mask_face(nbndvar), ivar

!-----Begin Exectuable Code

      If (timing_mpi) Then
         time1 = mpi_wtime()
      End If

      If (iface == 1) Then
        interp_mask_face(:) = interp_mask_facex(:)
      Elseif (iface == 2) Then
        interp_mask_face(:) = interp_mask_facey(:)
      Elseif (iface == 3) Then
        interp_mask_face(:) = interp_mask_facez(:)
      End If

      Do ivar = 1, nbndvar

      If (any(int_gcell_on_fc(1:ndim,ivar))) Then
         
      If (interp_mask_face(ivar) < 20) Then

      If (interp_mask_face(ivar) == 0) Then
!--------Injection - zeroth order prolongation 
         Call amr_1blk_fc_prol_inject(                                 & 
              recv,ia,ib,ja,jb,ka,kb,                                  & 
              idest,ioff,joff,koff,mype,iface,ivar)
         
      Elseif (interp_mask_face(ivar) == 1) Then
!--------Tri-linear interpolation 
         Call amr_1blk_fc_prol_linear(                                 & 
              recv,ia,ib,ja,jb,ka,kb,                                  & 
              idest,ioff,joff,koff,mype,iface,ivar)
         
      Else If (interp_mask_face(ivar) > 1) Then
!--------General order
         Call amr_1blk_fc_prol_genorder(                               & 
              recv,ia,ib,ja,jb,ka,kb,                                  & 
              idest,ioff,joff,koff,mype,iface,ivar,                    & 
              interp_mask_face(ivar))

      End If  ! End If (interp_mask_face(ivar) == 0)

      Elseif (interp_mask_face(ivar) == 20) Then

!--------User defined
         Call amr_1blk_fc_prol_user()

      End If  ! End If (interp_mask_face(ivar) < 20)

      End If  ! End If (any(int_gcell_on_fc(1:ndim,ivar)))

      End Do  ! End Do ivar = 1, nbndvar


      If (timing_mpi) Then
              timer_amr_1blk_fc_prol_gen =                             & 
                                timer_amr_1blk_fc_prol_gen             & 
                                + mpi_wtime() - time1
      End If

      Return
      End Subroutine amr_1blk_fc_prol_gen_fun
