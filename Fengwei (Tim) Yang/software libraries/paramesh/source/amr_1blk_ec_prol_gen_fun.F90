!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_ec_prol_gen_fun
!! NAME
!!
!!   amr_1blk_ec_prol_gen_fun
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_ec_prol_gen_fun (recv, ia, ib, ja, jb, ka, kb, idest
!!                                  ioff, joff, koff, mype,
!!                                  iedge_dir)
!!   Call amr_1blk_ec_prol_gen_fun (real, 
!!                                  integer, integer, integer, integer,
!!                                  integer, integer, integer, integer,
!!                                  integer, integer, integer, integer,
!!                                  integer, integer)
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
!!   Integer, Intent(in) :: iedge_dir
!!     Edge to apply prolongation to.
!!     iedge_dir = 1 -> unk_e_x
!!     iedge_dir = 2 -> unk_e_y
!!     iedge_dir = 3 -> unk_e_z
!!
!! INCLUDES
!!   
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
!!    amr_1blk_ec_prol_linear 
!!    amr_1blk_ec_prol_genorder
!!    amr_1blk_ec_prol_user
!!
!! RETURNS
!!
!!   Upon return prolonged data is placed in unk_e_x(y,z)1.
!!
!! DESCRIPTION
!!
!!   This routine is a wrapper routine which calls the functions
!!   which actually prolong the data for UNK_E_X(Y,Z). 
!!
!! AUTHORS
!!
!!   Written by Peter MacNeice with modification be Kevin Olson for
!!   high order ploynomial interpolation.
!!
!!***

      Subroutine amr_1blk_ec_prol_gen_fun(recv,   &
             ia,ib,ja,jb,ka,kb,idest,             & 
             ioff,joff,koff,mype,iedge_dir)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use prolong_arrays

      Use paramesh_interfaces, only :             & 
           amr_1blk_ec_prol_linear,               & 
           amr_1blk_ec_prol_genorder,             & 
           amr_1blk_ec_prol_user

      Implicit None

!-----Input/Output Variables
      real,    intent(inout) :: recv(:,:,:,:)
      integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in)    :: ioff,joff,koff,mype,iedge_dir

!-----Include Statements
      Include 'mpif.h'

!-----Local Variables
      Integer          :: ivar
      Double Precision :: time1
      
!-----Begin Executable Code

      If (timing_mpi) Then
         time1 = mpi_wtime()
      End If  ! End If (timing_mpi)

      Do ivar = 1,nbndvare

      If (any(int_gcell_on_ec(1:ndim,ivar))) Then

      If (interp_mask_ec(ivar) < 20) Then

      If (interp_mask_ec(ivar) == 1) Then
!--------standard linear interpolation
         Call amr_1blk_ec_prol_linear                        & 
              (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  & 
               mype,ivar,iedge_dir)
      Else
!--------High order Lagrange polynomial interpolation
         Call amr_1blk_ec_prol_genorder                      & 
              (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  & 
               mype,ivar,iedge_dir,interp_mask_ec(ivar))
      End If  ! End If (inter_mask_ec(ivar == 1)

      Else  ! If (interp_mask_ec(ivar) < 20)

         Call amr_1blk_ec_prol_user()

      End If  ! End If (interp_mask_ec(ivar) < 20)

      End If  ! End If (any(int_gcell_on_ec(1:ndim,ivar)))

      End Do  ! End Do ivar = 1,nbndvare


      If (timing_mpi) Then
              timer_amr_1blk_ec_prol_gen =                 & 
                               timer_amr_1blk_ec_prol_gen  & 
                                + mpi_wtime() - time1
      End If

      Return
      End Subroutine amr_1blk_ec_prol_gen_fun





