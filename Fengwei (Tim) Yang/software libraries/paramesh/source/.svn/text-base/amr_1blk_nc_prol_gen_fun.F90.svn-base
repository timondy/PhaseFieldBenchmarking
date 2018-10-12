!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_nc_prol_gen_fun
!! NAME
!!
!!   amr_1blk_nc_prol_gen_fun
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_nc_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest, 
!!                                 ioff,joff,koff,mype)
!!   Call amr_1blk_nc_prol_gen_fun(real array, 
!!                                 integer,integer,integer,integer,integer,integer,
!!                                  integer,
!!                                 integer,integer,integer,integer)
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
!!   amr_1blk_nc_prol_linear
!!   amr_1blk_nc_prol_genorder
!!   amr_1blk_nc_prol_user
!!
!! RETURNS
!!
!!   Upon return prolonged data is placed in unk_n1.
!!
!! DESCRIPTION
!!
!!   This routine is a wrapper routine which call the functions which actually prolong
!!   the data for UNK_N.  The order of interpolation which is used is determined
!!   by the data stored in the array interp_mask_nc.  Default is linear interpolation.
!! 
!!   Note: before using this routine in your program, make sure that the
!!   routine prolong_face_fun_init has been called.
!!
!! AUTHORS
!!
!! Written by Peter MacNeice December 2000 with modifications by Kevin Olson for
!! high order interpolation.
!!
!!***

      Subroutine amr_1blk_nc_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,      &
                                          idest,ioff,joff,koff,mype)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use prolong_arrays

      Use paramesh_interfaces, Only : amr_1blk_nc_prol_linear,         & 
                                      amr_1blk_nc_prol_genorder,       & 
                                      amr_1blk_nc_prol_user

      Implicit None

!-----Input/Output variables
      Real,    Intent(inout) :: recv(:,:,:,:)
      Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      Integer, Intent(in)    :: ioff,joff,koff,mype

!-----Include statements
      include 'mpif.h'

!-----Local variables
      Integer :: ivar
      Double Precision :: time1
      
!-----Begin Executable Code
      If (timing_mpi) then
         time1 = mpi_wtime()
      End If

      Do ivar=1,nvarcorn

      If (int_gcell_on_nc(ivar)) Then

      If (interp_mask_nc(ivar) < 20) Then

      If (interp_mask_nc(ivar) == 1) Then
!--------linear interpolation
         Call amr_1blk_nc_prol_linear(recv,ia,ib,ja,jb,ka,kb,idest,    & 
             ioff,joff,koff,mype,ivar)
      Else
!--------general, high order, interpolation routine
         Call amr_1blk_nc_prol_genorder(recv,ia,ib,ja,jb,ka,kb,idest,  & 
             ioff,joff,koff,mype,ivar,interp_mask_nc(ivar))
      End If  ! End If (interp_mask_nc(ivar) == 1)

      Else
         
!--------User defined interpolation
         call amr_1blk_nc_prol_user()

      End If  !  End If (interp_mask_nc(ivar) < 20)

      End If  !  End If (int_gcell_on_nc(ivar)

      End Do  !  End Do ivar=1,nvarcorn

      If (timing_mpi) Then
              timer_amr_1blk_nc_prol_gen =                             & 
                                timer_amr_1blk_nc_prol_gen             & 
                                + mpi_wtime() - time1
      End If

      Return
      End Subroutine amr_1blk_nc_prol_gen_fun





