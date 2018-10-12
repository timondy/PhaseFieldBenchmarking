!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_linear
!! NAME
!!  
!!   amr_1blk_cc_prol_linear
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_linear(recv,ia,ib,ja,jb,ka,kb,
!!                                idest,ioff,joff,koff,mype,ivar)
!!   Call amr_1blk_cc_prol_linear(real,
!!                                integer, integer, integer, integer,
!!                                integer, integer, integer, integer,
!!                                integer, integer, integer, integer)
!!
!!
!! ARGUMENTS
!!
!!  Real,    intent(inout) :: recv(:,:,:,:)
!!    Data array holding the data extracted from unk which will be prolonged
!!    and placed into the unk1 array.
!!
!!  Integer, Intent(in) :: ia,ib,ja,jb,ka,kb
!!    Integers which control the limits into unk1 where the prolonged data
!!    will be placed.
!!
!!  Integer, Intent(in) :: idest,ioff,joff,koff,mype
!!    Idest controls which 'layer' into which the prolonged data will be
!!    placed in unk1.  ioff, joff and koff are offsets.  mype is is the
!!    local processor id.
!!
!!  Integer, intent(in) :: ivar
!!    ivar is the variable number in unk which is prolonged.
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
!!   prolong_arrays
!!   
!! CALLS
!! 
!!   amr_abort
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit the data in recv is interpolated
!!   and placed into the unk1 array.
!!
!! DESCRIPTION
!!
!!   This routine takes data from the array recv, originally extracted 
!!   from the solution array unk, and performs a prolongation operation 
!!   on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
!!   The data in recv is from a parent block and the
!!   result of the prolongation operation is written directly into one
!!   layer of the working block array unk1(...,idest).
!!   The position of the child within the parent block is specified by 
!!   the ioff, joff and koff arguments.
!!
!!   This particular prolongation is simple linear interpolation. It can
!!   be used for blocks with an even or odd number of grid cells.
!!
!!   It is applied to all UNK variables whose corresponding element
!!   of interp_mask is set to 1.
!!
!!   Conservative prolongation. Special treatment for the  cells immediately
!!   adjacent to a boundary 
!!   (ie i=nguard,nguard+1,iu_bnd1-nguard,iu_bnd1-nguard+1
!!   and likewise for j and k indeces) 
!!   if using an even number of grid cells per block along that axis. 
!!   No special treatment is required when the number of cells is odd.
!!
!!   Note: before using this routine in your program, make sure that the
!!   routine prolong_unk_fun_init has been called.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 1997
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_cc_prol_linear               & 
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  & 
         mype,ivar)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use prolong_arrays

      Implicit None

!-----Input/Output Variables
      Real,    Intent(inout) :: recv(:,:,:,:)
      Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb
      Integer, Intent(in)    :: idest,ioff,joff,koff,mype
      Integer, Intent(in)    :: ivar

!-----Local variables
      Real    :: dx,dy,dz,cx,cy,cz
      Integer :: icl,icu,jcl,jcu,kcl,kcu,i_ind,j_ind,k_ind
      Integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p

!-----Begin exectable code

      If (prol_init.ne.100) Then
       Write(*,*) 'PARAMESH ERROR !'
       Write(*,*) 'Error : prolong_gen_unk_fun. ',  & 
             'You must call amr_initialize ',       & 
             'before you can use this routine!'
       Call amr_abort
      End If  ! End If (prol_init.ne.100)


! Set the bounds on the loop controlling the interpolation.
      icl=ia
      icu=ib
      jcl=ja
      jcu=jb
      kcl=ka
      kcu=kb

      i_ind = 1
      j_ind = 1
      k_ind = 1
      If (ioff.gt.0) i_ind = 2
      If (joff.gt.0) j_ind = 2
      If (koff.gt.0) k_ind = 2

!-----Interpolation loop.
      Do k=kcl,kcu
           k1 = prol_indexz(1,k,k_ind)
           k1p= prol_indexz(2,k,k_ind)
           dz = prol_dz(k)
           cz = 1.-dz
           Do j=jcl,jcu
                 j1 = prol_indexy(1,j,j_ind)
                 j1p= prol_indexy(2,j,j_ind)
                 dy = prol_dy(j)
                 cy = 1.-dy
                 Do i=icl,icu
                       i1 = prol_indexx(1,i,i_ind)
                       i1p= prol_indexx(2,i,i_ind)
                       dx = prol_dx(i)
                       cx = 1.-dx

!----------------------compute interpolated values at location (i,j,k)
                       unk1(ivar,i,j,k,idest) =                    & 
                              dz*( dy*( dx*recv(ivar,i1,j1,k1) +   & 
                              cx*recv(ivar,i1p,j1,k1))  +          & 
                              cy*( dx*recv(ivar,i1,j1p,k1) +       & 
                              cx*recv(ivar,i1p,j1p,k1) ) ) +       & 
                              cz*( dy*( dx*recv(ivar,i1,j1,k1p) +  & 
                              cx*recv(ivar,i1p,j1,k1p))  +         & 
                              cy*( dx*recv(ivar,i1,j1p,k1p) +      & 
                              cx*recv(ivar,i1p,j1p,k1p) ) )


                  End Do  ! End Do i=icl,icu
           End Do  ! End Do j=jcl,jcu
      End Do  ! End Do k=kcl,kcu


      Return
      End Subroutine amr_1blk_cc_prol_linear
