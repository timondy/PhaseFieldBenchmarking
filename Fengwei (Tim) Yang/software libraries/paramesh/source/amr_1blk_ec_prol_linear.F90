!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_ec_prol_linear
!! NAME
!!  
!!   amr_1blk_ec_prol_linear
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_ec_prol_linear(recv,ia,ib,ja,jb,ka,kb,
!!                                idest,ioff,joff,koff,mype,ivar,
!!                                iedge_dir)
!!   Call amr_1blk_ec_prol_linear(real,
!!                                integer, integer, integer, integer,
!!                                integer, integer, integer, integer,
!!                                integer, integer, integer, integer,
!!                                integer)
!!
!!
!! ARGUMENTS
!!
!!  Real,    intent(inout) :: recv(:,:,:,:)
!!    Data array holding the data extracted from unk_e_x(y,z) which will be prolonged
!!    and placed into the unk_e_x(y,z)1 array.
!!
!!  Integer, Intent(in) :: ia,ib,ja,jb,ka,kb
!!    Integers which control the limits into unk_e_z(y,z)1 where the prolonged data
!!    will be placed.
!!
!!  Integer, Intent(in) :: idest,ioff,joff,koff,mype
!!    Idest controls which 'layer' into which the prolonged data will be
!!    placed in unk_e_x(y,z)1.  ioff, joff and koff are offsets.  mype is is the
!!    local processor id.
!!
!!  Integer, Intent(in) :: ivar
!!    ivar is the variable number in unk_e_x(y,z) which is prolonged.
!!
!!  Integer, Intent(in) :: iedge_dir
!!    Which edge to apply interpolation to
!!    if iedge_dir = 1, apply interpolation unk_e_x
!!    if iedge_dir = 2, apply interpolation unk_e_y
!!    if iedge_dir = 3, apply interpolation unk_e_z
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
!!   and placed into the unk_e_x(y,z)1 array.
!!
!! DESCRIPTION
!!
!!   This routine takes data from the array recv, originally extracted 
!!   from the solution array unk_e_(y,z), and performs a prolongation operation 
!!   on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
!!   The data in recv is from a parent block and the
!!   result of the prolongation operation is written directly into one
!!   layer of the working block array unk_e_x(y,z)1(...,idest).
!!   The position of the child within the parent block is specified by 
!!   the ioff, joff and koff arguments.
!!
!!   This particular prolongation is simple linear interpolation. It can
!!   be used for blocks with an even or odd number of grid cells.
!!
!!   It is applied to all UNK_E_X(Y,Z) variables whose corresponding element
!!   of interp_mask_ec is set to 1.
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

      Subroutine amr_1blk_ec_prol_linear                               & 
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,               & 
            mype,ivar,iedge_dir)

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
      Integer, Intent(in)    :: ivar, iedge_dir

!-----Local arrays and Variables
      real    :: dx,dy,dz,cx,cy,cz
      integer :: icl,icu,jcl,jcu,kcl,kcu,i_ind,j_ind,k_ind
      integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p

!-----Begin Exectuable Code

      If (prol_init.ne.100) Then
         Write(*,*) 'PARAMESH ERROR !'
         Write(*,*) 'Error : prolong_face_fun. ',                      & 
              'You must call amr_prolong_face_fun_init ',              & 
              'before you can use this routine!'
         Call amr_abort
      End If

!-----Set the bounds on the loop controlling the interpolation.
      icl=ia
      icu=ib
      jcl=ja
      jcu=jb
      kcl=ka
      kcu=kb


      i_ind = 1
      j_ind = 1
      k_ind = 1
      If (ioff > 0) i_ind = 2
      If (joff > 0) j_ind = 2
      If (koff > 0) k_ind = 2

!-----Interpolation loop.

!-----Note that the range of indeces used in the facevar plane differs
!-----depending on the value of iface_off. This assumes that the face values
!-----corresponding to index 0 (ie nguard faces to the left of the block
!-----boundary) are never needed, when iface_off=-1. 

!-----Interpolate to unk_e_x1
      If (iedge_dir == 1) Then

        Do k=kcl,kcu+iface_off
             k1 = prol_f_indexz(1,k,k_ind)
             k1p= prol_f_indexz(2,k,k_ind)
             dz = prol_f_dz(k)
             cz = 1.-dz
             Do j=jcl,jcu+iface_off
                   j1 = prol_f_indexy(1,j,j_ind)
                   j1p= prol_f_indexy(2,j,j_ind)
                   dy = prol_f_dy(j)
                   cy = 1.-dy
                   Do i=icl,icu
                         i1 = prol_indexx(1,i,i_ind)
                         i1p= prol_indexx(2,i,i_ind)
                         dx = prol_dx(i)
                         cx = 1.-dx

!------------------------compute interpolated values at location (i,j,k)
                         unk_e_x1(ivar,i,j,k,idest) =                  & 
                              dz*( dy*( dx*recv(ivar,i1,j1,k1) +       & 
                              cx*recv(ivar,i1p,j1,k1))  +              & 
                              cy*( dx*recv(ivar,i1,j1p,k1) +           & 
                              cx*recv(ivar,i1p,j1p,k1) ) ) +           & 
                              cz*( dy*( dx*recv(ivar,i1,j1,k1p) +      & 
                              cx*recv(ivar,i1p,j1,k1p))  +             & 
                              cy*( dx*recv(ivar,i1,j1p,k1p) +          & 
                              cx*recv(ivar,i1p,j1p,k1p) ) )

                   End Do  ! End Do i=icl,icu
             End Do  ! End Do j=jcl,jcu+iface_off
        End Do  ! End Do k=kcl,kcu+iface_off

!-----Interpolate to unk_e_y1
      Elseif (iedge_dir == 2) Then

        Do k=kcl,kcu+iface_off
             k1 = prol_f_indexz(1,k,k_ind)
             k1p= prol_f_indexz(2,k,k_ind)
             dz = prol_f_dz(k)
             cz = 1.-dz
             Do j=jcl,jcu
                   j1 = prol_indexy(1,j,j_ind)
                   j1p= prol_indexy(2,j,j_ind)
                   dy = prol_dy(j)
                   cy = 1.-dy
                   Do i=icl,icu+iface_off
                         i1 = prol_f_indexx(1,i,i_ind)
                         i1p= prol_f_indexx(2,i,i_ind)
                         dx = prol_f_dx(i)
                         cx = 1.-dx

!------------------------compute interpolated values at location (i,j,k)
                         unk_e_y1(ivar,i,j,k,idest) =                  & 
                              dz*( dy*( dx*recv(ivar,i1,j1,k1) +       & 
                              cx*recv(ivar,i1p,j1,k1))  +              & 
                              cy*( dx*recv(ivar,i1,j1p,k1) +           & 
                              cx*recv(ivar,i1p,j1p,k1) ) ) +           & 
                              cz*( dy*( dx*recv(ivar,i1,j1,k1p) +      & 
                              cx*recv(ivar,i1p,j1,k1p))  +             & 
                              cy*( dx*recv(ivar,i1,j1p,k1p) +          & 
                              cx*recv(ivar,i1p,j1p,k1p) ) )

                   End Do  ! End Do i=icl,icu+iface_off
              End Do  ! End Do j=jcl,jcu
        End Do  ! End Do k=kcl,kcu+kface_off

!-----Interpolate to unk_e_z1
      Elseif (iedge_dir == 3) Then

        Do k=kcl,kcu
             k1 = prol_indexz(1,k,k_ind)
             k1p= prol_indexz(2,k,k_ind)
             dz = prol_dz(k)
             cz = 1.-dz
             Do j=jcl,jcu+iface_off
                   j1 = prol_f_indexy(1,j,j_ind)
                   j1p= prol_f_indexy(2,j,j_ind)
                   dy = prol_f_dy(j)
                   cy = 1.-dy
                   Do i=icl,icu+iface_off
                         i1 = prol_f_indexx(1,i,i_ind)
                         i1p= prol_f_indexx(2,i,i_ind)
                         dx = prol_f_dx(i)
                         cx = 1.-dx

!------------------------compute interpolated values at location (i,j,k)
                         unk_e_z1(ivar,i,j,k,idest) =                  & 
                              dz*( dy*( dx*recv(ivar,i1,j1,k1) +       & 
                              cx*recv(ivar,i1p,j1,k1))  +              & 
                              cy*( dx*recv(ivar,i1,j1p,k1) +           & 
                              cx*recv(ivar,i1p,j1p,k1) ) ) +           & 
                              cz*( dy*( dx*recv(ivar,i1,j1,k1p) +      & 
                              cx*recv(ivar,i1p,j1,k1p))  +             & 
                              cy*( dx*recv(ivar,i1,j1p,k1p) +          & 
                              cx*recv(ivar,i1p,j1p,k1p) ) )

                   End Do  ! End Do i=icl,icu+iface_off
             End Do  ! End Do j=jcl,jcu+iface_off
        End Do  ! End Do k=kcl,kcu

      End If  ! End If (iedge_dir == 1)


      Return
      End Subroutine amr_1blk_ec_prol_linear
