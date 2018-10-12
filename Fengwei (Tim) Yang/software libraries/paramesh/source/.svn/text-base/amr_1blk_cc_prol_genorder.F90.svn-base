!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_genorder
!! NAME
!!
!!   amr_1blk_cc_prol_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_genoder (recv, ia, ib, ja, jb, ka, kb,  &
!!                                  idest, ioff, joff, koff,       &
!!                                  mype, ivar, order)
!!
!!   Call amr_1blk_cc_prol_genorder (real, integer, integer, 
!!                                   integer, integer, integer, integer,
!!                                   integer, integer, integer, integer,
!!                                   integer, integer, integer)
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
!!  Integer, intent(in) :: ivar,order
!!    ivar is the varible number in unk which is prolonged.
!!    order is the order of the Lagrange polynomial to use.
!!   
!!
!! INCLUDE
!!   
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!!
!!   No calls made to other PARAMESH routines.
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit the data in recv is interpolated
!!   and placed into the unk1 array.
!!
!! DESCRIPTION
!!
!!  This routine takes data from the array recv, originally extracted
!!  from the solution array unk, and performs a prolongation operation
!!  on it, between the bounds ranges ia to ib, ja to jb, and ka to kb.
!!  The data in recv is from a parent block and the
!!  result of the prolongation operation is written directly into one
!!  layer of the working block array unk1(...,idest).
!!  The position of the child within the parent block is specified by
!!  the ioff, joff and koff arguments.
!!
!!  This particular prolongation uses a more general interpolation proceedure
!!  than some of the other routines provided with PARAMESH.  Lagrange interpolation 
!!  is used where Lagrange polynomials are fit through the available data points.
!!  Interpolation is done using Lagrange interpolation.  Any 'order' of interpolation 
!!  can be selected for any variable (as described below). 
!!  The interpolations are performed first in the 'X' direction.  `Y' 
!!  interpolations follow, but use the interpolated data from the `X' sweep.  
!!  The 'Z' sweep is similarly performed.
!!
!!  Here the variable 'order' selects the highest order of the Lagrange polynomial
!!  to use.  Since the interpolation scheme is general, one can select
!!  different orders of interpolation for different variables as.  For instance,
!!  if,
!!  interp_mask(1) = 0
!!  interp_mask(2) = 1
!!  interp_mask(3) = 2
!!  then variable 1 will be prolongated used simple direct injection, variable 2
!!  will be prolongated using linear interpolation and variable 3 will be
!!  prolongated using quadratic interpolation.
!!
!!  Finally, the `order' of interpolation must be equal or less than nguard.
!!  This ensures that enough guardcells space is available to compute
!!  the interpolation weights for the polynomial fits.
!!
!!  It is applied to all UNK variables whose corresponding element
!!  of interp_mask is set to 0.
!!
!!  Finally, the `order' of interpolation must be equal or less than nguard.
!!  This ensures that enough guardcells space is available to compute
!!  the interpolation weights for the polynomial fits.
!!
!!  It is applied to all UNK variables whose corresponding element
!!  of interp_mask is set to 0.
!!
!!  NOTE: This routine may not be as effcient as some of the other, similar
!!        routine provided for prolongation. So, if you don't need the
!!        flexibility of this routine, you might want to consider using another
!!        or writing another yourself.
!!
!!  NOTE2:  This routine does NOT guarantee conservative prologation at
!!          refinement jumps.  This is described in the documentation.
!!
!! AUTHORS
!!
!!   Written by Kevin Olson,  March 2002 and based on similar routines
!!   by Peter MacNeice.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_cc_prol_genorder            & 
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
         mype,ivar,order)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree

      Implicit None

!-----Input/Output Arguments
      Real,    Intent(inout) :: recv(:,:,:,:)
      Integer, Intent(in) :: ia,ib,ja,jb,ka,kb
      Integer, Intent(in) :: idest,ioff,joff,koff,mype
      Integer, Intent(in) :: ivar,order

!-----Local arrays and variables
!-----maxorder is the maximum order of the Lagrange polynomial allowed                   
!-----4 was chosen to make this routine more memory efficient                           
!-----but the routine will function if maxorder is changed and then                      
!-----PARAMESH is recompiled.                                                          
      Integer, Parameter :: largei = 100
      Integer, Parameter :: maxorder = 4

      Real :: weight_right
      Real :: weight_left
      Real :: tempy, tempx
      Real, Save, Allocatable    :: weightx(:,:,:)
      Real, Save, Allocatable    :: weighty(:,:,:)
      Real, Save, Allocatable    :: weightz(:,:,:)

      Integer :: i,j,k,ii
      Integer :: offi,offj,offk
      Integer :: iorder
      Integer :: icmin,jcmin,kcmin
      Integer :: ifmin,ifmax,jfmin,jfmax,kfmin,kfmax
      Integer :: ipar, jpar, kpar
      Integer, Save, Allocatable :: imina(:,:),   & 
                                    imaxa(:,:)
      Integer, Save, Allocatable :: jmina(:,:),   & 
                                    jmaxa(:,:)
      Integer, Save, Allocatable :: kmina(:,:),   & 
                                    kmaxa(:,:)
      
      Logical, Save :: first_call = .True.
 
!-----Begin Exectuable Code Section

!-----Computations done on first call to this routine
      If (first_call) Then

!--------For X direction

         Allocate (imina(iu_bnd1,0:maxorder))
         Allocate (imaxa(iu_bnd1,0:maxorder))
         Allocate (weightx(0:iu_bnd1,iu_bnd1,0:maxorder))

         Allocate (jmina(ju_bnd1,0:maxorder))
         Allocate (jmaxa(ju_bnd1,0:maxorder))
         Allocate (weighty(0:ju_bnd1,ju_bnd1,0:maxorder))

         Allocate (kmina(ku_bnd1,0:maxorder))
         Allocate (kmaxa(ku_bnd1,0:maxorder))
         Allocate (weightz(0:ku_bnd1,ku_bnd1,0:maxorder))

         first_call = .False.

         Do iorder = 0,maxorder
         
         i = ((1-nguard-1+largei)/2 +                               & 
               nguard - largei/2 ) + 1 
         Do ii = 1,iu_bnd1
            
            If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.  & 
                (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) Then
               
!--------------right point
               
               If (ii < nguard + nxb/2 + 1) Then
                  imina(ii,iorder) = i
                  imaxa(ii,iorder) = i + iorder
               Else
                  imina(ii,iorder) = i - iorder + 1
                  imaxa(ii,iorder) = i + 1
               End If  ! End If (ii < nguard + nxb/2 + 1)

               Do ipar = imina(ii,iorder),imaxa(ii,iorder)
                  weight_right = 1.
                  Do jpar = imina(ii,iorder),imaxa(ii,iorder)
                     If (jpar.ne.ipar) Then
                        weight_right =                              & 
                         weight_right*(.25-(jpar-i))/(ipar-jpar)
                     End If  ! End If (jpar.ne.ipar)
                  End Do  ! End Do jpar = imina(ii,iorder),imaxa(ii,iorder)
                  weightx(ipar,ii,iorder) = weight_right
               End Do  ! End Do ipar = imina(ii,iorder),imaxa(ii,iorder)
!--------------update parent index
               i = i + 1
               
            Else

!--------------left point
               
               If (ii < nguard + nxb/2 + 1) Then
                  imina(ii,iorder) = i - 1
                  imaxa(ii,iorder) = i - 1 + iorder
               Else
                  imina(ii,iorder) = i - iorder
                  imaxa(ii,iorder) = i
               End If  ! End If (ii < nguard + nxb/2 + 1)

               Do ipar = imina(ii,iorder),imaxa(ii,iorder)
                  weight_left = 1.
                  Do jpar = imina(ii,iorder),imaxa(ii,iorder)
                     If (jpar.ne.ipar) Then
                        weight_left =                                & 
                         weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     End If  ! End If (jpar.ne.ipar)
                  End Do  ! End Do jpar = imina(ii,iorder),imaxa(ii,iorder)
                  weightx(ipar,ii,iorder) = weight_left
               End Do  ! End Do ipar = imina(ii,iorder),imaxa(ii,iorder)
            
            End If  ! End If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) 
                    !          .or.  
                    !         (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0))

         End Do  ! Do ii = 1,iu_bnd1

!--------For Y Direction

         If (ndim >= 2) Then

         i = ((1-nguard-1+largei)/2 +                               & 
               nguard - largei/2 ) + 1 
         Do ii = 1,ju_bnd1
            
            If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.  & 
                (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) Then
               
!--------------right point
               
               If (ii < nguard + nyb/2 + 1) Then
                  jmina(ii,iorder) = i
                  jmaxa(ii,iorder) = i + iorder
               Else
                  jmina(ii,iorder) = i - iorder + 1
                  jmaxa(ii,iorder) = i + 1
               End If  ! End  If (ii < nguard + nyb/2 + 1)

               Do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weight_right = 1.
                  Do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                     If (jpar.ne.ipar) Then
                        weight_right =                               & 
                         weight_right*(.25-(jpar-i))/(ipar-jpar)
                     End If  ! End If (jpar.ne.ipar)
                  End Do  ! End Do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weighty(ipar,ii,iorder) = weight_right
               End Do  ! End Do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
                                ! update parent index
               i = i + 1
               
            Else

!--------------left point
               
               If (ii < nguard + nyb/2 + 1) Then
                  jmina(ii,iorder) = i - 1
                  jmaxa(ii,iorder) = i - 1 + iorder
               Else
                  jmina(ii,iorder) = i - iorder
                  jmaxa(ii,iorder) = i
               End If  ! End If (ii < nguard + nyb/2 + 1)

               Do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weight_left = 1.
                  Do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                     If (jpar.ne.ipar) Then
                        weight_left =                                & 
                         weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     End If  ! End If (jpar.ne.ipar)
                  End Do  ! End Do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weighty(ipar,ii,iorder) = weight_left
               End Do  ! End Do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
            
            End If  ! End If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) 
                    !          .or.
                    !         (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0))

         End Do  ! Do ii = 1,ju_bnd1

         End If  ! If (ndim >= 2)

!--------For Z Direction

         If (ndim == 3) Then

         i = ((1-nguard-1+largei)/2 +                                & 
               nguard - largei/2 ) + 1 
         Do ii = 1,ku_bnd1
            
            If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.   & 
                (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) Then
               
!--------------right point
               
               If (ii < nguard + nzb/2 + 1) Then
                  kmina(ii,iorder) = i
                  kmaxa(ii,iorder) = i + iorder
               Else
                  kmina(ii,iorder) = i - iorder + 1
                  kmaxa(ii,iorder) = i + 1
               End If  ! End If (ii < nguard + nzb/2 + 1)

               Do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weight_right = 1.
                  Do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                     If (jpar.ne.ipar) Then
                        weight_right =                               & 
                         weight_right*(.25-(jpar-i))/(ipar-jpar)
                     End If  ! End If (jpar.ne.ipar)
                  End Do  ! End Do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weightz(ipar,ii,iorder) = weight_right
               End Do  ! Do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
                                ! update parent index
               i = i + 1
               
            Else

!--------------left point
               
               If (ii < nguard + nzb/2 + 1) Then
                  kmina(ii,iorder) = i - 1
                  kmaxa(ii,iorder) = i - 1 + iorder
               Else
                  kmina(ii,iorder) = i - iorder
                  kmaxa(ii,iorder) = i
               End If  ! End If (ii < nguard + nzb/2 + 1)

               Do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weight_left = 1.
                  Do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                     If (jpar.ne.ipar) Then
                        weight_left =                                & 
                         weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     End If  ! End If (jpar.ne.ipar)
                  End Do  ! End Do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weightz(ipar,ii,iorder) = weight_left
               End Do  ! End Do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
            
            End If  ! End If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) 
                    !          .or.  
                    !         (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0))

         End Do  ! End Do ii = 1,ku_bnd1

         End If  ! End If (ndim == 3)

         End Do  ! End Do iorder = 0,maxorder

      End If  ! End If (first_call)



!-----Set the bounds on the loop controlling the interpolation.

      ifmin=ia
      ifmax=ib
      jfmin=ja
      jfmax=jb
      kfmin=ka
      kfmax=kb


      offi = 0
      offj = 0
      offk = 0
      If (ioff > 0) offi = nxb/2
      If (joff > 0) offj = nyb*k2d/2
      If (koff > 0) offk = nzb*k3d/2

      kcmin = ((kfmin-nguard-1+largei)/2 +       & 
                      nguard - largei/2 )*k3d +  & 
                      1 + offk
      jcmin = ((jfmin-nguard-1+largei)/2 +       & 
                      nguard - largei/2 )*k2d +  & 
                      1 + offj
      icmin = ((ifmin-nguard-1+largei)/2 +       & 
                      nguard - largei/2 ) +      & 
                      1 + offi


!-----Main Interpolation loops.

      Do k = kfmin,kfmax
      Do j = jfmin,jfmax
      Do i = ifmin,ifmax
         
         unk1(ivar,i,j,k,idest) = 0.

         If (ndim == 3) Then
               
            Do kpar = kmina(k,order),kmaxa(k,order)
            tempy = 0.
            Do jpar = jmina(j,order),jmaxa(j,order)
            tempx = 0.
            Do ipar = imina(i,order),imaxa(i,order)
               tempx = tempx +                                         & 
                    weightx(ipar,i,order)*                             & 
                    recv(ivar,ipar+offi,jpar+offj,kpar+offk)
            End Do  ! End Do ipar = imina(i,order),imaxa(i,order)
               tempy = tempy +                                         & 
                    weighty(jpar,j,order)*tempx
            End Do  ! End Do jpar = jmina(j,order),jmaxa(j,order)
               unk1(ivar,i,j,k,idest) = unk1(ivar,i,j,k,idest) +       & 
                    weightz(kpar,k,order)*tempy
            End Do  ! End Do kpar = kmina(k,order),kmaxa(k,order)

         Elseif (ndim == 2) Then

            kpar = 1
            Do jpar = jmina(j,order),jmaxa(j,order)
            Do ipar = imina(i,order),imaxa(i,order)
               unk1(ivar,i,j,k,idest) = unk1(ivar,i,j,k,idest) +       & 
                    weightx(ipar,i,order)*                             & 
                    weighty(jpar,j,order)*                             & 
                    recv(ivar,ipar+offi,jpar+offj,kpar+offk)
            End Do  ! End Do ipar = imina(i,order),imaxa(i,order)
            End Do  ! End Do jpar = jmina(j,order),jmaxa(j,order)

         Elseif (ndim == 1) Then

            kpar = 1
            jpar = 1
            Do ipar = imina(i,order),imaxa(i,order)
               unk1(ivar,i,j,k,idest) = unk1(ivar,i,j,k,idest) +       & 
                    weightx(ipar,i,order)*                             & 
                    recv(ivar,ipar+offi,jpar+offj,kpar+offk)
            End Do  ! End Do ipar = imina(i,order),imaxa(i,order)

         End If  ! End If (ndim == 3)

      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do j = jfmin,jfmax
      End Do  ! End Do k = kfmin,kfmax


      Return
      End Subroutine amr_1blk_cc_prol_genorder

