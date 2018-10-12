!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_nc_prol_genorder
!! NAME
!!
!!   amr_1blk_nc_prol_genorder
!!
!! SYNOPSIS
!!
!!    Call amr_1blk_nc_prol_genorder 
!!        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, 
!!         mype,ivar,order)
!!    Call amr_1blk_nc_prol_genorder 
!!        (real array,integer,integer,integer,integer,integer,integer,
!!          integer,integer,integer,integer, 
!!         integer,integer,integer)
!!
!! ARGUMENTS
!!
!!  Real,    Intent(inout) :: recv(:,:,:,:)                                             
!!    Data array holding the data extracted from unk which will be prolonged            
!!    and placed into the unk_n1 array.                                                 
!!                                                                                      
!!  Integer, Intent(in) :: ia,ib,ja,jb,ka,kb                                            
!!    Integers which control the limits into unk_n1 where the prolonged data            
!!    will be placed.                                                                   
!!                                                                                      
!!  Integer, Intent(in) :: idest,ioff,joff,koff,mype                                    
!!    Idest controls which 'layer' into which the prolonged data will be                
!!    placed in unk_n1.  ioff, joff and koff are offsets.  mype is is the               
!!    local processor id.                                                               
!!                                                                                      
!!  Integer, Intent(in) :: ivar,order                                                   
!!    ivar is the variable number in unk_n which is prolonged.                           
!!    order is the order of the Lagrange polynomial to use.                             
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
!!
!! CALLS
!!     
!!   No call made to other PARAMESH routines.
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
!!  This particular prolongation uses a more general interpolation proceedure then
!!  some of the other routines provided with PARAMESH.  Any 'order' if interpolation
!!  can be selected for any variable (as described below).
!!  It does this by explicitly computing the necessary Taylor expansions out 
!!  to the specified order.  The interpolations are performed first in the `x' 
!!  direction.  `Y' interapolations follow, but use the interpolated
!!  data from the `x' sweep.  The 'Z' sweep is similarly performed.  
!!
!!  To select the `order' (we use the term order here loosely) of interpolation 
!!  the array interp_mask must have data in it that is >= 0.  
!!  Since the interpolation scheme is general, one can select
!!  different orders of interpolation for different variables as.  For instance,
!!  if,
!!  interp_mask(1) = 0
!!  interp_mask(2) = 1
!!  interp_mask(3) = 2
!!  then variable 1 will be prolongated used simple direct injection, variable 2
!!  will be prolongated using linear interpolation and variable 3 will be prolongated
!!  using quadratic interpolation.
!!
!!  Finally, the `order' of interpolation must be equal or less than nguard.
!!
!!  It is applied to all UNK variables whose corresponding element
!!  of interp_mask is set to 0.
!!
!!  NOTE: This routine may not be as effcient as some of the other, similar routines
!!        provided for prolongation. So, if you don't need the flexibility 
!!        of this routine, you might want to consider using another or writing 
!!        another yourself.
!!
!!  NOTE2:  This routine does NOT guarantee conservative prologation at refinement
!!          jumps.  This is described in the documentation.
!!
!! AUTHORS
!!
!!   Written :     Kevin Olson,  March 2002 and based on similar routines 
!!                 by Peter MacNeice.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_nc_prol_genorder                             & 
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,                  & 
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
      Integer,Parameter :: maxorder = 4

      Real,Save :: weight(2,-maxorder-maxorder/2:                      & 
                             maxorder+maxorder/2,0:maxorder)
      Real :: f_intx(iu_bnd1+1,                                        & 
                     ju_bnd1+k2d,                                      & 
                     ku_bnd1+k3d)
      Real :: f_inty(iu_bnd1+1,                                        & 
                     ju_bnd1+k2d,                                      & 
                     ku_bnd1+k3d)

      Integer :: i,j,k
      Integer :: offi,offj,offk
      Integer :: ii,jj,kk,iorder
      Integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax
      Integer :: ifmin,ifmax,jfmin,jfmax,kfmin,kfmax
      Integer :: imin, imax, jmin, jmax, kmin, kmax
      Integer,Save :: iminh(2,0:maxorder),imaxh(2,0:maxorder)
      Integer :: ipar, jpar, kpar
      Integer :: iw, jw, kw
      Integer,Parameter :: largei = 100      

      Logical,Save :: first_call = .true.
      
!-----Begin Executable Code

      If (first_call) Then
         first_call = .False.

         Do iorder = 0,maxorder

!-----------HALF CELL TO THE RIGHT of CELL FACE

            iminh(1,iorder) = 0
            imaxh(1,iorder) = iorder

            Do ipar = iminh(1,iorder),imaxh(1,iorder)
               weight(1,ipar,iorder) = 1.
               Do jpar = iminh(1,iorder),imaxh(1,iorder)
                  If (jpar.ne.ipar) Then
                     weight(1,ipar,iorder) =                           & 
                          weight(1,ipar,iorder)*                       & 
                          (.5-jpar)/(ipar-jpar)
                  End If
               End Do
            End Do

            iminh(2,iorder) = 1 - iorder
            imaxh(2,iorder) = 1

            Do ipar = iminh(2,iorder),imaxh(2,iorder)
               weight(2,ipar,iorder) = 1.
               Do jpar = iminh(2,iorder),imaxh(2,iorder)
                  If (jpar.ne.ipar) Then
                     weight(2,ipar,iorder) =                           & 
                          weight(2,ipar,iorder)*                       & 
                          (.5-jpar)/(ipar-jpar)
                  End If
               End Do
            End Do

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

      kcmin = ((kfmin-nguard-1+largei)/2 +                             & 
                      nguard - largei/2 )*k3d +                        & 
                      1 + offk
      kcmax = ((kfmax-nguard-1+largei)/2 +                             & 
                      nguard - largei/2 )*k3d +                        & 
                      1 + offk
      jcmin = ((jfmin-nguard-1+largei)/2 +                             & 
                      nguard - largei/2 )*k2d +                        & 
                      1 + offj
      jcmax = ((jfmax-nguard-1+largei)/2 +                             & 
                      nguard - largei/2 )*k2d +                        & 
                      1 + offj
      icmin = ((ifmin-nguard-1+largei)/2 +                             & 
                      nguard - largei/2 ) +                            & 
                      1 + offi
      icmax = ((ifmax-nguard-1+largei)/2 +                             & 
                      nguard - largei/2 ) +                            & 
                      1 + offi





!-----Main Interpolation loop.






!-----Interpolate in x direction

      If (ndim >= 1) Then


      f_intx(:,:,:) = 0.

      kmin = kcmin-(nguard+2)*k3d
      If (kmin < 1) kmin = 1
      kmax = kcmax+(nguard+2)*k3d 
      If (kmax > nguard*2+nzb+k3d) kmax = nguard*2+nzb + k3d
      jmin = jcmin-(nguard+2)*k2d
      If (jmin < 1) jmin = 1
      jmax = jcmax+(nguard+2)*k2d
      If (jmax > nguard*2+nyb+k2d) jmax = nguard*2+nyb + k2d

      Do k = kmin,kmax
      Do j = jmin,jmax

         ! 1) now interpolate to half points

         ! starting parent index
         i = icmin
         
         Do ii = ifmin,ifmax

            If ((mod(ii,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.     & 
                (mod(ii,2)  ==  0 .and. mod(nguard,2) .ne. 0)) Then
                                       ! this point is on one of parent's points
                                       ! and Does not need to be interpolated
               f_intx(ii,j,k) = recv(ivar,i,j,k)
            Else

               If (ii < nguard + nxb/2) Then
                  iw = 1
               Else
                  iw = 2
               End If

               imin = iminh(iw,order) + i
               imax = imaxh(iw,order) + i

               Do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +                    & 
                       weight(iw,ipar-i,order)*recv(ivar,ipar,j,k)
               End Do
               ! update parent index
               i = i + 1
            End If  ! End If ((mod(ii,2) .ne. 0 ...

         End Do   ! End Do ii = ifmin, if
      End Do   ! End Do j = jmin,jmax
      End Do   ! End Do k = kmin,kmax

      If (ndim == 1) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_n1(ivar,i,j,k,idest) = f_intx(i,j,k)
               End Do
            End Do
         End Do
      End If

      End If  ! End If (ndim >= 1)

!-----Interpolate in y direction

      If (ndim >= 2) Then
      

      f_inty(:,:,:) = 0.

      kmin = kcmin-(nguard+2)*k3d
      If (kmin < 1) kmin = 1
      kmax = kcmax+(nguard+2)*k3d
      If (kmax > nguard*2+nzb + k3d) kmax = nguard*2+nzb + k3d

      Do k = kmin,kmax
      Do i = ifmin,ifmax

         ! 1) interpolate to half points
 
         ! starting parent index
         j = jcmin
            
         Do jj = jfmin,jfmax
            
            If ((mod(jj,2) .ne. 0 .and. mod(nguard,2)  == 0) .or.      & 
                (mod(jj,2)  ==  0 .and. mod(nguard,2) .ne.  0)) Then
               f_inty(i,jj,k) = f_intx(i,j,k)
            Else

               If (jj < nguard + nyb/2) Then
                  jw = 1
               Else
                  jw = 2
               End If

               jmin = iminh(jw,order) + j
               jmax = imaxh(jw,order) + j

               Do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +             & 
                    weight(jw,jpar-j,order)*f_intx(i,jpar,k)
               End Do
               ! update parent index
               j = j + 1

            End If  ! End If ((mod(jj,2) .ne. 0 ...
            
         End Do  ! End Do jj = jfmin,jfmax
      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do k = kmin,kmax

      If (ndim == 2) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_n1(ivar,i,j,k,idest) = f_inty(i,j,k)
               End Do
            End Do
         End Do
      End If

      End If  ! End If (ndim >= 2)

!-----Interpolate in z direction

      If (ndim == 3) Then


      Do j = jfmin,jfmax
      Do i = ifmin,ifmax

         ! 1) interpolate to half points

         ! starting parent index
         k = kcmin

         Do kk = kfmin,kfmax
            
            unk_n1(ivar,i,j,kk,idest) = 0.
            If ((mod(kk,2) .ne. 0 .and. mod(nguard,2)  == 0) .or.      & 
                (mod(kk,2)  ==  0 .and. mod(nguard,2) .ne.  0)) Then
               unk_n1(ivar,i,j,kk,idest) = f_inty(i,j,k)
            Else

               If (kk < nguard + nzb/2) Then
                  kw = 1
               Else
                  kw = 2
               End If

               kmin = iminh(kw,order) + k
               kmax = imaxh(kw,order) + k

               Do kpar = kmin,kmax
                  unk_n1(ivar,i,j,kk,idest) =                          & 
                       unk_n1(ivar,i,j,kk,idest) +                     & 
                       weight(kw,kpar-k,order)*f_inty(i,j,kpar)
               End Do
               ! update parent index
               k = k + 1
            End If  ! End If ((mod(kk,2) .ne. 0 ...
            
         End Do  ! End Do kk = kfmin,kfmax
      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do j = jfmin,jfmax

      End If  ! End If (ndim == 3)

      
      End Subroutine amr_1blk_nc_prol_genorder






