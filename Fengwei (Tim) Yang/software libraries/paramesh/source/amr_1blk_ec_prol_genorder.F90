!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_ec_prol_genorder
!! NAME
!!
!!   amr_1blk_ec_prol_genorder
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_ec_prol_genoder (recv, ia, ib, ja, jb, ka, kb,  &
!!                                  idest, ioff, joff, koff,       &
!!                                  mype, ivar, iedge_dir, order)
!!
!!   Call amr_1blk_cc_prol_genorder (real, integer, integer, 
!!                                   integer, integer, integer, integer,
!!                                   integer, integer, integer, integer,
!!                                   integer, integer, integer, integer)
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
!!  Integer, intent(in) :: ivar,iedge_dir,order
!!    ivar is the varible number in unk which is prolonged.
!!    iedge_dir is the coordinate direction interpolation is
!!    being applied to (1 -> x, 2 -> y, 3 -> z).
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
!!   and placed into the unk_e_x(y,z)1 array.
!!
!! DESCRIPTION
!!
!!  This routine takes data from the array recv, originally extracted 
!!  from the solution array(s) unk_e_x(y,z), and performs a prolongation  
!!  operation on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
!!  The data in recv is from a parent block and the result of the prolongation 
!!  operation is written directly into one layer of the working block array 
!!  unk_e_x(y,z)1(...,idest).  The position of the child within the parent block 
!!  is specified by the ioff, joff and koff arguments.
!!
!!  This particular prolongation uses a more general interpolation proceedure than
!!  some of the other routines provided with PARAMESH.  Lagrange interpolation is used
!!  where Lagrange polynomials are fit through the available data points.
!!  Any 'order' of interpolation can be selected for any variable (as described below).
!!  The interpolations are performed first in the `x' direction.  `Y' 
!!  interapolations follow, but use the interpolated data from the `x' sweep.  
!!  The 'Z' sweep is similarly performed.  
!!
!!  Here the variable 'order' selects the highest order of the Lagrange polynomial
!!  to use.  Since the interpolation scheme is general, one can select
!!  different orders of interpolation for different variables as.  For instance,
!!  if,
!!  interp_mask_ec(1) = 0
!!  interp_mask_ec(2) = 1
!!  interp_mask_ec(3) = 2
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
!! Written :     Kevin Olson,  March 2002 and based on similar routines 
!!               by Peter MacNeice.
!!***


!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_ec_prol_genorder             & 
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  & 
         mype,ivar,iedge_dir,order)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree

      Implicit None

!-----Input/Output Variables
      Real,    Intent(inout) :: recv(:,:,:,:)
      Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb
      Integer, Intent(in)    :: idest,ioff,joff,koff,mype
      Integer, Intent(in)    :: ivar,iedge_dir,order


!-----Local Variables
!-----maxorder is the maximum order of the Lagrange polynomial allowed
!-----4 was chosen to make this routine more memory efficient
!-----but the routine will function if maxorder is changed and then
!-----PARAMESH is recompiled.
      Integer, Parameter :: maxorder = 4
      Integer, Parameter :: largei = 100

      Real, Save :: weight_right(2,-maxorder-maxorder/2:           & 
                                   maxorder+maxorder/2,0:maxorder)
      Real, Save :: weight_left(2,-maxorder-maxorder/2:            & 
                                   maxorder+maxorder/2,0:maxorder)
      Real, Save :: weight_half(2,-maxorder-maxorder/2:            &   
                                   maxorder+maxorder/2,0:maxorder)
      Real :: f_intx(iu_bnd1+1,                                    & 
                     ju_bnd1+k2d,                                  & 
                     ku_bnd1+k3d)
      Real :: f_inty(iu_bnd1+1,                                    &  
                     ju_bnd1+k2d,                                  & 
                     ku_bnd1+k3d)
      Real :: f_intz(iu_bnd1+1,                                    & 
                     ju_bnd1+k2d,                                  & 
                     ku_bnd1+k3d)

      Integer :: i,j,k
      Integer :: offi,offj,offk
      Integer :: ii,jj,kk,iorder
      Integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax
      Integer :: ifmin,ifmax,jfmin,jfmax,kfmin,kfmax
      Integer :: ipar, jpar, kpar
      Integer :: imin, jmin, kmin, imax, jmax, kmax
      Integer :: iw, jw, kw
      Integer, Save :: iminl(2,0:maxorder)
      Integer, Save :: imaxl(2,0:maxorder)
      Integer, Save :: iminr(2,0:maxorder)
      Integer, Save :: imaxr(2,0:maxorder)
      Integer, Save :: iminh(2,0:maxorder)
      Integer, Save :: imaxh(2,0:maxorder)

      Logical, Save :: first_call = .true.
      
!-----Begin Executable Code

      If (first_call) Then
         first_call = .False.

         Do iorder = 0,maxorder

!-----------LEFT

            iminl(1,iorder) = -1
            imaxl(1,iorder) = -1 + iorder

            Do ipar = iminl(1,iorder),imaxl(1,iorder)
               weight_left(1,ipar,iorder) = 1.
               Do jpar = iminl(1,iorder),imaxl(1,iorder)
                  If (jpar.ne.ipar) Then
                     weight_left(1,ipar,iorder) =         & 
                          weight_left(1,ipar,iorder)*     & 
                          (-.25-jpar)/(ipar-jpar)
                  End If  ! End If (jpar.ne.ipar)
               End Do  ! End Do jpar = iminl(1,iorder),imaxl(1,iorder)
            End Do  ! End Do ipar = imin1(1,iorder),imaxl(1,iorder)

            iminl(2,iorder) = -iorder
            imaxl(2,iorder) = 0

            Do ipar = iminl(2,iorder),imaxl(2,iorder)
               weight_left(2,ipar,iorder) = 1.
               Do jpar = iminl(2,iorder),imaxl(2,iorder)
                  If (jpar.ne.ipar) Then
                     weight_left(2,ipar,iorder) =         & 
                          weight_left(2,ipar,iorder)*     & 
                          (-.25-jpar)/(ipar-jpar)
                  End If  ! End If (jpar.ne.ipar)
               End Do  ! End Do jpar = iminl(2,iorder),imaxl(2,iorder)
            End Do  ! End Do ipar = iminl(2,iorder),imaxl(2,iorder)

!-----------RIGHT

            iminr(1,iorder) = 0
            imaxr(1,iorder) = iorder

            Do ipar = iminr(1,iorder),imaxr(1,iorder)
               weight_right(1,ipar,iorder) = 1.
               Do jpar = iminr(1,iorder),imaxr(1,iorder)
                  If (jpar.ne.ipar) Then
                     weight_right(1,ipar,iorder) =        & 
                          weight_right(1,ipar,iorder)*    & 
                          (.25-jpar)/(ipar-jpar)
                  End If  ! End If (jpar.ne.ipar)
               End Do  ! End Do jpar = iminr(1,iorder),imaxr(1,iorder)
            End Do  ! End Do ipar = iminr(1,iorder),imaxr(1,iorder)

            iminr(2,iorder) = 1 - iorder
            imaxr(2,iorder) = 1

            Do ipar = iminr(2,iorder),imaxr(2,iorder)
               weight_right(2,ipar,iorder) = 1.
               Do jpar = iminr(2,iorder),imaxr(2,iorder)
                  If (jpar.ne.ipar) then
                     weight_right(2,ipar,iorder) =        & 
                          weight_right(2,ipar,iorder)*    & 
                          (.25-jpar)/(ipar-jpar)
                  End If  ! End If (jpar.ne.ipar)
               End Do  ! End Do jpar = iminr(2,iorder),imaxr(2,iorder)
            End Do  ! End Do ipar = iminr(2,iorder),imaxr(2,iorder)

!-----------HALF CELL TO THE RIGHT of CELL FACE

            iminh(1,iorder) = 0
            imaxh(1,iorder) = iorder

            Do ipar = iminh(1,iorder),imaxh(1,iorder)
               weight_half(1,ipar,iorder) = 1.
               Do jpar = iminh(1,iorder),imaxh(1,iorder)
                  If (jpar.ne.ipar) Then
                     weight_half(1,ipar,iorder) =        & 
                          weight_half(1,ipar,iorder)*    & 
                          (.5-jpar)/(ipar-jpar)
                  End If  ! End If (jpar.ne.ipar)
               End Do  ! End Do jpar = iminh(1,iorder),imaxh(1,iorder)
            End Do  ! End Do ipar = iminh(1,iorder),imaxh(1,iorder)

            iminh(2,iorder) = 1 - iorder
            imaxh(2,iorder) = 1

            Do ipar = iminh(2,iorder),imaxh(2,iorder)
               weight_half(2,ipar,iorder) = 1.
               Do jpar = iminh(2,iorder),imaxh(2,iorder)
                  If (jpar.ne.ipar) Then
                     weight_half(2,ipar,iorder) =        & 
                          weight_half(2,ipar,iorder)*    & 
                          (.5-jpar)/(ipar-jpar)
                  End If  ! End If (jpar.ne.ipar)
               End Do  ! End Do jpar = iminh(2,iorder),imaxh(2,iorder)
            End Do  ! End Do ipar = iminh(2,iorder),imaxh(2,iorder)

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
      if (ioff > 0) offi = nxb/2
      if (joff > 0) offj = nyb*k2d/2
      if (koff > 0) offk = nzb*k3d/2

      kcmin = ((kfmin-nguard-1+largei)/2 +       & 
                      nguard - largei/2 )*k3d +  & 
                      1 + offk
      kcmax = ((kfmax-nguard-1+largei)/2 +       & 
                      nguard - largei/2 )*k3d +  & 
                      1 + offk
      jcmin = ((jfmin-nguard-1+largei)/2 +       & 
                      nguard - largei/2 )*k2d +  & 
                      1 + offj
      jcmax = ((jfmax-nguard-1+largei)/2 +       & 
                      nguard - largei/2 )*k2d +  & 
                      1 + offj
      icmin = ((ifmin-nguard-1+largei)/2 +       &    
                      nguard - largei/2 ) +      & 
                      1 + offi
      icmax = ((ifmax-nguard-1+largei)/2 +       & 
                      nguard - largei/2 ) +      & 
                      1 + offi





!-----Main Interpolation loops.






!-----Start interpolation of unk_e_y


      If (iedge_dir == 2) Then 


!-----Interpolate unk_e_y in x direction


      if (ndim >= 1) then


      f_intx(:,:,:) = 0.

      kmin = kcmin-(nguard+2)*k3d
      If (kmin < 1) kmin = 1
      kmax = kcmax+(nguard+2)*k3d + k3d
      If (kmax > nguard*2+nzb +k3d) kmax = nguard*2+nzb +k3d
      jmin = jcmin-(nguard+2)*k2d
      If (jmin < 1) jmin = 1
      jmax = jcmax+(nguard+2)*k2d
      If (jmax > nguard*2+nyb) jmax = nguard*2+nyb

      Do k = kmin,kmax
      Do j = jmin,jmax

         ! 1) now interpolate to half points

         ! starting parent index
         i = icmin
         
         Do ii = ifmin,ifmax

            If ((mod(ii,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.   & 
                (mod(ii,2)  ==  0 .and. mod(nguard,2) .ne. 0)) Then
               ! this point is on parent's face
               ! and does not need to be interpolated
               f_intx(ii,j,k) = recv(ivar,i,j,k)
            Else
               If (ii < nguard + nxb/2) then
                  iw = 1
               Else
                  iw = 2
               End If  ! End If (ii < nguard + nxb/2)

               imin = i + iminh(iw,order)
               imax = i + imaxh(iw,order)

               Do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +                    & 
                     weight_half(iw,ipar-i,order)*recv(ivar,ipar,j,k)
               End Do
               ! update parent index
               i = i + 1
            End If  ! End If ((mod(ii,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.
                    !         (mod(ii,2)  ==  0 .and. mod(nguard,2) .ne. 0))

         End Do  ! End Do ii = ifmin,ifmax
      End Do  ! End Do j = jmin,jmax
      End Do  ! End Do k = kmin,kmax

      If (ndim == 1) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_e_y1(ivar,i,j,k,idest) = f_intx(i,j,k)
               End Do  ! End Do i = ifmin,ifmax
            End Do  ! End Do j = jfmin,jfmax
         End Do  ! End Do k = kfmin,kfmax
      End If  ! End If (ndim == 1)

      End If  ! End If (ndim >= 1)


!-----Interpolate unk_e_y in y direction


      If (ndim >= 2) Then
      

      f_inty(:,:,:) = 0.

      kmin = kcmin-(nguard+2)*k3d
      If (kmin < 1) kmin = 1
      kmax = kcmax+(nguard+2)*k3d + k3d
      If (kmax > nguard*2+nzb+k3d) kmax = nguard*2+nzb + k3d

      Do k = kmin,kmax
      Do i = ifmin,ifmax

         ! 1) interpolate to half points

         ! starting parent index
         j = jcmin
            
         Do jj = jfmin,jfmax
            
            If ((mod(jj,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.   & 
                (mod(jj,2)  ==  0 .and. mod(nguard,2)  ==  0)) Then

                                ! right point

               If (jj < nguard + nyb/2) Then
                  jw = 1
               Else
                  jw = 2
               End If  ! End If (jj < nguard + nyb/2)

               jmin = iminr(jw,order) + j
               jmax = imaxr(jw,order) + j

               Do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +                   & 
                       weight_right(jw,jpar-j,order)*f_intx(i,jpar,k)
               End Do  ! End Do jpar = jmin,jmax
               ! update parent index
               j = j + 1

            Else

                                ! left point

               If (jj < nguard + nyb/2) Then
                  jw = 1
               Else
                  jw = 2
               End If  ! End If (jj < nguard + nyb/2)
               
               jmin = iminl(jw,order) + j
               jmax = imaxl(jw,order) + j
               
               Do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +                  & 
                       weight_left(jw,jpar-j,order)*f_intx(i,jpar,k)
               End Do  ! End Do jpar = jmin,jmax

            End If  ! End If ((mod(jj,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
                    !         (mod(jj,2)  ==  0 .and. mod(nguard,2)  ==  0))
            
         End Do  ! End Do jj = jfmin,jfmax
      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do k = kmin,kmax

      If (ndim == 2) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_e_y1(ivar,i,j,k,idest) = f_inty(i,j,k)
               End Do  ! End Do i = ifmin,ifmax
            End Do  ! End Do j = jfmin,jfmax
         End Do  ! End Do k = kfmin,kfmax
      End If  ! End If (ndim == 2)

      End If  ! End If (ndim >= 2)

      
!-----Interpolate unk_e_y in z direction


      If (ndim == 3) Then


      Do j = jfmin,jfmax
      Do i = ifmin,ifmax

         ! 1) interpolate to half points

         ! starting parent index
         k = kcmin

         Do kk = kfmin,kfmax
            
            unk_e_y1(ivar,i,j,kk,idest) = 0.
            
            If ((mod(kk,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.   & 
                (mod(kk,2)  ==  0 .and. mod(nguard,2) .ne. 0)) Then
               ! this point is on parent's face
               ! and does not need to be interpolated
               unk_e_y1(ivar,i,j,kk,idest) = f_inty(i,j,k)
            Else
               If (kk < nguard + nzb/2) Then
                  kw = 1
               Else
                  kw = 2
               End If

               kmin = iminh(kw,order) + k
               kmax = imaxh(kw,order) + k

               Do kpar = kmin,kmax
                  unk_e_y1(ivar,i,j,kk,idest) =                      & 
                       unk_e_y1(ivar,i,j,kk,idest) +                 & 
                       weight_half(kw,kpar-k,order)*f_inty(i,j,kpar)
               End Do  ! End Do kpar = kmin,kmax
               ! update parent index
               k = k + 1
            End If  ! End If ((mod(kk,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or. 
                    !         (mod(kk,2)  ==  0 .and. mod(nguard,2) .ne. 0))

         End Do  ! End Do kk = kfmin,kfmax
      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do j = jfmin,jfmax

      End If  ! End If (ndim == 3)


!-----Start interpolation of unk_e_x


      Elseif (iedge_dir == 1) Then 


!-----Interpolate unk_e_x in y direction


      If (ndim >= 1) Then


      f_inty(:,:,:) = 0.

      kmin = kcmin-(nguard+2)*k3d
      If (kmin < 1) kmin = 1
      kmax = kcmax+(nguard+2)*k3d + k3d
      If (kmax > nguard*2+nzb+k3d) kmax = nguard*2+nzb+k3d
      imin = icmin-(nguard+2)
      If (imin < 1) imin = 1
      imax = icmax+(nguard+2)
      If (imax > nguard*2+nxb) imax = nguard*2+nxb

      Do k = kmin,kmax
      Do i = imin,imax

         ! 1) interpolate to half points

         ! starting parent index
         j = jcmin
         
         Do jj = jfmin,jfmax

            If ((mod(jj,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.  & 
                (mod(jj,2)  ==  0 .and. mod(nguard,2) .ne. 0)) Then
               ! this point is on the parent's face
               ! and does not need to be interpolated
               f_inty(i,jj,k) = recv(ivar,i,j,k)
            Else
               If (jj < nguard + nyb/2) Then
                  jw = 1
               Else
                  jw = 2
               End If  ! End If (jj < nguard + nyb/2)

               jmin = iminh(jw,order) + j
               jmax = imaxh(jw,order) + j

               Do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +                    & 
                      weight_half(jw,jpar-j,order)*recv(ivar,i,jpar,k)
               End Do  ! End Do jpar = jmin,jmax
               ! update parent index
               j = j + 1
            End If  ! End If ((mod(jj,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.
                    !         (mod(jj,2)  ==  0 .and. mod(nguard,2) .ne. 0))

         End Do  ! End Do jj = jfmin,jfmax
      End Do  ! End Do i = imin,imax
      End Do  ! End Do k = kmin,kmax

      If (ndim == 1) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_e_x1(ivar,i,j,k,idest) = f_inty(i,j,k)
               End Do  ! End Do i = ifmin,ifmax
            End Do  ! End Do j = jfmin,jfmax
         End Do  ! End Do k = kfmin,kfmax
      End If  ! End If (ndim == 1)

      End If  ! End If (ndim >= 1)


!-----Interpolate unk_e_x in x direction


      If (ndim >= 2) Then
      

      f_intx(:,:,:) = 0.

      kmin = kcmin-(nguard+2)*k3d
      If (kmin < 1) kmin = 1
      kmax = kcmax+(nguard+2)*k3d + k3d
      If (kmax > nguard*2+nzb+k3d) kmax = nguard*2+nzb+k3d

      Do k = kmin,kmax
      Do j = jfmin,jfmax

         ! 1) interpolate to half points

         ! starting parent index
         i = icmin
            
         Do ii = ifmin,ifmax
            
            If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.   & 
                (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) Then

                                ! right point

               If (ii < nguard + nxb/2) Then
                  iw = 1
               Else
                  iw = 2
               End If  ! End If (ii < nguard + nxb/2)

               imin = iminr(iw,order) + i
               imax = imaxr(iw,order) + i

               Do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +                  & 
                       weight_right(iw,ipar-i,order)*f_inty(ipar,j,k)
               End Do  ! End Do ipar = imin,imax
               ! update parent index
               i = i + 1

            Else

                                ! left point

               If (ii < nguard + nxb/2) Then
                  iw = 1
               Else
                  iw = 2
               End If  ! End If (ii < nguard + nxb/2)

               imin = iminl(iw,order) + i
               imax = imaxl(iw,order) + i

               Do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +                 & 
                       weight_left(iw,ipar-i,order)*f_inty(ipar,j,k)
               End Do  ! End Do ipar = imin,imax

            End If  ! End If ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
                    !         (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0))
            
         End Do  ! End Do ii = ifmin,ifmax
      End Do  ! End Do j = jfmin,jfmax
      End Do  ! End Do k = kmin,kmax

      If (ndim == 2) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_e_x1(ivar,i,j,k,idest) = f_intx(i,j,k)
               End Do  ! End Do i = ifmin,ifmax
            End Do  ! End Do j = jfmin,jfmax
         End Do  ! End Do k = kfmin,kfmax
      End If  ! End If (ndim == 2)

      End If  ! End If (ndim >= 2)

      
!-----Interpolate unk_e_x in z direction


      If (ndim == 3) Then


      Do j = jfmin,jfmax
      Do i = ifmin,ifmax

         ! 2) now interpolate to half points

         ! starting parent index
         k = kcmin
         
         Do kk = kfmin,kfmax

            unk_e_x1(ivar,i,j,kk,idest) = 0.

            If ((mod(kk,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.  & 
                (mod(kk,2)  ==  0 .and. mod(nguard,2) .ne. 0)) then
               ! this point is on the parent's face
               ! and does not need to be interpolated
               unk_e_x1(ivar,i,j,kk,idest) = f_intx(i,j,k)
            Else

               If (kk < nguard + nzb/2) Then
                  kw = 1
               Else
                  kw = 2
               End If  ! End If (kk < nguard + nzb/2)

               kmin = iminh(kw,order) + k
               kmax = imaxh(kw,order) + k

               Do kpar = kmin,kmax
                  unk_e_x1(ivar,i,j,kk,idest) =                      & 
                       unk_e_x1(ivar,i,j,kk,idest) +                 & 
                       weight_half(kw,kpar-k,order)*f_intx(i,j,kpar)
               End Do  ! End Do kpar = kmin,kmax
               ! update parent index
               k = k + 1
            End If  ! End If ((mod(kk,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or. 
                    !         (mod(kk,2)  ==  0 .and. mod(nguard,2) .ne. 0))

         End Do  ! End Do kk = kfmin,kfmax
      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do j = jfmin,jfmax

      End If  ! End If (ndim == 3)


!-----Begin interpolation of unk_e_z

      
      Elseif (iedge_dir == 3) Then 


!-----Interpolate unk_e_z in x direction


      If (ndim >= 1) Then


      f_intx(:,:,:) = 0.

      kmin = kcmin-(nguard+2)*k3d
      If (kmin < 1) kmin = 1
      kmax = kcmax+(nguard+2)*k3d
      If (kmax > nguard*2+nzb) kmax = nguard*2+nzb
      jmin = jcmin-(nguard+2)*k2d
      If (jmin < 1) jmin = 1
      jmax = jcmax+(nguard+2)*k2d + k2d
      If (jmax > nguard*2+nyb+k2d) jmax = nguard*2+nyb+k2d

      Do k = kmin,kmax
      Do j = jmin,jmax

         ! 1) now interpolate to half points

         ! starting parent index
         i = icmin

         Do ii = ifmin,ifmax

            If ((mod(ii,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.   & 
                (mod(ii,2)  ==  0 .and. mod(nguard,2) .ne. 0)) Then
               ! this point is on parent's face
               ! and does not need to be interpolated
               f_intx(ii,j,k) = recv(ivar,i,j,k)
            Else
               If (ii < nguard + nxb/2) Then
                  iw = 1
               Else
                  iw = 2
               End If  ! End If (ii < nguard + nxb/2)
               
               imin = iminh(iw,order) + i
               imax = imaxh(iw,order) + i

               Do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +                    & 
                       weight_half(iw,ipar-i,order)*recv(ivar,ipar,j,k)
               End Do  ! End Do ipar = imin,imax
               ! update parent index
               i = i + 1
            End If  ! End If ((mod(ii,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.
                    !         (mod(ii,2)  ==  0 .and. mod(nguard,2) .ne. 0))

         End Do  ! End Do ii = ifmin,ifmax
      End Do  ! End Do j = jmin,jmax
      End Do  ! End Do k = kmin,kmax

      If (ndim == 1) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_e_z1(ivar,i,j,k,idest) = f_intx(i,j,k)
               End Do  ! End Do i = ifmin,ifmax
            End Do  ! End Do j = jfmin,jfmax
         End Do  ! End Do k = kfmin,kfmax
      End If  ! End If (ndim == 1)

      End If  ! End If (ndim >= 1)


!-----Interpolate unk_e_z in z direction


      If (ndim >= 2) Then
      

      f_intz(:,:,:) = 0.

      jmin = jcmin-(nguard+2)*k2d
      If (jmin < 1) jmin = 1
      jmax = jcmax+(nguard+2)*k2d + k2d
      If (jmax > nguard*2+nyb+k2d) jmax = nguard*2+nyb+k2d

      Do j = jmin,jmax
      Do i = ifmin,ifmax

         ! 1) interpolate to half points

         ! starting parent index
         k = kcmin
            
         Do kk = kfmin,kfmax
            
            If ((mod(kk,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.   & 
                (mod(kk,2)  ==  0 .and. mod(nguard,2)  ==  0)) Then

                                ! right point

               If (kk < nguard + nzb/2) Then
                  kw = 1
               Else
                  kw = 2
               End If  ! End If (kk < nguard + nzb/2)

               kmin = iminr(kw,order) + k
               kmax = imaxr(kw,order) + k

               Do kpar = kmin,kmax
                  f_intz(i,j,kk) = f_intz(i,j,kk) +                   & 
                       weight_right(kw,kpar-k,order)*f_intx(i,j,kpar)
               End Do  ! End Do kpar = kmin,kmax
               ! update parent index
               k = k + 1

            Else
                                ! left point

               If (kk < nguard + nzb/2) Then
                  kw = 1
               Else
                  kw = 2
               End If  ! End If (kk < nguard + nzb/2)

               kmin = iminl(kw,order) + k
               kmax = imaxl(kw,order) + k

               Do kpar = kmin,kmax
                  f_intz(i,j,kk) = f_intz(i,j,kk) +                  & 
                       weight_left(kw,kpar-k,order)*f_intx(i,j,kpar)
               End Do  ! End Do kpar = kmin,kmax

            End If  ! End If ((mod(kk,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
                    !         (mod(kk,2)  ==  0 .and. mod(nguard,2)  ==  0))
            
         End Do  ! End Do kk = kfmin,kfmax
      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do j = jmin,jmax

      If (ndim == 2) Then
         Do k = kfmin,kfmax
            Do j = jfmin,jfmax
               Do i = ifmin,ifmax
                  unk_e_z1(ivar,i,j,k,idest) = f_intz(i,j,k)
               End Do  ! End Do i = ifmin,ifmax
            End Do  ! End Do j = jfmin,jfmax
         End Do  ! End Do k = kfmin,kfmax
      End If  ! End If (ndim == 2)

      End If  ! End If (ndim >= 2)


      
!-----Interpolate unk_e_z in y direction



      If (ndim == 3) Then


      Do k = kfmin,kfmax
      Do i = ifmin,ifmax

         ! 1) interpolate to half points

         ! starting parent index
         j = jcmin
         
         Do jj = jfmin,jfmax

            unk_e_z1(ivar,i,jj,k,idest) = 0.

            If ((mod(jj,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.   & 
                (mod(jj,2)  ==  0 .and. mod(nguard,2) .ne. 0)) Then
               ! this point is on parent's face
               ! and does not need to be interpolated
               unk_e_z1(ivar,i,jj,k,idest) = f_intz(i,j,k) 
            Else
               If (jj < nguard + nyb/2) Then
                  jw = 1
               Else
                  jw = 2
               End If  ! End If (jj < nguard + nyb/2)

               jmin = iminh(jw,order) + j
               jmax = imaxh(jw,order) + j

               Do jpar = jmin,jmax
                  unk_e_z1(ivar,i,jj,k,idest) =                      & 
                       unk_e_z1(ivar,i,jj,k,idest) +                 & 
                       weight_half(jw,jpar-j,order)*f_intz(i,jpar,k)
               End Do  ! End Do jpar = jmin,jmax
               ! update parent index
               j = j + 1
            End If  ! End If ((mod(jj,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.
                    !         (mod(jj,2)  ==  0 .and. mod(nguard,2) .ne. 0))

         End Do  ! End Do jj = jfmin,jfmax
      End Do  ! End Do i = ifmin,ifmax
      End Do  ! End Do k = kfmin,kfmax

      End If  ! End If (ndim == 3)

      End If  ! End If (iedge_dir == 2)
      
      End Subroutine amr_1blk_ec_prol_genorder






