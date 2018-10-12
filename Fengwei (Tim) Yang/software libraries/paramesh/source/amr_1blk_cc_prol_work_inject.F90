!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_work_inject
!! NAME
!! 
!!   amr_1blk_cc_prol_work_inject
!!
!! SYNOPSIS
!!
!!      Call amr_1blk_cc_prol_work_inject(recvt,ia,ib,ja,jb,ka,kb,  
!!                                        idest,ioff,joff,koff,mype)
!!      Call amr_1blk_cc_prol_work_inject(real, 
!!                                        integer, integer, integer, integer,
!!                                        integer, integer, integer, integer,
!!                                        integer, integer, integer)
!!
!! ARGUMENTS
!!
!!  Real,    intent(inout) :: recvt(:,:,:,:)                                            
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
!! INCLUDES
!!
!!   No includes
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
!!   amr_abort
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit the data in recv is interpolated    
!!   and placed into the work1 array. 
!! 
!! DESCRIPTION
!!
!!   This routine takes data from the array recvtw1, originally extracted 
!!   from the workspace array work on some block, 
!!   and performs a prolongation operation on it, between the bounds ranges 
!!   ia to ib, ja to jb, and ka to kb. The data in recvtw1 is from a parent 
!!   block and the result of the prolongation operation is returned in
!!   the working block `work' array work1.
!!   The position of the child within the 
!!   parent block is specified by the ioff, joff and koff arguments.
!!
!!   This particular prolongation is simple injection.
!!
!!   Note: before using this routine in your program, make sure that the
!!   routine prolong_fun_init has been called.
!!
!! AUTHOURS
!!
!!   Written :     Peter MacNeice          January 1997
!!
!!***

      subroutine amr_1blk_cc_prol_work_inject(recvt,   & 
             ia,ib,ja,jb,ka,kb,                        & 
             idest,ioff,joff,koff,mype)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use prolong_arrays

      Implicit None

!-----Input/Output Variables
      Real,    Intent(inout) :: recvt(:,:,:)
      Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb
      Integer, Intent(in)    :: idest,ioff,joff,koff,mype

!-----Local variables
      Real    :: dx,dy,dz,cx,cy,cz
      Integer :: icl,icu,jcl,jcu,kcl,kcu
      Integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p
      Integer :: offi,offj,offk

      Integer, Parameter :: largei = 100

!-----Begin Executable Code

      If (prolw_init.ne.100) Then
       Write(*,*) 'PARAMESH ERROR !'
       Write(*,*) 'Error : prolong_work_fun. ',      & 
             'You must call amr_prolong_fun_init ',  & 
             'before you can use this routine!'
       Call amr_abort
      End If  ! End If (prolw_init.ne.100)


!-------Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb


        offi = 0
        offj = 0
        offk = 0
        If (ioff > 0) offi = nxb/2
        If (joff > 0) offj = nyb*k2d/2
        If (koff > 0) offk = nzb*k3d/2

!-------Interpolation loop.
        Do k=kcl,kcu
             k1 = ((k-nguard_work-1+largei)/2 +                      & 
                      nguard_work - largei/2 )*k3d + 1 +offk
             k1p= k1
             dz = 1.
             cz = 0.
             Do j=jcl,jcu
                   j1 = ((j-nguard_work-1+largei)/2 +                & 
                            nguard_work - largei/2 )*k2d + 1 + offj
                   j1p= j1
                   dy = 1.
                   cy = 0.
                   Do i=icl,icu
                         i1 = (i-nguard_work-1+largei)/2 +           & 
                                 nguard_work - largei/2 + 1 + offi
                         i1p = i1
                         dx = 1.
                         cx = 0.

!------------------------compute interpolated values at location (i,j,k)
                         work1(i,j,k,idest) =                        & 
                               dz*( dy*( dx*recvt(i1,j1,k1) +        & 
                               cx*recvt(i1p,j1,k1))  +               & 
                               cy*( dx*recvt(i1,j1p,k1) +            & 
                               cx*recvt(i1p,j1p,k1) ) ) +            & 
                               cz*( dy*( dx*recvt(i1,j1,k1p) +       & 
                               cx*recvt(i1p,j1,k1p))  +              & 
                               cy*( dx*recvt(i1,j1p,k1p) +           & 
                               cx*recvt(i1p,j1p,k1p) ) )



                  End Do  ! End Do i=icl,icu
             End Do  ! End Do j=jcl,jcu
        End Do  ! End Do k=kcl,jcl


      Return
      End Subroutine amr_1blk_cc_prol_work_inject
