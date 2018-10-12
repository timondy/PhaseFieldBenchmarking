!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_prol_work_linear
!! NAME
!!
!!   amr_1blk_cc_prol_work_linear
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_work_linear(recvt,ia,ib,ja,jb,ka,kb,idest,                   
!!                                     ioff,joff,koff,mype)                            
!!   Call amr_1blk_cc_prol_inject (real,                                                
!!                                 integer, integer, integer, integer,                  
!!                                 integer, integer, integer, integer,                  
!!                                 integer, integer, integer)                  
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
!!   and placed into the unk1 array.                    
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
!!   This particular prolongation is simple linear interpolation. It can
!!   be used for blocks with an even or odd number of grid cells.
!!
!!   Conservative prolongation. Special treatment for the  cells immediately
!!   adjacent to a boundary (ie i=nguard_work,nguard_work+1,
!!   iu_bnd-nguard_work,iu_bnd-nguard_work+1 and likewise for j and k indeces)
!!   if using an even number of grid cells  per block along that axis. 
!!   No special treatment is required when the number
!!   of cells is odd.
!!
!!   Note: before using this routine in your program, make sure that the
!!   routine prolong_fun_init has been called.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 1997
!!
!!***

      Subroutine amr_1blk_cc_prol_work_linear(recvt,  & 
             ia,ib,ja,jb,ka,kb,                       & 
             idest,ioff,joff,koff,mype)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use prolong_arrays

      implicit none

!-----Input/Output Variables
      real,    intent(inout) :: recvt(:,:,:)
      integer, intent(in)    :: ia,ib,ja,jb,ka,kb
      integer, intent(in)    :: idest,ioff,joff,koff,mype

!-----Local Variables
      real    :: dx,dy,dz,cx,cy,cz
      integer :: icl,icu,jcl,jcu,kcl,kcu,i_ind,j_ind,k_ind
      integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p

!-----Begin Exectuable Code

      If (prolw_init.ne.100) Then
       Write(*,*) 'PARAMESH ERROR !'
       Write(*,*) 'Error : prolong_work_fun. ', & 
             'You must call amr_prolong_fun_init ', & 
             'before you can use this routine!'
       Call amr_abort
      End If  ! End If (prolw_init.ne.100)


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
      if(ioff.gt.0) i_ind = 2
      if(joff.gt.0) j_ind = 2
      if(koff.gt.0) k_ind = 2

!-------Interpolation loop.
        Do k=kcl,kcu
             k1 = prolw_indexz(1,k,k_ind)
             k1p= prolw_indexz(2,k,k_ind)
             dz = prolw_dz(k)
             cz = 1.-dz
             Do j=jcl,jcu
                   j1 = prolw_indexy(1,j,j_ind)
                   j1p= prolw_indexy(2,j,j_ind)
                   dy = prolw_dy(j)
                   cy = 1.-dy
                   Do i=icl,icu
                         i1 = prolw_indexx(1,i,i_ind)
                         i1p= prolw_indexx(2,i,i_ind)
                         dx = prolw_dx(i)
                         cx = 1.-dx

!------------------------compute interpolated values at location (i,j,k)
                         work1(i,j,k,idest) =                    & 
                                dz*( dy*( dx*recvt(i1,j1,k1) +   & 
                                cx*recvt(i1p,j1,k1))  +          & 
                                cy*( dx*recvt(i1,j1p,k1) +       & 
                                cx*recvt(i1p,j1p,k1) ) ) +       & 
                                cz*( dy*( dx*recvt(i1,j1,k1p) +  & 
                                cx*recvt(i1p,j1,k1p))  +         & 
                                cy*( dx*recvt(i1,j1p,k1p) +      & 
                                cx*recvt(i1p,j1p,k1p) ) )


                  End Do  ! End Do i=icl,icu
             End Do  ! End Do j=jcl,jcu
        End Do  ! End Do k=kcl,kcu


      Return
      End Subroutine amr_1blk_cc_prol_work_linear
