!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
subroutine amr_1blk_bcset(mype,ibc,lb,pe, & 
     &    idest,iopt,iface,jface,kface,surrblks)

! Arguments:
!      mype             local processor
!      ibc              the integer specifying the particular boundary
!                        condition to be imposed
!      lb               block number of selected block
!      pe               processor on which block lb is located
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      iface            a selector setting designating whether the guardcells
!                        to be set are at the left, center or right sections
!                        of the i index range, eg
!                           iface = -1      left end
!                                 =  0      middle
!                                 = +1      right. For example, if iface=-1,
!                        the i index applied when filling unk will run
!                        from 1:nguard, if iface=0 from 1+nguard:nxb+nguard,
!                        and if iface=+1 from nxb+nguard+1:nxb+2*nguard.
!      jface            a selector setting designating whether the guardcells
!                        to be set are at the left, center or right sections
!                        of the j index range.
!      kface            a selector setting designating whether the guardcells
!                        to be set are at the left, center or right sections
!                        of the k index range.
!
!
! Written :     Peter MacNeice          August 1998
! Modified:     Peter MacNeice          January 2001
! Modified:     J. Green                April 2008
!
! CEG notes: have now put in correct Neumann boundary conditions for the 
!            12 surrounding edges.  The corners have not been fully tested as 
!            are currently unused until it goes to a 27-point stencil.
!------------------------------------------------------------------------
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  use solution_parameters
  implicit none
  include 'mpif.h'
  
  integer, intent(in) :: mype,ibc,lb,pe
  integer, intent(in) :: idest,iopt,iface,jface,kface
  integer, intent(in) :: surrblks(:,:,:,:)
  integer :: i,j,k,nguard0,offset
  integer :: i1=0, i2=0, j1=0, j2=0, k1=0, k2=0, switch
  double precision dx, dy, dz
  
  
  !---------------------------------------------------------------------------
  ! Section to be modified by user
  
  ! Which boundary condition has been specified? ibc is the value
  ! of NEIGH for the current block face. If ibc is less than or equal
  ! to -20 you are on an external boundary. This if test should 
  ! conditionally execute the appropriate code for each of the different 
  ! boundary conditions you wish to impose.
  ! The example here serves to indicate the index ranges which need to be set.
  nguard0 = nguard*npgs
  !if(jface.eq.-1.and.kface.eq.-1)then
  !   !axis_node_count=axis_node_count+1
  !end if
  !if(bnd_box(1,2,lb).eq.nucleate_y)then
  !   ! lower y bnd matches
  !   if(bnd_box(1,3,lb).eq.nucleate_z)then
  !      ! lower z bnd matches
  !      ! y and z co-ords match so we're on the x axis
  !      axis_node_count=axis_node_count+1
  !   end if
  !end if
  !-------------------------

! CEG  :  Case labelling is as follows
!
!	kface=-1		kface=0			kface=1
!	
!	1e1  |  2f  |  3e1	1e2  |  2g  |  3e2	1e3  |  2h  |  3e3	
!	------------------	------------------	------------------
!	 1b  |  2d  |  3b	 1c  |  -   |  3c	 1d  |  2e  |  3d
!	------------------	------------------	------------------
!	1a1  |  2a  |  3a1	1a2  |  2b  |  3a2	1a3  |  2c  |  3a3

  if(ibc.eq.-21 ) then

     switch = 9*kface +3*jface+iface
!	 -7 |  -6   |  -5	  2  |  3  |   4	 11 |  12  |  13	
!	------------------	------------------	------------------
!	 -10 |  -9  |  -8	 -1  |  0   |   1	 8  |   9  |  10
!	------------------	------------------	------------------
!	 -13 | -12  |  -11	 -4  |  -3  |  -2	 5  |   6  |  7

     !-------------------------
     ! Boundary condition 2
     !
     !--------------------
     !--------------------        
     ! Do cell centered data
     if(nvar.gt.0) then
        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        dz = bsize(3,lb)/real(nzb)
        if(iface.eq.-1)then
           i1=1
           i2=nguard
        elseif(iface.eq.0)then
           i1 = nguard+1
           i2 = nxb+nguard
        elseif(iface.eq.1)then
           i1 = nxb+nguard+1
           i2 = nxb+2*nguard
        else
           print *,"Serious error in boundary allocations iface != -1,0,1"
        end if
        if(jface.eq.-1)then
           j1=1
           j2=nguard
        elseif(jface.eq.0)then
           j1 = nguard+1
           j2 = nyb+nguard
        elseif(jface.eq.1)then
           j1 = nyb+nguard+1
           j2 = nyb+2*nguard
        else
           print *,"Serious error in boundary allocations jface != -1,0,1"
        end if
        if(kface.eq.-1)then
           k1=1
           k2=nguard
        elseif(kface.eq.0)then
           k1 = k3d*nguard+1
           k2 = nzb+k3d*nguard
        elseif(kface.eq.1)then
           k1 = nxb+nguard+1
           k2 = nzb+2*nguard
        else
           print *,"Serious error in boundary allocations kface != -1,0,1"
        end if
        do k = k1,k2
           do j = j1,j2
              do i = i1,i2
                 if(iopt.eq.1)then
                    !unk1(,i,j,k,idest) = 0.0

                    select case (switch)
                       case(-13)
! 1 a1 - extras added by CEG 13/4/10 - was just ,k, in James version
                             if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+1-j,k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+1-j,2*nguard+1-k,idest) 
                             endif

                       case(-4)
!1 a2
                             if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+1-j,k,idest) 
                             endif

                       case(5)
!1 a3
                             if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20 .and. surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif
                       case(-10)
! 1 b
                             if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,2*nguard+1-k,idest)
                             endif
                       case(-1)
! 1 c
                             unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest)

                       case(8)
! 1 d
                             if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,2*nguard+2*nzb+1-k,idest)
                             endif
                       case(-7)
!1 e1
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,1)>-20 .and. surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             endif

                       case(2)
!1 e2
                             if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest)
                             else if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+2*nyb+1-j,k,idest)
                             endif
                       case(11)
!1 e3
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif

                       case(-12)
! 2 a
                             if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,2*nguard+1-k,idest)
                             endif
                       case(-3)
! 2 b
                             unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest)
! 2 c
                       case(6)
                             if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest)
                             endif

                        case(-9)
! 2 d
                             unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest)
                        case(9)
! 2 e
                             unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest)
                        case(-6)
! 2 f
                             if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest)
                             endif
                        case(3)
! 2 g
                             unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest)
                        case(12)
! 2 h
                             if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest)
                             endif

                        case(-11)
! 3 a1
                             if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,1)>-20 .and. surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20 .and. surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+1-j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,2*nguard+1-k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+1-j,2*nguard+1-k,idest) 
                             endif
                        case(-2)
! 3 a2
                             if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+1-j,k,idest) 
                             endif
                        case(7)
! 3 a3
                             if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20 .and. surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif

                        case(-8)
! 3 b
                             if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,2*nguard+1-k,idest)
!                                print *,'what 3b case are we doing?', ibc,surrblks(1,:,:,:)
                             endif
!          if (lb.eq.4.or.lb.eq.4100) print *,'3b blk',lb, unk1(:,i,j,k,idest)
                        case(1)
! 3 c
                             unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest)
                        case(10)
! 3 d
                             if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,2*nguard+2*nzb+1-k,idest)
                             endif
! 3 e - extras added by CEG 13/4/10 - was just ,k, in James version
                        case(-5)
! 3 e1
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,1)>-20 .and. surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,2,3,2)>-20 .and. surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             endif

                        case(4)
! 3 e2
                             if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest)
                             else if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest)
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,k,idest)
                             endif

                        case(13)
! 3 e3
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,2,3,2)>-20 .and. surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,2,3)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif

                      end select



                    !if(iface.eq.-1)then
                    !   unk1(:,i,j,k,idest) = unk1(:,2*nguard+1-i,j,k,idest)
                    !else if(iface.eq.1) then
                    !   unk1(:,i,j,k,idest) = unk1(:,2*nguard+2*nxb+1-i,j,k,idest)
                    !else if(jface.eq.-1)then
                    !   unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+1-j,k,idest)
                    !else if(jface.eq.1) then
                    !   unk1(:,i,j,k,idest) = unk1(:,i,2*nguard+2*nyb+1-j,k,idest)
                    !else if(kface.eq.-1)then
                    !   unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+1-k,idest)
                    !else if(kface.eq.1) then
                    !   unk1(:,i,j,k,idest) = unk1(:,i,j,2*nguard+2*nzb+1-k,idest)
                    !end if
                 else if((iopt.eq.2).or.(iopt.eq.3).or.(iopt.eq.4).or.(iopt.eq.5))then
!                    print *,'iopt is ',iopt
!CEG definitely used - iopt=3 occurs pretty soon
                    select case (switch)
                       case(-13)
! 1 a - extras added by CEG 13/4/10 - was just ,k, in James version
                             if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+1-j,k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+1-j,2*nguard+1-k,idest) 
                             endif
                       case(-4)
!1 a2
                             if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+1-j,k,idest) 
                             endif
                       case(5)
!1 a3
                             if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20 .and. surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif
                       case(-10)
! 1 b
                             if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,2*nguard+1-k,idest)
                             endif
                       case(-1)
                             work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest)
                       case(8)
! 1 d
                             if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,2*nguard+2*nzb+1-k,idest)
                             endif
                       case(-7)
!1 e1
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,1)>-20 .and. surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             endif
                       case(2)
!1 e2
                             if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest)
                             else if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+2*nyb+1-j,k,idest)
                             endif
                       case(11)
!1 e3
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20 .and. surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,1,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+1-i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif
                       case(-12)
! 2 a
                             if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,2*nguard+1-k,idest)
                             endif
                       case(-3)
                             work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest)
                       case(6)
! 2 c
                             if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest)
                             endif
                       case(-9)
                             work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest)
                       case(9)
                             work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest)
                       case(-6)
! 2 f
                             if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest)
                             endif
                       case(3)
                             work1(i,j,k,idest) = work1(i,j,k,idest)
                       case(12)
! 2 h
                             if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest)
                             endif
                       case(-11)
! 3 a1
                             if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,1)>-20 .and. surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20 .and. surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+1-j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,2*nguard+1-k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+1-j,2*nguard+1-k,idest) 
                             endif
                       case(-2)
! 3 a2
                             if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+1-j,k,idest) 
                             endif
                       case(7)
! 3 a3
                             if (surrblks(1,2,1,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20 .and. surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,1,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif
                       case(-8)
! 3 b
                             if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest)
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,2*nguard+1-k,idest)
                             endif
                       case(1)
                             work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest)
                       case(10)
! 3 d
                             if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest)
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,2*nguard+2*nzb+1-k,idest)
                             endif
!                             work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,2*nguard+2*nzb+1-k,idest)
                       case(-5)
! 3 e1
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,1)>-20 .and. surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,2,3,2)>-20 .and. surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,2*nguard+1-k,idest) 
                             else if (surrblks(1,2,2,1)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,2*nguard+1-k,idest) 
                             endif
                       case(4)
! 3 e2
                             if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest)
                             else if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest)
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,k,idest)
                             endif
                       case(13)
! 3 e3
                             if (surrblks(1,2,3,2)>-20 .and. surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,k,idest) 
                             else if (surrblks(1,2,2,3)>-20 .and. surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,2,3,2)>-20 .and. surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,3,2)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,j,2*nguard+2*nzb+1-k,idest) 
                             else if (surrblks(1,2,2,3)>-20) then
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,k,idest) 
                             else if (surrblks(1,3,2,2)>-20) then
                                work1(i,j,k,idest) = work1(i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             else
                                work1(i,j,k,idest) = work1(2*nguard+2*nxb+1-i,2*nguard+2*nyb+1-j,2*nguard+2*nzb+1-k,idest) 
                             endif
                      end select

                 else
                    print *,"Unhandled boundary condition call with iopt=",iopt,"ibc=",ibc
                 end if
              enddo
           enddo
        enddo
        
        
           
     endif                          ! end of nvar if test
     

     !-------------------------
  elseif(ibc.eq.-22.or.ibc.eq.-23 ) then
     print *,"Unhandled boundary condition call with iopt=",iopt,"ibc=",ibc
     print *,"removed at TriPhaseField revision 28"
     !-------------------------
  endif                            ! end of test of bc flag
  !-------------------------
  
  
  
  ! End of Section to be modified by user
  !---------------------------------------------------------------------------
  
  return
end subroutine amr_1blk_bcset
