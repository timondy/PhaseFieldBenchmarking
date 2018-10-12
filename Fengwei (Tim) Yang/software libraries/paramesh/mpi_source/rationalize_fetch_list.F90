!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

       subroutine rationalize_fetch_list (fetch_list,                  &
                                          istart,                      &
                                          iend,                        &
                                          nptsneigh)

! This routine analyses the parts of a block which are required,
! ensuring that a large enough part of the block is transmitted
! to satisfy the guardcell filling operation.


       use paramesh_dimensions
       use physicaldata
       use tree
       use mpi_morton

       implicit none

       integer, intent(in) :: istart, iend, nptsneigh
       integer, intent(inout) :: fetch_list(3,nptsneigh)

       integer :: fetch_list_old(3,istart:iend)
       integer :: imask(27)
       integer :: imask_loc(1)
       integer :: i,ib
       logical :: lfullblock

! initialize mask array
       imask = 0

       if(istart.eq.iend) return

       do i = istart,iend
          fetch_list_old(3,i) = fetch_list(3,i)
          if (fetch_list(3,i) > 27) fetch_list(3,i) =  & 
     &                                 fetch_list(3,i) - 27
       end do
       
! mark mask array
       lfullblock = .false.
       do ib = istart,iend
         imask(fetch_list(3,ib)) = 1
         if(fetch_list(3,ib).eq.14) lfullblock = .true.
       enddo
       if(lfullblock) then
         fetch_list(3,istart:iend)=14
         go to 2
       endif

! combine corner pairs to edge centers.
       if(ndim.ge.1) then 
       imask(14) = max(imask(14),(imask(13)+imask(15))/2)

! if full block is selected through combination of sides, 
! set sides mask to zero.
       imask(13) = max(0,imask(13)-imask(14))
       imask(15) = max(0,imask(15)-imask(14))
       endif

       if(ndim.ge.2) then 
       imask(11) = max(imask(11),(imask(10)+imask(12))/2)
       imask(13) = max(imask(13),(imask(10)+imask(16))/2)
       imask(15) = max(imask(15),(imask(12)+imask(18))/2)
       imask(17) = max(imask(17),(imask(16)+imask(18))/2)

! if any faces are formed zero out the edges responsible so
! they can be ignored below
       imask(10) = max(0,imask(10)-max(imask(11),imask(13)))
       imask(12) = max(0,imask(12)-max(imask(11),imask(15)))
       imask(16) = max(0,imask(16)-max(imask(13),imask(17)))
       imask(18) = max(0,imask(18)-max(imask(15),imask(17)))

       endif


! combine corner pairs to edge centers.
       if(ndim.eq.3) then

       imask(2)  = max(imask(2 ),(imask(1)+imask(3))/2)
       imask(4)  = max(imask(4 ),(imask(1)+imask(7))/2)
       imask(6)  = max(imask(6 ),(imask(3)+imask(9))/2)
       imask(8)  = max(imask(8 ),(imask(7)+imask(9))/2)
       imask(10)  = max(imask(10 ),(imask(1)+imask(19))/2)
       imask(12)  = max(imask(12 ),(imask(3)+imask(21))/2)
       imask(16)  = max(imask(16 ),(imask(7)+imask(25))/2)
       imask(18)  = max(imask(18 ),(imask(9)+imask(27))/2)
       imask(20) = max(imask(20),(imask(19)+imask(21))/2)
       imask(22) = max(imask(22),(imask(19)+imask(25))/2)
       imask(24) = max(imask(24),(imask(21)+imask(27))/2)
       imask(26) = max(imask(26),(imask(25)+imask(27))/2)

! if any edges are formed zero out the corners responsible so
! they can be ignored below
       imask(1) = max(0,imask(1) - max(0,imask(2),imask(4),imask(10)))
       imask(3) = max(0,imask(3) - max(0,imask(2),imask(6),imask(12)))
       imask(7) = max(0,imask(7) - max(0,imask(4),imask(8),imask(16)))
       imask(9) = max(0,imask(9) - max(0,imask(6),imask(8),imask(18)))
       imask(19) = max(0,imask(19) -  & 
     &                           max(0,imask(10),imask(20),imask(22)))
       imask(21) = max(0,imask(21) - & 
     &                           max(0,imask(12),imask(20),imask(24)))
       imask(25) = max(0,imask(25) -  & 
     &                           max(0,imask(16),imask(22),imask(26)))
       imask(27) = max(0,imask(27) - & 
     &                           max(0,imask(18),imask(24),imask(26)))

       endif

! combine edges to make full faces
       if(ndim.eq.3) then
       imask(5)  = max( imask(5 ), & 
     &              max( 0, & 
     &        min(1,(imask(2)+imask(4)+imask(6)+imask(8))/2) ) )
       imask(23)  = max( imask(23), & 
     &              max(0,min(1, & 
     &              (imask(20)+imask(22)+imask(24)+imask(26))/2) ) )
       imask(13)  = max( imask(13), & 
     &              max(0,min(1, & 
     &              (imask(4)+imask(10)+imask(16)+imask(22))/2) ) )
       imask(15)  = max( imask(15), & 
     &              max(0,min(1, & 
     &              (imask(6)+imask(12)+imask(18)+imask(24))/2) ) )
       imask(11)  = max( imask(11), & 
     &              max(0,min(1, & 
     &              (imask(2)+imask(10)+imask(12)+imask(20))/2) ) )
       imask(17)  = max( imask(17), & 
     &              max(0,min(1, & 
     &             (imask(8)+imask(16)+imask(18)+imask(26))/2) ) )


       imask(2) = max(0,imask(2) - max(0,imask(5),imask(11)))
       imask(4) = max(0,imask(4) - max(0,imask(5),imask(13)))
       imask(6) = max(0,imask(6) - max(0,imask(5),imask(15)))
       imask(8) = max(0,imask(8) - max(0,imask(5),imask(17)))
       imask(20) = max(0,imask(20) - max(0,imask(23),imask(11)))
       imask(22) = max(0,imask(22) - max(0,imask(23),imask(13)))
       imask(24) = max(0,imask(24) - max(0,imask(23),imask(15)))
       imask(26) = max(0,imask(26) - max(0,imask(23),imask(17)))
       imask(10) = max(0,imask(10) - max(0,imask(13),imask(11)))
       imask(12) = max(0,imask(12) - max(0,imask(15),imask(11)))
       imask(16) = max(0,imask(16) - max(0,imask(13),imask(17)))
       imask(18) = max(0,imask(18) - max(0,imask(15),imask(17)))

! combine faces to make full block
       imask(14)  = max( imask(14 ), & 
     &            max(0,min(1, (imask(5)+imask(23)+imask(13)+imask(15) & 
     &                     +imask(11)+imask(17))/2 ) ) )

       endif

! If 1D calculation
       if(ndim.eq.1) then 
         if(imask(14).eq.1) fetch_list(3,istart:iend)=14
         go to 2
       endif

       if(ndim.ge.2) then 

       if(sum(imask).gt.1) then
         imask_loc(1)=14
       else
         imask_loc = maxloc(imask)
       endif
       fetch_list(3,istart:iend)=imask_loc(1)

       endif

 2     continue
       
       if (any(fetch_list_old(3,istart:iend) > 27)) then
          fetch_list(3,istart:iend) = fetch_list(3,istart:iend) + & 
     &                                   27
       end if

       return
       end subroutine rationalize_fetch_list

