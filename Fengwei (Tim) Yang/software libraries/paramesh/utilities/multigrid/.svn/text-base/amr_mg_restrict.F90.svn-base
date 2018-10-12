!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_mg_restrict (nprocs, mype2, level)

      use tree
      use paramesh_interfaces, ONLY: amr_restrict
      
      implicit none
      
      integer, intent(in) :: nprocs, mype2, level
      integer             :: lb
      
!               Call the PARAMESH restriction routine.

      call amr_restrict (mype2, 2, 0, .false.)
      
      do lb = 1, lnblocks
         
         if (lrefine(lb) > level-1)  nodetype(lb) = -1
         if (lrefine(lb) == level-1) nodetype(lb) = 1
         if (lrefine(lb) == level-2 .and. nodetype(lb) /= 1)  & 
     &        nodetype(lb) = 2
      
      end do

      call amr_get_new_nodetypes (nprocs, mype2, level-1)

      return
      end subroutine amr_mg_restrict
