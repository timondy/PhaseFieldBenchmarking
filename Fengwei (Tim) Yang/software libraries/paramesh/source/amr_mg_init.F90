!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_mg_init()

      use tree
      use amr_mg_common

      logical, save :: first = .true.

      if (first) then
         allocate(nodetype_old(maxblocks_tr))
         first = .false.
      end if

      nodetype_old(1:lnblocks) = nodetype(1:lnblocks)

      call amr_mg_morton_process()

      return
      end subroutine amr_mg_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

      subroutine pf_mg_init()

      use tree
      use amr_mg_common

      logical, save :: first = .true.

      if (first) then
         allocate(nodetype_old(maxblocks_tr))
         first = .false.
      end if

      nodetype_old(1:lnblocks) = nodetype(1:lnblocks)

!      call amr_mg_morton_process()
      call pf_mg_morton_process()

      return
      end subroutine pf_mg_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

! Added by CEG to clean up memory at end

      subroutine amr_mg_end()

      use tree
      use amr_mg_common

      deallocate(nodetype_old)

      call amr_mg_morton_process_end()

      return
      end subroutine amr_mg_end
