!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_mg_morton_process()

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_mpi_interfaces, only : mpi_morton_bnd, & 
     &                                    mpi_morton_bnd_prolong, & 
     &                                    mpi_morton_bnd_fluxcon, & 
     &                                    mpi_morton_bnd_restrict, & 
     &                                    mpi_setup

      implicit none

      include 'mpif.h'

      integer :: level, max_level, itemp, ierr, lb
      integer :: tag_offset
      logical :: lec, lnc, lfulltree
      logical, save :: first = .true.
      integer, save :: nprocs, mype
      integer :: nodetype_old(maxblocks_tr)
      logical :: newchild_old(maxblocks_tr)

      if (grid_changed == 1) then

      if (first) then
         call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
         call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)
         first = .false.
      end if

      nodetype_old(1:lnblocks) = nodetype(1:lnblocks)
      newchild_old(1:lnblocks) = newchild(1:lnblocks)

      max_level = 0
      do lb = 1,lnblocks
         max_level = max(max_level,lrefine(lb))
      end do
      call MPI_ALLREDUCE(max_level,itemp,1,MPI_INTEGER, & 
     &                   MPI_MAX,MPI_COMM_WORLD,ierr)
      max_level = itemp

      do level = max_level, 1, -1

         newchild(:) = .FALSE.
         do lb = 1,lnblocks

!         if (nodetype(lb) == 1 .and. lrefine(lb) == level)
!     &        nodetype(lb) = -1
!         if (nodetype(lb) == 2 .and. lrefine(lb) == level-1) then
!            nodetype(lb) = 1
!         end if
!         if (nodetype(lb) == 3 .and. lrefine(lb) == level-2)
!     &        nodetype(lb) = 2

            nodetype(lb) = nodetype_old(lb)
            if (lrefine(lb) > level) & 
     &           nodetype(lb) = -1
            if (lrefine(lb) == level) then
               nodetype(lb) = 1
               newchild(lb) = .TRUE.
            end if
            if (lrefine(lb) == level-1 .and. nodetype(lb) /= 1) & 
     &           nodetype(lb) = 2

         end do

         if (first) then
            call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
            call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)
            first = .false.
         end if
!---------------------------
! call setup routines in preparation for calling
! all the mpi_morton_bnd_XXX routines.

! Find the coordinate ranges
         call mpi_amr_global_domain_limits
         call mpi_setup(mype,nprocs)

!---------------------------

         tag_offset = 100
         call mpi_morton_bnd(mype,nprocs,tag_offset)
         call mpi_amr_write_guard_comm_mg(nprocs,level)

         tag_offset = 100
         call mpi_morton_bnd_prolong(mype,nprocs,tag_offset)
         call mpi_amr_write_prol_comm_mg(nprocs,level)

         tag_offset = 100
         call mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)
         call mpi_amr_write_flux_comm_mg(nprocs,level)

         lec = .false.
         lnc = .false.
         if(nvaredge.gt.0) lec = .true.
         if(nvarcorn.gt.0) lnc = .true.
         tag_offset = 100
         lfulltree = .false.
         call mpi_morton_bnd_restrict(mype,nprocs, & 
     &                                tag_offset)
         call mpi_amr_write_restrict_comm_mg(nprocs,level)

       end do

       nodetype(1:lnblocks) = nodetype_old(1:lnblocks)
       newchild(1:lnblocks) = newchild_old(1:lnblocks)

       call amr_morton_process()

       grid_changed = 0

       end if

       end subroutine amr_mg_morton_process
