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


      subroutine mesh_test(mype)





!------------------------------------------------------------------------
!
! This routine tests the mesh to see if there is a jump of more
! than one refinement level anywhere in the grid. It does it
! by looping over all parents of leaf nodes, getting the list of blocks
! surrounding them, and checking that the appropriate
! parts of this list exist.
!
! Arguments:
!      mype             local processor
!
!
! Written :     Peter MacNeice          August 1998
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use paramesh_mpi_interfaces, only : mpi_morton_bnd

      implicit none

      include 'mpif.h'

#ifdef TIMINGS
#include "timer.fh"
#endif

      integer, intent(in) :: mype

!-------------------------


! local arrays

      real :: xtest,ytest,ztest

      logical ldiag
      integer psurround(3,3,3,3)
      integer psurrblks(3,3,3,3)
      integer :: tag_offset
      integer :: nprocs,lb,i,j,k,ierror,max_no_of_blocks

!-------------------------

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)

      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)
      call MPI_ALLREDUCE(lnblocks, & 
     &                   max_no_of_blocks, & 
     &                   1, & 
     &                   MPI_INTEGER, & 
     &                   MPI_MAX, & 
     &                   MPI_COMM_WORLD, & 
     &                   ierror)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Loop over blocks.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.2) then


        ldiag = .true.
        psurround(:,:,2-k2d:2+k2d,2-k3d:2+k3d) = & 
     &       surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb)
 
       if (ndim == 1) then
       do i = 1,3,2
        if(psurround(1,i,2,2).gt.-20.and.psurround(1,i,2,2).lt.0) then
         write(*,*) 'Bad grid : 2 level jump near block  (', & 
     &        mype,lb,') ','  surrounding block (',i,',2,2)'
        endif
       enddo
       end if

       if (ndim == 2) then
       do j = 2-k2d,2+k2d
       do i = 1,3
        if(i.ne.2.or.j.ne.2) then
        if(psurround(1,i,j,2).gt.-20.and.psurround(1,i,j,2).lt.0) then
         write(*,*) 'Bad grid : 2 level jump near block  (', & 
     &        mype,lb,') ','  surrounding block (',i,j,',2)'
        endif
        endif
       enddo
       enddo
       end if

       if (ndim == 3) then
       do k = 2-k3d,2+k3d
       do j = 2-k2d,2+k2d
       do i = 1,3
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then
        if(psurround(1,i,j,k).gt.-20.and.psurround(1,i,j,k).lt.0) then
         write(*,*) 'Bad grid : 2 level jump near block  (', & 
     &        mype,lb,') ','  surrounding block (',i,j,k,')'
        endif
        endif
       enddo
       enddo
       enddo 
       end if


      endif
      enddo
      endif


      Call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! Check consistency of bounding box and coordinate info.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
        xtest = (coord(1,lb)-bnd_box(1,1,lb))* & 
     &          (coord(1,lb)-bnd_box(2,1,lb))
        if(xtest.ge.0.) write(*,*) 'coord and bnd_box are ', & 
     &'inconsistent in x direction for block ',lb,' proc ',mype
        if (ndim >= 2) then
        ytest = (coord(2,lb)-bnd_box(1,2,lb))* & 
     &          (coord(2,lb)-bnd_box(2,2,lb))
        if(ytest.ge.0.) write(*,*) 'coord and bnd_box are ', & 
     &'inconsistent in y direction for block ',lb,' proc ',mype
        end if
        if (ndim == 3) then
        ztest = (coord(3,lb)-bnd_box(1,3,lb))* & 
     &          (coord(3,lb)-bnd_box(2,3,lb))
        if(ztest.ge.0.) write(*,*) 'coord and bnd_box are ', & 
     &'inconsistent in z direction for block ',lb,' proc ',mype
        end if
      enddo
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD,ierror)

      return
      end subroutine mesh_test
