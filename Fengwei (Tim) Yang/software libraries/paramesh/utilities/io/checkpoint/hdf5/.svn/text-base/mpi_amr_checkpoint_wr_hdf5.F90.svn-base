!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* utilities/io/checkpoint/hdf5/amr_checkpoint_wr_hdf5
!! NAME
!!
!!   amr_checkpoint_wr_hdf5
!!
!! SYNOPSIS
!!
!!   call amr_checkpoint_wr_hdf5(file_num)
!!   call amr_checkpoint_wr_hdf5(file_num, l_with_guardcells,  
!!                               user_attr1, user_attr2, user_attr3, 
!!                               user_attr4, user_attr5)
!!
!!   call amr_checkpoint_wr_hdf5(integer, optional logical, optional char*80,
!!                               optional real, optional real, optional real, 
!!                               optional real, optional real)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: file_num
!!     An integer number which will be appended to the end of the file name.
!!
!!   optional, logical, intent(in) :: l_with_guardcells
!!     If true, then guardcells are included in the checkpoint file.  Otherwise 
!!     (the default) they are not included.
!!   
!!   optional, real, intent(in) :: user_attr1(2,3,4,5)
!!     Arguments which allow the user to add some extra information to the file.  
!!     Currently only 5 real numbers can be added.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   io
!!
!! CALLS
!!
!!   write_blocks_hdf5_r8 or write_blocks_hdf5_r4
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a checkpoint file has been written.
!!
!! DESCRIPTION
!! 
!!  Subroutine to checkpoint PARAMESH runs in parallel using the HDF5 library.
!!  It writes out the tree data structure and data stored in PARAMESH blocks.
!!
!!  The files produced will have names of the form 
!!  'paramesh_chk_######.hdf5'. where '######' is the file_num argument passed into
!!  this routine.
!!
!!  The routine gets called directly by the driver routine 'amr_checkpoint_wr' if
!!  the user selects HDF5 format as the means of outputting the file.
!!
!!  The routine is not, by default, included in the PARAMESH default setup.  In order to
!!  use it an install script in utilities/io/checkpoint/hdf5 must first be run and
!!  PARAMESH must be recompiled.
!!
!! AUTHORS
!!
!!   Kevin Olson (2004)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f, ordering
#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_checkpoint_wr_hdf5 (file_num,  & 
     &                                   l_with_guardcells, & 
     &                                   user_attr_1, & 
     &                                   user_attr_2, & 
     &                                   user_attr_3, & 
     &                                   user_attr_4, & 
     &                                   user_attr_5)


      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use io

!-----------------------------!
! Start variable declarations !
!-----------------------------!

      implicit none

      include 'mpif.h'

      integer, intent(in)           :: file_num
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(in)    :: user_attr_1, & 
     &                                 user_attr_2, & 
     &                                 user_attr_3, & 
     &                                 user_attr_4, & 
     &                                 user_attr_5


      integer,dimension (:), allocatable :: n_to_left
      integer,dimension (:), allocatable :: glnblocks
      integer :: gid(nfaces+1+nchild,maxblocks_tr)
      integer :: nvar_chk_cc,nvar_chk_fc,nvar_chk_ec,nvar_chk_nc
      integer :: tot_blocks
      integer :: mype, nprocs
      integer :: il0, iu0, jl0, ju0, kl0, ku0
      integer :: ierr, nguard0, i, j, block_no
      integer :: lnblocks_wr, tot_blocks_wr, max_lnblocks, ngid
      integer, allocatable :: icheckp_on_cc(:)
      integer, allocatable :: icheckp_on_fc(:,:)
      integer, allocatable :: icheckp_on_ec(:,:)
      integer, allocatable :: icheckp_on_nc(:)
      logical :: l_with_guardcells2
      character (len=80) :: filename
      character (len=6)  :: fnum_string
      real :: user_attr_1_value, & 
     &        user_attr_2_value, & 
     &        user_attr_3_value, & 
     &        user_attr_4_value, & 
     &        user_attr_5_value
      integer :: iorder(4)
      integer :: ordering(1,2,3,4)
      
!---------------------------!
! End variable declarations !
!---------------------------!

      allocate(icheckp_on_cc(nvar))
      allocate(icheckp_on_fc(3,nfacevar))
      allocate(icheckp_on_ec(3,nvaredge))
      allocate(icheckp_on_nc(nvarcorn))

      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

      if (present(l_with_guardcells)) then
         l_with_guardcells2 = l_with_guardcells
      else
         l_with_guardcells2 = .false.
      end if

      nguard0 = nguard*npgs

      if(ndim.lt.3) then
         bnd_box(1,3,:) = 0.
         bnd_box(2,3,:) = 1.
         coord(3,:) = .5*(bnd_box(2,3,:)+bnd_box(1,3,:))
      endif
      if(ndim.lt.2) then
         bnd_box(1,2,:) = 0.
         bnd_box(2,2,:) = 1.
         coord(2,:) = .5*(bnd_box(2,2,:)+bnd_box(1,2,:))
      endif

      lnblocks_wr = lnblocks
      if (lnblocks_wr == 0) lnblocks_wr = lnblocks + 1

      call MPI_ALLREDUCE(lnblocks, tot_blocks, 1, MPI_INTEGER, & 
     &                   MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(lnblocks_wr, tot_blocks_wr, 1, MPI_INTEGER, & 
     &                   MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(lnblocks_wr, max_lnblocks, 1, MPI_INTEGER, & 
     &                   MPI_MAX, MPI_COMM_WORLD, ierr)

! COMPUTE TOTAL NO. OF BLOCKS STORED TO THE 'LEFT' OF THIS PROCESSOR

      if(allocated(n_to_left)) deallocate( n_to_left )
      allocate ( n_to_left(0:nprocs-1) )

      if(allocated(glnblocks)) deallocate( glnblocks )
      allocate ( glnblocks(0:nprocs-1) )

! CEG modified
!      glnblocks(mype) = lnblocks_wr
      call MPI_Allgather(lnblocks_wr, 1,MPI_INTEGER, & 
     &                   glnblocks,1,MPI_INTEGER, & 
     &                   MPI_COMM_WORLD,ierr)
      n_to_left = glnblocks

      do i = nprocs-1,1,-1
         n_to_left(i) = n_to_left(i-1)
      end do

      n_to_left(0) = 0
      do i = 2,nprocs-1
         n_to_left(i) = n_to_left(i) + n_to_left(i-1)
      end do

#ifdef DEBUG
      write(*,*) 'pe ',mype,' n_to_left ',n_to_left
      write(*,*) 'pe ',mype,' tot_blocks ',tot_blocks
#endif /* DEBUG */

! COMPUTE GLOBAL INDIRECT ADDRESSES FOR TREE DATA (gid)

      do block_no = 1,lnblocks

         ngid = 0
         do j = 1,nfaces
            ngid = ngid + 1
            if (neigh(1,j,block_no).gt.0) then
               gid(ngid,block_no) = neigh(1,j,block_no) + & 
     &              n_to_left(neigh(2,j,block_no))
            else
               gid(ngid,block_no) = neigh(1,j,block_no)
            end if
         end do
         
         ngid = ngid + 1
         if (parent(1,block_no).gt.0) then
            gid(ngid,block_no) = parent(1,block_no) + & 
     &           n_to_left(parent(2,block_no))
         else
            gid(ngid,block_no) = parent(1,block_no)
         end if

         do j = 1,nchild
            ngid = ngid + 1
            if (child(1,j,block_no).gt.0) then
               gid(ngid,block_no) = child(1,j,block_no) + & 
     &              n_to_left(child(2,j,block_no))
            else
               gid(ngid,block_no) = child(1,j,block_no)
            end if
         end do
         
      end do
      
      ngid = nfaces + 1 + nchild

      write (fnum_string, '(i6.6)') file_num
      filename = trim(output_dir) //  & 
     &           'paramesh_chk_' //  & 
     &           fnum_string //  & 
     &           '.hdf5'

! set limits on data arrays
      il0 = nguard0
      iu0 = nxb+nguard0
      jl0 = nguard0*k2d
      ju0 = nyb+nguard0*k2d
      kl0 = nguard0*k3d
      ku0 = nzb+nguard0*k3d
      if (.not.no_permanent_guardcells) then
         if(l_with_guardcells2) then
            il0 = 0
            iu0 = nxb+2*nguard0
            jl0 = k2d-k2d
            ju0 = nyb+(2*nguard0)*k2d
            kl0 = k3d-k3d
            ku0 = nzb+(2*nguard0)*k3d
         endif
      end if                    ! no_permanent_guardcells

      nvar_chk_cc =  0
      do i=1,nvar
        icheckp_on_cc(i) = 0
        if(checkp_on_cc(i)) then
           nvar_chk_cc = nvar_chk_cc + 1
           icheckp_on_cc(i) = 1
        end if
      enddo
      nvar_chk_fc =  0
      do i=1,nfacevar
         icheckp_on_fc(:,i) = 0
         if(checkp_on_fc(1,i)) then
           nvar_chk_fc = nvar_chk_fc + 1
           icheckp_on_fc(:,i) = 1
        end if
      enddo
      nvar_chk_ec =  0
      do i=1,nvaredge
        icheckp_on_ec(:,i) = 0
        if(checkp_on_ec(1,i)) then
           nvar_chk_ec = nvar_chk_ec + 1
           icheckp_on_ec(:,i) = 1
        end if
      enddo
      nvar_chk_nc =  0
      do i=1,nvarcorn
        icheckp_on_nc(i) = 0
        if(checkp_on_nc(i)) then
           nvar_chk_nc = nvar_chk_nc + 1
           icheckp_on_nc(i) = 1
        end if
      enddo

! handle user attributes

      user_attr_1_value = 0.
      if (present(user_attr_1)) then
         user_attr_1_value = user_attr_1
      end if

      user_attr_2_value = 0.
      if (present(user_attr_2)) then
         user_attr_2_value = user_attr_2
      end if
         
      user_attr_3_value = 0.
      if (present(user_attr_3)) then
         user_attr_3_value = user_attr_3
      end if
         
      user_attr_4_value = 0.
      if (present(user_attr_4)) then
         user_attr_4_value = user_attr_4
      end if
         
      user_attr_5_value = 0.
      if (present(user_attr_5)) then
         user_attr_5_value = user_attr_5
      end if

      do i = 1,4
         iorder(i) = size(ordering,dim=i)
      end do

         

#ifdef REAL8

      call write_blocks_hdf5_r8(filename, & 
     &                          tot_blocks, & 
     &                          tot_blocks_wr, & 
     &                          max_lnblocks, & 
     &                          lnblocks_wr, & 
     &                          n_to_left(mype), & 
     &                          mdim, & 
     &                          ndim, & 
     &                          ngid, & 
     &                          mflags, & 
     &                          lrefine, & 
     &                          nodetype, & 
     &                          which_child, & 
     &                          gid, & 
     &                          bflags, & 
     &                          coord, & 
     &                          bnd_box, & 
     &                          work_block, & 
     &                          unk, & 
     &                          nvar, & 
     &                          nvar_chk_cc, & 
     &                          icheckp_on_cc, & 
     &                          facevarx, facevary, facevarz, & 
     &                          nbndvar, & 
     &                          nvar_chk_fc, & 
     &                          icheckp_on_fc, & 
     &                          unk_e_x, unk_e_y, unk_e_z, & 
     &                          nbndvare, & 
     &                          nvar_chk_ec, & 
     &                          icheckp_on_ec, & 
     &                          unk_n, & 
     &                          nbndvarc, & 
     &                          nvar_chk_nc, & 
     &                          icheckp_on_nc, & 
     &                          iu_bnd, ju_bnd, ku_bnd, & 
     &                          il0, iu0, & 
     &                          jl0, ju0, & 
     &                          kl0, ku0, & 
     &                          user_attr_1_value, & 
     &                          user_attr_2_value, & 
     &                          user_attr_3_value, & 
     &                          user_attr_4_value, & 
     &                          user_attr_5_value, &
     &                          iorder)

#else

      call write_blocks_hdf5_r4(filename, & 
     &                          tot_blocks, & 
     &                          tot_blocks_wr, & 
     &                          max_lnblocks, & 
     &                          lnblocks_wr, & 
     &                          n_to_left(mype), & 
     &                          mdim, & 
     &                          ndim, & 
     &                          ngid, & 
     &                          mflags, & 
     &                          lrefine, & 
     &                          nodetype, & 
     &                          which_child, & 
     &                          gid, & 
     &                          bflags, & 
     &                          coord, & 
     &                          bnd_box, & 
     &                          work_block, & 
     &                          unk, & 
     &                          nvar, & 
     &                          nvar_chk_cc, & 
     &                          icheckp_on_cc, & 
     &                          facevarx, facevary, facevarz, & 
     &                          nbndvar, & 
     &                          nvar_chk_fc, & 
     &                          icheckp_on_fc, & 
     &                          unk_e_x, unk_e_y, unk_e_z, & 
     &                          nbndvare, & 
     &                          nvar_chk_ec, & 
     &                          icheckp_on_ec, & 
     &                          unk_n, & 
     &                          nbndvarc, & 
     &                          nvar_chk_nc, & 
     &                          icheckp_on_nc, & 
     &                          iu_bnd, ju_bnd, ku_bnd, & 
     &                          il0, iu0, & 
     &                          jl0, ju0, & 
     &                          kl0, ku0, & 
     &                          user_attr_1_value, & 
     &                          user_attr_2_value, & 
     &                          user_attr_3_value, & 
     &                          user_attr_4_value, & 
     &                          user_attr_5_value, &
     &                          iorder)

#endif
      
      deallocate(icheckp_on_cc)
      deallocate(icheckp_on_fc)
      deallocate(icheckp_on_ec)
      deallocate(icheckp_on_nc)

      return
      end subroutine amr_checkpoint_wr_hdf5


