!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* utilities/io/checkpoint/hdf5/amr_checkpoint_re_hdf5
!! NAME
!!
!!   amr_checkpoint_re_hdf5
!!
!! SYNOPSIS
!!
!!   call amr_checkpoint_re_hdf5(file_num)
!!   call amr_checkpoint_re_hdf5(file_num, l_with_guardcells,  
!!                               user_attr1, user_attr2, user_attr3, 
!!                               user_attr4, user_attr5)
!!
!!   call amr_checkpoint_re_hdf5(integer, optional logical, optional char*80,
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
!!   mpi_morton
!!
!! CALLS
!!
!!   read_blocks_hdf5_r8 or read_blocks_hdf5_r4
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a checkpoint file has been read in.
!!
!! DESCRIPTION
!! 
!!  Subroutine to read PARAMESH checkpoint files in parallel which have been 
!!  produced using the routine amr_checkpoint_wr_hdf5.  It reads in the tree data 
!!  structure and data stored in PARAMESH blocks.
!!
!!  The files read in must have names of the form 
!!  'paramesh_chk_######.hdf5'. where '######' is the file_num argument passed into
!!  this routine.
!!
!!  The routine gets called directly by the driver routine 'amr_checkpoint_re' if
!!  the user selects HDF5 format as the data formet for the file.
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

      subroutine amr_checkpoint_re_hdf5 (file_num, & 
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
#ifdef SAVE_MORTS
      use mpi_morton
#endif /*  SAVE_MORTS */

      use paramesh_interfaces, only : amr_morton_order, & 
     &                                amr_guardcell
      use paramesh_mpi_interfaces, only : mpi_amr_global_domain_limits, & 
     &                                    mpi_amr_boundary_block_info

!---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'

      integer, intent(in) :: file_num
      logical, optional, intent(in)  :: l_with_guardcells
      real, optional, intent(out) :: user_attr_1, & 
     &                               user_attr_2, & 
     &                               user_attr_3, & 
     &                               user_attr_4, & 
     &                               user_attr_5

      integer,dimension (:), allocatable :: n_to_left
      integer :: gid(nfaces+1+nchild,maxblocks_tr)
      integer :: proc, icount
      integer :: nvar_chk_cc,nvar_chk_fc,nvar_chk_ec,nvar_chk_nc
      integer :: tot_blocks
      integer :: mype, nprocs
      integer :: il0, iu0, jl0, ju0, kl0, ku0
      integer :: ierr, nguard0, i, j, block_no, lb
      integer :: ngid, lnblocks_old
      integer, allocatable :: icheckp_on_cc(:)
      integer, allocatable :: icheckp_on_fc(:,:)
      integer, allocatable :: icheckp_on_ec(:,:)
      integer, allocatable :: icheckp_on_nc(:)
      logical :: l_with_guardcells2, l_move_solution
      double precision :: time1
      real :: user_attr_1_value, & 
     &        user_attr_2_value, & 
     &        user_attr_3_value, & 
     &        user_attr_4_value, & 
     &        user_attr_5_value
#ifdef SAVE_MORTS
      integer :: lb
      integer :: mort_neigh(6,3,3,3)
      real    :: xmin,ymin,zmin,xmax,ymax,zmax
#endif /* SAVE_MORTS */

      character (len=80) :: filename
      character (len=6)  :: fnum_string

      integer :: iorder(4)
      integer :: ordering(1,2,3,4)

!---------------------------------------------------------------------------

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
            ju0 = nyb+(2*nguard0*k2d)
            kl0 = k3d-k3d
            ku0 = nzb+(2*nguard0*k3d)
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

      ngid = nfaces + 1 + nchild

      do i = 1, 4
         iorder(i) = size(ordering,dim=i)
      end do

#ifdef REAL8

      call read_blocks_hdf5_r8(filename, & 
     &                         tot_blocks, & 
     &                         lnblocks, & 
     &                         mdim, & 
     &                         ndim, & 
     &                         ngid, & 
     &                         mflags, & 
     &                         lrefine, & 
     &                         nodetype, & 
     &                         which_child, & 
     &                         gid, & 
     &                         bflags, & 
     &                         coord, & 
     &                         bnd_box, & 
     &                         work_block, & 
     &                         unk, & 
     &                         nvar, & 
     &                         nvar_chk_cc, & 
     &                         icheckp_on_cc, & 
     &                         facevarx, facevary, facevarz, & 
     &                         nbndvar, & 
     &                         nvar_chk_fc, & 
     &                         icheckp_on_fc, & 
     &                         unk_e_x, unk_e_y, unk_e_z, & 
     &                         nbndvare, & 
     &                         nvar_chk_ec, & 
     &                         icheckp_on_ec, & 
     &                         unk_n, & 
     &                         nbndvarc, & 
     &                         nvar_chk_nc, & 
     &                         icheckp_on_nc, & 
     &                         iu_bnd, ju_bnd, ku_bnd, & 
     &                         il0, iu0, & 
     &                         jl0, ju0, & 
     &                         kl0, ku0, & 
     &                         user_attr_1_value, & 
     &                         user_attr_2_value, & 
     &                         user_attr_3_value, & 
     &                         user_attr_4_value, & 
     &                         user_attr_5_value, &
     &                         iorder)

#else

      call read_blocks_hdf5_r4(filename, & 
     &                         tot_blocks, & 
     &                         lnblocks, & 
     &                         mdim, & 
     &                         ndim, & 
     &                         ngid, & 
     &                         mflags, & 
     &                         lrefine, & 
     &                         nodetype, & 
     &                         which_child, & 
     &                         gid, & 
     &                         bflags, & 
     &                         coord, & 
     &                         bnd_box, & 
     &                         work_block, & 
     &                         unk, & 
     &                         nvar, & 
     &                         nvar_chk_cc, & 
     &                         icheckp_on_cc, & 
     &                         facevarx, facevary, facevarz, & 
     &                         nbndvar, & 
     &                         nvar_chk_fc, & 
     &                         icheckp_on_fc, & 
     &                         unk_e_x, unk_e_y, unk_e_z, & 
     &                         nbndvare, & 
     &                         nvar_chk_ec, & 
     &                         icheckp_on_ec, & 
     &                         unk_n, & 
     &                         nbndvarc, & 
     &                         nvar_chk_nc, & 
     &                         icheckp_on_nc, & 
     &                         iu_bnd, ju_bnd, ku_bnd, & 
     &                         il0, iu0, & 
     &                         jl0, ju0, & 
     &                         kl0, ku0, & 
     &                         user_attr_1_value, & 
     &                         user_attr_2_value, & 
     &                         user_attr_3_value, & 
     &                         user_attr_4_value, & 
     &                         user_attr_5_value, &
     &                         iorder)

#endif

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

      do block_no = 1,lnblocks
         bsize(:,block_no) = bnd_box(2,:,block_no)- & 
     &        bnd_box(1,:,block_no)
      enddo

! COMPUTE TREE DATA FROM gid

      if(allocated(n_to_left)) deallocate( n_to_left )
      allocate ( n_to_left(0:nprocs-1) )

      proc = 0
      icount = 0
      n_to_left(:) = 0
      do while (icount < tot_blocks)
         if (proc < mype) then
            n_to_left(mype) = n_to_left(mype) + 1
         end if
         proc = proc + 1
         if (proc > nprocs-1) proc = 0
         icount = icount + 1
      end do

! CEG modified
      sendmsg = n_to_left(mype)
      call MPI_Allgather(sendmsg, 1,MPI_INTEGER, & 
     &                   n_to_left,1,MPI_INTEGER, & 
     &                   MPI_COMM_WORLD,ierr)

      do block_no = 1,lnblocks

! neighbor data
         ngid = 0
         do j = 1,nfaces
            ngid = ngid + 1
            if (gid(ngid,block_no).gt.0) then
               proc = 0
               do while (n_to_left(proc) < gid(ngid,block_no))
                  proc = proc + 1
                  if (proc >= nprocs) exit
               end do
               proc = proc - 1
               neigh(2,j,block_no) = proc 
               neigh(1,j,block_no) = gid(ngid,block_no) - & 
     &                               n_to_left(proc)
            else
               neigh(1,j,block_no) = gid(ngid,block_no)
               neigh(2,j,block_no) = gid(ngid,block_no)
            end if
         end do


! parent data
         ngid = ngid + 1
         if (gid(ngid,block_no).gt.0) then
            proc = 0
            do while (n_to_left(proc) < gid(ngid,block_no))
               proc = proc + 1
               if (proc >= nprocs) exit
            end do
            proc = proc - 1
            parent(2,block_no) = proc
            parent(1,block_no) = gid(ngid,block_no) - & 
     &                           n_to_left(proc)
         else
            parent(1,block_no) = gid(ngid,block_no)
            parent(2,block_no) = gid(ngid,block_no)
         end if

! children data
         do j = 1,nchild
            ngid = ngid + 1
            if (gid(ngid,block_no).gt.0) then
               proc = 0
               do while (n_to_left(proc) < gid(ngid,block_no))
                  proc = proc + 1
                  if (proc >= nprocs) exit
               end do
               proc = proc - 1
               child(2,j,block_no) = proc
               child(1,j,block_no) = gid(ngid,block_no) - & 
     &                               n_to_left(proc)
            else
               child(1,j,block_no) = gid(ngid,block_no)
               child(2,j,block_no) = gid(ngid,block_no)
            end if
         end do

      end do

! Now reorder blocks such that they are better balanced
! NOTE: this assumes that the total number of blocks is > nprocs

! NOTE: We cannot do a morton ordering here if l_with_guardcells is defined since
!       amr_redist_blks do not move the guardcells.

      if (.not.l_with_guardcells2) then  

      lnblocks_old = lnblocks
      l_move_solution = .true.
      call amr_morton_order (lnblocks_old,nprocs,mype, & 
     &                       l_move_solution)

      write(*,*) 'after amr_morton_order : pe lnb ',mype,lnblocks

      end if

!---------------------------------------------
! compute grid_xmax, etc
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call mpi_amr_global_domain_limits
      if (timing_mpi) then
      timer_amr_global_domain_limits = & 
     &     timer_amr_global_domain_limits + mpi_wtime() - time1
      endif

#ifdef SAVE_MORTS
! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

      write(*,*) 'checkre xmin etc ',xmin,ymin,zmin,xmax,ymax,zmax
      write(*,*) 'checkre lperiodicxtc ', & 
     &                        lperiodicx,lperiodicy,lperiodicz
      do lb = 1,lnblocks
        call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                        lperiodicx,lperiodicy,lperiodicz, & 
     &                        coord(:,lb),bsize(:,lb),ndim, & 
     &                        lrefine(lb),lrefine_max,mort_neigh)
        surr_morts(:,:,:,:,lb) = & 
     &                    mort_neigh(:,:,2-k2d:2+k2d,2-k3d:2+k3d)
      write(*,*) 'check_re set up surr_morts for initial blk ',lb,mype, & 
     &            ' surr_morts(6,:,:,1,lb) ',surr_morts(6,:,:,1,lb)
      enddo
#endif /* SAVE_MORTS */

      call amr_morton_process()

!
! St up an array of cell sizes for each grid refinement level.
! These can be used to minimize variation due to roundoff, but
! should ONLY be used with a uniformly spaced grid.
      level_cell_sizes = 0.
      level_cell_sizes(1,1) = (grid_xmax-grid_xmin)/real(nxb)
      if(ndim.gt.1) & 
     &  level_cell_sizes(2,1) = (grid_ymax-grid_ymin)/real(nyb)
      if(ndim.eq.3) & 
     &  level_cell_sizes(3,1) = (grid_zmax-grid_zmin)/real(nzb)
      do i=2,lrefine_max
        level_cell_sizes(1:ndim,i) = .5*level_cell_sizes(1:ndim,i-1)
      enddo
!---------------------------------------------

!
! mark grid as changed
      grid_changed = 1
      grid_analysed_mpi = 1


! Now make sure guardcell information is up to date

      if (.not.l_with_guardcells2)  & 
     &     call amr_guardcell(mype,1,nguard)

      call mpi_amr_boundary_block_info(mype,nprocs)

      deallocate(n_to_left)

      deallocate(icheckp_on_cc)
      deallocate(icheckp_on_fc)
      deallocate(icheckp_on_ec)
      deallocate(icheckp_on_nc)

! handle user attributes

      if (present(user_attr_1)) then
         user_attr_1 = user_attr_1_value
      end if

      if (present(user_attr_2)) then
         user_attr_2 = user_attr_2_value
      end if
         
      if (present(user_attr_3)) then
         user_attr_3 = user_attr_3_value
      end if
         
      if (present(user_attr_4)) then
         user_attr_4 = user_attr_4_value
      end if
         
      if (present(user_attr_5)) then
         user_attr_5 = user_attr_5_value
      end if

      return
      end subroutine amr_checkpoint_re_hdf5
