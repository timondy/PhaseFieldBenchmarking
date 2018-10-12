!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f, ordering
#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_plotfile_chombo (file_num)


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

      integer, intent(in) :: file_num

      integer,dimension (:), allocatable :: n_to_left
      integer,dimension (:,:), allocatable :: n_to_left_level
      integer,dimension (:), allocatable :: n_to_left_level2
      integer,dimension (:), allocatable :: glnblocks
      integer :: nvar_chk_cc,nvar_chk_fc,nvar_chk_ec,nvar_chk_nc
      integer :: num_components, num_levels
      integer :: tot_blocks
      integer :: mype, nprocs
      integer :: il0, iu0, jl0, ju0, kl0, ku0
      integer :: ierr, nguard0, i, j, block_no
      integer :: lnblocks_wr, tot_blocks_wr, max_lnblocks 
      integer :: icount, icount2
      integer :: minlevel, maxlevel
      integer :: tmp_int
      integer, dimension(:), allocatable :: no_at_level
      integer, allocatable :: icheckp_on_cc(:)
      integer, allocatable :: icheckp_on_fc(:,:)
      integer, allocatable :: icheckp_on_ec(:,:)
      integer, allocatable :: icheckp_on_nc(:)
      character (len=80) :: filename
      character (len=4)  :: fnum_string
      character (len=20),allocatable :: compNames(:)

      integer :: ordering(1,2,3,4)
      integer :: iorder(4)
      integer :: sendmsg

!---------------------------!
! End variable declarations !
!---------------------------!

      allocate(icheckp_on_cc(nvar))
      allocate(icheckp_on_fc(3,nfacevar))
      allocate(icheckp_on_ec(3,nvaredge))
      allocate(icheckp_on_nc(nvarcorn))

      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

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

! number of refinement levels

      maxlevel = -1
      do i = 1, lnblocks
         if (lrefine(i) > maxlevel) maxlevel = lrefine(i)
      end do

      call MPI_ALLREDUCE (minlevel, tmp_int, 1, & 
     &     MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      minlevel = tmp_int

      call MPI_ALLREDUCE (maxlevel, tmp_int, 1, & 
     &     MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      maxlevel = tmp_int

      num_levels = maxlevel

      allocate(no_at_level(num_levels))
      
      no_at_level(:) = 0
      do i = 1, lnblocks
         no_at_level(lrefine(i)) = no_at_level(lrefine(i)) + 1
      end do

      do i = 1,num_levels
         call MPI_ALLREDUCE (no_at_level(i), tmp_int, 1, & 
     &        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         no_at_level(i) = tmp_int
      end do


! COMPUTE TOTAL NO. OF BLOCKS STORED TO THE 'LEFT' OF THIS PROCESSOR

      if(allocated(n_to_left)) deallocate( n_to_left )
      allocate ( n_to_left(0:nprocs-1) )

      if(allocated(n_to_left_level)) deallocate( n_to_left_level )
      allocate ( n_to_left_level(0:nprocs-1,num_levels) )

      if(allocated(n_to_left_level2)) deallocate( n_to_left_level2 )
      allocate ( n_to_left_level2(num_levels) )

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

      n_to_left_level(:,:) = 0
      do i = 1, lnblocks
         n_to_left_level(mype,lrefine(i)) =  & 
     &        n_to_left_level(mype,lrefine(i)) + 1
      end do
      do i = 1, num_levels
! CEG modified
         sendmsg = n_to_left_level(mype,i)
         call MPI_Allgather(sendmsg, 1,MPI_INTEGER, & 
     &                      n_to_left_level(:,i),1,MPI_INTEGER, & 
     &                      MPI_COMM_WORLD,ierr)
         do j = nprocs-1,1,-1
            n_to_left_level(j,i) = n_to_left_level(j-1,i)
         end do
         n_to_left_level(0,i) = 0
         do j = 2,nprocs-1
            n_to_left_level(j,i) = n_to_left_level(j,i) +  & 
     &                             n_to_left_level(j-1,i)
         end do
      end do

      n_to_left_level2(:) = n_to_left_level(mype,:)

      if (sum(n_to_left_level2(:)) .ne. n_to_left(mype)) then
         print *,' ERROR sum inconsistent '
      end if

! CREATE output file name

      write (fnum_string, '(i4.4)') file_num
      filename = trim(output_dir) //  & 
     &           'paramesh_chombo_' //  & 
     &           fnum_string //  & 
     &           '.hdf5'

! set limits on data arrays
      il0 = nguard0
      iu0 = nxb+nguard0
      jl0 = nguard0*k2d
      ju0 = nyb+nguard0*k2d
      kl0 = nguard0*k3d
      ku0 = nzb+nguard0*k3d

      num_components = 0
      nvar_chk_cc =  0
      do i=1,nvar
        icheckp_on_cc(i) = 0
        if(checkp_on_cc(i)) then
           nvar_chk_cc = nvar_chk_cc + 1
           num_components = num_components + 1
           icheckp_on_cc(i) = 1
        end if
      enddo

      nvar_chk_fc =  0
      do i=1,nfacevar
        icheckp_on_fc(:,i) = 0
        if(checkp_on_fc(1,i)) then
           nvar_chk_fc = nvar_chk_fc + 1
           icheckp_on_fc(:,i) = 1
           num_components = num_components + 3
        end if
      enddo

      nvar_chk_ec =  0
      do i=1,nvaredge
        icheckp_on_ec(:,i) = 0
        if(checkp_on_ec(1,i)) then
           nvar_chk_ec = nvar_chk_ec + 1
           icheckp_on_ec(:,i) = 1
           num_components = num_components + 3
        end if
      enddo

      nvar_chk_nc =  0
      do i=1,nvarcorn
        icheckp_on_nc(i) = 0
        if(checkp_on_nc(i)) then
           nvar_chk_nc = nvar_chk_nc + 1
           icheckp_on_nc(i) = 1
           num_components = num_components + 1
        end if
      enddo


! component names (need to changes this to allow user to define)

      allocate(compNames(num_components))

      icount = 0

      icount2 = 0
      do i=1,nvar
        if(checkp_on_cc(i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'unk_' // fnum_string
        end if
      enddo

      icount2 = 0
      do i=1,nfacevar
        if(checkp_on_fc(1,i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'facevarx_' // fnum_string
        end if
      end do

      icount2 = 0
      do i=1,nfacevar
        if(checkp_on_fc(2,i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'facevary_' // fnum_string
        end if
      end do

      icount2 = 0
      do i=1,nfacevar
        if(checkp_on_fc(3,i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'facevarz_' // fnum_string
        end if
      enddo

      icount2 = 0
      do i=1,nvaredge
        if(checkp_on_ec(1,i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'unk_e_x_' // fnum_string
        end if
      enddo

      icount2 = 0
      do i=1,nvaredge
        if(checkp_on_ec(2,i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'unk_e_y_' // fnum_string
        end if
      enddo

      icount2 = 0
      do i=1,nvaredge
        if(checkp_on_ec(3,i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'unk_e_z_' // fnum_string
        end if
      enddo

      icount2 = 0
      do i=1,nvarcorn
        if(checkp_on_nc(i)) then
           icount = icount + 1
           icount2 = icount2 + 1
           write (fnum_string, '(i4.4)') icount2
           compNames(icount) = 'unk_n_' // fnum_string
        end if
      enddo
      
      do i = 1,4
         iorder(i) = size(ordering,dim=i)
      end do

#ifdef REAL8

      call write_blocks_chombo_r8(filename, & 
     &                            num_components, & 
     &                            num_levels, & 
     &                            compNames, & 
     &                            ndim, & 
     &                            nxb, nyb, nzb, & 
     &                            no_at_level, & 
     &                          tot_blocks, & 
     &                          tot_blocks_wr, & 
     &                          max_lnblocks, & 
     &                          lnblocks_wr, & 
     &                          n_to_left(mype), & 
     &                          n_to_left_level2, & 
     &                          mdim, & 
     &                          mflags, & 
     &                          lrefine, & 
     &                          nodetype, & 
     &                          which_child, & 
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
     &                          iorder)

#else

      call write_blocks_chombo_r4(filename, & 
     &                            num_components, & 
     &                            num_levels, & 
     &                            compNames, & 
     &                            ndim, & 
     &                            nxb, nyb, nzb, & 
     &                            no_at_level, & 
     &                          tot_blocks, & 
     &                          tot_blocks_wr, & 
     &                          max_lnblocks, & 
     &                          lnblocks_wr, & 
     &                          n_to_left(mype), & 
     &                          n_to_left_level2, & 
     &                          mdim, & 
     &                          mflags, & 
     &                          lrefine, & 
     &                          nodetype, & 
     &                          which_child, & 
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
     &                          iorder)

#endif

      deallocate(compNames)
      deallocate(no_at_level)
      deallocate(n_to_left)
      deallocate(n_to_left_level)
      deallocate(n_to_left_level2)
      deallocate(glnblocks)
      deallocate(icheckp_on_cc)
      deallocate(icheckp_on_fc)
      deallocate(icheckp_on_ec)
      deallocate(icheckp_on_nc)
      
      return
      end subroutine amr_plotfile_chombo

