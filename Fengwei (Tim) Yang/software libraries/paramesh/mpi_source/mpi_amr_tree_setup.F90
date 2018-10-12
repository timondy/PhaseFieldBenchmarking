#define ORIG
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

!#define DEBUG

      subroutine mpi_amr_tree_setup(mype,nprocs,tag_offset)


!------------------------------------------------------------------------
!
! This routine is a pre-analysis operation which precedes
! a section of code which requires tree info from remote blocks.
!
!
! Written :     Peter MacNeice          July 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!      tag_offset     defines the last tag number used for an mpi message.
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      Use paramesh_comm_data

      implicit none

      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: tag_offset

      include 'mpif.h'

#ifdef ORIG
      real, dimension (:), allocatable :: send_buf
      real, dimension (:), allocatable :: recv_buf
#endif /* ORIG */
      integer :: buffer_dim,buf_dim_bytes
      integer :: max_blks_sent
      integer :: len_surr_blks
      integer :: itemp
      integer :: no_of_bytes_per_real
      integer :: no_of_bytes_per_integer
      integer :: no_of_bytes_per_logical
      integer :: ierr, ierror

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      itemp = max(sum(commatrix_send),sum(commatrix_recv))
      call MPI_ALLREDUCE (itemp, & 
     &                    max_blks_sent, & 
     &                    1, & 
     &                    MPI_INTEGER, & 
     &                    MPI_MAX, & 
     &                    MPI_COMM_WORLD, & 
     &                    ierror)

!--------------------------------------------------------------
! 
! Establish tree info which is needed about remote blocks

#ifdef DEBUG
      write(*,*)  & 
     &    'pe ',mype,' entered tree setup : max_blks_sent ', & 
     &           max_blks_sent,' offset ',tag_offset
#endif /* DEBUG */

!------
! If we change this buffer size, remember to make the same change in
! mpi_amr_comm_setup
!------
!
! Set up buffer size info for tree messages
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
#ifdef SAVE_MORTS & 
     &                + 6*3*(1+2*k2d)*(1+2*k3d)
#endif 
      buffer_dim = (32+16+len_surr_blks)*max_blks_sent

! this is not a portable computation of buf_dim_bytes
      call mpi_type_size(MPI_INTEGER,no_of_bytes_per_integer,ierr)
      call mpi_type_size(amr_mpi_real, & 
     &                   no_of_bytes_per_real   ,ierr)
      call mpi_type_size(MPI_LOGICAL,no_of_bytes_per_logical,ierr)
#ifdef DEBUG
      write(*,*) 'no_of_bytes_per_integer = ', & 
     &               no_of_bytes_per_integer
      write(*,*) 'no_of_bytes_per_real = ', & 
     &               no_of_bytes_per_real
#endif /* DEBUG */
      buf_dim_bytes =  max_blks_sent* & 
     &    (12*no_of_bytes_per_real  & 
     &     + (16+19+len_surr_blks)*no_of_bytes_per_integer & 
     &     + 1*no_of_bytes_per_logical)

#ifdef REAL8
      buf_dim_bytes =  max_blks_sent* & 
     &    (12*no_of_bytes_per_real  & 
     &     + (16+19+len_surr_blks)*no_of_bytes_per_real & 
     &     + 1*no_of_bytes_per_real)
#endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : buffer_dim ', & 
     &           buffer_dim
#endif /* DEBUG */

      if(allocated(send_buf)) deallocate(send_buf)
      allocate( send_buf(buffer_dim))
      if(allocated(recv_buf)) deallocate(recv_buf)
      allocate( recv_buf(buffer_dim))

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : allocated '
#endif /* DEBUG */

!
! Pack the tree info to be sent
      call mpi_pack_tree_info(mype,nprocs,buf_dim_bytes, & 
     &                     buffer_dim,send_buf)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : exited mpi_pack_blocks'
#endif /* DEBUG */

#ifdef DEBUG
      write(*,*)  & 
     &'pe ',mype,' tree setup : entering mpi_xchange_blocks'
#endif /* DEBUG */

!
! Exchange the tree info to be sent
      call mpi_xchange_blocks(mype,nprocs,tag_offset, & 
     &                        buffer_dim,send_buf, &
     &                        buffer_dim,recv_buf)


#ifdef DEBUG
      write(*,*)  & 
     &'pe ',mype,' tree setup : exited mpi_xchange_blocks'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */

!
! Unpack the tree info which has been received
#ifdef DEBUG
      write(*,*)  & 
     & 'pe ',mype,' tree setup : entering mpi_unpack_treeinfo'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */
      call mpi_unpack_tree_info(mype,nprocs,buf_dim_bytes, & 
     &                       buffer_dim,recv_buf)


#ifdef DEBUG
      write(*,*)  & 
     & 'pe ',mype,' tree setup : exited mpi_unpack_tree_info'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */


      if(allocated(send_buf)) deallocate(send_buf)
      if(allocated(recv_buf)) deallocate(recv_buf)


#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : deallocated '
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */

      return
      end subroutine mpi_amr_tree_setup
