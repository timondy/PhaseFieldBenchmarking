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

      subroutine mpi_unpack_fluxes(mype, & 
     &                             buf_dim,R_buffer,flux_dir)

!------------------------------------------------------------------------
!
! This subroutine unpacks all fluxes which are to be received on mype.
! It further stores the local (receiving) block id, the neighboring remote
! (sending) block id, and the local guard block id into the array laddress
! which is to be used in the subroutine mpi_1blk_guardcell
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      nprocs         number of processors
!      buf_dim        dimension of buffer
!      R_buffer       receive buffer 
!      ir_buf         starting & ending indices of buffer
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_put_flux_buffer

      implicit none

      include 'mpif.h'

      integer, intent(in) :: mype,buf_dim
      real,    intent(inout) ::  R_buffer(buf_dim)
      integer, optional, intent(in) :: flux_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lblk, lnumb, lb
      integer :: index
      integer :: ierrorcode,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#define DEBUG

      lnumb = sum(commatrix_recv(:))
      if(lnumb.gt.maxblocks_alloc) then
            call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      index = 1

 
      do lblk=1,lnumb
        lb = lblk + strt_buffer - 1

                                  ! unpack all arrays from buffer into lb
#ifdef DEBUG
        write(*,*) 'pe ',mype,' lblk ',lblk,' unpacking starting ', & 
     &        ' at index ',index
#endif /* DEBUG */
        call mpi_put_flux_buffer(mype, & 
     &         lb,index,buf_dim,R_buffer,flux_dir)
#ifdef DEBUG
        write(*,*) 'pe ',mype,' lblk ',lblk,' unpacked into ',lb
#endif /* DEBUG */

      enddo


#ifdef DEBUG
      if (index .ne. buf_dim) then
      print *,' '
      print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
      print *,' Rbuffer is too big by ',size(R_buffer,dim=1)-index,mype
      print *,' index = ',index,' buf_dim = ',buf_dim,mype
      print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
      print *,' '
      end if
#endif

      return
      end subroutine mpi_unpack_fluxes
