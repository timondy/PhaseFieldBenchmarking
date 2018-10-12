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

!#define DEBUGX
!#define DEBUGY
!#define DEBUG


        subroutine comm_start(MaxProcs,nprocs,mype)

        use paramesh_mpi_interfaces, only : mpi_array_allocate

        implicit none

        integer, intent(out) :: nprocs,mype
        integer, intent(in)  :: MaxProcs

        include 'mpif.h'
        integer :: ierr
        logical :: initialized

        call mpi_initialized(initialized, ierr)
        if (.not.initialized) then
          call mpi_init(ierr)
#ifdef AUTOPACK
          call AP_INIT
#endif /* AUTOPACK */
        end if

        Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
        Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)

        call mpi_array_allocate(nprocs)

        return
        end subroutine comm_start

        subroutine comm_finish()
        use paramesh_mpi_interfaces, only : mpi_array_deallocate

        implicit none

        integer :: ierr

        call mpi_array_deallocate
#ifdef AUTOPACK
          call AP_FINALIZE
#endif /* AUTOPACK */
        call mpi_finalize(ierr)

        return
        end subroutine comm_finish





        subroutine comm_logical_or_to_all(target,source)
        implicit none
        logical, intent(in)  :: source
        logical, intent(out) :: target
        integer :: nred, ierr
        common/comm_lcommon/ reduce_datain,reduce_dataout
        logical :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_logical_allreduce(reduce_datain(1),reduce_dataout(1), & 
     &                     nred, & 
     &                     MPI_LOGICAL, & 
     &                     MPI_LOR,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)
        return
        end subroutine comm_logical_or_to_all

        subroutine comm_real_sum_to_all(target,source)
        Use paramesh_comm_data
        implicit none
        real, intent(in)  :: source
        real, intent(out) :: target
        integer :: nred, ierr
        common/comm_rcommon/ reduce_datain,reduce_dataout
        real :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_real_allreduce(reduce_datain(1),reduce_dataout(1), & 
     &                     nred, & 
     &                     amr_mpi_real, & 
     &                     MPI_SUM,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)
        return
        end subroutine comm_real_sum_to_all


        subroutine comm_dble_sum_to_all(target,source)
        implicit none
        double precision, intent(in)  :: source
        double precision, intent(out) :: target
        integer :: nred, ierr
        common/comm_dcommon/ reduce_datain,reduce_dataout
        double precision :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_dble_allreduce(reduce_datain(1),reduce_dataout(1), & 
     &                          nred, & 
     &                          MPI_DOUBLE_PRECISION, & 
     &                          MPI_SUM,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)
        return
        end subroutine comm_dble_sum_to_all





        subroutine comm_real_min_to_all(target,source)
        Use paramesh_comm_data
        implicit none
        real, intent(in)  :: source
        real, intent(out) :: target
        integer :: nred, ierr
        common/comm_rcommon/ reduce_datain,reduce_dataout
        real :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_real_allreduce(reduce_datain(1),reduce_dataout(1), & 
     &                     nred, & 
     &                     amr_mpi_real, & 
     &                     MPI_MIN,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)

        return
        end subroutine comm_real_min_to_all





        subroutine comm_real_max_to_all(target,source)
        Use paramesh_comm_data
        implicit none
        real, intent(in)  :: source
        real, intent(out) :: target
        integer :: nred, ierr
        common/comm_rcommon/ reduce_datain,reduce_dataout
        real :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_real_allreduce(reduce_datain(1),reduce_dataout(1), & 
     &                     nred, & 
     &                     amr_mpi_real, & 
     &                     MPI_MAX,MPI_COMM_WORLD,ierr)


        target = reduce_dataout(1)

        return
        end subroutine comm_real_max_to_all






        subroutine comm_int_sum_to_all(target,source)
        implicit none
        integer, intent(in)  :: source
        integer, intent(out) :: target
        integer :: nred, ierr
        common/comm_icommon/ ireduce_datain,ireduce_dataout
        integer :: ireduce_datain(1),ireduce_dataout(1)

        include 'mpif.h'

        nred = 1
        ireduce_datain(1) = source
        call mpi_int_allreduce(ireduce_datain(1),ireduce_dataout(1), & 
     &                  nred, & 
     &                  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        target = ireduce_dataout(1)

        return
        end subroutine comm_int_sum_to_all



        subroutine comm_int_min_to_all(target,source)
        implicit none
        integer, intent(in)  :: source
        integer, intent(out) :: target
        integer :: nred, ierr
        common/comm_icommon/ ireduce_datain,ireduce_dataout
        integer :: ireduce_datain(1),ireduce_dataout(1)

        include 'mpif.h'

        nred = 1
        ireduce_datain(1) = source
        call mpi_int_allreduce(ireduce_datain(1),ireduce_dataout(1), & 
     &                  nred, & 
     &                  MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

        target = ireduce_dataout(1)

        return
        end subroutine comm_int_min_to_all






        subroutine comm_int_max_to_all(target,source)
        implicit none
        integer, intent(in)  :: source
        integer, intent(out) :: target
        integer :: nred, ierr
        common/comm_icommon/ ireduce_datain,ireduce_dataout
        integer :: ireduce_datain(1),ireduce_dataout(1)

        include 'mpif.h'

        nred = 1
        ireduce_datain(1) = source
        call mpi_int_allreduce(ireduce_datain(1),ireduce_dataout(1), & 
     &                  nred, & 
     &                  MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        target = ireduce_dataout(1)


        return
        end subroutine comm_int_max_to_all


      subroutine mpi_setup(mype,nprocs)


!------------------------------------------------------------------------
!
! This routine sets up some variables which are needed to control
! the MPI message passing. 
! This routine should be called every time the grid is rebuilt.
!
!
! Written :     Peter MacNeice and Michael Gehmeyr          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use constants

      implicit none

      integer, intent(in) :: mype,nprocs

      include 'mpif.h'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lb
      logical :: l_periodx,l_periody,l_periodz
      integer :: ierr, ierror

      logical :: lreduce_datain(1),lreduce_dataout(1)
      real    :: eps,accuracy

      if (spherical_pm) then
        accuracy = 100./10.**precision(accuracy)
        if (accuracy > 1.0e-10) then
           eps = 1.e-10
        else
           eps = accuracy
        end if
      end if

!
! Set up error handler
       call mpi_errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN, & 
     &                         ierr)

! Detect if grid is periodic in any direction
       l_periodx = .false.
       l_periody = .false.
       l_periodz = .false.
       do lb = 1,lnblocks
         if( (bnd_box(1,1,lb)-.1*bsize(1,lb).lt.grid_xmin) & 
     &                           .and. & 
     &       (neigh(1,1,lb).gt.0) ) & 
     &        l_periodx = .true.
         if( (bnd_box(1,2,lb)-.1*bsize(2,lb).lt.grid_ymin) & 
     &       .and. (neigh(1,3,lb).gt.0) .and. (ndim.ge.2) ) & 
     &        l_periody = .true.
         if( (bnd_box(1,3,lb)-.1*bsize(3,lb).lt.grid_zmin) & 
     &       .and. (neigh(1,5,lb).gt.0) .and. (ndim.eq.3) ) & 
     &        l_periodz = .true.
       enddo

       if(spherical_pm) then
         l_periody = .false.
         if(abs(grid_zmin).lt.eps.and.abs(grid_zmax-pi).lt.eps) & 
     &        l_periodz = .true.
       endif

       lreduce_datain(1) = l_periodx
       call mpi_logical_allreduce(lreduce_datain(1),lreduce_dataout(1), & 
     &                 1, MPI_LOGICAL, & 
     &                 MPI_LOR,MPI_COMM_WORLD,ierror)
       l_periodx = lreduce_dataout(1)
       lreduce_datain(1) = l_periody
       call mpi_logical_allreduce(lreduce_datain(1),lreduce_dataout(1), & 
     &                 1, MPI_LOGICAL, & 
     &                 MPI_LOR,MPI_COMM_WORLD,ierror)
       l_periody = lreduce_dataout(1)
       lreduce_datain(1) = l_periodz
       call mpi_logical_allreduce(lreduce_datain(1),lreduce_dataout(1), & 
     &                 1, MPI_LOGICAL, & 
     &                 MPI_LOR,MPI_COMM_WORLD,ierror)
       l_periodz = lreduce_dataout(1)


        lperiodicx = l_periodx
        lperiodicy = l_periody
        lperiodicz = l_periodz

!       write(*,*) 'pe ',mype,' lperiodicx ',lperiodicx
!       write(*,*) 'pe ',mype,' lperiodicy ',lperiodicy
!       write(*,*) 'pe ',mype,' lperiodicz ',lperiodicz


! Initialize logical array which records the type of data in the last
! message which was packed.
        l_datapacked = .false.

      return
      end subroutine mpi_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine pf_mpi_setup(mype,nprocs)
      use mpi_morton
      implicit none
      integer, intent(in) :: mype,nprocs
                                          
!      lperiodicx = .false.
!      lperiodicy = .false.
!      lperiodicz = .false.
      l_datapacked = .false.
      return
      end subroutine pf_mpi_setup              


#ifdef MOVED
      subroutine send_real(rv,n,idest,itag,icom,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for real arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      rv             real valued array
!      n              dimension of rv
!      idest          MPI destination of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      Use paramesh_comm_data
      implicit none
      include 'mpif.h'
    
      real rv(*)
      integer n,idest,itag,icom,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_SEND(rv,n,amr_mpi_real, & 
     &                   idest,itag,icom,ierr)

      return
      end subroutine send_real

      subroutine send_int(iv,n,idest,itag,icom,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for integer arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      iv             integer valued array
!      n              dimension of iv
!      idest          MPI destination of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
     
      integer iv(*)
      integer n,idest,itag,icom,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_SEND(iv,n,MPI_INTEGER, & 
     &                   idest,itag,icom,ierr)

      return
      end subroutine send_int

      subroutine recv_real(rv,n,isource,itag,icom,istatus,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for integer arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      rv             real valued array
!      n              dimension of rv
!      isource        MPI source of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      istatus        MPI status flag
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      Use paramesh_comm_data
      implicit none
      include 'mpif.h'
     
      real rv(*)
      integer n,isource,itag,icom,ierr
      integer istatus(MPI_STATUS_SIZE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_RECV(rv,n,amr_mpi_real, & 
     &                   isource,itag,icom,istatus,ierr)

      return
      end subroutine recv_real

      subroutine recv_int(iv,n,isource,itag,icom,istatus,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for integer arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      iv             integer valued array
!      n              dimension of iv
!      isource        MPI source of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      istatus        MPI status flag
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
     
      integer iv(*)
      integer n,isource,itag,icom,ierr
      integer istatus(MPI_STATUS_SIZE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_RECV(iv,n,MPI_INTEGER, & 
     &                   isource,itag,icom,istatus,ierr)

      return
      end subroutine recv_int


#endif /*  MOVED */

      subroutine mpi_array_allocate(nprocs)


!------------------------------------------------------------------------
!
! Some of the arrays which are used by the routines which manage MPI
! communications need to be dimensioned with sizes determined at
! run time. These arrays are declared as allocatable in  the
! module mpi_morton, and the actual allocation is performed here.
!
! Written :     Peter MacNeice          June 2000
!------------------------------------------------------------------------
!
! Arguments:
!      nprocs         number of processors
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton


      integer, intent(in) :: nprocs


!------------------------------------------------------------------------

      allocate( pe_source(1:nprocs)             )
      allocate( morton_limits(6,1:2,1:2,1:nprocs) )
      allocate( commatrix_send(1:nprocs)        )
      allocate( commatrix_recv(1:nprocs)        )
      allocate( ir_buf(2,1:nprocs)              )
      allocate( is_buf(2,1:nprocs)              )


      return
      end subroutine mpi_array_allocate



      subroutine mpi_array_deallocate


!------------------------------------------------------------------------
!
! Some of the arrays which are used by the routines which manage MPI
! communications need to be dimensioned with sizes determined at
! run time. These arrays are declared as allocatable in the
! module mpi_morton. This routine deallocates them for a clean
! exit from the program.
!
! Written :     Peter MacNeice          June 2000
!------------------------------------------------------------------------
!
! Arguments:
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

!------------------------------------------------------------------------

      if(allocated(pe_source)) deallocate( pe_source     )
      if(allocated(morton_limits)) deallocate( morton_limits )
      if(allocated(commatrix_send)) deallocate( commatrix_send)
      if(allocated(commatrix_recv)) deallocate( commatrix_recv)
      if(allocated(ir_buf)) deallocate( ir_buf        )
      if(allocated(is_buf)) deallocate( is_buf        )


      return
      end subroutine mpi_array_deallocate


      subroutine morton_number( & 
     &          x0,y0,z0,bbsize,ndim,lrefine_max,lrefine, & 
     &          mort)

!---------------------------------------------------------------------
!
! This subroutine computes a morton number for a given set of
! coordinates and refinement level.
!
! Written :     December 2002     Peter MacNeice
! Modified for new morton numbering scheme: March 2003  Kevin Olson
!---------------------------------------------------------------------

      integer,intent(in ) :: lrefine,lrefine_max,ndim
      real,intent(in)     :: bbsize(3),x0,y0,z0
      integer,intent(out) :: mort(6)
      logical, save :: first = .true.
      integer, save :: nbits,nbits2
!------------------------------------------------------------

      integer, parameter :: ikind = selected_int_kind(18)
      integer (kind=selected_int_kind(18)) :: ix, iy, iz
      integer :: ipos2, j, ipos, nbitshft, nshifts, k

      if (first) then
       first = .false.
       nbits = bit_size(mort(1))
       if (nbits == 32) then
          nbits = nbits-2
       elseif (nbits == 64) then
          nbits = nbits-4
       end if
       nbits2 = bit_size(ix)-4
       if (ndim == 1) nbits2 = nbits
      end if

! compute morton number
      ix = int(x0/bbsize(1),kind=ikind)
      if (ndim >= 2) then
         iy = int(y0/bbsize(2),kind=ikind)
      else
         iy = 0
      end if
      if (ndim == 3) then
         iz = int(z0/bbsize(3),kind=ikind)
      else
         iz = 0
      end if
     
      mort(:) = 0

      ipos2 = 0
      j = 6
      do while (ipos2 < nbits2-1)

      ipos = 0
      do while (ipos <= nbits-ndim .and. ipos2 <= nbits2-1)

         if (btest(ix,ipos2)) then
            mort(j) = ibset(mort(j),ipos)
         else
            mort(j) = ibclr(mort(j),ipos)
         end if
         ipos = ipos + 1
         if (ndim >= 2) then
            if (btest(iy,ipos2)) then
               mort(j) = ibset(mort(j),ipos)
            else
               mort(j) = ibclr(mort(j),ipos)
            end if
            ipos = ipos + 1
         end if
         if (ndim == 3) then
            if (btest(iz,ipos2)) then
               mort(j) = ibset(mort(j),ipos)
            else
               mort(j) = ibclr(mort(j),ipos)
            end if
            ipos = ipos + 1
         end if

         ipos2 = ipos2 + 1

      end do

      j = j - 1

      end do

! now shift bits to the left by max_levels - level

      nbitshft = ndim*(lrefine_max-lrefine)
      if(nbitshft.gt.0) then

        nshifts = 0
        do while (nbitshft > nbits)
           nbitshft = nbitshft - nbits
           nshifts = nshifts + 1
        end do

        do k = 1,nshifts
        mort(1) = ishft(mort(1),nbits)
        do j = 1,5
        ipos2 = 0
        do ipos = 0,nbits-1
           call mvbits (mort(j+1),ipos,1,mort(j),ipos2)
           mort(j+1) = ibclr(mort(j+1),ipos)
           ipos2 = ipos2 + 1
        end do
        mort(j+1) = ishft(mort(j+1),nbits)
        end do
        end do

        mort(1) = ishft(mort(1),nbitshft)
        do j = 1,5
        ipos2 = 0
        do ipos = nbits-1-(nbitshft-1),nbits-1
           call mvbits (mort(j+1),ipos,1,mort(j),ipos2)
           mort(j+1) = ibclr(mort(j+1),ipos)
           ipos2 = ipos2 + 1
        end do
        mort(j+1) = ishft(mort(j+1),nbitshft)
        end do

      end if

      return
      end subroutine morton_number


      subroutine boundary_locator(x0,y0,z0,lboundary,ibc)

!------------------------------------------------------------------------
!
! This routine detects whether a given coordinate is on a particular
! boundary.
! The boundary regions are assumed to be volumes defined by their 
! bounding boxes, and each has a boundary condition index associated
! with it.
!
! Written :     Peter MacNeice         June 2000
!------------------------------------------------------------------------
!
! Arguments:
!
!------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      logical,intent(out) :: lboundary
      integer,intent(out) :: ibc
      real,intent(in)     :: x0,y0,z0

      real    :: xtest, ytest, ztest
      integer :: ibnd 
  

!------------------------------------
 
      lboundary = .false.
      do ibnd = 1, nboundaries

        xtest = (x0 - boundary_box(1,1,ibnd))* & 
     &          (x0 - boundary_box(2,1,ibnd))
        ytest = (y0 - boundary_box(1,2,ibnd))* & 
     &          (y0 - boundary_box(2,2,ibnd))
        ztest = (z0 - boundary_box(1,3,ibnd))* & 
     &          (z0 - boundary_box(2,3,ibnd))
        if(ndim.lt.3) ztest = -1.
        if(ndim.lt.2) ytest = -1.
        if(xtest.le.0. .and. ytest.le.0. .and. ztest.le.0.) then
          lboundary = .true.
          ibc       = boundary_index(ibnd)
        endif

      enddo

      return
      end subroutine boundary_locator



      subroutine mpi_xchange_blocks(mype,nprocs, tag_offset, & 
     &                              buf_dim_send, S_buffer, &
     &                              buf_dim_recv, R_buffer)

!------------------------------------------------------------------------
!
! This routine uses the general hand-shaking scheme for message passing 
! between pe's. 
! sends and receives buffers containing all requested block data.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!      buf_dim_send        dimension of send buffers
!      buf_dim_recv        dimension of recv buffers
!      tag_offset     offset for MPI tag
!      S_buffer       send buffer
!      R_buffer       recv buffer
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton
      Use paramesh_comm_data
 
      implicit none

      include 'mpif.h'

      integer, intent(in)    :: mype,nprocs,buf_dim_send,buf_dim_recv
      integer, intent(inout) :: tag_offset
      real,    intent(inout) :: S_buffer(buf_dim_send)
      real,    intent(inout) :: R_buffer(buf_dim_recv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: i,j
      integer :: istrt,ilast,isize
      integer :: isrc,idest, itag
      integer :: recvrequest(nprocs)
      integer :: recvstatus(MPI_STATUS_SIZE,nprocs)
      integer :: sendvrequest(nprocs)
      integer :: ierrorcode,ierr
      integer :: ij, ji


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!      R_buffer = 0

! Removed by CEG
!      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if (commatrix_send(mype+1) > 0  & 
     &    .or.  & 
     &    commatrix_recv(mype+1) > 0) then
        write(*,*) 'Paramesh error :  error in xchange : pe ',mype, & 
     &      ' diagonal element of commatrix is non-zero ', & 
     &      commatrix_recv(mype+1), & 
     &      commatrix_send(mype+1)
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

      ij = 0
      ji = 0

      do i = 1,nprocs
        isrc = i-1
        idest= mype
        itag = isrc*nprocs + idest+1 + tag_offset
                                 ! receive to pe=j
        if(commatrix_recv(i).gt.0) then
          ij = ij+1
          istrt = ir_buf(1,i)
          ilast = ir_buf(2,i)
          isize = ilast-istrt+1

          if (amr_error_checking) then
            if(isize.gt.buf_dim_recv) then
              write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 1:', & 
     &                    ' message is bigger than buffer space  ' & 
     &           ,' isrc ',isrc,' idest ',idest,' i ',i, & 
     &           ' istrt ',istrt,' ilast ',ilast,' isize ',isize, & 
     &           ' buf_dim_recv ',buf_dim_recv
                 call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
            endif
          endif
          call Mpi_Irecv(R_buffer(istrt),isize, & 
     &                     amr_mpi_real, & 
     &                     isrc ,itag,MPI_COMM_WORLD, & 
     &                     recvrequest(ij),ierr)
        endif
      enddo

!  print *,commatrix_send(mype), commatrix_recv(mype)

! CEG - This was two distinct but identical loops - first being mype+1->noprocs, second being 1->mype

#ifdef OLDWAY
      do j = mype+1,nprocs
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
				! send from mype=i
         if(commatrix_send(j).gt.0) then
            ji = ji+1
            istrt = is_buf(1,j)
            ilast = is_buf(2,j)
            isize = ilast-istrt+1

            if (amr_error_checking) then
             if(isize.gt.buf_dim_send) then
                 write(*,*) 'PARAMESH ERROR : ', & 
     &                   ' message is bigger than buffer space  '
                 call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
             endif
             endif
! CEG changed from Ssend to Isend
             call MPI_Isend(S_buffer(istrt),isize, & 
!             call MPI_Ssend(S_buffer(istrt),isize, & 
     &                      amr_mpi_real, & 
     &                      idest,itag,MPI_COMM_WORLD, &
     &                      sendvrequest(ji), &
     &                      ierr)
        endif
      enddo

      do j = 1,mype
#endif
!CEG switched last line to this
      do j = 1,nprocs
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
				! send from mype=i
         if(commatrix_send(j).gt.0) then
            ji = ji+1
            istrt = is_buf(1,j)
            ilast = is_buf(2,j)
            isize = ilast-istrt+1

            if (amr_error_checking) then
               if(isize.gt.buf_dim_send) then
                  write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 2:', & 
     &                   ' message is bigger than buffer space  '
                  call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
               endif
            endif
! CEG changed from Ssend to Isend
            call MPI_Isend(S_buffer(istrt),isize, & 
!            call MPI_Ssend(S_buffer(istrt),isize, & 
     &                      amr_mpi_real, & 
     &                      idest,itag,MPI_COMM_WORLD, &
     &                      sendvrequest(ji),  &
     &                      ierr)
        endif
      enddo

! CEG Ssend / Isend hanging - could it be I had the wait for receives before the wait for sends?

      if(ij.gt.0) & 
     &   call MPI_Waitall(ij,recvrequest,recvstatus, & 
     &                    ierrorcode)
      if(ji.gt.0) & 
     &   call MPI_Waitall(ji,sendvrequest, MPI_STATUSES_IGNORE, & 
     &                    ierrorcode)

! reset offset with largest tag
      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset          

      return
      end subroutine mpi_xchange_blocks



      subroutine mpi_xchange_tree_info(mype,nprocs, tag_offset, & 
     &                              buf_dim, S_buffer, R_buffer)

!------------------------------------------------------------------------
!
! This routine uses the general hand-shaking scheme for message passing 
! between pe's. 
! sends and receives buffers containing all requested block data.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!      buf_dim        dimension of buffers
!      tag_offset     offset for MPI tag
!      S_buffer       send buffer
!      R_buffer       recv buffer
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton
      Use paramesh_comm_data

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype,nprocs,buf_dim
      integer, intent(inout) ::  tag_offset
      real,    intent(inout) :: S_buffer(buf_dim), R_buffer(buf_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: i,j, ij, ji
      integer :: istrt,ilast,isize
      integer :: isrc,idest, itag
      integer :: ierr,errcode
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:),  allocatable :: sendvrequest
      integer,dimension (:,:),allocatable :: recvstatus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if(allocated(recvrequest)) deallocate( recvrequest )
      allocate ( recvrequest(nprocs) )
      if(allocated(recvstatus)) deallocate( recvstatus )
      allocate ( recvstatus(MPI_STATUS_SIZE,nprocs) )

!!      R_buffer = 0

      ij = 0
      ji = 0

      do i = 1,nprocs
         isrc = i-1
         idest= mype
         itag = isrc*nprocs + idest+1 + tag_offset
                                 ! receive to pe=j
         if(commatrix_recv(i).gt.0) then
            ij = ij+1
            istrt = ir_buf(1,i)
            ilast = ir_buf(2,i)
            isize = ilast-istrt+1

            if (amr_error_checking) then
             if(isize.gt.buf_dim) then
               write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 3:', & 
     &                   ' message is bigger than buffer space  '
                 call mpi_abort(MPI_COMM_WORLD,errcode,ierr)
             endif
             end if

             call Mpi_Irecv(R_buffer(istrt),isize,amr_mpi_real, & 
     &           isrc ,itag,MPI_COMM_WORLD,recvrequest(ij),ierr)
        endif
      enddo

      do j = 1,nprocs
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
				! send from mype=i
         if(commatrix_send(j).gt.0) then
            ji = ji+1
            istrt = is_buf(1,j)
            ilast = is_buf(2,j)
            isize = ilast-istrt+1

            if (amr_error_checking) then
             if(isize.gt.buf_dim) then
               write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 4:', & 
     &                   ' message is bigger than buffer space  '
                 call mpi_abort(MPI_COMM_WORLD,errcode,ierr)
             endif
             end if

! CEG changed to Isend from Ssend
             call MPI_Isend(S_buffer(istrt),isize,amr_mpi_real, & 
!             call MPI_Ssend(S_buffer(istrt),isize,amr_mpi_real, & 
     &           idest,itag,MPI_COMM_WORLD, &
     &           sendvrequest(ji), ierr)
!     &           ierr)
        endif
      enddo

! reset offset with largest tag
      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset          


      if(ij.gt.0) & 
     &   call MPI_Waitall(ij,recvrequest,recvstatus, & 
     &                    ierr)
      deallocate( recvrequest )
      deallocate( recvstatus )

      if(ji.gt.0) & 
     &   call MPI_Waitall(ji,sendvrequest,MPI_STATUSES_IGNORE, & 
     &                    ierr)

      deallocate( sendvrequest )

      return
      end subroutine mpi_xchange_tree_info

!--------------------------------------------------------------------

