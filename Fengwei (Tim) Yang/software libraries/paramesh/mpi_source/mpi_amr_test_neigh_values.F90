!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_test_neigh_values


! This routine tests the values in the neigh array for consistency.
! If the neighbors across a given face do not agree for a particular
! pair of blocks it reports an error and its location.

! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use constants
      Use paramesh_comm_data

      use paramesh_interfaces

      implicit none

      include 'mpif.h'

! local amr variables
      integer :: nprocs,mype
      integer :: ierror
      integer :: status1(MPI_STATUS_SIZE)
      logical :: lflag

!---------------------------------------------------------------
!
        integer,dimension (:),  allocatable :: glnblocks
        integer,dimension (:,:,:,:),  allocatable :: neight
        real,dimension (:,:,:,:),  allocatable :: bndboxt


        integer :: lb,iproc,k,isize,nrecv,lnblocks_max
        integer :: jface,jr,rem_blk,rem_pe,nerrors
        real    :: eps,accuracy

!---------------------------------------------------------------
        Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
        Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)

        accuracy = 100./10.**precision(accuracy)
        eps = pi*accuracy

#ifdef DEBUG
        if(mype.eq.0) write(*,*) 'Started amr_test_neigh_values'
        do lb=1,lnblocks
          write(*,*) 'blk/pe ',lb,mype,' neigh ',neigh(:,:,lb)
        enddo
#endif /* DEBUG */

! Initialize some message buffers


! broadcast number of blocks on each processor
        if(allocated(glnblocks)) deallocate( glnblocks )
        allocate ( glnblocks(0:nprocs-1) )
! CEG modofied
!        glnblocks(mype) = lnblocks
        call MPI_ALLGATHER(lnblocks, 1,MPI_INTEGER, & 
     &                   glnblocks,1,MPI_INTEGER, & 
     &                   MPI_COMM_WORLD,ierror)
#ifdef DEBUG
        write(*,*) 'amr_test_neigh_values : pe ',mype, & 
     &             ' glnblocks ',glnblocks
#endif /* DEBUG */

        lnblocks_max = maxval(glnblocks)

        if(mype.eq.0) then
          if(allocated(neight)) deallocate(neight)
          allocate(neight(2,6,lnblocks_max,0:nprocs-1))
          if(allocated(bndboxt)) deallocate(bndboxt)
          allocate(bndboxt(2,3,lnblocks_max,0:nprocs-1))

          neight = -1
          neight(:,:,:lnblocks,0) = neigh(:,:,:lnblocks)
          bndboxt = -1.0
          bndboxt(:,:,:lnblocks,0) = bnd_box(:,:,:lnblocks)
        endif


        do iproc=1,nprocs-1

        
!------------
! CEG removed
!        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        if(mype.eq.0) then

! non-blocking receive of neigh data from processor iproc.
          nrecv  = 0
          isize  = 2*6*glnblocks(iproc)
          if(isize.gt.0) then
            nrecv = 1
            call MPI_int_RECV(neight(1,1,1,iproc), & 
     &                       isize, & 
     &                       MPI_INTEGER, & 
     &                       iproc, & 
     &                       iproc, & 
     &                       MPI_COMM_WORLD, & 
     &                       ierror)

          endif
        endif

        if(mype.eq.iproc) then
          isize = 2*6*glnblocks(iproc)
! blocking send of neigh data from processor iproc.
          if(isize.gt.0) & 
     &      call MPI_int_SSEND(neigh(1,1,1), & 
     &                    isize, & 
     &                    MPI_INTEGER, & 
     &                    0, & 
     &                    iproc, & 
     &                    MPI_COMM_WORLD, & 
     &                    ierror)

        endif

#ifdef DEBUG
        if(mype.eq.0) & 
     &   write(*,*) 'a amr_test_neigh_values : pe ',mype, & 
     &             ' nrecv ',nrecv
#endif /* DEBUG */
!------------


        if(spherical_pm) then
! CEG removed
!        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        if(mype.eq.0) then

! non-blocking receive of bnd_box data from processor iproc.
          nrecv  = 0
          isize  = 2*3*glnblocks(iproc)
          if(isize.gt.0) then
            nrecv = 1
            call MPI_RECV(bndboxt(1,1,1,iproc), & 
     &                       isize, & 
     &                       amr_mpi_real, & 
     &                       iproc, & 
     &                       iproc, & 
     &                       MPI_COMM_WORLD, & 
     &                       ierror)

          endif
        endif

        if(mype.eq.iproc) then
          isize = 2*3*glnblocks(iproc)
! blocking send of neigh data from processor iproc.
          if(isize.gt.0) & 
     &      call MPI_real_SSEND(bnd_box(1,1,1), & 
     &                    isize, & 
     &                    amr_mpi_real, & 
     &                    0, & 
     &                    iproc, & 
     &                    MPI_COMM_WORLD, & 
     &                    ierror)

        endif

#ifdef DEBUG
        if(mype.eq.0) & 
     &   write(*,*) 'b amr_test_neigh_values : pe ',mype, & 
     &             ' nrecv ',nrecv
#endif /* DEBUG */
!------------
        endif





        enddo
! CEG removed
!        Call MPI_BARRIER(MPI_COMM_WORLD, ierror)

!
! Test neigh values
        if(mype.eq.0) then
#ifdef DEBUG
        do iproc=0,nprocs-1
        do lb=1,glnblocks(iproc)
          write(*,*) 'blk/pe ',lb,iproc,' neight ',neight(:,:,lb,iproc)
        enddo
        enddo
#endif /* DEBUG */
        nerrors = 0
        do iproc=0,nprocs-1
          do lb=1,glnblocks(iproc)
            do jface=1,2*ndim
              rem_blk = neight(1,jface,lb,iproc)
              rem_pe  = neight(2,jface,lb,iproc)
              if(jface.eq.1) jr = 2
              if(jface.eq.2) jr = 1
              if(jface.eq.3) jr = 4
              if(jface.eq.4) jr = 3
              if(jface.eq.5) jr = 6
              if(jface.eq.6) jr = 5

         if(jface.eq.3.and.abs(bndboxt(1,2,lb,iproc)).lt.eps)    jr=3
         if(jface.eq.4.and.abs(bndboxt(2,2,lb,iproc)-pi).lt.eps) jr=4

              if(rem_blk.gt.0) then
              if(neight(1,jr,rem_blk,rem_pe).gt.0) then
              if(neight(1,jr,rem_blk,rem_pe).ne.lb.or. & 
     &           neight(2,jr,rem_blk,rem_pe).ne.iproc) then

                write(*,*) 'Neighbor inconsistency : block ',lb,iproc, & 
     &                     ' face ',jface,' has neighbor ', & 
     &                     neight(:,jface,lb,iproc), & 
     &                     ' but that neighbor thinks its neighbor is ', & 
     &                     neight(:,jr,rem_blk,rem_pe)
                nerrors = nerrors + 1
              endif
              endif
              endif

            enddo
          enddo
        enddo
        write(*,*) 'A total of ',nerrors,' neigh inconsistencies found.'
        endif
        Call MPI_BARRIER(MPI_COMM_WORLD, ierror)


! deallocate temporary arrays
        if(mype.eq.0.and.allocated(neight)) deallocate(neight)
        if(mype.eq.0.and.allocated(bndboxt)) deallocate(bndboxt)


        return
      end subroutine amr_test_neigh_values
