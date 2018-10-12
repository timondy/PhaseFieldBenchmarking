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


      subroutine mpi_amr_write_guard_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr
      logical,save        :: first = .true.

      if(.not.allocated(pe_source_guard_mg))  & 
     &   allocate(pe_source_guard_mg(1:nprocs,1:lrefine_max), & 
     &            stat = istat) 


!!!

      if (.not.allocated(to_be_levels))  & 
     &     allocate(to_be_levels(1:lrefine_max))

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_received_guard_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_received_guard_mg)
      endif
      if(allocated(to_be_received)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &       to_be_received_guard_mg(3,size(to_be_received,2), & 
     &                               size(to_be_received,3)), & 
     &                               stat = istat)
      else
        write(*,*)'Paramesh error: write_guard_mg ', & 
     &            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_received_guard_mg(:,:,:) =  & 
     &       to_be_received(:,:,:)

!!!

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_sent_guard_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_sent_guard_mg)
       first = .false.
      end if
      if(allocated(to_be_sent)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &         to_be_sent_guard_mg(3,size(to_be_sent,2), & 
     &                             size(to_be_sent,3)), & 
     &                             stat = istat) 
      else
        write(*,*)'Paramesh error: write_guard ', & 
     &            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_sent_guard_mg(:,:,:) =  & 
     &        to_be_sent(:,:,:)

!!!

      if(.not.allocated(commatrix_guard_mg))  & 
     &      allocate(commatrix_guard_mg(1:nprocs,2,1:lrefine_max), & 
     &               stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_write_guard_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      pe_source_guard_mg(:,level) = pe_source(:)
      commatrix_guard_mg(:,:,level) = 0

      do i = 1,nprocs
         commatrix_guard_mg(i,1,level) = commatrix_recv(i)
         commatrix_guard_mg(i,2,level) = commatrix_send(i)
      enddo

      if (.not.allocated(strt_guard_mg)) then
         allocate(strt_guard_mg(1:lrefine_max))
         allocate(max_no_to_send_guard_mg(1:lrefine_max))
      end if
      strt_guard_mg(level) = strt_buffer
      max_no_to_send_guard_mg(level) = max_no_to_send

      if (.not.allocated(laddress_guard_mg)) then
        allocate(laddress_guard_mg(1:2,1:maxblocks_alloc,1:lrefine_max))
      end if
      laddress_guard_mg(:,:,level) = laddress(:,:)

      return
      end subroutine mpi_amr_write_guard_comm_mg




      subroutine mpi_amr_read_guard_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr

      pe_source_guard(:) = pe_source_guard_mg(:,level)
      commatrix_guard(:,:) = 0

      do i = 1,nprocs
        commatrix_guard(i,1) = commatrix_guard_mg(i,1,level)
        commatrix_guard(i,2) = commatrix_guard_mg(i,2,level)
      enddo

      strt_guard = strt_guard_mg(level)
      max_no_to_send_guard = max_no_to_send_guard_mg(level)

      laddress_guard(:,:) = laddress_guard_mg(:,:,level)


      if(max_no_to_send_guard.gt.0) then

      if(allocated(to_be_received_guard))  & 
     &                 deallocate(to_be_received_guard)
      allocate(to_be_received_guard(3,size(to_be_levels(level) % & 
     &                               to_be_received_guard_mg(:,:,:),2), & 
     &                          size(to_be_levels(level) %  & 
     &                               to_be_received_guard_mg(:,:,:),3)), & 
     &                        stat = istat) 
      to_be_received_guard(:,:,:) = to_be_levels(level) %  & 
     &                        to_be_received_guard_mg(:,:,:)

      if(allocated(to_be_sent_guard)) deallocate(to_be_sent_guard)
      allocate(                                                        & 
        to_be_sent_guard(3,                                            &
              size(to_be_levels(level)%to_be_sent_guard_mg(:,:,:),2),  & 
              size(to_be_levels(level)%to_be_sent_guard_mg(:,:,:),3)), & 
              stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_read_guard_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      to_be_sent_guard(1:3,:,:)= & 
       to_be_levels(level)%to_be_sent_guard_mg(1:3,:,:) 

      endif

      mpi_pattern_id = 10

      call amr_1blk_guardcell_reset

      return
      end subroutine mpi_amr_read_guard_comm_mg




!!!!!!!!!!!!!!
! Prolongation

      subroutine mpi_amr_write_prol_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr
      logical, save       :: first = .true.

      if(.not.allocated(pe_source_prol_mg))  & 
     &   allocate(pe_source_prol_mg(1:nprocs,1:lrefine_max), & 
     &            stat = istat) 


!!!

      if (.not.allocated(to_be_levels)) & 
     &     allocate(to_be_levels(1:lrefine_max))

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_received_prol_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_received_prol_mg)
      end if
      if(allocated(to_be_received)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &       to_be_received_prol_mg(3,size(to_be_received,2), & 
     &                               size(to_be_received,3)), & 
     &                               stat = istat)
      else
        write(*,*)'Paramesh error: write_prol_mg ', & 
     &            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_received_prol_mg(:,:,:) =  & 
     &       to_be_received(:,:,:)

!!!

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_sent_prol_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_sent_prol_mg)
       first = .false.
      end if
      if(allocated(to_be_sent)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &         to_be_sent_prol_mg(3,size(to_be_sent,2), & 
     &                             size(to_be_sent,3)), & 
     &                             stat = istat) 
      else
        write(*,*)'Paramesh error: write_prol ', & 
     &            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_sent_prol_mg(:,:,:) =  & 
     &        to_be_sent(:,:,:)

!!!

      if(.not.allocated(commatrix_prol_mg))  & 
     &      allocate(commatrix_prol_mg(1:nprocs,2,1:lrefine_max), & 
     &               stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_write_prol_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      pe_source_prol_mg(:,level) = pe_source(:)
      commatrix_prol_mg(:,:,level) = 0

      do i = 1,nprocs
         commatrix_prol_mg(i,1,level) = commatrix_recv(i)
         commatrix_prol_mg(i,2,level) = commatrix_send(i)
      enddo

      if (.not.allocated(strt_prol_mg)) then
         allocate(strt_prol_mg(1:lrefine_max))
         allocate(max_no_to_send_prol_mg(1:lrefine_max))
      end if
      strt_prol_mg(level) = strt_buffer
      max_no_to_send_prol_mg(level) = max_no_to_send

      if (.not.allocated(laddress_prol_mg)) then
        allocate(laddress_prol_mg(1:2,1:maxblocks_alloc,1:lrefine_max))
      end if
      laddress_prol_mg(:,:,level) = laddress(:,:)

      return
      end subroutine mpi_amr_write_prol_comm_mg




      subroutine mpi_amr_read_prol_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr

      pe_source_prol(:) = pe_source_prol_mg(:,level)
      commatrix_prol(:,:) = 0

      do i = 1,nprocs
        commatrix_prol(i,1) = commatrix_prol_mg(i,1,level)
        commatrix_prol(i,2) = commatrix_prol_mg(i,2,level)
      enddo

      strt_prol = strt_prol_mg(level)
      max_no_to_send_prol = max_no_to_send_prol_mg(level)

      laddress_prol(:,:) = laddress_prol_mg(:,:,level)


      if(max_no_to_send_prol.gt.0) then

      if(allocated(to_be_received_prol))  & 
     &                 deallocate(to_be_received_prol)
      allocate(to_be_received_prol(3,size(to_be_levels(level) % & 
     &                               to_be_received_prol_mg(:,:,:),2), & 
     &                          size(to_be_levels(level) %  & 
     &                               to_be_received_prol_mg(:,:,:),3)), & 
     &                        stat = istat) 
      to_be_received_prol(:,:,:) = to_be_levels(level) %  & 
     &                        to_be_received_prol_mg(:,:,:)

      if(allocated(to_be_sent_prol)) deallocate(to_be_sent_prol)
      allocate(                                                        & 
       to_be_sent_prol(3,                                              &
               size(to_be_levels(level)%to_be_sent_prol_mg(:,:,:),2),  &
               size(to_be_levels(level)%to_be_sent_prol_mg(:,:,:),3)), & 
               stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_read_prol_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif


      to_be_sent_prol(1:3,:,:) =                                       &
       to_be_levels(level)%to_be_sent_prol_mg(1:3,:,:)

      endif

      mpi_pattern_id = 20

      return
      end subroutine mpi_amr_read_prol_comm_mg






!!!!!!!!!!!!!!!!!!!
! Flux Conservation

      subroutine mpi_amr_write_flux_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr
      logical, save       :: first = .true.

      if(.not.allocated(pe_source_flux_mg))  & 
     &   allocate(pe_source_flux_mg(1:nprocs,1:lrefine_max), & 
     &            stat = istat) 


!!!

      if (.not.allocated(to_be_levels)) & 
     &     allocate(to_be_levels(1:lrefine_max))

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_received_flux_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_received_flux_mg)
      end if
      if(allocated(to_be_received)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &       to_be_received_flux_mg(3,size(to_be_received,2), & 
     &                               size(to_be_received,3)), & 
     &                               stat = istat)
      else
        write(*,*)'Paramesh error: write_flux_mg ', & 
     &            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_received_flux_mg(:,:,:) =  & 
     &       to_be_received(:,:,:)

!!!

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_sent_flux_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_sent_flux_mg)
       first = .false.
      end if
      if(allocated(to_be_sent)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &         to_be_sent_flux_mg(3,size(to_be_sent,2), & 
     &                             size(to_be_sent,3)), & 
     &                             stat = istat) 
      else
        write(*,*)'Paramesh error: write_flux ', & 
     &            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_sent_flux_mg(:,:,:) =  & 
     &        to_be_sent(:,:,:)

!!!

      if(.not.allocated(commatrix_flux_mg))  & 
     &      allocate(commatrix_flux_mg(1:nprocs,2,1:lrefine_max), & 
     &               stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_write_flux_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      pe_source_flux_mg(:,level) = pe_source(:)
      commatrix_flux_mg(:,:,level) = 0

      do i = 1,nprocs
         commatrix_flux_mg(i,1,level) = commatrix_recv(i)
         commatrix_flux_mg(i,2,level) = commatrix_send(i)
      enddo

      if (.not.allocated(strt_flux_mg)) then
         allocate(strt_flux_mg(1:lrefine_max))
         allocate(max_no_to_send_flux_mg(1:lrefine_max))
      end if
      strt_flux_mg(level) = strt_buffer
      max_no_to_send_flux_mg(level) = max_no_to_send

      if (.not.allocated(laddress_flux_mg)) then
        allocate(laddress_flux_mg(1:2,1:maxblocks_alloc,1:lrefine_max))
      end if
      laddress_flux_mg(:,:,level) = laddress(:,:)

      return
      end subroutine mpi_amr_write_flux_comm_mg




      subroutine mpi_amr_read_flux_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr

      pe_source_flux(:) = pe_source_flux_mg(:,level)
      commatrix_flux(:,:) = 0

      do i = 1,nprocs
        commatrix_flux(i,1) = commatrix_flux_mg(i,1,level)
        commatrix_flux(i,2) = commatrix_flux_mg(i,2,level)
      enddo

      strt_flux = strt_flux_mg(level)
      max_no_to_send_flux = max_no_to_send_flux_mg(level)

      laddress_flux(:,:) = laddress_flux_mg(:,:,level)


      if(max_no_to_send_flux.gt.0) then

      if(allocated(to_be_received_flux))  & 
     &                 deallocate(to_be_received_flux)
      allocate(to_be_received_flux(3,size(to_be_levels(level) % & 
     &                               to_be_received_flux_mg(:,:,:),2), & 
     &                          size(to_be_levels(level) %  & 
     &                               to_be_received_flux_mg(:,:,:),3)), & 
     &                        stat = istat) 
      to_be_received_flux(:,:,:) = to_be_levels(level) %  & 
     &                        to_be_received_flux_mg(:,:,:)

      if(allocated(to_be_sent_flux)) deallocate(to_be_sent_flux)
      allocate(                                                        & 
        to_be_sent_flux(3,                                             &
               size(to_be_levels(level)%to_be_sent_flux_mg(:,:,:),2),  & 
               size(to_be_levels(level)%to_be_sent_flux_mg(:,:,:),3)), & 
               stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_read_flux_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

       to_be_sent_flux(:,:,:) =                                        &
         to_be_levels(level)%to_be_sent_flux_mg(:,:,:)

      endif

      mpi_pattern_id = 30

      return
      end subroutine mpi_amr_read_flux_comm_mg





!!!!!!!!!!!!!
! Restriction


      subroutine mpi_amr_write_restrict_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr
      logical, save       :: first = .true.

      if(.not.allocated(pe_source_restrict_mg))  & 
     &   allocate(pe_source_restrict_mg(1:nprocs,1:lrefine_max), & 
     &            stat = istat) 


!!!

      if (.not.allocated(to_be_levels)) & 
     &     allocate(to_be_levels(1:lrefine_max))

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_received_restrict_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_received_restrict_mg)
      end if
      if(allocated(to_be_received)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &       to_be_received_restrict_mg(3,size(to_be_received,2), & 
     &                               size(to_be_received,3)), & 
     &                               stat = istat)
      else
        write(*,*)'Paramesh error: write_restrict_mg ', & 
     &            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_received_restrict_mg(:,:,:) =  & 
     &       to_be_received(:,:,:)

!!!

      if (.not.first) then
       if(associated(to_be_levels(level) % to_be_sent_restrict_mg))  & 
     &  deallocate(to_be_levels(level) % to_be_sent_restrict_mg)
       first = .false.
      end if
      if(allocated(to_be_sent)) then
        allocate( & 
     &       to_be_levels(level) %  & 
     &         to_be_sent_restrict_mg(3,size(to_be_sent,2), & 
     &                             size(to_be_sent,3)), & 
     &                             stat = istat) 
      else
        write(*,*)'Paramesh error: write_restrict ', & 
     &            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_levels(level) % to_be_sent_restrict_mg(:,:,:) =  & 
     &        to_be_sent(:,:,:)

!!!

      if(.not.allocated(commatrix_restrict_mg))  & 
     &      allocate(commatrix_restrict_mg(1:nprocs,2,1:lrefine_max), & 
     &               stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_write_restrict_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      pe_source_restrict_mg(:,level) = pe_source(:)
      commatrix_restrict_mg(:,:,level) = 0

      do i = 1,nprocs
         commatrix_restrict_mg(i,1,level) = commatrix_recv(i)
         commatrix_restrict_mg(i,2,level) = commatrix_send(i)
      enddo

      if (.not.allocated(strt_restrict_mg)) then
         allocate(strt_restrict_mg(1:lrefine_max))
         allocate(max_no_to_send_restrict_mg(1:lrefine_max))
      end if
      strt_restrict_mg(level) = strt_buffer
      max_no_to_send_restrict_mg(level) = max_no_to_send

      if (.not.allocated(laddress_restrict_mg)) then
       allocate & 
     & (laddress_restrict_mg(1:2,1:maxblocks_alloc,1:lrefine_max))
      end if
      laddress_restrict_mg(:,:,level) = laddress(:,:)

      return
      end subroutine mpi_amr_write_restrict_comm_mg




      subroutine mpi_amr_read_restrict_comm_mg(nprocs,level)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use amr_mg_common

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs, level

      integer             :: i,istat,ierrorcode,ierr

      pe_source_restrict(:) = pe_source_restrict_mg(:,level)
      commatrix_restrict(:,:) = 0

      do i = 1,nprocs
        commatrix_restrict(i,1) = commatrix_restrict_mg(i,1,level)
        commatrix_restrict(i,2) = commatrix_restrict_mg(i,2,level)
      enddo

      strt_restrict = strt_restrict_mg(level)
      max_no_to_send_restrict = max_no_to_send_restrict_mg(level)

      laddress_restrict(:,:) = laddress_restrict_mg(:,:,level)


      if(max_no_to_send_restrict.gt.0) then

      if(allocated(to_be_received_restrict))  & 
     &                 deallocate(to_be_received_restrict)
      allocate(to_be_received_restrict(3,size(to_be_levels(level) % & 
     &                             to_be_received_restrict_mg(:,:,:),2), & 
     &                          size(to_be_levels(level) %  & 
     &                            to_be_received_restrict_mg(:,:,:),3)), & 
     &                        stat = istat) 
      to_be_received_restrict(:,:,:) = to_be_levels(level) %  & 
     &                        to_be_received_restrict_mg(:,:,:)

      if(allocated(to_be_sent_restrict)) deallocate(to_be_sent_restrict)
      allocate(                                                        & 
        to_be_sent_restrict(3,                                         &
           size(to_be_levels(level)%to_be_sent_restrict_mg(:,:,:),2),  &
           size(to_be_levels(level)%to_be_sent_restrict_mg(:,:,:),3)), &
           stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_read_restrict_comm_mg: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      to_be_sent_restrict(1:3,:,:) =                                   &
       to_be_levels(level)%to_be_sent_restrict_mg(1:3,:,:)

      endif

      mpi_pattern_id = 40

      return
      end subroutine mpi_amr_read_restrict_comm_mg

