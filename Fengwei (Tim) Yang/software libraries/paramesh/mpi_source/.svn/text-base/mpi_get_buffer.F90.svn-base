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

      subroutine mpi_get_buffer(mype,lb,dtype,iopt,offset, & 
     &                          lcc,lfc,lec,lnc, & 
     &                          buffer_size,S_buffer, & 
     &                          nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine packs the block lb of the local arrays work,unk,facevar
! into the buffer S_buffer which is to be sent from mype.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor id 
!      lb             local block id to be packed for sending
!      iopt           option setting for unk or work array
!      offset         offset for buffer index
!      lcc            logical switch controlling whether unk or work data
!                     is packed
!      lfc            logical switch controlling whether facevar data
!                     is packed
!      lec            logical switch controlling whether unk_e_? data
!                     is packed
!      lnc            logical switch controlling whether unk_n data
!                     is packed
!      buffer_size    size of send buffer
!      S_buffer       send buffer
!
! new code
!      dtype          type of message to be added to buffer
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: dtype

      integer, intent(in)    :: lb,mype,iopt,buffer_size
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      real,    intent(inout) :: S_buffer(buffer_size)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0
      integer :: nguard_work0

      integer :: index,ierrorcode,ierr
      integer :: ilimit(1)
      integer :: i, j, k
      integer :: ia, ib, ja, jb, ka, kb
      integer :: n, ii
      integer :: vtype
      integer :: invar,ivar,ivar_next

#ifdef DEBUG
      logical :: lprint
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      ilimit = shape(S_buffer)

      if(lb.gt.maxblocks_alloc) then
        write(*,*) 'ERROR : mpi_get_buffer pe ',mype, & 
     &        ' putting blk ',lb,' into Sbuf'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif


! take over incoming index offset for block lb to be sent to remote pe

      index = offset          


! store block location as first info in buffer
      S_buffer(index) = lb
      S_buffer(index+1) = mype

#ifdef DEBUG
      lprint = .false.
      if(lcc) lprint = .true.
      if(lfc) lprint = .true.
      if(lec) lprint = .true.
      if(lnc) lprint = .true.
      if(iopt.gt.1) lprint = .true.
      if(.not.lprint) & 
     & write(*,*) 'pe ',mype,' in get_buf with no flags on! - ', & 
     &    ' while packing block ',lb, & 
     &    ' message type ',dtype
#endif /* DEBUG */


! new code
      S_buffer(index+2) = dtype
      index = index + 3


      if(iopt.gt.1) then


! pack the work array for block lb

        vtype = 0
        call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
        do k = ka , kb
          do j = ja , jb
            do i = ia , ib
              S_buffer(index) = work(i,j,k,lb,iopt-1)
              index  = index + 1
            enddo
          enddo
        enddo

      else                     ! iopt iftest


! pack the unk array for block lb

        invar = nvar
        if (lcc.and.lguard_in_progress)  & 
     &         invar = ngcell_on_cc
        if (lcc.and.invar.gt.0) then


          vtype = 1
          call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
          do k = ka , kb
            do j = ja , jb
              do i = ia , ib
                do ivar = 1,invar
                  ivar_next = gcell_on_cc_pointer(ivar)
                  if (no_permanent_guardcells) then
                    S_buffer(index+ivar-1) = gt_unk(ivar_next,i,j,k,lb)
                  else
                    S_buffer(index+ivar-1) = unk(ivar_next,i,j,k,lb)
                  end if
                enddo
              index = index + invar
            enddo
          enddo
        enddo

      endif                    ! lcc


! logical switch controling whether facevar data is being packed

      invar = nfacevar
      if (lfc.and.lguard_in_progress)  & 
     &         invar = maxval(ngcell_on_fc(1:ndim))
      if (lfc.and.invar.gt.0) then

! pack the facevarx array for block lb


      vtype = 2
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_fc_pointer(1,ivar)
          if (no_permanent_guardcells) then
          S_buffer(index+ivar-1) = gt_facevarx(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = facevarx(ivar_next,i,j,k,lb)
          end if
        enddo
        index  = index + invar
      enddo
      enddo
      enddo



! pack the facevary array for block lb

        if(ndim.ge.2) then


      vtype = 3
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_fc_pointer(2,ivar)
          if (no_permanent_guardcells) then
          S_buffer(index+ivar-1) = gt_facevary(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = facevary(ivar_next,i,j,k,lb)
          end if
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


        endif


! pack the facevarz array for block lb

        if(ndim.eq.3.or.l2p5d.eq.1) then


      vtype = 4
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_fc_pointer(3,ivar)
          if (no_permanent_guardcells) then
          S_buffer(index+ivar-1) = gt_facevarz(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = facevarz(ivar_next,i,j,k,lb)
          end if
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


        endif


      endif                    ! lfc

      invar = nvaredge
      if (lec.and.lguard_in_progress)  & 
     &         invar = maxval(ngcell_on_ec(1:3))
      if (lec .and. ndim.gt.1.and.invar.gt.0) then

! pack the unk_e_x array for block lb

      vtype = 5
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_ec_pointer(1,ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_e_x(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_e_x(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


      if (ndim > 1) then
! pack the unk_e_y array for block lb


      vtype = 6
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_ec_pointer(2,ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_e_y(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_e_y(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo

      if ((ndim == 3).or.(l2p5d == 1)) then
! pack the unk_e_z array for block lb


      vtype = 7
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_ec_pointer(3,ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_e_z(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_e_z(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo

      end if
      end if

      endif                    ! lec


! pack the unk_n array for block lb
      invar = nvarcorn
      if (lnc.and.lguard_in_progress)  & 
     &         invar = ngcell_on_nc
      if (lnc .and.invar.gt.0) then


      vtype = 8
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)

      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_nc_pointer(ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_n(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_n(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


      endif                    ! lnc

      endif                    ! end of iopt iftest

! Add tree info to buffer

      do i = 1,mdim
         S_buffer(index) = coord(i,lb)
         index = index + 1
      end do

      do i = 1,mdim
         S_buffer(index) = bsize(i,lb)
         index = index + 1
      end do

      do i = 1,2
         do j = 1,mdim
            S_buffer(index) = bnd_box(i,j,lb)
            index = index + 1
         end do
      end do

      do i = 1,2
         S_buffer(index) = real(parent(i,lb))
         index = index + 1
      end do
 
      do i = 1,2
         do j = 1,mchild
            S_buffer(index) = child(i,j,lb)
            index = index + 1
         end do
      end do

      S_buffer(index) = real(which_child(lb))
      index = index + 1

      if (newchild(lb)) then
         S_buffer(index) = real(1)
      else
         S_buffer(index) = real(0)
      end if
      index = index + 1

      do i = 1,2
         do j = 1,mfaces
            S_buffer(index) = real(neigh(i,j,lb))
            index = index + 1
         end do
      end do

      S_buffer(index) = real(lrefine(lb))
      index = index + 1

      S_buffer(index) = real(nodetype(lb))
      index = index + 1

      do ii = 1,3
         do k = 1,1+2*k3d
            do j = 1,1+2*k2d
               do i = 1,3
                  S_buffer(index) = & 
     &                 real(surr_blks(ii,i,j,k,lb))
                  index = index + 1
               end do
            end do
         end do
      end do

#ifdef SAVE_MORTS
      do ii = 1,6
         do k = 1,1+2*k3d
            do j = 1,1+2*k2d
               do i = 1,3
                  S_buffer(index) = & 
     &                 real(surr_morts(ii,i,j,k,lb))
                  index = index + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

      S_buffer(index) = real(empty(lb))
      index = index + 1

      if (ldtcomplete(lb)) then
         S_buffer(index) = real(1)
      else
         S_buffer(index) = real(0)
      end if
      index = index + 1

! overwrite outgoing index offset

      offset = index 

      return
      end subroutine mpi_get_buffer

!------------------------------------------------------------------------

      subroutine mpi_get_Sbuffer_size(mype,lb,dtype,iopt,offset, & 
     &                                lcc,lfc,lec,lnc, & 
     &                                nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine computes the buffer size for messages from a single proc.
!
!
! Written :     Kevin Olson January 2007
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor id 
!      lb             local block id to be packed for sending
!      iopt           option setting for unk or work array
!      offset         offset for buffer index
!      lcc            logical switch controlling whether unk or work data
!                     is packed
!      lfc            logical switch controlling whether facevar data
!                     is packed
!      lec            logical switch controlling whether unk_e_? data
!                     is packed
!      lnc            logical switch controlling whether unk_n data
!                     is packed
!
! new code
!      dtype          type of message to be added to buffer
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: dtype

      integer, intent(in)    :: lb,mype,iopt
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0
      integer :: nguard_work0

      integer :: index,ierrorcode,ierr
      integer :: ilimit(1)
      integer :: i, j, k
      integer :: ia, ib, ja, jb, ka, kb
      integer :: n, ii
      integer :: vtype
      integer :: invar,ivar,ivar_next

#ifdef DEBUG
      logical :: lprint
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      if(lb.gt.maxblocks_alloc) then
        write(*,*) 'ERROR : mpi_get_buffer pe ',mype, & 
     &        ' putting blk ',lb,' into Sbuf'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif


! take over incoming index offset for block lb to be sent to remote pe

      index = offset          

#ifdef DEBUG
      lprint = .false.
      if(lcc) lprint = .true.
      if(lfc) lprint = .true.
      if(lec) lprint = .true.
      if(lnc) lprint = .true.
      if(iopt.gt.1) lprint = .true.
      if(.not.lprint) & 
     & write(*,*) 'pe ',mype,' in get_buf with no flags on! - ', & 
     &    ' while packing block ',lb, & 
     &    ' message type ',dtype
#endif /* DEBUG */


      index = index + 3

      if(iopt.gt.1) then


! pack the work array for block lb

      vtype = 0
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + 1
      enddo
      enddo
      enddo

      else                     ! iopt iftest


! pack the unk array for block lb

      invar = nvar
      if (lcc.and.lguard_in_progress)  & 
     &         invar = ngcell_on_cc
      if (lcc.and.invar.gt.0) then


      vtype = 1
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index = index + invar
      enddo
      enddo
      enddo

      endif                    ! lcc


! logical switch controling whether facevar data is being packed

      invar = nfacevar
      if (lfc.and.lguard_in_progress)  & 
     &         invar = maxval(ngcell_on_fc(1:ndim))
      if (lfc.and.invar.gt.0) then

! pack the facevarx array for block lb


      vtype = 2
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + invar
      enddo
      enddo
      enddo



! pack the facevary array for block lb

        if(ndim.ge.2) then


      vtype = 3
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + invar
      enddo
      enddo
      enddo


        endif


! pack the facevarz array for block lb

        if(ndim.eq.3.or.l2p5d.eq.1) then


      vtype = 4
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + invar
      enddo
      enddo
      enddo


        endif


      endif                    ! lfc

      invar = nvaredge
      if (lec.and.lguard_in_progress)  & 
     &         invar = maxval(ngcell_on_ec(1:3))
      if (lec .and. ndim.gt.1.and.invar.gt.0) then

! pack the unk_e_x array for block lb

      vtype = 5
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + invar
      enddo
      enddo
      enddo


      if (ndim > 1) then
! pack the unk_e_y array for block lb


      vtype = 6
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + invar
      enddo
      enddo
      enddo

      if ((ndim == 3).or.(l2p5d == 1)) then
! pack the unk_e_z array for block lb


      vtype = 7
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + invar
      enddo
      enddo
      enddo

      end if
      end if

      endif                    ! lec


! pack the unk_n array for block lb
      invar = nvarcorn
      if (lnc.and.lguard_in_progress)  & 
     &         invar = ngcell_on_nc
      if (lnc .and.invar.gt.0) then


      vtype = 8
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)

      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        index  = index + invar
      enddo
      enddo
      enddo


      endif                    ! lnc

      endif                    ! end of iopt iftest

! Add tree info to buffer

      do i = 1,mdim
         index = index + 1
      end do

      do i = 1,mdim
         index = index + 1
      end do

      do i = 1,2
         do j = 1,mdim
            index = index + 1
         end do
      end do

      do i = 1,2
         index = index + 1
      end do
 
      do i = 1,2
         do j = 1,mchild
            index = index + 1
         end do
      end do

      index = index + 1

      index = index + 1

      do i = 1,2
         do j = 1,mfaces
            index = index + 1
         end do
      end do

      index = index + 1

      index = index + 1

      do ii = 1,3
         do k = 1,1+2*k3d
            do j = 1,1+2*k2d
               do i = 1,3
                  index = index + 1
               end do
            end do
         end do
      end do

#ifdef SAVE_MORTS
      do ii = 1,6
         do k = 1,1+2*k3d
            do j = 1,1+2*k2d
               do i = 1,3
                  index = index + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

      index = index + 1

      index = index + 1

! overwrite outgoing index offset

      offset = index 

      return
      end subroutine mpi_get_Sbuffer_size




!-------------------------------------------------------------------------------------

      subroutine pf_mpi_get_buffer(mype,lb,dtype,iopt,offset, & 
     &                          lcc,lfc,lec,lnc, & 
     &                          buffer_size,S_buffer, & 
     &                          nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine packs the block lb of the local arrays work,unk,facevar
! into the buffer S_buffer which is to be sent from mype.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor id 
!      lb             local block id to be packed for sending
!      iopt           option setting for unk or work array
!      offset         offset for buffer index
!      lcc            logical switch controlling whether unk or work data
!                     is packed
!      lfc            logical switch controlling whether facevar data
!                     is packed
!      lec            logical switch controlling whether unk_e_? data
!                     is packed
!      lnc            logical switch controlling whether unk_n data
!                     is packed
!      buffer_size    size of send buffer
!      S_buffer       send buffer
!
! new code
!      dtype          type of message to be added to buffer
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: dtype

      integer, intent(in)    :: lb,mype,iopt,buffer_size
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      real,    intent(inout) :: S_buffer(buffer_size)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0
      integer :: nguard_work0

      integer :: index,ierrorcode,ierr
      integer :: ilimit(1)
      integer :: i, j, k
      integer :: ia, ib, ja, jb, ka, kb
      integer :: n, ii
      integer :: vtype
      integer :: invar,ivar,ivar_next

#ifdef DEBUG
      logical :: lprint
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      ilimit = shape(S_buffer)

      if(lb.gt.maxblocks_alloc) then
        write(*,*) 'ERROR : mpi_get_buffer pe ',mype, & 
     &        ' putting blk ',lb,' into Sbuf'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif


! take over incoming index offset for block lb to be sent to remote pe

      index = offset          


! store block location as first info in buffer
      S_buffer(index) = lb
      S_buffer(index+1) = mype

#ifdef DEBUG
      lprint = .false.
      if(lcc) lprint = .true.
      if(lfc) lprint = .true.
      if(lec) lprint = .true.
      if(lnc) lprint = .true.
      if(iopt.gt.1) lprint = .true.
      if(.not.lprint) & 
     & write(*,*) 'pe ',mype,' in get_buf with no flags on! - ', & 
     &    ' while packing block ',lb, & 
     &    ' message type ',dtype
#endif /* DEBUG */


! new code
      S_buffer(index+2) = dtype
      index = index + 3


      if(iopt.gt.1) then


! pack the work array for block lb

        vtype = 0
        call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
        do k = ka , kb
          do j = ja , jb
            do i = ia , ib
              S_buffer(index) = work(i,j,k,lb,iopt-1)
              index  = index + 1
            enddo
          enddo
        enddo

      else                     ! iopt iftest


! pack the unk array for block lb

        invar = nvar
        if (lcc.and.lguard_in_progress)  & 
     &         invar = ngcell_on_cc
        if (lcc.and.invar.gt.0) then


          vtype = 1
          call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
          do k = ka , kb
            do j = ja , jb
              do i = ia , ib
                do ivar = 1,invar
                  ivar_next = gcell_on_cc_pointer(ivar)
                  if (no_permanent_guardcells) then
                    S_buffer(index+ivar-1) = gt_unk(ivar_next,i,j,k,lb)
                  else
                    S_buffer(index+ivar-1) = unk(ivar_next,i,j,k,lb)
                  end if
                enddo
              index = index + invar
            enddo
          enddo
        enddo

      endif                    ! lcc


! logical switch controling whether facevar data is being packed

      invar = nfacevar
      if (lfc.and.lguard_in_progress)  & 
     &         invar = maxval(ngcell_on_fc(1:ndim))
      if (lfc.and.invar.gt.0) then

! pack the facevarx array for block lb


      vtype = 2
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_fc_pointer(1,ivar)
          if (no_permanent_guardcells) then
          S_buffer(index+ivar-1) = gt_facevarx(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = facevarx(ivar_next,i,j,k,lb)
          end if
        enddo
        index  = index + invar
      enddo
      enddo
      enddo



! pack the facevary array for block lb

        if(ndim.ge.2) then


      vtype = 3
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_fc_pointer(2,ivar)
          if (no_permanent_guardcells) then
          S_buffer(index+ivar-1) = gt_facevary(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = facevary(ivar_next,i,j,k,lb)
          end if
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


        endif


! pack the facevarz array for block lb

        if(ndim.eq.3.or.l2p5d.eq.1) then


      vtype = 4
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_fc_pointer(3,ivar)
          if (no_permanent_guardcells) then
          S_buffer(index+ivar-1) = gt_facevarz(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = facevarz(ivar_next,i,j,k,lb)
          end if
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


        endif


      endif                    ! lfc

      invar = nvaredge
      if (lec.and.lguard_in_progress)  & 
     &         invar = maxval(ngcell_on_ec(1:3))
      if (lec .and. ndim.gt.1.and.invar.gt.0) then

! pack the unk_e_x array for block lb

      vtype = 5
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_ec_pointer(1,ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_e_x(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_e_x(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


      if (ndim > 1) then
! pack the unk_e_y array for block lb


      vtype = 6
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_ec_pointer(2,ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_e_y(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_e_y(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo

      if ((ndim == 3).or.(l2p5d == 1)) then
! pack the unk_e_z array for block lb


      vtype = 7
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_ec_pointer(3,ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_e_z(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_e_z(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo

      end if
      end if

      endif                    ! lec


! pack the unk_n array for block lb
      invar = nvarcorn
      if (lnc.and.lguard_in_progress)  & 
     &         invar = ngcell_on_nc
      if (lnc .and.invar.gt.0) then


      vtype = 8
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)

      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do ivar = 1,invar
          ivar_next = gcell_on_nc_pointer(ivar)
          if (no_permanent_guardcells) then
        if(index.gt.ilimit(1)) then
           write(*,*) 'Error : pe ',mype,' index ',index, & 
     & ' ilimit ',ilimit,'i,j,k,lb ', i,j,k,lb & 
     & ,' dtype ',dtype,' vtype ',vtype
           call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
          S_buffer(index+ivar-1) = gt_unk_n(ivar_next,i,j,k,lb)
          else
          S_buffer(index+ivar-1) = unk_n(ivar_next,i,j,k,lb)
          endif
        enddo
        index  = index + invar
      enddo
      enddo
      enddo


      endif                    ! lnc

      endif                    ! end of iopt iftest

! Add tree info to buffer

      do i = 1,mdim
         S_buffer(index) = coord(i,lb)
         index = index + 1
      end do

      do i = 1,mdim
         S_buffer(index) = bsize(i,lb)
         index = index + 1
      end do

      do i = 1,2
         do j = 1,mdim
            S_buffer(index) = bnd_box(i,j,lb)
            index = index + 1
         end do
      end do

      do i = 1,2
         S_buffer(index) = real(parent(i,lb))
         index = index + 1
      end do
 
      do i = 1,2
         do j = 1,mchild
            S_buffer(index) = child(i,j,lb)
            index = index + 1
         end do
      end do

      S_buffer(index) = real(which_child(lb))
      index = index + 1

      if (newchild(lb)) then
         S_buffer(index) = real(1)
      else
         S_buffer(index) = real(0)
      end if
      index = index + 1

      do i = 1,2
         do j = 1,mfaces
            S_buffer(index) = real(neigh(i,j,lb))
            index = index + 1
         end do
      end do

      S_buffer(index) = real(lrefine(lb))
      index = index + 1

      S_buffer(index) = real(nodetype(lb))
      index = index + 1

      do ii = 1,3
         do k = 1,1+2*k3d
            do j = 1,1+2*k2d
               do i = 1,3
                  S_buffer(index) = & 
     &                 real(surr_blks(ii,i,j,k,lb))
                  index = index + 1
               end do
            end do
         end do
      end do

#ifdef SAVE_MORTS
      do ii = 1,6
         do k = 1,1+2*k3d
            do j = 1,1+2*k2d
               do i = 1,3
                  S_buffer(index) = & 
     &                 real(surr_morts(ii,i,j,k,lb))
                  index = index + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

      S_buffer(index) = real(empty(lb))
      index = index + 1

! CEG removed this as this is for fluxes, unset in PF code
!      if (ldtcomplete(lb)) then
!         S_buffer(index) = real(1)
!      else
         S_buffer(index) = real(0)
!      end if
      index = index + 1

! overwrite outgoing index offset

      offset = index 

      return
      end subroutine pf_mpi_get_buffer


