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

      subroutine mpi_get_flux_buffer(mype,lb,dtype,offset, & 
     &                          buffer_size,S_buffer,flux_dir)

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
!      offset         offset for buffer index
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

      integer, intent(in)    :: lb,mype,buffer_size
      integer, intent(inout) :: offset
      real,    intent(inout) :: S_buffer(buffer_size)
      integer, optional, intent(in) :: flux_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0 
      integer :: nguard_work0 

      integer :: index,ierrorcode,ierr
      integer :: ia,ib,ja,jb,ka,kb
      integer :: ia0,ib0,ja0,jb0,ka0,kb0
      integer :: ilimit(1)
      integer :: vtype, flux_dirt
      integer :: i, j, k, n, ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      ilimit = shape(S_buffer)

      if(lb.gt.maxblocks_alloc) then
        write(*,*) 'ERROR : mpi_get_buffer pe ',mype, & 
     &        ' putting blk ',lb,' into Sbuf'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      end if


! take over incoming index offset for block lb to be sent to remote pe

      index = offset          
      

! store block location as first info in buffer
      S_buffer(index) = lb
      S_buffer(index+1) = mype

      S_buffer(index+2) = dtype
      index = index + 3


! set vtype as though appropriate for cell-centered data
      vtype = 1
      call mpi_set_message_limits(dtype, & 
     &                            ia0,ib0,ja0,jb0,ka0,kb0,vtype)

! pack the flux_x array for block lb

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if (flux_dirt == 1 .or. flux_dirt == 0) then

      if(dtype.eq.13.or.dtype.eq.15.or.dtype.eq.14) then

      if(dtype.eq.13) then
        ia = 1
        ib = 1
      elseif(dtype.eq.15) then
        ia = 2
        ib = 2
      elseif(dtype.eq.14) then
        ia = 1
        ib = 2
      endif
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nfluxes
          S_buffer(index) = flux_x(n,i,j,k,lb)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif

      end if

! pack the flux_y array for block lb

      if(ndim.ge.2) then
      if (flux_dirt == 2 .or. flux_dirt == 0) then
      if(dtype.eq.11.or.dtype.eq.17.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if(dtype.eq.11) then
        ja = 1
        jb = 1
      elseif(dtype.eq.17) then
        ja = 2
        jb = 2
      elseif(dtype.eq.14) then
        ja = 1
        jb = 2
      endif
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nfluxes
          S_buffer(index) = flux_y(n,i,j,k,lb)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo


      endif
      endif
      endif


! pack the flux_z array for block lb

      if(ndim.eq.3) then
      if (flux_dirt == 3 .or. flux_dirt == 0) then
      if(dtype.eq.5.or.dtype.eq.23.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if(dtype.eq.5) then
        ka = 1
        kb = 1
      elseif(dtype.eq.23) then
        ka = 2
        kb = 2
      elseif(dtype.eq.14) then
        ka = 1
        kb = 2
      endif
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nfluxes
          S_buffer(index) = flux_z(n,i,j,k,lb)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif
      endif

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
      end subroutine mpi_get_flux_buffer

!------------------------------------------------------------------------

      subroutine mpi_get_Sbuffer_size_fluxes(mype,lb,dtype,offset, & 
     &                                       flux_dir)

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
!      offset         offset for buffer index
!      buffer_size    size of send buffer
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

      integer, intent(in)    :: lb,mype
      integer, intent(inout) :: offset
      integer, optional, intent(in) :: flux_dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: index,ierrorcode,ierr
      integer :: ia,ib,ja,jb,ka,kb
      integer :: ia0,ib0,ja0,jb0,ka0,kb0
      integer :: ilimit(1)
      integer :: vtype, flux_dirt
      integer :: i, j, k, n, ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      end if

! take over incoming index offset for block lb to be sent to remote pe

      index = offset          
      
! store block location as first info in buffer

      index = index + 3

! set vtype as though appropriate for cell-centered data
      vtype = 1
      call mpi_set_message_limits(dtype, & 
     &                            ia0,ib0,ja0,jb0,ka0,kb0,vtype)

! pack the flux_x array for block lb

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if (flux_dirt == 1 .or. flux_dirt == 0) then

      if(dtype.eq.13.or.dtype.eq.15.or.dtype.eq.14) then

      if(dtype.eq.13) then
        ia = 1
        ib = 1
      elseif(dtype.eq.15) then
        ia = 2
        ib = 2
      elseif(dtype.eq.14) then
        ia = 1
        ib = 2
      endif
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif

      end if

! pack the flux_y array for block lb

      if(ndim.ge.2) then
      if (flux_dirt == 2 .or. flux_dirt == 0) then
      if(dtype.eq.11.or.dtype.eq.17.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if(dtype.eq.11) then
        ja = 1
        jb = 1
      elseif(dtype.eq.17) then
        ja = 2
        jb = 2
      elseif(dtype.eq.14) then
        ja = 1
        jb = 2
      endif
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo


      endif
      endif
      endif


! pack the flux_z array for block lb

      if(ndim.eq.3) then
      if (flux_dirt == 3 .or. flux_dirt == 0) then
      if(dtype.eq.5.or.dtype.eq.23.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if(dtype.eq.5) then
        ka = 1
        kb = 1
      elseif(dtype.eq.23) then
        ka = 2
        kb = 2
      elseif(dtype.eq.14) then
        ka = 1
        kb = 2
      endif
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif
      endif

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
      end subroutine mpi_get_Sbuffer_size_fluxes
