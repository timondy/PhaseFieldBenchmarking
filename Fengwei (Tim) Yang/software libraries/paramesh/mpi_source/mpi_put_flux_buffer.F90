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


      subroutine mpi_put_flux_buffer(mype,lb,offset, & 
     &                          buffer_size,R_buffer,flux_dir)

!------------------------------------------------------------------------
!
! This subroutine unpacks the buffer R_buffer which has been received on 
! this processor into the flux arrays for the local block lb.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      lb             guard block id to be unpacked from receiving
!      offset         offset for buffer index
!      buffer_size    size of receive buffer
!      R_buffer       receive buffer
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'


      integer, intent(in)    :: mype,lb,buffer_size
      integer, intent(inout) :: offset
      real,    intent(inout) :: R_buffer(buffer_size)
      integer, optional, intent(in) :: flux_dir


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0 
      integer :: nguard_work0 
      integer :: index
      integer :: ia,ib,ja,jb,ka,kb
      integer :: ia0,ib0,ja0,jb0,ka0,kb0
      integer :: dtype,vtype,flux_dirt
      integer :: i, j, k, n, ii, jj, iii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

! take over incoming index offset for block lb to be sent to remote pe

      index = offset          


      dtype = anint(R_buffer(index+2))
      index = index+3

      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      end if

      vtype = 1
      call mpi_set_message_limits(dtype, & 
     &                            ia0,ib0,ja0,jb0,ka0,kb0,vtype)

! unpack the flux_x array for block lb
      if (flux_dirt == 1 .or. flux_dirt == 0) then
      if(dtype.eq.13.or.dtype.eq.15.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

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

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif

! pack the facevary array for block lb

      if(ndim.ge.2) then
      if(flux_dirt == 2 .or. flux_dirt == 0) then
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

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif
      endif

! pack the facevarz array for block lb

      if(ndim.eq.3) then
      if(flux_dirt == 3 .or. flux_dirt == 0) then
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


      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif
      endif

! unpack tree info

      do ii = 1,mdim
         coord(ii,lb) = R_buffer(index)
         index = index + 1
      end do

      do ii = 1,mdim
         bsize(ii,lb) = R_buffer(index)
         index = index + 1
      end do

      do ii = 1,2
         do jj = 1,mdim
            bnd_box(ii,jj,lb) = R_buffer(index)
            index = index + 1
         end do
      end do

      do ii = 1,2
         parent(ii,lb) = int(R_buffer(index))
         index = index + 1
      end do

      do ii = 1,2
         do jj = 1,mchild
            child(ii,jj,lb) = R_buffer(index)
            index = index + 1
         end do
      end do

      which_child(lb) = int(R_buffer(index))
      index = index + 1

      if (int(R_buffer(index)) == 1) then
         newchild(lb) = .true.
      else
         newchild(lb) = .false.
      end if
      index = index + 1

      do ii = 1,2
         do jj = 1,mfaces
            neigh(ii,jj,lb) = int(R_buffer(index))
            index = index + 1
         end do
      end do

      lrefine(lb) = int(R_buffer(index))
      index = index + 1

      nodetype(lb) = int(R_buffer(index))
      index = index + 1

      do ii = 1,3
         do k = 1,1+2*k3d
            do jj = 1,1+2*k2d
               do iii = 1,3
                  surr_blks(ii,iii,jj,k,lb) = & 
     &                 int(R_buffer(index))
                  index = index + 1
               end do
            end do
         end do
      end do


#ifdef SAVE_MORTS
      do ii = 1,6
         do k = 1,1+2*k3d
            do jj = 1,1+2*k2d
               do iii = 1,3
                  surr_morts(ii,iii,jj,k,lb) = & 
     &                 int(R_buffer(index))
                  index = index + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

      empty(lb) = int(R_buffer(index))
      index = index + 1

      if (int(R_buffer(index)) == 1) then
         ldtcomplete(lb) = .true.
      else
         ldtcomplete(lb) = .false.
      end if
      index = index + 1

! overwrite outgoing index

      offset = index

      return
      end subroutine mpi_put_flux_buffer

!------------------------------------------------------------------------

      subroutine mpi_get_Rbuffer_size_fluxes(lb,dtype,offset, & 
     &                                       flux_dir)

!------------------------------------------------------------------------
!
! This subroutine unpacks the buffer R_buffer which has been received on 
! this processor into the flux arrays for the local block lb.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      lb             guard block id to be unpacked from receiving
!      offset         offset for buffer index
!      buffer_size    size of receive buffer
!      R_buffer       receive buffer
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'


      integer, intent(in)    :: lb, dtype
      integer, intent(inout) :: offset
      integer, optional, intent(in) :: flux_dir


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: index
      integer :: ia,ib,ja,jb,ka,kb
      integer :: ia0,ib0,ja0,jb0,ka0,kb0
      integer :: vtype,flux_dirt
      integer :: i, j, k, n, ii, jj, iii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! take over incoming index offset for block lb to be sent to remote pe

      index = offset          

      index = index+3

      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      end if

      vtype = 1
      call mpi_set_message_limits(dtype, & 
     &                            ia0,ib0,ja0,jb0,ka0,kb0,vtype)

! unpack the flux_x array for block lb
      if (flux_dirt == 1 .or. flux_dirt == 0) then
      if(dtype.eq.13.or.dtype.eq.15.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

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

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif

! pack the facevary array for block lb

      if(ndim.ge.2) then
      if(flux_dirt == 2 .or. flux_dirt == 0) then
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

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif
      endif

! pack the facevarz array for block lb

      if(ndim.eq.3) then
      if(flux_dirt == 3 .or. flux_dirt == 0) then
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


      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nfluxes
          index  = index + 1
        enddo
      enddo
      enddo
      enddo

      endif
      endif
      endif

! unpack tree info

      do ii = 1,mdim
         index = index + 1
      end do

      do ii = 1,mdim
         index = index + 1
      end do

      do ii = 1,2
         do jj = 1,mdim
            index = index + 1
         end do
      end do

      do ii = 1,2
         index = index + 1
      end do

      do ii = 1,2
         do jj = 1,mchild
            index = index + 1
         end do
      end do

      index = index + 1

      index = index + 1

      do ii = 1,2
         do jj = 1,mfaces
            index = index + 1
         end do
      end do

      index = index + 1

      index = index + 1

      do ii = 1,3
         do k = 1,1+2*k3d
            do jj = 1,1+2*k2d
               do iii = 1,3
                  index = index + 1
               end do
            end do
         end do
      end do


#ifdef SAVE_MORTS
      do ii = 1,6
         do k = 1,1+2*k3d
            do jj = 1,1+2*k2d
               do iii = 1,3
                  index = index + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

      index = index + 1

      index = index + 1

! overwrite outgoing index

      offset = index

      return
      end subroutine mpi_get_Rbuffer_size_fluxes
