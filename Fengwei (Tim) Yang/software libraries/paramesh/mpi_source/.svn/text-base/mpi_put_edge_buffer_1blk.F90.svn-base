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


      subroutine mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)

!------------------------------------------------------------------------
!
! This subroutine unpacks the buffer R_buffer which has been received on 
! this processor into the bedge_faceX_X arrays for the local block lb.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      lb             guard block id to be unpacked from receiving
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton

      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: lb,remote_pe,remote_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0 
      integer :: nguard_work0
      integer :: index, index0
      integer :: ia,ib,ja,jb,ka,kb
      integer :: ia0,ib0,ja0,jb0,ka0,kb0
      integer :: dtype,vtype
      integer :: i, j, k, n, ii, jj, iii
      integer :: mype, ierr
      logical :: lfound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)

! take over incoming index offset for block lb to be sent to remote pe
      
      call amr_mpi_find_blk_in_buffer(mype,remote_block, & 
     &     remote_pe,1,dtype,index0,lfound)
      vtype = 8
      call mpi_set_message_limits(dtype, & 
     &                            ia0,ib0,ja0,jb0,ka0,kb0,vtype)

      index = index0 + 1

! unpack the bedge_facex_y and bedge_facex_z arrays for block lb
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
        do n=1,nedges
          recvarx1e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.13) then
        recvarz2e(1:nedges,1+nguard0,ja:jb,1) = & 
     &      recvarx1e(1:nedges,1,ja:jb,1+nguard0*k3d)
        recvarz2e(1:nedges,1+nguard0,ja:jb,1+k3d) = & 
     &      recvarx1e(1:nedges,1,ja:jb,nzb+(1+nguard0)*k3d)
      elseif(dtype.eq.15) then
        recvarz2e(1:nedges,nxb+1+nguard0,ja:jb,1) = & 
     &      recvarx1e(1:nedges,2,ja:jb,1+nguard0*k3d)
        recvarz2e(1:nedges,nxb+1+nguard0,ja:jb,1+k3d) = & 
     &      recvarx1e(1:nedges,2,ja:jb,nzb+(1+nguard0)*k3d)
      endif

      if(ndim.eq.3.or.l2p5d.eq.1) then
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nedges
          recvarx2e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.13) then
        recvary2e(1:nedges,1+nguard0,1,ka:kb) = & 
     &      recvarx2e(1:nedges,1,1+nguard0*k2d,ka:kb)
        recvary2e(1:nedges,1+nguard0,2,ka:kb) = & 
     &      recvarx2e(1:nedges,1,nyb+(1+nguard0)*k2d,ka:kb)
      elseif(dtype.eq.15) then
        recvary2e(1:nedges,nxb+1+nguard0,1,ka:kb) = & 
     &      recvarx2e(1:nedges,2,1+nguard0*k2d,ka:kb)
        recvary2e(1:nedges,nxb+1+nguard0,2,ka:kb) = & 
     &      recvarx2e(1:nedges,2,nyb+(1+nguard0)*k2d,ka:kb)
      endif

      endif

      endif


! pack the bedge_facey_x and bedge_facey_z arrays for block lb

      if(ndim.ge.2) then
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
        do n=1,nedges
          recvary1e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.11) then
        recvarz1e(1:nedges,ia:ib,1+nguard0*k2d,1) = & 
     &      recvary1e(1:nedges,ia:ib,1,1+nguard0*k3d)
        recvarz1e(1:nedges,ia:ib,1+nguard0*k2d,1+k3d) = & 
     &      recvary1e(1:nedges,ia:ib,1,nzb+(1+nguard0)*k3d)
      elseif(dtype.eq.17) then
        recvarz1e(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,1) = & 
     &      recvary1e(1:nedges,ia:ib,2,1+nguard0*k3d)
        recvarz1e(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,1+k3d) = & 
     &      recvary1e(1:nedges,ia:ib,2,nzb+(1+nguard0)*k3d)
      endif

      if(ndim.eq.3.or.l2p5d.eq.1) then
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nedges
          recvary2e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.11) then
         recvarx2e(1:nedges,1,1+nguard0*k2d,ka:kb) = & 
     &        recvary2e(1:nedges,1+nguard0,1,ka:kb)
         recvarx2e(1:nedges,2,1+nguard0*k2d,ka:kb) = & 
     &        recvary2e(1:nedges,nxb+1+nguard0,1,ka:kb)
      elseif(dtype.eq.17) then
         recvarx2e(1:nedges,1,nyb+(1+nguard0)*k2d,ka:kb) = & 
     &        recvary2e(1:nedges,1+nguard0,2,ka:kb)
         recvarx2e(1:nedges,2,nyb+(1+nguard0)*k2d,ka:kb) = & 
     &        recvary2e(1:nedges,nxb+1+nguard0,2,ka:kb)
      endif

      endif


      endif
      endif

! pack the bedge_facez_x array for block lb

      if(ndim.eq.3) then
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
        kb = 1+k3d
      endif


      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nedges
          recvarz1e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.5) then
        recvary1e(1:nedges,ia:ib,1,1+nguard0*k3d) = & 
     &        recvarz1e(1:nedges,ia:ib,1+nguard0*k2d,1)
        recvary1e(1:nedges,ia:ib,2,1+nguard0*k3d) = & 
     &        recvarz1e(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,1)
      elseif(dtype.eq.23) then
        recvary1e(1:nedges,ia:ib,1,nzb+(1+nguard0)*k3d) = & 
     &        recvarz1e(1:nedges,ia:ib,1+nguard0*k2d,2)
        recvary1e(1:nedges,ia:ib,2,nzb+(1+nguard0)*k3d) = & 
     &        recvarz1e(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,2)
      endif

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nedges
          recvarz2e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.5) then
        recvarx1e(1:nedges,1,ja:jb,1+nguard0*k3d) = & 
     &        recvarz2e(1:nedges,1+nguard0,ja:jb,1)
        recvarx1e(1:nedges,2,ja:jb,1+nguard0*k3d) = & 
     &        recvarz2e(1:nedges,nxb+1+nguard0,ja:jb,1)
      elseif(dtype.eq.23) then
        recvarx1e(1:nedges,1,ja:jb,nzb+(1+nguard0)*k3d) = & 
     &        recvarz2e(1:nedges,1+nguard0,ja:jb,2)
        recvarx1e(1:nedges,2,ja:jb,nzb+(1+nguard0)*k3d) = & 
     &        recvarz2e(1:nedges,nxb+1+nguard0,ja:jb,2)
      endif

      endif
      endif


! Now unpack single edges


! first x edges
      if(ndim.eq.3) then
      if(dtype.eq.2.or.dtype.eq.8.or.dtype.eq.20.or.dtype.eq.26) then

      ia = ia0
      ib = ib0
      j = 1+nguard0
      k = 1

      if(dtype.eq.8) then
        j = nyb+k2d*(1+nguard0)
      elseif(dtype.eq.20) then
        k = 2
      elseif(dtype.eq.26) then
        k = 2
        j = nyb+k2d*(1+nguard0)
      endif
      do i = ia , ib
        do n=1,nedges
          recvarz1e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo

      ia = ia0
      ib = ib0
      j = 1
      k = 1+nguard0*k3d

      if(dtype.eq.8) then
        j = 2
      elseif(dtype.eq.20) then
        k = nzb+k3d*(1+nguard0)
      elseif(dtype.eq.26) then
        k = nzb+k3d*(1+nguard0)
        j = 2
      endif
      do i = ia , ib
        do n=1,nedges
          recvary1e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo


      endif
      endif


! now y edges

      if(ndim.eq.3) then
      if(dtype.eq.4.or.dtype.eq.6.or.dtype.eq.22.or.dtype.eq.24) then

      ja = ja0
      jb = jb0
      i = 1
      k = 1+nguard0*k3d

      if(dtype.eq.6) then
        i = 2
      elseif(dtype.eq.22) then
        k = nzb+k3d*(1+nguard0)
      elseif(dtype.eq.24) then
        k = nzb+k3d*(1+nguard0)
        i = 2
      endif
      do j = ja , jb
        do n=1,nedges
          recvarx1e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo

      ja = ja0
      jb = jb0
      i = 1+nguard0
      k = 1

      if(dtype.eq.6) then
        i = nxb+1+nguard0
      elseif(dtype.eq.22) then
        k = 2
      elseif(dtype.eq.24) then
        i = nxb+1+nguard0
        k = 2
      endif
      do j = ja , jb
        do n=1,nedges
          recvarz2e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo


      endif
      endif


! finally z edges

      if(ndim.eq.3.or.l2p5d.eq.1) then
      if(dtype.eq.10.or.dtype.eq.12.or.dtype.eq.16.or.dtype.eq.18) then

      ka = ka0
      kb = kb0
      i = 1
      j = 1+nguard0*k2d

      if(dtype.eq.12) then
        i = 2
      elseif(dtype.eq.16) then
        j = nyb+k2d*(1+nguard0)
      elseif(dtype.eq.18) then
        i = 2
        j = nyb+k2d*(1+nguard0)
      endif
      do k = ka , kb
        do n=1,nedges
          recvarx2e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo

      ka = ka0
      kb = kb0
      i = 1+nguard0
      j = 1

      if(dtype.eq.12) then
        i = nxb+1+nguard0
      elseif(dtype.eq.16) then
        j = 2
      elseif(dtype.eq.18) then
        i = nxb+1+nguard0
        j = 2
      endif
      do k = ka , kb
        do n=1,nedges
          recvary2e(n,i,j,k) = temprecv_buf(index)
          index  = index + 1
        enddo
      enddo


      endif
      endif


      return
      end subroutine mpi_put_edge_buffer_1blk
