!#define DEBUG
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


      subroutine mpi_put_buffer(lb,ioptw,offset, & 
     &                          lcc,lfc,lec,lnc, & 
     &                          buffer_size,R_buffer, & 
     &                          nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine unpacks the buffer R_buffer which has been received on 
! this processor into the local guard block lb of the arrays work,unk,
! and facevar.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      lb             guard block id to be unpacked from receiving
!      ioptw          option setting for work array
!      offset         offset for buffer index
!      lcc            logical switch controlling whether unk or work data
!                     is unpacked
!      lfc            logical switch controlling whether facevar data
!                     is unpacked
!      lec            logical switch controlling whether unk_e_? or work data
!                     is unpacked
!      lnc            logical switch controlling whether unk_n data
!                     is unpacked
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

      integer, intent(in)    :: lb,ioptw,buffer_size
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      real,    intent(inout) :: R_buffer(buffer_size)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0 
      integer :: nguard_work0 
      integer :: index
      integer :: dtype,vtype
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ii, jj, k, iii, invar
      integer :: ngcell_on_fc0,ngcell_on_ec0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DEBUG

      integer :: mype, ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
#endif /* DEBUG */


      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

! take over incoming index offset for block lb to be sent to remote pe

      index = offset          

       dtype = anint(R_buffer(index+2))
       index = index+3

! pack the work array for block lb
      if(ioptw.gt.1) then

      vtype = 0
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + (kb-ka+1)*(jb-ja+1)*(ib-ia+1)


      else                           ! ioptw



#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index on entry ',index
#endif /* DEBUG */


! logical switch controls whether unk data are being packed
      invar = nvar
!      if (lcc.and.lguard_in_progress.and.ngcell_on_cc.gt.0)
      if (lcc.and.lguard_in_progress) & 
     &         invar = ngcell_on_cc
      if (lcc) then

! pack the unk array for block lb

      vtype = 1
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for cc ',index, & 
     &         ' invar ',invar,' ia ib ja jb ka kb ',ia,ib,ja,jb,ka,kb, & 
     &         ' dtype ',dtype
#endif /* DEBUG */

      endif                       ! lcc



! logical switch controls whether facevar data are being packed
      invar = nfacevar
      ngcell_on_fc0 = maxval(ngcell_on_fc(1:ndim))
      if (lfc.and.lguard_in_progress) & 
     &         invar = ngcell_on_fc0
      if (lfc) then

! pack the facevarx array for block lb

      vtype = 2
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

! pack the facevary array for block lb

      if(ndim.ge.2) then

      vtype = 3
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      endif

! pack the facevarz array for block lb

      if(ndim.eq.3.or.l2p5d.eq.1) then

      vtype = 4
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      endif

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for fc ',index
#endif /* DEBUG */

      endif                           ! lfc

! logical switch controls whether unk_n data are being packed
      invar = nvaredge
      ngcell_on_ec0 = maxval(ngcell_on_ec(1:ndim))
      if (lec.and.lguard_in_progress) & 
     &         invar = ngcell_on_ec0
      if (lec .and. ndim.gt.0) then

! pack the unk_e_x array for block lb

      vtype = 5
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      if (ndim > 1) then
! pack the unk_e_y array for block lb

      vtype = 6
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      if (ndim == 3) then
! pack the unk_e_z array for block lb

      vtype = 7
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      endif
      end if
#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for ec ',index
#endif /* DEBUG */
      endif                       ! lec


! logical switch controls whether unk_n data are being packed
      invar = nvarcorn
      if (lnc.and.lguard_in_progress) & 
     &         invar = ngcell_on_nc
      if (lnc) then

! pack the unk_n array for block lb

      vtype = 8
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for nc ',index
#endif /* DEBUG */
      endif                       ! lnc

      endif                       ! ioptw


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

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype, & 
     &                  ' tree info unpacked into block ',lb
#endif /* DEBUG */

! overwrite outgoing index

      offset = index

      return
      end subroutine mpi_put_buffer

!------------------------------------------------------------------------

      subroutine mpi_get_Rbuffer_size(lb,dtype,ioptw,offset, & 
     &                                lcc,lfc,lec,lnc, & 
     &                                nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine unpacks the buffer R_buffer which has been received on 
! this processor into the local guard block lb of the arrays work,unk,
! and facevar.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      lb             guard block id to be unpacked from receiving
!      ioptw          option setting for work array
!      offset         offset for buffer index
!      lcc            logical switch controlling whether unk or work data
!                     is unpacked
!      lfc            logical switch controlling whether facevar data
!                     is unpacked
!      lec            logical switch controlling whether unk_e_? or work data
!                     is unpacked
!      lnc            logical switch controlling whether unk_n data
!                     is unpacked
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

      integer, intent(in)    :: lb,dtype,ioptw
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer :: nguard0 
      integer :: nguard_work0 
      integer :: index
      integer :: vtype
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ii, jj, k, iii, invar
      integer :: ngcell_on_fc0,ngcell_on_ec0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DEBUG

      integer :: mype, ierr

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
#endif /* DEBUG */


      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

! take over incoming index offset for block lb to be sent to remote pe

      index = offset          
      index = index+3

! pack the work array for block lb
      if(ioptw.gt.1) then

      vtype = 0
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + (kb-ka+1)*(jb-ja+1)*(ib-ia+1)


      else                           ! ioptw



#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index on entry ',index
#endif /* DEBUG */


! logical switch controls whether unk data are being packed
      invar = nvar
      if (lcc.and.lguard_in_progress) & 
     &         invar = ngcell_on_cc
      if (lcc) then

! pack the unk array for block lb

      vtype = 1
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for cc ',index, & 
     &         ' invar ',invar,' ia ib ja jb ka kb ',ia,ib,ja,jb,ka,kb, & 
     &         ' dtype ',dtype
#endif /* DEBUG */

      endif                       ! lcc



! logical switch controls whether facevar data are being packed
      invar = nfacevar
      ngcell_on_fc0 = maxval(ngcell_on_fc(1:ndim))
      if (lfc.and.lguard_in_progress) & 
     &         invar = ngcell_on_fc0
      if (lfc) then

! pack the facevarx array for block lb

      vtype = 2
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

! pack the facevary array for block lb

      if(ndim.ge.2) then

      vtype = 3
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      endif

! pack the facevarz array for block lb

      if(ndim.eq.3.or.l2p5d.eq.1) then

      vtype = 4
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      endif

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for fc ',index
#endif /* DEBUG */

      endif                           ! lfc

! logical switch controls whether unk_n data are being packed
      invar = nvaredge
      ngcell_on_ec0 = maxval(ngcell_on_ec(1:ndim))
      if (lec.and.lguard_in_progress) & 
     &         invar = ngcell_on_ec0
      if (lec .and. ndim.gt.0) then

! pack the unk_e_x array for block lb

      vtype = 5
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      if (ndim > 1) then
! pack the unk_e_y array for block lb

      vtype = 6
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      if (ndim == 3) then
! pack the unk_e_z array for block lb

      vtype = 7
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

      endif
      end if
#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for ec ',index
#endif /* DEBUG */
      endif                       ! lec


! logical switch controls whether unk_n data are being packed
      invar = nvarcorn
      if (lnc.and.lguard_in_progress) & 
     &         invar = ngcell_on_nc
      if (lnc) then

! pack the unk_n array for block lb

      vtype = 8
      call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                            nlayersx,nlayersy,nlayersz)
      index = index + invar*(kb-ka+1)*(jb-ja+1)*(ib-ia+1)

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype,' index update for nc ',index
#endif /* DEBUG */
      endif                       ! lnc

      endif                       ! ioptw


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

#ifdef DEBUG
      write(*,*) 'put_buffer : pe ',mype, & 
     &                  ' tree info unpacked into block ',lb
#endif /* DEBUG */

! overwrite outgoing index

      offset = index

      return
      end subroutine mpi_get_Rbuffer_size
