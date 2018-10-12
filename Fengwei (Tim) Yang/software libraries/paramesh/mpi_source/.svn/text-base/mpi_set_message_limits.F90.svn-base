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

      subroutine mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                                  nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine sets up the range of indeces to be fetched from a
! remote block, when only a part of a block is to be fetched.
!
!
! Written :     Peter MacNeice          April 2001
!------------------------------------------------------------------------
!
! Arguments:
!      dtype          sets message type, ie what section of the block is
!                     required. Must be a number between 1 and 27.
!      ia             lower range for i index.
!      ib             upper range for i index.
!      ja             lower range for j index.
!      jb             upper range for j index.
!      ka             lower range for k index.
!      kb             upper range for k index.
!      vtype          sets variable type.
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree

      use mpi_morton

      implicit none

      include 'mpif.h'

      integer, intent(in)  ::  dtype,vtype
      integer, intent(out) ::  ia,ib,ja,jb,ka,kb
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

! local variables

      integer :: ierrorcode,ierr
      integer :: nlayers0x, nlayers0y, nlayers0z
      integer :: i, nguard0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(vtype.gt.8) then
        write(*,*) 'Paramesh error: vtype too large.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

! set i to message type
      i = dtype

      if(vtype.ne.0) then
         if (.not.present(nlayersx)) then
            nlayers0x = nguard
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard
         else
            nlayers0z = nlayersz
         end if
         if (i > 27) then
            i = i - 27
            nlayers0x = min(nlayers0x*2+1,nguard*2)
            nlayers0y = min(nlayers0y*2+1,nguard*2)
            nlayers0z = min(nlayers0z*2+1,nguard*2)
         end if
      else
         if (.not.present(nlayersx)) then
            nlayers0x = nguard_work
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard_work
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard_work
         else
            nlayers0z = nlayersz
         end if
         if (i > 27) then
            i = i - 27
            nlayers0x = min(nlayers0x*2+1,nguard_work*2)
            nlayers0y = min(nlayers0y*2+1,nguard_work*2)
            nlayers0z = min(nlayers0z*2+1,nguard_work*2)
         end if
      endif

      if(vtype.eq.0) then
         nguard0 = nguard_work
      else
         nguard0 = nguard
      endif

!
! define starting end ending indices for the send and recv buffers

! set x index extent
        if(mod(i,3).eq.1) then
          ia = 1
          ib = nlayers0x
        elseif(mod(i,3).eq.2) then
          ia = 1
          ib = nxb
        elseif(mod(i,3).eq.0) then
          ia = nxb - nlayers0x + 1
          ib = nxb
        endif
! set y index extent
        if(mod((i-1)/3,3).eq.0) then
          ja = 1
          jb = (nlayers0y-1)*k2d+1
        elseif(mod((i-1)/3,3).eq.1) then
          ja = 1
          jb = nyb
        elseif(mod((i-1)/3,3).eq.2) then
          ja = (nyb - nlayers0y)*k2d + 1
          jb = nyb
        endif
! set z index extent
        if(i.le.9) then
          ka = 1
          kb = (nlayers0z-1)*k3d+1
        elseif(i.ge.10.and.i.le.18) then
          ka = 1
          kb = nzb
        elseif(i.ge.19) then
          ka = (nzb - nlayers0z)*k3d + 1
          kb = nzb
        endif

!  Now adjust bounds appropriately for the variable being
!  considered
      select case(vtype)
        case(0)                    ! work variable

        case(1)                    ! unk

        case(2)                    ! facevarx
          ib = ib+1
        case(3)                    ! facevary
          jb = jb+k2d
        case(4)                    ! facevarz
          kb = kb+k3d
        case(5)                    ! unk_e_x
          jb = jb+k2d
          kb = kb+k3d
        case(6)                    ! unk_e_y
          ib = ib+1
          kb = kb+k3d
        case(7)                    ! unk_e_z
          ib = ib+1
          jb = jb+k2d
        case(8)                    ! unk_n
          ib = ib+1
          jb = jb+k2d
          kb = kb+k3d
      end select


      if (.not.no_permanent_guardcells) then
! If permanent guardcell storage is allocated then we need to offset
! the data to be fetched by the number of guardcells at the left block
! boundaries.
          ia = ia + nguard0
          ib = ib + nguard0
          ja = ja + nguard0*k2d
          jb = jb + nguard0*k2d
          ka = ka + nguard0*k3d
          kb = kb + nguard0*k3d
       end if

      return
      end subroutine mpi_set_message_limits
