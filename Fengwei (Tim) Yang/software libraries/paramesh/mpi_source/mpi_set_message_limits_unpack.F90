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

      subroutine mpi_set_message_limits_unpack & 
     &                 (dtype,ia,ib,ja,jb,ka,kb,vtype)

!------------------------------------------------------------------------
!
! This subroutine sets up the range of indeces to be set when
! reading a message directly out of R_buffer into unk1 or work1.
!
!
! Written :     Peter MacNeice          April 2001
!------------------------------------------------------------------------
!
! Arguments:
!      dtype          sets message type, ie what section of the block is
!                     required. Must be a number between 1 and 27. Note
!                     this index refers to the portion of the working
!                     array unk1, unk_n1, etc which is to be filled.
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


      integer, intent(in)  ::  dtype,vtype
      integer, intent(out) ::  ia,ib,ja,jb,ka,kb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(vtype.eq.0) then
         nguard0 = nguard_work
      else
         nguard0 = nguard
      endif

!
! define starting end ending indices for the send and recv buffers


! set i to message type
        i = dtype

! set x index extent
        if(mod(i,3).eq.1) then
          ia = 1
          ib = nguard0
        elseif(mod(i,3).eq.2) then
          ia = 1+nguard0
          ib = nxb+nguard0
        elseif(mod(i,3).eq.0) then
          ia = nxb + nguard0 + 1
          ib = nxb + 2*nguard0
        endif
! set y index extent
        if(mod((i-1)/3,3).eq.0) then
          ja = 1
          jb = (nguard0-1)*k2d+1
        elseif(mod((i-1)/3,3).eq.1) then
          ja = 1+nguard0*k2d
          jb = nyb+nguard0*k2d
        elseif(mod((i-1)/3,3).eq.2) then
          ja = (nyb + nguard0)*k2d + 1
          jb = nyb + 2*nguard0*k2d
        endif
! set z index extent
        if(i.le.9) then
          ka = 1
          kb = (nguard0-1)*k3d+1
        elseif(i.ge.10.and.i.le.18) then
          ka = 1+nguard0*k3d
          kb = nzb+nguard0*k3d
        elseif(i.ge.19) then
          ka = (nzb + nguard0)*k3d + 1
          kb = nzb + 2*nguard0*k3d
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

      return
      end subroutine mpi_set_message_limits_unpack
