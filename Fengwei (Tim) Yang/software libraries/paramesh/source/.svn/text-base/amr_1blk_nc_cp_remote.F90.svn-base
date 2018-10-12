!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_nc_cp_remote
!! NAME
!!   
!!   amr_1blk_nc_cp_remote
!!
!! SYNOPSIS
!!
!!    Call amr_1blk_nc_cp_remote(mype,remote_pe,remote_block, 
!!          idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1, 
!!          ip3,jp3,kp3,nblk_ind)
!!    Call amr_1blk_nc_cp_remote(integer,integer,integer,
!!          integer,integer,integer,integer,integer,integer,integer,
!!           integer,integer,integer,integer,integer,integer,
!!          integer,integer,integer,integer)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in) :: mype local processor
!!   Integer, Intent(in) :: remote_pe remote processor
!!   Integer, Intent(in) :: remote_block local block id of the block to be copied from
!!                                       the remote processor
!!   Integer, Intent(in) :: idest selects the storage space in data_1blk.fh which is to
!!                                be used in this call. If the leaf node is having its
!!                                guardcells filled then set this to 1, if its parent
!!                                is being filled set it to 2.
!!   Integer, Intent(in) :: id lower limit of index range of points in x direction
!!                             on destination block
!!   Integer, Intent(in) :: jd lower limit of index range of points in y direction
!!                             on destination block
!!   Integer, Intent(in) :: kd lower limit of index range of points in z direction
!!                             on destination block
!!   Integer, Intent(in) :: is lower limit of index range of points in x direction
!!                             on source block
!!   Integer, Intent(in) :: js lower limit of index range of points in y direction
!!                             on source block
!!   Integer, Intent(in) :: ks lower limit of index range of points in z direction
!!                             on source block
!!   Integer, Intent(in) :: ilays no. of mesh points in x direction to be copied
!!   Integer, Intent(in) :: jlays no. of mesh points in y direction to be copied
!!   Integer, Intent(in) :: klays no. of mesh points in z direction to be copied
!!   Integer, Intent(in) :: ip1 offset added to index range defined by (id,ilays)
!!                              0 if guardcells are at lower end of i index
!!                              1 if guardcells are at upper end of i index
!!   Integer, Intent(in) :: jp1 offset added to index range defined by (jd,jlays)
!!                              0 if guardcells are at lower end of j index
!!                              1 if guardcells are at upper end of j index
!!   Integer, Intent(in) :: kp1 offset added to index range defined by (kd,klays)
!!                              0 if guardcells are at lower end of k index
!!                              1 if guardcells are at upper end of k index
!!   Integer, Intent(in) :: ip3 array index adjustment
!!   Integer, Intent(in) :: jp3 array index adjustment
!!   Integer, Intent(in) :: kp3 array index adjustment
!!   Integer, Intent(in) :: nblk_ind block index is doing fine to course
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   mpi_morton
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_mpi_find_blk_in_buffer
!!   mpi_set_message_limits
!!   
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine copies guard cell information for cell corner
!!   data to layer idest of unk_n, from the appropriate
!!   corner data of the neighboring block.
!!
!! AUTHOR
!!
!! Written :     Peter MacNeice          December 2000
!!
!!***

#include "paramesh_preprocessor.fh"
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
!#define DEBUG

      subroutine amr_1blk_nc_cp_remote(mype,remote_pe,remote_block,    & 
         idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,        & 
         ip3,jp3,kp3,nblk_ind)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits

      Implicit None

!-----Input/Output Arguments
      Integer, Intent(in) :: mype,remote_pe,remote_block
      Integer, Intent(in) :: idest,id,jd,kd,is,js,ks
      Integer, Intent(in) :: ilays,jlays,klays
      Integer, Intent(in) :: ip1,jp1,kp1,ip3,jp3,kp3
      Integer, Intent(in) :: nblk_ind

!-----Local arays and variables
      Integer :: il,jl,kl,id1,jd1,kd1,is1,js1,ks1
      Integer :: ill,jll,kll
      Integer :: index, dtype 
      Integer :: ii, jj, kk, i, j, k
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: ivar, ivar_next
      Logical :: lfound
      Integer :: vtype

!-----Begin Executable Code

!-----Adjust index ranges
      il = ilays - ip3
      jl = jlays*k2d - jp3*k2d
      kl = klays*k3d - kp3*k3d

      id1 = id + ip1
      jd1 = jd + jp1*k2d
      kd1 = kd + kp1*k3d
      is1 = is + ip1
      js1 = js + jp1*k2d
      ks1 = ks + kp1*k3d

      If (remote_block <= lnblocks.and.remote_pe == mype) Then

       if (no_permanent_guardcells) Then

       unk_n1(1:nvarcorn,id1:id1+il,jd1:jd1+jl,                        & 
                                        kd1:kd1+kl,idest)              & 
          =  gt_unk_n(1:nvarcorn,is1:is1+il,js1:js1+jl,                & 
                                 ks1:ks1+kl,                           & 
                                 remote_block)
       Else ! no_permanent_guardcells

       unk_n1(1:nvarcorn,id1:id1+il,jd1:jd1+jl,                        & 
                                        kd1:kd1+kl,idest)              & 
          =  unk_n(1:nvarcorn,is1:is1+il,js1:js1+jl,                   & 
                              ks1:ks1+kl,                              & 
                              remote_block)

       End If ! End If (no_permanent_guardcells)

!--
      Else                          ! otherwise block is remote
!--

        Call amr_mpi_find_blk_in_buffer(mype,remote_block,             & 
                              remote_pe,idest,dtype,index,lfound)

#ifdef DEBUG
         Write(*,*) 'pe ',mype,' find blk in buff ',remote_block,      & 
                        remote_pe,' lfound ',lfound,' dtype ',dtype,   & 
                        ' index ',index,' original address ',          & 
                        laddress(:,remote_block)
#endif /* DEBUG */

!-------If this routine is executing a copy to fill guardcells of a
!-------leaf blocks^s parent, and the remote block is not found, Then
!-------it is assumed that it is not in the list of buffered remote blocks
!-------because it is not really needed. Therefore in this case we
!-------return without copying anything.
        If (idest == 2.and.(.Not.lfound)) Return

!-------starting index if cell-centered data is also included in recv_buf
        If (l_datapacked(2)) index =                                   & 
                     index + ngcell_on_cc*message_size_cc(dtype)
        If (l_datapacked(3)) index =                                   & 
                            index + maxval(ngcell_on_fc(1:ndim))       & 
                                    *message_size_fc(dtype)
        If (l_datapacked(4)) index =                                   & 
                            index + maxval(ngcell_on_ec(1:3))          & 
                                    *message_size_ec(dtype)

        ill = ilays
        jll = jlays
        kll = klays

        vtype = 8
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,                    & 
                     ill,jll,kll)

        kk = kd1
        Do k = ka,kb
        jj = jd1
        Do j = ja,jb
        ii = id1
        Do i = ia,ib
          If (k >= ks1 .and. k <= ks1 + kl) Then
          If (j >= js1 .and. j <= js1 + jl) Then
          If (i >= is1 .and. i <= is1 + il) Then

          Do ivar=1,ngcell_on_nc
           ivar_next = gcell_on_nc_pointer(ivar)

           unk_n1(ivar_next,ii,jj,kk,idest) =                          & 
                     temprecv_buf(index+ivar)
          End Do  ! End Do ivar=1,ngcell_on_nc

          End If  ! End If (i >= is1 .and. i <= is1 + il)
          End If  ! End If (j >= js1 .and. j <= js1 + jl)
          End If  ! End If (k >= ks1 .and. k <= ks1 + kl)

          If (i >= is1 .and. i <= is1 + il) ii = ii + 1
          index = index+ngcell_on_nc

#ifdef DEBUG
         Write(*,*) 'pe ',mype,' unpacked recv for unk_n1 ',i,j,k,     & 
            unk_n1(1:nvarcorn,ii,jj,kk,idest),' index ',index+1
#endif /* DEBUG */

        End Do  ! End Do i = ia,ib
        If (j >= js1 .and. j <= js1 + jl) jj = jj + 1
        End Do  ! End Do j = ja,jb
        If (k >= ks1 .and. k <= ks1 + kl) kk = kk + 1
        End Do  ! End Do k = ka,kb

      End If  ! End If (remote_block <= lnblocks.and.remote_pe == mype)

      Return
      End Subroutine amr_1blk_nc_cp_remote
