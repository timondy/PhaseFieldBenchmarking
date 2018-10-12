!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_ec_cp_remote
!! NAME
!!
!!   amr_1blk_ec_cp_remote
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_ec_cp_remote(mype,remote_pe,remote_block,  
!!           idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1, 
!!           ip2,jp2,kp2,ip3,jp3,kp3,iface,nblk_ind)
!!   Call amr_1blk_ex_cp_remote(All integers)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in) :: mype             
!!     local processor
!!
!!   Integer, Intent(in) :: remote_pe        
!!     remote processor
!!
!!   Integer, Intent(in) :: remote_block     
!!     local block id of the block to be copied from the remote processor
!!
!!   Integer, Intent(in) :: idest           
!!     selects the storage space in data_1blk.fh which is to
!!     be used in this call. If the leaf node is having its
!!     guardcells filled then set this to 1, if its parent
!!      is being filled set it to 2.
!!
!!   Integer, Intent(in) :: id               
!!     lower limit of index range of points in x direction
!!     on destination block
!!
!!   Integer, Intent(in) :: jd               
!!     lower limit of index range of points in y direction
!!     on destination block
!!
!!   Integer, Intent(in) :: kd               
!!     lower limit of index range of points in z direction
!!     on destination block
!!   
!!   Integer, Intent(in) :: is               
!!     lower limit of index range of points in x direction
!!     on source block
!!
!!   Integer, Intent(in) :: js               
!!     lower limit of index range of points in y direction
!!     on source block
!!
!!   Integer, Intent(in) :: ks               
!!     lower limit of index range of points in z direction
!!     on source block
!!
!!   Integer, Intent(in) :: ilay             
!!     no. of mesh points in x direction to be copied
!!
!!   Integer, Intent(in) :: jlay             
!!     no. of mesh points in y direction to be copied
!!
!!   Integer, Intent(in) :: klay             
!!     no. of mesh points in z direction to be copied
!!
!!   Integer, Intent(in) :: ip1              
!!     i index offset dependent on the face being treated,
!!      1 if at low or high x range, 0 otherwise.
!!
!!   Integer, Intent(in) :: jp1              
!!     j index offset dependent on the face being treated,
!!      1 if at low or high y range, 0 otherwise.
!!
!!   Integer, Intent(in) :: kp1              
!!     k index offset dependent on the face being treated,
!!      1 if at low or high z range, 0 otherwise.
!!
!!   Integer, Intent(in) :: ip2              
!!     extend range in i coord for unk_e_x by this amount
!!      must be set to either 1 or 0 
!!      (should be 1 for face 2, otherwise 0.)
!!
!!   Integer, Intent(in) :: jp2              
!!     extend range in j coord for unk_e_y by this amount
!!      must be set to either 1 or 0
!!      (should be 1 for face 4, otherwise 0.)
!!
!!   Integer, Intent(in) :: kp2              
!!     extend range in k coord for unk_e_z by this amount
!!       must be set to either 1 or 0
!!      (should be 1 for face 6, otherwise 0.)
!!
!!   Integer, Intent(in) :: ip3
!!     adjustment for range in i coord
!!
!!   Integer, Intent(in) :: jp3
!!     adjustment for range in j coord
!!
!!   Integer, Intent(in) :: kp3
!!     adjustment for range in k coord
!!
!!   Integer, Intent(in) :: iface            
!!     set between 1 and 6 if working on a block face,
!!     otherwise set to 0.
!!
!!   Integer, Intent(in) :: nblk_ind
!!     A block index which is used when l_f_to_c is set to TRUE.
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
!!   Nothing explicitly returned.  Data is copied into edge centered arrays
!!   from previously comminicated buffers which contain information from
!!   off processor blocks.
!!
!! DESCRIPTION
!!
!!   This routine copies guard cell information for cell edge centered
!!   data to layer idest of unk_e_x, unk_e_y and unk_e_z, from the appropriate
!!   edge data of the neighboring block.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          December 2000
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_ec_cp_remote(mype,remote_pe,remote_block,  & 
         idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,      & 
         ip2,jp2,kp2,ip3,jp3,kp3,iface,nblk_ind)


!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits

      Implicit None

!-----Input/Output Variables
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3,iface
      integer, intent(in) :: nblk_ind

!-----Local Variables
      integer :: il,jl,kl
      integer :: ill,jll,kll
      integer :: il1, jl1, kl1, id1, jd1, kd1
      integer :: is1, js1, ks1, index
      integer :: ii, jj, kk, i, j, k
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ivar, ivar_next
      integer :: dtype
      integer :: vtype

      logical :: lfound

!-----Begin Executable Code

      If (ndim > 1) Then

!-----Adjust index ranges
      il = ilays - ip3
      jl = jlays*k2d - jp3*k2d
      kl = klays*k3d - kp3*k3d

      il1 = il-1
      jl1 = (jl-1)*k2d
      kl1 = (kl-1)*k3d

      id1 = id + ip1
      jd1 = jd + jp1*k2d
      kd1 = kd + kp1*k3d
      is1 = is + ip1
      js1 = js + jp1*k2d
      ks1 = ks + kp1*k3d

      If (remote_block <= lnblocks .and. remote_pe == mype) Then

       If (no_permanent_guardcells) Then

       unk_e_x1(1:nbndvare,                              & 
                id:id+il1+ip2,                           & 
                jd1:jd1+jl,                              & 
                kd1:kd1+kl,                              & 
                idest)                                   & 
          =  gt_unk_e_x(1:nbndvare,                      & 
                        is:is+il1+ip2,                   & 
                        js1:js1+jl,                      & 
                        ks1:ks1+kl,remote_block)

       unk_e_y1(1:nbndvare,                              & 
                id1:id1+il,                              & 
                jd:jd+jl1+jp2*k2d,                       & 
                kd1:kd1+kl,                              & 
                idest)                                   &  
          =  gt_unk_e_y(1:nbndvare,                      & 
                        is1:is1+il,                      & 
                        js:js+jl1+jp2*k2d,               & 
                        ks1:ks1+kl,remote_block)

       If (ndim == 3) Then
       unk_e_z1(1:nbndvare,                              & 
                id1:id1+il,                              & 
                jd1:jd1+jl,                              & 
                kd:kd+kl1+kp2*k3d,                       & 
                idest) & 
          =  gt_unk_e_z(1:nbndvare,                      & 
                        is1:is1+il,                      & 
                        js1:js1+jl,                      & 
                        ks:ks+kl1+kp2*k3d,remote_block)
       End If  ! End If (ndim == 3)

      Else ! If (no_permanent_guardcells)

       unk_e_x1(1:nbndvare,                              & 
                id:id+il1+ip2,                           & 
                jd1:jd1+jl,                              & 
                kd1:kd1+kl,                              & 
                idest)                                   & 
          =  unk_e_x(1:nbndvare,                         & 
                     is:is+il1+ip2,                      & 
                     js1:js1+jl,                         & 
                     ks1:ks1+kl,remote_block)

       unk_e_y1(1:nbndvare,                              & 
                id1:id1+il,                              & 
                jd:jd+jl1+jp2*k2d,                       & 
                kd1:kd1+kl,                              & 
                idest)                                   & 
          =  unk_e_y(1:nbndvare,                         & 
                     is1:is1+il,                         & 
                     js:js+jl1+jp2*k2d,                  & 
                     ks1:ks1+kl,remote_block)

       If (ndim == 3) Then
       unk_e_z1(1:nbndvare,                              & 
                id1:id1+il,                              & 
                jd1:jd1+jl,                              & 
                kd:kd+kl1+kp2*k3d,                       & 
                idest)                                   & 
          =  unk_e_z(1:nbndvare,                         & 
                     is1:is1+il,                         & 
                     js1:js1+jl,                         & 
                     ks:ks+kl1+kp2*k3d,remote_block)
       End If  ! End If (ndim == 3)

      End If ! If (no_permanent_guardcells)

      Else  ! If (remote_block <= lnblocks .and. remote_pe == mype)
            ! In other words, this block is remote

        Call amr_mpi_find_blk_in_buffer(mype,remote_block,        & 
                                        remote_pe,idest,dtype,    &
                                        index,lfound)

!-------If this routine is executing a copy to fill guardcells of a
!-------leaf blocks^s parent, and the remote block is not found, then
!-------it is assumed that it is not in the list of buffered remote blocks
!-------because it is not really needed. Therefore in this case we
!-------return without copying anything.

        If (idest == 2 .and. (.Not.lfound)) Return

!-------starting index if cell-centered data is also included in recv_buf
        If (l_datapacked(2)) index =                                   & 
                             index + ngcell_on_cc*message_size_cc(dtype)
        If (l_datapacked(3)) index =                                   & 
                             index + maxval(ngcell_on_fc(1:ndim))      & 
                                         *message_size_fc(dtype)

        ill = ilays
        jll = jlays
        kll = klays

        vtype = 5
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,                    & 
                     ill,jll,kll)

        kk = kd1
        Do k = ka,kb
        jj = jd1
        Do j = ja,jb
        ii = id
        Do i = ia,ib
          If (k >= ks1 .and. k <= ks1 + kl) Then
          If (j >= js1 .and. j <= js1 + jl) Then
          If (i >= is  .and. i <= is + il1 + ip2) Then

          Do ivar=1,ngcell_on_ec(1)
            ivar_next = gcell_on_ec_pointer(1,ivar)
             unk_e_x1(ivar_next,ii,jj,kk,idest) =                      & 
                      temprecv_buf(index+ivar)
          End Do  ! End Do ivar=1,ngcell_on_ec(1)

          endif  ! End If (i >= is  .and. i <= is + il1 + ip2)
          endif  ! End If (j >= js1 .and. j <= js1 + jl)
          endif  ! End If (k >= ks1 .and. k <= ks1 + kl)
          if (i >= is .and. i <= is + il1 + ip2) ii = ii + 1
          index = index+ngcell_on_ec(1)
        End Do  ! End Do i = ia,ib
        if (j >= js1 .and. j <= js1 + jl) jj = jj + 1
        End Do  ! End Do j = ja,jb
        if (k >= ks1 .and. k <= ks1 + kl) kk = kk + 1
        End Do  ! End Do k = ka,kb

        If (ndim >= 2) Then

        vtype = 6
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,                    & 
                     ill,jll,kll)

        kk = kd1
        Do k = ka,kb
        jj = jd
        Do j = ja,jb
        ii = id1
        Do i = ia,ib
          If (k >= ks1 .and. k <= ks1 + kl) Then
          If (j >= js  .and. j <= js  + jl1 + jp2*k2d) Then
          If (i >= is1 .and. i <= is1 + il) Then

          Do ivar=1,ngcell_on_ec(2)
            ivar_next = gcell_on_ec_pointer(2,ivar)
             unk_e_y1(ivar_next,ii,jj,kk,idest) =                      & 
                      temprecv_buf(index+ivar)
          End Do  ! Do ivar=1,ngcell_on_ec(2)

          End If  ! End If (i >= is1 .and. i <= is1 + il)
          End If  ! End If (j >= js  .and. j <= js  + jl1 + jp2*k2d)
          End If  ! End If (k >= ks1 .and. k <= ks1 + kl)
          If (i >= is1 .and.  i <= is1 + il) ii = ii + 1
          index = index+ngcell_on_ec(2)
        End Do  ! End Do k = ia,ib
        If (j >= js  .and.  j <= js + jl1 + jp2*k2d) jj = jj + 1
        End Do  ! End Do j = ja,jb
        if (k >= ks1 .and. k <= ks1 + kl) kk = kk + 1
        End Do  ! End Do k = ka,kb
        End If  ! End If (ndim >= 2)
 
        If (ndim == 3 .or. l2p5d == 1) Then

        vtype =7 
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,                    & 
                     ill,jll,kll)

        kk = kd
        Do k = ka,kb
        jj = jd1
        Do j = ja,jb
        ii = id1
        Do i = ia,ib
          If (k >= ks  .and.  k <= ks + kl1 + kp2*k3d) Then
          If (j >= js1 .and.  j <= js1 + jl) Then
          If (i >= is1 .and.  i <= is1 + il) Then

          Do ivar=1,ngcell_on_ec(3)
            ivar_next = gcell_on_ec_pointer(3,ivar)
            unk_e_z1(ivar_next,ii,jj,kk,idest) =                       & 
                      temprecv_buf(index+ivar)
          End Do  ! End Do ivar=1,ngcell_on_ec(3)

          End If  ! End If (i >= is1 .and.  i <= is1 + il)
          End If  ! End If (j >= js1 .and.  j <= js1 + jl)
          End If  ! End If (k >= ks  .and.  k <= ks + kl1 + kp2*k3d)
          If (i >= is1 .and.  i <= is1 + il) ii = ii + 1
          index = index+ngcell_on_ec(3)
        End Do  ! End Do i = ia,ib
        if (j >= js1 .and.  j <= js1 + jl) jj = jj + 1
        End Do  ! End Do j = ja,jb
        if (k >= ks .and. k <= ks + kl1 + kp2*k3d) kk = kk + 1
        End Do  ! End Do k = ka,kb
        End If  ! End If (ndim == 3 .or. l2p5d == 1)

        End If  ! End If (remote_block <= lnblocks .and. remote_pe == mype)

#ifdef XFORCE_CONSISTENCY_AT_SRL_INTERFACES
!------This is not yet complete. Need to cater for cases where 4 blocks share
!------A common edge.

       If (iface == 1) Then

        Do ivar=1,ngcell_on_ec(2)
         ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d,    & 
                  1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(      & 
         unk_e_y1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d,    & 
                             1+nguard*k3d:nzb+(1+nguard)*k3d,idest)  & 
         + recvy(ivar_next,nxb+1,1:nyb,1:nzb+k3d) )
        End Do  ! End Do ivar=1,ngcell_on_ec(2)

        If (ndim == 3) Then
        Do ivar=1,ngcell_on_ec(3)
         ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,1+nguard,1+nguard*k2d:nyb+(1+nguard)*k2d, & 
                  1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(           & 
         unk_e_z1(ivar_next,1+nguard,1+nguard*k2d:nyb+(1+nguard)*k2d, & 
                  1+nguard*k3d:nzb+nguard*k3d,idest)                  & 
         + recvz(ivar_next,nxb+1,1:nyb+k2d,1:nzb) )
        End Do  ! End Do ivar=1,ngcell_on_ec(3)
        End If  ! End If (ndim == 3)

       Elseif (iface == 2) Then

        Do ivar=1,ngcell_on_ec(2)
         ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
                  1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(       & 
         unk_e_y1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
                  1+nguard*k3d:nzb+(1+nguard)*k3d,idest)              & 
         + recvy(ivar_next,1,1:nyb,1:nzb+k3d) )
        End Do  ! End Do ivar=1,ngcell_on_ec(2)

        If (ndim == 3) Then
        Do ivar=1,ngcell_on_ec(3)
         ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,nxb+1+nguard,                             & 
                  1+nguard*k2d:nyb+(1+nguard)*k2d,                    & 
                  1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(           & 
         unk_e_z1(ivar_next,nxb+1+nguard,                             & 
                  1+nguard*k2d:nyb+(1+nguard)*k2d,                    & 
                  1+nguard*k3d:nzb+nguard*k3d,idest)                  &
        + recvz(ivar_next,1,1:nyb+k2d,1:nzb) )
        End Do  ! End Do ivar=1,ngcell_on_ec(3)
        End If  ! End If (ndim == 3)

       Elseif (iface == 3) Then

        Do ivar=1,ngcell_on_ec(1)
         ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d,         & 
                  1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(       & 
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d,         & 
                  1+nguard*k3d:nzb+(1+nguard)*k3d,idest)              & 
         + recvx(ivar_next,1:nxb,nyb+k2d,1:nzb+k3d) )
        End Do  ! End Do ivar=1,ngcell_on_ec(1)

        If (ndim == 3) Then
        Do ivar=1,ngcell_on_ec(3)
         ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,1+nguard*k2d,       & 
                  1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(           & 
         unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,1+nguard*k2d,       & 
                  1+nguard*k3d:nzb+nguard*k3d,idest)                  & 
         + recvz(ivar_next,1:nxb+1,1,1:nzb) )
        End Do  ! End Do ivar=1,ngcell_on_ec(3)
        End If  ! End If (ndim == 3)

       Elseif (iface == 4) Then

        Do ivar=1,ngcell_on_ec(1)
         ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,nyb+(1+nguard)*k2d,   & 
                  1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(       & 
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,nyb+(1+nguard)*k2d,   & 
                  1+nguard*k3d:nzb+(1+nguard)*k3d,idest)              & 
         + recvx(ivar_next,1:nxb,1,1:nzb+k3d) )
        End Do  ! End Do ivar=1,ngcell_on_ec(1)

        If (ndim == 3) Then
        Do ivar=1,ngcell_on_ec(3)
         ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,nyb+(1+nguard)*k2d, & 
                  1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(           & 
         unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,nyb+(1+nguard)*k2d, & 
                  1+nguard*k3d:nzb+nguard*k3d,idest)                  & 
         + recvz(ivar_next,1:nxb+1,1,1:nzb) )
        End Do
        End If  ! End If (ndim == 3)

       Elseif (iface == 5 .and. ndim == 3) Then

        do ivar=1,ngcell_on_ec(1)
         ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,                      & 
                  1+nguard*k2d:nyb+(1+nguard)*k2d,                    & 
                  1+nguard*k3d,idest) = .5*(                          & 
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,                      & 
                  1+nguard*k2d:nyb+(1+nguard)*k2d,                    & 
                  1+nguard*k3d,idest)                                 & 
         + recvx(ivar_next,1:nxb,1:nyb+k2d,nzb) )
        End Do  ! End do ivar=1,ngcell_on_ec(1)

        Do ivar=1,ngcell_on_ec(2)
         ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,                    & 
                  1+nguard*k2d:nyb+nguard*k2d,                        & 
                  1+nguard*k3d,idest) = .5*(                          & 
         unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,                    & 
                  1+nguard*k2d:nyb+nguard*k2d,                        & 
                  1+nguard*k3d,idest)                                 & 
         + recvy(ivar_next,1:nxb+1,1:nyb,nzb) )
        End Do  ! End Do ivar=1,ngcell_on_ec(2)

       Elseif(iface == 6 .and. ndim == 3) Then

        Do ivar=1,ngcell_on_ec(1)
         ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,                      & 
                  1+nguard*k2d:nyb+(1+nguard)*k2d,                    & 
                  nzb+(1+nguard)*k3d,idest) = .5*(                    & 
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,                      & 
                  1+nguard*k2d:nyb+(1+nguard)*k2d,                    &  
                  nzb+(1+nguard)*k3d,idest)                           & 
         + recvx(ivar_next,1:nxb,1:nyb+k2d,1) )
        End Do  ! End Do ivar=1,ngcell_on_ec(1)

        Do ivar=1,ngcell_on_ec(2)
         ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,                    & 
                  1+nguard*k2d:nyb+nguard*k2d,                        & 
                  nzb+(1+nguard)*k3d,idest) = .5*(                    & 
         unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,                    & 
                  1+nguard*k2d:nyb+nguard*k2d,                        & 
                  nzb+(1+nguard)*k3d,idest)                           & 
         + recvy(ivar_next,1:nxb+1,1:nyb,1) )
        End Do  ! End Do ivar=1,ngcell_on_ec(2)

       End If  ! End If (iface == 1)
#endif /* FORCE_CONSISTENCY_AT_SRL_INTERFACES */

      End If ! End If (ndim > 1)

      Return
      End Subroutine amr_1blk_ec_cp_remote


