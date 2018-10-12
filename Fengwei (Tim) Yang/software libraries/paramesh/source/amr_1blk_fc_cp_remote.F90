!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_fc_cp_remote
!! NAME
!!
!!   amr_1blk_fc_cp_remote
!!
!! SYNOPSIS
!! 
!!  Call amr_1blk_fc_cp_remote(mype, remote_pe, remote_block, idest,  
!!                             iopt, id, jd, kd, is, js, ks, ilays,   
!!                             jlays, klays, iface, nblk_ind, ipolar
!!  Call amr_1blk_fc_cp_remote(integer, integer, integer, integer,    
!!                             integer, integer, integer, integer,    
!!                             integer, integer, integer, integer,    
!!                             integer, integer, integer, integer,           
!!                             integer automatic array)
!!
!! ARGUMENTS
!!  
!!  Integer, intent(in) :: mype, remote_pe, remote_block
!!    mype           The local calling processor.
!!    remote_pe      The remote processor to fetch data from.
!!    remote_block   The remote block that data is fetched from.
!!  Integer, Intent(in) :: idest, iopt, id, jd, kd, is, js, ks
!!    idest           selects the storage space in data_1blk.fh which is to
!!                    be used in this call. If the leaf node is having its
!!                    guardcells filled then set this to 1, if its parent
!!                    is being filled set it to 2.
!!    iopt            a switch to control which data source is to be used
!!                    iopt=1 will use 'unk'
!!                    iopt>=2 will use 'work'
!!    id              lower limit of index range of points in x direction
!!                    on destination block
!!    jd              lower limit of index range of points in y direction
!!                    on destination block
!!    kd              lower limit of index range of points in z direction
!!                    on destination block
!!    is              lower limit of index range of points in x direction
!!                    on source block
!!    js              lower limit of index range of points in y direction
!!                    on source block
!!    ks              lower limit of index range of points in z direction
!!                    on source block
!!  Integer, Intent(in) :: ilays, jlays, klays, iface, nblk_ind
!!    ilays           no. of mesh points in x direction to be copied
!!    jlays           no. of mesh points in y direction to be copied
!!    klays           no. of mesh points in z direction to be copied
!!    iface           contains the block face on input. If this is
!!                    set between 1 to 6 the facevar variables on this
!                     face are averaged with those on a neighbor
!!                    at the same refinement level. If iface=0 this
!!                    averaging is not done.
!!    nblk_ind        index, running from 1-27 denoting location of 
!!                    neighbor block
!!  Integer, intent(in) :: ipolar(:)
!!    ipolar          switch used to apply special code at poles in 
!!                    spherical coordinates
!!
!! INCLUDES
!! 
!!  paramesh_preprocessor.fh
!!  mpif.h
!!
!! USES
!!
!!  paramesh_dimensions
!!  physicaldata
!!  tree
!!  workspace
!!  mpi_morton
!!  paramesh_interfaces
!!  paramesh_mpi_interfaces
!!  timings
!!
!! CALLS
!!
!!  amr_mpi_find_blk_in_buffer
!!  mpi_set_message_limits
!!
!! RETURNS
!!  
!!  Upon exit, data from a previously communicated block is copied into
!!  the appropriate on-processor buffer.
!!
!! DESCRIPTION
!!
!!   This routine copies guard cell information for cell face centered
!!   data to face iface in block
!!   idest, from the appropriate face of the neighboring block, assuming
!!   that the neighboring block is on a different processor.
!!   It can be easily edited to alter the data pattern required for schemes
!!   of different order.
!!
!! AUTHORS
!!
!!  Peter MacNeice, July 1998.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_fc_cp_remote(mype,remote_pe,remote_block,    & 
         idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,        & 
         ip2,jp2,kp2,iface,nblk_ind,ipolar)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits

      Implicit None

!-----Input/Output Variables
      Integer, Intent(in) :: mype,remote_pe,remote_block
      Integer, Intent(in) :: idest,id,jd,kd,is,js,ks
      Integer, Intent(in) :: ilays,jlays,klays
      Integer, Intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,iface
      Integer, Intent(in) :: nblk_ind
      Integer, Intent(in) :: ipolar(:)

!-----Local Arrays and Variables
      Real, Allocatable :: recvxf(:,:,:,:),                            & 
                           recvyf(:,:,:,:),                            & 
                           recvzf(:,:,:,:)

      Integer :: il,jl,kl,id1,jd1,kd1,is1,js1,ks1,js2
      Integer :: ilo,ihi,jlo,jhi,klo,khi
      Integer :: ill,jll,kll
      Integer :: index0, ii, jj, kk, i, j, k, jbface
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: ivar, ivar_next
      Integer :: nguard0, jstride
      Integer :: dtype
      Integer :: vtype

      Logical, Save :: first_call = .true.
      Logical :: lfound

!-----Begin Exectuable Code

      nguard0 = nguard*npgs

      If (force_consistency) Then

      If (no_permanent_guardcells) Then
       Allocate(recvxf(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,          & 
                       kl_bnd:ku_bnd))
       Allocate(recvyf(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,        & 
                       kl_bnd:ku_bnd))
       Allocate(recvzf(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,            & 
                       kl_bnd:ku_bnd+k3d))
      Else
       Allocate(recvxf(nbndvar,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd))
       Allocate(recvyf(nbndvar,il_bnd:iu_bnd,1:1+k2d,kl_bnd:ku_bnd))
       Allocate(recvzf(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:1+k3d))
      End If  ! End If (no_permanent_guardcells)

      recvxf = 0.
      recvyf = 0.
      recvzf = 0.

      End If  ! End If (force_consistency)


!-----Adjust index ranges
      il = ilays-1 
      jl = (jlays-1)*k2d 
      kl = (klays-1)*k3d 

      id1 = id + ip1
      jd1 = jd + jp1*k2d
      kd1 = kd + kp1*k3d
      is1 = is + ip1
      js1 = js + jp1*k2d
      ks1 = ks + kp1*k3d
      ja  = js
      jb  = js + jl 
      js2 = js1 + jl + jp2
      jstride = 1

      If (spherical_pm) Then
      If (lsingular_line) Then
      If (ipolar(1) == -1.and.jd <= nguard) Then
        jstride = -1
        If (no_permanent_guardcells) Then
          ja = nguard-jd+1
          jb = ja - jl
          js1 = nguard-(jd1-2)
        Else
          ja = 2*nguard
          jb = ja -nguard + k2d
          js1 = 2*(nguard+1) - jd1
          js2 = js1 + jl*jstride
        End If  ! End If (no_permanent_guardcells)
      Elseif(ipolar(2) == +1.and.jd > nyb+nguard) Then
        jstride = -1
        If (no_permanent_guardcells) Then
          ja = (nyb+1)-(jd-(nyb+nguard))
          jb = ja - jl
          js1 = (nyb+1)-(jd1-(nyb+nguard+1))
        Else
          ja = (nyb+nguard)-(jd-(nyb+nguard+1))
          jb = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          js1 = (nyb+nguard+1)-(jd1-(nyb+nguard+1))
          js2 = js1 + jl*jstride
        End If  ! End If (no_permanent_guardcells)
      End If  ! End If (ipolar(1) == -1.and.jd <= nguard)
      End If  ! End If (lsingular_line)
      End If  ! End If (spherical_pm)

      If (remote_block <= lnblocks .and. remote_pe == mype) Then

       If (no_permanent_guardcells) Then

       facevarx1(1:nbndvar,id1:id1+il+ip2,jd:jd+jl,kd:kd+kl,idest)     & 
          = gt_facevarx(1:nbndvar,is1:is1+il+ip2,                      & 
                                  ja:jb:jstride,                       & 
                                  ks:ks+kl,                            & 
                                  remote_block)

       If (ndim >= 2) Then
         facevary1(1:nbndvar,id:id+il,jd1:jd1+jl+jp2,kd:kd+kl,idest)   & 
            = gt_facevary(1:nbndvar,is:is+il,                          & 
                                    js1:js1+(jl+jp2)*jstride:jstride,  & 
                                    ks:ks+kl,                          & 
                                    remote_block)
       End If  ! End If (ndim >= 2)

       If (ndim == 3) Then
         facevarz1(1:nbndvar,id:id+il,jd:jd+jl,kd1:kd1+kl+kp2,idest)   & 
            = gt_facevarz(1:nbndvar,is:is+il,                          & 
                                    ja:jb:jstride,                     & 
                                    ks1:ks1+kl+kp2,                    & 
                                    remote_block)
       End If  ! End If (ndim == 3)

       If (force_consistency) Then
        recvxf(:,1,:,:) = gt_facevarx(:,1,:,:,remote_block)
        recvxf(:,2,:,:) = gt_facevarx(:,nxb+1,:,:,remote_block)
        If (ndim >= 2) Then
        recvyf(:,:,1,:) = gt_facevary(:,:,1,:,remote_block)
        recvyf(:,:,1+k2d,:) = gt_facevary(:,:,nyb+k2d,:,remote_block)
        End If  ! End If (ndim >= 2)
        If (ndim == 3) Then
        recvzf(:,:,:,1) = gt_facevarz(:,:,:,1,remote_block)
        recvzf(:,:,:,1+k3d) = gt_facevarz(:,:,:,nzb+k3d,remote_block)
        End If  ! End If (ndim == 3)
       End If  ! End If (force_consistency)

      Else !  If (no_permanent_guardcells)

       facevarx1(1:nbndvar,id1:id1+il+ip2,jd:jd+jl,kd:kd+kl,idest)     & 
          = facevarx(1:nbndvar,is1:is1+il+ip2,ja:jb:jstride,           & 
                               ks:ks+kl,                               & 
                               remote_block)

       If (ndim >= 2) Then
         facevary1(1:nbndvar,id:id+il,jd1:jd1+jl+jp2,kd:kd+kl,idest)   & 
            = facevary(1:nbndvar,is:is+il,                             & 
                                 js1:js2:jstride,                      & 
                                 ks:ks+kl,                             & 
                                 remote_block)
       End If  ! End If (ndim >= 2)

       If (ndim == 3) Then
         facevarz1(1:nbndvar,id:id+il,jd:jd+jl,kd1:kd1+kl+kp2,idest)   & 
           = facevarz(1:nbndvar,is:is+il,ja:jb:jstride,                &
                                ks1:ks1+kl+kp2,                        & 
                                remote_block)
       End If  ! If (ndim == 3)

       If (force_consistency) Then
        recvxf(:,1,:,:) = gt_facevarx(:,1,:,:,remote_block)
        recvxf(:,2,:,:) = gt_facevarx(:,2,:,:,remote_block)
        If (ndim >= 2) Then
        recvyf(:,:,1,:) = gt_facevary(:,:,1,:,remote_block)
        recvyf(:,:,1+k2d,:) = gt_facevary(:,:,1+k2d,:,remote_block)
        End If  ! End If (ndim >= 2)
        If (ndim == 3) Then
        recvzf(:,:,:,1) = gt_facevarz(:,:,:,1,remote_block)
        recvzf(:,:,:,1+k3d) = gt_facevarz(:,:,:,1+k3d,remote_block)
        End If  ! End If (ndim == 3)
       End If  ! End If (force_consistency)

      End If !  End If (no_permanent_guardcells)

!-----The block is remote
      Else  ! If (remote_block <= lnblocks .and. remote_pe == mype)

        call amr_mpi_find_blk_in_buffer(mype,remote_block,             & 
                       remote_pe,idest,dtype,index0,lfound)

!-------If this routine is executing a copy to fill guardcells of a
!-------leaf blocks^s parent, and the remote block is not found, then
!-------it is assumed that it is not in the list of buffered remote blocks
!-------because it is not really needed. Therefore in this case we
!-------return without copying anything.
        If (idest == 2 .and. (.not.lfound)) Then
         If (force_consistency) Then
          Deallocate(recvxf)
          Deallocate(recvyf)
          Deallocate(recvzf)
         End If  ! End If (force_consistency)
         Return
        End If  ! End If (idest == 2 .and. (.not.lfound))

!-------starting index if cell-centered data is also included in recv_buf
        If (l_datapacked(2)) index0 = index0 +                         &
                             ngcell_on_cc*message_size_cc(dtype)

        ill = ilays
        jll = jlays
        kll = klays

        vtype = 2
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,                    & 
                     ill,jll,kll)

        ilo = is1
        ihi = is1+il+ip2
        If (iface == 1) ihi = nxb+1+nguard0
        If (iface == 2) ilo = 1+nguard0


        kk = kd
        Do k = ka,kb

        jj = jd
        jstride = 1
        js2 = js
        js1 = js+jl

        If (spherical_pm) Then
        If (lsingular_line) Then
        If (ipolar(1) == -1 .and. jd <= nguard) then
        jstride = -1
        If (no_permanent_guardcells) Then
          js2 = jd + jl + (nguard - 2*(jd+jl)) +1
          js1 = jd      + (nguard - 2* jd    ) +1
          jj  = jd + jl
        Else
          js2 = jd + jl + 2*(nguard - (jd+jl)) +1
          js1 = jd      + 2*(nguard -  jd    ) +1
          jj  = jd + jl
        End If  ! End If (no_permanent_guardcells)
        Elseif (ipolar(2) == +1.and.jd > nyb+nguard) Then
        jstride = -1
        If (no_permanent_guardcells) Then
          js1 = nyb - ( jd    -(nyb+nguard+1))
          js2 = nyb - ((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        Else
          js1 = (nyb+nguard)-( jd    -(nyb+nguard+1))
          js2 = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        End If  ! End If (no_permanent_guardcells)
        End If  ! End If (ipolar(1) == -1 .and. jd <= nguard)
        End If  ! End If (lsingular_line)
        End If  ! End If (spherical_pm)

        Do j = ja,jb
        ii = id1
        If (iface == 2) ii = id1-1
        Do i = ia,ib
          If (k >= ks .and. k <= ks + kl) Then
          If (j >= js2 .and. j <= js1) Then
          If (i >= ilo .and. i <= ihi) Then

          Do ivar=1,ngcell_on_fc(1)

          ivar_next = gcell_on_fc_pointer(1,ivar)

          If (i == 1+nguard0 .and. iface == 2) Then
            If (force_consistency) Then
             recvxf(ivar_next,1,j,k) = temprecv_buf(index0+ivar)
            End If  ! End If (force_consistency)
          Else If (i == nxb+1+nguard0 .and. iface == 1) Then
            If (force_consistency) Then 
             recvxf(ivar_next,2,j,k) = temprecv_buf(index0+ivar)
            End If  ! End If (force_consistency)
          Else
            facevarx1(ivar_next,ii,jj,kk,idest) =                      &
                           temprecv_buf(index0+ivar)
          End If  ! End If (If(i == 1+nguard0 .and. iface == 2)

          End Do  ! End Do ivar=1,ngcell_on_fc(1)

          End If  ! End If (i >= ilo .and. i <= ihi)
          End If  ! End If (j >= js2 .and. j <= js1)
          End If  ! End If (k >= ks .and. k <= ks + kl)

          if (i >= ilo .and. i <= ihi) ii = ii + 1
          index0 = index0+ngcell_on_fc(1)

        End Do  ! End Do i = ia,ib
        If (j >= js2 .and. j <= js1) jj = jj + jstride
        End Do  ! End Do j = ja,jb
        If (k >= ks .and. k <= ks + kl) kk = kk + 1
        End Do  ! End Do k = ka,kb

        if(ndim >= 2) Then

        vtype = 3 
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,                    & 
                     ill,jll,kll)


! reset js1, js2 to values set on entry
        js1 = js + jp1*k2d
        js2 = js1 + jl + jp2

        jlo = js1
        jhi = js1+jl+jp2
        If (iface == 3) jhi = nyb+1+nguard0*k2d
        If (iface == 4) jlo = 1 + nguard0*k2d

        kk = kd
        Do k = ka,kb

        jj = jd1
        If (iface == 4) jj = jd1-1
        jstride = 1

        If (spherical_pm) Then
        If (lsingular_line) Then
        If (ipolar(1) == -1 .and. jd <= nguard) Then
        jstride = -1
        If (no_permanent_guardcells) Then 
          jlo = jd + jl + (nguard+2 - 2*(jd+jl))
          jhi = jd      + (nguard+2 - 2* jd    )
          jj  = jd + jl
        Else
          jlo = jd + jl + 2*(nguard+1 - (jd+jl))
          jhi = jd      + 2*(nguard+1 -  jd    ) 
          jj  = jd + jl
        End If  ! End If (no_permanent_guardcells)
        Elseif(ipolar(2) == +1.and.jd > nyb+nguard) Then
        jstride = -1
        If (no_permanent_guardcells) Then
          jhi = nyb - ( jd    -(nyb+nguard+1))
          jlo = nyb - ((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl + 1
        Else
          jhi = (nyb+nguard)-( jd    -(nyb+nguard+1))
          jlo = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl + 1
        End If  ! End If (no_permanent_guardcells)
        End If  ! End If (ipolar(1) == -1 .and. jd <= nguard)
        End If  ! End If (lsingular_line)
        End If  ! End If (spherical_pm)

        Do j = ja,jb
        ii = id
        Do i = ia,ib
          If (k >= ks .and. k <= ks + kl) Then
          If (j >= jlo .and. j <= jhi) Then
          If (i >= is .and. i <= is + il) Then

          Do ivar=1,ngcell_on_fc(2)

          ivar_next = gcell_on_fc_pointer(2,ivar)

          If (j == 1+nguard0*k2d .and. iface == 4) then
            If (force_consistency) Then
            recvyf(ivar_next,i,1,k) = temprecv_buf(index0+ivar)
            End If  ! End If (force_consistency)
          Elseif (j == nyb+1+nguard0*k2d .and. iface == 3) Then
             If (force_consistency) Then
             recvyf(ivar_next,i,1+k2d,k) = temprecv_buf(index0+ivar)
            End If  ! End (force_consistency)
          Else
            facevary1(ivar_next,ii,jj,kk,idest) =                      & 
                                    temprecv_buf(index0+ivar)
          End If  ! End If (j == 1+nguard0*k2d .and. iface == 4)

           End Do  ! End Do ivar=1,ngcell_on_fc(2)

           End If  ! End If (i >= is .and. i <= is + il)
           End If  ! End If (j >= jlo .and. j <= jhi)
           End If  ! End If (k >= ks .and. k <= ks + kl)
           If (i >= is .and. i <= is + il) ii = ii + 1
           index0 = index0+ngcell_on_fc(2)
        End Do  ! End Do i = ia,ib
        If (j >= jlo .and. j <= jhi) jj = jj + jstride
        End Do  ! End Do j = ja,jb
        if (k >= ks .and. k <= ks + kl) kk = kk + 1
        End Do  ! End Do k = ka,kb

        End If  ! End If (ndim >= 2)

        if(ndim == 3) then

        vtype = 4 
        Call mpi_set_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,                    & 
                     ill,jll,kll)

        klo = ks1
        khi = ks1+kl+kp2
        If (iface == 5) khi = nzb+1+nguard0*k3d
        If (iface == 6) klo = 1+nguard0*npgs

        kk = kd1
        If(iface == 6) kk = kd1-1
        Do k = ka,kb

        jj = jd
        jstride = 1

        js2 = js
        js1 = js+jl

        If (spherical_pm) Then
        If (lsingular_line) Then
        If (ipolar(1) == -1 .and. jd <= nguard) Then
        jstride = -1
        If (no_permanent_guardcells) Then
          js2 = jd + jl + (nguard - 2*(jd+jl)) +1
          js1 = jd      + (nguard - 2* jd    ) +1
          jj  = jd + jl
        Else
          js2 = jd + jl + 2*(nguard - (jd+jl)) +1
          js1 = jd      + 2*(nguard -  jd    ) +1
          jj  = jd + jl
        End If  ! End If (no_permanent_guardcells)
        Elseif(ipolar(2) == +1 .and. jd > nyb+nguard) Then
        jstride = -1
        If (no_permanent_guardcells) Then
          js1 = nyb - ( jd    -(nyb+nguard+1))
          js2 = nyb - ((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        Else
          js1 = (nyb+nguard)-( jd    -(nyb+nguard+1))
          js2 = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        End If  ! End If (no_permanent_guardcells)
        End If  ! End If (ipolar(1) == -1 .and. jd <= nguard)
        End If  ! End If (lsingular_line)
        End If  ! End If (spherical_pm)

        Do j = ja,jb
        ii = id
        Do i = ia,ib
          If (k >= klo .and. k <= khi) Then
          If (j >= js2 .and. j <= js1) Then
          If (i >= is .and. i <= is + il) Then

          Do ivar=1,ngcell_on_fc(3)
          ivar_next = gcell_on_fc_pointer(3,ivar)

          If (k == 1+nguard0*k3d .and. iface == 6) Then
            If (force_consistency) Then
             recvzf(ivar_next,i,j,1) = temprecv_buf(index0+ivar)
            End If  ! End If (force_consistency)
          Elseif(k == nzb+1+nguard0*k3d .and. iface == 5) Then
            If (force_consistency) Then
             recvzf(ivar_next,i,j,1+k3d) = temprecv_buf(index0+ivar)
            End If  ! End If (force_consistency)
          Else
            facevarz1(ivar_next,ii,jj,kk,idest) =                    &
                       temprecv_buf(index0+ivar)
          End If  ! End If (k == 1+nguard0*k3d .and. iface == 6)

          End Do  ! End do ivar=1,ngcell_on_fc(3)

          End If  ! End If (i >= is .and. i <= is + il)
          End If  ! End If (j >= js2 .and. j <= js1)
          End If  ! End If (k >= klo .and. k <= khi)

          if (i >= is .and. i <= is + il) ii = ii + 1
          index0 = index0+ngcell_on_fc(3)

        End Do  ! End Do i = ia,ib
        If (j >= js2 .and. j <= js1) jj = jj + jstride
        End Do  ! End Do j = ja,jb
        If (k >= klo .and. k <= khi) kk = kk + 1
        End Do  ! End Do k = ka,kb

        End If  ! End If (ndim == 3)


      End If  ! End If (remote_block <= lnblocks .and. remote_pe == mype)

!-----make sure srl shared faces end up with the mean value from each
!-----block^s data 
      If (force_consistency) Then

       If (iface == 1) Then

       Do ivar=1,ngcell_on_fc(1)
         ivar_next = gcell_on_fc_pointer(1,ivar)

         facevarx1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d,     & 
                             1+nguard*k3d:nzb+nguard*k3d,idest) =      & 
         .5*(                                                          & 
         facevarx1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d,     & 
                             1+nguard*k3d:nzb+nguard*k3d,idest) +      & 
           recvxf(ivar_next,2,1+nguard0*k2d:nyb+nguard0*k2d,           & 
                              1+nguard0*k3d:nzb+nguard0*k3d) ) 
       End Do  ! End Do ivar=1,ngcell_on_fc(1)

       Elseif (iface == 2) Then

       Do ivar=1,ngcell_on_fc(1)
         ivar_next = gcell_on_fc_pointer(1,ivar)

         facevarx1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
                            1+nguard*k3d:nzb+nguard*k3d,idest) =       & 
         .5*(                                                          & 
           recvxf(ivar_next,1,1+nguard0*k2d:nyb+nguard0*k2d,           & 
                              1+nguard0*k3d:nzb+nguard0*k3d) +         & 
         facevarx1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
                             1+nguard*k3d:nzb+nguard*k3d,idest) )      
       End Do  ! End Do ivar=1,ngcell_on_fc(1)

       Elseif (iface == 3) Then

       jbface = 2
       If (spherical_pm) then
         If (ipolar(1) == -1) jbface = 1
       End If  ! End If (spherical_pm)
       Do ivar=1,ngcell_on_fc(2)
         ivar_next = gcell_on_fc_pointer(2,ivar)

         facevary1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d,         & 
                             1+nguard*k3d:nzb+nguard*k3d,idest) =      & 
         .5*(                                                          & 
         facevary1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d,         & 
                             1+nguard*k3d:nzb+nguard*k3d,idest) +      & 
           recvyf(ivar_next,1+nguard0:nxb+nguard0,jbface,              & 
                           1+nguard0*k3d:nzb+nguard0*k3d) )
       End Do  ! End Do ivar=1,ngcell_on_fc(2)

       Elseif (iface == 4) Then

       jbface = 1
       If (spherical_pm) then
         If(ipolar(2) == 1) jbface = 2
       End If  ! End If (spherical_pm)
       Do ivar=1,ngcell_on_fc(2)
         ivar_next = gcell_on_fc_pointer(2,ivar)

         facevary1(ivar_next,1+nguard:nxb+nguard,nyb+k2d+nguard*k2d,   & 
                             1+nguard*k3d:nzb+nguard*k3d,idest) =      & 
         .5*(                                                          & 
           recvyf(ivar_next,1+nguard0:nxb+nguard0,jbface,              & 
                           1+nguard0*k3d:nzb+nguard0*k3d)  +           & 
           facevary1(ivar_next,1+nguard:nxb+nguard,nyb+k2d+nguard*k2d, & 
                             1+nguard*k3d:nzb+nguard*k3d,idest) ) 
       End Do  ! End Do ivar=1,ngcell_on_fc(2)

       Elseif (iface == 5) Then

       Do ivar=1,ngcell_on_fc(3)
         ivar_next = gcell_on_fc_pointer(3,ivar)

         facevarz1(ivar_next,1+nguard:nxb+nguard,                      & 
                   1+nguard*k2d:nyb+nguard*k2d,1+nguard*k3d,idest) =   & 
         .5*(                                                          & 
         facevarz1(ivar_next,1+nguard:nxb+nguard,                      & 
                   1+nguard*k2d:nyb+nguard*k2d,1+nguard*k3d,idest) +   & 
           recvzf(ivar_next,1+nguard0:nxb+nguard0,                     & 
                           1+nguard0*k2d:nyb+nguard0*k2d,1+k3d) )
       End Do  ! End Do ivar=1,ngcell_on_fc(3)

       Elseif(iface == 6) Then

       Do ivar=1,ngcell_on_fc(3)
         ivar_next = gcell_on_fc_pointer(3,ivar)

         facevarz1(ivar_next,1+nguard:nxb+nguard,                      & 
             1+nguard*k2d:nyb+nguard*k2d,nzb+k3d+nguard*k3d,idest) =   & 
         .5*(                                                          & 
           recvzf(ivar_next,1+nguard0:nxb+nguard0,                     &  
                           1+nguard0*k2d:nyb+nguard0*k2d,1) +          & 
           facevarz1(ivar_next,1+nguard:nxb+nguard,                    & 
             1+nguard*k2d:nyb+nguard*k2d,nzb+k3d+nguard*k3d,idest) )   

       End Do  ! End Do ivar=1,ngcell_on_fc(3)
       End If  ! End If (iface == 1)

      End If  ! End If (force_consistency)

      If (force_consistency) Then
       Deallocate(recvxf)
       Deallocate(recvyf)
       Deallocate(recvzf)
      End If  ! End If (force_consistency)

      Return
      End Subroutine amr_1blk_fc_cp_remote
