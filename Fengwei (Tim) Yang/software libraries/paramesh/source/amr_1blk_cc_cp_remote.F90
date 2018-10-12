!#define DEBUG

!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_cc_cp_remote
!! NAME
!!
!!   amr_1blk_cc_cp_remote
!!
!! SYNOPSIS
!! 
!!  Call amr_1blk_cc_cp_remote(mype, remote_pe, remote_block, idest,  &
!!                             iopt, id, jd, kd, is, js, ks, ilays,   &
!!                             jlays, klays, nblk_ind, ipolar)
!!  Call amr_1blk_cc_cp_remote(integer, integer, integer, integer,    &
!!                             integer, integer, integer, integer,    &
!!                             integer, integer, integer, integer,    &
!!                             integer, integer, integer,             &
!!                             integer automatic array)
!!
!! ARGUMENTS
!!  
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
!!  Integer, Intent(in) :: ilays, jlays, klays, nblk_ind
!!    ilays           no. of mesh points in x direction to be copied
!!    jlays           no. of mesh points in y direction to be copied
!!    klays           no. of mesh points in z direction to be copied
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
!! This routine copies guard cell information to face iface in layer
!! idest of the working block, from the appropriate face of the neighboring 
!! block, assuming that the neighboring block is on a different processor.
!!
!! AUTHORS
!!
!!  Peter MacNeice, July 1998.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      Subroutine amr_1blk_cc_cp_remote(mype,remote_pe,remote_block,  & 
                                       idest,iopt,id,jd,kd,is,js,ks, &
                                       ilays,jlays,klays,            & 
                                       nblk_ind,ipolar)


!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use mpi_morton
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits
      Use timings

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, Intent(in) :: mype,remote_pe,remote_block
      Integer, Intent(in) :: idest,iopt,id,jd,kd,is,js,ks
      Integer, Intent(in) :: ilays,jlays,klays,nblk_ind
      Integer, Intent(in) :: ipolar(:)


!-----Local variables and arrays
      Integer :: il,jl,kl,iopt0
      Integer :: ill,jll,kll
      Integer :: ia, ib, ja, jb, ka, kb, jstride, js0, jsl
      Integer :: js1,js2
      Integer :: i, j, k, ii, jj, kk, ivar, ivar_next
      Integer :: indx
      Double Precision :: time1
      Integer :: ierr,dtype
      Integer :: vtype
      Logical :: lfound

!-----Begin Executable code.

!-----Adjust index ranges

      il = ilays-1
      jl = (jlays-1)*k2d
      kl = (klays-1)*k3d

      If (remote_block <= lnblocks .and. remote_pe == mype) Then

      If (timing_mpix) Then
         time1 = mpi_wtime()
      End If

      jstride = 1

!-------------------------
!-----Compute indeces for unk1, facevarx1(yz), unk_e_x1(yz), unk_n1
      If (iopt == 1) Then
!-------------------------

         ia = is
         ib = is+il
         ja = js
         jb = js+jl
         ka = ks
         kb = ks+kl

!-------Handle poles in spherical coordinates
        If(spherical_pm) Then
        If(lsingular_line) Then
        If(ipolar(1) == -1 .and. jd <= nguard) Then
          jstride = -1
          If (no_permanent_guardcells) Then
            ja = nguard-jd+1
            jb = ja - jl
          Else
            ja = 2*nguard     
            jb = ja -nguard     
          Endif  ! End If (no_pemanent_guardcells)
        Elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) Then
          jstride = -1
          if (no_permanent_guardcells) Then
            ja = (nyb+1)-(jd-(nyb+nguard))
            jb = ja - jl
          Else
            ja = (nyb+nguard)-(jd-(nyb+nguard+1))
            jb = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          End If  ! End If (no_permanent_guardcells)
        End If  ! End If (ipolar(1)
        End If  ! End If (lsingular_line)
        End If  ! End If (spherical_pm)

! Copy complete remote block into a buffer block called recv.
        If (no_permanent_guardcells) Then      

        Do ivar=1, nvar
          If (int_gcell_on_cc(ivar)) Then
            unk1(ivar,id:id+il,jd:jd+jl,kd:kd+kl,idest) = & 
              gt_unk(ivar,ia:ib,ja:jb:jstride,ka:kb,remote_block)
          End If  ! End If (int_gcell_on_cc(ivar))
        End Do  ! End Do ivar = 1, nvar

        Else ! If (no_permanent_guardcells)

        Do ivar=1, nvar
          If (int_gcell_on_cc(ivar)) Then
            unk1(ivar,id:id+il,jd:jd+jl,kd:kd+kl,idest) = & 
              unk(ivar,ia:ib,ja:jb:jstride,ka:kb,remote_block)
          End If  ! End If (int_gcell_on_cc(ivar))
        End Do  ! End Do ivar = 1, nvar

        End If ! End If (no_permanent_guardcells)

!-------------------------
!-----Compute indeces for work1 arrays
      Else If (iopt.ge.2) Then
!-------------------------
        jstride = 1
        js0 = js
        jsl = js0 + jl

!-------Handle poles in spherical coordinates
        If (spherical_pm) Then
        If (lsingular_line) Then
        If (ipolar(1) == -1 .and. jd <= nguard_work) Then
          jstride = -1
          If (no_permanent_guardcells) Then
            js0 = nguard_work-jd+1
            jsl = js0 - jl
          Else
            js0 = 2*nguard_work
            jsl = js0 -nguard_work
          End If  ! End If (no_permanent_guardcells)
        Elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) Then
          jstride = -1
          If (no_permanent_guardcells) Then
            js0 = (nyb+1)-(jd-(nyb+nguard_work))
            jsl = js0 - jl
          Else
            js0 = (nyb+nguard_work)-(jd-(nyb+nguard_work+1))
            jsl = (nyb+nguard_work)-((jd+jl)-(nyb+nguard_work+1))
          End If  ! End If (no_permanent_guardcells)
        End If  ! End If (ipolar == -1 .and. jd <= nguard_work)
        End If  ! End If (lsingular_line)
        End If  ! End If (spherical_pm)

        iopt0 = iopt-1
! Copy complete remote block into a buffer block called recvw.

        work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) =                      & 
           work(is:is+il,js0:jsl:jstride,ks:ks+kl,remote_block,iopt0)
!-------------------------
      End If  ! End If (iopt == 1)
!-------------------------

      If (timing_mpix) Then
         timer_amr_1blk_cc_cp_remote(1) =                              & 
           timer_amr_1blk_cc_cp_remote(1) + mpi_wtime() - time1
      Else
         timer_amr_1blk_cc_cp_remote(1) = -999.
      End If  ! End If (timing_mpix)

!-----Start code to fetch block data if off-processor
      Else ! If (remote_block <= lnblocks .and. remote_pe == mype)

        If (timing_mpix) Then
           time1 = mpi_wtime()
        End If  ! End If (timing_mpix)

        Call amr_mpi_find_blk_in_buffer(mype,remote_block,             & 
                              remote_pe,idest,dtype,indx,lfound)

        If (timing_mpix) Then
         timer_amr_1blk_cc_cp_remote(2) =                             & 
           timer_amr_1blk_cc_cp_remote(2) + mpi_wtime() - time1
         time1 = mpi_wtime()
        Else
         timer_amr_1blk_cc_cp_remote(2) = -999.
        End If  ! End If (timing_mpix)

!-------If this routine is executing a copy to fill guardcells of a
!-------leaf blocks^s parent, and the remote block is not found, then
!-------it is assumed that it is not in the list of buffered remote blocks
!-------because it is not really needed. Therefore in this case we
!-------return without copying anything.
        If (idest == 2 .and. (.Not.lfound)) Then
           Return
        End If  ! End If (idest == 2 .and. (.Not.lfound))


        If (iopt == 1) Then

          vtype = 1
          ill = ilays
          jll = jlays
          kll = klays

          Call mpi_set_message_limits(                 & 
                       dtype,ia,ib,ja,jb,ka,kb,vtype,  & 
                       ill,jll,kll)

          kk = kd
          Do k = ka,kb
           jj = jd
           jstride = 1
           js2 = js
           js1 = js+jl

!---------Handle poles in spherical coordinates
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
          ElseIf(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
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
            If (k >= ks .and. k <= ks + kl) Then
            If (j >= js2.and. j <= js1)     Then
            If (i >= is .and. i <= is + il) Then

            Do ivar = 1, ngcell_on_cc
              ivar_next = gcell_on_cc_pointer(ivar)
              unk1(ivar_next,ii,jj,kk,idest) =                         & 
                       temprecv_buf(indx+ivar)
            End Do  !  End Do ivar=1,ngcell_on_cc

            End If  ! End If (i >= is .and. i <= is + il)
            End If  ! End If (j >= js2.and. j <= js1)
            End If  ! End If (k >= ks .and. k <= ks + kl)
            if (i >= is .and. i <= is + il) ii = ii + 1
            indx = indx+ngcell_on_cc

          End Do  !  End Do i = ia,ib

          If (j >= js2 .and. j <= js1) jj = jj + jstride

          End Do  !  End Do j = ja,jb

          If (k >= ks .and. k <= ks + kl) kk = kk + 1

          End Do  !  End Do k = ka,kb

        Elseif (iopt > 1) Then

          vtype = 0
          ill = ilays
          jll = jlays
          kll = klays

          Call mpi_set_message_limits(                  & 
                       dtype,ia,ib,ja,jb,ka,kb,vtype,   & 
                       ill,jll,kll)
          kk = kd
          Do k = ka,kb

          jj = jd
          jstride = 1
          js2 = js
          js1 = js+jl

          If (spherical_pm) Then
          If (lsingular_line) Then
          If (ipolar(1).eq.-1.and.jd.le.nguard_work) Then
            jstride = -1
            If (no_permanent_guardcells) Then
              js2 = jd + jl + (nguard_work - 2*(jd+jl)) +1
              js1 = jd      + (nguard_work - 2* jd    ) +1
              jj  = jd + jl
            Else
              js2 = jd + jl + 2*(nguard_work - (jd+jl)) +1
              js1 = jd      + 2*(nguard_work -  jd    ) +1
              jj  = jd + jl
            End If  ! End If (no_permanent_guardcells)
          Elseif(ipolar(2) == +1 .and. jd > nyb+nguard_work) Then
            jstride = -1
            If (no_permanent_guardcells) Then
              js1 = nyb - ( jd    -(nyb+nguard_work+1))
              js2 = nyb - ((jd+jl)-(nyb+nguard_work+1))
              jj  = jd + jl
            Else
              js1 = (nyb+nguard_work)-( jd    -(nyb+nguard_work+1))
              js2 = (nyb+nguard_work)-((jd+jl)-(nyb+nguard_work+1))
              jj  = jd + jl
            End If  ! End If (no_permanent_guardcells)
          End If  ! End If (ipolar(1).eq.-1.and.jd.le.nguard_work)
          End If  ! End If (lsingular_line)
          End If  ! End If (spherical_pm)

          Do j = ja,jb
          ii = id

          Do i = ia,ib
            If (k >= ks .and. k <= ks + kl) Then
            If (j >= js2.and. j <= js1)     Then
            If (i >= is .and. i <= is + il) Then
               work1(ii,jj,kk,idest) =                                 & 
                       temprecv_buf(indx+1)
            End If  ! End If (i >= is .and. i <= is + il)
            End If  ! End If (j >= js2.and. j <= js1) 
            End If  ! End If (k >= ks .and. k <= ks + kl)
            If (i >= is .and. i <= is + il) ii = ii + 1
            indx = indx+1
          End Do  ! End Do i = ia,ib
          If (j >= js2 .and. j <= js1) jj = jj + jstride
          End Do  ! End Do j = ja,jb
          If (k >= ks .and. k <= ks + kl) kk = kk + 1
          End Do  ! End Do I = Ka,Kb

        End If  ! End If (iopt == 1)
      
        If (timing_mpix) Then
         timer_amr_1blk_cc_cp_remote(3) =                              & 
           timer_amr_1blk_cc_cp_remote(3) + mpi_wtime() - time1
        Else
         timer_amr_1blk_cc_cp_remote(3) = -999.
        End If

      End If  ! End If (remote_block <= lnblocks .and. remote_pe == mype)

      Return
      End Subroutine amr_1blk_cc_cp_remote





