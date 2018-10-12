!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source amr_perm_to_1blk
!! NAME
!!
!!   amr_perm_to_1blk
!!
!! SYNOPSIS
!!
!!   Call amr_perm_to_1blk(lcc,lfc,lec,lnc,lb,pe,iopt,idest)
!!   Call amr_perm_to_1blk(logical,logical,logical,logical,
!!                         integer,integer,integer,integer)
!!
!! ARGUMENTS
!!
!!   Logical, Intent(in) :: lcc   copies cell centered data if true
!!   Logical, Intent(in) :: lfc   copies cell face-centered data if true
!!   Logical, Intent(in) :: lec   copies cell edge-centered data if true
!!   Logical, Intent(in) :: lnc   copies cell corner data if true
!!   Integer, Intent(in) :: lb    block from which data is to be copied
!!   Integer, Intent(in) :: pe    processor from which data is to be copied
!!   Integer, Intent(in) :: iopt  data structure to be copied
!!   Integer, Intent(in) :: idest sets value for last dimension index
!!                                  in the 1-blk data arrays!!   
!! INCLUDES
!!
!!   paramesh_preprocessorh.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   workspace
!!   mpi_morton
!!   paramesh_mpi_interfaces
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   mpi_set_message_limits
!!   amr_mpi_find_blk_in_buffer
!!
!! RETURNS
!!
!!   Nothing Returned.  
!!   
!! DESCRIPTION
!!
!!   This routine copies data to the 1-block working arrays with guardcells
!!   from the permanent data arrays, which may or may not have permanent
!!   guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
!!   in physicaldata.fh.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          February 1999
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_perm_to_1blk( lcc,lfc,lec,lnc,lb,pe,iopt,idest)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use mpi_morton
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguemtns
      Integer, intent(in) ::  lb,pe,iopt,idest
      Logical, intent(in) ::  lcc,lfc,lec,lnc


!-----Local variables and arrays.
      Integer :: iopt0
      Integer :: nguard0,nguard_work0
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: i, j, k, ivar, ivar_next
      Integer :: vtype,dtype,rem_blk,rem_pe,mype,ierr
      Integer :: index,index0
      Logical :: lfound

!-----Begin executable code.

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      If (lb > lnblocks) Then
           Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
           rem_blk = lb
           rem_pe  = mype
           Call amr_mpi_find_blk_in_buffer(mype,rem_blk,               & 
                             rem_pe,idest,dtype,index0,lfound)
           If (.Not.lfound) Then 
             Write(*,*)                                                & 
              'perm to 1blk reporting blk not found',                  & 
              ' mype=',mype,' looked for ',lb,                         & 
              ' where lnblocks=',lnblocks,                             & 
              ' strt_buffer=',strt_buffer,' last_buffer=',last_buffer, & 
              ' laddress ',laddress(:,strt_buffer:last_buffer)
             Call amr_abort()    ! remove this abort after testing
           End If
      End If  ! End If (lb > lnblocks)

!-----Put block lb's data into the data_1blk.fh datastructures, with the
!-----appropriate guardcell padding.
      If (iopt == 1) Then

        If (lcc) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nvar
              If (int_gcell_on_cc(ivar)) Then
                 unk1(ivar,                                            &
                      1+nguard:nxb+nguard,                             &
                      1+nguard*k2d:nyb+nguard*k2d,                     & 
                      1+nguard*k3d:nzb+nguard*k3d,                     &
                      idest) =                                         & 
                 unk(ivar,                                             &
                     1+nguard0:nxb+nguard0,                            &
                     1+nguard0*k2d:nyb+nguard0*k2d,                    & 
                     1+nguard0*k3d:nzb+nguard0*k3d,                    &
                     lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

             vtype = 1
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)
             index = index0

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_cc
                 ivar_next = gcell_on_cc_pointer(ivar)
                 unk1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
               End Do
               index = index + ngcell_on_cc
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

        End If  ! End If  (lcc)

        If (lfc) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nfacevar
              If (int_gcell_on_fc(1,ivar)) Then
                 facevarx1(ivar,                                       &
                           1+nguard:nxb+nguard+1,                      & 
                           1+nguard*k2d:nyb+nguard*k2d,                & 
                           1+nguard*k3d:nzb+nguard*k3d,                &
                           idest) =                                    & 
                 facevarx(ivar,                                        &
                          1+nguard0:nxb+nguard0+1,                     & 
                          1+nguard0*k2d:nyb+nguard0*k2d,               & 
                          1+nguard0*k3d:nzb+nguard0*k3d,lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

!------------starting index if cell-centered data is also included in recv_buf
             index = index0 + ngcell_on_cc*message_size_cc(dtype)

             vtype = 2
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_fc(1)
                 ivar_next = gcell_on_fc_pointer(1,ivar)
                 facevarx1(ivar_next,i,j,k,idest) =                    & 
                          temprecv_buf(index+ivar)
               End Do
               index = index+ngcell_on_fc(1)
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

          If (ndim >= 2) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nfacevar
              If (int_gcell_on_fc(2,ivar)) Then
                 facevary1(ivar,                                       &
                           1+nguard:nxb+nguard,                        & 
                           1+nguard*k2d:nyb+(nguard+1)*k2d,            & 
                           1+nguard*k3d:nzb+nguard*k3d,                &
                           idest) =                                    & 
                 facevary(ivar,                                        &
                          1+nguard0:nxb+nguard0,                       & 
                          1+nguard0*k2d:nyb+(nguard0+1)*k2d,           & 
                          1+nguard0*k3d:nzb+nguard0*k3d,               &
                          lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

             vtype = 3
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_fc(2)
                 ivar_next = gcell_on_fc_pointer(2,ivar)
                 facevary1(ivar_next,i,j,k,idest) =                    & 
                          temprecv_buf(index+ivar)
               End Do
               index = index + ngcell_on_fc(2)
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

          End If  ! end If (ndim >= 2)

          If (ndim == 3) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nfacevar
              If (int_gcell_on_fc(3,ivar)) Then
                 facevarz1(ivar,                                       &
                           1+nguard:nxb+nguard,                        & 
                           1+nguard*k2d:nyb+nguard*k2d,                & 
                           1+nguard*k3d:nzb+(nguard+1)*k3d,            &
                           idest) =                                    & 
                  facevarz(ivar,                                       &
                           1+nguard0:nxb+nguard0,                      & 
                           1+nguard0*k2d:nyb+nguard0*k2d,              & 
                           1+nguard0*k3d:nzb+(nguard0+1)*k3d,          &
                           lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

             vtype = 4
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_fc(3)
                 ivar_next = gcell_on_fc_pointer(3,ivar)
                 facevarz1(ivar_next,i,j,k,idest) =                    & 
                          temprecv_buf(index+ivar)
               End Do
               index = index + ngcell_on_fc(3)
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

          End If  ! End If (ndim == 3)

        End If  ! End If (lfc)

        If (ndim > 1) Then
          If (lec) Then

            If (lb <= lnblocks) Then

            Do ivar=1,nvaredge
              If (int_gcell_on_ec(1,ivar)) Then
                 unk_e_x1(ivar,                                        &
                          1+nguard:nxb+nguard,                         & 
                          1+nguard*k2d:nyb+nguard*k2d+k2d,             & 
                          1+nguard*k3d:nzb+nguard*k3d+k3d,             &
                          idest) =                                     & 
                 unk_e_x(ivar,                                         &
                         1+nguard0:nxb+nguard0,                        & 
                         1+nguard0*k2d:nyb+nguard0*k2d+k2d,            & 
                         1+nguard0*k3d:nzb+nguard0*k3d+k3d,            &
                         lb)
              End If
            End Do

            ElseIf (lb > lnblocks) Then

!-----------starting index if cell-centered data is also included in recv_buf
            index = index0 + ngcell_on_cc*message_size_cc(dtype)       & 
                           + maxval(ngcell_on_fc(1:ndim))*             & 
                                   message_size_fc(dtype)

            vtype = 5
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)
  
            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_ec(1)
                ivar_next = gcell_on_ec_pointer(1,ivar)
                unk_e_x1(ivar_next,i,j,k,idest) =                      & 
                         temprecv_buf(index+ivar)
              End Do
              index = index + ngcell_on_ec(1)
            End Do
            End Do
            End Do

            End If  ! End If (lb <= lnblocks)

            If (lb <= lnblocks) Then

            Do ivar=1,nvaredge
              If (int_gcell_on_ec(2,ivar)) Then
                 unk_e_y1(ivar,                                        &
                          1+nguard:nxb+nguard+1,                       &  
                          1+nguard*k2d:nyb+nguard*k2d,                 & 
                          1+nguard*k3d:nzb+(nguard+1)*k3d,             &
                          idest) =                                     & 
                 unk_e_y(ivar,                                         &
                         1+nguard0:nxb+nguard0+1,                      & 
                         1+nguard0*k2d:nyb+nguard0*k2d,                & 
                         1+nguard0*k3d:nzb+(nguard0+1)*k3d,            &
                         lb)
              End If
            End Do

            ElseIf (lb > lnblocks) Then

            vtype = 6
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_ec(2)
                ivar_next = gcell_on_ec_pointer(2,ivar)
                unk_e_y1(ivar_next,i,j,k,idest) =                      & 
                         temprecv_buf(index+ivar)
              End Do
              index = index + ngcell_on_ec(2)
            End Do
            End Do
            End Do

            End If  ! End If (lb <= lnblocks)

            If (ndim == 3) Then

            If (lb <= lnblocks) Then

            Do ivar=1,nvaredge
              If (int_gcell_on_ec(3,ivar)) Then
                 unk_e_z1(ivar,                                        & 
                          1+nguard:nxb+nguard+1,                       & 
                          1+nguard*k2d:nyb+(nguard+1)*k2d,             & 
                          1+nguard*k3d:nzb+nguard*k3d,                 &
                          idest) =                                     & 
                 unk_e_z(ivar,                                         &
                         1+nguard0:nxb+nguard0+1,                      & 
                         1+nguard0*k2d:nyb+(nguard0+1)*k2d,            & 
                         1+nguard0*k3d:nzb+nguard0*k3d,                &
                         lb)
              End If
            End Do

            ElseIf (lb > lnblocks) Then

            vtype = 7
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_ec(3)
                ivar_next = gcell_on_ec_pointer(3,ivar)
                unk_e_z1(ivar_next,i,j,k,idest) =                      & 
                        temprecv_buf(index+ivar)
              End Do
              index = index + ngcell_on_ec(3)
            End Do
            End Do
            End Do

            End If  ! End If (ln <= lnblocks)

            End If  ! End If (ndim == 3)

          End If  ! End If (lec)
        End If  ! End If (ndim > 1)

        If (lnc) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nvarcorn
              If (int_gcell_on_nc(ivar)) Then
                 unk_n1(ivar,                                          &
                        1+nguard:nxb+nguard+1,                         & 
                        1+nguard*k2d:nyb+(nguard+1)*k2d,               & 
                        1+nguard*k3d:nzb+(nguard+1)*k3d,               &
                        idest) =                                       & 
                 unk_n(ivar,                                           &
                       1+nguard0:nxb+nguard0+1,                        & 
                       1+nguard0*k2d:nyb+(nguard0+1)*k2d,              & 
                       1+nguard0*k3d:nzb+(nguard0+1)*k3d,              &
                       lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

!-----------starting index if cell-centered data is also included in recv_buf
            index = index0 + ngcell_on_cc*message_size_cc(dtype)       & 
                           + maxval(ngcell_on_fc(1:ndim))*             & 
                                     message_size_fc(dtype)            & 
                           + maxval(ngcell_on_ec(1:ndim))*             & 
                                     message_size_ec(dtype)

            vtype = 8
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_nc
                ivar_next = gcell_on_nc_pointer(ivar)
                unk_n1(ivar_next,i,j,k,idest) =                        &
                          temprecv_buf(index+ivar)
              End Do
              index = index+ngcell_on_nc
            End Do
            End Do
            End Do

          End If  ! End If (lb <= lnblocks)

        End If  ! End If (lnc)

      Else If (iopt >= 2) Then

          iopt0 = iopt-1

            If (lb <= lnblocks) Then

               work1(1+nguard_work:nxb+nguard_work,                       & 
                  1+nguard_work*k2d:nyb+nguard_work*k2d,               & 
                  1+nguard_work*k3d:nzb+nguard_work*k3d,               &
                  idest) =                                             & 
                          work(1+nguard_work0:nxb+nguard_work0,                      & 
                 1+nguard_work0*k2d:nyb+nguard_work0*k2d,              & 
                 1+nguard_work0*k3d:nzb+nguard_work0*k3d,              &
                 lb,                                                   &
                 iopt0)

            Else If (lb > lnblocks) Then
 
               vtype = 0
               index = index0
               Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

               If (no_permanent_guardcells) Then
                  ia = ia + nguard_work
                  ib = ib + nguard_work
                  ja = ja + nguard_work*k2d
                  jb = jb + nguard_work*k2d
                  ka = ka + nguard_work*k3d
                  kb = kb + nguard_work*k3d
               End If

               Do k = ka,kb
                  Do j = ja,jb
                     Do i = ia,ib
                        work1(i,j,k,idest) =                                     & 
                            temprecv_buf(index+1)
                        index = index + 1
                     End Do
                  End Do
               End Do

             End If  ! End If (lb <= lnblocks)

      End If  ! End If (iopt == 1)

      Return
      End Subroutine amr_perm_to_1blk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine pf_perm_to_1blk( lcc,lfc,lec,lnc,lb,pe,iopt, idest)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use mpi_morton
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguemtns
      Integer, intent(in) ::  lb,pe,iopt, idest
      Logical, intent(in) ::  lcc,lfc,lec,lnc


!-----Local variables and arrays.
      Integer :: iopt0
      Integer :: nguard0,nguard_work0
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: i, j, k, ivar, ivar_next
      Integer :: vtype,dtype,rem_blk,rem_pe,mype,ierr
      Integer :: index,index0
      Logical :: lfound

!-----Begin executable code.

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      If (lb > lnblocks) Then
           Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
           rem_blk = lb
           rem_pe  = mype
           Call amr_mpi_find_blk_in_buffer(mype,rem_blk,               & 
                             rem_pe,idest,dtype,index0,lfound)
           If (.Not.lfound) Then 
             Write(*,*)                                                & 
              'perm to 1blk reporting blk not found',                  & 
              ' mype=',mype,' looked for ',lb,                         & 
              ' where lnblocks=',lnblocks,                             & 
              ' strt_buffer=',strt_buffer,' last_buffer=',last_buffer, & 
              ' laddress ',laddress(:,strt_buffer:last_buffer)
             Call amr_abort()    ! remove this abort after testing
           End If
      End If  ! End If (lb > lnblocks)

!-----Put block lb's data into the data_1blk.fh datastructures, with the
!-----appropriate guardcell padding.
      If (iopt == 1) Then

        If (lcc) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nvar
              If (int_gcell_on_cc(ivar)) Then
                 unk1(ivar,                                            &
                      1+nguard:nxb+nguard,                             &
                      1+nguard*k2d:nyb+nguard*k2d,                     & 
                      1+nguard*k3d:nzb+nguard*k3d,                     &
                      idest) =                                         & 
                 unk(ivar,                                             &
                     1+nguard0:nxb+nguard0,                            &
                     1+nguard0*k2d:nyb+nguard0*k2d,                    & 
                     1+nguard0*k3d:nzb+nguard0*k3d,                    &
                     lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

             vtype = 1
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)
             index = index0

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_cc
                 ivar_next = gcell_on_cc_pointer(ivar)
                 unk1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
               End Do
               index = index + ngcell_on_cc
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

        End If  ! End If  (lcc)

        If (lfc) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nfacevar
              If (int_gcell_on_fc(1,ivar)) Then
                 facevarx1(ivar,                                       &
                           1+nguard:nxb+nguard+1,                      & 
                           1+nguard*k2d:nyb+nguard*k2d,                & 
                           1+nguard*k3d:nzb+nguard*k3d,                &
                           idest) =                                    & 
                 facevarx(ivar,                                        &
                          1+nguard0:nxb+nguard0+1,                     & 
                          1+nguard0*k2d:nyb+nguard0*k2d,               & 
                          1+nguard0*k3d:nzb+nguard0*k3d,lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

!------------starting index if cell-centered data is also included in recv_buf
             index = index0 + ngcell_on_cc*message_size_cc(dtype)

             vtype = 2
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_fc(1)
                 ivar_next = gcell_on_fc_pointer(1,ivar)
                 facevarx1(ivar_next,i,j,k,idest) =                    & 
                          temprecv_buf(index+ivar)
               End Do
               index = index+ngcell_on_fc(1)
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

          If (ndim >= 2) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nfacevar
              If (int_gcell_on_fc(2,ivar)) Then
                 facevary1(ivar,                                       &
                           1+nguard:nxb+nguard,                        & 
                           1+nguard*k2d:nyb+(nguard+1)*k2d,            & 
                           1+nguard*k3d:nzb+nguard*k3d,                &
                           idest) =                                    & 
                 facevary(ivar,                                        &
                          1+nguard0:nxb+nguard0,                       & 
                          1+nguard0*k2d:nyb+(nguard0+1)*k2d,           & 
                          1+nguard0*k3d:nzb+nguard0*k3d,               &
                          lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

             vtype = 3
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_fc(2)
                 ivar_next = gcell_on_fc_pointer(2,ivar)
                 facevary1(ivar_next,i,j,k,idest) =                    & 
                          temprecv_buf(index+ivar)
               End Do
               index = index + ngcell_on_fc(2)
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

          End If  ! end If (ndim >= 2)

          If (ndim == 3) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nfacevar
              If (int_gcell_on_fc(3,ivar)) Then
                 facevarz1(ivar,                                       &
                           1+nguard:nxb+nguard,                        & 
                           1+nguard*k2d:nyb+nguard*k2d,                & 
                           1+nguard*k3d:nzb+(nguard+1)*k3d,            &
                           idest) =                                    & 
                  facevarz(ivar,                                       &
                           1+nguard0:nxb+nguard0,                      & 
                           1+nguard0*k2d:nyb+nguard0*k2d,              & 
                           1+nguard0*k3d:nzb+(nguard0+1)*k3d,          &
                           lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

             vtype = 4
             Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

             If (no_permanent_guardcells) Then
                ia = ia + nguard
                ib = ib + nguard
                ja = ja + nguard*k2d
                jb = jb + nguard*k2d
                ka = ka + nguard*k3d
                kb = kb + nguard*k3d
             End If

             Do k = ka,kb
             Do j = ja,jb
             Do i = ia,ib
               Do ivar=1,ngcell_on_fc(3)
                 ivar_next = gcell_on_fc_pointer(3,ivar)
                 facevarz1(ivar_next,i,j,k,idest) =                    & 
                          temprecv_buf(index+ivar)
               End Do
               index = index + ngcell_on_fc(3)
             End Do
             End Do
             End Do

          End If  ! End If (lb <= lnblocks)

          End If  ! End If (ndim == 3)

        End If  ! End If (lfc)

        If (ndim > 1) Then
          If (lec) Then

            If (lb <= lnblocks) Then

            Do ivar=1,nvaredge
              If (int_gcell_on_ec(1,ivar)) Then
                 unk_e_x1(ivar,                                        &
                          1+nguard:nxb+nguard,                         & 
                          1+nguard*k2d:nyb+nguard*k2d+k2d,             & 
                          1+nguard*k3d:nzb+nguard*k3d+k3d,             &
                          idest) =                                     & 
                 unk_e_x(ivar,                                         &
                         1+nguard0:nxb+nguard0,                        & 
                         1+nguard0*k2d:nyb+nguard0*k2d+k2d,            & 
                         1+nguard0*k3d:nzb+nguard0*k3d+k3d,            &
                         lb)
              End If
            End Do

            ElseIf (lb > lnblocks) Then

!-----------starting index if cell-centered data is also included in recv_buf
            index = index0 + ngcell_on_cc*message_size_cc(dtype)       & 
                           + maxval(ngcell_on_fc(1:ndim))*             & 
                                   message_size_fc(dtype)

            vtype = 5
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)
  
            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_ec(1)
                ivar_next = gcell_on_ec_pointer(1,ivar)
                unk_e_x1(ivar_next,i,j,k,idest) =                      & 
                         temprecv_buf(index+ivar)
              End Do
              index = index + ngcell_on_ec(1)
            End Do
            End Do
            End Do

            End If  ! End If (lb <= lnblocks)

            If (lb <= lnblocks) Then

            Do ivar=1,nvaredge
              If (int_gcell_on_ec(2,ivar)) Then
                 unk_e_y1(ivar,                                        &
                          1+nguard:nxb+nguard+1,                       &  
                          1+nguard*k2d:nyb+nguard*k2d,                 & 
                          1+nguard*k3d:nzb+(nguard+1)*k3d,             &
                          idest) =                                     & 
                 unk_e_y(ivar,                                         &
                         1+nguard0:nxb+nguard0+1,                      & 
                         1+nguard0*k2d:nyb+nguard0*k2d,                & 
                         1+nguard0*k3d:nzb+(nguard0+1)*k3d,            &
                         lb)
              End If
            End Do

            ElseIf (lb > lnblocks) Then

            vtype = 6
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_ec(2)
                ivar_next = gcell_on_ec_pointer(2,ivar)
                unk_e_y1(ivar_next,i,j,k,idest) =                      & 
                         temprecv_buf(index+ivar)
              End Do
              index = index + ngcell_on_ec(2)
            End Do
            End Do
            End Do

            End If  ! End If (lb <= lnblocks)

            If (ndim == 3) Then

            If (lb <= lnblocks) Then

            Do ivar=1,nvaredge
              If (int_gcell_on_ec(3,ivar)) Then
                 unk_e_z1(ivar,                                        & 
                          1+nguard:nxb+nguard+1,                       & 
                          1+nguard*k2d:nyb+(nguard+1)*k2d,             & 
                          1+nguard*k3d:nzb+nguard*k3d,                 &
                          idest) =                                     & 
                 unk_e_z(ivar,                                         &
                         1+nguard0:nxb+nguard0+1,                      & 
                         1+nguard0*k2d:nyb+(nguard0+1)*k2d,            & 
                         1+nguard0*k3d:nzb+nguard0*k3d,                &
                         lb)
              End If
            End Do

            ElseIf (lb > lnblocks) Then

            vtype = 7
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_ec(3)
                ivar_next = gcell_on_ec_pointer(3,ivar)
                unk_e_z1(ivar_next,i,j,k,idest) =                      & 
                        temprecv_buf(index+ivar)
              End Do
              index = index + ngcell_on_ec(3)
            End Do
            End Do
            End Do

            End If  ! End If (ln <= lnblocks)

            End If  ! End If (ndim == 3)

          End If  ! End If (lec)
        End If  ! End If (ndim > 1)

        If (lnc) Then

          If (lb <= lnblocks) Then

            Do ivar=1,nvarcorn
              If (int_gcell_on_nc(ivar)) Then
                 unk_n1(ivar,                                          &
                        1+nguard:nxb+nguard+1,                         & 
                        1+nguard*k2d:nyb+(nguard+1)*k2d,               & 
                        1+nguard*k3d:nzb+(nguard+1)*k3d,               &
                        idest) =                                       & 
                 unk_n(ivar,                                           &
                       1+nguard0:nxb+nguard0+1,                        & 
                       1+nguard0*k2d:nyb+(nguard0+1)*k2d,              & 
                       1+nguard0*k3d:nzb+(nguard0+1)*k3d,              &
                       lb)
              End If
            End Do

          ElseIf (lb > lnblocks) Then

!-----------starting index if cell-centered data is also included in recv_buf
            index = index0 + ngcell_on_cc*message_size_cc(dtype)       & 
                           + maxval(ngcell_on_fc(1:ndim))*             & 
                                     message_size_fc(dtype)            & 
                           + maxval(ngcell_on_ec(1:ndim))*             & 
                                     message_size_ec(dtype)

            vtype = 8
            Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

            If (no_permanent_guardcells) Then
               ia = ia + nguard
               ib = ib + nguard
               ja = ja + nguard*k2d
               jb = jb + nguard*k2d
               ka = ka + nguard*k3d
               kb = kb + nguard*k3d
            End If

            Do k = ka,kb
            Do j = ja,jb
            Do i = ia,ib
              Do ivar=1,ngcell_on_nc
                ivar_next = gcell_on_nc_pointer(ivar)
                unk_n1(ivar_next,i,j,k,idest) =                        &
                          temprecv_buf(index+ivar)
              End Do
              index = index+ngcell_on_nc
            End Do
            End Do
            End Do

          End If  ! End If (lb <= lnblocks)

        End If  ! End If (lnc)

      Else If (iopt >= 2) Then

          iopt0 = iopt-1
!         do iopt0 = iopt-1, nvar_work

            If (lb <= lnblocks) Then
!if (lb.eq.105.or.lb.eq.121) print *, 'perm :lb,iopt0,work', lb, iopt0, work(2,2,1,lb,iopt0)
               work1(1+nguard_work:nxb+nguard_work,                       & 
                  1+nguard_work*k2d:nyb+nguard_work*k2d,               & 
                  1+nguard_work*k3d:nzb+nguard_work*k3d,               &
                  idest) =                                             & 
!                  iopt0) =                                             & !was idest
                          work(1+nguard_work0:nxb+nguard_work0,                      & 
                 1+nguard_work0*k2d:nyb+nguard_work0*k2d,              & 
                 1+nguard_work0*k3d:nzb+nguard_work0*k3d,              &
                 lb,                                                   &
                 iopt0)

            Else If (lb > lnblocks) Then
 
               vtype = 0
               index = index0
               Call mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype)

               If (no_permanent_guardcells) Then
                  ia = ia + nguard_work
                  ib = ib + nguard_work
                  ja = ja + nguard_work*k2d
                  jb = jb + nguard_work*k2d
                  ka = ka + nguard_work*k3d
                  kb = kb + nguard_work*k3d
               End If

               Do k = ka,kb
                  Do j = ja,jb
                     Do i = ia,ib
                        work1(i,j,k,idest) =                                     &
!                        work1(i,j,k,iopt0) =                                     &! was idest 
                            temprecv_buf(index+1)
                        index = index + 1
                     End Do
                  End Do
               End Do

             End If  ! End If (lb <= lnblocks)
! CEG end iopt0 do
!         End Do

      End If  ! End If (iopt == 1)

      Return
      End Subroutine pf_perm_to_1blk

