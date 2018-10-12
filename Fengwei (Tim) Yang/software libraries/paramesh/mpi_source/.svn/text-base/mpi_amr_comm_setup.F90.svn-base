!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_comm_setup
!! NAME
!!
!!   mpi_amr_comm_setup
!! 
!! SYNOPSIS
!!
!!   Call mpi_amr_comm_setup(mype,nprocs,
!!                           lguard,lprolong,
!!                           lflux,ledge,lrestrict,lfulltree,
!!                           iopt,lcc,lfc,lec,lnc,tag_offset,
!!                           nlayersx,nlayersy,nlayersz,
!!                           flux_dir)
!!   Call mpi_amr_comm_setup(mype,nprocs,
!!                           lguard,lprolong,
!!                           lflux,ledge,lrestrict,lfulltree,
!!                           iopt,lcc,lfc,lec,lnc,tag_offset)
!!
!!   Call mpi_amr_comm_setup(integer, integer,
!!                           logical, logical,
!!                           logical, logical, logical, logical,
!!                           integer, logical, logical, logical, logical, integer
!!                           optional integer, optional integer, optional, integer,
!!                           optional integer)
!!
!! ARGUMENTS      
!!
!!    integer, intent(in) :: mype,nprocs,iopt
!!      mype   -> The local processor number.
!!      nprocs -> The total number of processors being used.
!!       
!!   logical, intent(in)  :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
!!      lguard    -> If true, info for guardcell filling will be
!!                    communicated.
!!      lprolong  -> If true, info for prolongation will be
!!                    communicated.
!!      lflux     -> If true, info for flux conservation will be
!!                    communicated.
!!      ledge     -> If true, info for edge averaging will be
!!                    communicated.
!!      lrestrict -> If true, info for restriction will be
!!                    communicated.
!!      lfulltree -> If true, info for restricting data up the entire
!!                    tree will be communicated.
!!
!!   integer, intent(in) :: iopt
!!      iopt   -> Controls which arrays are operated on:
!!                If iopt = 1 then 'unk', 'facevarx(y,z)', 'unk_e_x(y,z)'
!!                and 'unk_n'
!!                If iopt >= 2 then 'work'.
!!
!!   logical, intent(in) :: lcc,lfc,lec,lnc
!!      lcc -> A logical switch controlling whether unk or work data
!!              is filled.
!!      lfc -> A logical switch controlling whether facevar data
!!              is filled.
!!      lec -> A logical switch controlling whether unk_e_? data
!!              is filled.
!!      lnc -> A logical switch controlling whether unk_n data
!!              is filled.
!!
!!   integer, intent(inout) :: tag_offset
!!      tag_offset -> Defines the last tag number used for an mpi message.
!!                    This can be almost anything, but 100 is usually a 
!!                    good choice.
!!
!!   integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
!!     Optional integer arguments specifying how many guardcells to 
!!      exchange in each coordinate direction.  If these are not 
!!      specified then all guardcells are filled.
!!
!!   integer, intent(in), optional :: flux_dir
!!     Optional argument specifying which direction to operate on 
!!       for flux conservation. 
!!       If flux_dir = 1, operate on 'x' direction.
!!       If flux_dir = 2, operate on 'y' direction.
!!       If flux_dir = 3, operate on 'z' direction.
!!     If not specified, all directions are operated on.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   tree
!!   timings
!!   mpi_morton
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   mpi_pack_blocks
!!   mpi_Sbuffer_size
!!   mpi_unpack_blocks
!!   mpi_Rbuffer_size
!!   mpi_pack_edges
!!   mpi_unpack_edges
!!   mpi_pack_fluxes
!!   mpi_unpack_fluxes
!!   mpi_pack_tree_info
!!   mpi_unpack_tree_info
!!   mpi_amr_read_guard_comm
!!   mpi_amr_read_prol_comm
!!   mpi_amr_read_flux_comm
!!   mpi_amr_read_restrict_comm
!!   mpi_set_message_sizes
!!   mpi_xchange_blocks
!!
!! RETURNS
!!
!!   Upon exit, all necessary block data has been communicated to the local
!!   Calling processor for the operation selected through the specified
!!   arguments.
!!
!! DESCRIPTION
!!
!!   This routine is the pre-analysis operation which precedes a code block 
!!   in which interprocessor communications may be required.  What is does
!!   is set up and communicate the necessary block information for performing
!!   various communications tasks.  Those tasks are:
!!     1) guardcell filling
!!     2) prolongation
!!     3) restriction
!!     4) flux conservation
!!     5) edge averaging
!!
!! AUTHORS
!!
!!   Peter MacNeice (June 2000) with modifications by Kevin Olson for
!!   directional guardcell filling and flux conservation.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_amr_comm_setup(mype,nprocs,                       & 
                                    lguard,lprolong,                   & 
                                    lflux,ledge,lrestrict,lfulltree,   & 
                                    iopt,lcc,lfc,lec,lnc,tag_offset,   & 
                                    nlayersx,nlayersy,nlayersz,        & 
                                    flux_dir)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use tree
      Use timings
      Use mpi_morton
      Use paramesh_mpi_interfaces, Only :                              & 
                                      mpi_pack_blocks,                 & 
                                      mpi_Sbuffer_size,                & 
                                      mpi_unpack_blocks,               & 
                                      mpi_Rbuffer_size,                & 
                                      mpi_pack_edges,                  & 
                                      mpi_unpack_edges,                & 
                                      mpi_pack_fluxes,                 & 
                                      mpi_unpack_fluxes,               & 
                                      mpi_pack_tree_info,              & 
                                      mpi_unpack_tree_info,            & 
                                      mpi_amr_read_guard_comm,         & 
                                      mpi_amr_read_prol_comm,          & 
                                      mpi_amr_read_flux_comm,          & 
                                      mpi_amr_read_restrict_comm,      & 
                                      mpi_set_message_sizes,           & 
                                      mpi_xchange_blocks

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in)    :: mype,nprocs,iopt
      Integer, Intent(inout) :: tag_offset
      Logical, Intent(in)    :: lcc,lfc,lec,lnc,lfulltree
      Logical, Intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict
      Integer, Intent(in), Optional :: nlayersx,nlayersy,nlayersz
      Integer, Intent(in), Optional :: flux_dir

!-----Local variables and arrays.
      Real, Save, Dimension (:), allocatable :: send_buf
      Real, Save, Dimension (:), allocatable :: recv_buf
      Integer :: buffer_dim_send, buffer_dim_recv
      Integer :: buffer_dim, buf_dim_bytes
      Integer :: bufdim
      Integer :: max_blks_sent
      Integer :: len_surr_blks
      Integer :: itemp
      Integer :: offset_tree
      Integer :: nlayerstx, nlayersty, nlayerstz
      Integer :: flux_dirt
      Integer :: ierror
      Integer :: ii,jj
      Logical, Save :: first_Call = .True.

!-----Begin executable code.
      If (iopt == 1) Then

!-----install user selection for guardcell variables on reset defaults.
      If (lguard) Then
!        print *,'gcell_on_cc',gcell_on_cc
!CEG removed
!        int_gcell_on_cc = gcell_on_cc
        int_gcell_on_fc = gcell_on_fc
        int_gcell_on_ec = gcell_on_ec
        int_gcell_on_nc = gcell_on_nc
        lguard_in_progress = .True.

        If (nvar > 0) Then
          ngcell_on_cc = 0
          Do ii = 1,nvar
            If (int_gcell_on_cc(ii)) Then
              ngcell_on_cc =  ngcell_on_cc + 1
              gcell_on_cc_pointer(ngcell_on_cc) = ii
            End If
          End Do
        End If

        If (nfacevar > 0) Then
          ngcell_on_fc = 0
          Do ii = 1,nfacevar
          Do jj = 1,3
            If (int_gcell_on_fc(jj,ii)) Then
              ngcell_on_fc(jj) =  ngcell_on_fc(jj) + 1
              gcell_on_fc_pointer(jj,ngcell_on_fc(jj)) = ii
            End If
          End Do
          End Do
        End If

        If (nvaredge > 0) Then
          ngcell_on_ec = 0
          Do ii = 1,nvaredge
          Do jj = 1,3
            If (int_gcell_on_ec(jj,ii)) Then
              ngcell_on_ec(jj) =  ngcell_on_ec(jj) + 1
              gcell_on_ec_pointer(jj,ngcell_on_ec(jj)) = ii
            End If
          End Do
          End Do
        End If

        If (nvarcorn > 0) Then
          ngcell_on_nc = 0
          Do ii = 1,nvarcorn
            If (int_gcell_on_nc(ii)) Then
              ngcell_on_nc =  ngcell_on_nc + 1
              gcell_on_nc_pointer(ngcell_on_nc) = ii
            End If
          End Do
        End If

      Else                                

!CEG removed
!        int_gcell_on_cc(:) = .True.
        int_gcell_on_fc(:,:) = .True.
        int_gcell_on_ec(:,:) = .True.
        int_gcell_on_nc(:) = .True.
        lguard_in_progress = .False.
        ngcell_on_cc = nvar
        Do ii=1,nvar
          gcell_on_cc_pointer(ii) = ii
        End Do
        ngcell_on_fc = nfacevar
        Do ii=1,nfacevar
          gcell_on_fc_pointer(:,ii) = ii
        End Do
        ngcell_on_ec = nvaredge
        Do ii=1,nvaredge
          gcell_on_ec_pointer(:,ii) = ii
        End Do
        ngcell_on_nc = nvarcorn
        Do ii=1,nvarcorn
          gcell_on_nc_pointer(ii) = ii
        End Do

      End If  ! End If (lguard)

      If (Present(nlayersx)) Then
         nlayerstx = nlayersx
      Else
         nlayerstx = nguard
      End If

      If (Present(nlayersy)) Then
         nlayersty = nlayersy
      Else
         nlayersty = nguard
      End If

      If (Present(nlayersz)) Then
         nlayerstz = nlayersz
      Else
         nlayerstz = nguard
      End If

      Else

      If (Present(nlayersx)) Then
         nlayerstx = nlayersx
      Else
         nlayerstx = nguard_work
      End If

      If (Present(nlayersy)) Then
         nlayersty = nlayersy
      Else
         nlayersty = nguard_work
      End If

      If (Present(nlayersz)) Then
         nlayerstz = nlayersz
      Else
         nlayerstz = nguard_work
      End If

      End If  ! End If (iopt == 1)

      If (lguard) Then
         If (nxb/nguard < 2) nlayerstx = min(nlayerstx+1,   nguard)
         If (nyb/nguard < 2) nlayersty = min(nlayersty+k2d, nguard)
         If (nzb/nguard < 2) nlayerstz = min(nlayerstz+k3d, nguard)
      End If
      
      If (Present(flux_dir)) Then
         flux_dirt = flux_dir
      Else
         flux_dirt = 0
      End If

      Call mpi_set_message_sizes(iopt,nlayerstx,nlayersty,nlayerstz)

      If (lguard.and.(.not.lrestrict) .or. lfulltree ) Then

         Call mpi_amr_read_guard_comm(nprocs)

      ElseIf (lprolong) Then

         Call mpi_amr_read_prol_comm(nprocs)

      ElseIf ((lflux.or.ledge).and.(.not.lrestrict)) Then

        Call mpi_amr_read_flux_comm(nprocs)

      ElseIf (lrestrict) Then

        Call mpi_amr_read_restrict_comm(nprocs)

      End If

!-----If lrestrict is true Then we only need tree information.
      If ( (.not.lrestrict) .or. (lrestrict.and.lflux)                 & 
                            .or. (lrestrict.and.ledge)                 & 
                            .or. (lrestrict.and.lguard) ) Then

      If (lguard.or.lprolong) Then

        Call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_send,offset_tree,             & 
                              .True., .False., .False., flux_dir,      &
                              nlayerstx,nlayersty,nlayerstz)

        Call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_recv,                         & 
                              .True.,.False., .False., flux_dir,       &
                              nlayerstx,nlayersty,nlayerstz)

      ElseIf (lflux) Then

        Call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_send,offset_tree,             & 
                              .False., .True., .False., flux_dir,      &
                              nlayerstx,nlayersty,nlayerstz)

        Call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_recv,                         & 
                              .False.,.True., .False., flux_dir,       &
                              nlayerstx,nlayersty,nlayerstz)

      ElseIf (ledge) Then

        Call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_send,offset_tree,             & 
                              .False., .False., .True., flux_dir,      &
                              nlayerstx,nlayersty,nlayerstz)

        Call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_recv,                         & 
                              .False.,.False., .True., flux_dir,       &
                              nlayerstx,nlayersty,nlayerstz)

      End If  ! End If (lguard.or.lprolong)

!-----Set up buffer size info, including space for tree data
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
      offset_tree = 32+16+len_surr_blks

      If (Allocated(send_buf)) Deallocate(send_buf)
      Allocate( send_buf(buffer_dim_send))
      If (Allocated(temprecv_buf)) Deallocate(temprecv_buf)
      Allocate( temprecv_buf(buffer_dim_recv))

      If (lguard.or.lprolong) Then

        Call mpi_pack_blocks(mype,nprocs,iopt,lcc,lfc,lec,lnc,         & 
                             buffer_dim_send,send_buf,offset_tree,     & 
                             nlayerstx,nlayersty,nlayerstz)

      ElseIf (lflux) Then

        Call mpi_pack_fluxes(mype,nprocs,buffer_dim_send,send_buf,     & 
                             offset_tree,flux_dir)

      ElseIf (ledge) Then

        Call mpi_pack_edges(mype,nprocs,buffer_dim_send,send_buf,      & 
                            offset_tree)

      End If  ! End If (lguard.or.lprolong)

      Call mpi_xchange_blocks(mype,nprocs,tag_offset,                  & 
                              buffer_dim_send,send_buf,                &
                              buffer_dim_recv,temprecv_buf)


      If (lguard.or.lprolong) Then

         Call mpi_unpack_blocks(mype,iopt,lcc,lfc,lec,lnc,             & 
                                buffer_dim_recv,temprecv_buf,          & 
                                nlayerstx,nlayersty,nlayerstz)

      ElseIf (lflux) Then
         
         Call mpi_unpack_fluxes(mype,buffer_dim_recv,temprecv_buf,     & 
                                flux_dir)

      ElseIf (ledge) Then

         Call mpi_unpack_edges(mype,buffer_dim_recv,temprecv_buf)

      End If  ! End If (lguard.or.lprolong)

      If (Allocated(send_buf)) Deallocate(send_buf)

      Else

! CEG moved this from out of the IF as the AllReduce is very costly and only needed here
        itemp = max(sum(commatrix_send), sum(commatrix_recv))
        Call MPI_ALLREDUCE (itemp,                                       & 
                          max_blks_sent,                               & 
                          1,                                           & 
                          MPI_INTEGER,                                 & 
                          MPI_MAX,                                     & 
                          MPI_COMM_WORLD,                              & 
                          ierror)
!  write (6,*) 'CEG - Allreduce in mpi_amr_comm_setup  : Faked at 40 : max_blks_sent=', itemp
!      max_blks_sent=40

!-----Establish tree info which is needed about remote blocks
!-----Set up buffer size info for tree messages
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
      buffer_dim = (32+16+len_surr_blks)*max_blks_sent

      If (Allocated(send_buf)) Deallocate(send_buf)
      Allocate( send_buf(buffer_dim))
      If (Allocated(recv_buf)) Deallocate(recv_buf)
      Allocate( recv_buf(buffer_dim))

!-----Pack the tree info to be sent
      Call mpi_pack_tree_info(mype,nprocs,buf_dim_bytes,               & 
                              buffer_dim,send_buf)

!-----Exchange the tree info to be sent
      Call mpi_xchange_blocks(mype,nprocs,tag_offset,                  & 
                              buffer_dim, send_buf,                    &
                              buffer_dim, recv_buf)


!-----Unpack the tree info which has been received
      Call mpi_unpack_tree_info(mype,nprocs,buf_dim_bytes,             & 
                                buffer_dim,recv_buf)


      If (Allocated(send_buf)) Deallocate(send_buf)
      If (Allocated(recv_buf)) Deallocate(recv_buf)

      End If  ! End If ( (.not.lrestrict) .or. (lrestrict.and.lflux)

!-----Store laddress pe limits to help optimize searching
      If (last_buffer > strt_buffer.and.nprocs > 1) Then

         ladd_strt(1:nprocs-1) = last_buffer
         ladd_strt(0         ) = strt_buffer
         Do jj= last_buffer,strt_buffer,-1
           ladd_strt(laddress(2,jj)) = jj
         End Do
         Do jj= nprocs-2,0,-1
           ladd_strt(jj) = minval(ladd_strt(jj:jj+1))
         End Do

         ladd_end(0:nprocs-2) = strt_buffer
         ladd_end(nprocs-1  ) = last_buffer
         Do jj= strt_buffer,last_buffer
           ladd_end(laddress(2,jj)) = jj
         End Do
         Do jj= 1,nprocs-1
           ladd_end(jj) = maxval(ladd_end(jj-1:jj))
         End Do

      Else

        ladd_strt = strt_buffer
        ladd_end  = strt_buffer

      End If

      Return
      End Subroutine mpi_amr_comm_setup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine pf_mpi_amr_comm_setup(mype,nprocs,                       & 
                                    lguard,lprolong,                   & 
                                    lflux,ledge,lrestrict,lfulltree,   & 
                                    iopt,lcc,lfc,lec,lnc,tag_offset,   & 
                                    mglevel)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use tree
      Use timings
      Use mpi_morton
      Use paramesh_mpi_interfaces, Only :                              & 
                                      mpi_pack_blocks,                 & 
                                      mpi_Sbuffer_size,                & 
                                      mpi_unpack_blocks,               & 
                                      mpi_Rbuffer_size,                & 
                                      mpi_pack_edges,                  & 
                                      mpi_unpack_edges,                & 
                                      mpi_pack_fluxes,                 & 
                                      mpi_unpack_fluxes,               & 
                                      mpi_pack_tree_info,              & 
                                      mpi_unpack_tree_info,            & 
                                      mpi_amr_read_guard_comm,         & 
                                      mpi_amr_read_prol_comm,          & 
                                      mpi_amr_read_flux_comm,          & 
                                      mpi_amr_read_restrict_comm,      & 
                                      mpi_set_message_sizes,           & 
                                      mpi_xchange_blocks

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in)    :: mype,nprocs,iopt
      Integer, Intent(inout) :: tag_offset
      Logical, Intent(in)    :: lcc,lfc,lec,lnc,lfulltree
      Logical, Intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict
      Integer, Intent(in)    :: mglevel

!-----Local variables and arrays.
      Real, Save, Dimension (:), allocatable :: send_buf
      Real, Save, Dimension (:), allocatable :: recv_buf
      Integer :: buffer_dim_send, buffer_dim_recv
      Integer :: buffer_dim, buf_dim_bytes
      Integer :: bufdim
      Integer :: max_blks_sent
      Integer :: len_surr_blks
!      Integer :: itemp
      Integer :: offset_tree
      Integer :: nlayerstx, nlayersty, nlayerstz
      Integer :: flux_dirt
      Integer :: ierror
      Integer :: ii,jj
      Logical, Save :: first_Call = .True.

!-----Begin executable code.
      If (iopt == 1) Then
!-----install user selection for guardcell variables on reset defaults.
        int_gcell_on_cc = gcell_on_cc
        int_gcell_on_fc = gcell_on_fc
        int_gcell_on_ec = gcell_on_ec
        int_gcell_on_nc = gcell_on_nc
        lguard_in_progress = .True.

        ngcell_on_cc = 0
        Do ii = 1,nvar
          If (int_gcell_on_cc(ii)) Then
            ngcell_on_cc =  ngcell_on_cc + 1
            gcell_on_cc_pointer(ngcell_on_cc) = ii
          End If
        End Do

        nlayerstx = nguard
        nlayersty = nguard
        nlayerstz = nguard
      Else
        nlayerstx = nguard_work
        nlayersty = nguard_work
        nlayerstz = nguard_work
      End If  ! End If (iopt == 1)

      Call pf_mpi_set_message_sizes(iopt)
      
      If (.not.lrestrict) Then
         Call mpi_amr_read_guard_comm(nprocs)
      Else
         Call mpi_amr_read_restrict_comm(nprocs)
      End If

      Call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_send,offset_tree,             & 
                              .True., .False., .False.)
                              
      Call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_recv,                         & 
                              .True.,.False., .False.)

!-----Set up buffer size info, including space for tree data
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
      offset_tree = 32+16+len_surr_blks

!      If (Allocated(send_buf)) Deallocate(send_buf)
      Allocate( send_buf(buffer_dim_send))
      If (Allocated(temprecv_buf)) Deallocate(temprecv_buf)
      Allocate( temprecv_buf(buffer_dim_recv))

      Call pf_mpi_pack_blocks(mype,nprocs,iopt,lcc,lfc,lec,lnc,         & 
                             buffer_dim_send,send_buf,offset_tree,     & 
                             nlayerstx,nlayersty,nlayerstz)

      Call mpi_xchange_blocks(mype,nprocs,tag_offset,                  & 
                              buffer_dim_send,send_buf,                &
                              buffer_dim_recv,temprecv_buf)

      Call mpi_unpack_blocks(mype,iopt,lcc,lfc,lec,lnc,             & 
                                buffer_dim_recv,temprecv_buf,          & 
                                nlayerstx,nlayersty,nlayerstz)

!      If (Allocated(send_buf)) 
      Deallocate(send_buf)

!-----Store laddress pe limits to help optimize searching
      If (last_buffer > strt_buffer.and.nprocs > 1) Then

         ladd_strt(1:nprocs-1) = last_buffer
         ladd_strt(0         ) = strt_buffer
         Do jj= last_buffer,strt_buffer,-1
           ladd_strt(laddress(2,jj)) = jj
         End Do
         Do jj= nprocs-2,0,-1
           ladd_strt(jj) = minval(ladd_strt(jj:jj+1))
         End Do

         ladd_end(0:nprocs-2) = strt_buffer
         ladd_end(nprocs-1  ) = last_buffer
         Do jj= strt_buffer,last_buffer
           ladd_end(laddress(2,jj)) = jj
         End Do
         Do jj= 1,nprocs-1
           ladd_end(jj) = maxval(ladd_end(jj-1:jj))
         End Do

      Else

         ladd_strt = strt_buffer
         ladd_end  = strt_buffer

      End If

      Return
      End Subroutine pf_mpi_amr_comm_setup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! New routine added by CEG to clean up Adaptive MG memory

      subroutine mpi_amr_comm_free()
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use tree
      Use timings
      Use mpi_morton

      deallocate(temprecv_buf)
!      deallocate()
!      deallocate()
!      deallocate()
!      deallocate()
!      deallocate()

      Return
      End Subroutine mpi_amr_comm_free

