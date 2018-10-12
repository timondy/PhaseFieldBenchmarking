!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_edge_average_udt
!! NAME
!!
!!   amr_edge_average_udt
!!
!! SYNOPSIS
!!
!!   call amr_edge_average_udt(mype)
!!
!!   call amr_edge_average_udt(integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype
!!     The calling processors id.
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
!!   tree
!!   mpi_morton
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_restrict_edge_data
!!   mpi_amr_comm_setup
!!   mpi_put_edge_buffer_1blk
!!
!!
!! RETURNS
!!
!!   Nothing returned.  Upon return the edge data at refinement jumps is
!!   properly averaged.
!!
!! DESCRIPTION
!! 
!!   This routine gets cell edge-based data at block boundaries from 
!!   neighbors who are parents of leaf blocks. 
!!
!!   The data structure used to store and pass this data is defined
!!   in the module 'physicaldata'.
!!
!!   This version is called when uniform timesteps are being used across
!!   the blocks in the computation.  It is called by the wrapper routine
!!   'amr_edge_average'.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (July 1997).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_edge_average_udt(mype)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_interfaces, only : amr_restrict_edge_data
      Use paramesh_mpi_interfaces, only : mpi_amr_comm_setup,          &  
                                          mpi_put_edge_buffer_1blk

      Implicit None

!-----Include statements.
      include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in)  ::  mype

!-----Local arrays and variables.
      Integer :: nguard0
      Integer :: ng_off
      Integer :: kup,klo,kup1
      Integer :: remote_pe,remote_block
      Integer :: remote_pe2,remote_block2
      Integer :: cnodetype
      Integer :: tag_offset,nprocs,iopt,ierr
      Integer :: lb, jf, iblk
      Integer :: i, j, k, ia0, ib0, ja0, jb0, ka0, kb0
      Integer :: dtype, vtype, index, index0
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: n
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      Logical :: lcc,lfc,lec,lnc
      Logical :: lfound

!-----Begin executable code.
      nguard0 = nguard*npgs
      ng_off = nguard0+iface_off

      klo  = 1+k3d*nguard0
      kup  = 1+k3d*(nzb+nguard0-1)
      kup1 = k3d+nzb+k3d*nguard0

      If (ndim >= 2) Then

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      tag_offset = 100
      tag_offset = 100
      lcc = .False.
      lfc = .False.
      lec = .False.
      lnc = .False.
      iopt = 1

!-----Note, both ledge and lrestrict are true so that the fluxes
!-----are acquired which are needed in the restriction operation.
      lguard    = .False.
      lprolong  = .False.
      lflux     = .False.
      ledge     = .True.
      lrestrict = .True.
      lfulltree = .False.
      Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,             & 
                              lflux,ledge,lrestrict,lfulltree,         & 
                              iopt,lcc,lfc,lec,lnc,tag_offset)

!-----all leaf blocks provide reduced boundary edge data to their parents
      Call amr_restrict_edge_data(mype)

      tag_offset = 100

      lguard    = .False.
      lprolong  = .False.
      lflux     = .False.
      ledge     = .True.
      lrestrict = .False.
      lfulltree = .False.
      Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,             & 
                              lflux,ledge,lrestrict,lfulltree,         & 
                              iopt,lcc,lfc,lec,lnc,tag_offset)

!-----cycle through the grid blocks on this processor
      If (lnblocks > 0) Then
      Do lb = 1,lnblocks

!-----Is this a leaf block and not at the original refinement level ?
      If (nodetype(lb) == 1) Then

!-------Cycle over the blocks faces
        Do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)
          remote_pe2 = neigh(2,jf,lb)
          remote_block2  = neigh(1,jf,lb)
          cnodetype = 0
          lfound = .False.


          If (remote_block > 0) Then

!---------(remote_block,remote_pe) may be a local block, a remote block,
!---------or it may not exist.
!---------If it is a local block Then check its nodetype.
!---------If it is found in the list of remote blocks stored in buffer space
!---------then check its nodetype.
!---------If it is not found in either of these places, Then set its nodetype
!---------to 0.
         If (remote_pe2.ne.mype) Then
            lfound = .False.
            Do iblk = strt_buffer,last_buffer
               If (remote_block2 == laddress(1,iblk).And.              & 
                    remote_pe2  == laddress(2,iblk) ) Then
                  remote_block2 = iblk
                  remote_pe2    = mype
                  lfound = .True.
               End If
            End Do
         ElseIf (remote_pe2 == mype) Then
            lfound = .True.
         End If

!--------Is the neighbor to this face a parent of a leaf block?
         If (lfound) Then
            cnodetype = nodetype(remote_block2)
         End If
         
         End If  ! End If (remote_block > 0)

         If (cnodetype == 2) Then

! If yes Then copy the appropriate layer from its boundary variable data 

          If (remote_pe == mype .And. remote_block <= lnblocks) Then

            If (jf == 1) Then

              bedge_facex_y(:,1,:,:,lb) =                              &
                 bedge_facex_y(:,2,:,:,remote_block)
              If ((ndim == 3).Or.(l2p5d == 1)) Then
                bedge_facex_z(:,1,:,:,lb) =                            &
                   bedge_facex_z(:,2,:,:,remote_block)
              End If

            ElseIf (jf == 2) Then

              bedge_facex_y(:,2,:,:,lb) =                              &
                 bedge_facex_y(:,1,:,:,remote_block)
              If ((ndim == 3).Or.(l2p5d == 1)) Then
                bedge_facex_z(:,2,:,:,lb) =                            &
                   bedge_facex_z(:,1,:,:,remote_block)
              End If

            ElseIf (jf == 3) Then

              If ((ndim == 3).Or.(l2p5d == 1)) Then
                bedge_facey_z(:,:,1,:,lb) =                            &
                   bedge_facey_z(:,:,2,:,remote_block)
              End If
              bedge_facey_x(:,:,1,:,lb) =                              &
                 bedge_facey_x(:,:,2,:,remote_block)

            ElseIf (jf == 4) Then
                 
              If ((ndim == 3).Or.(l2p5d == 1)) Then
                 bedge_facey_z(:,:,2,:,lb) =                           &
                    bedge_facey_z(:,:,1,:,remote_block)
              End If
              bedge_facey_x(:,:,2,:,lb) =                              &
                 bedge_facey_x(:,:,1,:,remote_block)
              
            ElseIf (jf == 5) Then
              
              bedge_facez_x(:,:,:,1,lb) =                              &
                 bedge_facez_x(:,:,:,2,remote_block)
              bedge_facez_y(:,:,:,1,lb) =                              &
                 bedge_facez_y(:,:,:,2,remote_block)
              
            ElseIf (jf == 6) Then
              
              bedge_facez_x(:,:,:,2,lb) =                              &
                 bedge_facez_x(:,:,:,1,remote_block)
              bedge_facez_y(:,:,:,2,lb) =                              &
                 bedge_facez_y(:,:,:,1,remote_block)
              
            End If  ! End If (jf == 1)

           Else

            Call mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)

            If (jf == 1) Then
               bedge_facex_y(:,1,:,:,lb) = recvarx1e(:,2,:,:)
            ElseIf (jf == 2) Then
               bedge_facex_y(:,2,:,:,lb) = recvarx1e(:,1,:,:)
            End If

            If (jf == 1) Then
               If ((ndim == 3).Or.(l2p5d == 1)) Then
                  bedge_facex_z(:,1,:,:,lb) = recvarx2e(:,2,:,:)
               End If
            ElseIf (jf == 2) Then
               If ((ndim == 3).Or.(l2p5d == 1)) Then
                  bedge_facex_z(:,2,:,:,lb) = recvarx2e(:,1,:,:)
               End If
            End If

            If (jf == 3) Then
               bedge_facey_x(:,:,1,:,lb) = recvary1e(:,:,2,:)
            ElseIf (jf == 4) Then
               bedge_facey_x(:,:,2,:,lb) = recvary1e(:,:,1,:)
            End If

            If (jf == 3) Then
               If ((ndim == 3).Or.(l2p5d == 1)) Then
                  bedge_facey_z(:,:,1,:,lb) = recvary2e(:,:,2,:)
               End If
            ElseIf (jf == 4) Then
               If ((ndim == 3).Or.(l2p5d == 1)) Then
                  bedge_facey_z(:,:,2,:,lb) = recvary2e(:,:,1,:)
               End If
            End If

            If (jf == 5) Then
               bedge_facez_x(:,:,:,1,lb) = recvarz1e(:,:,:,2)
            ElseIf (jf == 6) Then
               bedge_facez_x(:,:,:,2,lb) = recvarz1e(:,:,:,1)
            End If

            If (jf == 5) Then
               bedge_facez_y(:,:,:,1,lb) = recvarz2e(:,:,:,2)
            ElseIf (jf == 6) Then
               bedge_facez_y(:,:,:,2,lb) = recvarz2e(:,:,:,1)
            End If

         End If  ! End If (remote_pe == mype .And. remote_block <= lnblocks)

!--------make common variables on an edge consistent
         If (jf == 1) Then

            bedge_facey_z(:,1+nguard0,1,klo:kup,lb) =                  & 
                 bedge_facex_z(:,1,1+nguard0*k2d,klo:kup,lb)
            
            bedge_facey_z(:,1+nguard0,2,klo:kup,lb) =                  & 
                 bedge_facex_z(:,1,k2d+nyb+nguard0*k2d,klo:kup,lb)
            
            If ((ndim == 3).Or.(l2p5d == 1)) Then
               bedge_facez_y(:,1+nguard0,                              & 
                    1+nguard0*k2d:nyb+nguard0*k2d,1,lb)                & 
                    = bedge_facex_y(:,1,                               & 
                    1+nguard0*k2d:nyb+nguard0*k2d,klo,lb)
               
               If (ndim == 3)                                          & 
                    bedge_facez_y(:,1+nguard0,                         & 
                    1+nguard0*k2d:nyb+nguard0*k2d,2,lb)                & 
                    = bedge_facex_y(:,1,                               & 
                    1+nguard0*k2d:nyb+nguard0*k2d,kup1,lb)
            End If

         ElseIf (jf == 2) Then

            bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb) =              & 
                 bedge_facex_z(:,2,1+nguard0*k2d,klo:kup,lb)

            bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb) =              & 
                 bedge_facex_z(:,2,k2d+nyb+nguard0*k2d,                & 
                                klo:kup,lb)

            If ((ndim == 3).Or.(l2p5d == 1)) Then
               bedge_facez_y(:,1+nxb+nguard0,                          & 
                    1+nguard0*k2d:nyb+nguard0*k2d,                     & 
                    1,lb)= & 
                    bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,   & 
                    klo,lb)
               
               If (ndim == 3)                                          & 
                    bedge_facez_y(:,1+nxb+nguard0,                     & 
                    1+nguard0*k2d:nyb+nguard0*k2d,                     & 
                    2,lb)=                                             & 
                    bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,   & 
                    kup1,lb)

            End If

         ElseIf (jf == 3) Then

            bedge_facex_z(:,1,1+nguard0*k2d,klo:kup,lb) =              & 
                 bedge_facey_z(:,1+nguard0,1,klo:kup,lb)
            
            bedge_facex_z(:,2,1+nguard0*k2d,klo:kup,lb) =              & 
                 bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb)
            
            If ((ndim == 3).Or.(l2p5d == 1)) Then
               bedge_facez_x(:,1+nguard0:nxb+nguard0,                  & 
                    1+nguard0*k2d,1,lb)=                               & 
                    bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)
               
               If (ndim == 3)                                          & 
                    bedge_facez_x(:,1+nguard0:nxb+nguard0,             & 
                    1+nguard0*k2d,2,lb)=                               & 
                    bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)
            End If
            
         ElseIf (jf == 4) Then

            bedge_facex_z(:,1,k2d+nyb+nguard0*k2d,klo:kup,lb) =        & 
                 bedge_facey_z(:,1+nguard0,2,klo:kup,lb)
            
            bedge_facex_z(:,2,k2d+nyb+nguard0*k2d,klo:kup,lb) =        & 
                 bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb)
            
            If ((ndim == 3).Or.(l2p5d == 1)) Then
               bedge_facez_x(:,1+nguard0:nxb+nguard0,                  & 
                    k2d+nyb+nguard0*k2d,                               & 
                    1,lb)=                                             & 
                    bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)
               
               If (ndim == 3)                                          & 
                    bedge_facez_x(:,1+nguard0:nxb+nguard0,             & 
                    k2d+nyb+nguard0*k2d,                               & 
                    2,lb)=                                             & 
                    bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)
            End If
            
         ElseIf (jf == 5) Then

            bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)=           & 
                 bedge_facez_x(:,1+nguard0:nxb+nguard0,                & 
                 1+nguard0*k2d,                                        & 
                 1,lb)
            
            bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)=           & 
                 bedge_facez_x(:,1+nguard0:nxb+nguard0,                & 
                 k2d+nyb+nguard0*k2d,1,lb)
            
            bedge_facex_y(:,1,1+nguard0*k2d:nyb+nguard0*k2d,           & 
                 klo,lb)=                                              & 
                 bedge_facez_y(:,1+nguard0,                            & 
                 1+nguard0*k2d:nyb+nguard0*k2d,                        & 
                 1,lb)
            
            bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,           & 
                 klo,lb)=                                              & 
                 bedge_facez_y(:,1+nxb+nguard0,                        & 
                 1+nguard0*k2d:nyb+nguard0*k2d                         & 
                 ,1,lb)
            
         ElseIf (jf == 6) Then
            
            bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)=          & 
                 bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0*k2d,  & 
                 2,lb)
            
            bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)=          & 
                 bedge_facez_x(:,1+nguard0:nxb+nguard0,                & 
                 k2d+nyb+nguard0*k2d,                                  & 
                 2,lb)
            
            bedge_facex_y(:,1,1+nguard0*k2d:nyb+nguard0*k2d,           & 
                 kup1,lb)=                                             & 
                 bedge_facez_y(:,1+nguard0,                            & 
                 1+nguard0*k2d:nyb+nguard0*k2d,2,lb)
            
            bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,           & 
                 kup1,lb)=                                             & 
                 bedge_facez_y(:,1+nxb+nguard0,                        & 
                 1+nguard0*k2d:nyb+nguard0*k2d,                        & 
                 2,lb)
            
        End If  ! End If (jf == 1)

        End If  ! End If (cnodetype == 2)

        End Do  ! End Do jf = 1,nfaces

      End If ! End If (nodetype(lb) == 1)
      End Do ! End Do lb = 1,lnblocks
      End If ! End If (lnblocks > 0)

      End If ! If (ndim >=2)

      Return
      End Subroutine amr_edge_average_udt
