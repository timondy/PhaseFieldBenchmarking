!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_edge_diagonal_check
!! NAME
!!
!!   amr_edge_diagonal_check (mype)
!!
!! SYNOPSIS
!!
!!   Call amr_edge_diagonal_check (mype)
!!
!!   Call amr_edge_diagonal_check (integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype
!!     The Calling processors id.
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
!!   amr_mpi_find_blk_in_buffer
!!   mpi_set_message_limits
!!   mpi_put_edge_buffer_1blk
!!
!! RETURNS
!!
!!   Nothing returned. 
!!
!! DESCRIPTION
!! 
!!  This routine checks to see if the diagonal block between two
!!  leaf-neighbors at the same refinement level as the current block,
!!  is refined. If it is then the edge-based variables along the edge
!!  shared with that diagonal block is given the edge values
!!  form the refined diagonal block, to insure conservation properties.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (October 1997).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_edge_diagonal_check(mype)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_set_message_limits,      & 
                                          mpi_put_edge_buffer_1blk

      Implicit None 

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, intent(in)  ::  mype

!-----Local arrays and variables.
      Integer :: nguard0
      Integer :: klo,kup
      Integer :: jlo,jup
      Integer :: ilo,iup
      Integer :: remote_pe,remote_block
      Integer :: remote_pe2,remote_block2
!      Integer :: mark_edge(12,maxblocks)
      Integer,allocatable :: mark_edge(:,:)
      Integer :: i, ie, iblk, lb, k, j, ierrcode, ierr
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: ia0, ib0, ja0, jb0, ka0, kb0
      Integer :: dtype, vtype
      Integer :: index, index0, n
      Real,Allocatable :: receive(:,:)
      Logical :: lfound

!-----Begin executable code.
! CEG allocate memory
  allocate(mark_edge(12,maxblocks))
      nguard0 = nguard*npgs
      klo=1+nguard0*k3d
      kup=klo*k3d+nzb
      jlo=1+nguard0*k2d
      jup=jlo*k2d+nyb
      ilo=1+nguard0
      iup=ilo+nxb

      Allocate(receive(nedges,maxdim+2*nguard0))

      If (ndim >= 2) Then

!-----Initialize array marking edges for diagonal patching.
      mark_edge(:,:) = 0

      Do i = 1,no_of_diagonal_edges
        ie   = edge_mark(6,1,i)
        iblk = edge_mark(6,2,i)
        mark_edge(ie,iblk) = i
      End Do

!-----Loop over the blocks on this processor.
      If (lnblocks > 0) Then
      Do lb=1,lnblocks

!-----Is this a leaf block which has finished its current timestep?
      If ((nodetype(lb) == 1.And..not.var_dt) .Or.                     & 
          (nodetype(lb) == 1.And.ldtcomplete(lb))) Then

!-----Any edges on this block which are still marked need a diagonal patch.
!-----Note that in the data copies below, we can always assume that a
!-----neighbor block exists, since the edge would not have been marked
!-----earlier if that was not so.

!------Loop over the edges on this block.
       Do ie=1,nbedges

       If (mark_edge(ie,lb) >= 1) Then

        lfound = .False.
        remote_block  = edge_mark(6,3,mark_edge(ie,lb))
        remote_pe     = edge_mark(6,4,mark_edge(ie,lb))
        remote_block2 = edge_mark(6,3,mark_edge(ie,lb))
        remote_pe2    = edge_mark(6,4,mark_edge(ie,lb))

!--------(remote_block,remote_pe) may be a local block or a remote block.
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

!-----The edge data on the neighboring faces can be assumed to have been 
!-----averaged correctly from the refined diagonal blocks.
      If (remote_pe == mype .And. remote_block <= lnblocks) Then

!------Now copy over the edge data from one of the neighbors.
       If (ie == 1) Then                    ! edge: x low edge, y low edge
         Do k=klo,kup-k3d
           bedge_facex_z(:,1,1+nguard0*k2d,k,lb) =                     &
             bedge_facex_z(:,2,jup,k,remote_block)
         End Do
         bedge_facey_z(:,1+nguard0,1,klo:kup-k3d,lb)=                  & 
                      bedge_facex_z(:,1,1+nguard0*k2d,klo:kup-k3d,lb)

       ElseIf (ie == 2) Then               ! edge: x low edge, y high edge
         Do k=klo,kup-k3d
           bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,k,lb) =               & 
              bedge_facex_z(:,2,jlo,k,remote_block)
         End Do
         bedge_facey_z(:,1+nguard0,2,klo:kup-k3d,lb)=                  & 
           bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)

       ElseIf (ie == 3) Then               ! edge: x high edge, y low edge
         Do k=klo,kup-k3d
           bedge_facex_z(:,2,1+nguard0*k2d,k,lb) =                     &
              bedge_facex_z(:,1,jup,k,remote_block)
         End Do
         bedge_facey_z(:,1+nguard0+nxb,1,klo:kup-k3d,lb)=              & 
           bedge_facex_z(:,2,1+nguard0*k2d,klo:kup-k3d,lb)

       ElseIf (ie == 4) Then               ! edge: x high edge, y high edge
         Do k=klo,kup-k3d
           bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,k,lb) =               & 
              bedge_facex_z(:,1,jlo,k,remote_block)
         End Do
         bedge_facey_z(:,1+nguard0+nxb,2,klo:kup-k3d,lb)=              & 
           bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)

       ElseIf (ie == 5) Then                ! edge: y low edge, z low edge
         Do i=ilo,iup-1
           bedge_facey_x(:,i,1,klo,lb) =                               &
              bedge_facey_x(:,i,2,kup,remote_block)
         End Do
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k3d,1,lb)=                & 
                      bedge_facey_x(:,ilo:iup-1,1,klo,lb)

       ElseIf (ie == 6) Then                ! edge: y high edge, z low edge
         Do i=ilo,iup-1
           bedge_facey_x(:,i,2,klo,lb) =                               &
              bedge_facey_x(:,i,1,kup,remote_block)
         End Do
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,1,lb)=          & 
                      bedge_facey_x(:,ilo:iup-1,2,klo,lb)

       ElseIf (ie == 7) Then                ! edge: y low edge, z high edge
         Do i=ilo,iup-1
           bedge_facey_x(:,i,1,kup,lb) =                               &
              bedge_facey_x(:,i,2,klo,remote_block)
         End Do
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k2d,2,lb)=                & 
                      bedge_facey_x(:,ilo:iup-1,1,kup,lb)

       ElseIf (ie == 8) Then                ! edge: y high edge, z high edge
         Do i=ilo,iup-1
           bedge_facey_x(:,i,2,kup,lb) =                               &
              bedge_facey_x(:,i,1,klo,remote_block)
         End Do
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,2,lb)=          & 
                      bedge_facey_x(:,ilo:iup-1,2,kup,lb)

       ElseIf (ie == 9) Then                ! edge: x low edge, z low edge
         Do j=jlo,jup-k2d
           bedge_facex_y(:,1,j,klo,lb) =                               &
              bedge_facex_y(:,2,j,kup,remote_block)
         End Do
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)=                  & 
                      bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)

       ElseIf (ie == 10) Then                ! edge: x low edge, z high edge
         Do j=jlo,jup-k2d
           bedge_facex_y(:,1,j,kup,lb) =                               &
              bedge_facex_y(:,2,j,klo,remote_block)
         End Do
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)=                  & 
                      bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)

       ElseIf (ie == 11) Then                ! edge: x high edge, z low edge
         Do j=jlo,jup-k2d
           bedge_facex_y(:,2,j,klo,lb) =                               &
              bedge_facex_y(:,1,j,kup,remote_block)
         End Do
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)=              & 
                      bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)

       ElseIf (ie == 12) Then                ! edge: x high edge, z high edge
         Do j=jlo,jup-k2d
           bedge_facex_y(:,2,j,kup,lb) =                               &
              bedge_facex_y(:,1,j,klo,remote_block)
         End Do
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)=              & 
                      bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)
       End If

      Else                      ! If (remote_pe

         Call mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)

         If (ie == 1) Then       ! edge: x low edge, y low edge
            Do k=klo,kup-k3d
               bedge_facex_z(:,1,jlo,k,lb)=                            & 
                    recvarx2e(:,2,jup,k)
            End Do
            bedge_facey_z(:,ilo,1,klo:kup-k3d,lb)=                     & 
                 bedge_facex_z(:,1,jlo,klo:kup-k3d,lb)
            
         ElseIf (ie == 2) Then   ! edge: x low edge, y high edge
            Do k=klo,kup-k3d
               bedge_facex_z(:,1,jup,k,lb)=                            & 
                    recvarx2e(:,2,jlo,k)
            End Do
            bedge_facey_z(:,ilo,2,klo:kup-k3d,lb)=                     & 
                 bedge_facex_z(:,1,jup,klo:kup-k3d,lb)
            
         ElseIf (ie == 3) Then   ! edge: x high edge, y low edge
            Do k=klo,kup-k3d
               bedge_facex_z(:,2,jlo,k,lb)=                            & 
                    recvarx2e(:,1,jup,k)
            End Do
            bedge_facey_z(:,iup,1,klo:kup-k3d,lb)=                     & 
                 bedge_facex_z(:,2,jlo,klo:kup-k3d,lb)
            
         ElseIf (ie == 4) Then   ! edge: x high edge, y high edge
            Do k=klo,kup-k3d
               bedge_facex_z(:,2,jup,k,lb)=                            & 
                    recvarx2e(:,1,jlo,k)
            End Do
            bedge_facey_z(:,iup,2,klo:kup-k3d,lb)=                     & 
                 bedge_facex_z(:,2,jup,klo:kup-k3d,lb)
         ElseIf (ie == 5) Then   ! edge: y low edge, z low edge
            Do i=ilo,iup-1
               bedge_facey_x(:,i,1,klo,lb)= recvary1e(:,i,2,kup)
            End Do
            bedge_facez_x(:,ilo:iup-1,jlo,1,lb)=                       & 
                 bedge_facey_x(:,ilo:iup-1,1,klo,lb)
            
         ElseIf (ie == 6) Then   ! edge: y high edge, z low edge
            Do i=ilo,iup-1
               bedge_facey_x(:,i,2,klo,lb)= recvary1e(:,i,1,kup)
            End Do
            bedge_facez_x(:,ilo:iup-1,jup,1,lb)=                       & 
                 bedge_facey_x(:,ilo:iup-1,2,klo,lb)
            
         ElseIf (ie == 7) Then   ! edge: y low edge, z high edge
            Do i=ilo,iup-1
               bedge_facey_x(:,i,1,kup,lb)= recvary1e(:,i,2,klo)
            End Do
            bedge_facez_x(:,ilo:iup-1,jlo,2,lb)=                       & 
                 bedge_facey_x(:,ilo:iup-1,1,kup,lb)
            
         ElseIf (ie == 8) Then   ! edge: y high edge, z high edge
            Do i=ilo,iup-1
               bedge_facey_x(:,i,2,kup,lb)= recvary1e(:,i,1,klo)
            End Do
            bedge_facez_x(:,ilo:iup-1,jup,2,lb)=                       & 
                 bedge_facey_x(:,ilo:iup-1,2,kup,lb)
            
         ElseIf (ie == 9) Then   ! edge: x low edge, z low edge
            Do j=jlo,jup-k2d
               bedge_facex_y(:,1,j,klo,lb)= recvarx1e(:,2,j,kup)
            End Do
            bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)=               & 
                 bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)
            
         ElseIf (ie == 10) Then  ! edge: x low edge, z high edge
            Do j=jlo,jup-k2d
               bedge_facex_y(:,1,j,kup,lb)= recvarx1e(:,2,j,klo)
            End Do
            bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)=               & 
                 bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)
            
         ElseIf (ie == 11) Then  ! edge: x high edge, z low edge
            Do j=jlo,jup-k2d
               bedge_facex_y(:,2,j,klo,lb)= recvarx1e(:,1,j,kup)
            End Do
            bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)=           & 
                 bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)
            
         ElseIf (ie == 12) Then  ! edge: x high edge, z high edge
            Do j=jlo,jup-k2d
               bedge_facex_y(:,2,j,kup,lb)= recvarx1e(:,1,j,klo)
            End Do
            bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)=           & 
                 bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)
            
         End If  ! End If (ie == 1)

      End If  ! End If (remote_pe == mype .And. remote_block <= lnblocks)
      End If  ! End If (mark_edge(ie,lb) >= 1)
      End Do  ! End Do ie=1,nbedges
      End If  ! End If If ((nodetype(lb) == 1.And..not.var_dt) .Or. ...

      End Do  ! End Do lb = 1, lnblocks
      End If  ! End If (lnblocks > 0)

      End If  ! If (ndim >= 2)

      Deallocate(receive)
! CEG deallocate memory
  deallocate(mark_edge)

      Return
      End Subroutine amr_edge_diagonal_check
