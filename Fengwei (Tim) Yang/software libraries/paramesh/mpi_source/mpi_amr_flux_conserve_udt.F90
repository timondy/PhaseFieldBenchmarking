!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/amr_flux_conserve_udt
!! NAME
!!
!!   amr_flux_conserve_udt
!!
!! SYNOPSIS
!!
!!   call amr_flux_conserve_udt (mype)
!!   call amr_flux_conserve_udt (mype, flux_dir)
!!
!!   call amr_flux_conserve_udt (integer, optional integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   optional, integer, intent(in) :: flux_dir
!!     Option integer which selects which coordinate direction to apply
!!     the flux conservation operation to:
!!     If flux_dir = 1 -> x direction
!!        flux_dir = 2 -> y direction
!!        flux_dir = 3 -> z direction
!!     If this argument is not specified, then the default behaviour is
!!     to operate on all directions.  Using this argument can be useful
!!     for Strang-split schemes to improve performance.
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
!!
!! CALLS
!! 
!!   amr_restrict_bnd_data
!!   mpi_amr_comm_setup
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit the fluxes stored in the arrays
!!   flux_x, flux_y, or flux_z are corrected at jumps in refinement.  This
!!   is either an averaging proceedure or a sum as selected by the user
!!   by adjusting the logical variables which control this behaviour
!!   that are read at runtime from the file 'amr_runtime_parameters'.
!!
!! DESCRIPTION
!!
!!   This routine gets block boundary data from neighbors who are
!!   parents of leaf blocks. This is required in flux conserving schemes
!!   where the coarser block needs to use the same fluxes and mean pressures
!!   as will be used on the finer blocks across their shared boundary.
!!
!!   This version is called when uniform timesteps are being used across
!!   the blocks in the computation.
!!
!! NOTES
!!  
!!   This routine is NOT user callable.  The user interacts with this routine
!!   by calling 'amr_flux_conserve' which, in turn, calls this routine.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson for 
!!   directional flux correction (i.e. flux_dir option).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_flux_conserve_udt(mype,flux_dir)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_interfaces, only : amr_restrict_bnd_data,           & 
                                      amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_amr_comm_setup,          & 
                                          mpi_set_message_limits

      Implicit None

!-----Include statements.
      include 'mpif.h'

!-----Input/Output arguments.
      integer, optional, intent(in)  ::  flux_dir
      integer, intent(in)  ::  mype

!-----Local arrays and variables.
      Integer :: remote_pe,remote_block
      Integer :: remote_pe2,remote_block2
      Integer,Save :: anodetype(1)
      Integer,Save :: cnodetype
      Integer :: tag_offset,nprocs,ierr
      Integer :: iopt, lb, jf, iblk
      Integer :: face_min, face_max, flux_dirt
      Integer :: i, j, k, ia0, ib0, ja0, jb0, ka0, kb0
      Integer :: dtype, vtype, index, index0
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: n
      Logical :: lfound
      Logical :: lcc ,lfc,lec,lnc
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree

!-----Begin executable code.
      If (present(flux_dir)) Then
         flux_dirt = flux_dir
      Else
         flux_dirt = 0
      End If

      If (flux_dirt == 1) Then
         face_min = 1
         face_max = 2
      ElseIf (flux_dirt == 2) Then
         face_min = 3
         face_max = 4
      ElseIf (flux_dirt == 3) Then
         face_min = 5
         face_max = 6
      Else
         face_min = 1
         face_max = nfaces
      End If

      lcc = .False.
      lfc = .False.
      lec = .False.
      lnc = .False.
      iopt = 1

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      tag_offset = 100

!-----Note, both lflux and lrestrict are true so that the fluxes
!-----are acquired which are needed in the restriction operation.
      lguard    = .False.
      lprolong  = .False.
      lflux     = .True.
      ledge     = .False.
      lrestrict = .True.
      lfulltree = .False.
      Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,             & 
                              lflux,ledge,lrestrict,lfulltree,         & 
                              iopt,lcc,lfc,lec,lnc,tag_offset,         & 
                              flux_dir=flux_dirt)

!-----all leaf blocks provide reduced boundary data to their parents
      Call amr_restrict_bnd_data(mype,flux_dirt)

      tag_offset = 100

      lguard    = .False.
      lprolong  = .False.
      lflux     = .True.
      ledge     = .False.
      lrestrict = .False.
      lfulltree = .False.
      Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,             & 
                              lflux,ledge,lrestrict,lfulltree,         & 
                              iopt,lcc,lfc,lec,lnc,tag_offset,         & 
                              flux_dir=flux_dirt)

!-----cycle through the grid blocks on this processor
      If (lnblocks > 0) Then
      Do lb = 1,lnblocks

!-----Is this a leaf block and not at the original refinement level ?
      If (nodetype(lb) == 1) Then

!------Cycle over the blocks faces
       Do jf = face_min, face_max

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)
          remote_block2 = remote_block
          remote_pe2    = remote_pe
          cnodetype = 0

          If (remote_block > 0) Then

!---------(remote_block,remote_pe) may be a local block, a remote block,
!---------or it may not exist.
!---------If it is a local block Then check its nodetype.
!---------If it is found in the list of remote blocks stored in buffer space
!---------Then check its nodetype.
!---------If it is not found in either of these places, Then set its nodetype
!---------to 0.
          If (remote_pe2.ne.mype) Then
            lfound = .False.
            Do iblk = strt_buffer,last_buffer
              If (remote_block2 == laddress(1,iblk).And.               & 
                  remote_pe2  == laddress(2,iblk) ) Then
                remote_block2 = iblk
                remote_pe2    = mype
                lfound = .True.
              End If
            End Do
          ElseIf (remote_pe2 == mype) Then
            lfound = .True.
          End If

! Is the neighbor to this face a parent of a leaf block?
          If (lfound) Then
            cnodetype = nodetype(remote_block2)
          End If

          End If  ! End If (remote_pe2.ne.mype)

          If (cnodetype == 2) Then

!-----------If yes then copy the appropriate layer from its boundary 
!-----------variable data 
            If (remote_pe == mype .And. remote_block <= lnblocks) Then

            If (jf == 1) Then
               flux_x(1:nfluxes,1,:,:,lb) =                            &
                 flux_x(1:nfluxes,2,:,:,remote_block)
            ElseIf (jf == 2) Then
               flux_x(1:nfluxes,2,:,:,lb) =                            &
                 flux_x(1:nfluxes,1,:,:,remote_block)
            ElseIf (jf == 3) Then
               flux_y(1:nfluxes,:,1,:,lb) =                            &
                 flux_y(1:nfluxes,:,2,:,remote_block)
            ElseIf (jf == 4) Then
               flux_y(1:nfluxes,:,2,:,lb) =                            &
                 flux_y(1:nfluxes,:,1,:,remote_block)
            ElseIf (jf == 5) Then
               flux_z(1:nfluxes,:,:,1,lb) =                            &
                 flux_z(1:nfluxes,:,:,2,remote_block)
            ElseIf (jf == 6) Then
               flux_z(1:nfluxes,:,:,2,lb) =                            &
                 flux_z(1:nfluxes,:,:,1,remote_block)
            End If  ! End If (jf == 1)

            Else ! If (remote_pe

            Call amr_mpi_find_blk_in_buffer(mype,remote_block,         & 
                                            remote_pe,1,dtype,index0,  &
                                            lfound)
            vtype = 1
            Call mpi_set_message_limits(dtype,                         & 
                                        ia0,ib0,ja0,jb0,ka0,kb0,vtype)

            index = index0 + 1

            If (dtype == 13.Or.dtype == 15.Or.dtype == 14) Then

            ia = ia0
            ib = ib0
            ja = ja0
            jb = jb0
            ka = ka0
            kb = kb0

            If (dtype == 13) Then
               ia = 1
               ib = 1
            ElseIf (dtype == 15) Then
               ia = 2
               ib = 2
            ElseIf (dtype == 14) Then
               ia = 1
               ib = 2
            End If

            If (jf == 1 .Or. jf == 2) Then

            Do k = ka,kb
               Do j = ja,jb
                  Do i = ia,ib
                     Do n=1,nfluxes
                        recvarxf(n,i,j,k) = temprecv_buf(index)
                        index  = index + 1
                     End Do
                  End Do
               End Do
            End Do

            Else

               If (flux_dirt == 0) Then
               index = index + nfluxes*(ib-ia+1)*(jb-ja+1)*(kb-ka+1)
               End If

            End If

            End If  ! End If (dtype == 13.Or.dtype == 15.Or.dtype == 14)

            If (ndim >= 2) Then
               If (dtype == 11.Or.dtype == 17.Or.dtype == 14) Then

                  ia = ia0
                  ib = ib0
                  ja = ja0
                  jb = jb0
                  ka = ka0
                  kb = kb0
                  
                  If (dtype == 11) Then
                     ja = 1
                     jb = 1
                  ElseIf (dtype == 17) Then
                     ja = 2
                     jb = 2
                  ElseIf (dtype == 14) Then
                     ja = 1
                     jb = 2
                  End If

                  If (jf == 3 .Or. jf == 4) Then

                  Do k = ka,kb
                     Do j = ja,jb
                        Do i = ia,ib
                           Do n=1,nfluxes
                              recvaryf(n,i,j,k) =                      & 
                                   temprecv_buf(index)
                              index  = index + 1
                           End Do
                        End Do
                     End Do
                  End Do
                  
                  Else
                     
                  If (flux_dirt == 0) Then
                  index = index + nfluxes*(ib-ia+1)*(jb-ja+1)*         & 
                                          (kb-ka+1)
                  End If

                  End If  ! End If (jf == 3 .Or. jf == 4)

               End If  ! End If (dtype == 11.Or.dtype == 17.Or.dtype == 14)
            End If  ! End If (ndim >= 2)

            If (ndim == 3) Then
               If (dtype == 5.Or.dtype == 23.Or.dtype == 14) Then
                  
                  ia = ia0
                  ib = ib0
                  ja = ja0
                  jb = jb0
                  ka = ka0
                  kb = kb0
                  
                  If (dtype == 5) Then
                     ka = 1
                     kb = 1
                  ElseIf (dtype == 23) Then
                     ka = 2
                     kb = 2
                  ElseIf (dtype == 14) Then
                     ka = 1
                     kb = 2
                  End If
                  
                  If (jf == 5 .Or. jf == 6) Then

                  Do k = ka,kb
                     Do j = ja,jb
                        Do i = ia,ib
                           Do n=1,nfluxes
                              recvarzf(n,i,j,k) =                      & 
                                   temprecv_buf(index)
                              index  = index + 1
                           End Do
                        End Do
                     End Do
                  End Do
                  
                  Else

                  If (flux_dirt == 0) Then
                  index = index + nfluxes*(ib-ia+1)*(jb-ja+1)*         & 
                                          (kb-ka+1)
                  End If

                  End If  ! End If (jf == 5 .Or. jf == 6)

               End If  ! End If (dtype == 5.Or.dtype == 23.Or.dtype == 14)
            End If  ! End If (ndim == 3)

            If (jf == 1) Then
               flux_x(1:nfluxes,1,:,:,lb) = recvarxf(1:nfluxes,2,:,:)
            ElseIf (jf == 2) Then
               flux_x(1:nfluxes,2,:,:,lb) = recvarxf(1:nfluxes,1,:,:)
            ElseIf (jf == 3) Then
               flux_y(1:nfluxes,:,1,:,lb) = recvaryf(1:nfluxes,:,2,:)
            ElseIf (jf == 4) Then
               flux_y(1:nfluxes,:,2,:,lb) = recvaryf(1:nfluxes,:,1,:)
            ElseIf (jf == 5) Then
               flux_z(1:nfluxes,:,:,1,lb) = recvarzf(1:nfluxes,:,:,2)
            ElseIf (jf == 6) Then
               flux_z(1:nfluxes,:,:,2,lb) = recvarzf(1:nfluxes,:,:,1)
            End If

            End If ! End If (remote_pe == mype .And. remote_block <= lnblocks)

          End If  ! End If (cnodetype == 2)

       End Do  ! End Do jf = face_min, face_max

      End If  ! End If (nodetype(lb) == 1)
      End Do  ! End Do lb = 1,lnblocks  
      End If  ! End If (lnblocks > 0)

      Return
      End Subroutine amr_flux_conserve_udt
