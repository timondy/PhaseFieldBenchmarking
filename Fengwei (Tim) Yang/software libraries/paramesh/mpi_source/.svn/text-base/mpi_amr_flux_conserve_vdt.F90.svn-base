!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/amr_flux_conserve_vdt
!! NAME
!!
!!   amr_flux_conserve_vdt
!!
!! SYNOPSIS
!!
!!   call amr_flux_conserve_vdt (mype, nsub)
!!
!!   call amr_flux_conserve_vdt (integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype          
!!     The calling processor number.
!!
!!   integer, intent(in) :: nsub          
!!     The current time subcycle. If this is 1 then this info is used to 
!!     reset the temporary boundary flux arrays to 0.
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
!!   amr_restrict_bnd_data_vdt
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
!!   The data structure used to store and pass this data is defined
!!   in the include file 'block_boundary_data.fh' which can be included
!!   in 'physicaldata.fh'.
!!
!!   This version is used when variable timesteps are allowed across the
!!   blocks in the computation.
!!
!! NOTES
!!  
!!   This routine is NOT user callable.  The user interacts with this routine
!!   by calling 'amr_flux_conserve' which, in turn, calls this routine.
!!
!!   THIS ROUTINE HAS NOT BEEN TESTED - July 14,2000
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) 
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_flux_conserve_vdt(mype,nsub)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_interfaces, only : amr_restrict_bnd_data_vdt,       & 
                                      amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpi_amr_comm_setup,          & 
                                          mpi_set_message_limits

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output statements.
      Integer, intent(in)  ::  mype,nsub

!-----Local arrays and variables.
      Integer :: remote_pe,remote_block
      Integer :: remote_pe2,remote_block2
      Integer,Save ::  anodetype(1)
      Integer,Save :: cnodetype
      Integer :: tag_offset,nprocs,ierr
      Integer :: iopt
      Integer :: lb, lcycle, jf, iblk, iface
      Integer :: i, j, k, ia0, ib0, ja0, jb0, ka0, kb0
      Integer :: dtype, vtype, index, index0
      Integer :: ia, ib, ja, jb, ka, kb
      Integer :: n
      Real    :: phase0, phase1
      Logical :: lfound
      Logical :: lcc ,lfc,lec,lnc
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree

!-----Begin executable code.
      If (var_dt) Then

      lcc = .False.
      lfc = .False.
      lec = .False.
      lnc = .False.
      iopt = 1

      If (lnblocks > 0) Then
      Do lb = 1,lnblocks

!-----Is this a parent of at least one leaf block ?
      If (nodetype(lb) == 2) Then

!-------Set timestep phases for the current block, and for the next finer level.
        lcycle = loc_cycle(lrefine(lb))
        phase0 = phase_dt(lrefine(lb))
        phase1 = phase_dt(lrefine(lb)+1)

!-------At start of the current blocks timestep zero out the arrays used to 
!-------accumulate boundary fluxes from its children.
        If (lcycle == 1) Then
           ttflux_x(:,:,:,:,lb) = 0.
           If (ndim >= 2) ttflux_y(:,:,:,:,lb) = 0.
           If (ndim == 3) ttflux_z(:,:,:,:,lb) = 0.
        End If

      End If  ! End If (nodetype(lb) == 2)
      End Do  ! End Do lb = 1, lnblocks
      End If  ! End If (lnblocks > 0)

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      tag_offset = 100

!-----Note, both lflux and lrestrict are true so that the fluxes
!-----are acquired which are needed in the restriction operation.
      lguard    = .False.
      lprolong  = .False.
      lflux     = .True.
      ledge     = .False.
      lrestrict = .True.
      lrestrict = .False.
      Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,             & 
                              lflux,ledge,lrestrict,lfulltree,         & 
                              iopt,lcc,lfc,lec,lnc,tag_offset)

!-----Leaf blocks which have completed their timestep provide reduced 
!-----boundary data to their parents.
!-----Fluxes are accumulated in the ttflux_ arrays.
      Call amr_restrict_bnd_data_vdt(mype)

!-----Parents who have completed their timestep and border a leaf block
!-----update their fluxes.
      Do lb = 1,lnblocks

!-----Is this a parent block of at least one leaf node?
      If ((nodetype(lb) == 2).and.ldtcomplete(lb)) Then

!-------If yes Then cycle through its neighbors.
        Do iface=1,nfaces

!---------If this neighbor is a leaf block or an external boundary Then 
!---------replace fluxes with restricted fluxes.
          cnodetype = 1
          If (neigh(1,iface,lb) >= 1) Then

          remote_pe    = neigh(2,iface,lb)
          remote_block = neigh(1,iface,lb)

!---------If (remote_block,remote_pe) is not a local block Then it must have a
!---------local copy available in the buffer space at the end of the local
!---------block list.
          If (remote_pe.ne.mype) Then
            Do iblk = strt_buffer,last_buffer
              If (remote_block == laddress(1,iblk).and.                & 
                  remote_pe  == laddress(2,iblk) ) Then
                remote_block = iblk
                remote_pe    = mype
              End If
            End Do
          End If

          cnodetype = nodetype(remote_block)

          End If  ! End If (neigh(1,iface,lb) >= 1)

          If (cnodetype == 1) Then
            If (iface == 1) flux_x(:,1,:,:,lb)=ttflux_x(:,1,:,:,lb)
            If (iface == 2) flux_x(:,2,:,:,lb)=ttflux_x(:,2,:,:,lb)
            If (iface == 3) flux_y(:,:,1,:,lb)=ttflux_y(:,:,1,:,lb)
            If (iface == 4) flux_y(:,:,2,:,lb)=ttflux_y(:,:,2,:,lb)
            If (iface == 5) flux_z(:,:,:,1,lb)=ttflux_z(:,:,:,1,lb)
            If (iface == 6) flux_z(:,:,:,2,lb)=ttflux_z(:,:,:,2,lb)
          End If

        End Do  ! End Do iface=1,nfaces
      End If  ! End If ((nodetype(lb) == 2).and.ldtcomplete(lb))
      End Do  ! End Do lb = 1, lnblocks

      tag_offset = 100

      lguard    = .False.
      lprolong  = .False.
      lflux     = .True.
      ledge     = .False.
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

!-----Has this block completed its timestep?
      If (ldtcomplete(lb)) Then

!------Cycle over the blocks faces
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
!---------Then check its nodetype.
!---------If it is not found in either of these places, Then set its nodetype
!---------to 0.
          If (remote_pe2.ne.mype) Then
             lfound = .False.
             Do iblk = strt_buffer,last_buffer
                If (remote_block2 == laddress(1,iblk).and.             & 
                    remote_pe2  == laddress(2,iblk) ) Then
                   remote_block2 = iblk
                   remote_pe2    = mype
                   lfound = .True.
                End If
             End Do
          ElseIf (remote_pe2 == mype) Then
             lfound = .True.
          End If

!---------Is the neighbor to this face a parent of a leaf block?
          If (lfound) Then
             cnodetype = nodetype(remote_block2)
          End If

          End If  ! End If (remote_block > 0)

          If (cnodetype == 2) Then

!-----------If yes then copy the appropriate layer from its boundary 
!-----------variable data 
            If (remote_pe == mype .and. remote_block <= lnblocks) Then

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
            End If  

            Else ! If (remote_pe

            Call amr_mpi_find_blk_in_buffer(mype,remote_block,         & 
                                            remote_pe,1,dtype,index0,  &
                                            lfound)
            vtype = 1
            Call mpi_set_message_limits(dtype,                         & 
                                        ia0,ib0,ja0,jb0,ka0,kb0,vtype)

            index = index0 + 1

            If (dtype == 13.or.dtype == 15.or.dtype == 14) Then

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

            If (jf == 1 .or. jf == 2) Then

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

            End If  ! End If (jf == 1 .or. jf == 2)

            End If  ! End If (dtype == 13.or.dtype == 15.or.dtype == 14)

            If (ndim >= 2) Then
               If (dtype == 11.or.dtype == 17.or.dtype == 14) Then

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

                  If (jf == 3 .or. jf == 4) Then

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

                  End If  ! End If (jf == 3 .or. jf == 4)

               End If  ! End If (dtype == 11.or.dtype == 17.or.dtype == 14)
            End If  ! End If (ndim >= 2)

            If (ndim == 3) Then
               If (dtype == 5.or.dtype == 23.or.dtype == 14) Then

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

                  If (jf == 5 .or. jf == 6) Then
                                       
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

                  End If  ! End If (jf == 5 .or. jf == 6)

               End If  ! End If (dtype == 5.or.dtype == 23.or.dtype == 14)
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

            End If  ! End If (remote_pe == mype .and. remote_block <= lnblocks)

          End If  ! End If (cnodetype == 2)

       End Do  ! End Do jf = 1,nfaces

      End If  ! End If (ldtcomplete(lb))

      End If  ! End If (nodetype(lb) == 1)
      End Do  ! End Do lb = 1, lnblocks
      End If  ! End If (lnblocks > 0)

      End If  ! End If (var_dt)

      Return
      End Subroutine amr_flux_conserve_vdt
