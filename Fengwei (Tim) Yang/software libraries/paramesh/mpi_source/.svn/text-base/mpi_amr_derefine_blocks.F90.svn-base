!---------------------------------------------------------------------- 
! PARAMESH - an adaptive mesh library.  
! Copyright (C) 2003 
! 
! Use of the PARAMESH software is governed by the terms of the 
! usage agreement which can be found in the file 
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory. 
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_derefine_blocks 
!! NAME 
!! 
!!   mpi_amr_derefine_blocks 
!! 
!! SYNOPSIS 
!! 
!!   Call mpi_amr_derefine_blocks (lnblocks_old, mype) 
!! 
!!   Call mpi_amr_derefine_blocks (integer, integer) 
!! 
!! ARGUMENTS 
!! 
!!   integer, intent(in) :: lnblocks_old, mype 
!!   lnbllocks_old -> The number of local blocks before derefinement. 
!!   mype -> The local processor id. 
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
!!   constants 
!! 
!! CALLS 
!! 
!! RETURNS 
!! 
!!   Nothing returned. 
!! 
!! DESCRIPTION 
!! 
!!   This routine handles derefinement to the mesh.  It removes blocks 
!!   which are marked for derinement. 
!! 
!! AUTHORS 
!! 
!!   Kevin Olson, 1997 
!! 
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz] !!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_derefine_blocks(lnblocks_old, mype)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use constants

      Implicit None

!-----Input/Output arguments.
      Integer, intent(inout) :: lnblocks_old
      Integer, intent(in) :: mype

!-----Include statements.
      Include 'mpif.h'

!-----Local arrays and variables
      Real :: eps,accuracy
!      Integer :: new_loc(maxblocks_tr)
      Integer,allocatable :: new_loc(:)
      Integer :: i,j,k,jsend
      Integer :: lnblocks2
! CEG rewritten to have all the arrays as allocated rather than on the stack
!      Integer :: neight(2,mfaces,maxblocks_tr)
!      Integer :: childt(2,mchild,maxblocks_tr)
!      Integer :: parentt(2,maxblocks_tr)
!      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
!      Integer :: reqr(maxblocks_tr)
!      Integer :: nodetype_chi(nchild,maxblocks_tr)
      Integer,allocatable :: neight(:,:,:)
      Integer,allocatable :: childt(:,:,:)
      Integer,allocatable :: parentt(:,:)
      Integer,allocatable :: statr(:,:)
      Integer,allocatable :: reqr(:)
      Integer,allocatable :: nodetype_chi(:,:)
      Integer :: ierr,nsend,nrecv
      Integer :: jr0
      Logical :: lsend

!-----Begin Executable code.

      accuracy = 10./10.**precision(accuracy)
      eps = accuracy

!-----remove blocks marked for derefinement by packing the data 
!  print *,maxblocks_tr
      allocate(new_loc(1:maxblocks_tr))
      allocate(neight(2, mfaces, 1:maxblocks_tr))
      allocate(childt(2, mchild, 1:maxblocks_tr))
      allocate(parentt(2, 1:maxblocks_tr))
      allocate(statr(MPI_STATUS_SIZE, 1:maxblocks_tr))
      allocate(reqr(1:maxblocks_tr))
      allocate(nodetype_chi(nchild, 1:maxblocks_tr))
!      print *, i, maxblocks_tr
!      Do i = 1, maxblocks_tr
!         new_loc(i) = -1
!      End Do

      new_loc(:) = -1

!-----Compute new_loc, new_loc marks where each block will end up after the !-----derefinement is Do ne
      k = 1
      Do i = 1,lnblocks
         If (.not.derefine(i)) Then
            new_loc(i) = k
            k = k + 1
          End If
      End Do

!  print *,'mpi_amr_derefine_blocks.F90', lnblocks, maxblocks_tr !-----reconnect all pointers
      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
             nrecv = nrecv + 1
             Call MPI_IRECV(parentt(1,i),1,MPI_INTEGER, &
                  parent(2,i),i,MPI_COMM_WORLD, &
                  reqr(nrecv),ierr)
           Else
             parentt(1,i) = new_loc(parent(1,i))
           End If
         End If
       End Do
       
       nsend = 0
       Do i = 1,lnblocks
         Do j = 1,nchild
           If (child(1,j,i) > 0) Then
             If (child(2,j,i).ne.mype) Then
               nsend = nsend + 1
               Call MPI_SSEND (new_loc(i),1,MPI_INTEGER, &
                    child(2,j,i),child(1,j,i),MPI_COMM_WORLD, &
                    ierr)
             End If
           End If
         End Do
       End Do

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      nrecv = 0
      Do i = 1,lnblocks
        Do j = 1,nchild
          If (child(1,j,i) > 0) Then
            If (child(2,j,i).ne.mype) Then
              nrecv = nrecv + 1
              Call MPI_IRECV(childt(1,j,i),1,MPI_INTEGER, &
                   child(2,j,i),child(1,j,i),MPI_COMM_WORLD, &
                   reqr(nrecv),ierr)
            Else
              childt(1,j,i) = new_loc(child(1,j,i))
            End If
          End If
        End Do
       End Do
       
       nsend = 0
       Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
             nsend = nsend + 1
             Call MPI_SSEND (new_loc(i),1,MPI_INTEGER, &
                  parent(2,i),i,MPI_COMM_WORLD, &
                  ierr)
           End If
         End If
       End Do

       If (nrecv > 0) Then
         Call MPI_WAITALL(nrecv,reqr,statr,ierr)
       End If

       Do j = 1,nfaces

         If (mod(j,2) == 0) Then
            jsend = j - 1
         Else
            jsend = j + 1
         End If
            
         nrecv = 0
         Do i = 1,lnblocks
            If (neigh(1,j,i) > 0) Then
               If (neigh(2,j,i).ne.mype) Then
                  nrecv = nrecv + 1
                  Call MPI_IRECV(neight(1,j,i),1,MPI_INTEGER, &
                       neigh(2,j,i),neigh(1,j,i),MPI_COMM_WORLD, &
                       reqr(nrecv),ierr)
               Else
                  neight(1,j,i) = new_loc(neigh(1,j,i))
               End If
            End If
         End Do
      
         nsend = 0
         Do i = 1,lnblocks

         lsend=.True.

         If (spherical_pm) Then
          If (j == 3.and.abs(bnd_box(2,2,i)-pi) < eps) lsend=.False.
          If (j == 4.and.abs(bnd_box(1,2,i)) < eps) lsend=.False.
         End If

         If (lsend) Then
           If (neigh(1,jsend,i) > 0) Then
             If (neigh(2,jsend,i).ne.mype) Then
               nsend = nsend + 1
               Call MPI_SSEND (new_loc(i),1,MPI_INTEGER, &
                    neigh(2,jsend,i),i,MPI_COMM_WORLD, &
                    ierr)
             End If
           End If
         End If

         If (spherical_pm) Then
           lsend = .True.
           jr0 = jsend
           If (j == 3.and.abs(bnd_box(1,2,i)) < eps) jr0 = 3
           If (j == 4.and.abs(bnd_box(1,2,i)) < eps) lsend = .False.
           If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) jr0 = 4
           If (j == 3.and.abs(bnd_box(2,2,i)-pi) < eps) lsend = .False.
           If (abs(bnd_box(1,2,i)) < eps.and.  &
               abs(bnd_box(2,2,i)-pi) < eps) Then
             write(*,*) 'both poles in blk ',i,abs(bnd_box(1,2,i)), &
               abs(bnd_box(2,2,i)-pi)
             lsend = .True.
           End If

           If (lsend.and.jr0 == j) Then
             If (neigh(1,jr0,i) > 0) Then
               If (neigh(2,jr0,i).ne.mype) Then
                  nsend = nsend + 1
                  Call MPI_SSEND (new_loc(i),1,MPI_INTEGER, &
                       neigh(2,jr0,i),i,MPI_COMM_WORLD, &
                       ierr)
               End If
             End If
            End If
         End If ! End If (spherical_pm)

         End Do ! End Do i = 1,lnblocks

         If (nrecv > 0) Then
            Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         End If

      End Do ! End Do j = 1,nfaces

      Do i = 1,lnblocks_old
        If (new_loc(i).ne.i.and.new_loc(i) > 0) Then
          If (nvar > 0) unk(:,:,:,:,new_loc(i)) = unk(:,:,:,:,i)
          If (nfacevar > 0) Then
             facevarx(:,:,:,:,new_loc(i)) = facevarx(:,:,:,:,i)
             facevary(:,:,:,:,new_loc(i)) = facevary(:,:,:,:,i)
             facevarz(:,:,:,:,new_loc(i)) = facevarz(:,:,:,:,i)
          End If
          If (nvaredge > 0) Then
             unk_e_x(:,:,:,:,new_loc(i)) = unk_e_x(:,:,:,:,i)
             unk_e_y(:,:,:,:,new_loc(i)) = unk_e_y(:,:,:,:,i)
             unk_e_z(:,:,:,:,new_loc(i)) = unk_e_z(:,:,:,:,i)
          End If
          If (nvarcorn > 0) unk_n(:,:,:,:,new_loc(i)) = &
                                             unk_n(:,:,:,:,i)
        End If
      End Do

      parent(1,1:lnblocks) = parentt(1,1:lnblocks)
      child(1,:,1:lnblocks) = childt(1,:,1:lnblocks)
      neigh(1,:,1:lnblocks) = neight(1,:,1:lnblocks)

      k = 1
      lnblocks2 = lnblocks
      Do i = 1,lnblocks
         
         If (.not.derefine(i)) Then
            
            If (k.ne.i) Then
               Do j = 1,nchild
                  child(1,j,k) = child(1,j,i)
                  child(2,j,k) = child(2,j,i)
               End Do
               parent(1,k) = parent(1,i)
               parent(2,k) = parent(2,i)
               Do j = 1,nfaces
                  neigh(1,j,k) = neigh(1,j,i)
                  neigh(2,j,k) = neigh(2,j,i)
               End Do
               Do j = 1,ndim
                  coord(j,k) = coord(j,i)
                  bnd_box(1,j,k) = bnd_box(1,j,i)
                  bnd_box(2,j,k) = bnd_box(2,j,i)
               End Do
               bsize(:,k) = bsize(:,i)
               newchild(k) = newchild(i)
               which_child(k) = which_child(i)
               lrefine(k) = lrefine(i)
               bflags(:,k) = bflags(:,i)
               work_block(k) = work_block(i)
               If (empty_cells) Then
                  empty(k) = empty(i)
               End If
               
            End If

            k = k + 1
            
         Else
            
            lnblocks2 = lnblocks2 - 1
            lnblocks_old = lnblocks_old - 1
            
         End If ! End If (.not.derefine(i))
         
      End Do ! End Do i = 1,lnblocks

!-----overwrite old locations
      Do i = lnblocks2+1,lnblocks
         
         derefine(i) = .FALSE.
         Do j = 1,nchild
            child(1,j,i) = -1
            child(2,j,i) = -1
         End Do
         parent(1,i) = -1
         parent(2,i) = -1
         Do j = 1,nfaces
            neigh(1,j,i) = -1
            neigh(2,j,i) = -1
         End Do
         Do j = 1,ndim
            coord(j,i) = -1.
            bnd_box(1,j,i) = -1.
            bnd_box(2,j,i) = -1.
         End Do
         bsize(:,i) = -1.
         nodetype(i) = -1
         which_child(i) = -1
         newchild(i) = .FALSE.
         lrefine(i) = -1
         bflags(:,i) = -1
         work_block(i) = 0.
         If (empty_cells) Then
            empty(i) = 0
         End If
         
      End Do ! End Do i = lnblocks2+1,lnblocks
      
      lnblocks = lnblocks2

!-----reset node types
      Do i = 1,lnblocks
         nodetype(i) = 3
         If (child(1,1,i) <= -1) Then
            nodetype(i) = 1
         End If
      End Do
      nrecv = 0
      Do i = 1,lnblocks
         Do j = 1,nchild
            nodetype_chi(j,i) = -1
            If (child(1,j,i) > -1) Then
            If (child(2,j,i).ne.mype) Then
               nrecv = nrecv + 1
               Call MPI_IRECV(nodetype_chi(j,i), &
                              1, &
                              MPI_INTEGER, &
                              child(2,j,i), &
                              child(1,j,i), &
                              MPI_COMM_WORLD, &
                              reqr(nrecv), &
                              ierr)
            Else
               nodetype_chi(j,i) = nodetype(child(1,j,i))
            End If
            End If
         End Do
      End Do ! End Do i = 1, lnblocks

      nsend = 0
      Do i = 1,lnblocks !--------send nodetype to your parent
         If (parent(1,i) >= 1) Then
         If (parent(2,i).ne.mype) Then
            nsend = nsend + 1
            Call MPI_SSEND(nodetype(i), &
                           1, &
                           MPI_INTEGER, &
                           parent(2,i), &
                           i, &
                           MPI_COMM_WORLD, &
                           ierr)
         End If
         End If
      End Do ! End Do i = 1, lnblocks

      If (nrecv > 0) Then
         Call MPI_WAITALL (nrecv, reqr, statr, ierr)
      End If

      Do i = 1,lnblocks
         Do j = 1,nchild
            If (nodetype_chi(j,i) == 1) nodetype(i) = 2
         End Do
      End Do

!-----reset derefine flags
      Do i = 1,maxblocks_tr
         derefine(i) = .FALSE.
      End Do

      deallocate(new_loc)
      deallocate(neight)
      deallocate(childt)
      deallocate(parentt)
      deallocate(statr)
      deallocate(reqr)
      deallocate(nodetype_chi)

      Return
      End Subroutine amr_derefine_blocks
