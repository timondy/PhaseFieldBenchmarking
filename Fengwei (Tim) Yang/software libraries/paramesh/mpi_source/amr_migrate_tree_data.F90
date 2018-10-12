!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/amr_migrate_tree_data
!! NAME
!!
!!   amr_migrate_tree_data
!!
!! SYNOPSIS
!!
!!   call amr_migrate_tree_data (new_loc,nprocs,mype)
!!   call amr_migrate_tree_data (integer array, integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(inout) :: new_loc(:,:)
!!     new locations (processor and location in morton list) to migrate data to
!!   integer, intent(in) :: nprocs, mype
!!     number of procs. and proc. id.
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
!!   io
!!   paramesh_comm_data
!!
!! CALLS
!!
!!   fill_old_loc
!! 
!! RETURNS
!!
!!   Nothing returned
!!
!! DESCRIPTION
!!
!!   This subroutine moves tree data and reconnects all pointers according to the 
!!   new locations which are passed in in the argument 'new_loc'.  
!!   'new_loc' is constructed during the morton sort and ordering by work 
!!   before it is passed to this routine.
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2001).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_migrate_tree_data (new_loc,nprocs,mype)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_comm_data
      Use paramesh_interfaces, only : fill_old_loc

      Implicit None

      Include 'mpif.h'

!-----Input/Output variables.
      Integer, Intent(inout) :: new_loc(:,:)
      Integer, Intent(in)    :: nprocs,mype

!-----Local variables and arrays
      Integer, Parameter :: buf_size = mdim+mdim+2*mdim
      Integer, Parameter :: ibuf_size = 2*mfaces+2*mchild+2+5+10
      Real    :: buffer(buf_size)
!      Real    :: buffert(buf_size,maxblocks)
      Integer :: ibuffer(ibuf_size)
!      Integer :: ibuffert(ibuf_size,maxblocks)
!      Integer :: neight(2,mfaces,maxblocks_tr)
!      Integer :: childt(2,mchild,maxblocks_tr)
!      Integer :: parentt(2,maxblocks_tr)
      Integer :: i,j,k,jj
!      Integer :: old_loc(2,maxblocks_tr)
!      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
!      Integer :: reqr(maxblocks_tr)
      Integer :: ierr,nsend,nrecv
!      Logical :: newchildt(maxblocks_tr)
! CEG rewritten to have all the arrays as allocated rather than on the stack
      Integer,allocatable :: neight(:,:,:)
      Integer,allocatable :: childt(:,:,:)
      Integer,allocatable :: parentt(:,:)
      Integer,allocatable :: old_loc(:,:)
      Integer,allocatable :: statr(:,:)
      Integer,allocatable :: reqr(:)
      Logical,allocatable :: newchildt(:)
      Real,allocatable :: buffert(:,:)
      Integer,allocatable :: ibuffert(:,:)

! CEG allocates memory
      allocate(old_loc(2, maxblocks_tr))
      allocate(neight(2, mfaces, maxblocks_tr))
      allocate(childt(2, mchild, maxblocks_tr))
      allocate(parentt(2, maxblocks_tr))
      allocate(statr(MPI_STATUS_SIZE, maxblocks_tr))
      allocate(reqr(maxblocks_tr))
      allocate(newchildt(maxblocks_tr))
      allocate(buffert(buf_size,maxblocks))
      allocate(ibuffert(ibuf_size,maxblocks))


!-----Begin executable code.

      Call fill_old_loc(new_loc,old_loc,nprocs,mype)

!-----count no. of new blocks
      new_lnblocks = 0
      Do i = 1,maxblocks
        If (old_loc(1,i) > -1) Then
          new_lnblocks = new_lnblocks + 1
        End If
      End Do

!-----update pointers to parents, children and neighbors
      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
             nrecv = nrecv + 1
             Call MPI_IRECV(parentt(1,i),2,MPI_INTEGER,                & 
                  parent(2,i),i,MPI_COMM_WORLD,                        & 
                  reqr(nrecv),ierr)
           Else
             parentt(:,i) = new_loc(:,parent(1,i))
           End If
         End If
       End Do
       
       nsend = 0
       Do i = 1,lnblocks
         Do j = 1,nchild
           If (child(1,j,i) > 0) Then
             If (child(2,j,i).ne.mype) Then
               ! parent is sending to all its children
               nsend = nsend + 1
               Call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER,             & 
                    child(2,j,i),child(1,j,i),MPI_COMM_WORLD,          & 
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
              Call MPI_IRECV(childt(1,j,i),2,MPI_INTEGER,              & 
                   child(2,j,i),child(1,j,i),MPI_COMM_WORLD,           & 
                   reqr(nrecv),ierr)
            Else
              childt(:,j,i) = new_loc(:,child(1,j,i))
            End If
          End If
        End Do
       End Do
       
       nsend = 0
       Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
!------------child is sending to its parent
             nsend = nsend + 1
             Call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER,               & 
                  parent(2,i),i,MPI_COMM_WORLD,                        & 
                  ierr)
           End If
         End If
       End Do

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      nrecv = 0
      Do i = 1,lnblocks
        Do j = 1,nfaces
          If (neigh(1,j,i) > 0) Then
            If (neigh(2,j,i).ne.mype) Then
              nrecv = nrecv + 1
              Call MPI_IRECV(neight(1,j,i),2,MPI_INTEGER,              & 
                   neigh(2,j,i),neigh(1,j,i),MPI_COMM_WORLD,           & 
                   reqr(nrecv),ierr)
           Else
             neight(:,j,i) = new_loc(:,neigh(1,j,i))
            End If
          End If
        End Do
      End Do

      nsend = 0
      Do i = 1,lnblocks
        Do j = 1,nfaces
          If (neigh(1,j,i) > 0) Then
            If (neigh(2,j,i).ne.mype) Then
              nsend = nsend + 1
              Call MPI_int_SSEND (new_loc(1,i),2,MPI_INTEGER,          & 
                   neigh(2,j,i),i,MPI_COMM_WORLD,                      & 
                   ierr)
            End If
          End If
        End Do
      End Do

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      parent(:,1:lnblocks) = parentt(:,1:lnblocks)
      child(:,:,1:lnblocks) = childt(:,:,1:lnblocks)
      neigh(:,:,1:lnblocks) = neight(:,:,1:lnblocks)

!-----initialize temp buffer array
      Do i = 1,maxblocks
         buffert(:,i) = -1.
         ibuffert(:,i) = -1
         newchildt(i) = .FALSE.
      End Do

      nrecv = 0
      Do i = 1,new_lnblocks
        If (old_loc(2,i).ne.mype) Then
          nrecv = nrecv + 1
          Call MPI_IRECV(buffert(1,i),buf_size,                        & 
               amr_mpi_real,                                           & 
               old_loc(2,i),i,MPI_COMM_WORLD,                          & 
               reqr(nrecv),ierr)
          nrecv = nrecv + 1
          Call MPI_IRECV(ibuffert(1,i),ibuf_size,MPI_INTEGER,          & 
               old_loc(2,i),i+maxblocks_tr,MPI_COMM_WORLD,             & 
               reqr(nrecv),ierr)
          nrecv = nrecv + 1
          Call MPI_IRECV(newchildt(i),1,MPI_LOGICAL,                   & 
               old_loc(2,i),i+2*maxblocks_tr,MPI_COMM_WORLD,           & 
               reqr(nrecv),ierr)
        End If
      End Do

      nsend = 0

      Do i = 1,lnblocks

!-------pack buffer for sending
        k = 0
        Do j = 1,mdim
          k = k + 1
          buffer(k) = coord(j,i)
        End Do
        Do j = 1,mdim
          Do jj = 1,2
            k = k + 1
            buffer(k) = bnd_box(jj,j,i)
          End Do
        End Do
        Do j = 1,mdim
          k = k + 1
          buffer(k) = bsize(j,i)
        End Do

        k = 0
        Do j = 1,mchild
          Do jj = 1,2
            k = k + 1
            ibuffer(k) = child(jj,j,i)
          End Do
        End Do
        Do j = 1,mfaces
          Do jj = 1,2
            k = k + 1
            ibuffer(k) = neigh(jj,j,i)
          End Do
        End Do
        Do j = 1,2
          k = k + 1
          ibuffer(k) = parent(j,i)
        End Do
        k = k + 1
        ibuffer(k) = lrefine(i)
        k = k + 1
        ibuffer(k) = nodetype(i)
        k = k + 1
        ibuffer(k) = which_child(i)
        k = k + 1
        ibuffer(k) = empty(i)
        k = k + 1
        ibuffer(k) = work_block(i)
        Do j=1,mflags
        k = k + 1
        ibuffer(k) = bflags(j,i)
        End Do

        If (new_loc(2,i).ne.mype) Then
          nsend = nsend + 1
          Call MPI_SSEND(buffer(1),buf_size,amr_mpi_real,              & 
               new_loc(2,i),new_loc(1,i),                              & 
               MPI_COMM_WORLD,ierr)
          nsend = nsend + 1
          Call MPI_SSEND(ibuffer(1),ibuf_size,MPI_INTEGER,             & 
               new_loc(2,i),new_loc(1,i)+maxblocks_tr,                 & 
               MPI_COMM_WORLD,ierr)
          nsend = nsend + 1
          Call MPI_SSEND(newchild(i),1,MPI_LOGICAL,                    & 
               new_loc(2,i),new_loc(1,i)+2*maxblocks_tr,               & 
               MPI_COMM_WORLD,ierr)
        Else
          buffert(1:buf_size,new_loc(1,i)) = buffer(1:buf_size)
          ibuffert(1:ibuf_size,new_loc(1,i)) = ibuffer(1:ibuf_size)
          newchildt(new_loc(1,i)) = newchild(i)
        End If

      End Do  ! End Do i = 1,lnblocks

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

!-----now unpack the buffer
      Do i = 1,new_lnblocks

        k = 0
        Do j = 1,mdim
          k = k + 1
          coord(j,i) = buffert(k,i)
        End Do
        Do j = 1,mdim
          Do jj = 1,2
            k = k + 1
            bnd_box(jj,j,i) = buffert(k,i)
          End Do
        End Do
        Do j = 1,mdim
          k = k + 1
          bsize(j,i) = buffert(k,i)
        End Do
        k = 0
        Do j = 1,mchild
          Do jj = 1,2
            k = k + 1
            child(jj,j,i) = ibuffert(k,i)
          End Do
        End Do
        Do j = 1,mfaces
          Do jj = 1,2
            k = k + 1
            neigh(jj,j,i) = ibuffert(k,i)
          End Do
        End Do
        Do j = 1,2
          k = k + 1
          parent(j,i) = ibuffert(k,i)
        End Do
        k = k + 1
        lrefine(i) = ibuffert(k,i)
        k = k + 1
        nodetype(i) = ibuffert(k,i)
        k = k + 1
        which_child(i) = ibuffert(k,i)
        k = k + 1
        empty(i) = max(ibuffert(k,i),0)     ! empty must be either 1 or 0.
                                            ! It cannot be -1.
        k = k + 1
        work_block(i) = ibuffert(k,i)
        Do j=1,mflags
        k = k + 1
        bflags(j,i) = ibuffert(k,i)
        End Do

        newchild(i) = newchildt(i)

      End Do  ! End Do i = 1,new_lnblocks


! CEG deallocate memory
      deallocate(old_loc)
      deallocate(neight)
      deallocate(childt)
      deallocate(parentt)
      deallocate(statr)
      deallocate(reqr)
      deallocate(newchildt)
      deallocate(buffert)
      deallocate(ibuffert)


      Return
      End Subroutine amr_migrate_tree_data


