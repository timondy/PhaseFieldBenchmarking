!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/amr_sort_by_work
!! NAME
!!
!!   amr_sort_by_work
!!
!! SYNOPSIS
!!
!!   call amr_sort_by_work (new_loc, nprocs, mype)
!!
!!   call amr_compute_morton(integer, integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(inout) :: new_loc(:,:)
!!    The new locations which blocks are to migrate to.  new_loc(1,lb) indicates the
!!    location internal to processor and new_loc(2,lb) indicates which processor
!!    to migrate block 'lb' to in order to acheive load balancing.
!!
!!   integer, intent(in)    :: nprocs
!!    The number of processors.
!!
!!   integer, intent(in)    :: mype
!!     The processor id of the calling processor.
!!
!! INCLUDES
!!
!!   mpif.h
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   io
!!   paramesh_interfaces
!!   paramesh_comm_data
!!
!! CALLS
!! 
!!   fill_old_loc
!!    
!! RETURNS
!!
!!   Returns, in a list (new_loc) the locations where the ordered blocks are to
!!   migrate to.
!!
!! DESCRIPTION
!!
!!   Subroutine to compute where blocks should migrate to in order to acheive 
!!   work load balancing.  The morton ordered list is scanned by this routine
!!   and the total amount of 'work' (as indicated by the user supplied work weighting
!!   values for each block stored in the array 'work_block').  This information
!!   is Then used to compute where to 'cut' the morton ordered list into a 
!!   number of pieces equal to the number of processors.  Each piece has roughly the
!!   the same total amout of work.  This routine does not move the blocks from
!!   their current locations in the list, but only computes where they are to
!!   move to and returns this information in the array 'new_loc'.
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2004).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_sort_by_work (new_loc,nprocs,mype)

!-----Use statements.
      use paramesh_dimensions
      Use physicaldata
      Use tree
      Use paramesh_interfaces, only : fill_old_loc
      Use io
      Use paramesh_comm_data

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      integer, intent(inout) :: new_loc(:,:)
      integer, intent(in)    :: nprocs,mype

!-----Local variables and arrays.
      Real    :: work_per_proc,work_left,loc_work,tot_work
!      Real    :: work(maxblocks_tr),workt(maxblocks_tr)
      Real,allocatable    :: work(:),workt(:)
      Real    :: reduce_datain(1),reduce_dataout(1)
      Real    :: rdatain(1),rdataout(1)
      Integer :: ierr,errorcode
      Integer :: lnblocks2
      Integer :: i,j
      Integer :: pidt=-1,lidt,lid_old
      Integer :: left,right
      Integer :: nsend,nrecv
!      Integer :: reqr(2*maxblocks_tr)
!      Integer :: pid(maxblocks_tr),lid(maxblocks_tr),lid2(maxblocks_tr)
!      Integer :: stat (MPI_STATUS_SIZE)
!      Integer :: new_loc_temp(2,maxblocks_tr)
!      Integer :: old_loc(2,maxblocks_tr)
!      Integer :: statr(MPI_STATUS_SIZE,2*maxblocks_tr)
      Integer,allocatable :: reqr(:)
      Integer,allocatable :: pid(:),lid(:),lid2(:)
      Integer,allocatable :: stat (:)
      Integer,allocatable :: new_loc_temp(:,:)
      Integer,allocatable :: old_loc(:,:)
      Integer,allocatable :: statr(:,:)
      Logical :: repeat,repeatt
      Logical :: lreduce_datain(1),lreduce_dataout(1)

!-----Begin executable code.
! CEG alloaction
  allocate(work(maxblocks_tr))
  allocate(workt(maxblocks_tr))
  allocate(reqr(2*maxblocks_tr))
  allocate(pid(maxblocks_tr))
  allocate(lid(maxblocks_tr))
  allocate(lid2(maxblocks_tr))
  allocate(stat (MPI_STATUS_SIZE))
  allocate(new_loc_temp(2,maxblocks_tr))
  allocate(old_loc(2,maxblocks_tr))
  allocate(statr(MPI_STATUS_SIZE,2*maxblocks_tr))


!-----Initialize work arrary and temp work array.
      Do i = 1,maxblocks_tr
         work(i) = -1.
         workt(i) = -1.
      End Do 

!-----Assign values to work array.
      Do i = 1,maxblocks_tr
         work(i) = 0.
         work(i) = work_block(i)
      End Do 

!-----move work to temp work array workt
      Call fill_old_loc(new_loc,old_loc,nprocs,mype)

      lnblocks2 = 0
      Do i = 1,maxblocks_tr
        If (old_loc(1,i) > -1) Then
          lnblocks2 = lnblocks2 + 1
        End If
      End Do 

      nrecv = 0
      Do i = 1,lnblocks2
        If (old_loc(2,i).ne.mype) Then
          nrecv = nrecv + 1
          Call MPI_IRECV(workt(i),1,amr_mpi_real,                      & 
               old_loc(2,i),old_loc(1,i),MPI_COMM_WORLD,               & 
               reqr(nrecv),ierr)
        End If
      End Do 

      nsend = 0
      Do i = 1,lnblocks
        If (new_loc(2,i).ne.mype) Then
          nsend = nsend + 1
          Call MPI_SSEND(work(i),1,amr_mpi_real,                       & 
               new_loc(2,i),i,MPI_COMM_WORLD,ierr)
        Else
          workt(new_loc(1,i)) = work(i)
        End If
      End Do 
         
      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      Do i = 1,lnblocks2
         work(i) = workt(i)
      End Do 

!-----SUM total work within each processosr
      If (lnblocks2 > 0) Then
        workt(1) = work(1)
      Else
        workt(1) = 0
      End If
      Do i = 2,lnblocks2
         workt(i) = workt(i-1) + work(i)
      End Do 

!-----SUM work across processors
      If (lnblocks2 > 0) Then
        loc_work = workt(lnblocks2)
      Else
        loc_work = 0
      End If
      reduce_datain(1) = loc_work
!  print *, 'CEG - Allreduce in amr_sort_by_work'
      Call MPI_ALLREDUCE (reduce_datain(1),reduce_dataout(1),          & 
           1,amr_mpi_real,                                             & 
           MPI_SUM,MPI_COMM_WORLD,ierr)
      tot_work = reduce_dataout(1)

!-----Compute work per processor
      work_per_proc = tot_work/nprocs

!-----Compute final work by looking left
      work_left = 0.

      rdatain(1) = loc_work
      Call MPI_SCAN (rdatain(1),rdataout(1),1,                         & 
                     amr_mpi_real,MPI_SUM,MPI_COMM_WORLD,              & 
                     ierr)
      work_left = rdataout(1)
      work_left = work_left - loc_work
      Do j = 1,lnblocks2
         workt(j) = workt(j) + work_left
      End Do 

!-----compute processor ids
      Do i = 1,maxblocks_tr
         pid(i) = 0
         lid(i) = 0
      End Do 

      Do i = 1,lnblocks2
         pid(i) = int((workt(i)-1.)/work_per_proc)
         If (pid(i) < 0) pid(i) = 0
         If (pid(i) > nprocs-1) pid(i) = nprocs-1
      End Do 

!-----compute local ids
      lid(1) = 1
      Do i = 2,lnblocks2
         lid(i) = lid(i-1) + 1
         If (pid(i-1) < pid(i)) lid(i) = 1  ! start a new group
      End Do 

      Do i = 1,maxblocks_tr
         lid2(i) = lid(i)
      End Do 

      left = mype - 1
      right = mype + 1
      If (mype == 0) left = MPI_PROC_NULL
      If (mype == nprocs-1) right = MPI_PROC_NULL

      Call MPI_SENDRECV                                                & 
           (pid(lnblocks2),1,MPI_INTEGER,right,1,                      & 
            pidt,          1,MPI_INTEGER,left, 1,                      & 
            MPI_COMM_WORLD,stat,ierr)

      lidt = 0
 27   lid_old = lidt ! lid_old stores last fetched value of lid to left

      lidt = 0
      Call MPI_SENDRECV                                                & 
           (lid(lnblocks2),1,MPI_INTEGER,right,1,                      & 
            lidt,          1,MPI_INTEGER,left, 1,                      & 
            MPI_COMM_WORLD,stat,ierr)

      Do j = 1,lnblocks2
         If (pidt == pid(j)) Then ! if pidt (which was fetched)
                                  ! equals local pid Then the list
                                  ! has been split across processors
               
            lid(j) = lid2(j) + lidt
            
         End If
      End Do 
      
      repeat = .FALSE.
      If (lidt.ne.lid_old) repeat = .TRUE.
      lreduce_datain(1) = repeat
!  print *, 'CEG - Allreduce in amr_sort_by_work'
      Call MPI_ALLREDUCE(lreduce_datain(1),lreduce_dataout(1),         & 
                         1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
      repeatt = lreduce_dataout(1)
      If (repeatt) go to 27

!-----now reorder according to new pid and lid numbers
      nrecv = 0
      Do i = 1,lnblocks
         If (new_loc(2,i).ne.mype) Then
            nrecv = nrecv + 1
            Call MPI_IRECV(new_loc_temp(1,i),1,MPI_INTEGER,            & 
                 new_loc(2,i),new_loc(1,i),MPI_COMM_WORLD,             & 
                 reqr(nrecv),ierr)
         Else
            new_loc_temp(1,i) = lid(new_loc(1,i))
         End If
      End Do 

      nsend = 0
      Do i = 1,lnblocks2
        If (old_loc(2,i).ne.mype) Then
           nsend = nsend + 1
           Call MPI_SSEND(lid(i),1,MPI_INTEGER,                        & 
                old_loc(2,i),i,                                        & 
                MPI_COMM_WORLD,ierr)
        End If
      End Do 

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      nrecv = 0
      Do i = 1,lnblocks
         If (new_loc(2,i).ne.mype) Then
            nrecv = nrecv + 1
            Call MPI_IRECV(new_loc_temp(2,i),1,MPI_INTEGER,            &  
                 new_loc(2,i),new_loc(1,i),                            & 
                 MPI_COMM_WORLD,                                       & 
                 reqr(nrecv),ierr)
         Else
            new_loc_temp(2,i) = pid(new_loc(1,i))
         End If
      End Do 
      
      nsend = 0
      Do i = 1,lnblocks2
         If (old_loc(2,i).ne.mype) Then
           nsend = nsend + 1
           Call MPI_SSEND(pid(i),1,MPI_INTEGER,                        & 
                old_loc(2,i),i,                                        & 
                MPI_COMM_WORLD,ierr)
        End If
      End Do 
         
      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      Do i = 1,lnblocks

         new_loc(:,i) = new_loc_temp(:,i)

         If (new_loc(1,i) > maxblocks) Then
            Open (unit=30, & 
                 file=amr_log_file,                                    & 
                 position='append',                                    & 
                 status='unknown',                                     & 
                 form='formatted')
            Write(30,*) 'PARAMESH ERROR !'
            Write(30,*) mype,i
            Write(30,*) ' new_loc(1 = ',new_loc(1,i),new_loc(2,i)
            Write(30,*) 'New block location exceeds MAXBLOCKS limit'
            Write(30,*) 'Suggestion: increase MAXBLOCKS or modify',    & 
                 ' refinement criteria'
            Close(30)
            Call MPI_ABORT (MPI_COMM_WORLD,errorcode,ierr)
         End If         
         
         If (new_loc(2,i) > nprocs.or.new_loc(2,i) < 0) Then
            Open (unit=30,                                             & 
                 file=amr_log_file,                                    & 
                 position='append',                                    & 
                 status='unknown',                                     & 
                 form='formatted')
            Write(30,*) 'PARAMESH ERROR !'
            Write(30,*) mype,i
            Write(30,*) 'New block location out of bounds',new_loc(2,i)
            Write(30,*) 'Suggestion: increase MAXBLOCKS or modify', & 
                        ' refinement criteria'
            Close(30)
            Call MPI_ABORT (MPI_COMM_WORLD,errorcode,ierr)
         End If         

      End Do 
      new_loc(:,lnblocks+1:maxblocks_tr) = -1

! CEG free
  deallocate(reqr)
  deallocate(pid)
  deallocate(lid)
  deallocate(lid2)
  deallocate(stat)
  deallocate(new_loc_temp)
  deallocate(old_loc)
  deallocate(statr)
  deallocate(work)
  deallocate(workt)


      Return
      End Subroutine amr_sort_by_work

