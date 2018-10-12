!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine process_fetch_list(fetch_list,                        &
                                    istack,                            &
                                    mype,                              &
                                    nprocs,                            &
                                    n_to_left,                         &
                                    tag_offset)


      Use tree
      Use paramesh_Dimensions
      Use mpi_morton
      Use paramesh_mpi_interfaces, only : compress_fetch_list 

      Implicit None

      Include 'mpif.h'

      Integer, Intent(inout), Dimension(:,:) :: fetch_list
      Integer, Intent(in)    :: istack, mype, nprocs
      Integer, Intent(in)    :: n_to_left(0:nprocs-1)
      Integer, Intent(inout) :: tag_offset

      Integer :: no_of_remote_neighs, max_no_to_be_received
      Integer,Dimension (:),  allocatable :: recvrequest
      Integer,Dimension (:,:),allocatable :: recvstatus
      Integer :: i, j, k, ierror, ierrorcode, iprocs, ii, jj, jm1, jp
      Integer :: ll, kk, isrc, idest, itag, isize
      integer :: ierr
      Integer, allocatable :: reqR(:), reqS(:)!, stat(:,:)

!-----Compress the list of possible off processor blocks by eliminating 
!-----redundant entries (this routine also sorts the list)

      no_of_remote_neighs = 0
      If (istack > 0) Then
       Call compress_fetch_list(fetch_list,                            &
                                istack,                                &
                                no_of_remote_neighs,                   &
                                mype,                                  &
                                nprocs,                                &
                                n_to_left)
      End If

      Do i = 1, no_of_remote_neighs
        fetch_list(2,i) = fetch_list(2,i) + 1
      End Do

!-----Construct commatrix_recv (the number of blocks to receive from each
!-----processor)

      commatrix_send(:) = 0
      commatrix_recv(:) = 0
      do i = 1,no_of_remote_neighs
         commatrix_recv(fetch_list(2,i)) =                             &
            commatrix_recv(fetch_list(2,i)) + 1
      End Do
      
!-----Constrcut commatrix_send (the number of blocks to send to each other processor) 
!-----by providing the complete commatrix_recv to all processors

! CEG removed to just transpose
!      Call MPI_ALLTOALL (commatrix_recv,1,MPI_INTEGER,                 & 
!                         commatrix_send,1,MPI_INTEGER,                 & 
!                         MPI_COMM_WORLD,ierror)

      allocate(reqS(nprocs-1))
      allocate(reqR(nprocs-1))
!      allocate(stat(MPI_STATUS_SIZE,nprocs))
!      allocate(stat(MPI_STATUS_SIZE,nprocs))
      reqS(:) = 0
      reqR(:) = 0


! CEG Swapped the order of the two MPI calls over and combined them into one loop, along with iprocs calculation
      iprocs = 0
      do i = 1, nprocs
            iprocs = iprocs + min(1,commatrix_recv(i))
      end do

      ii = 0
      do i = 1, nprocs
         if (i-1.eq.mype) then
           commatrix_send(i) = commatrix_recv(i)
         else
            ii = ii + 1
!         call MPI_Recv(commatrix_send(i), 1, MPI_INTEGER, i-1, i-1, MPI_COMM_WORLD, stat(1,i), ierr)
            call MPI_IRecv(commatrix_send(i), 1, MPI_INTEGER, i-1, i-1, MPI_COMM_WORLD, reqR(ii), ierr)
            call MPI_ISend(commatrix_recv(i), 1, MPI_INTEGER, i-1, mype, MPI_COMM_WORLD, reqS(ii), ierr)
!            print*,'reqS(i)',mype,i,reqS(ii), reqR(ii)
         endif
      end do

!-----Compute the maximum no. of blocks which any processor
!-----is going to receive.
       max_no_to_be_received = max(1,iprocs)

!-----Evaluate smallest guard block starting index over all pe
!-----store this into variable strt_buffer which is Used in amr_1blk_guardcell

      last_buffer = maxblocks_alloc

! CEG Moved from higher up to allow maximum (albeit small) Comp/comm overlap
      iprocs = 0
      if (nprocs.gt.1) then
         call MPI_Waitall(nprocs-1, reqR, MPI_STATUSES_IGNORE, ierr)

!CEG Finish overlapping communication with computation
         call MPI_Waitall(nprocs-1, reqS, MPI_STATUSES_IGNORE, ierr)
!         call MPI_Waitall(nprocs-1, reqS, stat, ierr)
      endif
      deallocate(reqR)
      deallocate(reqS)

      k = last_buffer
      Do i = 0,nprocs-1
         k = k - commatrix_recv(i+1)
      End Do
      strt_buffer = k 

      If (strt_buffer <= lnblocks) Then
        Write(*,*)  & 
        'ERROR in process_fetch_list : guard block starting index',    & 
        strt_buffer,' not larger than lnblocks',lnblocks,              & 
        ' processor no. ',mype,' maxblocks_alloc ',                    & 
        maxblocks_alloc
        Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierror)
      End If

      If (Allocated(to_be_sent))     Deallocate(to_be_sent)
      If (Allocated(to_be_received)) Deallocate(to_be_received)

!------Compute the maximum no. of blocks which any processor
!------is going to receive.

       do i = 1,nprocs
          iprocs = iprocs + min(1,commatrix_send(i))
       End Do
       max_no_to_send = max(1,iprocs)

!-----Dynamically allocate memory to store the lists of blocks to be
!-----sent and received.

      Allocate ( to_be_sent(3,                                         & 
                            max(1,maxval(commatrix_send)),             & 
                            max(1,max_no_to_send) ) )
      Allocate ( to_be_received(3,                                     & 
                                max(1,maxval(commatrix_recv)),         & 
                                max(1,max_no_to_be_received) ) )

!-----Construct arrays to_be_sent and to_be_received which contain
!-----the lists of blocks to be packaged.

      to_be_sent = -1
      to_be_received = -1
      laddress = 0

!-----Set up the array to_be_received on each processor
      If (no_of_remote_neighs > 0) Then

         jp = 1
         ii = 0
         jm1 = 1
         Do jj = 1, no_of_remote_neighs
            If (jj == 1 .or.                                          &
            fetch_list(2,jj) == fetch_list(2,jm1)) Then
               ii = ii + 1
            Else
               jp = jp + 1
               ii = 1
            End If
            to_be_received(:,ii,jp) = fetch_list(:,jj)
            jm1 = jj
         End Do

         jj = no_of_remote_neighs
         laddress(1,strt_buffer:strt_buffer+jj-1) = fetch_list(1,1:jj)
         laddress(2,strt_buffer:strt_buffer+jj-1) = fetch_list(2,1:jj)-1
         
      End If

      If (Allocated(recvrequest)) Deallocate( recvrequest )
      Allocate ( recvrequest(nprocs) )
      If (Allocated(recvstatus)) Deallocate( recvstatus )
      Allocate ( recvstatus(MPI_STATUS_SIZE, nprocs) )

! Post receives 
      kk = 0
      Do i = 1,nprocs
         
         isrc = i-1
         idest= mype
         itag = isrc*nprocs + idest+1 + tag_offset

                                ! receive to pe=j
         If (commatrix_send(i).gt.0) Then
            kk = kk+1
            isize = 3*commatrix_send(i)
            call MPI_IRECV(to_be_sent(1,1,kk),isize,                   & 
                 MPI_INTEGER,isrc ,itag,MPI_COMM_WORLD,                & 
                 recvrequest(kk),ierror)
         End If
      End Do

! Post sends

      ll = 0
      Do j = 1,nprocs
         
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
         
                                ! send from mype=i
         If (commatrix_recv(j).gt.0) Then
            ll = ll+1
            isize = 3*commatrix_recv(j)
            Call MPI_SSEND(to_be_received(1,1,ll),isize,MPI_INTEGER,   & 
                           idest,itag,MPI_COMM_WORLD,ierror)
         End If
      End Do

      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset

      If (kk.gt.0) Then
         Call MPI_Waitall(kk,recvrequest,recvstatus,ierror)
      End If
     
      If (Allocated(recvrequest)) Deallocate( recvrequest )
      If (Allocated(recvstatus)) Deallocate( recvstatus )


      End Subroutine process_fetch_list

