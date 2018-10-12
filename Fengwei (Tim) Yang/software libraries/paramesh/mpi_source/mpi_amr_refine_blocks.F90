!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_refine_blocks
!! NAME
!!
!!   amr_refine_blocks
!! 
!! SYNOPSIS
!!
!!   Call amr_refine_blocks (nprocs, mype)
!!
!!   Call amr_refine_blocks (integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: nprocs
!!     The number of processors being used.
!!
!!   integer, intent(in) :: mype     
!!     The Calling processor number.
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
!!   paramesh_interfaces
!!   mpi_morton
!!   timings
!!   constants
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   morton_neighbors
!!
!! RETURNS
!!
!!   This routine does not return anything.  Upon return new child blocks have been
!!   created and their locations in the tree data structure have been established.
!!
!! DESCRIPTION
!!
!!   The routine takes care of adding blocks to local list of blocks based 
!!   the setting of the flags set in the logical array 'refine(:)'.  If 
!!   refine(lb) is set to .true. for block lb, then the children of that block
!!   are created and added to the end of the list of on-processor blocks.
!!   The tree data is also updated by this routine to reflect these new blocks.
!!   The routine does NOT redistribute any data between processors.  That is handled
!!   by the routine amr_redist_blk.
!!
!!   Note that this routine is Called by the routine 'amr_refine_derefine' and should
!!   not need to be Called directly by a user's application.
!!
!! AUTHORS
!!
!!   Kevin Olson (original version 1996, MPI version 1997,1998)
!!   Modifications by Mike Zingale and Jonathan Dursi (U. of Chicago), 1999
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_refine_blocks (nprocs,mype)

      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use timings
      Use constants
      Use paramesh_mpi_interfaces, only :  morton_neighbors

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output statements.
      Integer, Intent(in)    :: nprocs,mype

!-----Local array and variables.
      Real :: hx,hy,hz
      Real :: hh(mdim)      
      Real :: eps,accuracy
      Integer :: i,j,k,jj,ipar,ipar_proc,iproc
      Integer :: ineigh,ineigh_proc
      Integer :: lnblocks2,ichi
      Integer :: ierr,jr,n3
      Integer :: nsend,nrecv
      Integer :: istart,iend,istart2,iend2,istep
      Integer :: errorcode
      Integer,allocatable :: nsend_to_proc(:), nrecv_pack(:)
      Integer,allocatable :: reqr(:)
      Integer,allocatable :: kk(:)
      Integer,allocatable :: nodetype_chi(:,:)
      Integer,allocatable :: statr(:,:)
      Integer,allocatable :: tneigh(:,:,:)
!      Integer :: nsend_to_proc(0:16384), nrecv_pack(0:16384)
!      Integer :: reqr(maxblocks_tr)
!      Integer :: kk(maxblocks_tr)
!      Integer :: nodetype_chi(nchild,maxblocks_tr)
!      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
!      Integer :: tneigh(2,mfaces,maxblocks_tr)
      Integer :: idummy_array(1)
      Integer :: iii,ii,jjj,kkk
      Integer :: jr0,kr,isweep,nsweep
      Logical :: lrecv
      Logical :: lsend
!      Logical :: refine_neigh(maxblocks_tr)
      Logical,allocatable :: refine_neigh(:)

!-----Begin executable code.
! CEG allocate memory
  allocate(nsend_to_proc(0:16384))
  allocate(nrecv_pack(0:16384))
  allocate(reqr(maxblocks_tr))
  allocate(kk(maxblocks_tr))
  allocate(nodetype_chi(nchild,maxblocks_tr))
  allocate(statr(MPI_STATUS_SIZE,maxblocks_tr))
  allocate(tneigh(2,mfaces,maxblocks_tr))
  allocate(refine_neigh(maxblocks_tr))

      accuracy = 10./10.**precision(accuracy)
      eps = accuracy

!-----refine blocks, create their children, the children's parents,
!-----turn the children on, and turn the parents off
      Do i = 1,maxblocks_tr
         newchild(i) = .False.
      End Do

      lnblocks2 = lnblocks

      Do i = 1,lnblocks

         If (refine(i)) Then ! refine block 'i' on this processor

            Do j = 1,nchild ! create 'nchild' new child blocks for block 'i'

               lnblocks2 = lnblocks2 + 1

               if(lnblocks2 > maxblocks_tr) Then
                  Open(unit=30,file='amr_log',status='unknown',        & 
                       position='append')
                  Write(30,*) 'PARAMESH ERROR !'
                  Write(30,*) 'Too many blocks created! '
                  Write(30,*) 'Increase MAXBLOCKS_TR and rerun! '
                  Close(30)
                  Write(*,*) 'PARAMESH ERROR !'
                  Write(*,*) 'Too many blocks created! '
                  Write(*,*) 'Increase MAXBLOCKS_TR and rerun! '
                  Call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
               End If

               child(1,j,i) = lnblocks2 ! child j's on-processor id
               child(2,j,i) = mype      ! child j's processor no.

               lrefine(lnblocks2) = lrefine(i) + 1 ! refinement level of child
               newchild(lnblocks2) = .True. ! this is a new child

               parent(1,lnblocks2) = i ! child j's parent
               parent(2,lnblocks2) = mype ! child j's parent's processor

               bnd_box(:,:,lnblocks2) = bnd_box(:,:,i)

!--------------ordering of children in space is morton ordered !
               hh(:) = bsize(:,i)/4.
               If (j == 1) Then
                  hx = -1.
                  hy = -1.
                  hz = -1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               Else If (j == 2) Then
                  hx = 1.
                  hy = -1.
                  hz = -1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               Else If (j == 3) Then
                  hx = -1.
                  hy = 1.
                  hz = -1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               Else If (j == 4) Then
                  hx = 1.
                  hy = 1.
                  hz = -1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               Else If (j == 5) Then
                  hx = -1.
                  hy = -1.
                  hz = 1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               Else If (j == 6) Then
                  hx = 1.
                  hy = -1.
                  hz = 1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               Else If (j == 7) Then
                  hx = -1.
                  hy = 1.
                  hz = 1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               Else If (j == 8) Then
                  hx = 1.
                  hy = 1.
                  hz = 1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               End If

               If (ndim < 3) bnd_box(:,3,lnblocks2) = bnd_box(:,3,i)
               If (ndim < 2) bnd_box(:,2,lnblocks2) = bnd_box(:,2,i)

               which_child(lnblocks2) = j
                     
               coord(1,lnblocks2) = coord(1,i) + hx*hh(1)
               If (ndim >= 2) Then
                  coord(2,lnblocks2) = coord(2,i) + hy*hh(2)
               End If
               If (ndim == 3) Then
                  coord(3,lnblocks2) = coord(3,i) + hz*hh(3)
               End If

               bsize(:,lnblocks2) = bsize(:,i)/2.

               bflags(:,lnblocks2) = bflags(:,i)

               If (empty_cells) Then
                  If (empty(i) == 1) empty(lnblocks2)=1
               End If

            End Do  ! End Do j = 1, nchild

            If (empty_cells) Then
               If (empty(i) == 1) empty(i)=0
            End If

         End If  ! End If (refine(i))

       End Do  ! End Do i = 1, lnblocks

!------Connect neighbors of newly created sub-blocks
!
!                    4              6
!                    |              |
!  in x-y plane  1 - i - 2 ; in z   i
!                    |              |
!                    3              5
!
!-----In some cases, old data in the parts of the neigh array
!-----about to be Used by the new children can confUse the algorithm
!-----which searches for off processor neighbors. In particular, if
!-----an external boundary is not properly detected, Then a spurious
!-----neighbor link can be set up in section XXX below.
!-----Reinitializing these parts of neigh avoids the confusion.
      neigh(:,:,lnblocks+1:lnblocks2) = -1
      child(:,:,lnblocks+1:lnblocks2) = -1


!-----connect with siblings (which at this point are all on processor)
      Do i = 1,lnblocks
         If (refine(i)) Then
            Do j = 1,nchild ! cycle through children
               ichi = child(1,j,i)
               Do k = 1,nfaces
                  neigh(2,k,ichi) = mype
               End Do
            End Do

!-----------connect in x direction
            Do j = 2,nchild,2
               ichi = child(1,j,i)
               k = j - 1
               ineigh = child(1,k,i)
               neigh(1,1,ichi) = ineigh
            End Do
            Do j = 1,nchild-1,2
               ichi = child(1,j,i)
               k = j + 1
               ineigh = child(1,k,i)
               neigh(1,2,ichi) = ineigh
            End Do

!-----------connect in y direction
            If (ndim >= 2) Then
            Do j = 3,4
               ichi = child(1,j,i)
               k = j - 2
               ineigh = child(1,k,i)
               neigh(1,3,ichi) = ineigh
            End Do
            If (ndim == 3) Then
            Do j = 7,8
               ichi = child(1,j,i)
               k = j - 2
               ineigh = child(1,k,i)
               neigh(1,3,ichi) = ineigh
            End Do
            End If
            Do j = 1,2
               ichi = child(1,j,i)
               k = j + 2
               ineigh = child(1,k,i)
               neigh(1,4,ichi) = ineigh
            End Do
            If (ndim == 3) Then
            Do j = 5,6
               ichi = child(1,j,i)
               k = j + 2
               ineigh = child(1,k,i)
               neigh(1,4,ichi) = ineigh
            End Do
            End If
            End If  ! End If (ndim == 3)

!-----------connect in z direction
            If (ndim == 3) Then
            Do j = 5,8
               ichi = child(1,j,i)
               k = j - 4
               ineigh = child(1,k,i)
               neigh(1,5,ichi) = ineigh
            End Do
            Do j = 1,4
               ichi = child(1,j,i)
               k = j + 4
               ineigh = child(1,k,i)
               neigh(1,6,ichi) = ineigh
            End Do
            End If

         End If  ! End If (refine(i))
      End Do  ! End Do i = 1, lnblocks

!-----Now connect with off-processor neighbors by looking at neighbors of parent

!-----In order to establish pointers between children of newly refined blocks,
!-----each parent block must know if its neighbor is refined. In the next 
!-----section the refine flags are passed between neighbors. In a regular 
!-----cartesian grid this is a simple send-recv operation across every face. 
!-----For face 3 say (ie j=3) each block sends its refine value to its face 3 
!-----neighbor. Each block also issues a receive from its face 4 neighbor. 
!-----In this way sends and receives are automatiCally matched.
!-----In non-cartesian geometry with singular lines, eg a pole in spherical 
!-----coordinates, block faces touching the pole, touch their neighbor with 
!-----the same face number, eg face 3 touches face 3. Therefore in the 
!-----send/recv example above, polar blocks must issue a send but also a 
!-----receive for their face 3 neighbor, not their face 4 neighbor.
!-----To make this work, we have written the following section in two sweeps. 
!-----The first sweep operates as in the cartesian example, except blocks 
!-----with polar faces Do not issue any sends or receives for their polar 
!-----faces. On the second sweep, only the polar faces issue sends and receives.

      tneigh(:,:,:) = 0

!-----Send refine flags to neighbors
      Do j = 1,nfaces

      nsweep = 1

      If (spherical_pm) Then
         If (j == 3.or.j == 4) nsweep = 2
      End If

      Do isweep = 1,nsweep

         If (j == 1) Then
            jr = 2
            istart = 1
            iend = nchild-1
            istart2 = 2
            iend2 = nchild
            istep = 2
         ElseIf (j == 2) Then
            jr = 1
            istart = 2
            iend = nchild
            istart2 = 1
            iend2 = nchild-1
            istep = 2
         ElseIf (j == 3) Then
            jr = 4               
            istart = 1
            iend = 2
            istart2 = 3
            iend2 = 4
            istep = 1
         ElseIf (j == 4) Then
            jr = 3
            istart = 3
            iend = 4
            istart2 = 1
            iend2 = 2
            istep = 1
         ElseIf (j == 5) Then
            jr = 6
            istart = 1
            iend = 4
            istart2 = 5
            iend2 = 8
            istep = 1
         ElseIf (j == 6) Then
            jr = 5
            istart = 5
            iend = 8
            istart2 = 1
            iend2 = 4
            istep = 1
         End If

         refine_neigh(:) = .False.
         refine(lnblocks+1:maxblocks_tr) = .False.

         nrecv = 0
         Do i = 1,lnblocks

!--------Do receive from opposite index face, unless opposite face is
!--------a theta face in spherical coords
         lrecv=.True.

         If (isweep == 1) Then

         If (spherical_pm) Then
!---------In first sweep switch off any receives which might be generated at 
!---------polar faces
          If (j == 3.and.abs(bnd_box(2,2,i)-pi) < eps) lrecv=.False.
          If (j == 4.and.abs(bnd_box(1,2,i)) < eps) lrecv=.False.
         End If

         If (lrecv) Then

            If (neigh(1,jr,i) > 0) Then
               If (neigh(2,jr,i).ne.mype) Then
                  nrecv = nrecv + 1
                  Call MPI_IRECV(refine_neigh(i),                      & 
                                 1,                                    & 
                                 MPI_LOGICAL,                          & 
                                 neigh(2,jr,i),                        & 
                                 neigh(1,jr,i),                        &
                                 MPI_COMM_WORLD,                       & 
                                 reqr(nrecv),                          & 
                                 ierr)
               End If
            End If

         End If  ! End If (lrecv)

         End If  ! End If (isweep == 1)

         If (spherical_pm) Then
         If (isweep == 2) Then
!--------This section executes receives for polar faces only.
            lrecv = .False.
            If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) lrecv=.True.
            If (j == 3.and.abs(bnd_box(1,2,i)) < eps) lrecv=.True.
            If (abs(bnd_box(1,2,i)) < eps.and.                            &  
                abs(bnd_box(2,2,i)-pi) < eps) Then
             lrecv=.True.
            End If
            If (lrecv) Then
            If (neigh(1,j,i) > 0) Then
               If (neigh(2,j,i).ne.mype) Then
                  nrecv = nrecv + 1
                  Call MPI_IRECV(refine_neigh(i),                      & 
                                 1,                                    & 
                                 MPI_LOGICAL,                          & 
                                 neigh(2,j,i),                         & 
                                 neigh(1,j,i),                         &
                                 MPI_COMM_WORLD,                       & 
                                 reqr(nrecv),                          & 
                                 ierr)
               End If
            End If  ! End If (neigh(1,j,i) > 0)
            End If  ! End If (lrecv)
           End If  ! End If (isweep == 2)
         End If  ! End If (spherical_pm)
         End Do  ! End Do i = 1, lnblocks

         nsend = 0
         Do i = 1,lnblocks

         lsend = .True.

         If (isweep == 1) Then
            If (spherical_pm) Then
!-----------In first sweep switch off any sends which might be generated at 
!-----------polar faces
            If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) lsend=.False.
            If (j == 3.and.abs(bnd_box(1,2,i)) < eps) lsend=.False.
            End If

            If (lsend) Then
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            If (ineigh >= 1) Then
               If (ineigh_proc.ne.mype) Then
                  nsend = nsend + 1
                  Call MPI_SSEND(refine(i),                            & 
                                 1,                                    & 
                                 MPI_LOGICAL,                          & 
                                 ineigh_proc,                          & 
                                 i,                                    &
                                 MPI_COMM_WORLD,                       & 
                                 ierr)
               Else
                  If (lsend) refine_neigh(ineigh) = refine(i)
               End If
            End If
            End If  ! End If (lsend)
         End If  ! End If (isweep == 1)

         If (isweep == 2) Then
            If (spherical_pm) Then
            lsend = .False.
!-----------In second sweep switch on sends associated with polar faces only
            If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) lsend=.True.
            If (j == 3.and.abs(bnd_box(1,2,i)) < eps) lsend=.True.
            End If

            If (lsend) Then
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            If (ineigh >= 1) Then
               If (ineigh_proc.ne.mype) Then
                  nsend = nsend + 1
                  Call MPI_SSEND(refine(i),                            & 
                                 1,                                    & 
                                 MPI_LOGICAL,                          & 
                                 ineigh_proc,                          & 
                                 i,                                    &
                                 MPI_COMM_WORLD,                       & 
                                 ierr)
               Else
                  If (lsend) refine_neigh(ineigh) = refine(i)
               End If
            End If  ! End If (ineigh >= 1)
            End If  ! End If (lsend)
         End If  ! End If (isweep == 2)

         End Do  ! End Do i = 1, lnblocks

         If (nrecv > 0) Then
            Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         End If
 100     Continue

         Do jj = istart,iend,istep
            
         If (j == 1) Then
            k = jj+1         ! child no. of neighbor which lies on border
         ElseIf (j == 2) Then
            k = jj-1
         ElseIf (j == 3) Then
            k = jj+2
         ElseIf (j == 4) Then
            k = jj-2
         ElseIf (j == 5) Then
            k = jj+4
         ElseIf (j == 6) Then
            k = jj-4
         End If

         nrecv = 0
         Do i = 1,lnblocks
         kr = k
            If (refine(i)) Then

              If (isweep == 1) Then
              lrecv = .True.

              If (spherical_pm) Then
!--------------In first sweep switch off any receives which might be generated 
!--------------at polar faces
               If (j == 3.and.abs(bnd_box(1,2,i)) < eps) lrecv=.False.
               If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) lrecv=.False.
              End If

              ElseIf (isweep == 2) Then
              lrecv = .False.

              If (spherical_pm) Then
!--------------In first sweep switch off any receives which might be 
!--------------generated at polar faces
               If (j == 3.and.abs(bnd_box(1,2,i)) < eps) lrecv=.True.
               If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) lrecv=.True.
              End If

              End If  ! End If (isweep == 1)

              If (lrecv) Then

               ineigh = neigh(1,j,i) ! neighbor of parent
               ineigh_proc = neigh(2,j,i) ! neighbor of parent
               If (ineigh > 0) Then
                                       ! 'i' on left side
                  ichi = child(1,jj,i) ! child 'jj' of parent 'i'
!-----------------fetch child no. k of ineigh from ineigh_proc
!-----------------and store in neighbor 1 of child jj (ichi)
                  If (ineigh_proc.ne.mype) Then
                     nrecv = nrecv + 1
                     Call MPI_IRECV(neigh(1,j,ichi),                   & 
                                    2,                                 & 
                                    MPI_INTEGER,                       & 
                                    neigh(2,j,i),                      & 
                                    neigh(1,j,i),                      & 
                                    MPI_COMM_WORLD,                    & 
                                    reqr(nrecv),                       & 
                                    ierr)
                  Else

                  If (spherical_pm) Then
                     If (isweep == 2) Then
!--------------------if this child is at the pole Then modify the index of 
!--------------------the child which should send back.
!--------------------(ie if spherical pole block Then k=jj in next line)
                     If ( j == 3                                      & 
                         .and.abs(bnd_box(1,2,i)) < eps) Then
                       kr = jj
                     End If
                     If ( j == 4                                       & 
                         .and.abs(bnd_box(2,2,i)-pi) < eps) Then
                       kr = jj
                     End If
                     End If  ! End If (isweep == 2)
                   End If  ! End If (spherical_pm)
                   neigh(:,j,ichi) = child(:,kr,ineigh)
                  End If  ! End If (ineigh_proc.ne.mype)
               End If  ! End If (ineigh > 0)
              End If  ! End If (lrecv)
            End If  ! End If (lrefine(i)
         End Do  ! End Do i = 1, lnblocks

         nsend = 0
         Do i = 1,lnblocks
            kr = k

            If (refine_neigh(i)) Then

            If (isweep == 1) Then

            If (spherical_pm) Then
             lsend = .True.
!------------In first sweep switch off sends associated with polar faces
             If (j == 4.and.abs(bnd_box(1,2,i)) < eps) lsend=.False.
             If (j == 3.and.abs(bnd_box(2,2,i)-pi) < eps) lsend=.False.
            End If

            If (lsend) Then

               ineigh = neigh(1,jr,i) ! right neighbor of parent
               ineigh_proc = neigh(2,jr,i)
               If (ineigh > 0) Then
                  ! send child no. jj of block i to neigh on right
                  If (ineigh_proc.ne.mype) Then
                     nsend = nsend + 1
                     Call MPI_SSEND(child(1,kr,i),                     & 
                          2,                                           & 
                          MPI_INTEGER,                                 & 
                          neigh(2,jr,i),                               & 
                          i,                                           & 
                          MPI_COMM_WORLD,                              & 
                          ierr)
                  End If
               End If  ! End If (ineigh > 0)

             End If  ! End If (lsend)
             End If  ! End If (isweep == 1)


             If (isweep == 2) Then

             If (spherical_pm) Then
              lsend = .False.
!-------------In second sweep switch on sends associated with polar faces only
              If (j == 3.and.abs(bnd_box(1,2,i)) < eps) lsend=.True.
              If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) lsend=.True.

              If (lsend) Then

               jr0 = j
               ineigh = neigh(1,jr0,i) ! right neighbor of parent
               ineigh_proc = neigh(2,jr0,i)
               If (ineigh > 0 .and. j == jr0) Then
!-----------------send child no. jj of block i to neigh on right
                  If (ineigh_proc.ne.mype) Then
                     nsend = nsend + 1

!--------------------if spherical pole block Then kr=jj in next line
                     If (j == 3                                        & 
                         .and.abs(bnd_box(1,2,i)) < eps) Then
                       kr = jj
                     End If
                     If (j == 4                                        & 
                         .and.abs(bnd_box(2,2,i)-pi) < eps) Then
                       kr = jj
                     End If
                     Call MPI_SSEND(child(1,kr,i),                     & 
                          2,                                           & 
                          MPI_INTEGER,                                 & 
                          neigh(2,jr0,i),                              & 
                          i,                                           & 
                          MPI_COMM_WORLD,                              & 
                          ierr)
                  End If  ! End If (ineigh_proc.ne.mype)
               End If  ! End If (ineigh > 0 .and. j == jr0)
              End If  ! End If (lsend)
             End If  ! End If (spherical_pm)
            End If  ! End If (isweep == 2)
            End If  ! End If (refine_neigh(i)
         End Do  ! End Do i = 1, lnblocks
         
         If (nrecv > 0) Then
            Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         End If

         End Do  ! End Do jj = istart,iend,istep

         If (j == 3.and.istart == 1.and.ndim == 3) Then
            istart = 5
            iend = 6
            istart2 = 7
            iend2 = 8
            go to 100
         ElseIf (j == 4.and.istart == 3.and.ndim == 3) Then
            istart = 7
            iend = 8
            istart2 = 5
            iend2 = 6
            go to 100
         End If

!--------count no. of sends on each proc. to all other procs
         nsend_to_proc(:) = 0
         Do i = lnblocks+1,lnblocks2
            If (newchild(i)) Then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               If (ineigh >= 1) Then
                  If (ineigh_proc.ne.mype) Then
                     nsend_to_proc(ineigh_proc) =                      & 
                          nsend_to_proc(ineigh_proc) + 1
                  End If
               End If
            End If
         End Do

!--------collect data for `this' proc from other procs so that
!--------the total no. of receives to post can be computed
!--------(Change made by M. Zingale and J. Dursi)
         nrecv = 0
!  print *, 'CEG - Allreduce in mpi_amr_refine_blocks'
         Call MPI_ALLREDUCE(nsend_to_proc, nrecv_pack, nprocs,         & 
              MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         nrecv = nrecv_pack(MyPE)

         Do i = 1,nrecv
            Call MPI_IRECV(kk(i),                                      & 
                           1,                                          & 
                           MPI_INTEGER,                                & 
                           MPI_ANY_SOURCE,                             & 
                           MPI_ANY_TAG,                                & 
                           MPI_COMM_WORLD,                             & 
                           reqr(i),                                    & 
                           ierr)
         End Do  

!--------NOW new children send to their new neighbors so that they can be set
         nsend = 0
         Do i = lnblocks+1,lnblocks2
            If (newchild(i)) Then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               If (ineigh >= 1) Then
                  If (ineigh_proc.ne.mype) Then
                     nsend = nsend + 1
                     idummy_array(1) = ineigh
                     Call MPI_SSEND(idummy_array(1),                   & 
                                    1,                                 & 
                                    MPI_INTEGER,                       & 
                                    ineigh_proc,                       & 
                                    i,                                 & 
                                    MPI_COMM_WORLD,                    & 
                                    ierr)
                  Else
                     tneigh(1,jr,ineigh) = i
                     tneigh(2,jr,ineigh) = mype
                  End If
               End If
            End If
         End Do

!--------it seems to work only if tneigh is used ????
         If (nrecv > 0) Then
            Call MPI_WAITALL(nrecv,reqr,statr,ierr)
            Do i = 1,nrecv
               tneigh(1,jr,kk(i)) = statr(MPI_TAG,i)
               tneigh(2,jr,kk(i)) = statr(MPI_SOURCE,i)
            End Do
         End If
!CEG tried removing but it stopped it all working
         Call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! NEEDED ????

      End Do  ! End Do isweep = 1,nsweep

      End Do  ! End Do j = 1, nfaces

                     
      Do i = 1,lnblocks2
         Do j = 1,nfaces
            If (tneigh(1,j,i) > 0.and..not.newchild(i).and.            & 
                tneigh(1,j,i).ne.neigh(1,j,i)) Then
               neigh(1,j,i) = tneigh(1,j,i)
               neigh(2,j,i) = tneigh(2,j,i)
          
            End If
         End Do
      End Do
      

      lnblocks = lnblocks2

!-----reset node types
      Do i = 1,lnblocks
         If (newchild(i)) nodetype(i) = 1
         If (refine(i)) Then
            nodetype(i) = 2
         End If
      End Do

      nrecv = 0
      nodetype_chi(:,1:lnblocks) = 0
      Do i = 1,lnblocks
         Do j = 1,nchild
            If (child(1,j,i) > 0) Then
            If (child(2,j,i).ne.mype) Then
               nrecv = nrecv + 1
               Call MPI_IRECV(nodetype_chi(j,i),                       & 
                              1,                                       & 
                              MPI_INTEGER,                             & 
                              child(2,j,i),                            & 
                              child(1,j,i),                            & 
                              MPI_COMM_WORLD,                          & 
                              reqr(nrecv),                             & 
                              ierr)
            Else
               nodetype_chi(j,i) = nodetype(child(1,j,i))
            End If
            End If
         End Do
      End Do  

      nsend = 0
      Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
         If (parent(2,i).ne.mype) Then
            nsend = nsend + 1
            Call MPI_SSEND(nodetype(i),                                & 
                           1,                                          & 
                           MPI_INTEGER,                                & 
                           parent(2,i),                                &  
                           i,                                          & 
                           MPI_COMM_WORLD,                             & 
                           ierr)
         End If
         End If
      End Do

      If (nrecv > 0) Then
         Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If
      Do i = 1,lnblocks
         n3 = 0
         Do j = 1,nchild
            If (nodetype_chi(j,i).ne.1) Then
               n3 = n3 + 1
            End If
            If (n3 == nchild) nodetype(i) = 3
         End Do
      End Do 

!-----Now set neighbor pointers of new children if they lie on a boundary
      Do i = 1,lnblocks

         If (newchild(i)) Then

!-----------fetch i's parent
            ipar = parent(1,i)
            ipar_proc = parent(2,i)

            Do j = 1,nfaces

               If (neigh(1,j,i) <= -1) Then ! this neighbor may be on a border

!-----------------fetch i's parent's neighbor j
!-----------------Here we know that all new children on the same processor as
!-----------------the parent so we don't need any communications !!!

                  ineigh = neigh(1,j,ipar)
               
!-----------------if parent's neighbor is lt -1 Then i's neighbor is
!-----------------also on the domain border and is set to the parent's
!-----------------value

                  If (ineigh <= -20) neigh(1,j,i) = ineigh

               End If

            End Do

         End If  ! End If (newchild(i))

      End Do  ! End Do i = 1, lnblocks

!-----reset refine flags
      Do i = 1,maxblocks_tr
         refine(i) = .False.
      End Do

! CEG deallocate memory
  deallocate(nsend_to_proc)
  deallocate(nrecv_pack)
  deallocate(reqr)
  deallocate(kk)
  deallocate(nodetype_chi)
  deallocate(statr)
  deallocate(tneigh)
  deallocate(refine_neigh)
      Return
      End Subroutine amr_refine_blocks

