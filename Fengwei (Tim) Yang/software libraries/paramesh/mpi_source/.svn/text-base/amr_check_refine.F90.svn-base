!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_check_refine
!! NAME
!!
!!   amr_check_refine
!! 
!! SYNOPSIS
!!
!!   call amr_refine_blocks (nprocs, mype, icontinue)
!!
!!   call amr_refine_blocks (integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: nprocs
!!     The number of processors being used.
!!
!!   integer, intent(in) :: mype     
!!     The calling processor number.
!!
!!   integer, intent(in) :: icontinue     
!!     A variable controling algorithm behavior.
!!
!! INCLUDES
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
!!   morton_neighbors
!!
!! RETURNS
!!
!!   This routine does not return anything.  Upon return additional refine flags
!!   may be set to .true. to ensure a valid mesh.
!!
!! DESCRIPTION
!!
!!   This routine performs checks on the list of refine flags to make sure that
!!   a valid mesh results when the new, child blocks are added.  It makes sure that
!!   a jump in refinement of no more than one level is ensured.
!!  
!!   The routine is called only my amr_refine_derefine and a user's application should
!!   not need to call this routine directly.
!!
!! AUTHORS
!!
!!   Kevin Olson (1999)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_check_refine(nprocs,mype,icontinue)

      use paramesh_dimensions
      use physicaldata
      use tree
      use constants

      implicit none

      integer, intent(in)    :: nprocs,mype
      integer, intent(in)    :: icontinue

      include 'mpif.h'

! local variables ----------------------------------------------------------

      integer :: i,j,ipar,ipar_proc
      integer :: ineigh,ineigh_proc
      integer :: ix, iy, iz, ix_n, iy_n, iz_n
      integer :: ierr,jr
      integer :: nsend,nrecv
!      integer :: reqr(maxblocks_tr)
      integer,allocatable :: reqr(:)
!      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer,allocatable :: statr(:,:)

      logical :: repeat,repeat_t
!      logical :: ref_testt(2,1+k2d,1+k3d,maxblocks_tr)
!      logical :: ref_testts(2,1+k2d,1+k3d,maxblocks_tr)
!      logical :: ref_test(2,1+k2d,1+k3d,maxblocks_tr)
      logical,allocatable :: ref_testt(:,:,:,:)
      logical,allocatable :: ref_testts(:,:,:,:)
      logical,allocatable :: ref_test(:,:,:,:)
      logical :: lreduce_datain(1), lreduce_dataout(1)

      real, parameter :: eps = 1.e-10
      integer :: jr0
      logical :: lrecv,lsend

! --------------------------------------------------------------------------
! CEG allocate memory properly
  allocate(statr(MPI_STATUS_SIZE,maxblocks_tr))
  allocate(reqr(maxblocks_tr))
  allocate(ref_testt(2,1+k2d,1+k3d,maxblocks_tr))
  allocate(ref_testts(2,1+k2d,1+k3d,maxblocks_tr))
  allocate(ref_test(2,1+k2d,1+k3d,maxblocks_tr))

#ifdef DEBUG_FLOW_TRACE
         write(*,*) 'pe ',mype,' entered amr_check_refine'
#endif /* DEBUG_FLOW_TRACE */


! Turn off refine flags of non-leaf blocks
      do i = 1,lnblocks
         if (nodetype(i).gt.1.and.refine(i)) refine(i) = .false.
      end do


! If no leaf blocks were marked for refinement then return.
 21   if(icontinue.eq.0) return


!!!!!!!!!!      
! CHECK FOR neighboring blocks which differ by more than 1 level of refinement
      
      do i = 1,lnblocks
        do iz = 1,1+k3d
          do iy = 1,1+k2d
            do ix = 1,2
              ref_test(ix,iy,iz,i) = .FALSE.
            end do
          end do
        end do
      end do

#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' step 1', & 
     &     ' refine(1:lnblocks) ',refine(1:lnblocks)
#endif /* DEBUG */
! Step 1 - set up an array, ref_test, which contains the 
!          logical refine flag for each child of every local block

      nrecv = 0
      do i = 1,lnblocks
         do j = 1,nchild
            if (child(1,j,i).gt.0) then

              if (j.eq.1) then
                 ix = 1
                 iy = 1
                 iz = 1
              elseif (j.eq.2) then
                 ix = 2
                 iy = 1
                 iz = 1
              elseif (j.eq.3) then
                 ix = 1
                 iy = 2
                 iz = 1
              elseif (j.eq.4) then
                 ix = 2
                 iy = 2
                 iz = 1
              elseif (j.eq.5) then
                 ix = 1
                 iy = 1
                 iz = 2
              elseif (j.eq.6) then
                 ix = 2
                 iy = 1
                 iz = 2
              elseif (j.eq.7) then
                 ix = 1
                 iy = 2
                 iz = 2
              elseif (j.eq.8) then
                 ix = 2
                 iy = 2
                 iz = 2
              endif

               if (child(2,j,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_logical_IRECV(ref_test(ix,iy,iz,i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           child(2,j,i), & 
     &                           child(1,j,i), & 
     &                           MPI_COMM_WORLD, & 
     &                           reqr(nrecv), & 
     &                           ierr)
               else
                  ref_test(ix,iy,iz,i) = refine(child(1,j,i))
               end if
            end if
         end do
      end do
           
      nsend = 0
      do i = 1,lnblocks
            ipar = parent(1,i)
            ipar_proc = parent(2,i)
            if (ipar.ge.1) then
               if (ipar_proc.ne.mype) then
                  nsend = nsend + 1
                  call MPI_logical_SSEND(refine(i), & 
     &                           1, & 
     &                           MPI_LOGICAL, & 
     &                           ipar_proc, & 
     &                           i, & 
     &                           MPI_COMM_WORLD, & 
     &                           ierr)
               end if
            end if
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if


! Step 2 - now cycle over all faces, and for each set up the 
!          array, ref_testt, which contains the logical refine flag
!          for each child of the neighbor block across this face.

      do j = 1,nfaces     

#ifdef DEBUG
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! in DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' step 2',' jface ', & 
     &            j
#endif /* DEBUG */
         do i = 1,lnblocks
            ref_testt(:,:,:,i) = .FALSE.

      if(spherical_pm) then
            ref_testts(:,:,:,i) = .FALSE.
      endif

         end do

         if (j.eq.1) then
            jr = 2
         elseif (j.eq.2) then
            jr = 1
         elseif (j.eq.3) then
            jr = 4
         elseif (j.eq.4) then
            jr = 3
         elseif (j.eq.5) then
            jr = 6
         elseif (j.eq.6) then
            jr = 5
         end if
        
         nrecv = 0
         do i = 1,lnblocks

! do receive from opposite index face, unless opposite face is
! a theta face in spherical coords
         lrecv=.true.

      if(spherical_pm) then
         if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
         if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
      endif

            if (lrecv) then

            if (neigh(1,jr,i).gt.0) then
               if (neigh(2,jr,i).ne.mype) then
                  nrecv = nrecv + 1
#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' blk ',i, & 
     &           ' nrecv ',nrecv, & 
     &           ' source pe ',neigh(2,jr,i), & 
     &           ' posting n recv testt tag ',neigh(1,jr,i),jr
#endif /* DEBUG */
                  call MPI_logical_IRECV(ref_testt(1,1,1,i), & 
     &                           2**ndim, & 
     &                           MPI_LOGICAL, & 
     &                           neigh(2,jr,i), & 
     &                           neigh(1,jr,i), & 
     &                           MPI_COMM_WORLD, & 
     &                           reqr(nrecv), & 
     &                           ierr)
               end if
            end if

            end if


      if(spherical_pm) then
! additional receives needed by polar blocks
       lrecv=.true.
       jr0 = jr
       if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) jr0 = 3
       if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
       if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) jr0 = 4
       if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
       if(abs(bnd_box(1,2,i)).lt.eps.and. & 
     &        abs(bnd_box(2,2,i)-pi).lt.eps) then
             lrecv=.true.
        endif
            if (lrecv .and. jr0.eq.j) then
            if (neigh(1,jr0,i).gt.0) then
               if (neigh(2,jr0,i).ne.mype) then
                  nrecv = nrecv + 1
#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' blk ',i, & 
     &           ' nrecv ',nrecv, & 
     &           ' source pe ',neigh(2,jr0,i), & 
     &           ' posting recv testts tag ',neigh(1,jr0,i),jr0
#endif /* DEBUG */

! Note the array ref_testts is larger than is really needed. More careful
! coding may reduce the size.
                  call MPI_logical_IRECV(ref_testts(1,1,1,i), & 
     &                           2**ndim, & 
     &                           MPI_LOGICAL, & 
     &                           neigh(2,jr0,i), & 
     &                           neigh(1,jr0,i), & 
     &                           MPI_COMM_WORLD, & 
     &                           reqr(nrecv), & 
     &                           ierr)
               end if
            end if
            end if
         end if  ! end if (spherical_pm

         end do
         
         nsend = 0
         do i = 1,lnblocks

           lsend = .true.

           if(lsend) then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               if (ineigh.ge.1) then
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1
#ifdef DEBUG
      write(*,*) 'amr_refine_block : proc ',mype,' blk ',i, & 
     &          ' nsend ',nsend, & 
     &          ' dest pe ',ineigh_proc, & 
     &          ' posting send ref_test tag ',i,' neigh() ',neigh(:,j,i)
#endif /* DEBUG */
                     call MPI_logical_SSEND(ref_test(1,1,1,i), & 
     &                              2**ndim, & 
     &                              MPI_LOGICAL, & 
     &                              ineigh_proc, & 
     &                              i, & 
     &                              MPI_COMM_WORLD, & 
     &                              ierr)
                  else


!                     ref_testt(:,:,:,ineigh) = ref_test(:,:,:,i)
!start of new test section
         lrecv=.true.
      if(spherical_pm) then
         if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
         if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
      endif
            if (lrecv) then
            if (neigh(1,j,i).gt.0) then
               if (neigh(2,j,i).eq.mype) then
                     ref_testt(:,:,:,neigh(1,j,i)) = ref_test(:,:,:,i)
               end if
            end if
            end if

      if(spherical_pm) then
! additional receives needed by polar blocks
       lrecv=.true.
       jr0 = jr
       if(j.eq.3.and.abs(bnd_box(1,2,i)).lt.eps) jr0 = 3
       if(j.eq.4.and.abs(bnd_box(1,2,i)).lt.eps) lrecv=.false.
       if(j.eq.4.and.abs(bnd_box(2,2,i)-pi).lt.eps) jr0 = 4
       if(j.eq.3.and.abs(bnd_box(2,2,i)-pi).lt.eps) lrecv=.false.
       if(abs(bnd_box(1,2,i)).lt.eps.and. & 
     &        abs(bnd_box(2,2,i)-pi).lt.eps) then
             lrecv=.true.
       endif
            if (lrecv .and. jr0.eq.j) then
            if (neigh(1,jr0,i).gt.0) then
               if (neigh(2,jr0,i).eq.mype) then
! Note the array ref_testts is larger than is really needed. More careful
! coding may reduce the size.
                  ref_testts(:,:,:,i) = ref_test(:,:,:,i)
               end if
            end if
            end if
      endif

!end of new test section


                  end if
               end if

           end if

         end do
#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' waiting jface ',j, & 
     &             ' testt nrecv ',nrecv,' nsend ',nsend
#endif /* DEBUG */
         
         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if

#ifdef DEBUG
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! in DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' step 3',' jface ', & 
     &            j
#endif /* DEBUG */

! Step 3 - if any elements of ref_testt are true, then make sure
!          that the elements of ref_test across face j are also
!          set to true. (Note: j here is the face as indexed by
!          the neighbor of block i.)
         
         do i = 1,lnblocks
           
            if ( any(ref_testt(1:2,:,:,i)) & 
     &         .or.(spherical_pm) & 
     &         )  then
! CEG added extra copy in for non-terminated spherical test
            if ( spherical_pm &
     &         )   write(*,*) 'spherical_pm check taken out of mpi_source/amr_check_refine.F90'
!.and. any(ref_testts(:,:,:,i))) & 

               
               if (j.eq.1) then
                  
                  ix = 2
                  ix_n = 1
                  do iz = 1,1+k3d
                     do iy = 1,1+k2d
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix_n,iy,iz,i)
                     end do
                  end do
                  
               elseif (j.eq.2) then
                  
                  ix = 1
                  ix_n = 2
                  do iz = 1,1+k3d
                     do iy = 1,1+k2d
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix_n,iy,iz,i)
                     end do
                  end do
              
               elseif (j.eq.3) then
                  
                  iy = 2
                  iy_n = 1          
                  do iz = 1,1+k3d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy_n,iz,i) & 
     &           .or. (spherical_pm .and. ref_testts(ix,iy,iz,i))
                     end do
                  end do
                  
               elseif (j.eq.4) then
                  
                  iy = 1
                  iy_n = 2
                  do iz = 1,1+k3d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy_n,iz,i) & 
     &           .or. (spherical_pm .and. ref_testts(ix,iy,iz,i))
                     end do
                  end do
                  
               elseif (j.eq.5) then
                  
                  iz = 2
                  iz_n = 1
                  do iy = 1,1+k2d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy,iz_n,i)
                     end do
                  end do
                  
               elseif (j.eq.6) then
                  
                  iz = 1
                  iz_n = 2
                  do iy = 1,1+k2d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) =  & 
     &                       ref_test(ix,iy,iz,i) .or.  & 
     &                       ref_testt(ix,iy,iz_n,i)
                     end do
                  end do
                  
               end if

            end if

            
         end do
      end do

#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,' step 4'
#endif /* DEBUG */
! Step 4      
! SET REFINE FLAGS BASED ON ref_test flags and repeat refinement process
! if necessary


      repeat = .FALSE.
      do i = 1,lnblocks

         do iz = 1,1+k3d
           do iy = 1,1+k2d
             do ix = 1,2
           
               if (ref_test(ix,iy,iz,i).and.nodetype(i).eq.1 & 
     &               .and..not.refine(i)) then
                 repeat = .TRUE.
                 refine(i) = .TRUE.
                 derefine(i) = .FALSE.
               end if

             end do
           end do
         end do
               
      end do

! cycle through all processors to see if any repeat, if so all repeat
! this should be done via a scan function

      lreduce_datain(1) = repeat
      call mpi_logical_allreduce (lreduce_datain,lreduce_dataout, & 
     &                    1,MPI_LOGICAL, & 
     &                    MPI_LOR,MPI_COMM_WORLD,ierr)
      repeat_t = lreduce_dataout(1)
      repeat = repeat_t

      if (repeat) then
         go to 21
      end if
#ifdef DEBUG_FLOW_TRACE
         write(*,*) 'pe ',mype,' exiting amr_check_refine'
#endif /* DEBUG_FLOW_TRACE */


  deallocate(statr)
  deallocate(reqr)
  deallocate(ref_testt)
  deallocate(ref_testts)
  deallocate(ref_test)

      return
      end subroutine amr_check_refine
