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

!#define DEBUG

      subroutine mpi_amr_singular_line (istrategy,nprocs)

!----------------------------------------------------------------------


! This routine manages the process of refinement in the vicinity of a singular line.
! The polar axis is an example of a singular line in spherical coordinates.
! A singular line is a line along a coordinate direction. Blocks which touch
! this line have zero face area there.

! We wish to ensure that we do not violate the 2 refinement level criterion at
! this line.

!----------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use mpi_morton
      use constants
      Use paramesh_comm_data
      use paramesh_interfaces, only : amr_q_sort_real

      implicit none

      integer, intent(in) :: istrategy,nprocs

      include 'mpif.h'


      real,   allocatable,dimension(:)   :: singular_line_locn
      integer,allocatable,dimension(:)   :: singular_line_level
      real,   allocatable,dimension(:)   :: gsingular_line_locn
      integer,allocatable,dimension(:)   :: gsingular_line_level
      integer,allocatable,dimension(:)   :: gstack
      integer,allocatable,dimension(:)   :: displs
      integer,allocatable,dimension(:)   :: indx
      real,   allocatable,dimension(:)   :: temp

      real    :: eps,accuracy
      real    :: test,test1,test2,cposition,rcoord
      integer :: istack,istack_tot,jstack,no_of_singular_pts,lb,ierror
      integer :: iloop,i
      integer :: mype, ierr
 
!----------------------------------------------------------------------
      accuracy = 10./10.**precision(accuracy)


      if(allocated(gstack)) deallocate(gstack)
      if(allocated(displs)) deallocate(displs)
      if(allocated(singular_line_locn)) deallocate(singular_line_locn)
      allocate(singular_line_locn(maxblocks))
      allocate(gstack(nprocs))
      allocate(displs(nprocs))

!----------------------------------------------------------------------

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)


#ifdef DEBUG
      write(*,*) '[',mype,' ] ', & 
     &     ' entered mpi_amr_singular_line: istrategy = ', & 
     &                istrategy
      call amr_flush(6)
#endif /* DEBUG */
      if(istrategy == 0) then
        return
      
      elseif(istrategy == 1) then


! Strategy 1 - simplest. Assume a uniform grid in the directions orthogonal
! to the line for the blocks touching the line.

!----------------------------------------------------------------------

! Construct a list of the block coordinate of any leaf blocks which are
! flagged for refinement and which touch the singular line.
! Then exchange these lists with other processors. If any other leaf blocks
! share the same singular coordinate with an entry in one of these lists,
! then reset the refine flag.

! Then repeat the same loop, but this time analyse the derefinement requests.


! loop through twice. On iloop=1, make adjustments to the refine array along
! the singular line. For iloop=2 make sure the derefine array is consistent.
      do iloop=1,2

      istack = 0

      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then

        if(spherical_pm) then
        eps = pi*accuracy
        if(iloop.eq.1) then
          if( refine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            istack = istack+1
            singular_line_locn(istack)    = coord(1,lb)     ! coord along radial singular line
            if(abs(bnd_box(2,2,lb)-pi) .lt. eps )            & ! if in southern hemi
     &           singular_line_locn(istack) = -coord(1,lb)  ! assign -ve sign to 
                                                            ! location
          endif
        elseif(iloop.eq.2) then
          if( .not.derefine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            istack = istack+1
            singular_line_locn(istack)    = coord(1,lb)     ! coord along radial singular line
            if(abs(bnd_box(2,2,lb)-pi) .lt. eps )            & ! if in southern hemi
     &           singular_line_locn(istack) = -coord(1,lb)  ! assign -ve sign to 
                                                            ! location
          endif
        endif
        endif

        if (cylindrical_pm) then
        eps = abs(grid_xmax-grid_xmin)*accuracy
        if(iloop.eq.1) then
          if( refine(lb) .and. & 
     &        (abs(bnd_box(1,1,lb)-0.) .lt. eps)) then
            istack = istack+1
            singular_line_locn(istack)    = coord(3,lb)     ! coord along axial singular line
          endif
        elseif(iloop.eq.2) then
          if( .not.derefine(lb) .and. & 
     &        (abs(bnd_box(1,1,lb)-0.) .lt. eps)) then
            istack = istack+1
            singular_line_locn(istack)    = coord(3,lb)     ! coord along axial singular line
          endif
        endif
        end if

      if (polar_pm) then
        eps = abs(grid_xmax-grid_xmin)*accuracy
        if(iloop.eq.1) then
          if( refine(lb) .and. & 
     &        abs(bnd_box(1,1,lb)-0.) .lt. eps ) then
            istack = istack+1
            singular_line_locn(istack)    = coord(1,lb)     ! coord at zero radius
          endif
        elseif(iloop.eq.2) then
          if( .not.derefine(lb) .and. & 
     &        abs(bnd_box(1,1,lb)-0.) .lt. eps ) then
            istack = istack+1
            singular_line_locn(istack)    = coord(1,lb)     ! coord at zero radius
          endif
        endif
        end if

      endif
      enddo


#ifdef DEBUG
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      write(*,*) '[',mype,' ] ', & 
     &      'mpi_amr_singular_line: istack = ',istack
      if(istack.gt.0) & 
     & write(*,*) '[',mype,' ] ', & 
     & 'mpi_amr_singular_line: pe  singular_line_locn ', & 
     &           singular_line_locn(1:istack)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */

! each proc gets the number of refining singular blocks on every other proc
      call mpi_allgather(istack,1,MPI_INTEGER,gstack,1,MPI_INTEGER, & 
     &                    MPI_COMM_WORLD,ierror)

      call comm_int_sum_to_all(istack_tot,istack)
#ifdef DEBUG
      write(*,*) '[',mype,' ] ', & 
     & 'mpi_amr_singular_line: ', & 
     &           ' gstack = ',gstack
      write(*,*) '[',mype,' ] ', & 
     &     'mpi_amr_singular_line:  ', & 
     &           ' istack_tot = ',istack_tot
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */

      if(istack_tot.gt.0) then


      if(allocated(indx)) deallocate(indx)
      allocate(indx(istack_tot))
      if(allocated(gsingular_line_locn)) deallocate(gsingular_line_locn)
      allocate(gsingular_line_locn(istack_tot))

      displs(1) = 0
      do i=2,nprocs
        displs(i) = displs(i-1) + gstack(i-1)
      enddo

      if(nprocs.gt.1) then
      call mpi_allgatherv(singular_line_locn,istack, & 
     &                    amr_mpi_real, & 
     &                    gsingular_line_locn,gstack, & 
     &                    displs, & 
     &                    amr_mpi_real, & 
     &                    MPI_COMM_WORLD,ierror)
      else
        if(istack.gt.0) & 
     &      gsingular_line_locn(1:istack) = singular_line_locn(1:istack) 
      endif

#ifdef DEBUG
      write(*,*) '[',mype,' ] ', & 
     &  'mpi_amr_singular_line: gsingular_line_locn ', & 
     &           gsingular_line_locn(1:istack_tot)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */


! locally sort gsingular_line_locn
      do i = 1, istack_tot
         indx(i) = i
      end do
      call amr_q_sort_real (gsingular_line_locn,istack_tot,ia = indx)

#ifdef DEBUG
      write(*,*) '[',mype,' ] ', & 
     &  'mpi_amr_singular_line: after sort  indx ', & 
     &                  indx(:istack_tot)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */
!
! apply permutation to gsingular_line_locn
      if(allocated(temp)) deallocate(temp)
      allocate(temp(istack_tot))
      temp(1:istack_tot) = gsingular_line_locn(1:istack_tot)
      do i=1,istack_tot
       gsingular_line_locn(i) = temp(indx(i))
      enddo

! compress list
       jstack = 0
       cposition = 0.
       do i = 1,istack_tot
         if( gsingular_line_locn(i).ne.cposition ) then
           jstack = jstack + 1
           gsingular_line_locn(jstack) = gsingular_line_locn(i)
           cposition = gsingular_line_locn(i)
         endif
       enddo
       no_of_singular_pts = jstack

#ifdef DEBUG
      write(*,*) '[',mype,' ] ', & 
     &        'mpi_amr_singular_line: after compress ', & 
     &          ' gsingular_line_locn ', & 
     &           gsingular_line_locn(:no_of_singular_pts)
      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif /* DEBUG */

! now for each singular block see if its singular_locn matches anything in this list
      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then
        if(spherical_pm) then
        if(iloop.eq.1) then
          if( .not.refine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            rcoord = coord(1,lb)
            if(abs(bnd_box(2,2,lb)-pi) .lt. eps ) rcoord = -rcoord
            do i=1,no_of_singular_pts
              test = abs(gsingular_line_locn(i)-rcoord)
              if(test.lt.eps) refine(lb)=.true.
            enddo
          endif
        elseif(iloop.eq.2) then
          if( derefine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            rcoord = coord(1,lb)
            if(abs(bnd_box(2,2,lb)-pi) .lt. eps ) rcoord = -rcoord
            do i=1,no_of_singular_pts
              test = abs(gsingular_line_locn(i)-rcoord)
              if(test.lt.eps) derefine(lb)=.false.
            enddo
          endif
        endif
        endif

      if (polar_pm) then
        if(iloop.eq.1) then
          if( .not.refine(lb) .and. & 
     &        abs(bnd_box(1,1,lb)-0.) .lt. eps ) then
            do i=1,no_of_singular_pts
              test = abs(gsingular_line_locn(i)-coord(1,lb))
              if(test.lt.eps) refine(lb)=.true.
            enddo
          endif
        elseif(iloop.eq.2) then
          if( derefine(lb) .and. & 
     &        abs(bnd_box(1,1,lb)-0.) .lt. eps ) then
            do i=1,no_of_singular_pts
              test = abs(gsingular_line_locn(i)-coord(1,lb))
              if(test.lt.eps) derefine(lb)=.false.
            enddo
          endif
        endif
        end if
      endif
      enddo


      endif                           ! end of istack_tot if test


      enddo                           ! end of iloop do loop


      call amr_flush(6)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!----------------------------------------------------------------------
      elseif(istrategy == 2) then
!----------------------------------------------------------------------
#ifdef NOTNOW
! Strategy 2 - more complex. Allow variation in resolution in the directions orthogonal
! to the line for the blocks touching the line.

!----------------------------------------------------------------------


! Construct a list of leaf block corners which touch the singular line

      if(allocated(gsingular_line_locn)) deallocate(gsingular_line_locn)
      if(allocated(gsingular_line_level))  & 
     &                                  deallocate(gsingular_line_level)
      allocate(gsingular_line_locn(2*maxblocks))
      allocate(gsingular_line_level(2*maxblocks))

! loop through twice. On iloop=1, make adjustments to the refine array along
! the singular line. For iloop=2 make sure the derefine array is consistent.
      do iloop=1,2
      istack = 0

      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then

        istack = 0
        if(spherical_pm) then
        eps = abs(grid_ymax-grid_ymin)*accuracy
        if(iloop.eq.1) then
          if( refine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            istack = istack+1
            singular_line_locn(istack:istack+1)    = bnd_box(:,1,lb)     ! coord along radial singular line
            singular_line_level(istack)    = lrefine(lb)+1         ! store refinement level for each node
            singular_line_level(istack+1)  = lrefine(lb)+1         ! store refinement level for each node
            istack = istack+1
          endif
        elseif(iloop.eq.2) then
          if( .not.derefine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            istack = istack+1
            singular_line_locn(istack:istack+1)    = bnd_box(:,1,lb)     ! coord along radial singular line
            singular_line_level(istack)    = lrefine(lb)+1         ! store refinement level for each node
            singular_line_level(istack+1)  = lrefine(lb)+1         ! store refinement level for each node
          endif
        endif
        endif

      if (cylindrical_pm) the
        eps = abs(grid_xmax-grid_xmin)*accuracy
        if(iloop.eq.1) then
          if( refine(lb) .and. & 
     &        (abs(bnd_box(1,1,lb)-0.) .lt. eps)) then
            istack = istack+1
            singular_line_locn(istack:istack+1)    = bnd_box(:,3,lb)     ! coord along radial singular line
            singular_line_level(istack)    = lrefine(lb)+1         ! store refinement level for each node
            singular_line_level(istack+1)  = lrefine(lb)+1         ! store refinement level for each node
            istack = istack+1
          endif
        elseif(iloop.eq.2) then
          if( .not.derefine(lb) .and. & 
     &        (abs(bnd_box(1,1,lb)-0.) .lt. eps)) then
            istack = istack+1
            singular_line_locn(istack:istack+1)    = bnd_box(:,3,lb)     ! coord along radial singular line
            singular_line_level(istack)    = lrefine(lb)+1         ! store refinement level for each node
            singular_line_level(istack+1)  = lrefine(lb)+1         ! store refinement level for each node
            istack = istack+1
          endif
        endif
        end if

      endif
      enddo

! each proc gets the number of singular nodes on every other proc
      call mpi_allgather(istack,1,MPI_INTEGER,gstack,1,MPI_INTEGER, & 
     &                    MPI_COMM_WORLD,ierror)

      call comm_int_sum_to_all(istack_tot,istack)

! if there are any singular nodes, provide their info to every processor
      if(istack_tot.gt.0) then

      if(allocated(gsingular_line_locn)) deallocate(gsingular_line_locn)
      allocate(gsingular_line_locn(istack_tot))
      if(allocated(gsingular_line_level))  & 
     &                                deallocate(gsingular_line_level)
      allocate(gsingular_line_level(istack_tot))

      displs(1) = 0
      do i=2,nprocs
        displs(i) = displs(i-1) + gstack(i-1)
      enddo
      call mpi_allgatherv(singular_line_locn,istack, & 
     &                    amr_mpi_real,	 & 
     &                    gsingular_line_locn,gstack, & 
     &                    displs, & 
     &                    amr_mpi_real,	 & 
     &                    MPI_COMM_WORLD,ierror)
    

      call mpi_allgatherv(singular_line_level,istack,MPI_INTEGER, & 
     &                    gsingular_line_level,gstack, & 
     &                    displs,MPI_INTEGER, & 
     &                    MPI_COMM_WORLD,ierror)



! locally sort gsingular_line_locn
      do i = 1, istack_tot
         gsingular_line_level(i) = i
      end do
      call amr_q_sort_real(gsingular_line_locn,istack_tot, & 
     &                                      ia = gsingular_line_level)

! compress list
! and retain the largest refinement level found for each singular node
       jstack = 0
       cposition = 0.
       do i = 1,istack_tot
         if( gsingular_line_locn(i).ne.cposition ) then
           jstack = jstack + 1
           gsingular_line_locn(jstack) = gsingular_line_locn(i)
           gsingular_line_level(jstack) = gsingular_line_level(i)
           cposition = gsingular_line_locn(i)
         else
           gsingular_line_level(jstack) = max(  & 
     &                                      gsingular_line_level(jstack) & 
     &                                     ,gsingular_line_level(i) )
         endif
       enddo
       no_of_singular_pts = jstack



! now for each singular block see if its singular_locn matches anything in this list
      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then
        if(spherical_pm) then
        if(iloop.eq.1) then
          if( .not.refine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            do i=1,no_of_singular_pts
              test1 = abs(gsingular_line_locn(i)-bnd_box(1,1,lb))
              test2 = abs(gsingular_line_locn(i)-bnd_box(2,1,lb))
              if(test1.lt.eps.or.test2.lt.eps) then
                if(lrefine(lb).lt.gsingular_line_level(i)-1) then
                  refine(lb)=.true.
                endif
              endif
            enddo
          endif
        elseif(iloop.eq.2) then
          if( derefine(lb) .and. & 
     &        (abs(bnd_box(1,2,lb)-0.) .lt. eps .or. & 
     &         abs(bnd_box(2,2,lb)-pi) .lt. eps ) ) then
            do i=1,no_of_singular_pts
              test1 = abs(gsingular_line_locn(i)-bnd_box(1,1,lb))
              test2 = abs(gsingular_line_locn(i)-bnd_box(2,1,lb))
              if(test1.lt.eps.or.test2.lt.eps) then
                if(lrefine(lb).lt.gsingular_line_level(i)) & 
     &                 derefine(lb)=.false.
              endif
            enddo
          endif
        endif
        endif
      endif
      enddo


      endif                           ! end of istack_tot if test

      enddo                           ! end of iloop do loop


      deallocate(gsingular_line_locn)
      deallocate(gsingular_line_level)

#endif /* NOTNOW */

! end of strategy 2
      endif                           ! istrategy if test
!----------------------------------------------------------------------

      if(allocated(singular_line_locn)) deallocate(singular_line_locn)
      if(allocated(singular_line_level)) deallocate(singular_line_level)
      if(allocated(gstack)) deallocate(gstack)
      if(allocated(displs)) deallocate(displs)
      if(allocated(temp)) deallocate(temp)

      return
      end subroutine mpi_amr_singular_line

