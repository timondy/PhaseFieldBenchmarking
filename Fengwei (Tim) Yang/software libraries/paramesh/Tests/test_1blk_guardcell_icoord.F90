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

!# define DEBUG

      program test_1blk_guardcell_icoord





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use workspace
      use io

      use paramesh_interfaces, only : comm_start, & 
     &                                amr_initialize, & 
     &                                amr_refine_derefine, & 
     &                                amr_1blk_copy_soln, & 
     &                                amr_guardcell, & 
     &                                amr_prolong, & 
     &                                amr_1blk_guardcell, & 
     &                                amr_1blk_guardcell_reset, & 
     &                                amr_1blk_restrict, & 
     &                                guardcell_test, & 
     &                                mesh_test, & 
     &                                amr_close

      use paramesh_mpi_interfaces, only :  & 
     &                                Mpi_morton_bnd, & 
     &                                mpi_amr_comm_setup, & 
     &                                mpi_amr_1blk_restrict, & 
     &                                mpi_morton_bnd_prolong

! Only required for programs in ./Tests
#include "test_defs.fh"

      include 'mpif.h'

      integer :: tag_offset,max_blks_sent
      integer nguard0
      integer nguard_work0
      integer :: three,four,five, six, lb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local amr variables
      integer nprocs,mype
      integer num_procs

      save mype

      logical, allocatable :: tnewchild(:)

!
! application specific variables

      real :: accuracy
      integer iopt,nlayers,icoord
      integer ierror_sum,ierror_tot
      logical lrefine_again,ltype2only
      logical lcc, lfc, lec, lnc, l_srl_only, ldiag

      logical :: lmpi,lnperm,ladvanceall
      integer :: errorcode,errcode, ierr

      logical :: lguard,lprolong,lflux,ledge,lrestrict
      logical :: lfulltree

      integer   ngmax , ngd
      logical, allocatable :: lmask(:,:,:)

      character (len=80) :: filename


      save ierror_sum,ierror_tot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call amr_initialize ()

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs
      ngmax = max(nguard,nguard_work)

      allocate(lmask(nxb+2*ngmax+1,nyb+(2*ngmax+1)*k2d, & 
     &     nzb+(2*ngmax+1)*k3d))

      allocate(tnewchild(maxblocks_tr))

      ierror_sum = 0
      ierror_tot = 0

      do iloop = 1,ndim
      icoord = iloop
      Write(*,*) 'starting loop : icoord = ',icoord
        
!      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)



       if(iloop.gt.1) then
!-------------------------------
!      call amr_initialize

! We would have liked to call amr_initialize here, but because it
! calls comm_start, we cannot call it inside a loop. So instead we
! reproduce the necessary lines from amr_initialize here.

! initialize tree data structure
        bsize(:,:) = -1.
        lrefine(:) = -1
        nodetype(:) = -1
        stay(:) = .TRUE.
        refine(:) = .FALSE.
        derefine(:) = .FALSE.
        parent(:,:) = -1
        child(:,:,:) = -1
        which_child(:) = -1
        coord(:,:) = -1.
        bnd_box(:,:,:) = -1.
        neigh(:,:,:) = -1
        empty(:) = 0
        bflags(:,:) = -1
        work_block(:) = 0.
        surr_blks(:,:,:,:,:) = -1

! initialize solution arrays
        if (nvar > 0) then
          unk(:,:,:,:,:) = 0.
        end if
        if (nvarcorn > 0) then
          unk_n(:,:,:,:,:) = 0.
        end if
        if (nvaredge > 0) then
           unk_e_x(:,:,:,:,:) = 0.
           unk_e_y(:,:,:,:,:) = 0.
           unk_e_z(:,:,:,:,:) = 0.
        end if
        if (nfacevar > 0) then
           facevarx(:,:,:,:,:) = 0.
           facevary(:,:,:,:,:) = 0.
           facevarz(:,:,:,:,:) = 0.
	end if	

        call amr_1blk_guardcell_reset

! Mark amr_gsurrounding_blks uncalled. This will be set to +1 if
! and when amr_gsurrounding_blks is called.
        gsurrblks_set = -1

        
!-------------------------------
       endif



#ifdef CONSERVE
      write(*,*) 'CONSERVE must not be defined for this test!'
      call amr_abort()
#endif

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      accuracy = 100./10.**precision(accuracy)

      print *,' nprocs = ',nprocs,mype

      if (diagonals) then
         write(*,*) 'diagonals on '
      end if

! set default value of dz and z0 to cater for 2D case.
      z0 = 0.
      dz = 0.


      ierror_sum = 0
      ierror_tot = 0

      iopt = 1
      nlayers = nguard
      if(mype.eq.0) write(*,*) 'nlayers = ',nlayers

!

! set a limit on the refinement level
      lrefine_max = 5
      lrefine_min = 1

      ax = 1.
      ay = 10.
      az = 100.
#ifdef TESTXDIR
      ax = 1.
      ay = 0.
      az = 0.
#endif
#ifdef TESTYDIR
      if (ndim >= 2) then
      ax = 0.
      ay = 1.
      az = 0.
      end if
#endif
#ifdef TESTZDIR
      if (ndim == 3) then
      ax = 0.
      ay = 0.
      az = 1.
      end if
#endif


! set the workspace array layer to be tested
      ioptw = 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up initial grid state.

      g_xmin = 0.
      g_xmax = 1.
      g_ymin = 0.
      g_ymax = 1.
      g_zmin = 0.
      g_zmax = 1.

! set up a single block covering the whole cubic domain
      lnblocks = 0
      if(mype.eq.0.) then
                lnblocks = 1
                bsize(:,1)=1.
                coord(:,1) = .5
                bnd_box(1,:,1) = g_xmin
                bnd_box(2,:,1) = g_xmax
                nodetype(1) = 1
                lrefine(1) = 1

                neigh(:,:,1) = -21

                refine(1)=.true.
      endif


      boundary_index = -21
! x boundaries
      boundary_box(1,2:3,1:2) = -1.e10
      boundary_box(2,2:3,1:2) =  1.e10
      boundary_box(1,1,1) = -1.e10
      boundary_box(2,1,1) = g_xmin
      boundary_box(1,1,2) = g_xmax
      boundary_box(2,1,2) = 1.e10
! y boundaries
      if(ndim.ge.2) then
      three = (2*k2d) + 1
      four  = three + k2d
      boundary_box(1,1,three:four) = -1.e10
      boundary_box(2,1,three:four) =  1.e10
      boundary_box(1,3,three:four) = -1.e10
      boundary_box(2,3,three:four) =  1.e10
      boundary_box(1,2,three) = -1.e10
      boundary_box(2,2,three) = g_ymin
      boundary_box(1,2,four) = g_ymax
      boundary_box(2,2,four) = 1.e10
      endif
! z boundaries
      if(ndim.eq.3) then
      five = (4*k3d) + 1
      six  = five + k3d
      boundary_box(1,1:2,five:six) = -1.e10
      boundary_box(2,1:2,five:six) =  1.e10
      boundary_box(1,3,five) = -1.e10
      boundary_box(2,3,five) = g_zmin
      boundary_box(1,3,six) = g_zmax
      boundary_box(2,3,six) = 1.e10
      endif
      
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!start test
! set the solution array to be the grid points x,y or z coordinates
        do l=1,lnblocks

      if(nodetype(l).eq.1 .or. advance_all_levels) then

       if (ndim.eq.3) dz = bsize(3,l)/real(nzb)
       if (ndim >= 2) dy = bsize(2,l)/real(nyb)
       dx = bsize(1,l)/real(nxb)
       do k=1+nguard0*k3d,nzb+nguard0*k3d
       do j=1+nguard0*k2d,nyb+nguard0*k2d
       do i=1+nguard0,nxb+nguard0
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard0)
       if (ndim >= 2) then
        y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard0)
       else
        y0 = 0.
       end if
       if(ndim == 3) then
        z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard0)
       else
        z0 = 0.
       end if
       value = ax*x0 + ay*y0 + az*z0
       do ivar=1,nvar
       unk(ivar,i,j,k,l) = value*real(ivar)
       enddo
       enddo
       enddo
       enddo

       if (nvar_work > 0) then
       do k=1+nguard_work0*k3d,nzb+nguard_work0*k3d
       do j=1+nguard_work0*k2d,nyb+nguard_work0*k2d
       do i=1+nguard_work0,nxb+nguard_work0
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard_work0)
       if (ndim >= 2) then
        y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard_work0)
       else
        y0 = 0.
       end if
       if(ndim == 3) then
        z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard_work0)
       else
        z0 = 0.
       end if
       value = ax*x0 + ay*y0 + az*z0
       work(i,j,k,l,ioptw-1) = value
       enddo
       enddo
       enddo
       end if


       if(nvarcorn.gt.0) then

       if (ndim.eq.3) dz = bsize(3,l)/real(nzb)
       if (ndim >= 2) dy = bsize(2,l)/real(nyb)
       dx = bsize(1,l)/real(nxb)
       do k=1+nguard0*k3d,nzb+nguard0*k3d+k3d
       do j=1+nguard0*k2d,nyb+nguard0*k2d+k2d
       do i=1+nguard0,nxb+nguard0+1
       x0 = coord(1,l)-.5*bsize(1,l)-dx+dx*real(i-nguard0)
       if (ndim >= 2) then
        y0 = coord(2,l)-.5*bsize(2,l)-dy+dy*real(j-nguard0)
       else
        y0 = 0.
       end if
       if(ndim == 3) then
        z0 = & 
     &       coord(3,l)-.5*bsize(3,l)-dz+dz*real(k-nguard0)
       else
        z0 = 0.
       end if
       value = ax*x0 + ay*y0 + az*z0
       do ivar=1,nvarcorn
       unk_n(ivar,i,j,k,l) = value*real(ivar)
       enddo
       enddo
       enddo
       enddo

       endif

      endif

      enddo
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up data in facevarx etc

      if(nfacevar.gt.0) then

        do l=1,lnblocks

        if(nodetype(l).eq.1 .or. advance_all_levels) then

              if(ndim == 3) dz = bsize(3,l)/real(nzb)
              if(ndim >= 2) dy = bsize(2,l)/real(nyb)
              dx = bsize(1,l)/real(nxb)
              if(mod(nxb,2).eq.1) then
                      if (ndim == 3) dz = bsize(3,l)/real(nzb-k3d)
                      if (ndim >= 2) dy = bsize(2,l)/real(nyb-k2d)
                      dx = bsize(1,l)/real(nxb-1)
              endif


              do k=1+nguard0*k3d,nzb+nguard0*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
                zk = z0 + dz*real(k-nguard0)
                  do j=1+nguard0*k2d,nyb+nguard0*k2d
                    if (ndim >= 2) y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                    yj = y0 + dy*real(j-nguard0)
                    do i=1+nguard0,nxb+nguard0+1
                      x0 = coord(1,l)-.5*bsize(1,l)-dx
                      xi = x0 + dx*real(i-nguard0)
                      do ivar=1,nbndvar
                        value = ax*xi+ay*yj+az*zk
                        facevarx(ivar,i,j,k,l)=value*real(ivar)
                      enddo
                    enddo
                  enddo
              enddo

              do k=1+nguard0*k3d,nzb+nguard0*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
                zk = z0 + dz*real(k-nguard0)
                do i=1+nguard0,nxb+nguard0
                  x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                  xi = x0 + dx*real(i-nguard0)
                  do j=1+nguard0*k2d,nyb+(nguard0+1)*k2d
                    if (ndim >= 2) y0 = coord(2,l)-.5*bsize(2,l)-dy
                    yj = y0 + dy*real(j-nguard0)
                    do ivar=1,nbndvar
                      value = ax*xi+ay*yj+az*zk
                      facevary(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo

              do j=1+nguard0*k2d,nyb+nguard0*k2d
                if (ndim >= 2) y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                yj = y0 + dy*real(j-nguard0)
                do i=1+nguard0,nxb+nguard0
                  x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                  xi = x0 + dx*real(i-nguard0)
                  do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d
                    if(ndim.eq.3) z0 =  & 
     &                       coord(3,l)-.5*bsize(3,l)-dz
                    zk = z0 + dz*real(k-nguard0)
                    do ivar=1,nbndvar
                      value = ax*xi+ay*yj+az*zk
                      facevarz(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo

        endif

      enddo

      endif

! set up data in unk_e_x[y][z] 

      if(nvaredge.gt.0) then

      do l=1,lnblocks

      if(nodetype(l).le.2 .or. advance_all_levels) then
 
              if (ndim.eq.3) dz = bsize(3,l)/real(nzb)
              if (ndim >= 2) dy = bsize(2,l)/real(nyb)
              dx = bsize(1,l)/real(nxb)
              if(mod(nxb,2).eq.1) then
                      if (ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                      if (ndim >= 2) dy = bsize(2,l)/real(nyb-k2d)
                      dx = bsize(1,l)/real(nxb-1)
              endif
 
              do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
                zk = z0 + dz*real(k-nguard0)
                  do j=1+nguard0*k2d,nyb+(nguard0+1)*k2d
                    if(ndim >= 2) y0 = coord(2,l)-.5*bsize(2,l)-dy
                    yj = y0 + dy*real(j-nguard0)
                    do i=1+nguard0,nxb+nguard0
                      x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                      xi = x0 + dx*real(i-nguard0)
                      do ivar=1,nvaredge
                        value = ax*xi+ay*yj+az*zk
                        unk_e_x(ivar,i,j,k,l)=value*real(ivar)
                      enddo
                    enddo
                  enddo
              enddo                                          

	if (ndim >= 2) then
              do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
                zk = z0 + dz*real(k-nguard0)
                do i=1+nguard0,nxb+nguard0+1
                  x0 = coord(1,l)-.5*bsize(1,l)-dx
                  xi = x0 + dx*real(i-nguard0)
                  do j=1+nguard0*k2d,nyb+nguard0*k2d
                    if (ndim >= 2) y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                    yj = y0 + dy*real(j-nguard0)
                    do ivar=1,nvaredge
                      value = ax*xi+ay*yj+az*zk
                      unk_e_y(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo
	end if
 
	if (ndim == 3) then
              do j=1+nguard0*k2d,nyb+(nguard0+1)*k2d
                if (ndim >= 2) y0 = coord(2,l)-.5*bsize(2,l)-dy
                yj = y0 + dy*real(j-nguard0)
                do i=1+nguard0,nxb+nguard0+1
                  x0 = coord(1,l)-.5*bsize(1,l)-dx
                  xi = x0 + dx*real(i-nguard0)
                  do k=1+nguard0*k3d,nzb+nguard0*k3d
                    if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
                    zk = z0 + dz*real(k-nguard0)
                    do ivar=1,nvaredge
                      value = ax*xi+ay*yj+az*zk
                      unk_e_z(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo
	end if	
 
        endif

      enddo

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loop_count=0
! Now cycle over blocks adjusting refinement of initial setup as required
        do while(loop_count.lt.3)


        refine(:) = .false.
      if(loop_count.lt.2) then
       refine(:) = .true.
      elseif(loop_count.eq.2) then
      if(ndim.eq.3) then
       do l=1,lnblocks
       if( coord(1,l).eq..125.and.coord(2,l).eq..125.and. & 
     &       coord(3,l).eq..125) refine(l)=.true.
       if( coord(1,l).eq..375.and.coord(2,l).eq..375.and. & 
     &       coord(3,l).eq..375) refine(l)=.true.
       if( coord(1,l).eq..625.and.coord(2,l).eq..875.and. & 
     &       coord(3,l).eq..875) refine(l)=.true.
       enddo
      elseif(ndim.eq.2) then
                do l=1,lnblocks
                if( coord(1,l).eq..125.and.coord(2,l).eq..125) & 
     &                 refine(l)=.true.
                if( coord(1,l).eq..375.and.coord(2,l).eq..375) & 
     &                 refine(l)=.true.
                if( coord(1,l).eq..625.and.coord(2,l).eq..875) & 
     &                 refine(l)=.true.
                enddo
      elseif(ndim==1) then
                do l=1,lnblocks
                if( coord(1,l).eq..125) & 
     &                 refine(l)=.true.
                if( coord(1,l).eq..375) & 
     &                 refine(l)=.true.
                if( coord(1,l).eq..625) & 
     &                 refine(l)=.true.
                enddo

      endif
      endif


! refine grid and apply morton reordering to grid blocks if necessary
      call amr_refine_derefine

      if (nvar > 0 .and. nvar_work > 0) then

      do lb = 1, lnblocks
      noff = (nguard_work - nguard)*npgs
      if(noff.ge.0) then
        work(il_bnd+noff:iu_bnd+noff,jl_bnd+noff*k2d:ju_bnd+noff*k2d, & 
     &       kl_bnd+noff*k3d:ku_bnd+noff*k3d,lb,ioptw-1) =  & 
     &   unk(1,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
      else
        work(ilw:iuw,jlw:juw,klw:kuw,lb,ioptw-1) = & 
     &   unk(1,il_bnd-noff:iu_bnd+noff, & 
     &         jl_bnd-noff*k2d:ju_bnd+noff*k2d, & 
     &         kl_bnd-noff*k3d:ku_bnd+noff*k3d,lb)
      endif
      end do

      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! a global prolongation call resets the newchild marker flags to false.
! Thus to test prolong for work we will need to restore this after the
! prolongation is applied to unk and facevar's
      tnewchild(:) = newchild(:)

      tag_offset = 100
      call mpi_morton_bnd_prolong & 
     &             (mype,nprocs,tag_offset)

#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited mpi_morton_bnd_prolong : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */


      iopt = 1
      nlayers = nguard
      call amr_prolong(mype,iopt,nlayers)
#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited amr_prolong : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */

#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'finished amr_prolong for unk : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'starting amr_prolong for work : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */


      if (nvar_work > 0) then
      newchild(:) = tnewchild(:)
      iopt = ioptw
      nlayers = nguard_work
      lcc = .true.
      lfc = .false.
      lec = .false.
      lnc = .false.
      call amr_prolong(mype,iopt,nlayers)
#ifdef DEBUG
      do lb = 1,lnblocks
       do j=jlw,juw
       write(*,*) 'work ',mype,lb,j,work(:,j,1,lb,iopt-1)
       enddo
      enddo

       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'finished amr_prolong for work : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */

       end if

#ifdef TEMP_TEST
!-----
       if (.not.advance_all_levels) then
       write(*,*) 'This part not tested '
       write(*,*) 'Should be added in other test programs '
       call mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc,.true.)
       end if
!-----
#endif /* TEMP_TEST */



      if (no_permanent_guardcells) then
! Store a copy of the current solution in gt_unk
      call amr_1blk_copy_soln(-1)
      end if

      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)

      iopt = 1


      if (.not.no_permanent_guardcells) then


! set guard cell data to zero to ensure proper test of guardcell
! set external guard cell data.
       call zero_guardcells(ioptw)


      iopt = 1
      nlayers = nguard
      lcc = .true.
      lfc = .true.
      lec = .true.
      lnc = .true. 
      tag_offset = 100
      write(*,*) 'pe ',mype ,' calling amr_guardcell iopt ',iopt
      call amr_guardcell(mype,iopt,nlayers)
      write(*,*) 'pe ',mype ,' exited amr_guardcell iopt ',iopt

      if (nvar_work > 0) then
      iopt = ioptw
      nlayers = nguard_work
      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(nvar.gt.0) lcc = .true.
      if(nfacevar.gt.0) lfc = .true.
      if(nvaredge.gt.0) lec = .true.
      if(nvarcorn.gt.0) lnc = .true.
      tag_offset = 100

      call amr_guardcell(mype,iopt,nlayers)
      end if

      else                      ! no_permanent_guardcells

      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(nvar.gt.0) lcc = .true.
      if(nfacevar.gt.0) lfc = .true.
      if(nvaredge.gt.0) lec = .true.
      if(nvarcorn.gt.0) lnc = .true.
      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs, & 
     &                        lguard,lprolong,lflux,ledge, & 
     &                        lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)

      end if                    ! no_permanent_guardcells

        loop_count=loop_count+1
!        write(*,*) 'proc loop_count ',mype,loop_count

        enddo
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        do ii=0,nprocs-1
                if(mype.eq.ii) then
                do l=1,lnblocks
                write(*,*) 'proc ',ii,' block ',l,' coord= ', & 
     &                  (coord(lcoord,l),lcoord=1,ndim), & 
     &                  ' size = ',bsize(1,l)
                write(*,*) 'proc ',ii,' block ',l,' parent= ', & 
     &                  parent(1,l),parent(2,l)
                write(*,*) 'proc ',ii,' block ',l,' nodety= ', & 
     &                  nodetype(l)
                enddo
                endif
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        enddo

      iopt = 1

       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

      if(mype.eq.0) write(*,*) 'Start of automatic testing.'

       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test of unk communications

      if (nvar > 0) then

      do ii=0,nprocs-1
      if(mype.eq.ii) then

      do l=1,lnblocks

      if (no_permanent_guardcells) then
      iopt = 1
      nlayers = nguard 
      lcc = .true.
      lfc = .false.
      lec = .false.
      lnc = .false.
      l_srl_only = .false.
!     icoord = 0
      ldiag = .false.
      ldiag = diagonals
      if(nodetype(l).eq.1 .or. advance_all_levels) then

      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)

      endif

      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
      unk1(:,:,:,:,1) = unk(:,:,:,:,l)
#endif
      end if                    ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or. advance_all_levels) then

      ilbnd=il_bnd1
      iubnd=iu_bnd1
      jlbnd=jl_bnd1
      jubnd=ju_bnd1
      klbnd=kl_bnd1
      kubnd=ku_bnd1

      if(neigh(1,1,l).le.-20) ilbnd=1+nguard
      if(neigh(1,2,l).le.-20) iubnd=nxb+nguard
      if (ndim >= 2) then
       if(neigh(1,3,l).le.-20) jlbnd=1+nguard
       if(neigh(1,4,l).le.-20) jubnd=nyb+nguard
      end if
      if(ndim.eq.3) then
       if(neigh(1,5,l).le.-20) klbnd=1+nguard
       if(neigh(1,6,l).le.-20) kubnd=nzb+nguard
      endif

      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+1)*k2d:nyb+2*ngd*k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          nyb+(ngd+1)*k2d:nyb+2*ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        endif
      endif 


      if (ndim.eq.3) dz = bsize(3,l)/real(nzb)
      if (ndim >= 2) dy = bsize(2,l)/real(nyb)
      dx = bsize(1,l)/real(nxb)
      if(mod(nxb,2).eq.1) then
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       if (ndim >= 2) dy = bsize(2,l)/real(nyb-k2d)
       dx = bsize(1,l)/real(nxb-1)
      endif

      do k=klbnd,kubnd
      if(ndim.eq.3) then
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       if(mod(nxb,2).eq.1) z0 = coord(3,l)-(.5*bsize(3,l)+dz)
      endif
      zk = z0 + dz*real(k-nguard)
       do j=jlbnd,jubnd
       if (ndim >= 2) then
        y0 = coord(2,l)-.5*(bsize(2,l)+dy)
        if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       end if
       yj = y0 + dy*real(j-nguard)
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)
       if(mod(nxb,2).eq.1) x0 =  & 
     &       coord(1,l)-(.5*bsize(1,l)+dx)
       xi = x0 + dx*real(i-nguard)
       value = ax*xi+ay*yj+az*zk
       do ivar=1,nvar
       v0 = value*real(ivar)
       if(abs(v0-unk1(ivar,i,j,k,1))>accuracy.and. & 
     &    lmask(i,j,k)) then
         write(*,998) ii,l,ivar,i,j,k,unk1(ivar,i,j,k,1),v0
         print *,' ay, yj = ',ay,yj
         ierror_sum = ierror_sum + 1
       endif
       enddo
       enddo
       enddo
      enddo

      endif

      enddo

      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo

      end if ! end if (nvar > 0) 


      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'unk test complete.'
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test of work communications

      if (nvar_work > 0) then

      iopt = ioptw

      lcc = .true.
      lfc = .false.
      lec = .false.
      lnc = .false. 
      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs, & 
     &                        lguard,lprolong,lflux,ledge, & 
     &                        lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)

      do ii=0,nprocs-1
      if(mype.eq.ii) then

      do l=1,lnblocks

      if (no_permanent_guardcells) then
      iopt = ioptw
      nlayers = nguard_work
      lcc = .true.
      lfc = .false.
      lec = .false.
      lnc = .false.
      l_srl_only = .false.
!     icoord = 0
      ldiag = .false.
      ldiag = diagonals

      if(nodetype(l).eq.1 .or. advance_all_levels) then

      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)

      endif

      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
      work1(:,:,:,1) = work(:,:,:,l,ioptw-1)
#endif
      end if                    ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or. advance_all_levels) then

      ilbnd=ilw1
      iubnd=iuw1
      jlbnd=jlw1
      jubnd=juw1
      klbnd=klw1
      kubnd=kuw1

      if(neigh(1,1,l).le.-20) ilbnd=1+nguard_work
      if(neigh(1,2,l).le.-20) iubnd=nxb+nguard_work
      if (ndim >= 2) then
       if(neigh(1,3,l).le.-20) jlbnd=1+nguard_work
       if(neigh(1,4,l).le.-20) jubnd=nyb+nguard_work
      end if
      if(ndim.eq.3) then 
       if(neigh(1,5,l).le.-20) klbnd=1+nguard_work
       if(neigh(1,6,l).le.-20) kubnd=nzb+nguard_work
      endif


      ngd = nguard_work
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+1)*k2d:nyb+2*ngd*k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          nyb+(ngd+1)*k2d:nyb+2*ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        endif
      endif 


      if( ndim.eq.3) dz = bsize(3,l)/real(nzb)
      if (ndim >= 2) dy = bsize(2,l)/real(nyb)
      dx = bsize(1,l)/real(nxb)
      if(mod(nxb,2).eq.1) then
       if( ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       if (ndim >= 2) dy = bsize(2,l)/real(nyb-k2d)
       dx = bsize(1,l)/real(nxb-1)
      endif

      do k=klbnd,kubnd
      if(ndim.eq.3) then
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       if(mod(nxb,2).eq.1) z0 = coord(3,l)-(.5*bsize(3,l)+dz)
      endif
      zk = z0 + dz*real(k-nguard_work)
       do j=jlbnd,jubnd
       if (ndim >= 2) then
        y0 = coord(2,l)-.5*(bsize(2,l)+dy)
        if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       end if
       yj = y0 + dy*real(j-nguard_work)
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)
       if(mod(nxb,2).eq.1) x0 =  & 
     &       coord(1,l)-(.5*bsize(1,l)+dx)
       xi = x0 + dx*real(i-nguard_work)
       value = ax*xi+ay*yj+az*zk
       if(abs(value-work1(i,j,k,1))>accuracy.and. & 
     &    lmask(i,j,k)) then
         write(*,997) ii,l,i,j,k,work1(i,j,k,1),value
         ierror_sum = ierror_sum + 1
       endif
       enddo
       enddo
      enddo

      endif

      enddo
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo

      end if ! end if (nvar_work > 0) 

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'work test complete.'
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if (nvarcorn > 0) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test of unk_n communications
      iopt = 1
      nlayers = nguard

      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .true.
      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs, & 
     &                        lguard,lprolong,lflux,ledge, & 
     &                        lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)
 
      do ii=0,nprocs-1
      if(mype.eq.ii) then
 
      do l=1,lnblocks

      if(nodetype(l).eq.1 .or. advance_all_levels) then
 
      if (no_permanent_guardcells) then
      iopt = 1
      nlayers = nguard
      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .true.
      l_srl_only = .false.
!     icoord = 0
      ldiag = .false.
      ldiag = diagonals
 
      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)
 
      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
      unk_n1(:,:,:,:,1) = unk_n(:,:,:,:,l)
#endif
      end if                    ! no_permanent_guardcells
 
 
        ilbnd=il_bnd1
        iubnd=iu_bnd1
        jlbnd=jl_bnd1
        jubnd=ju_bnd1
        klbnd=kl_bnd1
        kubnd=ku_bnd1
 
        if(neigh(1,1,l).le.-20) ilbnd=1+nguard
        if(neigh(1,2,l).le.-20) iubnd=nxb+nguard
        if (ndim >= 2) then
         if(neigh(1,3,l).le.-20) jlbnd=1+nguard
         if(neigh(1,4,l).le.-20) jubnd=nyb+nguard
        end if
        if(ndim.eq.3) then
          if(neigh(1,5,l).le.-20) klbnd=1+nguard
          if(neigh(1,6,l).le.-20) kubnd=nzb+nguard
        endif
 
      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd+1,1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+k3d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+2)*k3d:nzb+2*ngd*k3d+k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          nzb+(ngd+2)*k3d:nzb+2*ngd*k3d+k3d) = .true.
        endif
      endif 

 
        if (ndim.eq.3) dz = bsize(3,l)/real(nzb)
        if (ndim >= 2) dy = bsize(2,l)/real(nyb)
        dx = bsize(1,l)/real(nxb)
        if(mod(nxb,2).eq.1) then
                if (ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                if (ndim >= 2) dy = bsize(2,l)/real(nyb-1)
                dx = bsize(1,l)/real(nxb-1)
        endif
                                                                  
      do k=klbnd,kubnd+k3d
        if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
        zk = z0 + dz*real(k-nguard)
        do j=jlbnd,jubnd+k2d
          if (ndim >= 2) y0 = coord(2,l)-.5*bsize(2,l)-dy
          yj = y0 + dy*real(j-nguard)
          do i=ilbnd,iubnd+1
            x0 = coord(1,l)-.5*bsize(1,l)-dx
            xi = x0 + dx*real(i-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nvarcorn
              v0 = value*real(ivar)
              if(abs(unk_n1(ivar,i,j,k,1)-v0)>accuracy.and. & 
     &           lmask(i,j,k)) then
                      write(*,993) mype,l,ivar,i,j,k, & 
     &                             unk_n1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
                      endif
            enddo
          enddo
        enddo
      enddo
 
      endif

      enddo
 
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo
 
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'unk_n test complete.'
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end if ! end if (nvarcorn > 0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nfacevar > 0) then
 
#ifdef FACE
! test of guardcell for facevar's (if appropriate)
      if(mod(nxb,2).eq.0) then

      iopt = 1
      lcc = .false.
      lfc = .true.
      lec = .false.
      lnc = .false.
      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs, & 
     &                        lguard,lprolong,lflux,ledge, & 
     &                        lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)


! test values
        do ii=0,nprocs-1
        if(mype.eq.ii) then

        do l=1,lnblocks

      if (no_permanent_guardcells) then
      iopt = 1
      nlayers = nguard
      lcc = .false.
      lfc = .true.
      lec = .false.
      lnc = .false.
      l_srl_only = .false.
!     icoord = 0
      ldiag = .false.
      ldiag = diagonals

      if(nodetype(l).eq.1 .or. advance_all_levels) then

      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)

      endif

      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
      facevarx1(:,:,:,:,1) = facevarx(:,:,:,:,l)
      facevary1(:,:,:,:,1) = facevary(:,:,:,:,l)
      facevarz1(:,:,:,:,1) = facevarz(:,:,:,:,l)
#endif
      end if                    ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or. advance_all_levels) then

        ilbnd=il_bnd1
        iubnd=iu_bnd1
        jlbnd=jl_bnd1
        jubnd=ju_bnd1
        klbnd=kl_bnd1
        kubnd=ku_bnd1

        if(neigh(1,1,l).le.-20) ilbnd=1+nguard
        if(neigh(1,2,l).le.-20) iubnd=nxb+nguard
        if (ndim >= 2) then
         if(neigh(1,3,l).le.-20) jlbnd=1+nguard
         if(neigh(1,4,l).le.-20) jubnd=nyb+nguard
        end if
        if(ndim.eq.3) then
          if(neigh(1,5,l).le.-20) klbnd=1+nguard
          if(neigh(1,6,l).le.-20) kubnd=nzb+nguard
        endif

        ione = 1

        if (ndim.eq.3) dz = bsize(3,l)/real(nzb)
        if (ndim >= 2) dy = bsize(2,l)/real(nyb)
        dx = bsize(1,l)/real(nxb)
        if(mod(nxb,2).eq.1) then
                if (ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                if (ndim >= 2) dy = bsize(2,l)/real(nyb-k2d)
                dx = bsize(1,l)/real(nxb-1)
        endif

! first test facevarx
#ifdef FACEX

      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd+1,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+1)*k2d:nyb+2*ngd*k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          nyb+(ngd+1)*k2d:nyb+2*ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        endif
      endif 



      do k=klbnd,kubnd
        if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
        zk = z0 + dz*real(k-nguard)
        do j=jlbnd,jubnd
          if (ndim >= 2) y0 = coord(2,l)-.5*(bsize(2,l)+dy)
          yj = y0 + dy*real(j-nguard)
          do i=ilbnd,iubnd+ione
            x0 = coord(1,l)-.5*bsize(1,l)-dx
            xi = x0 + dx*real(i-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevarx1(ivar,i,j,k,1)-v0)>accuracy & 
     &                        .and.lmask(i,j,k)) then
                      write(*,996) mype,l,ivar,i,j,k, & 
     &                             facevarx1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
                      endif
            enddo
          enddo
        enddo
      enddo
#endif
! now test facevary
#ifdef FACEY

      if (ndim >= 2) then
      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd,1+ngd*k2d:nyb+ngd*k2d+1, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d+1, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d+1, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+1,:) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+1, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d+1, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d+1, & 
     &          nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        endif
      endif 


      do k=klbnd,kubnd
        if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
        zk = z0 + dz*real(k-nguard)
        do i=ilbnd,iubnd
          x0 = coord(1,l)-.5*(bsize(1,l)+dx)
          xi = x0 + dx*real(i-nguard)
          do j=jlbnd,jubnd+ione
            if (ndim >= 2) y0 = coord(2,l)-.5*bsize(2,l)-dy
            yj = y0 + dy*real(j-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevary1(ivar,i,j,k,1)-v0)>accuracy & 
     &                        .and.lmask(i,j,k)) then
                      write(*,995) mype,l,ivar,i,j,k, & 
     &                             facevary1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
              endif
            enddo
          enddo
        enddo
      enddo
      end if
#endif

! finally test facevarz
#ifdef FACEZ

      if (ndim == 3) then
      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+1) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+1) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+1) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+1)*k2d:nyb+2*ngd*k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+1) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          nyb+(ngd+1)*k2d:nyb+2*ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+1) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+2)*k3d:nzb+2*ngd*k3d+1) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          nzb+(ngd+1)*k3d:nzb+2*ngd*k3d+1) = .true.
        endif
      endif 


      do j=jlbnd,jubnd
        y0 = coord(2,l)-.5*(bsize(2,l)+dy)
        yj = y0 + dy*real(j-nguard)
        do i=ilbnd,iubnd
          x0 = coord(1,l)-.5*(bsize(1,l)+dx)
          xi = x0 + dx*real(i-nguard)
          do k=klbnd,kubnd+ione
            z0 = coord(3,l)-.5*bsize(3,l)-dz
            zk = z0 + dz*real(k-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevarz1(ivar,i,j,k,1)-v0)>accuracy & 
     &                        .and.lmask(i,j,k)) then
                     write(*,994) mype,l,ivar,i,j,k, & 
     &                             facevarz1(ivar,i,j,k,1),v0
                     ierror_sum = ierror_sum + 1
              endif
            enddo
          enddo
        enddo
      enddo
      end if
#endif


      endif

      enddo
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'face var. guard cell test complete.'
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      endif

#endif

      end if ! end if (nfacevar > 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nvaredge > 0) then

#ifdef EDGES
      if (ndim == 3) then
! test of guardcell for unk_e's (if appropriate)
      if(mod(nxb,2).eq.0) then
 
      iopt = 1
      lcc = .false.
      lfc = .false.
      lec = .true.
      lnc = .false.
      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs, & 
     &                        lguard,lprolong,lflux,ledge, & 
     &                        lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)
 
! test values
        do ii=0,nprocs-1
        if(mype.eq.ii) then
 
        do l=1,lnblocks
 

      if(nodetype(l).eq.1 .or. advance_all_levels) then

      if (no_permanent_guardcells) then
      iopt = 1
      nlayers = nguard
      lcc = .false.
      lfc = .false.
      lec = .true.
      lnc = .false.
      l_srl_only = .false.
!     icoord = 0
      ldiag = .false.
      ldiag = diagonals
      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)
 
      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
      if (nvaredge > 0) then
        unk_e_x1(:,:,:,:,1) = unk_e_x(:,:,:,:,l)
        unk_e_y1(:,:,:,:,1) = unk_e_y(:,:,:,:,l)
        unk_e_z1(:,:,:,:,1) = unk_e_z(:,:,:,:,l)
      end if
#endif
      end if                    ! no_permanent_guardcells

        ilbnd=il_bnd1
        iubnd=iu_bnd1
        jlbnd=jl_bnd1
        jubnd=ju_bnd1
        klbnd=kl_bnd1
        kubnd=ku_bnd1
 
        if(neigh(1,1,l).le.-20) ilbnd=1+nguard
        if(neigh(1,2,l).le.-20) iubnd=nxb+nguard
        if(neigh(1,3,l).le.-20) jlbnd=1+nguard
        if(neigh(1,4,l).le.-20) jubnd=nyb+nguard
        if(ndim.eq.3) then
          if(neigh(1,5,l).le.-20) klbnd=1+nguard
          if(neigh(1,6,l).le.-20) kubnd=nzb+nguard
        endif
 
        ione = 1
 
        if (ndim.eq.3) dz = bsize(3,l)/real(nzb)
        if (ndim >= 2) dy = bsize(2,l)/real(nyb)
        dx = bsize(1,l)/real(nxb)
        if(mod(nxb,2).eq.1) then
                if (ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                if (ndim >= 2) dy = bsize(2,l)/real(nyb-k2d)
                dx = bsize(1,l)/real(nxb-1)
        endif                                 

! first test unk_e_x
#ifdef UNKE_X
      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd,1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
          lmask(nxb+ngd+1:nxb+2*ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+2)*k3d:nzb+2*ngd*k3d+k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          nzb+(ngd+2)*k3d:nzb+2*ngd*k3d+k3d) = .true.
        endif
      endif 


      do k=klbnd,kubnd+ione
        if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
        zk = z0 + dz*real(k-nguard)
        do j=jlbnd,jubnd+ione
          y0 = coord(2,l)-.5*bsize(2,l)-dy
          yj = y0 + dy*real(j-nguard)
          do i=ilbnd,iubnd
            x0 = coord(1,l)-.5*(bsize(1,l)+dx)
            xi = x0 + dx*real(i-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nvaredge
              v0 = value*real(ivar)
              if(abs(unk_e_x1(ivar,i,j,k,1)-v0)>accuracy.and. & 
     &           lmask(i,j,k)) then
                      write(*,990) mype,l,ivar,i,j,k, & 
     &                             unk_e_x1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
                      endif
            enddo
          enddo
        enddo
      enddo
#endif
! now test unk_e_y
#ifdef UNKE_Y                     
      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd+1,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+1)*k2d:nyb+2*ngd*k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          nyb+(ngd+1)*k2d:nyb+2*ngd*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d+k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+2)*k3d:nzb+2*ngd*k3d+k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d, & 
     &          nzb+(ngd+2)*k3d:nzb+2*ngd*k3d+k3d) = .true.
        endif
      endif 


      do k=klbnd,kubnd+ione
        if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
        zk = z0 + dz*real(k-nguard)
        do i=ilbnd,iubnd+ione
          x0 = coord(1,l)-.5*bsize(1,l)-dx
          xi = x0 + dx*real(i-nguard)
          do j=jlbnd,jubnd
            y0 = coord(2,l)-.5*(bsize(2,l)+dy)
            yj = y0 + dy*real(j-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nvaredge
              v0 = value*real(ivar)
              if(abs(unk_e_y1(ivar,i,j,k,1)-v0)>accuracy.and. & 
     &           lmask(i,j,k)) then
                      write(*,991) mype,l,ivar,i,j,k, & 
     &                             unk_e_y1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
              endif
            enddo
          enddo
        enddo
      enddo
#endif
 
! finally test unk_e_z
#ifdef UNKE_Z                       
      ngd = nguard
      lmask(:,:,:) = .false.
      lmask(1+ngd:nxb+ngd+1,1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.

      if(icoord.eq.1) then
        if(ldiag) then
          lmask(1:ngd,:,:) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1,:,:) = .true.
        else
          lmask(1:ngd,1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &      1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(nxb+ngd+2:nxb+2*ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.2) then
        if(ldiag) then
          lmask(:,1:1+(ngd-1)*k2d,:) = .true.
          lmask(:,nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+k2d,:) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1:1+(ngd-1)*k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          nyb+(ngd+2)*k2d:nyb+2*ngd*k2d+k2d, & 
     &          1+ngd*k3d:nzb+ngd*k3d) = .true.
        endif
      elseif(icoord.eq.3) then
        if(ldiag) then
          lmask(:,:,1:1+(ngd-1)*k3d) = .true.
          lmask(:,:,nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        else
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          1:1+(ngd-1)*k3d) = .true.
          lmask(1+ngd:nxb+ngd+1, & 
     &          1+ngd*k2d:nyb+ngd*k2d+k2d, & 
     &          nzb+(ngd+1)*k3d:nzb+2*ngd*k3d) = .true.
        endif
      endif 


      do j=jlbnd,jubnd+ione
        y0 = coord(2,l)-.5*bsize(2,l)-dy
        yj = y0 + dy*real(j-nguard)
        do i=ilbnd,iubnd+ione
          x0 = coord(1,l)-.5*bsize(1,l)-dx
          xi = x0 + dx*real(i-nguard)
          do k=klbnd,kubnd
            z0 = coord(3,l)-.5*(bsize(3,l)+dz)
            zk = z0 + dz*real(k-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nvaredge
              v0 = value*real(ivar)
              if(abs(unk_e_z1(ivar,i,j,k,1)-v0)>accuracy.and. & 
     &           lmask(i,j,k)) then
                     write(*,992) mype,l,ivar,i,j,k, & 
     &                             unk_e_z1(ivar,i,j,k,1),v0
                     ierror_sum = ierror_sum + 1
              endif
            enddo
          enddo
        enddo
      enddo
#endif
 
      endif

      enddo
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo
 
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'unk_e. guard cell test complete.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)         
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      endif

      end if

      end if ! end if (nvaredge > 0.

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(mype.eq.0) write(*,*)  & 
     &   'End of automatic testing for icoord = ',icoord
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      enddo                                  ! end of icoord loop


      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call comm_int_sum_to_all(ierror_tot,ierror_sum)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      lmpi = .false.
      lnperm = .false.
      ladvanceall = .false.
          lmpi = .true.
          if (no_permanent_guardcells) then
          lnperm = .true.
          end if
          if (advance_all_levels) then
          ladvanceall = .true.
          end if
      if(mype.eq.0) then
      filename = trim(output_dir) // 'test.log'
      open(unit = 55,file=filename, & 
     &            status='unknown',position='append')
        write(*,*) ' '
        write(*,*) ' '
        if(ierror_tot.eq.0) then
          write(*,*) 'No errors detected - Test Successful '
          write(55,*) 'No errors detected - ', & 
     &                'Test of 1blk_guardcell_icoord Successful ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        else
          write(*,*) ierror_tot,' errors detected - Test failed '
          write(55,*) ierror_tot,' errors detected - ', & 
     &                'Test of blk_guardcell_icoord failed ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        endif
      close(unit=55)
      endif


      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'Start guardcell consistency check.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!      call guardcell_test(mype)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'Guardcell consistency check done.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'Start mesh check.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call mesh_test(mype)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'Mesh check done.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call report(mype,nprocs,ndim,no_permanent_guardcells, & 
     &            advance_all_levels,lmpi, & 
     &            ierror_tot,.true., & 
     &            '1blk_guardcell_icoord         ')

 2    call amr_close()

         
998      format('u:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
997      format('w:error proc ',i3,' block l= ',4(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
996      format('fx:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
995      format('fy:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
994      format('fz:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)

993      format('n:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
992      format('ez:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
991      format('ey:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
990      format('ex:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)

      deallocate(tnewchild)
      deallocate(lmask)

      end
