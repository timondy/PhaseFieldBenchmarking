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

!#define LOU_2D_PATTERN
!#define KEVINS_2D_PATTERN
#define STANDARD_TESTS
#define DEBUG

! This program refines a spherical grid with radial grid between 1.0 and 6.0.
! Once the grid has been refined, the solution is reset to a specific linear function
! of theta, which is discontinuous at theta=pi/2, and zero at theta=0 and pi.
! the function is
!         f  =  theta            if theta leq pi/4
!         f  =  pi/2 - theta     if pi/4 < theta < 3*pi/4
!         f  =  pi - theta       if theta geq 3*pi/4
!
! In this way we can use linear interpolation without conservation to automatically
! test the guardcell filling, including at refinement jumps, in the vicinity
! of the poles, without special fixes to the testing. The one exception is excluding
! automatic testing across the theta=pi/2 block boundaries, since this is where the
! solution is discontinuous. Because of this last issue, the refinement patterns used
! in this test should avoid refinement jumps across theta=pi/2.

! This program also checks the grid.
! It checks the grid by comparing it with a file which stores
! the correct grid.

! In order to create the correct grid files for later checking define 
! RECORD_GRID once you known everything is working, and execute the code once
! with nprocs=1. Then for subsequent runs undefine RECORD_GRID and
! then define CHECK_GRID.

!#define RECORD_GRID
!#define CHECK_GRID

      program test_singular_coord_sph_gcell




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use workspace
      use io
      use constants
      Use paramesh_comm_data

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
     &                                gtest_neigh_data, & 
     &                                amr_close

      use paramesh_mpi_interfaces, only :  & 
     &                                mpi_morton_bnd, & 
     &                                mpi_amr_comm_setup, & 
     &                                mpi_amr_1blk_restrict, & 
     &                                mpi_morton_bnd_prolong

! Only required for programs in ./Tests
#include "test_defs.fh"

      include 'mpif.h'

      integer :: tag_offset,max_blks_sent
      integer nguard0
      integer nguard_work0
      integer :: four,five, six
      real,parameter :: eps = 1.e-10
      real    :: accuracy
      real    :: polarfact
      logical :: lpolar
      real,external :: ftheta
      logical :: ltest,lfy,lfz
      logical,external :: test_switch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local amr variables
      integer nprocs,mype
      integer num_procs

      save mype

      logical, allocatable :: tnewchild(:)
      logical, allocatable :: rflags(:)

!
! application specific variables

      integer iopt,nlayers,icoord
      integer ierror_sum,ierror_tot
      logical lrefine_again,ltype2only
      logical lcc, lfc, lec, lnc, l_srl_only, ldiag

      logical :: lmpi,lnperm,ladvanceall
      integer :: errorcode, ierr, errcode

      logical :: lguard,lprolong,lflux,ledge,lrestrict
      logical :: lfulltree

      integer :: surrblks(2,3,3,3)
      save       surrblks

      save ierror_sum,ierror_tot

      integer :: error_grid
      integer :: lnblocks_tot
      integer :: gmaxblocks
      real,allocatable  :: gcoord(:,:)
      real    :: rcoord(3)
      integer,allocatable  :: glnblocks(:)
      real    :: test1,test2,test3

      character (len=80) filename
      character (len=6) filenumber

      integer :: maxrefine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      call amr_initialize 

       write(*,*)
       write(*,*) 'WARNING : Temporarily set N_EDGE_VAR = 0 to make ', & 
     &            'this test work past the mpi_morton_bnd_fluxcon call.'
       write(*,*) 'WARNING : Temporarily set N_EDGE_VAR = 0 to make ', & 
     &            'this test work past the mpi_morton_bnd_fluxcon call.'
       write(*,*)

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      error_grid = 0

      allocate(tnewchild(maxblocks_tr))
      allocate(rflags(maxblocks_tr))

      if(ndim.ne.3) then
         write(*,*) 'ndim=',ndim,' : Must set ndim=3 for this test!'
         go to 222
      endif
      if (.not.curvilinear) then
         write(*,*) 'CURVILINEAR must be defined for this test!'
         go to 222
      end if
      if (conserve) then
         write(*,*) 'CONSERVE must not be defined for this test!'
         go to 222
      end if

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

      rflags(:) = .true.

      ierror_sum = 0
      ierror_tot = 0

      iopt = 1
      nlayers = nguard
      if(mype.eq.0) write(*,*) 'nlayers = ',nlayers

!

! set a limit on the refinement level
      lrefine_max = 6
      lrefine_min = 1

      ax = 0.
      ay = 1.
      az = 0.
#ifdef TESTXDIR
      ax = 1.
      ay = 0.
      az = 0.
#endif
#ifdef TESTYDIR
      ax = 0.
      ay = 1.
      az = 0.
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

      if (spherical_pm) then
      g_xmin = 1.
      g_xmax = 6.
      g_ymin = 0.
      g_ymax = pi
      g_zmin = 0.
      g_zmax = 2.*pi
      end if

      if (polar_pm) then
      g_xmin = 0.
      g_xmax = 1.
      g_ymin = 0.
      g_ymax = 2.*pi
      g_zmin = 0.
      g_zmax = 1.
      end if

! For spherical coordinates which include poles, set up at least 2 blocks
! covering the whole grid, one for 0<phi<pi and the other for pi<phi<2pi.
! This allows us to establish the correct guardcell assignments
! at the poles in the theta direction.

      lnblocks = 0
      if(mype.eq.0.) then
                lnblocks = 1

      if (spherical_pm) then
                lnblocks = 2
             end if

                bnd_box(1,1,1:lnblocks) = g_xmin
                bnd_box(2,1,1:lnblocks) = g_xmax
                bnd_box(1,2,1:lnblocks) = g_ymin
                bnd_box(2,2,1:lnblocks) = g_ymax

      if (polar_pm) then
                bnd_box(1,3,1:lnblocks) = g_zmin
                bnd_box(2,3,1:lnblocks) = g_zmax
             end if

             if (spherical_pm) then
                bnd_box(1,3,1) = g_zmin
                bnd_box(2,3,1) = g_zmax/2.
                bnd_box(1,3,2) = g_zmax/2.
                bnd_box(2,3,2) = g_zmax
       write(*,*) 'gminmax ',g_xmin,g_xmax,g_ymin,g_ymax,g_zmin, & 
     &                    g_zmax
             end if

                coord(:,1:lnblocks) =  & 
     &             .5*(bnd_box(1,:,1:lnblocks)+bnd_box(2,:,1:lnblocks))
                bsize(:,1:lnblocks) =  & 
     &                (bnd_box(2,:,1:lnblocks)-bnd_box(1,:,1:lnblocks))
                nodetype(1:lnblocks) = 1
                lrefine(1:lnblocks) = 1

                if (spherical_pm) then
                neigh(:,1,1:2) = -60        ! radial boundary condition
                neigh(:,2,1:2) = -61        ! radial boundary condition
                neigh(1,3:4,1) = 2        ! elevation - needs pointer modification to connect neighbors at poles
                neigh(1,3:4,2) = 1        ! elevation - needs pointer modification to connect neighbors at poles
                neigh(2,3:4,1:2) = 0      ! elevation - needs pointer modification to connect neighbors at poles
                neigh(1,5:6,1) = 2        ! azimuthal periodicity
                neigh(2,5:6,1) = 0        ! azimuthal periodicity
                neigh(1,5:6,2) = 1        ! azimuthal periodicity
                neigh(2,5:6,2) = 0        ! azimuthal periodicity
                end if

                if (polar_pm) then
                 
                neigh(1,1,1) = -101       ! inner radial boundary condition - if the inner
                                          ! radius is at r=0 then this BC must acquire the correct
                                          ! data 
                neigh(2,1,1) = -101       ! inner radial boundary condition 

!                neigh(1,1,1) = 1         ! inner radial boundary condition
!                neigh(2,1,1) = 0         ! inner radial boundary condition - block 
!                                         !   should point to its own face but the 
!                                         !   surrounding blks and neigh flags will 
!                                         !   need to be treated with special care later
                neigh(1:2,2,1) = -51      ! outer radial boundary condition
                neigh(1,3:4,1) = 1        ! polar angle periodicity
                neigh(2,3:4,1) = 0        ! polar angle periodicity
                end if

                refine(1)=.true.
      endif

      if (spherical_pm) then
      boundary_index(1) = -60
      boundary_index(2) = -61
      boundary_index(3:6) = -1
! inner radial boundary
      boundary_box(1,1,1) = -1.e10
      boundary_box(1,2,1) = -1.e10   ! may be able to set to 0
      boundary_box(1,3,1) = -1.e10   ! may be able to set to 0
      boundary_box(2,1,1) =  g_xmin
      boundary_box(2,2,1) = 1.e10    ! may be able to set to pi
      boundary_box(2,3,1) = 1.e10    ! may be able to set to 2pi
! outer radial boundary
      boundary_box(1,1,2) = g_xmax
      boundary_box(1,2,2) = -1.e10   ! may be able to set to 0
      boundary_box(1,3,2) = -1.e10   ! may be able to set to 0
      boundary_box(2,1,2) = 1.e10
      boundary_box(2,2,2) = 1.e10    ! may be able to set to pi
      boundary_box(2,3,2) = 1.e10    ! may be able to set to pi
      end if

      if (polar_pm) then
      boundary_index(1) = -51
! outer radial boundary
      boundary_box(1,2:3,1) = -1.e10
      boundary_box(2,2:3,1) =  1.e10
      boundary_box(1,1,1) = g_xmax
      boundary_box(2,1,1) = 1.e10
      end if

#ifdef NOTNOW
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
      five = (4*k3d)+1
      six  = five + k3d
      boundary_box(1,1:2,five:six) = -1.e10
      boundary_box(2,1:2,five:six) =  1.e10
      boundary_box(1,3,five) = -1.e10
      boundary_box(2,3,five) = g_zmin
      boundary_box(1,3,six) = g_zmax
      boundary_box(2,3,six) = 1.e10
      endif
#endif /* NOTNOW */
      
      

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize solution to zero
      unk = 0.
      if (nfacevar > 0) then
         facevarx = 0.
         facevary = 0.
         facevarz = 0.
      end if
      if (nvaredge > 0) then
         unk_e_x = 0.
         unk_e_y = 0.
         unk_e_z = 0.
      end if
      unk_n = 0.

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do ii=0,nprocs-1
                if(mype.eq.ii) then
                write(*,*) '[',mype,'] lnblocks = ',lnblocks
                do l=1,lnblocks
                write(*,*) '[',mype,']  block ',l,' coord= ', & 
     &                  (coord(icoord,l),icoord=1,ndim), & 
     &                  ' size = ',bsize(1,l)
                write(*,*) 'proc ',ii,' block ',l,' parent= ', & 
     &                  parent(1,l),parent(2,l)
                write(*,*) 'proc ',ii,' block ',l,' nodety= ', & 
     &                  nodetype(l)
                write(*,*) 'proc ',ii,' block ',l,' neigh= ', & 
     &                  neigh(1,1:6,l)
!     .                  neigh(:,1:4,l)
                enddo
                endif
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        enddo

!      call amr_block_boundary_tecplot(10)
!      call output_tecplot_refmap(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loop_count=0
      loop_count_max=3
      maxrefine = 1
! Now cycle over blocks adjusting refinement of initial setup as required
        do while(loop_count.lt.loop_count_max)

        write(*,*) '[',mype,'] !!!!    loop_count start = ',loop_count
        write(*,*) '[',mype,'] !!!!    loop_count start = ',loop_count
        write(*,*) '[',mype,'] !!!!    loop_count start = ',loop_count
        maxrefine = maxrefine + 1

#ifdef STANDARD_TESTS
        refine(:) = .false.

        if(ndim.eq.3) then

      if (spherical_pm) then
          if(loop_count.eq.0) then
           refine(:) = .true.
          elseif(loop_count.eq.1) then
           refine(:) = .true.
          elseif(loop_count.eq.2) then
            do l=1,lnblocks
              test1 = abs(coord(1,l)-(1.+.375*5.))
              test2 = abs(coord(2,l)-pi*7./8.)
              test3 = abs(coord(3,l)-pi*3./8.)
              if( test1.lt.eps.and.test2.lt.eps.and.test3.lt.eps) & 
     &            refine(l)=.true.
              if(refine(l)) write(*,*) 'blk ',l,' refine ',refine(l), & 
     &          ' coord ',coord(:,l),1.+5.*.25,pi/4.,pi/4.
              test1 = abs(coord(1,l)-(1.+.625*5.))
              test2 = abs(coord(2,l)-pi*1./8.)
              test3 = abs(coord(3,l)-pi*1./8.)
              if( test1.lt.eps.and.test2.lt.eps.and.test3.lt.eps) & 
     &            refine(l)=.true.
              if(refine(l)) write(*,*) 'blk ',l,' refine ',refine(l), & 
     &          ' coord ',coord(:,l),1.+5.*.25,pi/4.,pi/4.
            enddo
          endif
          end if

          elseif(ndim.eq.2) then

      if (polar_pm) then
          write(*,*) 'USING POLAR COORDS'
          if(loop_count.lt.2) then
           refine(:) = .true.
          elseif(loop_count.eq.2) then
            do l=1,lnblocks
              test1 = abs(coord(1,l)-.125)
              test2 = abs(coord(2,l)-pi/4.)
              if( test1.lt.eps.and.test2.lt.eps) & 
     &            refine(l)=.true.
            enddo
          elseif(loop_count.eq.3) then
            do l=1,lnblocks
              test1 = abs(coord(1,l)-.0625)
              test2 = abs(coord(2,l)-pi/8.)
              if( test1.lt.eps.and.test2.lt.eps) & 
     &            refine(l)=.true.
            enddo
          endif
          endif 

        endif

#endif /*  STANDARD_TESTS */

! refine grid and apply morton reordering to grid blocks if necessary
#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'entering amr_refine_derefine : pe ',mype, & 
     &    ' loop_count = ',loop_count, & 
     &    ' refine(1:lnblocks) = ',refine(1:lnblocks)
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */

      call amr_refine_derefine
#ifdef DEBUG
!       call amr_flush(6)
!       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited amr_refine_derefine : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */

#ifdef RECORD_GRID
       if(nprocs.eq.1) then
          write(filenumber,"(I6.6)") loop_count
          open(unit = 56,file='grid_sph_correct_lc.'//filenumber, & 
     &            status='unknown')
          write(56,450) lnblocks
450       format(i5)
          do l=1,lnblocks
            write(56,451) l,coord(:,l)
451         format(i5,3(1x,1pe15.8))
          enddo
          close(unit=56)
       endif
#endif /* RECORD_GRID */
#ifdef CHECK_GRID
! First collect coord data and lnblocks onto pe 0.
        gmaxblocks = maxblocks*nprocs
        write(*,*) 'gmaxblocks ',gmaxblocks
        if(.not.allocated(gcoord)) allocate(gcoord(3,gmaxblocks))

! broadcast number of blocks on each processor
        if(.not.allocated(glnblocks)) allocate ( glnblocks(0:nprocs-1) )
        glnblocks(mype) = lnblocks
        call MPI_Allgather(glnblocks(mype), 1,MPI_INTEGER, & 
     &                   glnblocks,1,MPI_INTEGER, & 
     &                   MPI_COMM_WORLD,ierror)

        call comm_int_sum_to_all(lnblocks_tot,lnblocks)
        call mpi_gather(coord,3*maxblocks, & 
     &                     amr_mpi_real, & 
     &                     gcoord,3*maxblocks, & 
     &                     amr_mpi_real, & 
     &                     0,MPI_COMM_WORLD,ierror)

! Then compare with stored correct result
       if(mype.eq.0) then

#ifdef RECORD_GRIDX
          write(filenumber,"(I6.6)") loop_count
          open(unit = 56,file='grid_sph_correct_lc.'//filenumber, & 
     &            status='unknown')
          write(56,450) lnblocks_tot
450       format(i5)
          do i=0,nprocs-1
          ll0 = maxblocks*i
          do l=1,glnblocks(i)
            ll = ll0 + l
            write(56,451) ll,gcoord(:,ll)
451         format(i5,3(1x,1pe15.8))
          enddo
          enddo

          close(unit=56)
#endif /* RECORD_GRIDX */

          write(filenumber,"(I6.6)") loop_count
          open(unit = 56,file='grid_sph_correct_lc.'//filenumber, & 
     &            Status='old')
          read(56,460) rlnblocks_tot
          if(rlnblocks_tot.ne.lnblocks_tot) then
            
          endif
460       format(i5)
          do i=0,nprocs-1
          ll0 = maxblocks*i
          do l=1,glnblocks(i)
            ll = ll0 + l
            read(56,461) lr,rcoord
461         format(i5,3(1x,1pe15.8))
            test1 = abs((gcoord(1,ll)-rcoord(1))/rcoord(1))
            test2 = abs((gcoord(2,ll)-rcoord(2))/pi)
            test3 = abs((gcoord(3,ll)-rcoord(3))/pi)
            if(test1.gt.1.e-8.or.test2.gt.1.e-8.or.test3.gt.1.e-8) then
              write(*,*) 'grid error blk ',l,' ll ',ll,' pe ',i, & 
     &                   ' gcoord ',gcoord(:,ll),' rcoord ',rcoord, & 
     &                    test1,test2,test3,' loop_count ',loop_count
              error_grid = error_grid + 1
            endif
          enddo
          enddo
          close(unit=56)
       endif
#endif /* CHECK_GRID */

!      call amr_block_boundary_tecplot(loop_count)
!      call output_tecplot_refmap(loop_count)

      if (nvar > 0 .and. nvar_work > 0) then

! Copy unk into work so we can test prolong for iopt > 1
      noff = (nguard_work - nguard)*npgs
      if(noff.ge.0) then
        work(il_bnd+noff:iu_bnd+noff,jl_bnd+noff*k2d:ju_bnd+noff*k2d, & 
     &       kl_bnd+noff*k3d:ku_bnd+noff*k3d,:,ioptw-1) =  & 
     &   unk(1,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,:)
      else
        work(ilw:iuw,jlw:juw,klw:kuw,:,ioptw-1) =  & 
     &   unk(1,il_bnd-noff:iu_bnd+noff, & 
     &         jl_bnd-noff*k2d:ju_bnd+noff*k2d, & 
     &         kl_bnd-noff*k3d:ku_bnd+noff*k3d,:)
      endif

      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! a global prolongation call resets the newchild marker flags to false.
! Thus to test prolong for work we will need to restore this after the
! prolongation is applied to unk and facevar's
      tnewchild(:) = newchild(:)

      tag_offset = 100
      call mpi_morton_bnd_prolong & 
     &             (mype,nprocs,tag_offset)

!#ifdef DEBUG
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited mpi_morton_bnd_prolong : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
!#endif /* DEBUG */

      iopt = 1
      nlayers = nguard
       write(*,*) 'calling amr_prolong : pe ',mype, & 
     &            ' loop_count ',loop_count & 
     &           ,' lnblocks ',lnblocks
       if(mype.eq.2.and.lnblocks.eq.80) write(*,*) & 
     &      'MAIN: blk 23-2 coord ',coord(:,23)
      call amr_prolong(mype,iopt,nlayers)
       write(*,*) 'exited amr_prolong : pe ',mype, & 
     &            ' loop_count ',loop_count
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
!      call amr_block_boundary_tecplot(loop_count)
!      call output_tecplot_refmap(loop_count)

      end if

!#ifdef DEBUG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if(loop_count.eq.loop_count_max-1) then

! set the solution array to a linear function of the grid points x,y or z coordinates
        do l=1,lnblocks

      if(nodetype(l).eq.1 .or. advance_all_levels) then

       if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
       dy = bsize(2,l)/real(nyb)
       dx = bsize(1,l)/real(nxb)
       do k=1+nguard0*k3d,nzb+nguard0*k3d
       do j=1+nguard0,nyb+nguard0
       do i=1+nguard0,nxb+nguard0
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard0)
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard0)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard0)
       value = ftheta(y0,z0)
       do ivar=1,nvar
       unk(ivar,i,j,k,l) = value*real(ivar)
       enddo
       enddo
       enddo
       enddo

       if (nvar_work > 0) then
       do k=1+nguard_work0*k3d,nzb+nguard_work0*k3d
       do j=1+nguard_work0,nyb+nguard_work0
       do i=1+nguard_work0,nxb+nguard_work0
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard_work0)
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard_work0)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard_work0)
       value = ftheta(y0,z0)
       work(i,j,k,l,ioptw-1) = value
       enddo
       enddo
       enddo
       end if


       if(nvarcorn.gt.0) then

       if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
       dy = bsize(2,l)/real(nyb)
       dx = bsize(1,l)/real(nxb)
       do k=1+nguard0*k3d,nzb+nguard0*k3d+k3d
       do j=1+nguard0,nyb+nguard0+k2d
       do i=1+nguard0,nxb+nguard0+1
       x0 = coord(1,l)-.5*bsize(1,l)-dx+dx*real(i-nguard0)
       y0 = coord(2,l)-.5*bsize(2,l)-dy+dy*real(j-nguard0)
       if(ndim.eq.3) z0 = & 
     &       coord(3,l)-.5*bsize(3,l)-dz+dz*real(k-nguard0)
       value = ftheta(y0,z0)
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

              if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
              dy = bsize(2,l)/real(nyb)
              dx = bsize(1,l)/real(nxb)
              if(mod(nxb,2).eq.1) then
                      if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                      dy = bsize(2,l)/real(nyb-1)
                      dx = bsize(1,l)/real(nxb-1)
              endif


              do k=1+nguard0*k3d,nzb+nguard0*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
                zk = z0 + dz*real(k-nguard0)
                  do j=1+nguard0,nyb+nguard0
                    y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                    yj = y0 + dy*real(j-nguard0)
                    do i=1+nguard0,nxb+nguard0+1
                      x0 = coord(1,l)-.5*bsize(1,l)-dx
                      xi = x0 + dx*real(i-nguard0)
                      do ivar=1,nbndvar
                        value = ftheta(yj,zk)
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
                  do j=1+nguard0,nyb+nguard0+1
                    y0 = coord(2,l)-.5*bsize(2,l)-dy
                    yj = y0 + dy*real(j-nguard0)
                    do ivar=1,nbndvar
                      value = ftheta(yj,zk)
                      facevary(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo

              do j=1+nguard0,nyb+nguard0
                y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                yj = y0 + dy*real(j-nguard0)
                do i=1+nguard0,nxb+nguard0
                  x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                  xi = x0 + dx*real(i-nguard0)
                  do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d
                    if(ndim.eq.3) z0 =  & 
     &                       coord(3,l)-.5*bsize(3,l)-dz
                    zk = z0 + dz*real(k-nguard0)
                    do ivar=1,nbndvar
                      value = ftheta(yj,zk)
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
 
              if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
              dy = bsize(2,l)/real(nyb)
              dx = bsize(1,l)/real(nxb)
              if(mod(nxb,2).eq.1) then
                      if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                      dy = bsize(2,l)/real(nyb-1)
                      dx = bsize(1,l)/real(nxb-1)
              endif
 
              do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
                zk = z0 + dz*real(k-nguard0)
                  do j=1+nguard0,nyb+nguard0+1
                    y0 = coord(2,l)-.5*bsize(2,l)-dy
                    yj = y0 + dy*real(j-nguard0)
                    do i=1+nguard0,nxb+nguard0
                      x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                      xi = x0 + dx*real(i-nguard0)
                      do ivar=1,nvaredge
                        value = ftheta(yj,zk)
                        unk_e_x(ivar,i,j,k,l)=value*real(ivar)
                      enddo
                    enddo
                  enddo
              enddo                                          

              do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
                zk = z0 + dz*real(k-nguard0)
                do i=1+nguard0,nxb+nguard0+1
                  x0 = coord(1,l)-.5*bsize(1,l)-dx
                  xi = x0 + dx*real(i-nguard0)
                  do j=1+nguard0,nyb+nguard0
                    y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                    yj = y0 + dy*real(j-nguard0)
                    do ivar=1,nvaredge
                      value = ftheta(yj,zk)
                      unk_e_y(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo
 
              do j=1+nguard0,nyb+nguard0+1
                y0 = coord(2,l)-.5*bsize(2,l)-dy
                yj = y0 + dy*real(j-nguard0)
                do i=1+nguard0,nxb+nguard0+1
                  x0 = coord(1,l)-.5*bsize(1,l)-dx
                  xi = x0 + dx*real(i-nguard0)
                  do k=1+nguard0*k3d,nzb+nguard0*k3d
                    if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
                    zk = z0 + dz*real(k-nguard0)
                    do ivar=1,nvaredge
                      value = ftheta(yj,zk)
                      unk_e_z(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo
 
        endif

      enddo

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)



       endif                      ! loop_count_max end if
!------------------------------------------------


       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'finished amr_prolong for work : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
!#endif /* DEBUG */

      if (no_permanent_guardcells) then
! Store a copy of the current solution in gt_unk
      call amr_1blk_copy_soln(-1)
      end if


      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)

      iopt = 1


      if (.not.no_permanent_guardcells) then


! set guard cell data to zero to ensure proper test of guardcell
       call zero_guardcells(ioptw)


      iopt = 1
      nlayers = nguard
      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(nvar.gt.0) lcc = .true.
      if(nfacevar.gt.0) lfc = .true.
      if(nvaredge.gt.0) lec = .true.
      if(nvarcorn.gt.0) lnc = .true.
      tag_offset = 100
      write(*,*) 'calling amr_guardcell iopt=',iopt
      call amr_guardcell(mype,iopt,nlayers)

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
      write(*,*) 'calling amr_guardcell iopt=',iopt
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
        write(*,*) 'proc loop_count ',mype,loop_count

!        enddo
!        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        do ii=0,nprocs-1
                if(mype.eq.ii) then
                do l=1,lnblocks
                write(*,*) 'proc ',ii,' block ',l,' coord= ', & 
     &                  (coord(icoord,l),icoord=1,ndim), & 
     &                  ' size = ',bsize(1,l)
                write(*,*) 'proc ',ii,' block ',l,' parent= ', & 
     &                  parent(1,l),parent(2,l)
                write(*,*) 'proc ',ii,' block ',l,' nodety= ', & 
     &                  nodetype(l)
                write(*,*) 'proc ',ii,' block ',l,' neigh= ', & 
     &                  neigh(1,1:6,l)
                write(*,*) 'proc ',ii,' block ',l,' surr_blks= ', & 
     &                  surr_blks(1,:,:,:,l)
!                write(*,*) 'proc ',ii,' block ',l,' neigh= ',
!     .                  neigh(:,1:4,l)
                enddo
                endif
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        enddo



        enddo
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)



      iopt = 1

       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

      if (nvar > 0 .and. nvar_work > 0) then
! temporary testing section
      noff = (nguard_work - nguard)*npgs
      if(noff.ge.0) then
        work(il_bnd+noff:iu_bnd+noff,jl_bnd+noff*k2d:ju_bnd+noff*k2d, & 
     &       kl_bnd+noff*k3d:ku_bnd+noff*k3d,:,ioptw-1) = & 
     &   unk(1,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,:)
      else
        work(ilw:iuw,jlw:juw,klw:kuw,:,ioptw-1) =  & 
     &   unk(1,il_bnd-noff:iu_bnd+noff, & 
     &         jl_bnd-noff*k2d:ju_bnd+noff*k2d, & 
     &         kl_bnd-noff*k3d:ku_bnd+noff*k3d,:)
      endif
! end temporary testing section
      end if



      if(mype.eq.0) write(*,*) 'Start of automatic testing.'


!
! make sure interp_mask is set for linear interpolation
      do ivar=1,nvar
        if(interp_mask_unk(ivar).ne.1) then
          write(*,*) 'test_singular_coord: ERROR. interp_mask_unk ', & 
     &               'is not set to 1 ie linear interpolation.'
          call amr_abort()
        endif
      enddo


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
      icoord = 0
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


      if(nodetype(l).eq.1 .or.  & 
     &     (advance_all_levels.and.lrefine(l).gt.1)) then
!      if(nodetype(l).eq.1) then

      write(*,*) 'block ',l,' bnd_box ',bnd_box(:,:,l), & 
     &            ' lrefine ',lrefine(l),' parent ',parent(:,l)

      ilbnd=il_bnd1
      iubnd=iu_bnd1
      jlbnd=jl_bnd1
      jubnd=ju_bnd1
      klbnd=kl_bnd1
      kubnd=ku_bnd1

#ifdef NOTNOW
      if(neigh(1,1,l).le.-20) ilbnd=1+nguard
      if(neigh(1,2,l).le.-20) iubnd=nxb+nguard
      if(neigh(1,3,l).le.-20) jlbnd=1+nguard
      if(neigh(1,4,l).le.-20) jubnd=nyb+nguard
      if(ndim.eq.3) then
      if(neigh(1,5,l).le.-20) klbnd=1+nguard
      if(neigh(1,6,l).le.-20) kubnd=nzb+nguard
      endif
#endif /* NOTNOW */

! do not test guardcells at the theta=pi/2 interior block interface
!      if(abs(bnd_box(2,2,l)-.5*pi).lt.eps) jubnd = nyb+nguard
!      if(abs(bnd_box(1,2,l)-.5*pi).lt.eps) jlbnd = 1+nguard
! do not test guardcells at the phi=0 or 2*pi interior block interfaces either
!      if(abs(bnd_box(2,3,l)-2.*pi).lt.eps) kubnd = nzb+nguard
!      if(abs(bnd_box(1,3,l)).lt.eps) klbnd = 1+nguard

      if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
      dy = bsize(2,l)/real(nyb)
      dx = bsize(1,l)/real(nxb)
      if(mod(nxb,2).eq.1) then
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       dy = bsize(2,l)/real(nyb-1)
       dx = bsize(1,l)/real(nxb-1)
      endif

      do k=klbnd,kubnd
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       zk = z0 + dz*real(k-nguard)
      do j=jlbnd,jubnd
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)
       if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       yj = y0 + dy*real(j-nguard)
       do i=ilbnd,iubnd
       value = ftheta(yj,zk)
       lfy = .false.
       lfz = .false.
       ltest = test_switch(neigh(:,:,l),bnd_box(:,:,l),i,j,k,iopt, & 
     &                    lfy,lfz)
       if(ltest) then
       do ivar=1,nvar
       v0 = value*real(ivar)
       if (abs(v0-unk1(ivar,i,j,k,1)) > accuracy) then
         write(*,998) ii,l,i,j,k,unk1(ivar,i,j,k,1),v0
         ierror_sum = ierror_sum + 1
       endif
       enddo
       endif
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
      icoord = 0
      ldiag = .false.
      ldiag = diagonals

      if(nodetype(l).eq.1 .or.  & 
     &     (advance_all_levels.and.lrefine(l).gt.1)) then

      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)

      endif

      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
         work1(:,:,:,1) = work(:,:,:,l,ioptw-1)
#endif
      end if                    ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or.  & 
     &     (advance_all_levels.and.lrefine(l).gt.1)) then

      ilbnd=ilw1
      iubnd=iuw1
      jlbnd=jlw1
      jubnd=juw1
      klbnd=klw1
      kubnd=kuw1

#ifdef NOTNOW
      if(neigh(1,1,l).le.-20) ilbnd=1+nguard_work
      if(neigh(1,2,l).le.-20) iubnd=nxb+nguard_work
      if(neigh(1,3,l).le.-20) jlbnd=1+nguard_work
      if(neigh(1,4,l).le.-20) jubnd=nyb+nguard_work
      if(ndim.eq.3) then
      if(neigh(1,5,l).le.-20) klbnd=1+nguard_work
      if(neigh(1,6,l).le.-20) kubnd=nzb+nguard_work
      endif
#endif

! do not test guardcells at the theta=pi/2 interior block interface
!      if(abs(bnd_box(2,2,l)-.5*pi).lt.eps) jubnd = nyb+nguard_work
!      if(abs(bnd_box(1,2,l)-.5*pi).lt.eps) jlbnd = 1+nguard_work
! do not test guardcells at the phi=0 or 2*pi interior block interfaces either
!      if(abs(bnd_box(2,3,l)-2.*pi).lt.eps) kubnd = nzb+nguard_work
!      if(abs(bnd_box(1,3,l)).lt.eps) klbnd = 1+nguard_work

      if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
      dy = bsize(2,l)/real(nyb)
      dx = bsize(1,l)/real(nxb)
      if(mod(nxb,2).eq.1) then
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       dy = bsize(2,l)/real(nyb-1)
       dx = bsize(1,l)/real(nxb-1)
      endif

      do k=klbnd,kubnd
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       zk = z0 + dz*real(k-nguard_work)
      do j=jlbnd,jubnd
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)
       if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       yj = y0 + dy*real(j-nguard_work)
       do i=ilbnd,iubnd
       value = ftheta(yj,zk)
       lfy = .false.
       lfz = .false.
       ltest = test_switch(neigh(:,:,l),bnd_box(:,:,l),i,j,k,iopt, & 
     &                    lfy,lfz)
       if(ltest) then
       if(abs(work1(i,j,k,1)-value) > accuracy) then
         write(*,997) ii,l,i,j,k,work1(i,j,k,1),value
         ierror_sum = ierror_sum + 1
       endif
       endif
       enddo
       enddo
      enddo

      endif

      enddo
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo

      end if  ! end if (nvar_work > 0)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'work test complete.'
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#ifdef NOTNOW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! test of unk_n communications

      if (nvarcorn > 0) then

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

      if(nodetype(l).eq.1 .or.  & 
     &     (advance_all_levels.and.lrefine(l).gt.1)) then
      if(refine(l).eq.maxrefine) then
 
      if (no_permanent_guardcells) then
      iopt = 1
      nlayers = nguard
      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .true.
      l_srl_only = .false.
      icoord = 0
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
        if(neigh(1,3,l).le.-20) jlbnd=1+nguard
        if(neigh(1,4,l).le.-20) jubnd=nyb+nguard
        if(ndim.eq.3) then
          if(neigh(1,5,l).le.-20) klbnd=1+nguard
          if(neigh(1,6,l).le.-20) kubnd=nzb+nguard
        endif
 
 
        if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
        dy = bsize(2,l)/real(nyb)
        dx = bsize(1,l)/real(nxb)
        if(mod(nxb,2).eq.1) then
                if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                dy = bsize(2,l)/real(nyb-1)
                dx = bsize(1,l)/real(nxb-1)
        endif
                                                                  
      do k=klbnd,kubnd+k3d
        if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
        zk = z0 + dz*real(k-nguard)
        do j=jlbnd,jubnd+k2d
          y0 = coord(2,l)-.5*bsize(2,l)-dy
          yj = y0 + dy*real(j-nguard)
          do i=ilbnd,iubnd+1
            x0 = coord(1,l)-.5*bsize(1,l)-dx
            xi = x0 + dx*real(i-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nvarcorn
              v0 = value*real(ivar)
              if(abs(unk_n1(ivar,i,j,k,1)-v0) > accuracy) then
                      write(*,993) mype,l,ivar,i,j,k, & 
     &                             unk_n1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
                      endif
            enddo
          enddo
        enddo
      enddo
 
      endif
      endif

      enddo
 
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo
 
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'unk_n test complete.'
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      end if

#endif /* NOTNOW */
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef FACE
! test of guardcell for facevar's (if appropriate)

      if (nfacevar > 0) then

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
      icoord = 0
      ldiag = .false.
      ldiag = diagonals

      if(nodetype(l).eq.1 .or.  & 
     &     (advance_all_levels.and.lrefine(l).gt.1)) then

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


      if(nodetype(l).eq.1 .or.  & 
     &     (advance_all_levels.and.lrefine(l).gt.1)) then

        ilbnd=il_bnd1
        iubnd=iu_bnd1
        jlbnd=jl_bnd1
        jubnd=ju_bnd1
        klbnd=kl_bnd1
        kubnd=ku_bnd1

#ifdef NOTNOW
        if(neigh(1,1,l).le.-20) ilbnd=1+nguard
        if(neigh(1,2,l).le.-20) iubnd=nxb+nguard
        if(neigh(1,3,l).le.-20) jlbnd=1+nguard
        if(neigh(1,4,l).le.-20) jubnd=nyb+nguard
        if(ndim.eq.3) then
          if(neigh(1,5,l).le.-20) klbnd=1+nguard
          if(neigh(1,6,l).le.-20) kubnd=nzb+nguard
        endif
#endif

! do not test guardcells at the theta=pi/2 interior block interface
!      if(abs(bnd_box(2,2,l)-.5*pi).lt.eps) jubnd = nyb+nguard
!      if(abs(bnd_box(1,2,l)-.5*pi).lt.eps) jlbnd = 1+nguard
! do not test guardcells at the phi=0 or 2*pi interior block interfaces either
!      if(abs(bnd_box(2,3,l)-2.*pi).lt.eps) kubnd = nzb+nguard
!      if(abs(bnd_box(1,3,l)).lt.eps) klbnd = 1+nguard

        ione = 1

        if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
        dy = bsize(2,l)/real(nyb)
        dx = bsize(1,l)/real(nxb)
        if(mod(nxb,2).eq.1) then
                if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                dy = bsize(2,l)/real(nyb-1)
                dx = bsize(1,l)/real(nxb-1)
        endif

! first test facevarx
#ifdef FACEX
      do k=klbnd,kubnd
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       zk = z0 + dz*real(k-nguard)
      do j=jlbnd,jubnd
          y0 = coord(2,l)-.5*(bsize(2,l)+dy)
          yj = y0 + dy*real(j-nguard)
          do i=ilbnd,iubnd+ione
          lfy = .false.
          lfz = .false.
          ltest = test_switch(neigh(:,:,l),bnd_box(:,:,l),i,j,k,iopt, & 
     &                        lfy,lfz)
          if(ltest) then
            value = ftheta(yj,zk)
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevarx1(ivar,i,j,k,1)-v0)>accuracy) then
                      write(*,996) mype,l,ivar,i,j,k, & 
     &                             facevarx1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
                      endif
            enddo
          endif
          enddo
        enddo
      enddo
#endif
! now test facevary
#ifdef FACEY

      do k=klbnd,kubnd
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       zk = z0 + dz*real(k-nguard)
        do i=ilbnd,iubnd
          do j=jlbnd,jubnd+ione
            y0 = coord(2,l)-.5*bsize(2,l)-dy
            yj = y0 + dy*real(j-nguard)
            value = ftheta(yj,zk)
            lfy = .true.
            lfz = .false.
            ltest = test_switch(neigh(:,:,l),bnd_box(:,:,l),i,j,k,iopt, & 
     &                          lfy,lfz)
          if(ltest) then
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevary1(ivar,i,j,k,1)-v0)>accuracy) then
                      write(*,995) mype,l,ivar,i,j,k, & 
     &                             facevary1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
              endif
            enddo
          endif
          enddo
        enddo
      enddo
#endif

! finally test facevarz
#ifdef FACEZ
! do not test guardcells at the phi=0, pi or 2*pi interior block interfaces either
      if(abs(bnd_box(1,3,l)).lt.eps) klbnd = 2+nguard
      if(abs(bnd_box(2,3,l)-pi).lt.eps) kubnd = nzb+nguard-1
      if(abs(bnd_box(1,3,l)-pi).lt.eps) klbnd = 2+nguard
      if(abs(bnd_box(2,3,l)-2.*pi).lt.eps) kubnd = nzb+nguard-1

      if (ndim == 3) then
      do j=jlbnd,jubnd
        y0 = coord(2,l)-.5*(bsize(2,l)+dy)
        yj = y0 + dy*real(j-nguard)
        do i=ilbnd,iubnd
          do k=klbnd,kubnd+ione
            z0 = coord(3,l)-.5*bsize(3,l)-dz
            zk = z0 + dz*real(k-nguard)
            lfy = .false.
            lfz = .true.
            ltest = test_switch(neigh(:,:,l),bnd_box(:,:,l),i,j,k,iopt, & 
     &                          lfy,lfz)
            if(ltest) then
            value = ftheta(yj,zk)
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevarz1(ivar,i,j,k,1)-v0)>accuracy) then
                     write(*,994) mype,l,ivar,i,j,k, & 
     &                             facevarz1(ivar,i,j,k,1),v0
                     ierror_sum = ierror_sum + 1
              endif
            enddo
            endif
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

      end if

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef NOTNOW

#ifdef EDGES
      
      if (nvaredge > 0) then

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
 
      if(nodetype(l).eq.1 .or.  & 
     &     (advance_all_levels.and.lrefine(l).gt.1)) then
      if(refine(l).eq.maxrefine) then

      if (no_permanent_guardcells) then
      iopt = 1
      nlayers = nguard
      lcc = .false.
      lfc = .false.
      lec = .true.
      lnc = .false.
      l_srl_only = .false.
      icoord = 0
      ldiag = .false.
      ldiag = diagonals
      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)
 
      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
         unk_e_x1(:,:,:,:,1) = unk_e_x(:,:,:,:,l)
         unk_e_y1(:,:,:,:,1) = unk_e_y(:,:,:,:,l)
         unk_e_z1(:,:,:,:,1) = unk_e_z(:,:,:,:,l)
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
 
        if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
        dy = bsize(2,l)/real(nyb)
        dx = bsize(1,l)/real(nxb)
        if(mod(nxb,2).eq.1) then
                if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
                dy = bsize(2,l)/real(nyb-1)
                dx = bsize(1,l)/real(nxb-1)
        endif                                 

! first test unk_e_x
#ifdef UNKE_X
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
              if(abs(unk_e_x1(ivar,i,j,k,1)-v0)>accuracy) then
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
              if(abs(unk_e_y1(ivar,i,j,k,1)-v0)>accuracy) then
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
              if(abs(unk_e_z1(ivar,i,j,k,1)-v0)>accuracy) then
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
      endif

      enddo
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo
 
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'unk_e. guard cell test complete.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)         

      endif
      endif

      end if

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif /* NOTNOW */

      if(mype.eq.0) write(*,*) 'End of automatic testing.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! Test consistency of face variables
      if(nfacevar.gt.0) then

      if(mype.eq.0) then
         write(*,*) ' '
         write(*,*) ' '
         write(*,*) 'Testing routine for div B consistency check'
         write(*,*) 'ax = ',ax
         write(*,*) 'ay = ',ay
         write(*,*) 'az = ',az
         if(ndim.eq.2) then
           write(*,*) 'div B should be = ax + ay = ',ax+ay
         elseif(ndim.eq.3) then
           write(*,*) 'div B should be = ax + ay + az = ',ax+ay+az
         endif
      endif
      test_a = ax + ay + az*k3d      

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      istep = 0
        write(*,*) 'calling gtest'
! commented out gtest_neigh_data call until facevar data works in spherical coords
!      call gtest_neigh_data(mype,istep,test_a)
        write(*,*) 'exited gtest'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      endif                       ! end of nfacevar iftest


! Sum number of errors 
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        write(*,*) 'ca isum'
      call comm_int_sum_to_all(ierror_tot,ierror_sum)
        write(*,*) 'ex isum'

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
 222  if(mype.eq.0) then
      filename = trim(output_dir) // 'test.log'
      open(unit = 55,file=filename, & 
     &            status='unknown',position='append')
        write(*,*) ' '
        write(*,*) ' '
        if(ierror_tot.eq.0.and.error_grid.eq.0) then
          write(*,*) 'No errors detected - Test Successful '
          write(55,*) 'No errors detected - ', & 
     &                'Test of singular_coord_sph_gcell Successful ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        else
          write(*,*) ierror_tot,' errors detected - Test failed '
          write(55,*) ierror_tot,' guardcell errors detected - ', & 
     &                error_grid,' grid errors detected - ', & 
     &                'Test of singular_coord_sph_gcell failed ', & 
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

      write(*,*) 'proc ',mype,' final lnblocks ',lnblocks

      call report(mype,nprocs,ndim,no_permanent_guardcells, & 
     &            advance_all_levels,lmpi, & 
     &            ierror_tot,.true., & 
     &            'singular_coord_sph_gcell      ')

      call amr_close()

         
998      format('u:error proc ',i3,' block l= ',4(2x,i3),2x,f7.4,2x, & 
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
!#ifdef CHECK_GRID
        if(allocated(gcoord)) deallocate(gcoord)
        if(allocated(glnblocks)) deallocate(glnblocks)
!#endif /* CHECK_GRID */


      deallocate(tnewchild)
      deallocate(rflags)

       write(*,*)
       write(*,*) 'WARNING : Temporarily set N_EDGE_VAR = 0 to make ', & 
     &            'this test work past the mpi_morton_bnd_fluxcon call.'
       write(*,*) 'WARNING : Temporarily set N_EDGE_VAR = 0 to make ', & 
     &            'this test work past the mpi_morton_bnd_fluxcon call.'
       write(*,*)


      end program test_singular_coord_sph_gcell


      function test_switch(neighl,bbox,i,j,k,iopt, & 
     &                     lfy,lfz)
      use paramesh_dimensions
      use physicaldata
      use constants
      implicit none
      real :: bbox(2,3)
      integer :: neighl(2,6),iopt,ng0
      integer :: jlbnd,jubnd,klbnd,kubnd,i,j,k
      integer :: joff,koff
      logical :: test_switch,lfy,lfz
      real,parameter :: eps = 1.e-10

      ng0 = nguard
      if(iopt.gt.1) ng0 = nguard_work

      joff = 0
      koff = 0
      if(lfy) joff = 1
      if(lfz) koff = 1

      jlbnd = 1
      jubnd = nyb+ng0*2 + joff
      klbnd = 1
      kubnd = nzb+ng0*2 + koff

! do not test guardcells at the theta=pi/2 interior block interface
      if(abs(bbox(2,2)-.5*pi).lt.eps) jubnd = nyb+ng0 + joff
      if(abs(bbox(1,2)-.5*pi).lt.eps) jlbnd = 1+ng0
! do not test guardcells at the phi=0 or 2*pi interior block interfaces either
      if(abs(bbox(2,3)-2.*pi).lt.eps) kubnd = nzb+ng0 + koff
      if(abs(bbox(1,3)).lt.eps) klbnd = 1+ng0

      test_switch = .true.
      test_switch = .false.
      if( (j.lt.jlbnd.or.j.gt.jubnd)  .or. & 
     &    (k.lt.klbnd.or.k.gt.kubnd)  ) then
        test_switch = .false.
      endif

! also dont test at some cells on just inside edges or corners next to the
! phi = 0 and 2*pi faces.
      if( (abs(bbox(1,3)).lt.eps) .and. & 
     &    (neighl(1,1).eq.-1) .and.  & 
     &    (k.eq.1+ng0).and.(i.le.ng0) ) test_switch = .false.
      if( (abs(bbox(1,3)).lt.eps) .and. & 
     &    (neighl(1,2).eq.-1) .and.  & 
     &    (k.eq.1+ng0).and.(i.gt.nxb+ng0) ) test_switch = .false.
      if( (abs(bbox(2,3)-2.*pi).lt.eps) .and. & 
     &    (neighl(1,1).eq.-1) .and.  & 
     &    (k.eq.nzb+ng0).and.(i.le.ng0) ) test_switch = .false.
      if( (abs(bbox(2,3)-2.*pi).lt.eps) .and. & 
     &    (neighl(1,2).eq.-1) .and.  & 
     &    (k.eq.nzb+ng0).and.(i.gt.nxb+ng0)) test_switch = .false.

      end function test_switch
