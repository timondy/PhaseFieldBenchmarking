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

#define NWORDS (NX_B + 2*N_GUARD_CELLS)

!-------------------------------------------------------------------
!
! This test is designed to check the routine amr_1blk_bcset_example
! which is located in ./templates.
! Before testing, ./templates/amr_1blk_bcset_example.F should be
! moved into the Tests directory.
! The test works by setting up a single block on pe 0, and
! filling its guardcells. The results are compared with a pre-existing
! file 
!    test_bcset_example.answer
! which contains the correct answer.
!
!-------------------------------------------------------------------

      program test_bcset_example


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use workspace

      use paramesh_interfaces, only : comm_start, & 
     &                                amr_initialize, & 
     &                                amr_refine_derefine, & 
     &                                amr_1blk_copy_soln, & 
     &                                amr_guardcell, & 
     &                                amr_prolong, & 
     &                                amr_restrict, & 
     &                                amr_1blk_guardcell, & 
     &                                amr_1blk_guardcell_reset, & 
     &                                amr_1blk_restrict, & 
     &                                amr_flux_conserve, & 
     &                                amr_edge_average, & 
     &                                guardcell_test, & 
     &                                mesh_test, & 
     &                                amr_close, & 
     &                                amr_surrounding_blks

      use paramesh_mpi_interfaces, only :  & 
     &                                mpi_morton_bnd, & 
     &                                mpi_amr_comm_setup, & 
     &                                mpi_amr_tree_setup, & 
     &                                mpi_amr_1blk_restrict, & 
     &                                mpi_morton_bnd_prolong, & 
     &                                mpi_amr_local_surr_blks


! Only required for programs in ./Tests
#include "test_defs.fh"

      include 'mpif.h'

      integer :: tag_offset,max_blks_sent
      integer nguard0
      integer nguard_work0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local amr variables
      integer nprocs,mype
      integer num_procs

      save mype

      logical :: tnewchild(maxblocks_tr)
      logical :: lfullblock
      logical :: lerror_setup
!
! application specific variables

      integer remote_block,remote_pe
      integer iopt,nlayers,icoord
      integer ierror_sum,ierror_tot
      integer cnodetype
      logical lrefine_again,ltype2only
      logical rflags(maxblocks_tr)
      logical lcc, lfc, lec, lnc, l_srl_only, ldiag

      logical :: lmpi,lnperm,ladvanceall
      integer :: errorcode

      logical :: lguard,lprolong,lflux,ledge,lrestrict
      logical :: lfulltree

      integer :: surrblks(2,3,3,3)
      integer :: psurrblks(2,3,3,3)
      logical :: l_parent
      integer :: ierrorcode,ierr

      save ierror_sum,ierror_tot
      save cnodetype
      save remote_block,remote_pe

      real :: unk1_a(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                    kl_bnd1:ku_bnd1)
      real :: unk_n1_a(nvarcorn,il_bnd1:iu_bnd1+1, & 
     &                 jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d)
      real :: unk_e_x1_a(nvaredge,il_bnd1:iu_bnd1, & 
     &                 jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d)
      real :: unk_e_y1_a(nvaredge,il_bnd1:iu_bnd1+1, & 
     &                 jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1+k3d)
      real :: unk_e_z1_a(nvaredge,il_bnd1:iu_bnd1+1, & 
     &                 jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1)
      real :: facevarx1_a(nfacevar,il_bnd1:iu_bnd1+1, & 
     &                 jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
      real :: facevary1_a(nfacevar,il_bnd1:iu_bnd1, & 
     &                 jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1)
      real :: facevarz1_a(nfacevar,il_bnd1:iu_bnd1, & 
     &                 jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1+k3d)

      character(len=5) :: char5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call amr_initialize

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

#ifdef CONSERVE
      write(*,*) 'CONSERVE must not be defined for this test!'
      call amr_abort()
#endif

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if(nprocs.gt.1) then
        write(*,*) & 
     &       'Error : this test is designed for 1 processor only'
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        stop
      endif

      ierror_sum = 0
      ierror_tot = 0

      if(mype.eq.0) then
      print *,' nprocs = ',nprocs,mype


! set default value of dz and z0 to cater for 2D case.
      z0 = 0.
      dz = 0.


      iopt = 1
      nlayers = nguard
      if(mype.eq.0) write(*,*) 'nlayers = ',nlayers

!

! set a limit on the refinement level
      lrefine_max = 5

      ax = 1.
      ay = 10.
      az = 100.
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

      if(mype.eq.0) write(*,*) 'ax = ',ax,' ay = ',ay,' az = ',az

! set the workspace array layer to be tested
      ioptw = 3


      if(ioptw-1.gt.nvar_work) then
          write(*,*) 'ERROR: Too few work arrays'
          write(*,*) 'ERROR: Reset nvar_work.'
          call amr_abort
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up initial grid state.

      g_xmin = 0.
      g_xmax = 1.
      g_ymin = 0.
      g_ymax = 1.
      g_zmin = 0.
      g_zmax = 1.


!---------------------------------


! reflecting BCs
      ibc = -22
! zero gradient BCs
      ibc = -21


!---------------------------------

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
                neigh(:,:,1) = ibc
                refine(1)=.true.
      endif


      boundary_index = ibc
! x boundaries
      boundary_box(1,2:3,1:2) = -1.e10
      boundary_box(2,2:3,1:2) =  1.e10
      boundary_box(1,1,1) = -1.e10
      boundary_box(2,1,1) = g_xmin
      boundary_box(1,1,2) = g_xmax
      boundary_box(2,1,2) = 1.e10
! y boundaries
      if(ndim.ge.2) then
      boundary_box(1,1,3:4) = -1.e10
      boundary_box(2,1,3:4) =  1.e10
      boundary_box(1,3,3:4) = -1.e10
      boundary_box(2,3,3:4) =  1.e10
      boundary_box(1,2,3) = -1.e10
      boundary_box(2,2,3) = g_ymin
      boundary_box(1,2,4) = g_ymax
      boundary_box(2,2,4) = 1.e10
      endif
! z boundaries
      if(ndim.eq.3) then
      boundary_box(1,1:2,5:6) = -1.e10
      boundary_box(2,1:2,5:6) =  1.e10
      boundary_box(1,3,5) = -1.e10
      boundary_box(2,3,5) = g_zmin
      boundary_box(1,3,6) = g_zmax
      boundary_box(2,3,6) = 1.e10
      endif
      
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!start test
! set the solution array to be the grid points x,y or z coordinates
        do l=1,lnblocks

      if(nodetype(l).eq.1 .or. advance_all_levels) then

       if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
       dy = bsize(2,l)/real(nyb)
       dx = bsize(1,l)/real(nxb)

       do k=1,nzb+2*nguard0*k3d
       do j=1,nyb+2*nguard0
       do i=1,nxb+2*nguard0
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard0)
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard0)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard0)
       value = ax*x0 + ay*y0 + az*z0
       do ivar=1,nvar
       unk(ivar,i,j,k,l) = value*real(ivar)
       enddo
       enddo
       enddo
       enddo

       do k=1,nzb+2*nguard_work0*k3d
       do j=1,nyb+2*nguard_work0
       do i=1,nxb+2*nguard_work0
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)+dx*real(i-nguard_work0)
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)+dy*real(j-nguard_work0)
       if(ndim.eq.3) z0 =  & 
     &       coord(3,l)-.5*(bsize(3,l)+dz)+dz*real(k-nguard_work0)
       value = ax*x0 + ay*y0 + az*z0
       work(i,j,k,l,ioptw-1) = value
       enddo
       enddo
       enddo

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
       value = ax*x0 + ay*y0 + az*z0
       do ivar=1,nvarcorn
       unk_n(ivar,i,j,k,l) = value*real(ivar)
       enddo
       enddo
       enddo
       enddo           


      endif

      enddo
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up data in facevarx etc


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


              do k=1,nzb+2*nguard0*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
                zk = z0 + dz*real(k-nguard0)
                  do j=1,nyb+2*nguard0
                    y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                    yj = y0 + dy*real(j-nguard0)
                    do i=1,nxb+2*nguard0+1
                      x0 = coord(1,l)-.5*bsize(1,l)-dx
                      xi = x0 + dx*real(i-nguard0)
                      do ivar=1,nbndvar
                        value = ax*xi+ay*yj+az*zk
                        facevarx(ivar,i,j,k,l)=value*real(ivar)
                      enddo
                    enddo
                  enddo
              enddo

              do k=1,nzb+2*nguard0*k3d
                if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
                zk = z0 + dz*real(k-nguard0)
                do i=1,nxb+2*nguard0
                  x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                  xi = x0 + dx*real(i-nguard0)
                  do j=1,nyb+(2*nguard0+1)*k2d
                    y0 = coord(2,l)-.5*bsize(2,l)-dy
                    yj = y0 + dy*real(j-nguard0)
                    do ivar=1,nbndvar
                      value = ax*xi+ay*yj+az*zk
                      facevary(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo

              do j=1,nyb+2*nguard0
                y0 = coord(2,l)-.5*(bsize(2,l)+dy)
                yj = y0 + dy*real(j-nguard0)
                do i=1,nxb+2*nguard0
                  x0 = coord(1,l)-.5*(bsize(1,l)+dx)
                  xi = x0 + dx*real(i-nguard0)
                  do k=1,nzb+(2*nguard0+1)*k3d
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

! set up data in unk_e_x[y][z]
 
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
                        value = ax*xi+ay*yj+az*zk
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
                      value = ax*xi+ay*yj+az*zk
                      unk_e_y(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo

              if (ndim == 3) then
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
                      value = ax*xi+ay*yj+az*zk
                      unk_e_z(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo
              else
              unk_e_z(:,:,:,:,l) = 0.
              end if

        endif

      enddo

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Store a copy of the current solution in gt_unk
      if (no_permanent_guardcells) then
      call amr_1blk_copy_soln(-1)
      end if


      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)

      iopt = 1


      if (.not.no_permanent_guardcells) then

      write(*,*) 'permanent guardcells'


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
      write(*,*) 'permanent guardcells: calling amr_guardcell'
      call amr_guardcell(mype,iopt,nlayers)
      write(*,*) 'permanent guardcells: exited amr_guardcell'

      iopt = ioptw
      nlayers = nguard_work
      lcc = .true.
      lfc = .false.
      lec = .false.
      lnc = .false.
      tag_offset = 100
      call amr_guardcell(mype,iopt,nlayers)

      else                      ! no_permanent_guardcells

      write(*,*) 'no permanent guardcells'
      lcc = .true.
      lfc = .true.
      lec = .true.
      lnc = .true.
      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      write(*,*) 'permanent guardcells: calling mpi_amr_comm_setup'
      call mpi_amr_comm_setup(mype,nprocs, & 
     &                        lguard,lprolong,lflux,ledge, & 
     &                        lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)
      end if                    ! no_permanent_guardcells


        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

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
                enddo
                endif
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      iopt = 1



      if (no_permanent_guardcells) then
! Store a copy of the current solution in gt_unk
      call amr_1blk_copy_soln(-1)
      end if

      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)

      iopt = 1

      if (no_permanent_guardcells) then
      lcc = .true.
      lfc = .true.
      lec = .true.
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
      end if                    ! no_permanent_guardcells



      if(mype.eq.0) write(*,*) 'Start of automatic testing.'


      if(mype.eq.0) open(unit=10,file='bcset_examples.dat', & 
     &                   status='unknown')
      if(mype.eq.0) then
        write(10,160) nvar
        write(10,160) nfacevar
        write(10,160) nvaredge
        write(10,160) nvarcorn
        write(10,160) nxb
        write(10,160) nyb
        write(10,160) nzb
        write(10,160) ndim
        write(10,160) nguard
        write(10,160) nprocs
        write(10,160) lnblocks
        write(10,160) ibc
        write(10,161) ax
        write(10,161) ay
        write(10,161) az
160     format(i3)
161     format(f7.5)
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set block number
      l = 1


      if (no_permanent_guardcells) then
      iopt = 1
      nlayers = nguard 
      lcc = .true.
      lfc = .true.
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
      unk1(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1) =  & 
     & unk(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,l)
      unk_n1(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,1) =  & 
     & unk_n(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,l)
      unk_e_x1(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,1) =  & 
     & unk_e_x(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,l)
      unk_e_y1(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,1) =  & 
     & unk_e_y(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,l)
      unk_e_z1(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,1) =  & 
     & unk_e_z(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,l)
      facevarx1(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1) =  & 
     & facevarx(:,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd,l)
      facevary1(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,1) =  & 
     & facevary(:,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,l)
      facevarz1(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,1) =  & 
     & facevarz(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,l)
      end if                    ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or. advance_all_levels) then

      if(nvar.gt.0) then

      do ivar=1,nvar
      write(10,150) ivar
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1
        write(10,151) j,k,unk1(ivar,:,j,k,1)
      enddo
      enddo
      enddo

      endif                   ! end of nvar if test


150   format('unk1 ',i3)
#if NWORDS == 6
151   format(i3,2x,i3,6(1x,f9.5))
153   format(i3,2x,i3,7(1x,f9.5))
#endif
#if NWORDS == 8
151   format(i3,2x,i3,8(1x,f9.5))
153   format(i3,2x,i3,9(1x,f9.5))
#endif
#if NWORDS == 10
151   format(i3,2x,i3,10(1x,f9.5))
153   format(i3,2x,i3,11(1x,f9.5))
#endif

      if(nvarcorn.gt.0) then

      do ivar=1,nvarcorn
      write(10,152) ivar
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1+k2d
        write(10,153) j,k,unk_n1(ivar,:,j,k,1)
      enddo
      enddo
      enddo
152   format('unk_n1 ',i3)

      endif                   ! end of nvarcorn if test

      if(nfacevar.gt.0) then


      do ivar=1,nfacevar
      write(10,154) ivar
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1
        write(10,153) j,k,facevarx1(ivar,:,j,k,1)
      enddo
      enddo
      enddo
154   format('facevarx1 ',i3)

      do ivar=1,nfacevar
      write(10,155) ivar
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1+k2d
        write(10,151) j,k,facevary1(ivar,:,j,k,1)
      enddo
      enddo
      enddo
155   format('facevary1 ',i3)

      if (ndim == 3) then
      do ivar=1,nfacevar
      write(10,156) ivar
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1
        write(10,151) j,k,facevarz1(ivar,:,j,k,1)
      enddo
      enddo
      enddo
156   format('facevarz1 ',i3)
      end if

      endif                   ! end of nfacevar if test

      if (ndim > 1) then
      if(nvaredge.gt.0) then

      do ivar=1,nvaredge
      write(10,157) ivar
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1+k2d
        write(10,151) j,k,unk_e_x1(ivar,:,j,k,1)
      enddo
      enddo
      enddo
157   format('unk_e_x1 ',i3)

      do ivar=1,nvaredge
      write(10,158) ivar
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1
        write(10,153) j,k,unk_e_y1(ivar,:,j,k,1)
      enddo
      enddo
      enddo
158   format('unk_e_y1 ',i3)

      if (ndim == 3) then
      do ivar=1,nvaredge
      write(10,159) ivar
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1+k2d
        write(10,153) j,k,unk_e_z1(ivar,:,j,k,1)
      enddo
      enddo
      enddo
159   format('unk_e_z1 ',i3)
      end if
      endif                   ! end of nvaredge if test
      end if


      endif

      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      close(unit=10)


      if(mype.eq.0) open(unit=11,file='test_bcset_examples.answer', & 
     &                   status='old')
      lerror_setup = .false.
      if(mype.eq.0) then
        read(11,160) nvar_a
        if(nvar_a.ne.nvar) lerror_setup = .true.
        read(11,160) nfacevar_a
        if(nfacevar_a.ne.nfacevar) lerror_setup = .true.
        read(11,160) nvaredge_a
        if(nvaredge_a.ne.nvaredge) lerror_setup = .true.
        read(11,160) nvarcorn_a
        if(nvarcorn_a.ne.nvarcorn) lerror_setup = .true.
        read(11,160) nxb_a
        if(nxb_a.ne.nxb) lerror_setup = .true.
        read(11,160) nyb_a
        if(nyb_a.ne.nyb) lerror_setup = .true.
        read(11,160) nzb_a
        if(nzb_a.ne.nzb) lerror_setup = .true.
        read(11,160) ndim_a
        if(ndim_a.ne.ndim) lerror_setup = .true.
        read(11,160) nguard_a
        if(nguard_a.ne.nguard) lerror_setup = .true.
        read(11,160) nprocs_a
        if(nprocs_a.ne.nprocs) lerror_setup = .true.
        read(11,160) lnblocks_a
        if(lnblocks_a.ne.lnblocks) lerror_setup = .true.
        read(11,160) ibc_a
        if(ibc_a.ne.ibc) lerror_setup = .true.
        read(11,161) ax_a
        if(ax_a.ne.ax) lerror_setup = .true.
        read(11,161) ay_a
        if(ay_a.ne.ay) lerror_setup = .true.
        read(11,161) az_a
        if(az_a.ne.az) lerror_setup = .true.
      endif
      ierror_tot = 1
      if(lerror_setup) ierror_tot = 1

      if(.not.lerror_setup) then

      lb = 1

      if(nvar.gt.0) then

      do ivar=1,nvar
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1
        read(11,151) j0,k0,unk1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1
         if(unk1_a(ivar,i,j,k).ne.unk1(ivar,i,j,k,1)) then
           write(55,998) mype,lb,ivar,i,j,k,unk1_a(ivar,i,j,k), & 
     &                   unk1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo
170   format(a5)

      endif                   ! end of nvar if test

      if(nvarcorn.gt.0) then

      do ivar=1,nvarcorn
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1+k2d
        read(11,153) j0,k0,unk_n1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1+1
         if(unk_n1_a(ivar,i,j,k).ne.unk_n1(ivar,i,j,k,1)) then
           write(55,993) mype,lb,ivar,i,j,k,unk_n1_a(ivar,i,j,k), & 
     &                   unk_n1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo

      endif                   ! end of nvarcorn if test

      if(nfacevar.gt.0) then

      do ivar=1,nfacevar
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1
        read(11,153) j0,k0,facevarx1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1+1
         if(facevarx1_a(ivar,i,j,k).ne. & 
     &                  facevarx1(ivar,i,j,k,1)) then
           write(55,896) mype,lb,ivar,i,j,k,facevarx1_a(ivar,i,j,k), & 
     &                   facevarx1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo

      do ivar=1,nfacevar
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1+k2d
        read(11,151) j0,k0,facevary1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1
         if(facevary1_a(ivar,i,j,k).ne. & 
     &                  facevary1(ivar,i,j,k,1)) then
           write(55,895) mype,lb,ivar,i,j,k,facevary1_a(ivar,i,j,k), & 
     &                   facevary1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo

      if (ndim == 3) then
      do ivar=1,nfacevar
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1
        read(11,151) j0,k0,facevarz1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1
         if(facevarz1_a(ivar,i,j,k).ne. & 
     &                  facevarz1(ivar,i,j,k,1)) then
           write(55,894) mype,lb,ivar,i,j,k,facevarz1_a(ivar,i,j,k), & 
     &                   facevarz1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo
      end if

      endif                   ! end of nfacevar if test



      if (ndim > 1) then
      if(nvaredge.gt.0) then

      do ivar=1,nvaredge
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1+k2d
        read(11,151) j0,k0,unk_e_x1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1
         if(unk_e_x1_a(ivar,i,j,k).ne.unk_e_x1(ivar,i,j,k,1)) then
           write(55,990) mype,lb,ivar,i,j,k,unk_e_x1_a(ivar,i,j,k), & 
     &                   unk_e_x1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo

      do ivar=1,nvaredge
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1+k3d
      do j = jl_bnd1,ju_bnd1
        read(11,153) j0,k0,unk_e_y1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1+1
         if(unk_e_y1_a(ivar,i,j,k).ne.unk_e_y1(ivar,i,j,k,1)) then
           write(55,991) mype,lb,ivar,i,j,k,unk_e_y1_a(ivar,i,j,k), & 
     &                   unk_e_y1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo

      if (ndim == 3) then
      do ivar=1,nvaredge
      read(11,170) char5
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1+k2d
        read(11,153) j0,k0,unk_e_z1_a(ivar,:,j,k)
        do i = il_bnd1,iu_bnd1+1
         if(unk_e_z1_a(ivar,i,j,k).ne.unk_e_z1(ivar,i,j,k,1)) then
           write(55,992) mype,lb,ivar,i,j,k,unk_e_z1_a(ivar,i,j,k), & 
     &                   unk_e_z1(ivar,i,j,k,1)
           ierror_tot = ierror_tot + 1
         endif
        enddo
      enddo
      enddo
      enddo
      end if
      endif                   ! end of nvaredge if test
      end if

      endif                    ! lerror_setup if test

      close(unit=11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(mype.eq.0) write(*,*) 'End of automatic testing.'
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
      open(unit = 55,file='test.log', & 
     &            status='unknown',position='append')
        write(*,*) ' '
        write(*,*) ' '
        if(ierror_tot.eq.0) then
          write(*,*) 'No errors detected - Test Successful '
          write(55,*) 'No errors detected - ', & 
     &                'Test of bcset_example Successful ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        else
          write(*,*) ierror_tot,' errors detected - Test failed '
          if(lerror_setup) then
          write(55,*) ierror_tot,' errors detected - ', & 
     &        'Test of bcset_example failed ', & 
     &        'problem setup is not consistent with regression ', & 
     &        'data file - modify setup or use different ', & 
     &        'regression data file.'
          else
          write(55,*) ierror_tot,' errors detected - ', & 
     &                'Test of bcset_example failed ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
          endif
        endif
      close(unit=55)
      endif


      call report(mype,nprocs,ndim,no_permanent_guardcells, & 
     &            advance_all_levels,lmpi, & 
     &            ierror_tot,.true., & 
     &            'bcset_example                 ')

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_close()

         
998      format('u:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
997      format('w:error proc ',i3,' block l= ',4(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
996      format('fx:error proc ',i3,' block l= ',5(2x,i3),2x,f8.5,2x, & 
     &       f8.5)
995      format('fy:error proc ',i3,' block l= ',5(2x,i3),2x,f8.5,2x, & 
     &       f8.5)
994      format('fz:error proc ',i3,' block l= ',5(2x,i3),2x,f8.5,2x, & 
     &       f8.5)
993      format('n:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
992      format('ez:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
991      format('ey:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
990      format('ex:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)

896      format('facex:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4, & 
     &       2x,f7.4)
895      format('facey:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4, & 
     &       2x,f7.4)
894      format('facez:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4, & 
     &       2x,f7.4)

893      format('ex_y:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
892      format('ex_z:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
891      format('ey_x:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
890      format('ey_z:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
889      format('ez_x:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)
888      format('ez_y:error proc ',i3,' block l= ',5(2x,i3),2x,f7.4,2x, & 
     &       f7.4)


      end
