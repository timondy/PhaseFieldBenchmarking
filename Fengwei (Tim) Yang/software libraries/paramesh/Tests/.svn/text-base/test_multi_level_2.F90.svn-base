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

!!!#define ALL_ON_PE1
#define PERIODIC

      program test_multi_level_2

! This test a 3D grid of 1x8x8 blocks of size 4x4x4 with periodic
! boundary conditions.
! x range 0 - 1
! y range 0 - 8
! z range 0 - 8
! The point of this test is to check that multiple root blocks
! can be defined.
!
! Test set up to diagnose problem detected by Rick DeVore
!                                              March 5, 2002.



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
     &                                amr_restrict, & 
     &                                amr_1blk_guardcell, & 
     &                                amr_1blk_guardcell_reset, & 
     &                                amr_1blk_restrict, & 
     &                                amr_flux_conserve, & 
     &                                amr_edge_average, & 
     &                                guardcell_test, & 
     &                                mesh_test, & 
     &                                amr_close, & 
     &                                amr_surrounding_blks, & 
     &                                amr_morton_order

      use paramesh_mpi_interfaces, only :  & 
     &                                mpi_morton_bnd, & 
     &                                mpi_amr_comm_setup, & 
     &                                mpi_amr_tree_setup, & 
     &                                mpi_amr_1blk_restrict, & 
     &                                mpi_morton_bnd_prolong

      use mpi_morton

! Only required for programs in ./Tests
#include "test_defs.fh"

      include 'mpif.h'

      integer :: tag_offset,max_blks_sent
      integer nguard0
      integer nguard_work0
      integer :: three,four,five, six
      character (len=80) :: filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local amr variables
      integer nprocs,mype
      integer num_procs

      save mype

      logical, allocatable :: tnewchild(:)
      logical, allocatable :: rflags(:)
      logical :: lfullblock
      logical :: lrelevant

!
! application specific variables

      real :: accuracy
      integer remote_block,remote_pe
      integer iopt,nlayers,icoord
      integer ierror_sum,ierror_tot
      integer cnodetype
      logical lrefine_again,ltype2only
      logical lcc, lfc, lec, lnc, l_srl_only, ldiag

      logical :: lmpi,lnperm,ladvanceall
      integer :: errorcode

      logical :: lguard,lprolong,lflux,ledge,lrestrict
      logical :: lfulltree

      integer :: surrblks(2,3,3,3)
      integer :: psurrblks(2,3,3,3)
      integer :: ierrorcode,ierr
      integer :: lnblocks_old

      save ierror_sum,ierror_tot
      save cnodetype
      save remote_block,remote_pe

      logical :: llperiodicx,llperiodicy,llperiodicz
      logical :: lpbc1,lpbc2,lpbc3
      logical :: lsbc1l,lsbc1r
      logical :: lsbc2l,lsbc2r
      logical :: lsbc3l,lsbc3r
      logical :: lebc1l,lebc1r
      logical :: lebc2l,lebc2r
      logical :: lebc3l,lebc3r
      integer :: pe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call amr_initialize

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      allocate(tnewchild(maxblocks_tr))
      allocate(rflags(maxblocks_tr))

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      ierror_sum = 0
      ierror_tot = 0

!!!#ifdef N_DIM == 3
      if (ndim == 3) then

#ifdef CONSERVE
      write(*,*) 'CONSERVE must not be defined for this test!'
      call amr_abort()
#endif


      accuracy = 100./10.**precision(accuracy)

      if(mype.eq.0) then
      print *,' nprocs = ',nprocs,mype


      if(nfacevar.eq.0) then
       write(*,*) ' This test requires nfacevar be', & 
     &       ' set to a value greater than 0. '
       stop
      endif
      endif

! set default value of dz and z0 to cater for 2D case.
      z0 = 0.
      dz = 0.

      iopt = 1
      nlayers = nguard
      if(mype.eq.0) write(*,*) 'nlayers = ',nlayers

!

! set a limit on the refinement level
      lrefine_max = 5
      lrefine_min = 1

      ax = .01
      ay = .1
      az = 1.
#ifdef TESTXDIR
      ax = .01
      ay = 0.
      az = 0.
#endif
#ifdef TESTYDIR
      ax = 0.
      ay = .1
      az = 0.
#endif
#ifdef TESTZDIR
      ax = 0.
      ay = 0.
      az = 1.
#endif


! set the workspace array layer to be tested
      ioptw = 3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up initial grid state.

      g_xmin = 0.
      g_xmax = 1.
      g_ymin = 0.
      g_ymax = 8.
      g_zmin = 0.
      g_zmax = 8.
      g_xmin = -2.
      g_xmax = 2.
      g_ymin = -16.
      g_ymax = 16.
      g_zmin = -16.
      g_zmax = 16.

      llperiodicx = .true.
      llperiodicy = .true.
      llperiodicz = .true.

! set up a 1x8x8 in 3D covering the domain 0 < x < 1,
! 0 < y < 8.0, 0 < z < 8.0.

      lnblocks = 0


#ifdef ALL_ON_PE1
      if(mype.eq.0.) then
                
                lnblocksx = 1
                lnblocksy = 8
                lnblocksz = 8
                lnblocks = lnblocksx*lnblocksy*lnblocksz
                bsize(:,1:64)=1.*4.

                do lbz = 1,lnblocksz
                do lby = 1,lnblocksy
                do lbx = 1,lnblocksx

                lb = lbx + (lby-1)*lnblocksx  + (lbz-1)*lnblocksy 
                coord(1,lb) = g_xmin+.5*4.
                coord(2,lb) = g_ymin + (.5 + real(lby-1))*4.
                coord(3,lb) = g_zmin + (.5 + real(lbz-1))*4.
                bnd_box(1,1,lb) = g_xmin+real(lbx-1)*4.
                bnd_box(2,1,lb) = g_xmin+real(lbx)*4.
                bnd_box(1,2,lb) = g_ymin+real(lby-1)*4.
                bnd_box(2,2,lb) = g_ymin+real(lby)*4.
                bnd_box(1,3,lb) = g_zmin+real(lbz-1)*4.
                bnd_box(2,3,lb) = g_zmin+real(lbz)*4.
                nodetype(lb) = 1
                lrefine(lb) = 1
                refine(lb)=.false.

                neigh(1,1:2,lb) = lb
                neigh(1,3,lb) = lbx  & 
     &                 + mod(lby-2+10*lnblocksy,lnblocksy)*lnblocksx   & 
     &                 + (lbz-1)*lnblocksy
                neigh(1,4,lb) = lbx  & 
     &                 + mod(lby,lnblocksy)*lnblocksx   & 
     &                 + (lbz-1)*lnblocksy
                neigh(1,5,lb) = lbx  & 
     &                 + (lby-1)*lnblocksx & 
     &                 + mod(lbz-2+10*lnblocksz,lnblocksz)*lnblocksy  
                neigh(1,6,lb) = lbx  & 
     &                 + (lby-1)*lnblocksx & 
     &                 + mod(lbz,lnblocksz)*lnblocksy  
                neigh(2,1:6,lb) = mype
!                write(*,*) 'block lb ',lb,' lbx lby lbz ',
!     .                        lbx,lby,lbz,' neighs ',
!     .                       neigh(1,:,lb)
                enddo
                enddo
                enddo



      endif

#else /* ALL_ON_PE1 */

! Ricks setup

      lpbc1 = .true.
      lpbc2 = .true.
      lpbc3 = .true.
      npes = nprocs
      nblockx = 1
      nblocky = 8
      nblockz = 8
      xrigh = 2.
      xleft = -2.
      yrigh = 16.
      yleft = -16.
      zrigh = 16.
      zleft = -16.
!
!     initialize the grid.
!
          nblocks  = nblockx * nblocky * nblockz
          nblock   = -1
          lnblocks = 0
!
          nbp = (nblocks+npes-1) / npes
          nbp = min (max (nbp, 2), nblocks)
!
!     distribute the nblockx*nblocky*nblockz blocks across the processors.
!
          nx = nblockx
          ny = nblocky
          nz = nblockz
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                nblock        = nblock + 1
                lb            = mod (nblock, nbp) + 1
                pe            = nblock / nbp
!
                if (mype .eq. pe) then
                  lnblocks      = lnblocks + 1
                  nodetype(lb)  = 1
                  lrefine(lb)   = 1
!
                  neigh(1,1,lb) = mod (nblock-1, nbp) + 1
                  neigh(2,1,lb) = (nblock-1) / nbp
                  neigh(1,2,lb) = mod (nblock+1, nbp) + 1
                  neigh(2,2,lb) = (nblock+1) / nbp
                  if (i .eq. 1) then
                    if (lpbc1) then
                      neigh(1,1,lb) = mod (nblock+nx-1, nbp) + 1
                      neigh(2,1,lb) = (nblock+nx-1) / nbp
                    else
                      if (lsbc1l) then
                        neigh(1,1,lb) = -30
                      else if (lebc1l) then
                        neigh(1,1,lb) = -31
                      else
                        neigh(1,1,lb) = -32
                      end if
                      neigh(2,1,lb) = pe
                    end if
                  end if
                  if (i .eq. nx) then
                    if (lpbc1) then
                      neigh(1,2,lb) = mod (nblock-nx+1, nbp) + 1
                      neigh(2,2,lb) = (nblock-nx+1) / nbp
                    else
                      if (lsbc1r) then
                        neigh(1,2,lb) = -30
                      else if (lebc1r) then
                        neigh(1,2,lb) = -31
                      else
                        neigh(1,2,lb) = -32
                      end if
                      neigh(2,2,lb) = pe
                    end if
                  end if
!
                  neigh(1,3,lb) = mod (nblock-nx, nbp) + 1
                  neigh(2,3,lb) = (nblock-nx) / nbp
                  neigh(1,4,lb) = mod (nblock+nx, nbp) + 1
                  neigh(2,4,lb) = (nblock+nx) / nbp
                  if (j .eq. 1) then
                    if (lpbc2) then
                      neigh(1,3,lb) = mod (nblock+nx*(ny-1), nbp) + 1
                      neigh(2,3,lb) = (nblock+nx*(ny-1)) / nbp
                    else
                      if (lsbc2l) then
                        neigh(1,3,lb) = -30
                      else if (lebc2l) then
                        neigh(1,3,lb) = -31
                      else
                        neigh(1,3,lb) = -32
                      end if
                      neigh(2,3,lb) = pe
                    end if
                  end if
                  if (j .eq. ny) then
                    if (lpbc2) then
                      neigh(1,4,lb) = mod (nblock-nx*(ny-1), nbp) + 1
                      neigh(2,4,lb) = (nblock-nx*(ny-1)) / nbp
                    else
                      if (lsbc2r) then
                        neigh(1,4,lb) = -30
                      else if (lebc2r) then
                        neigh(1,4,lb) = -31
                      else
                        neigh(1,4,lb) = -32
                      end if
                      neigh(2,4,lb) = pe
                    end if
                  end if
!
                  neigh(1,5,lb) = mod (nblock-nx*ny, nbp) + 1
                  neigh(2,5,lb) = (nblock-nx*ny) / nbp
                  neigh(1,6,lb) = mod (nblock+nx*ny, nbp) + 1
                  neigh(2,6,lb) = (nblock+nx*ny) / nbp
                  if (k .eq. 1) then
                    if (lpbc3) then
                      neigh(1,5,lb) = mod (nblock+nx*ny*(nz-1), nbp) + 1
                      neigh(2,5,lb) = (nblock+nx*ny*(nz-1)) / nbp
                    else
                      if (lsbc3l) then
                        neigh(1,5,lb) = -30
                      else if (lebc3l) then
                        neigh(1,5,lb) = -31
                      else
                        neigh(1,5,lb) = -32
                      end if
                      neigh(2,5,lb) = pe
                    end if
                  end if
                  if (k .eq. nz) then
                    if (lpbc3) then
                      neigh(1,6,lb) = mod (nblock-nx*ny*(nz-1), nbp) + 1
                      neigh(2,6,lb) = (nblock-nx*ny*(nz-1)) / nbp
                    else
                      if (lsbc3r) then
                        neigh(1,6,lb) = -30
                      else if (lebc3r) then
                        neigh(1,6,lb) = -31
                      else
                        neigh(1,6,lb) = -32
                      end if
                      neigh(2,6,lb) = pe
                    end if
                  end if
!
                  bsize(1,lb)     = (xrigh - xleft) / float(nx)
                  bnd_box(1,1,lb) = xleft + float(i-1) * bsize(1,lb)
                  bnd_box(2,1,lb) = bnd_box(1,1,lb) + bsize(1,lb)
                  coord(1,lb)     = 0.5 * (bnd_box(1,1,lb) & 
     &                                  +  bnd_box(2,1,lb))
                  bsize(2,lb)     = (yrigh - yleft) / float(ny)
                  bnd_box(1,2,lb) = yleft + float(j-1) * bsize(2,lb)
                  bnd_box(2,2,lb) = bnd_box(1,2,lb) + bsize(2,lb)
                  coord(2,lb)     = 0.5 * (bnd_box(1,2,lb) & 
     &                                  +  bnd_box(2,2,lb))
                  bsize(3,lb)     = (zrigh - zleft) / float(nz)
                  bnd_box(1,3,lb) = zleft + float(k-1) * bsize(3,lb)
                  bnd_box(2,3,lb) = bnd_box(1,3,lb) + bsize(3,lb)
                  coord(3,lb)     = 0.5 * (bnd_box(1,3,lb) & 
     &                                  +  bnd_box(2,3,lb))
                end if
              end do
            end do
          end do

#endif /* ALL_ON_PE1 */


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
      three= 2*k2d + 1
      four = 3*k2d + 1
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

      print *,' DONE with initial grid setup ',mype
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call amr_reorder_grid()

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

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!start test
! set the solution array to be the grid points x,y or z coordinates
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
       value = ax*x0 + ay*y0 + az*z0
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
       value = ax*x0 + ay*y0 + az*z0
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
                  do j=1+nguard0,nyb+nguard0+1
                    y0 = coord(2,l)-.5*bsize(2,l)-dy
                    yj = y0 + dy*real(j-nguard0)
                    do ivar=1,nbndvar
                      value = ax*xi+ay*yj+az*zk
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

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loop_count=0
! Now cycle over blocks adjusting refinement of initial setup as required
!      do while(loop_count.lt.2)
      do while(loop_count.lt.1)


       derefine(:) = .false.
       refine(:) = .false.
!       do l=1,lnblocks
!         if( loop_count.eq.0.and.abs(coord(1,l)).lt..6 )
!     .                              refine(l)=.true.
!         if( loop_count.eq.1.and.abs(coord(1,l)).lt..3 )
!     .                              refine(l)=.true.
!       enddo



! refine grid and apply morton reordering to grid blocks if necessary
      call amr_refine_derefine

      if (nvar > 0 .and. nvar_work > 0) then

      noff = (nguard_work - nguard)*npgs
      if(noff.ge.0) then
        work(il_bnd+noff:iu_bnd+noff,jl_bnd+noff*k2d:ju_bnd+noff*k2d, & 
     &       kl_bnd+noff*k3d:ku_bnd+noff*k3d,:,ioptw-1) =  & 
     &   unk(1,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,:)
      else
        work(ilw:iuw,jlw:juw,klw:kuw,:,ioptw-1) = & 
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

      if (nvar_work > 0) then
      iopt = ioptw
      nlayers = nguard_work
      lcc = .true.
      lfc = .false.
      lec = .false.
      lnc = .false.
      tag_offset = 100
      call amr_guardcell(mype,iopt,nlayers)
      end if

      else                      ! no_permanent_guardcells

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

      endif ! no_permanent_guardcells

        loop_count=loop_count+1
        write(*,*) 'proc loop_count ',mype,loop_count

        enddo
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

#ifdef FLUX
! Set up flux variables for testing

      flux_x(:,:,:,:,:) = 0.
      flux_y(:,:,:,:,:) = 0.
      flux_z(:,:,:,:,:) = 0.

      do l=1,lnblocks

      if(nodetype(l).eq.1 .or. advance_all_levels) then

      do ivar=1,nfluxes

      flux_x(ivar,1,jl_bndi:ju_bndi,kl_bndi:ku_bndi,l) =  & 
     & facevarx(1,1+nguard0,jl_bndi:ju_bndi,kl_bndi:ku_bndi,l)

      flux_x(ivar,2,jl_bndi:ju_bndi,kl_bndi:ku_bndi,l) =  & 
     & facevarx(1,1+nguard0+nxb,jl_bndi:ju_bndi,kl_bndi:ku_bndi,l)

      flux_y(ivar,il_bndi:iu_bndi,1,kl_bndi:ku_bndi,l) =  & 
     & facevary(1,il_bndi:iu_bndi,1+nguard0*k2d,kl_bndi:ku_bndi,l)

      flux_y(ivar,il_bndi:iu_bndi,2,kl_bndi:ku_bndi,l) =  & 
     & facevary(1,il_bndi:iu_bndi,1+(nguard0+nyb)*k2d,kl_bndi:ku_bndi,l)

      if(ndim.eq.3) then

      flux_z(ivar,il_bndi:iu_bndi,jl_bndi:ju_bndi,1,l) =  & 
     & facevarz(1,il_bndi:iu_bndi,jl_bndi:ju_bndi,1+nguard0*k3d,l)

      flux_z(ivar,il_bndi:iu_bndi,jl_bndi:ju_bndi,2,l) =  & 
     & facevarz(1,il_bndi:iu_bndi,jl_bndi:ju_bndi,1+(nguard0+nzb)*k3d,l)

      endif

      enddo

      endif

      enddo



! get remote tree info needed for zeroing fluxes and edges.
      tag_offset = 100

      do l=1,lnblocks
      if(nodetype(l).eq.1) then

! zero out fluxes on those interior faces bordering finer blocks
      remote_block = neigh(1,1,l) 
      if(remote_block.gt.0) then
      remote_pe = neigh(2,1,l) 
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
            endif
          enddo
        endif
      cnodetype = nodetype(remote_block)
      if(cnodetype.eq.2) flux_x(:,1,:,:,l) = 0.
      endif
! 
      remote_block = neigh(1,2,l) 
      if(remote_block.gt.0) then
      remote_pe = neigh(2,2,l) 
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
            endif
          enddo
        endif
      cnodetype = nodetype(remote_block)
      if(cnodetype.eq.2) flux_x(:,2,:,:,l) = 0.
      endif
! 
      remote_block = neigh(1,3,l) 
      if(remote_block.gt.0) then
      remote_pe = neigh(2,3,l) 
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
            endif
          enddo
        endif
      cnodetype = nodetype(remote_block)
      if(cnodetype.eq.2) flux_y(:,:,1,:,l) = 0.
      endif
! 
      remote_block = neigh(1,4,l) 
      if(remote_block.gt.0) then
      remote_pe = neigh(2,4,l) 
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
            endif
          enddo
        endif
      cnodetype = nodetype(remote_block) 
      if(cnodetype.eq.2) flux_y(:,:,1+k2d,:,l) = 0.
      endif
! 
      if(ndim.eq.3) then

      remote_block = neigh(1,5,l) 
      if(remote_block.gt.0) then
      remote_pe = neigh(2,5,l) 
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
            endif
          enddo
        endif
      cnodetype = nodetype(remote_block)
      if(cnodetype.eq.2) flux_z(:,:,:,1,l) = 0.
      endif
! 
      remote_block = neigh(1,6,l) 
      if(remote_block.gt.0) then
      remote_pe = neigh(2,6,l) 
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
            endif
          enddo
        endif
      cnodetype = nodetype(remote_block)
      if(cnodetype.eq.2) flux_z(:,:,:,1+k3d,l) = 0.
      endif
 
      endif

      endif
      enddo

#endif /* FLUX */
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up data on block boundary edges

! initialize
        bedge_facex_y = 0.
        bedge_facex_z = 0.
        bedge_facey_x = 0.
        bedge_facey_z = 0.
        bedge_facez_x = 0.
        bedge_facez_y = 0.

!----------------------------------------------------



#if defined(FLUX) || defined(EDGE_TEST)

      if (.not.consv_flux_densities) then
      if(mype.eq.0) then
         write(*,*) ' '
         write(*,*) 'This test will fail because CONSV_FLUX_DENSITIES'
         write(*,*) 'is undefined.'
         write(*,*) 'Edit ../headers/physicaldata.fh.'
         write(*,*) ' '
      endif 
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_abort
      endif

#endif


      iopt = 1



      if (no_permanent_guardcells) then
! Store a copy of the current solution in gt_unk
      call amr_1blk_copy_soln(-1)
      endif                     ! no_permanent_guardcells


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
      endif                     ! no_permanent_guardcells


      if (nvar > 0 .and. nvar_work > 0) then
      noff = (nguard_work - nguard)*npgs
      if(noff.ge.0) then
        work(il_bnd+noff:iu_bnd+noff,jl_bnd+noff*k2d:ju_bnd+noff*k2d, & 
     &       kl_bnd+noff*k3d:ku_bnd+noff*k3d,:,ioptw-1) = & 
     &   unk(1,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,:)
      else
        work(ilw:iuw,jlw:juw,klw:kuw,:,ioptw-1) = & 
     &   unk(1,il_bnd-noff:iu_bnd+noff, & 
     &         jl_bnd-noff*k2d:ju_bnd+noff*k2d, & 
     &         kl_bnd-noff*k3d:ku_bnd+noff*k3d,:)
      endif
      end if


      if(mype.eq.0) write(*,*) 'Start of automatic testing.'
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
      if (diagonals) then
         ldiag = .true.
      end if

      if(nodetype(l).eq.1 .or. advance_all_levels) then

      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)

      endif

      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
      unk1(:,:,:,:,1) = unk(:,:,:,:,l)
#endif
      endif                     ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or. advance_all_levels) then

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

      do k=klbnd,kubnd
      if(ndim.eq.3) then
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       if(mod(nxb,2).eq.1) z0 = coord(3,l)-(.5*bsize(3,l)+dz)
      endif
      zk = z0 + dz*real(k-nguard)
      if(llperiodicz) then
       if(zk.lt.g_zmin) zk = g_zmax - (g_zmin - zk)
       if(zk.gt.g_zmax) zk = g_zmin + (zk - g_zmax)
      endif
       do j=jlbnd,jubnd
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)
       if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       yj = y0 + dy*real(j-nguard)
       if(llperiodicy) then
        if(yj.lt.g_ymin) yj = g_ymax - (g_ymin - yj)
        if(yj.gt.g_ymax) yj = g_ymin + (yj - g_ymax)
       endif
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)
       if(mod(nxb,2).eq.1) x0 =  & 
     &       coord(1,l)-(.5*bsize(1,l)+dx)
       xi = x0 + dx*real(i-nguard)
       if(llperiodicx) then
        if(xi.lt.g_xmin) xi = g_xmax - (g_xmin - xi)
        if(xi.gt.g_xmax) xi = g_xmin + (xi - g_xmax)
       endif
       value = ax*xi+ay*yj+az*zk
       do ivar=1,nvar
       v0 = value*real(ivar)
       if(abs(v0-unk1(ivar,i,j,k,1))>accuracy) then
         write(*,998) ii,l,ivar,i,j,k,unk1(ivar,i,j,k,1),v0
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
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


#ifndef PERIODIC

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

      if(nodetype(l).eq.1 .or. advance_all_levels) then

      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)

      endif

      else                      ! no_permanent_guardcells
      work1(:,:,:,1) = work(:,:,:,l,ioptw-1)
      endif                     ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or. advance_all_levels) then

      ilbnd=ilw1
      iubnd=iuw1
      jlbnd=jlw1
      jubnd=juw1
      klbnd=klw1
      kubnd=kuw1

      if(neigh(1,1,l).le.-20) ilbnd=1+nguard_work
      if(neigh(1,2,l).le.-20) iubnd=nxb+nguard_work
      if(neigh(1,3,l).le.-20) jlbnd=1+nguard_work
      if(neigh(1,4,l).le.-20) jubnd=nyb+nguard_work
      if(ndim.eq.3) then
      if(neigh(1,5,l).le.-20) klbnd=1+nguard_work
      if(neigh(1,6,l).le.-20) kubnd=nzb+nguard_work
      endif

      if(ndim.eq.3) dz = bsize(3,l)/real(nzb)
      dy = bsize(2,l)/real(nyb)
      dx = bsize(1,l)/real(nxb)
      if(mod(nxb,2).eq.1) then
       if(ndim.eq.3) dz = bsize(3,l)/real(nzb-k3d)
       dy = bsize(2,l)/real(nyb-1)
       dx = bsize(1,l)/real(nxb-1)
      endif

      do k=klbnd,kubnd
      if(ndim.eq.3) then
       z0 = coord(3,l)-.5*(bsize(3,l)+dz)
       if(mod(nxb,2).eq.1) z0 = coord(3,l)-(.5*bsize(3,l)+dz)
      endif
      zk = z0 + dz*real(k-nguard_work)
       do j=jlbnd,jubnd
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)
       if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       yj = y0 + dy*real(j-nguard_work)
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)
       if(mod(nxb,2).eq.1) x0 =  & 
     &       coord(1,l)-(.5*bsize(1,l)+dx)
       xi = x0 + dx*real(i-nguard_work)
       value = ax*xi+ay*yj+az*zk
       if(abs(value-work1(i,j,k,1))>accuracy) then
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
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


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

      if(nodetype(l).eq.1 .or. advance_all_levels) then
 
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
      unk_n1(:,:,:,:,1) = unk_n(:,:,:,:,l)
      endif                     ! no_permanent_guardcells
 
 
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
              if(abs(unk_n1(ivar,i,j,k,1)-v0)>accuracy) then
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
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)                            
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end if

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef FACE

      if (nfacevar > 0) then

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
      icoord = 0
      ldiag = .false.
      ldiag = diagonals

      if(nodetype(l).eq.1 .or. advance_all_levels) then

      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)

      endif

      else                      ! no_permanent_guardcells
      facevarx1(:,:,:,:,1) = facevarx(:,:,:,:,l)
      facevary1(:,:,:,:,1) = facevary(:,:,:,:,l)
      facevarz1(:,:,:,:,1) = facevarz(:,:,:,:,l)

      endif                     ! no_permanent_guardcells


      if(nodetype(l).eq.1 .or. advance_all_levels) then

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

! first test facevarx
#ifdef FACEX
      do k=klbnd,kubnd
        if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
        zk = z0 + dz*real(k-nguard)
        do j=jlbnd,jubnd
          y0 = coord(2,l)-.5*(bsize(2,l)+dy)
          yj = y0 + dy*real(j-nguard)
          do i=ilbnd,iubnd+ione
            x0 = coord(1,l)-.5*bsize(1,l)-dx
            xi = x0 + dx*real(i-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevarx1(ivar,i,j,k,1)-v0)>accuracy) then
                      write(*,896) mype,l,ivar,i,j,k, & 
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

      do k=klbnd,kubnd
        if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
        zk = z0 + dz*real(k-nguard)
        do i=ilbnd,iubnd
          x0 = coord(1,l)-.5*(bsize(1,l)+dx)
          xi = x0 + dx*real(i-nguard)
          do j=jlbnd,jubnd+ione
            y0 = coord(2,l)-.5*bsize(2,l)-dy
            yj = y0 + dy*real(j-nguard)
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevary1(ivar,i,j,k,1)-v0)>accuracy) then
                      write(*,895) mype,l,ivar,i,j,k, & 
     &                             facevary1(ivar,i,j,k,1),v0
                      ierror_sum = ierror_sum + 1
              endif
            enddo
          enddo
        enddo
      enddo
#endif

! finally test facevarz
#ifdef FACEZ

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
              if(abs(facevarz1(ivar,i,j,k,1)-v0)>accuracy) then
                     write(*,894) mype,l,ivar,i,j,k, & 
     &                             facevarz1(ivar,i,j,k,1),v0
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
      if(mype.eq.0) write(*,*) 'face var. guard cell test complete.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end if

#endif /* FACE */


#ifdef FLUX

! test of flux conservation at block boundaries (if appropriate)
! ensure flux conservation
      nsub=1
      call amr_flux_conserve(mype,nsub)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! test values 
      do ii=0,nprocs-1
      if(mype.eq.ii) then

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

! first test flux_x
#ifdef FLUXX

        do k=1+nguard0*k3d,nzb+nguard0*k3d
          if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
          zk = z0 + dz*real(k-nguard0)
          do j=1+nguard0,nyb+nguard0
            y0 = coord(2,l)-.5*(bsize(2,l)+dy)
            yj = y0 + dy*real(j-nguard0)
            do i=1,2
              xi = coord(1,l)+.5*bsize(1,l)*real(2*i-3)
              value = ax*xi+ay*yj+az*zk
              do ivar=1,nfluxes
                v0 = value
                if(abs(flux_x(ivar,i,j,k,l)-v0)>accuracy) then
                  write(*,996) mype,l,ivar,i,j,k, & 
     &                         flux_x(ivar,i,j,k,l),v0
                  ierror_sum = ierror_sum + 1
                endif
              enddo
            enddo
          enddo
        enddo
#endif

! now test flux_y
#ifdef FLUXY

        do k=1+nguard0*k3d,nzb+nguard0*k3d
          if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
          zk = z0 + dz*real(k-nguard0)
          do i=1+nguard0,nxb+nguard0
            x0 = coord(1,l)-.5*(bsize(1,l)+dx)
            xi = x0 + dx*real(i-nguard0)
            do j=1,2
              yj = coord(2,l)+.5*bsize(2,l)*real(2*j-3)
              value = ax*xi+ay*yj+az*zk
              do ivar=1,nfluxes
                v0 = value
                if(abs(flux_y(ivar,i,j,k,l)-v0)>accuracy) then
                  write(*,995) mype,l,ivar,i,j,k, & 
     &                flux_y(ivar,i,j,k,l),v0
                  ierror_sum = ierror_sum + 1
                endif
              enddo
            enddo
          enddo
        enddo
#endif

! finally test flux_z
#ifdef FLUXZ
      if(ndim.eq.3) then

        do j=1+nguard0,nyb+nguard0
          y0 = coord(2,l)-.5*(bsize(2,l)+dy)
          yj = y0 + dy*real(j-nguard0)
          do i=1+nguard0,nxb+nguard0
            x0 = coord(1,l)-.5*(bsize(1,l)+dx)
            xi = x0 + dx*real(i-nguard0)
            do k=1,2
              zk = coord(3,l)+.5*bsize(3,l)*real(2*k-3)
              value = ax*xi+ay*yj+az*zk
              do ivar=1,nfluxes
                v0 = value
                if(abs(flux_z(ivar,i,j,k,l)-v0)>accuracy) then
                  write(*,994) mype,l,ivar,i,j,k, & 
     &                    flux_z(ivar,i,j,k,l),v0
                  ierror_sum = ierror_sum + 1
                endif
              enddo
            enddo
          enddo
        enddo

        endif
#endif

      endif

      enddo
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'flux conservation test complete.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#endif  /* FLUX */

      else
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'Flux conservation test omitted.', & 
     &                         'Not appropriate with odd number ', & 
     &                         'of mesh points per block'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef EDGES

      if (nvaredge > 0) then

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
      icoord = 0
      ldiag = .false.
      ldiag = diagonals
      call amr_1blk_guardcell(mype,iopt,nlayers,l,mype, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag)
 
      else                      ! no_permanent_guardcells
      unk_e_x1(:,:,:,:,1) = unk_e_x(:,:,:,:,l)
      unk_e_y1(:,:,:,:,1) = unk_e_y(:,:,:,:,l)
      unk_e_z1(:,:,:,:,1) = unk_e_z(:,:,:,:,l)
 
      endif                     ! no_permanent_guardcells

 
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
 
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up data on block boundary edges

! initialize temporary block face arrays for edge-centered data
        bedge_facex_y = 0.
        bedge_facex_z = 0.
        bedge_facey_x = 0.
        bedge_facey_z = 0.
        bedge_facez_x = 0.
        bedge_facez_y = 0.

!----------------------------------------------------

! Zero out cell-edge data which will be modified in the call
! to amr_edge_average to ensure a meaningful test.

#ifdef EDGE_TEST

        if(mod(nxb,2).eq.0) then

        do l=1,lnblocks
        if(nodetype(l).eq.1) then


        ia = 1+nguard0
        ib = nxb+nguard0
        ja = 1+nguard0*k2d
        jb = nyb+nguard0*k2d
        ka = 1+nguard0*k3d
        kb = nzb+nguard0*k3d

! zero out edge values on those interior faces bordering finer blocks
        remote_block = neigh(1,1,l)
        if(remote_block.gt.0) then
        remote_pe = neigh(2,1,l)
        cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
        if(cnodetype.eq.2) then
          unk_e_y(:,ia,ja:jb,ka:kb+k3d,l) = 0.
          if(ndim.eq.3) unk_e_z(:,ia,ja:jb+k2d,ka:kb,l) = 0.
        endif
        endif
!
        remote_block = neigh(1,2,l)
        if(remote_block.gt.0) then
        remote_pe = neigh(2,2,l)
        cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
        if(cnodetype.eq.2) then
          unk_e_y(:,ib+1,ja:jb,ka:kb+k3d,l) = 0.
          if(ndim.eq.3) unk_e_z(:,ib+1,ja:jb+k2d,ka:kb,l) = 0.
        endif
        endif
!
        remote_block = neigh(1,3,l)
        if(remote_block.gt.0) then
        remote_pe = neigh(2,3,l)
        cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
        if(cnodetype.eq.2) then
          unk_e_x(:,ia:ib,ja,ka:kb+k3d,l) = 0.
          if(ndim.eq.3) unk_e_z(:,ia:ib+1,ja,ka:kb,l) = 0.
        endif
        endif

!
        remote_block = neigh(1,4,l)
        if(remote_block.gt.0) then
        remote_pe = neigh(2,4,l)
        cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
        if(cnodetype.eq.2) then
          unk_e_x(:,ia:ib,jb+k2d,ka:kb+k3d,l) = 0.
          if(ndim.eq.3) unk_e_z(:,ia:ib+1,jb+k2d,ka:kb,l) = 0.
        endif
        endif
!

      if(ndim.eq.3) then

        remote_block = neigh(1,5,l)
        if(remote_block.gt.0) then
        remote_pe = neigh(2,5,l)
        cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
        if(cnodetype.eq.2) then
          unk_e_x(:,ia:ib,ja:jb+k2d,ka,l) = 0.
          unk_e_y(:,ia:ib+1,ja:jb,ka,l) = 0.
        endif
        endif
!
        remote_block = neigh(1,6,l)
        if(remote_block.gt.0) then
        remote_pe = neigh(2,6,l)
        cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
        if(cnodetype.eq.2) then
          unk_e_x(:,ia:ib,ja:jb+k2d,kb+k3d,l) = 0.
          unk_e_y(:,ia:ib+1,ja:jb,kb+k3d,l) = 0.
        endif
        endif

        endif           ! end of ndim=3 iftest

!
! call surrounding blocks to get addresses of diagonal neighbors
      call MPI_ALLREDUCE(lnblocks, & 
     &                   max_no_of_blocks, & 
     &                   1, & 
     &                   MPI_INTEGER, & 
     &                   MPI_MAX, & 
     &                   MPI_COMM_WORLD, & 
     &                   ierror)
       surrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =       &
          surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,l)

! loop over diagonal neighbors and identify their locations
       do ie = 1 , 12

         cnodetype = -1

         if(ie.eq.1) then
            if (ndim == 3) then

           remote_block = surrblks(1,1,1,2)
           remote_pe    = surrblks(2,1,1,2)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             if(ndim.eq.3) unk_e_z(:,ia,ja,ka:kb,l) = 0.
           endif

         elseif(ie.eq.2) then

           remote_block = surrblks(1,1,3,2)
           remote_pe    = surrblks(2,1,3,2)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             if(ndim.eq.3) unk_e_z(:,ia,jb+k2d,ka:kb,l) = 0.
           endif

         elseif(ie.eq.3) then

           remote_block = surrblks(1,3,1,2)
           remote_pe    = surrblks(2,3,1,2)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             if(ndim.eq.3) unk_e_z(:,ib+1,ja,ka:kb,l) = 0.
           endif

         elseif(ie.eq.4) then

           remote_block = surrblks(1,3,3,2)
           remote_pe    = surrblks(2,3,3,2)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             if(ndim.eq.3) unk_e_z(:,ib+1,jb+k2d,ka:kb,l) = 0.
           endif

           end if

         elseif(ie.eq.5) then

           remote_block = surrblks(1,2,1,2-k3d)
           remote_pe    = surrblks(2,2,1,2-k3d)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_x(:,ia:ib,ja,ka,l) = 0.
           endif

         elseif(ie.eq.6) then

           remote_block = surrblks(1,2,3,2-k3d)
           remote_pe    = surrblks(2,2,3,2-k3d)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
        if(remote_block.gt.0) then
        cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_x(:,ia:ib,jb+k2d,ka,l) = 0.
           endif

         elseif(ie.eq.7) then

            if (ndim == 3) then
           remote_block = surrblks(1,2,1,3)
           remote_pe    = surrblks(2,2,1,3)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_x(:,ia:ib,ja,kb+k3d,l) = 0.
           endif

         elseif(ie.eq.8) then

           remote_block = surrblks(1,2,3,3)
           remote_pe    = surrblks(2,2,3,3)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_x(:,ia:ib,jb+k2d,kb+k3d,l) = 0.
           endif
           end if

         elseif(ie.eq.9) then
           remote_block = surrblks(1,1,2,2-k3d)
           remote_pe    = surrblks(2,1,2,2-k3d)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_y(:,ia,ja:jb,ka,l) = 0.
           endif

         elseif(ie.eq.10) then

           if (ndim == 3) then
           remote_block = surrblks(1,1,2,3)
           remote_pe    = surrblks(2,1,2,3)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_y(:,ia,ja:jb,kb+k3d,l) = 0.
           endif
           end if

         elseif(ie.eq.11) then

           remote_block = surrblks(1,3,2,2-k3d)
           remote_pe    = surrblks(2,3,2,2-k3d)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_y(:,ib+1,ja:jb,ka,l) = 0.
           endif

         elseif(ie.eq.12) then

           if (ndim == 3) then
           remote_block = surrblks(1,3,2,3)
           remote_pe    = surrblks(2,3,2,3)
           cnodetype = -1
        if(remote_pe.ne.mype) then
          do iblk = strt_buffer,last_buffer
            if( remote_block.eq.laddress(1,iblk) .and. & 
     &             remote_pe.eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              cnodetype = nodetype(remote_block)
            endif
          enddo
        else
          cnodetype = nodetype(remote_block)
        endif
           if(cnodetype.eq.2) then
             unk_e_y(:,ib+1,ja:jb,kb+k3d,l) = 0.
           endif
           end if

         endif                      ! end of ie if test

       enddo                      ! end of loop over edges index ie

        endif
      enddo

      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!--------


      if(mod(nxb,2).eq.0) then


! adjust cell-edge data to give consistent circulation integrals 
      nsub=1
      lfullblock = .true.
      call amr_edge_average(mype,lfullblock,nsub)

!
! test of data on block boundary edges
      do ii=0,nprocs-1
      if(mype.eq.ii) then


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
         if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
         zk = z0 + dz*real(k-nguard0)
         do j=1+nguard0*k2d,nyb+nguard0*k2d
           y0 = coord(2,l)-.5*(bsize(2,l)+dy)
           yj = y0 + dy*real(j-nguard0)
           do i=1,2
             xi=coord(1,l)+.5*bsize(1,l)*real(2*i-3)
             value = ax*xi+ay*yj+az*zk
             do ivar=1,nedges
               vv = value*real(ivar)
               if(abs(bedge_facex_y(ivar,i,j,k,l)-vv)>accuracy) then
                 write(*,893) & 
     &             mype,l,ivar,i,j,k,bedge_facex_y(ivar,i,j,k,l),vv
                 ierror_sum = ierror_sum + 1
               endif
             enddo
           enddo
         enddo
       enddo

       if(ndim.eq.3) then
         do k=1+nguard0*k3d,nzb+nguard0*k3d
           z0 = coord(3,l)-.5*(bsize(3,l)+dz)
           zk = z0 + dz*real(k-nguard0)
           do j=1+nguard0*k2d,nyb+(nguard0+1)*k2d
             y0 = coord(2,l)-.5*bsize(2,l)-dy
             yj = y0 + dy*real(j-nguard0)
             do i=1,2
               xi=coord(1,l)+.5*bsize(1,l)*real(2*i-3)
               value = ax*xi+ay*yj+az*zk
               do ivar=1,nedges
                 vv = value*real(ivar)
                 if(abs(bedge_facex_z(ivar,i,j,k,l)-vv)>accuracy) then
                   write(*,892) & 
     &               mype,l,ivar,i,j,k,bedge_facex_z(ivar,i,j,k,l),vv
                   ierror_sum = ierror_sum + 1
                 endif
               enddo
             enddo
           enddo
         enddo
       endif

       do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d
         if(ndim.eq.3) z0 = coord(3,l)-.5*bsize(3,l)-dz
         zk = z0 + dz*real(k-nguard0)
         do i=1+nguard0,nxb+nguard0
           x0 = coord(1,l)-.5*(bsize(1,l)+dx)
           xi = x0 + dx*real(i-nguard0)
           do j=1,2
             yj=coord(2,l)+.5*bsize(2,l)*real(2*j-3)
             value = ax*xi+ay*yj+az*zk
             do ivar=1,nedges
               vv = value*real(ivar)
               if(abs(bedge_facey_x(ivar,i,j,k,l)-vv)>accuracy) then
                 write(*,891) & 
     &         mype,l,ivar,i,j,k,bedge_facey_x(ivar,i,j,k,l),vv
                 ierror_sum = ierror_sum + 1
               endif
             enddo
           enddo
         enddo
       enddo


       if(ndim.eq.3) then

       do k=1+nguard0*k3d,nzb+nguard0*k3d
         z0 = coord(3,l)-.5*(bsize(3,l)+dz)
         zk = z0 + dz*real(k-nguard0)
         do i=1+nguard0,nxb+nguard0+1
           x0 = coord(1,l)-.5*bsize(1,l)-dx
           xi = x0 + dx*real(i-nguard0)
           do j=1,2
             yj=coord(2,l)+.5*bsize(2,l)*real(2*j-3)
             value = ax*xi+ay*yj+az*zk
             do ivar=1,nedges
               vv = value*real(ivar)
               if(abs(bedge_facey_z(ivar,i,j,k,l)-vv)>accuracy) then
                 write(*,890) & 
     &         mype,l,ivar,i,j,k,bedge_facey_z(ivar,i,j,k,l),vv
                 ierror_sum = ierror_sum + 1
               endif
             enddo
           enddo
         enddo
       enddo

       do j=1+nguard0*k2d,nyb+(nguard0+1)*k2d
         y0 = coord(2,l)-.5*bsize(2,l)-dy
         yj = y0 + dy*real(j-nguard0)
         do i=1+nguard0,nxb+nguard0
           x0 = coord(1,l)-.5*(bsize(1,l)+dx)
           xi = x0 + dx*real(i-nguard0)
           do k=1,2
             zk=coord(3,l)+.5*bsize(3,l)*real(2*k-3)
             value = ax*xi+ay*yj+az*zk
             do ivar=1,nedges
               vv = value*real(ivar)
               if(abs(bedge_facez_x(ivar,i,j,k,l)-vv)>accuracy) then
                 write(*,889) & 
     &         mype,l,ivar,i,j,k,bedge_facez_x(ivar,i,j,k,l),vv
                 ierror_sum = ierror_sum + 1
               endif
             enddo
           enddo
         enddo
       enddo

       do j=1+nguard0*k2d,nyb+nguard0*k2d
         y0 = coord(2,l)-.5*(bsize(2,l)+dy)
         yj = y0 + dy*real(j-nguard0)
         do i=1+nguard0,nxb+nguard0+1
           x0 = coord(1,l)-.5*bsize(1,l)-dx
           xi = x0 + dx*real(i-nguard0)
           do k=1,2
             zk=coord(3,l)+.5*bsize(3,l)*real(2*k-3)
             value = ax*xi+ay*yj+az*zk
             do ivar=1,nedges
               vv = value*real(ivar)
               if(abs(bedge_facez_y(ivar,i,j,k,l)-vv)>accuracy) then
                 write(*,888) & 
     &         mype,l,ivar,i,j,k,bedge_facez_y(ivar,i,j,k,l),vv
                 ierror_sum = ierror_sum + 1
               endif
             enddo
           enddo
         enddo
       enddo

              endif            ! end of ndim=3 if test


        endif

      enddo
      endif
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'edge averaging test complete.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      else
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'Edge average test omitted.', & 
     &                         'Not appropriate with odd number ', & 
     &                         'of mesh points per block'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif

#endif /* EDGE_TEST */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif /* PERIODIC */

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
      filename = trim(output_dir) // 'test.log'
      open(unit = 55,file=filename, & 
     &            status='unknown',position='append')
        write(*,*) ' '
        write(*,*) ' '
        if(ierror_tot.eq.0) then
          write(*,*) 'No errors detected - Test Successful '
          write(55,*) 'No errors detected - ', & 
     &                'Test of multi_level_2 Successful ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        else
          write(*,*) ierror_tot,' errors detected - Test failed '
          write(55,*) ierror_tot,' errors detected - ', & 
     &                'Test of multi_level_2 failed ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        endif
      close(55)
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
      lrelevant = .true.
      else
        ierror_tot = 0
        lrelevant = .false.
        write(*,*) 'N_DIM must be 3 for this test'
      end if

      call report(mype,nprocs,ndim,no_permanent_guardcells, & 
     &            advance_all_levels,lmpi, & 
     &            ierror_tot, & 
     &            lrelevant, & 
     &            'multi_level_2                 ')

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

      deallocate(tnewchild)
      deallocate(rflags)

      end
