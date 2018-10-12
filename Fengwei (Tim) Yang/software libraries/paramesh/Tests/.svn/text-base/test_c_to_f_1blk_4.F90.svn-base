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

!#define KEVINS_2D_PATTERN

      program test_c_to_f




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
      integer :: order
      integer :: three,four,five, six

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local amr variables
      integer nprocs,mype
      integer num_procs

      save mype

      logical, allocatable :: tnewchild(:)
      logical, allocatable :: rflags(:)

!
! application specific variables

      real :: accuracy 
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

      character (len=80) :: filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call amr_initialize 

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

! reset interpolation order(s)

      interp_mask_unk(:) = 4
      interp_mask_work(:) = 4
      if (nfacevar > 0) then
       interp_mask_facex(:) = 4
       interp_mask_facey(:) = 4
       interp_mask_facez(:) = 4
      end if
      if (nvaredge > 0) then
       interp_mask_ec(:) = 4
      end if
      if (nvarcorn > 0) then
       interp_mask_nc(:) = 4
      end if

      interp_mask_unk_res(:) = 4
      interp_mask_work_res(:) = 4
      if (nfacevar > 0) then
       interp_mask_facex_res(:) = 4
       interp_mask_facey_res(:) = 4
       interp_mask_facez_res(:) = 4
      end if
      if (nvaredge > 0) then
       interp_mask_ec_res(:) = 4
      end if
      if (nvarcorn > 0) then
       interp_mask_nc_res(:) = 4
      end if

      if (curvilinear_conserve) then
         print *,' resetting interp_mask for curvilinear_conserve '
         interp_mask_unk(:) = 1
         interp_mask_unk_res(:) = 1
         interp_mask_work(:) = 1
         interp_mask_work_res(:) = 1
         interp_mask_facex(:) = 1
         interp_mask_facex_res(:) = 1
         interp_mask_facey(:) = 1
         interp_mask_facey_res(:) = 1
         interp_mask_facez(:) = 1
         interp_mask_facez_res(:) = 1
         interp_mask_ec(:) = 1
         interp_mask_ec_res(:) = 1
         interp_mask_nc(:) = 1
         interp_mask_nc_res(:) = 1
      end if

!      if (nguard < 2) then
!         print *,' ERROR: nguard must be >= 2 for this test '
!         call amr_abort()
!      end if

      allocate(tnewchild(maxblocks_tr))
      allocate(rflags(maxblocks_tr))

#ifdef CONSERVE
      write(*,*) 'CONSERVE must not be defined for this test!'
      call amr_abort()
#endif

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      accuracy = 100000./10.**precision(accuracy)

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
      lrefine_max = 4
      lrefine_min = 1

!      if (.not. advance_all_levels) then
! if a refinement jump exists, a restriction is needed unless all levels
! are being advanced. This will not give the correct answer until we provide
! a second order accurate restriction routine. Setting the maximum refinement level
! to 3 will prevent refinement jumps.
!      lrefine_max = 3
!      write(*,*) 'WARNING !!! - NO REFINEMENT JUMPS IN THIS TEST'
!      write(*,*) 'WARNING !!! - NO REFINEMENT JUMPS IN THIS TEST'
!      write(*,*) 'WARNING !!! - NO REFINEMENT JUMPS IN THIS TEST'
!      end if

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
       do ivar=1,nvar
       value = ax*(x0**interp_mask_unk(ivar)) +  & 
     &         ay*(y0**interp_mask_unk(ivar)) +  & 
     &         az*(z0**interp_mask_unk(ivar))
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
       value = ax*(x0**interp_mask_work(ioptw-1)) +  & 
     &         ay*(y0**interp_mask_work(ioptw-1)) +  & 
     &         az*(z0**interp_mask_work(ioptw-1))
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
       do ivar=1,nvarcorn
       value = ax*x0**interp_mask_nc(ivar) +  & 
     &         ay*y0**interp_mask_nc(ivar) +  & 
     &         az*z0**interp_mask_nc(ivar)
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
                        value = ax*xi**interp_mask_facex(ivar)+ & 
     &                          ay*yj**interp_mask_facex(ivar)+ & 
     &                          az*zk**interp_mask_facex(ivar)
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
                      value = ax*xi**interp_mask_facey(ivar)+ & 
     &                        ay*yj**interp_mask_facey(ivar)+ & 
     &                        az*zk**interp_mask_facey(ivar)
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
                      value = ax*xi**interp_mask_facez(ivar)+ & 
     &                        ay*yj**interp_mask_facez(ivar)+ & 
     &                        az*zk**interp_mask_facez(ivar)
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
                        value = ax*xi**interp_mask_ec(ivar)+ & 
     &                          ay*yj**interp_mask_ec(ivar)+ & 
     &                          az*zk**interp_mask_ec(ivar)
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
                      value = ax*xi**interp_mask_ec(ivar)+ & 
     &                        ay*yj**interp_mask_ec(ivar)+ & 
     &                        az*zk**interp_mask_ec(ivar)
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
                      value = ax*xi**interp_mask_ec(ivar)+ & 
     &                        ay*yj**interp_mask_ec(ivar)+ & 
     &                        az*zk**interp_mask_ec(ivar)
                      unk_e_z(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo 

        endif

      enddo

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loop_count=0
! Now cycle over blocks adjusting refinement of initial setup as required
        do while(loop_count.lt.3)
!        do while(loop_count.lt.2)


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
!extra - to vary test
!       if( coord(1,l).eq..625.and.coord(2,l).eq..375.and.
!     .       coord(3,l).eq..375) refine(l)=.true.
       enddo
      elseif(ndim.eq.2) then
                do l=1,lnblocks
#ifdef KEVINS_2D_PATTERN
                if( (coord(1,l).gt..374.and.coord(1,l).lt..626) & 
     &           .and.(coord(2,l).gt..374.and.coord(2,l).lt..626) ) & 
     &                 refine(l)=.true.
#else
                if( coord(1,l).eq..125.and.coord(2,l).eq..125) & 
     &                 refine(l)=.true.
                if( coord(1,l).eq..375.and.coord(2,l).eq..375) & 
     &                 refine(l)=.true.
                if( coord(1,l).eq..625.and.coord(2,l).eq..875) & 
     &                 refine(l)=.true.
#endif /*  KEVINS_2D_PATTERN */
                enddo

      endif
      endif


! refine grid and apply morton reordering to grid blocks if necessary
#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'entering amr_refine_derefine : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */

      call amr_refine_derefine
#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited amr_refine_derefine : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */

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
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited mpi_morton_bnd_prolong : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

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

      else ! no_permanent_guardcells
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

      iopt = 1

       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

      if (nvar_work > 0) then
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
       do j=jlbnd,jubnd
       y0 = coord(2,l)-.5*(bsize(2,l)+dy)
       if(mod(nxb,2).eq.1) y0 = coord(2,l)-(.5*bsize(2,l)+dy)
       yj = y0 + dy*real(j-nguard)
       do i=ilbnd,iubnd
       x0 = coord(1,l)-.5*(bsize(1,l)+dx)
       if(mod(nxb,2).eq.1) x0 =  & 
     &       coord(1,l)-(.5*bsize(1,l)+dx)
       xi = x0 + dx*real(i-nguard)
       do ivar=1,nvar
       value = ax*(xi**interp_mask_unk(ivar))+ & 
     &         ay*(yj**interp_mask_unk(ivar))+ & 
     &         az*(zk**interp_mask_unk(ivar))
       v0 = value*real(ivar)
       if(abs(v0-unk1(ivar,i,j,k,1)) > accuracy) then
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
       value = ax*(xi**interp_mask_work(1))+ & 
     &         ay*(yj**interp_mask_work(1))+ & 
     &         az*(zk**interp_mask_work(1))
       if(abs(value-work1(i,j,k,1)) > accuracy) then
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
#ifndef NO_PERMANENT_GUARDCELLS
      unk_n1(:,:,:,:,1) = unk_n(:,:,:,:,l)
#endif
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
            do ivar=1,nvarcorn
              value = ax*xi**interp_mask_nc(ivar)+ & 
     &                ay*yj**interp_mask_nc(ivar)+ & 
     &                az*zk**interp_mask_nc(ivar)
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
#ifndef NO_PERMANENT_GUARDCELLS
      facevarx1(:,:,:,:,1) = facevarx(:,:,:,:,l)
      facevary1(:,:,:,:,1) = facevary(:,:,:,:,l)
      facevarz1(:,:,:,:,1) = facevarz(:,:,:,:,l)
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
            do ivar=1,nfacevar
              value = ax*xi**interp_mask_facex(ivar)+ & 
     &                ay*yj**interp_mask_facex(ivar)+ & 
     &                az*zk**interp_mask_facex(ivar)
              v0 = value*real(ivar)
              if(abs(facevarx1(ivar,i,j,k,1)-v0)>accuracy) then
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

      do k=klbnd,kubnd
        if(ndim.eq.3) z0 = coord(3,l)-.5*(bsize(3,l)+dz)
        zk = z0 + dz*real(k-nguard)
        do i=ilbnd,iubnd
          x0 = coord(1,l)-.5*(bsize(1,l)+dx)
          xi = x0 + dx*real(i-nguard)
          do j=jlbnd,jubnd+ione
            y0 = coord(2,l)-.5*bsize(2,l)-dy
            yj = y0 + dy*real(j-nguard)
            do ivar=1,nfacevar
              value = ax*xi**interp_mask_facey(ivar)+ & 
     &                ay*yj**interp_mask_facey(ivar)+ & 
     &                az*zk**interp_mask_facey(ivar)
              v0 = value*real(ivar)
              if(abs(facevary1(ivar,i,j,k,1)-v0)>accuracy) then
                      write(*,995) mype,l,ivar,i,j,k, & 
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
      if (ndim == 3) then
      do j=jlbnd,jubnd
        y0 = coord(2,l)-.5*(bsize(2,l)+dy)
        yj = y0 + dy*real(j-nguard)
        do i=ilbnd,iubnd
          x0 = coord(1,l)-.5*(bsize(1,l)+dx)
          xi = x0 + dx*real(i-nguard)
          do k=klbnd,kubnd+ione
            z0 = coord(3,l)-.5*bsize(3,l)-dz
            zk = z0 + dz*real(k-nguard)
            do ivar=1,nfacevar
              value = ax*xi**interp_mask_facez(ivar)+ & 
     &                ay*yj**interp_mask_facez(ivar)+ & 
     &                az*zk**interp_mask_facez(ivar)
              v0 = value*real(ivar)
              if(abs(facevarz1(ivar,i,j,k,1)-v0)>accuracy) then
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

      end if

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
#ifndef NO_PERMANENT_GUARDCELLS
      unk_e_x1(:,:,:,:,1) = unk_e_x(:,:,:,:,l)
      unk_e_y1(:,:,:,:,1) = unk_e_y(:,:,:,:,l)
      unk_e_z1(:,:,:,:,1) = unk_e_z(:,:,:,:,l)
#endif
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
            do ivar=1,nvaredge
              value = ax*xi**interp_mask_ec(ivar)+ & 
     &                ay*yj**interp_mask_ec(ivar)+ & 
     &                az*zk**interp_mask_ec(ivar)
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
            do ivar=1,nvaredge
              value = ax*xi**interp_mask_ec(ivar)+ & 
     &                ay*yj**interp_mask_ec(ivar)+ & 
     &                az*zk**interp_mask_ec(ivar)
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
            do ivar=1,nvaredge
              value = ax*xi**interp_mask_ec(ivar)+ & 
     &                ay*yj**interp_mask_ec(ivar)+ & 
     &                az*zk**interp_mask_ec(ivar)
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

      endif
      endif
      
      end if

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(mype.eq.0) write(*,*) 'End of automatic testing.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! Test consistency of face variables
      if(nfacevar.gt.0) then

      if(mype.eq.0) then
!         write(*,*) ' '
!         write(*,*) ' '
!         write(*,*) 'Testing routine for div B consistency check'
!         write(*,*) 'ax = ',ax
!         write(*,*) 'ay = ',ay
!         write(*,*) 'az = ',az
!         if(ndim.eq.2) then
!           write(*,*) 'div B should be = ax + ay = ',ax+ay
!         elseif(ndim.eq.3) then
!           write(*,*) 'div B should be = ax + ay + az = ',ax+ay+az
!         endif
      endif
      test_a = ax + ay + az*k3d      

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      istep = 0
!      call gtest_neigh_data(mype,istep,test_a)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      endif                       ! end of nfacevar iftest


! Sum number of errors 
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
     &                'Test of c_to_f_1blk_4 Successful ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        else
          write(*,*) ierror_tot,' errors detected - Test failed '
          write(55,*) ierror_tot,' errors detected - ', & 
     &                'Test of c_to_f_1blk_4 failed ', & 
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

      deallocate(tnewchild)
      deallocate(rflags)

      call report(mype,nprocs,ndim,no_permanent_guardcells, & 
     &            advance_all_levels,lmpi, & 
     &            ierror_tot,.true., & 
     &            'c_to_f_1blk_4                 ')


      call amr_close()

998      format('u:error proc ',i3,' block l= ',5(2x,i3),2x,f11.8,2x, & 
     &       f11.8)
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

      end
