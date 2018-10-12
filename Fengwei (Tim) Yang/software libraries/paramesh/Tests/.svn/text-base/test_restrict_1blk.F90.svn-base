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


      program test_restrict



! This program tests restriction. It sets up the same refinement
! pattern as in test_c_to_f, then zeroes out parent blocks of
! leaf modes before refilling them with a restriction operation.
! Only these leaf node parents are then tested.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use workspace
      use io
      use mpi_morton

      use paramesh_interfaces, only : comm_start, & 
     &                                amr_initialize, & 
     &                                amr_refine_derefine, & 
     &                                amr_1blk_copy_soln, & 
     &                                amr_guardcell, & 
     &                                amr_prolong, & 
     &                                amr_restrict, & 
     &                                amr_1blk_guardcell, & 
     &                                amr_1blk_guardcell_reset, & 
     &                                guardcell_test, & 
     &                                mesh_test, & 
     &                                amr_close

      use paramesh_mpi_interfaces, only :  & 
     &                                mpi_morton_bnd, & 
     &                                mpi_amr_comm_setup, & 
     &                                mpi_morton_bnd_prolong, & 
     &                                mpi_morton_bnd_restrict


! Only required for programs in ./Tests
#include "test_defs.fh"

      include 'mpif.h'

      integer :: tag_offset,max_blks_sent
      integer nguard0
      integer nguard_work0
      integer ndel
      real vdefault
      parameter(vdefault=-1.e5)
      integer :: three,four,five, six, lb

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
      integer jchild,cnodetype,remote_blk,remote_pe
      save    jchild,cnodetype


      integer iopt,nlayers,icoord
      integer ierror_sum,ierror_tot
      logical lrefine_again,ltype2only
      logical lcc, lfc, lec, lnc, l_srl_only, ldiag

      logical :: lmpi,lnperm,ladvanceall
      integer :: errorcode, ierr

      logical :: lguard,lprolong,lflux,ledge,lrestrict
      logical :: lfulltree


      save ierror_sum,ierror_tot

      character (len=80) :: filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call amr_initialize 

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs
      ndel = nguard_work0-nguard0

      allocate(tnewchild(maxblocks_tr))
      allocate(rflags(maxblocks_tr))

#ifdef CONSERVE
      write(*,*) 'CONSERVE must not be defined for this test!'
      call amr_abort()
#endif

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)


      accuracy = 1000./10.**precision(accuracy)


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

      endif                     ! advance_all_levels

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
                      do ivar=1,nfacevar
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
                    do ivar=1,nfacevar
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
                    do ivar=1,nfacevar
                      value = ax*xi+ay*yj+az*zk
                      facevarz(ivar,i,j,k,l)=value*real(ivar)
                    enddo
                  enddo
                enddo
              enddo


           endif                ! advance_all_levels

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
 
        endif

      enddo

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loop_count=0
! Now cycle over blocks adjusting refinement of initial setup as required
        do while(loop_count.lt.3)

        write(*,*) 'proc start loop_count ',mype,loop_count

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

      endif
      endif


! refine grid and apply morton reordering to grid blocks if necessary
      call amr_refine_derefine
        write(*,*) 'proc end refder loop_count ',mype,loop_count

! a global prolongation call resets the newchild marker flags to false.
! Thus to test prolong for work we will need to restore this after the
! prolongation is applied to unk and facevar's
      tnewchild(:) = newchild(:)

      tag_offset = 100
      call mpi_morton_bnd_prolong & 
     &             (mype,nprocs,tag_offset)
        write(*,*) 'proc end mortpr loop_count ',mype,loop_count

      iopt = 1
      nlayers = nguard
      call amr_prolong(mype,iopt,nlayers)
        write(*,*) 'proc end pr loop_count ',mype,loop_count
#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited amr_prolong : pe ',mype, & 
     &    ' loop_count = ',loop_count
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif /* DEBUG */

      if (no_permanent_guardcells) then
! Store a copy of the current solution in gt_unk
      call amr_1blk_copy_soln(-1)
      end if

      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)
        write(*,*) 'proc end mortbn loop_count ',mype,loop_count

      iopt = 1


      if (.not.no_permanent_guardcells) then


! set guard cell data to zero to ensure proper test of guardcell
! set external guard cell data.
!        call zero_guardcells(ioptw)


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

      else 

      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(nvar.gt.0) lcc = .true.
      if(nfacevar.gt.0) lfc = .true.
      if(nvaredge.gt.0) lec = .true.
      if(nvarcorn.gt.0) lnc = .true.
      tag_offset = 100
      lguard    = .false.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .true.
      lfulltree = .false.
        write(*,*) 'proc ent comm loop_count ',mype,loop_count
      call mpi_amr_comm_setup(mype,nprocs, & 
     &                        lguard,lprolong,lflux,ledge, & 
     &                        lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)

      endif

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

        if (nvar > 0 .and. nvar_work > 0) then

! Set up work arrays for testing 
        do lb = 1,lnblocks
        do i=1,nvar_work
          work(1+nguard0+ndel:nxb+nguard0+ndel, & 
     &         1+(nguard0+ndel)*k2d:nyb+(nguard0+ndel)*k2d, & 
     &         1+(nguard0+ndel)*k3d:nzb+(nguard0+ndel)*k3d,lb,i) = & 
     &       unk(1,1+nguard0:nxb+nguard0, & 
     &             1+nguard0*k2d:nyb+nguard0*k2d, & 
     &             1+nguard0*k3d:nzb+nguard0*k3d,lb)
        enddo
        enddo

        end if

! Zero out leaf block parents in preparation for test of restriction
        do l=1,lnblocks
        if(nodetype(l).eq.2) then
         if (nvar > 0) then
            unk(:,:,:,:,l) = 0.
         end if
         if (nvar_work > 0) then
            work(:,:,:,l,ioptw-1) = 0.
         end if
         if (nfacevar > 0) then
            facevarx(:,:,:,:,l) = 0.
            facevary(:,:,:,:,l) = 0.
            facevarz(:,:,:,:,l) = 0.
         end if
         if (nvarcorn > 0) then
            unk_n(:,:,:,:,l) = 0.
         end if
         if (nvaredge > 0) then
            unk_e_x(:,:,:,:,l) = 0.
            unk_e_y(:,:,:,:,l) = 0.
            unk_e_z(:,:,:,:,l) = 0.
         end if
        endif
        enddo

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


        tag_offset = 100
        lfulltree = .false.
        call mpi_morton_bnd_restrict & 
     &             (mype,nprocs,tag_offset)
        tag_offset = 100

! If any children are themselves parents, then mark their volume
! with a default data value
        do l=1,lnblocks
        if(nodetype(l).eq.2) then

        do ich = 1,nchild
         remote_blk = child(1,ich,l)
         remote_pe  = child(2,ich,l)
         if(remote_pe .ne. mype) then
           do iblk = strt_buffer,last_buffer
             if(remote_blk.eq.laddress(1,iblk) .and. & 
     &          remote_pe .eq.laddress(2,iblk) ) then
               remote_blk = iblk
             endif
           enddo
         endif
         cnodetype = nodetype(remote_blk)
         if(cnodetype.ne.1) then
         jchild = ich
! compute the offset in the parent block appropriate for this child
         ioff = mod(jchild-1,2)*nxb/2
         joff = mod((jchild-1)/2,2)*nyb/2
         koff = mod((jchild-1)/4,2)*nzb/2

         nxh = nxb - nxb/2
         nyh = nyb - nyb/2
         nzh = nzb - nzb/2

         write(*,*) 'pe ',mype,' blk ',l,' is a parent with leaf ' & 
     &     ,'siblings but child ',jchild,' is itself a parent ', & 
     &     ': mark it with default in index range ', & 
     &      1+nguard0+ioff,nguard0+ioff+nxh, & 
     &         1+(nguard0+joff)*k2d,1+(nguard0+joff+nyh-1)*k2d, & 
     &         1+(nguard0+koff)*k3d,1+(nguard0+koff+nzh-1)*k3d

         if (nvar > 0) then
         unk(:,1+nguard0+ioff:nguard0+ioff+nxh, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh-1)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh-1)*k3d, & 
     &         l) = 2.*vdefault
         end if

         if (nvarcorn > 0) then
         unk_n(:,1+nguard0+ioff:nguard0+ioff+nxh+1, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh)*k3d, & 
     &         l) = 2.*vdefault
         end if

         if (nvaredge > 0) then
         unk_e_x(:,1+nguard0+ioff:nguard0+ioff+nxh, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh)*k3d, & 
     &         l) = 2.*vdefault

         unk_e_y(:,1+nguard0+ioff:nguard0+ioff+nxh+1, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh-1)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh)*k3d, & 
     &         l) = 2.*vdefault

         unk_e_z(:,1+nguard0+ioff:nguard0+ioff+nxh+1, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh-1)*k3d, & 
     &         l) = 2.*vdefault
         end if

         if (nvar_work > 0) then
         work(1+nguard_work0+ioff:nguard_work0+ioff+nxh, & 
     &   1+(nguard_work0+joff)*k2d:1+(nguard_work0+joff+nyh-1)*k2d, & 
     &   1+(nguard_work0+koff)*k3d:1+(nguard_work0+koff+nzh-1)*k3d, & 
     &   l,ioptw-1) = 2.*vdefault
         end if

         if (nfacevar > 0) then
         facevarx(:,1+nguard0+ioff:nguard0+ioff+nxh+1, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh-1)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh-1)*k3d, & 
     &         l) = 2.*vdefault
         facevary(:,1+nguard0+ioff:nguard0+ioff+nxh, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh-1)*k3d, & 
     &         l) = 2.*vdefault
         facevarz(:,1+nguard0+ioff:nguard0+ioff+nxh, & 
     &         1+(nguard0+joff)*k2d:1+(nguard0+joff+nyh-1)*k2d, & 
     &         1+(nguard0+koff)*k3d:1+(nguard0+koff+nzh)*k3d, & 
     &         l) = 2.*vdefault
         end if

        endif
        enddo

        endif
        enddo

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        
! Restrict
        iopt = 1
        write(*,*) 'main : pe ',mype,' calling amr_restrict for ', & 
     &    ' iopt = ',iopt
        iempty = 0
        call amr_restrict(mype,iopt,iempty,.false.)

        write(*,*) 'main : pe ',mype,' exited amr_restrict for ', & 
     &    ' iopt = ',iopt


        if(mod(nxb,2).eq.0) then
! If using an odd grid size then restriction on WORK will have
! errors near external boundaries, because during the restriction
! operation guardcell is called to fill the guardcells of the
! working block array WORK1(...,1), and unlike UNK and FACEVARs,
! there is no boundary condition routine available to fill guardcells
! of WORK arrays on external boundaries. If using an even sized grid
! the default restriction operation does not require any guardcells.
        if (nvar_work > 0) then
        iopt = ioptw
        write(*,*) 'main : pe ',mype,' calling amr_restrict for ', & 
     &    ' iopt = ',iopt
        iempty = 0
        call amr_restrict(mype,iopt,iempty,.false.)
        write(*,*) 'main : pe ',mype,' exited amr_restrict for ', & 
     &    ' iopt = ',iopt
        endif
        end if

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      if(mype.eq.0) write(*,*) 'Start of automatic testing.'


      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test of unk communications

      if (nvar > 0) then

      do ii=0,nprocs-1
      if(mype.eq.ii) then

      do l=1,lnblocks
      if(nodetype(l).eq.2) then

      if (no_permanent_guardcells) then
      unk1(:,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &       1+nguard*k3d:nzb+nguard*k3d,1) =  & 
     &  unk(:,1:nxb,1:nyb,1:nzb,l)
      else                      ! no_permanent_guardcells
#ifndef NO_PERMANENT_GUARDCELLS
      unk1(:,:,:,:,1) = unk(:,:,:,:,l)
#endif
      endif                     ! no_permanent_guardcells


      ilbnd=1+nguard
      iubnd=nxb+nguard
      jlbnd=1+nguard*k2d
      jubnd=nyb+nguard*k2d
      klbnd=1+nguard*k3d
      kubnd=nzb+nguard*k3d

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
       value = ax*xi+ay*yj+az*zk
       do ivar=1,nvar
       v0 = value*real(ivar)
       if(abs(v0-unk1(ivar,i,j,k,1))>accuracy & 
     &          .and.unk1(ivar,i,j,k,1).gt.vdefault) then
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

      if (mod(nxb,2).eq.0) then


      do ii=0,nprocs-1
      if(mype.eq.ii) then

      do l=1,lnblocks
      if(nodetype(l).eq.2) then

      if (no_permanent_guardcells) then
      work1(1+nguard_work:nxb+nguard_work, & 
     &      1+nguard_work*k2d:nyb+nguard_work*k2d, & 
     &      1+nguard_work*k3d:nzb+nguard_work*k3d,1) =  & 
     &  work(1:nxb,1:nyb,1:nzb,l,ioptw-1)
      else
#ifndef NO_PERMANENT_GUARDCELLS
      work1(:,:,:,1) = work(:,:,:,l,ioptw-1)
#endif
      endif   ! no_permanent_guardcells



      ilbnd=1+nguard_work
      iubnd=nxb+nguard_work
      jlbnd=1+nguard_work*k2d
      jubnd=nyb+nguard_work*k2d
      klbnd=1+nguard_work*k3d
      kubnd=nzb+nguard_work*k3d

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
       if(abs(value-work1(i,j,k,1))>accuracy & 
     &          .and.work1(i,j,k,1).gt.vdefault) then
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

      end if  ! end if (nvar_work > 0)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(mype.eq.0) write(*,*) 'work test complete.'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      else
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if(mype.eq.0) write(*,*)  & 
     &           'work test skipped because inappropriate.'
       Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif

      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test of unk_n communications

      if (nvarcorn > 0) then

      iopt = 1
      nlayers = nguard
 
      do ii=0,nprocs-1
      if(mype.eq.ii) then
 
      do l=1,lnblocks
      if(nodetype(l).eq.2) then
 
      if (no_permanent_guardcells) then
      unk_n1(:,1+nguard:nxb+nguard+1,1+nguard*k2d:nyb+(nguard+1)*k2d, & 
     &       1+nguard*k3d:nzb+(nguard+1)*k3d,1) = & 
     &  unk_n(:,1:nxb+1,1:nyb+k2d,1:nzb+k3d,l)
      else
#ifndef NO_PERMANENT_GUARDCELLS
      unk_n1(:,:,:,:,1) = unk_n(:,:,:,:,l)
#endif
      end if !no_permanent_guardcells
 
 
      ilbnd=1+nguard
      iubnd=nxb+nguard
      jlbnd=1+nguard*k2d
      jubnd=nyb+nguard*k2d
      klbnd=1+nguard*k3d
      kubnd=nzb+nguard*k3d
 
 
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
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end if
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef FACE

      if (nfacevar > 0) then

! test of guardcell for facevar's (if appropriate)
      if(mod(nxb,2).eq.0) then


! test values
        do ii=0,nprocs-1
        if(mype.eq.ii) then

      do l=1,lnblocks
      if(nodetype(l).eq.2) then

      if (no_permanent_guardcells) then
      facevarx1(:,1+nguard:nxb+nguard+1,1+nguard*k2d:nyb+nguard*k2d, & 
     &       1+nguard*k3d:nzb+nguard*k3d,1) =  & 
     &  facevarx(:,1:nxb+1,1:nyb,1:nzb,l)
      facevary1(:,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d+k2d, & 
     &       1+nguard*k3d:nzb+nguard*k3d,1) =  & 
     &  facevary(:,1:nxb,1:nyb+k2d,1:nzb,l)
      facevarz1(:,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d, & 
     &       1+nguard*k3d:nzb+nguard*k3d+k3d,1) =  & 
     &  facevarz(:,1:nxb,1:nyb,1:nzb+k3d,l)
      else
#ifndef NO_PERMANENT_GUARDCELLS
      facevarx1(:,:,:,:,1) = facevarx(:,:,:,:,l)
      facevary1(:,:,:,:,1) = facevary(:,:,:,:,l)
      facevarz1(:,:,:,:,1) = facevarz(:,:,:,:,l)
#endif
      endif ! no_permanent_guardcells


        ilbnd=1+nguard
        iubnd=nxb+nguard
        jlbnd=1+nguard*k2d
        jubnd=nyb+nguard*k2d
        klbnd=1+nguard*k3d
        kubnd=nzb+nguard*k3d

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
              if(abs(facevarx1(ivar,i,j,k,1)-v0)>accuracy & 
     &          .and.facevarx1(ivar,i,j,k,1).gt.vdefault) then
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
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevary1(ivar,i,j,k,1)-v0)>accuracy & 
     &          .and.facevary1(ivar,i,j,k,1).gt.vdefault) then
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
            value = ax*xi+ay*yj+az*zk
            do ivar=1,nfacevar
              v0 = value*real(ivar)
              if(abs(facevarz1(ivar,i,j,k,1)-v0)>accuracy & 
     &          .and.facevarz1(ivar,i,j,k,1).gt.vdefault) then
!                     write(*,994) mype,l,ivar,i,j,k,
                     write(*,*) mype,l,ivar,i,j,k, & 
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
 
        do ii=0,nprocs-1
        if(mype.eq.ii) then
 
        do l=1,lnblocks
        if(nodetype(l).eq.2) then

        if (no_permanent_guardcells) then
        unk_e_x1(:,1+nguard:nxb+nguard, & 
     &             1+nguard*k2d:nyb+(nguard+1)*k2d, & 
     &             1+nguard*k3d:nzb+(nguard+1)*k3d,1) = & 
     &    unk_e_x(:,1:nxb,1:nyb+k2d,1:nzb+k3d,l)
        unk_e_y1(:,1+nguard:nxb+nguard+1, & 
     &             1+nguard*k2d:nyb+nguard*k2d, & 
     &             1+nguard*k3d:nzb+(nguard+1)*k3d,1) = & 
     &    unk_e_y(:,1:nxb+1,1:nyb,1:nzb+k3d,l)
        unk_e_z1(:,1+nguard:nxb+nguard+1, & 
     &             1+nguard*k2d:nyb+(nguard+1)*k2d, & 
     &             1+nguard*k3d:nzb+nguard*k3d,1) = & 
     &    unk_e_z(:,1:nxb+1,1:nyb+k2d,1:nzb,l)
      else
#ifndef NO_PERMANENT_GUARDCELLS
      unk_e_x1(:,:,:,:,1) = unk_e_x(:,:,:,:,l)
      unk_e_y1(:,:,:,:,1) = unk_e_y(:,:,:,:,l)
      unk_e_z1(:,:,:,:,1) = unk_e_z(:,:,:,:,l)
#endif
      endif ! no_permanent_guardcells


        ilbnd=il_bnd1+nguard
        iubnd=iu_bnd1-nguard
        jlbnd=jl_bnd1+nguard*k2d
        jubnd=ju_bnd1-nguard*k2d
        klbnd=kl_bnd1+nguard*k3d
        kubnd=ku_bnd1-nguard*k3d
 
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
      call amr_1blk_guardcell_reset
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)         

      endif
      end if

      end if

#endif
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
      filename = trim(output_dir) // 'test.log'
      open(unit = 55,file=filename, & 
     &            status='unknown',position='append')
        write(*,*) ' '
        write(*,*) ' '
        if(ierror_tot.eq.0) then
          write(*,*) 'No errors detected - Test Successful '
          write(55,*) 'No errors detected - ', & 
     &                'Test of restrict_1blk Successful ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi
        else
          write(*,*) ierror_tot,' errors detected - Test failed '
          write(55,*) ierror_tot,' errors detected - ', & 
     &                'Test of restrict_1blk failed ', & 
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
     &            'restrict_1blk                 ')


2      call amr_close()

         
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
      deallocate(rflags)

      end
