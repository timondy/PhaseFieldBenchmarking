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

      program test_prolong




!
! This program generates an initial block, fills it with data
! which is a linear function of the local coordinate, uniformly
! refines the grid twice, and tests the prolongation operation.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use workspace
      use io
      use paramesh_comm_data

      use paramesh_interfaces, only : comm_start, & 
     &                                amr_initialize, & 
     &                                amr_refine_derefine, & 
     &                                amr_1blk_copy_soln, & 
     &                                amr_guardcell, & 
     &                                amr_1blk_guardcell, & 
     &                                amr_prolong, & 
     &                                amr_close

      use paramesh_mpi_interfaces, only :  & 
     &                                mpi_morton_bnd, & 
     &                                mpi_amr_comm_setup, & 
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
      logical, allocatable :: rflags(:)

!
! application specific variables

      real :: accuracy
      integer :: iopt,nlayers,icoord,pe
      integer :: ierror_sum,ierror_tot
      integer :: ierr, errcode
      logical lrefine_again,ltype2only
      logical lcc, lfc, lec, lnc, l_srl_only, ldiag

      logical :: lmpi,lnperm,ladvanceall

      logical :: lguard,lprolong,lflux,ledge,lrestrict
      logical :: lfulltree

      integer :: surrblks(2,3,3,3)
      integer :: psurrblks(2,3,3,3)

      save ierror_sum,ierror_tot

      character (len=80) :: filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      call amr_initialize

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      allocate(tnewchild(maxblocks_tr))
      allocate(rflags(maxblocks_tr))

      if (conserve) then
      write(*,*) 'CONSERVE must not be define dfor this test!'
      call amr_abort()
      end if

      if (nfacevar == 0) then
      write(*,*) 'ERROR for this test: NFACEVAR IS ZERO'
      go to 2
      end if

      if (nfield_divf <= 1) then 
       print *,' ERROR for this test: nfield_divf is too small, = ',nfield_divf
       go to 2
      end if

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)


      accuracy = 100./10.**precision(accuracy)

      print *,' nprocs = ',nprocs,mype
      print *,' nvar = ',nvar
      print *,' nfacevar = ',nfacevar
      print *,' nvaredge = ',nvaredge
      print *,' nvarcorn = ',nvarcorn

      if (diagonals) then
         write(*,*) 'diagonals on '
      end if

! Set up which field to clean

      i_divf_fc_vars(1:ndim,1) = 1
      i_divf_fc_vars(1:ndim,2) = 2

      ! Use divergence cleaning algorithm
      Call prol_fc_clean_divb_init(2, i_divf_fc_vars)
      print *,' Using divergence cleaning algorithm '


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

                neigh(:,:,1) = -210

                refine(1)=.true.
      endif


      boundary_index = -210
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
      
      

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!start test
! set the solution array to be a dipole field for the face variables

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
                        Call Dipole_field(xi, yj, zk, Ex, Ey, Ez)
                        facevarx(ivar,i,j,k,l)=Ex
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
                       Call Dipole_field(xi, yj, zk, Ex, Ey, Ez)
                      facevary(ivar,i,j,k,l)=Ey
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
                      Call Dipole_field(xi, yj, zk, Ex, Ey, Ez)
                      facevarz(ivar,i,j,k,l)=Ez
                    enddo
                  enddo
                enddo
              enddo

        endif

      enddo

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      iminx = nguard+1
      imaxx = nguard+nxb+1
      jminx = nguard*k2d+1
      jmaxx = nguard*k2d+nyb
      kminx = nguard*k3d+1
      kmaxx = nguard*k3d+nzb

      iminy = nguard+1
      imaxy = nguard+nxb
      jminy = nguard*k2d+1
      jmaxy = nguard*k2d+nyb+k2d
      kminy = nguard*k3d+1
      kmaxy = nguard*k3d+nzb

      iminz = nguard+1
      imaxz = nguard+nxb
      jminz = nguard*k2d+1
      jmaxz = nguard*k2d+nyb
      kminz = nguard*k3d+1
      kmaxz = nguard*k3d+nzb+k3d

!-----Make sure the top level field's are divergence free
      do ivar = 1, nfacevar
      Call clean_field   &
        (facevarx(ivar,iminx:imaxx,jminx:jmaxx,kminx:kmaxx,1), & 
         facevary(ivar,iminy:imaxy,jminy:jmaxy,kminy:kmaxy,1), &
         facevarz(ivar,iminz:imaxz,jminz:jmaxz,kminz:kmaxz,1), &
         iminx,imaxx,jminx,jmaxx,kminx,kmaxx,  &
         iminy,imaxy,jminy,jmaxy,kminy,kmaxy,  &
         iminz,imaxz,jminz,jmaxz,kminz,kmaxz,  &
         ndim, nxb, nyb, nzb, bsize(:,1))
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loop_count=0
! Now cycle over blocks adjusting refinement of initial setup as required
      do while(loop_count < 3)

! refine grid and apply morton reordering to grid blocks if necessary
      if (loop_count < 2) then
         derefine(:) = .false.
         stay(:) = .false.
         refine(:) = .true.
      else ! refine lower left corner
         derefine(:) = .false.
         stay(:) = .false.
         refine(:) = .false.
         do lb = 1, lnblocks

         if (ndim == 2) then
         ! lower left corner
         if (neigh(1,1,lb) <= -20 .and. neigh(1,3,lb) <= -20) then
            refine(lb) = .true.
         end if
         ! upper right corner
         if (neigh(1,2,lb) <= -20 .and. neigh(1,4,lb) <= -20) then
            refine(lb) = .true.
         end if
         ! block to the lower left of center
         if (bnd_box(2,1,lb) == .5 .and. bnd_box(2,2,lb) == .5) then
            refine(lb) = .true.
         end if

         else ! ndim == 3

         ! lower left corner
         if (neigh(1,1,lb) <= -20 .and. neigh(1,3,lb) <= -20 .and.   &
             neigh(1,5,lb) <= -20) then
            refine(lb) = .true.
         end if
         ! upper right corner
         if (neigh(1,2,lb) <= -20 .and. neigh(1,4,lb) <= -20 .and.   &
             neigh(1,6,lb) <= -20) then
            refine(lb) = .true.
         end if
         ! block to the lower left of center
         if (bnd_box(2,1,lb) == .5 .and. bnd_box(2,2,lb) == .5 .and. &
             bnd_box(2,3,lb) == .5) then
            refine(lb) = .true.
         end if

         end if ! End if (ndim == 2)

         end do
      end if

      call amr_refine_derefine()

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! a global prolongation call resets the newchild marker flags to false.
! Thus to test prolong for work we will need to restore this after the
! prolongation is applied to unk and facevar's
      tnewchild(:) = newchild(:)

      tag_offset = 100
      call mpi_morton_bnd_prolong & 
     &             (mype,nprocs,tag_offset)

       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited mpi_morton_bnd_prolong : pe ',mype
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

      iopt = 1
      nlayers = nguard
       write(*,*) 'entered amr_prolong : pe ',mype
       call amr_prolong(mype,iopt,nlayers)
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)
       write(*,*) 'exited amr_prolong : pe ',mype
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, errcode)

      loop_count=loop_count+1

      enddo
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!-----Compute the maximum divergence of the fields stored in facevarx, y z

      divb_max = -1.e30
      do lb = 1, lnblocks

         if (nodetype(lb) == 1) then

            if(ndim == 3) dz = bsize(3,lb)/real(nzb)
            dy = bsize(2,lb)/real(nyb)
            dx = bsize(1,lb)/real(nxb)

            do ivar = 1, nfacevar
            do k = nguard0+k3d, nguard0*k3d+nzb
               do j = nguard0+k2d, nguard0*k2d+nyb
                   do i = nguard0+1,nguard0+nxb
                      divb = (facevarx(ivar,i+1,j,k,lb) -       &
                              facevarx(ivar,i,j,k,lb))/dx
                      if (ndim >= 2) then
                      divb = divb +                             &
                             (facevary(ivar,i,j+k2d,k,lb) -     &
                              facevary(ivar,i,j,k,lb))/dy
                      end if
                      if (ndim == 3) then
                      divb = divb +                             &
                             (facevarz(ivar,i,j,k+k3d,lb) -     &
                              facevarz(ivar,i,j,k,lb))/dz
                      end if

                      divb_max = max(divb_max,divb)

                      unk(1,i,j,k,lb) = divb

                    end do
                 end do
              end do
              end do

           end if

        end do

        Call MPI_ALLREDUCE(divb_max, divb_max_all, 1, amr_mpi_real,  &
                           MPI_MAX, MPI_COMM_WORLD, ierr)

        print *,' MAX DIVB ACROSS ALL PROCS = ',divb_max_all

        call amr_plotfile_chombo(0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 2    if(mype.eq.0) write(*,*) 'End of automatic testing.'
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
     &                'Test of prolong_1blk Successful ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi & 
     &         ,' nxb nyb nzb ',nxb,nyb,nzb,' nguard ',nguard
        else
          write(*,*) ierror_tot,' errors detected - Test failed '
          write(55,*) ierror_tot,' errors detected - ', & 
     &                'Test of prolong_1blk failed ', & 
     &                ': nprocs ',nprocs,' ndim ',ndim, & 
     &                ' noperm ',lnperm,' advanceall ',ladvanceall, & 
     &                ' mpi ',lmpi & 
     &         ,' nxb nyb nzb ',nxb,nyb,nzb,' nguard ',nguard
        endif
      close(unit=55)
      endif

      call report(mype,nprocs,ndim,no_permanent_guardcells, & 
     &            advance_all_levels,lmpi, & 
     &            ierror_tot,.true., & 
     &            'prolong_1blk                  ')

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call amr_close()

         
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


      
