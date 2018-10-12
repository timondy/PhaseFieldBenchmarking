#define ALL_CELLS
!#define BLOCK_BND


!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

        subroutine output_tecplot_refmap(fileno)

!
! This routine outputs a uniformly refined paramesh grid in 2D as 
! multiple zones.
!


!
! Written:       Peter MacNeice
!                December 2004
!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata
      use workspace

! include file defining the tree
      use tree
      use io
      Use paramesh_comm_data

      use paramesh_interfaces, only : & 
     &                         amr_1blk_guardcell, & 
     &                         amr_perm_to_1blk, & 
     &                         amr_restrict, & 
     &                         amr_1blk_copy_soln, & 
     &                         amr_guardcell

      integer ::   fileno

      include 'mpif.h'

      integer :: nprocs,mype   ,jjjj
      logical :: lcc,lfc,lec,lnc
      logical :: ldiag,l_srl_only
      logical :: lguard,lprolong,lrestrict,lflux,ledge,lfulltree
      integer :: tag_offset
      real :: rdata(4,nxb+1,nyb+k2d,nzb+k3d)
      real :: sdata(4,nxb+1,nyb+k2d,nzb+k3d)

      real, allocatable :: x(:,:,:)
      real, allocatable :: y(:,:,:)
      real, allocatable :: z(:,:,:)
      integer,allocatable :: glnblocks(:)

      integer :: lb_max,ib1,local_ib

      integer,save :: idataout(1),idatain(1)

      real,parameter :: dummy_real = 0.

      integer cnodetype
      save mype,cnodetype

      character(len=6) filenumber

      integer :: leaf_lnblocks

      integer :: ierror, ierr
      integer :: lnodetype(maxblocks)
      integer,allocatable :: request_arr(:)
      integer :: status1(MPI_STATUS_SIZE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!--------------------------------------------------
! Make sure every pe knows how many leaf blocks every other pe has.

        if(.not.allocated(glnblocks)) allocate(glnblocks(0:nprocs-1))

        leaf_lnblocks = 0
        do lb = 1,lnblocks
          if(nodetype(lb).eq.1) leaf_lnblocks = leaf_lnblocks + 1
        enddo
        glnblocks(:) = 0
        call mpi_allgather(leaf_lnblocks,1,MPI_INTEGER, & 
     &                     glnblocks,1,MPI_INTEGER, & 
     &                     MPI_COMM_WORLD,ierror)
        lb_max = maxval(glnblocks)

        izone = 0

!--------------------------------------------------
! allocate temporary storage arrays
        allocate(x(nxb+1,nyb+k2d,nzb+k3d))
        allocate(y(nxb+1,nyb+k2d,nzb+k3d))
        allocate(z(nxb+1,nyb+k2d,nzb+k3d))
        allocate(request_arr(1:nprocs-1))

!----------------------------------------------------
! pe 0 writes file header

        if(mype.eq.0) then

        write(filenumber,"(I6.6)") fileno
        iout = 10
        open(unit=iout,status='unknown', & 
     &       file=trim(output_dir)//'tecplot_refmap.'//filenumber)
        write(*,*) ' file ','tecplot_refmap.'//filenumber

        write(iout,*) 'TITLE="Tecplot singular line test"'
        if (spherical_pm) then
        write(iout,191)
        else
        write(iout,191)
        end if
 191    format('VARIABLES="R" "THETA" "PHI" "Level"')
 192    format('VARIABLES="X" "Y" "Z" "Level"')

         endif                       ! end of mype if test
!----------------------------------------------------

         Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!----------------------------------------------------

         local_ib = 0
         do ib1 = 1,lb_max

!----------------------------------------------------

         Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


         if(lnblocks.gt.0) then

! get next leaf node id
          local_ib = local_ib+1
          do while ( (nodetype(local_ib).ne.1) .and. & 
     &                 (local_ib.le.lnblocks ))
            local_ib = local_ib+1
          enddo

! if we have not exhausted the list of blocks
          if(local_ib.le.lnblocks) then

            ip = mype
            ib = local_ib

! Now the basic solution data
            dy = 1.
            dz = 1.
            if(ndim.eq.3) dz =  & 
     &                       (bnd_box(2,3,ib)-bnd_box(1,3,ib))/real(nzb)
            if(ndim.ge.2) dy =  & 
     &                       (bnd_box(2,2,ib)-bnd_box(1,2,ib))/real(nyb)
            dx = (bnd_box(2,1,ib)-bnd_box(1,1,ib))/real(nxb)
           
            do k=1,nzb+k3d
            do j=1,nyb+k2d
            do i=1,nxb+1
              x(i,j,k) = bnd_box(1,1,ib)+dx*real(i-1)
              y(i,j,k) = bnd_box(1,2,ib)+dy*real(j-1)
              z(i,j,k) = bnd_box(1,3,ib)+dz*real(k-1)
            enddo
            enddo
            enddo

            do k=1,nzb+k3d
            do j=1,nyb+k2d
            do i=1,nxb+1
              sdata(1,i,j,k) = x(i,j,k)
              sdata(2,i,j,k) = y(i,j,k)
              sdata(3,i,j,k) = z(i,j,k)
!              sdata(4,i,j,k) = real(lrefine(ib))
            enddo
            enddo
            enddo
            sdata(4,:,:,:) = 0.
            do k=1,nzb
            do j=1,nyb
            do i=1,nxb
              sdata(4,i,j,k) = unk(i,j,k,1,ib)
            enddo
            enddo
            enddo

!------
! end compute for local leaf block 
!------------------------------------------------------
          endif                          ! end of local_ib if test  
         endif                          ! end of lnblocks if test  


        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


!------------------------------------------------------
! Output block ib1 from each processor


        do iproc = 0,nprocs-1

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! pe 0 will need to receive data from remote pes. Post
! all the necessary non-blocking receives here.

          if(mype.eq.0.and.iproc.gt.0) then
            if(ib1.le.glnblocks(iproc)) then
              call mpi_irecv(rdata(1,1,1,1), & 
     &                      4*(nxb+1)*(nyb+k2d)*(nzb+k3d), & 
     &                      amr_mpi_real, & 
     &                      iproc,1+iproc*10, & 
     &                      MPI_COMM_WORLD,request_arr(iproc), & 
     &                       ierror)
            endif
          endif

        if(mype.gt.0) then
          if(mype.eq.iproc) then
            if(ib1.le.leaf_lnblocks) then
            call mpi_send(sdata(1,1,1,1), & 
     &                    4*(nxb+1)*(nyb+k2d)*(nzb+k3d), & 
     &                    amr_mpi_real, & 
     &                    0,1+mype*10, & 
     &                    MPI_COMM_WORLD,ierror)
            endif
          endif
        endif
        

        if(mype.eq.0) then
            izone = izone+1
            if(izone.eq.1) then
#ifdef ALL_CELLS
              write(iout,193) (nxb+1)*(nyb+k2d)*(nzb+k3d),nxb*nyb*nzb
#endif /* ALL_CELLS */
#ifdef BLOCK_BND
              write(iout,193) (1+1)*(1+k2d)*(1+k3d),1*1*1
#endif /* BLOCK_BND */
193         format('ZONE N=',i4,',E=',i4, & 
     &                ',F=FEPOINT, ET=BRICK')
!     .                ',F=FEPOINT, ET=QUADRILATERAL')
            elseif(izone.gt.1) then
#ifdef ALL_CELLS
              write(iout,194) (nxb+1)*(nyb+k2d)*(nzb+k3d),nxb*nyb*nzb
#endif /* ALL_CELLS */
#ifdef BLOCK_BND
              write(iout,194) (1+1)*(1+k2d)*(1+k3d),1*1*1
#endif /* BLOCK_BND */
194           format('ZONE N=',i4,',E=',i4, & 
     &    ',F=FEPOINT, ET=BRICK, D=(FECONNECT)')
!     .    ',F=FEPOINT, ET=QUADRILATERAL, D=(FECONNECT)')
            endif

            if(ib1.le.glnblocks(iproc)) then
            if(iproc.gt.0) then
              call mpi_wait(request_arr(iproc),status1,ierror)
            else
              rdata = sdata
            endif
            endif

!---        
! output this block
#ifdef ALL_CELLS
            do k=1,nzb+k3d
            do j=1,nyb+k2d
            do i=1,nxb+1
#endif /* ALL_CELLS */
#ifdef BLOCK_BND
            do k=1,nzb+k3d,1+(nzb-1)*k3d
            do j=1,nyb+k2d,nyb*k2d
            do i=1,nxb+1,nxb
#endif /* BLOCK_BND */
              write(iout,50) rdata(1:4,i,j,k)
            enddo
            enddo
            enddo
!---        
! Output node lists for each quadrilateral
            if(izone.eq.1) then
#ifdef ALL_CELLS
            nz0 = 1
            ny0 = 1
            nx0 = 1
            do ke=1,nzb
            do je=1,nyb
            do ie=1,nxb
#endif /* ALL_CELLS */
#ifdef BLOCK_BND
            nz0 = nzb
            ny0 = nyb
            nx0 = nxb
            do ke=1,nzb,nzb-k3d
            do je=1,nyb,nyb-k2d
            do ie=1,nxb,nxb-1
#endif /* BLOCK_BND */
                nz1 = nzb/nz0
                ny1 = nyb/ny0
                nx1 = nxb/nx0
                n1 = ie+(je-1)*(nx1+1)+(ke-1)*(nx1+1)*(ny1+k2d)
                n2 = n1+1
                n3 = n2+(nx1+1)
                n4 = n3-1
                n5 = ie+(je-1)*(nx1+1)+ke*(nx1+1)*(ny1+k2d)
                n6 = n5+1
                n7 = n6+(nx1+1)
                n8 = n7-1
                write(iout,55) n1,n2,n3,n4,n5,n6,n7,n8
              enddo
              enddo
              enddo
            endif
!---        

        endif                          ! end of mype if test


        enddo                        ! end of iproc do


        enddo                          ! end of ib1 do loop

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!------------------------------------------------------

50      format(4(1x,1pe13.6))
55      format(8(1x,i6))

        if(mype.eq.0) close(unit=iout)

        deallocate(x)
        deallocate(y)
        deallocate(z)
        deallocate(glnblocks)
        deallocate(request_arr)

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      return
      end subroutine output_tecplot_refmap

