       subroutine user_coord_transfm(lb,pe)
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
!-----------------------------------------------------------------------
! This routine applies any desired coordinate transformation to the
! underlying, uniformly spaced grid, for the grid block lb,pe.
!
! Written:      Peter MacNeice                  November 2002
!
!
!-----------------------------------------------------------------------
       use physicaldata
       use tree

       implicit none

       include 'mpif.h'

       integer, intent(in) :: lb,pe

       integer :: mype, ierr
       integer :: i,j,ibxl,ibxr,ibyl,ibyr,ibzl,ibzr
       real,save     :: bbox(2,3)
       integer,save  :: cneigh(2,6)


!-----------------------------------------------------------------------

! local variables

!-----------------------------------------------------------------------

          Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)

          cneigh(:,:) = neigh(:,:,lb)
          if(pe.eq.mype) then
            bbox(:,:) = bnd_box(:,:,lb)
            cneigh(1,1) = surr_blks(1,1,1+k2d,1+k3d,lb)
            cneigh(1,2) = surr_blks(1,3,1+k2d,1+k3d,lb)
            cneigh(1,3) = surr_blks(1,2,1,1+k3d,lb)
            cneigh(1,4) = surr_blks(1,2,1+2*k2d,1+k3d,lb)
            ibxl = cneigh(1,1)
            ibxr = cneigh(1,2)
            if(ndim.ge.2) then
              ibyl = cneigh(1,3)
              ibyr = cneigh(1,4)
            endif
            if(ndim.eq.3) then
              ibzl = cneigh(1,5)
              ibzr = cneigh(1,6)
            endif
          else
            lb0 = -1
            do llb = strt_buffer,last_buffer
              if(laddress(1,llb).eq.lb.and.laddress(1,llb).eq.pe) & 
     &          lb0 = llb
            enddo
            if(lb0.eq.-1) then
              write(*,*) 'User coord transfm: error : blk not found'
              call amr_abort()
            endif
            bbox(:,:) = bnd_box(:,:,lb0)
          call mpi_amr_get_bc_settings(lb,pe,mype, & 
     &                        ibxl,ibxr,ibyl,ibyr,ibzl,ibzr)
          endif


!------------------------------------------------------
! Capture transformation in amr curvilinear coordinates
          do i = il_bnd1,iu_bnd1+1
            cell_face_coord1(i) = ???
          enddo

          if(ndim.ge.2) then
          do j = jl_bnd1,ju_bnd1+k2d
            cell_face_coord2(j) = ???
          enddo
          endif

          if(ndim.eq.3) then
          do k = kl_bnd1,ku_bnd1+k3d
            cell_face_coord3(k) = ???
          enddo
          endif
!------------------------------------------------------


      return
      end subroutine user_coord_transfm
