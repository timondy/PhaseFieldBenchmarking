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


      subroutine gtest_neigh_data(mype,istep,test_a)





!------------------------------------------------------------------------
!
! This routine compares face centered data on neighboring blocks at
! all relevant points in space. If there is no refinement jump
! across the block boundary the field value on the chosen cell face
! for both blocks are printed. If there is a refinement jump then the
! 4 face values on the finer block are summed, and the sum is printed
! along with the corresponding face value on the coarse neighbor.

! When interpreting the log file entries produced when errors are
! detected, bear in mind that this routine works by looping over all 
! faces of leaf blocks,
! 1. for each leaf block it tests to see if faces 2, 4 and 6 have
!    a neighbor at their refinement level. If they do, then it
!    computes the divergence for each cell immediately inside each
!    of these iblock boundary faces, first using local data, then using 
!    data for the block boundary acquired from the neighbor. It reports
!    a divergence error if either of these tests exceeds a preset
!    threshold value. Faces 1, 3 and 5 are not tested since these
!    will be handled as faces 2, 4 and 6 of the appropriate neighbor
!    blocks.
! 2. wherever a refinement jump is detected, the divergence test is
!    performed by computing the divergence for each 2x2x2 block
!    of cells immediately interior to the face of the finer block,
!    and then repeating the test with a data value from the coarser
!    neighbor block corresponding to the block boundary data used 
!    in the first computation. If either exceed the preset threshold
!    an error message is inserted in the log file identifying the
!    lowest i,j,k indeces of the 2x2x2 block in question.
!    (Obviously in 2D this would be a 2x2 block).
!
!
!
! Arguments:
!      mype             local processor
!
!
! Written :     Peter MacNeice          August 1998
! Modified:     Peter MacNeice          January 2002
! 
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use paramesh_comm_data

      use paramesh_interfaces, only :  comm_int_sum_to_all
      use paramesh_mpi_interfaces, only :   & 
     &                                 mpi_amr_comm_setup, & 
     &                                 mpi_amr_get_remote_block_fvar

      implicit none

      include 'mpif.h'

#ifdef TIMINGS
#include "timer.fh"
#endif

      integer, intent(in)    ::  mype,istep
      real,    intent(in)    ::  test_a


      integer :: nguard0

!-------------------------

      integer remote_pe,remote_block

! local arrays
      real recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd)
      real recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &       kl_bnd:ku_bnd)
      real recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd+k3d)


      integer,save :: anodetype(1)
      integer   cnodetype
      save cnodetype

      integer par_neigh(2,mfaces)
      integer parent_blk,parent_pe
      save par_neigh

      character*6 pattern,pattern1
      integer     pattern_buf(5000)
      integer     lb_buf(5000),iprint,jprint,kprint,level_buf(5000)
      real        data_buf(4,5000)
      integer     ijk_buf(3,5000)
      save        pattern_buf,lb_buf,data_buf,level_buf
      save        iprint,jprint,kprint,ijk_buf


      logical :: lguard,lprolong,lflux,ledge,lrestrict
      logical :: lfulltree
      logical :: lcc,lfc,lec,lnc
      integer :: tag_offset,iblk,nprocs,icoord
      integer :: ii,jj,kk,ii0,jj0,kk0,ii1,jj1,kk1,iproc,ip
      integer :: i0,j0,k0,i1,j1,k1,i2,j2,k2
      integer :: ioff,joff,koff,jchild,idest,iface
      integer :: neigh_blk,neigh_pe,iface_max,lb,iopt
      integer :: ndel

      real :: bxsum,bysum,bzsum,eps,delxi,delyi,delzi,dx,dy,dz
      real :: bxsum1,bysum1,bzsum1
      real :: divbm,divbl,divbnl,divbnlmax,divbmax
      integer,dimension (:), allocatable ::  g_iprint
      integer ::  l_iprint(1)
      integer :: status(MPI_STATUS_SIZE)
      integer :: jsrc,jdest,itag,isize,ierr
      real, parameter :: divb_threshold = 5.e-10
 
!-------------------------

      nguard0 = nguard*npgs

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      divbmax = 0.
      divbnlmax = 0.

      iprint = 0
      jprint = 0
      kprint = 0
      lb_buf(:) = 0
      level_buf(:) = 0
      pattern_buf(:) = 0
      data_buf(:,:) = 0.

! Prepare communication needs.
! Call comm_setup for guardcell filling
      iopt = 1
      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(iopt.eq.1) then
        if(nvar.gt.0) lcc = .true.
        if(nfacevar.gt.0) lfc = .true.
        if(nvaredge.gt.0) lec = .true.
        if(nvarcorn.gt.0) lnc = .true.
      elseif(iopt.ge.2) then
        lcc = .true.
      endif

      call amr_1blk_copy_soln(-1)

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


! Loop over blocks.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then

        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        delxi = 1./dx
        delyi = 1./dy
        delzi = 1.
        dz = 0.
        if (ndim == 3) then
        dz = bsize(3,lb)/real(nzb)
        delzi = 1./dz
        end if
        eps = .01*dx
        eps = .5*dx

        if (curvilinear) then
           call amr_block_geometry(lb,mype)
        endif

!----------------
! Neighbors at same refinement level


! Loop over one face for each dimension
        if (ndim == 3) then
        iface_max = 6
        else
        iface_max = 4
        end if
      do iface = 2,iface_max,2

         neigh_blk = neigh(1,iface,lb)
         neigh_pe  = neigh(2,iface,lb)

          remote_block = neigh_blk
          remote_pe    = neigh_pe
! if (neigh_blk,neigh_pe) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.
          if(neigh_pe.ne.mype) then
            do iblk = strt_buffer,last_buffer
              if(neigh_blk.eq.laddress(1,iblk).and. & 
     &           neigh_pe .eq.laddress(2,iblk) ) then
                neigh_blk = iblk
                neigh_pe  = mype
              endif
            enddo
          endif

          idest = 1
          lcc = .false.
          lfc = .true.
          lec = .false.
          lnc = .false.

! Is the neighbor at the same refinement level?
      if(neigh_blk.gt.0) then

          cnodetype = nodetype(neigh_blk)

! Only consider neighbor if it is a leaf node.
          if(cnodetype.eq.1) then


! Copy complete remote block into a buffer block called recv1.
! Only consider faces 2, 4 and 6, to avoid redundant check.
         if(iface.eq.2) then

           icoord = 1
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)

           i0 = 1+nxb+nguard0
           i1 = nguard0+1

           do k0 = 1+nguard0*k3d,nzb+nguard0*k3d
           do j0 = 1+nguard0,nyb+nguard0

              j1 = j0
              k1 = k0

              if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
              endif

              if(ndim.eq.3) then
           divbl = & 
     & delxi*(facevarx(1,i0,j0,k0,lb)    -facevarx(1,i0-1,j0,k0,lb)) & 
     &+delyi*(facevary(1,i0-1,j0+1,k0,lb)-facevary(1,i0-1,j0,k0,lb)) & 
     &+delzi*(facevarz(1,i0-1,j0,k0+1,lb)-facevarz(1,i0-1,j0,k0,lb))
              else
           divbl = & 
     & delxi*(facevarx(1,i0,j0,k0,lb)    -facevarx(1,i0-1,j0,k0,lb)) & 
     &+delyi*(facevary(1,i0-1,j0+1,k0,lb)-facevary(1,i0-1,j0,k0,lb))
              endif


              if(ndim.eq.3) then
           divbnl = & 
     & delxi*(recvx(1,i1,j1,k1)          -facevarx(1,i0-1,j0,k0,lb)) & 
     &+delyi*(facevary(1,i0-1,j0+1,k0,lb)-facevary(1,i0-1,j0,k0,lb)) & 
     &+delzi*(facevarz(1,i0-1,j0,k0+1,lb)-facevarz(1,i0-1,j0,k0,lb))
              else
           divbnl = & 
     & delxi*(recvx(1,i1,j1,k1)          -facevarx(1,i0-1,j0,k0,lb)) & 
     &+delyi*(facevary(1,i0-1,j0+1,k0,lb)-facevary(1,i0-1,j0,k0,lb))
              endif

           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then

           iprint = iprint+1
           pattern_buf(iprint) = 1
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = facevarx(1,i0,j0,k0,lb)
           data_buf(2,iprint) = recvx(1,i1,j1,k1)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0

           endif

           enddo
           enddo

         elseif(iface.eq.4) then
           icoord = 2
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)

           j0 = 1+nyb+nguard0
           j1 = nguard0+1

           do k0 = 1+nguard0*k3d,nzb+nguard0*k3d
           do i0 = 1+nguard0,nxb+nguard0

              i1 = i0
              k1 = k0

              if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
              endif

              if(ndim.eq.3) then
           divbl = & 
     & delxi*(facevarx(1,i0+1,j0-1,k0,lb)-facevarx(1,i0,j0-1,k0,lb)) & 
     &+delyi*(facevary(1,i0,j0,k0,lb)    -facevary(1,i0,j0-1,k0,lb)) & 
     &+delzi*(facevarz(1,i0,j0-1,k0+1,lb)-facevarz(1,i0,j0-1,k0,lb))
              else
           divbl = & 
     & delxi*(facevarx(1,i0+1,j0-1,k0,lb)-facevarx(1,i0,j0-1,k0,lb)) & 
     &+delyi*(facevary(1,i0,j0,k0,lb)    -facevary(1,i0,j0-1,k0,lb))
              endif


              if(ndim.eq.3) then
           divbnl = & 
     & delxi*(facevarx(1,i0+1,j0-1,k0,lb)-facevarx(1,i0,j0-1,k0,lb)) & 
     &+delyi*(recvy(1,i1,j1,k1)          -facevary(1,i0,j0-1,k0,lb)) & 
     &+delzi*(facevarz(1,i0,j0-1,k0+1,lb)-facevarz(1,i0,j0-1,k0,lb))
              else
           divbnl =  & 
     & delxi*(facevarx(1,i0+1,j0-1,k0,lb)-facevarx(1,i0,j0-1,k0,lb)) & 
     &+delyi*(recvy(1,i1,j1,k1)          -facevary(1,i0,j0-1,k0,lb))
              endif

           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then

           iprint = iprint+1
           pattern_buf(iprint) = 2
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = facevary(1,i0,j0,k0,lb)
           data_buf(2,iprint) = recvy(1,i1,j1,k1)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0

           endif

           enddo
           enddo

         elseif(iface.eq.6.and.ndim==3) then

           icoord = 3
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)

           k0 = 1+nzb+nguard0
           k1 = nguard0+1


           do j0 = 1+nguard0,nyb+nguard0
           do i0 = 1+nguard0,nxb+nguard0

              i1 = i0
              j1 = j0

              if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
              endif

           divbl = & 
     & delxi*(facevarx(1,i0+1,j0,k0-1,lb)-facevarx(1,i0,j0,k0-1,lb)) & 
     &+delyi*(facevary(1,i0,j0+1,k0-1,lb)-facevary(1,i0,j0,k0-1,lb)) & 
     &+delzi*(facevarz(1,i0,j0,k0,lb)    -facevarz(1,i0,j0,k0-1,lb))


           divbnl = & 
     & delxi*(facevarx(1,i0+1,j0,k0-1,lb)-facevarx(1,i0,j0,k0-1,lb)) & 
     &+delyi*(facevary(1,i0,j0+1,k0-1,lb)-facevary(1,i0,j0,k0-1,lb)) & 
     &+delzi*(recvz(1,i1,j1,k1)          -facevarz(1,i0,j0,k0-1,lb))


           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then

           iprint = iprint+1
           pattern_buf(iprint) = 3
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = facevarz(1,i0,j0,k0,lb)
           data_buf(2,iprint) = recvz(1,i1,j1,k1)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0

           endif

           enddo
           enddo

         endif                          ! end of iface if test


         endif                          ! end of neigh leaf node test

        endif                           ! end of test of neigh refinement level

       enddo                            ! end of loop over block faces


!----------------
! Now deal with neighbors at coarser level


! Get parents neighbor and child info
        parent_blk = parent(1,lb)
        parent_pe  = parent(2,lb)
! if (parent_blk,parent_pe) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.
          if(parent_pe.ne.mype) then
            do iblk = strt_buffer,last_buffer
              if(parent_blk.eq.laddress(1,iblk).and. & 
     &           parent_pe .eq.laddress(2,iblk) ) then
                parent_blk = iblk
                parent_pe  = mype
              endif
            enddo
          endif

         par_neigh(:,:) = neigh(:,:,parent_blk)

! Proceed only if a parent exists
        if(parent_blk.gt.0) then

         jchild = which_child(lb)


! compute the offset in the parent block appropriate for this child
         ioff = mod(jchild-1,2)*nxb/2
         joff = mod((jchild-1)/2,2)*nyb/2
         koff = mod((jchild-1)/4,2)*nzb/2



! Loop over block faces
        if (ndim == 3) then
        iface_max = 6
        else
        iface_max = 4
        end if
        do iface = 1,iface_max


! compute the offset in the parent block appropriate for this child
! (do this inside the iface loop since these offsets may be adjusted)
         ioff = mod(jchild-1,2)*nxb/2
         joff = mod((jchild-1)/2,2)*nyb/2
         koff = mod((jchild-1)/4,2)*nzb/2

! Select the correct neighbor
        neigh_blk = par_neigh(1,iface)
        neigh_pe  = par_neigh(2,iface)

! if (neigh_blk,neigh_pe) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.
          remote_block = neigh_blk
          remote_pe    = neigh_pe
          if(neigh_pe.ne.mype) then
            do iblk = strt_buffer,last_buffer
              if(neigh_blk.eq.laddress(1,iblk).and. & 
     &           neigh_pe .eq.laddress(2,iblk) ) then
                neigh_blk = iblk
                neigh_pe  = mype
              endif
            enddo
          endif

          idest = 1
          lcc = .false.
          lfc = .true.
          lec = .false.
          lnc = .false.

! Is this face next to a coarser block ( and not an external boundary)?
!        if(neigh_blk.eq.-1) then
        if(neigh(1,iface,lb).eq.-1) then


! flip the offset in the appropriate axis
         if(iface.le.2) then
           ioff = mod(jchild,2)*nxb/2
         elseif(iface.eq.3.or.iface.eq.4) then
           joff = mod((jchild)/2,2)*nyb/2
         elseif(iface.gt.4) then
           koff = mod((jchild)/4,2)*nzb/2
         endif



! Copy complete remote block into a buffer block called recv1.
         if(iface.eq.1) then

           icoord = 1
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)
            i0 = nguard0+1
            i2 = nxb+nguard0+1

            do k0 = 1+nguard0*k3d,nzb+nguard0*k3d,2
            do j0 = 1+nguard0,nyb+nguard0*k2d,2

            k2 = ((k0-nguard0-1)/2+nguard0+koff)*k3d+1
            j2 = (j0-nguard0-1)/2+1+nguard0+joff

            j1 = j0 + k2d
            k1 = k0 + k3d

            if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
           endif

            bxsum = facevarx(1,i0,j0,k0,lb)+facevarx(1,i0,j1,k0,lb) & 
     &             +facevarx(1,i0,j0,k1,lb)+facevarx(1,i0,j1,k1,lb)

            bxsum1=facevarx(1,i0+2,j0,k0,lb)+facevarx(1,i0+2,j1,k0,lb) & 
     &            +facevarx(1,i0+2,j0,k1,lb)+facevarx(1,i0+2,j1,k1,lb)
            ii = i0
            ii1 = ii+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            kk0 = min(k0,k1)
            kk1 = kk0+k3d
            bysum =facevary(1,ii,jj0,kk0,lb)+facevary(1,ii1,jj0,kk0,lb) & 
     &            +facevary(1,ii,jj0,kk1,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1= & 
     &         facevary(1,ii,jj0+2,kk0,lb)+facevary(1,ii1,jj0+2,kk0,lb) & 
     &        +facevary(1,ii,jj0+2,kk1,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            if (ndim == 3) then
            bzsum =facevarz(1,ii,jj0,kk0,lb)+facevarz(1,ii1,jj0,kk0,lb) & 
     &            +facevarz(1,ii,jj1,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1= & 
     &         facevarz(1,ii,jj0,kk0+2,lb)+facevarz(1,ii1,jj0,kk0+2,lb) & 
     &        +facevarz(1,ii,jj1,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)
            else
            bzsum = 0.
            bzsum1= 0.
            endif

            divbl  = ( (bxsum1 - bxsum)*delxi  & 
     &               + (bysum1 - bysum)*delyi  & 
     &               + (bzsum1 - bzsum)*delzi )*.125
            divbnl = ( (bxsum1 - 4.*recvx(1,i2,j2,k2))*delxi  & 
     &               + (bysum1 - bysum)*delyi  & 
     &               + (bzsum1 - bzsum)*delzi )*.125

           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))


           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then
           iprint = iprint+1
           pattern_buf(iprint) = 4
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bxsum
           data_buf(2,iprint) = recvx(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0
           endif

           enddo
           enddo

         elseif(iface.eq.2) then

           icoord = 1
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)
            i0 = nxb+nguard0+1
            i2 = nguard0+1

            do k0 = 1+nguard0*k3d,nzb+nguard0*k3d,2
            do j0 = 1+nguard0,nyb+nguard0*k2d,2

            k2 = ((k0-nguard0-1)/2+nguard0+koff)*k3d+1
            j2 = (j0-nguard0-1)/2+1+nguard0+joff

            j1 = j0 + k2d
            k1 = k0 + k3d

            if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
           endif

            bxsum = facevarx(1,i0,j0,k0,lb)+facevarx(1,i0,j1,k0,lb) & 
     &             +facevarx(1,i0,j0,k1,lb)+facevarx(1,i0,j1,k1,lb)

            bxsum1=facevarx(1,i0-2,j0,k0,lb)+facevarx(1,i0-2,j1,k0,lb) & 
     &            +facevarx(1,i0-2,j0,k1,lb)+facevarx(1,i0-2,j1,k1,lb)
            ii = i0-2
            ii1 = ii+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            kk0 = min(k0,k1)
            kk1 = kk0+k3d
            bysum =facevary(1,ii,jj0,kk0,lb)+facevary(1,ii1,jj0,kk0,lb) & 
     &            +facevary(1,ii,jj0,kk1,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1= & 
     &         facevary(1,ii,jj0+2,kk0,lb)+facevary(1,ii1,jj0+2,kk0,lb) & 
     &        +facevary(1,ii,jj0+2,kk1,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            if (ndim == 3) then
            bzsum =facevarz(1,ii,jj0,kk0,lb)+facevarz(1,ii1,jj0,kk0,lb) & 
     &            +facevarz(1,ii,jj1,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1= & 
     &         facevarz(1,ii,jj0,kk0+2,lb)+facevarz(1,ii1,jj0,kk0+2,lb) & 
     &        +facevarz(1,ii,jj1,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)
            else
            bzsum = 0.
            bzsum1= 0.
            endif

            divbl  = ( (bxsum  - bxsum1)*delxi & 
     &               + (bysum1 - bysum)*delyi & 
     &               + (bzsum1 - bzsum)*delzi )*.125
            divbnl = ( (4.*recvx(1,i2,j2,k2) - bxsum1)*delxi  & 
     &               + (bysum1 - bysum)*delyi & 
     &               + (bzsum1 - bzsum)*delzi )*.125

           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then
           iprint = iprint+1
           pattern_buf(iprint) = 4
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bxsum
           data_buf(2,iprint) = recvx(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0
           endif

           enddo
           enddo

         elseif(iface.eq.3) then

           icoord = 2
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)
            j0 = nguard0+1
            j2 = nyb+nguard0+1

            do k0 = 1+nguard0*k3d,nzb+nguard0*k3d,2
            do i0 = 1+nguard0,nxb+nguard0,2

            k2 = ((k0-nguard0-1)/2+nguard0+koff)*k3d + 1
            i2 = (i0-nguard0-1)/2+1+nguard0+ioff

            i1 = i0 + 1
            k1 = k0 + k3d


            if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
           endif

            bysum = facevary(1,i0,j0,k0,lb)+facevary(1,i1,j0,k0,lb) & 
     &        +facevary(1,i0,j0,k1,lb)+facevary(1,i1,j0,k1,lb)

            bysum1= facevary(1,i0,j0+2,k0,lb)+facevary(1,i1,j0+2,k0,lb) & 
     &        +facevary(1,i0,j0+2,k1,lb)+facevary(1,i1,j0+2,k1,lb)
            jj = j0
            jj1 = jj+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            kk0 = min(k0,k1)
            kk1 = kk0+k3d
            bxsum =facevarx(1,ii0,jj,kk0,lb)+facevarx(1,ii0,jj1,kk0,lb) & 
     &            +facevarx(1,ii0,jj,kk1,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1= & 
     &         facevarx(1,ii0+2,jj,kk0,lb)+facevarx(1,ii0+2,jj1,kk0,lb) & 
     &        +facevarx(1,ii0+2,jj,kk1,lb)+facevarx(1,ii0+2,jj1,kk1,lb)
            if (ndim == 3) then
            bzsum =facevarz(1,ii0,jj,kk0,lb)+facevarz(1,ii0,jj1,kk0,lb) & 
     &            +facevarz(1,ii1,jj,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1= & 
     &         facevarz(1,ii0,jj,kk0+2,lb)+facevarz(1,ii0,jj1,kk0+2,lb) & 
     &        +facevarz(1,ii1,jj,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)
            else
            bzsum = 0.
            bzsum1= 0.
            end if


            divbl  = ( (bxsum1 - bxsum)*delxi & 
     &               + (bysum1 - bysum)*delyi & 
     &               + (bzsum1 - bzsum)*delzi )*.125
            divbnl = ( (bxsum1 - bxsum)*delxi & 
     &               + (bysum1 - 4.*recvy(1,i2,j2,k2))*delyi  & 
     &               + (bzsum1 - bzsum)*delzi )*.125


           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then
           iprint = iprint+1
           pattern_buf(iprint) = 5
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bysum
           data_buf(2,iprint) = recvy(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0
           endif

           enddo
           enddo

         elseif(iface.eq.4) then
           icoord = 2
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)

            j0 = nyb+nguard0+1
            j2 = nguard0+1

            do k0 = 1+nguard0*k3d,nzb+nguard0*k3d,2
            do i0 = 1+nguard0,nxb+nguard0,2

            k2 = ((k0-nguard0-1)/2+nguard0+koff)*k3d+1
            i2 = (i0-nguard0-1)/2+1+nguard0+ioff

            i1 = i0 + 1
            k1 = k0 + k3d

            if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
           endif

            bysum = facevary(1,i0,j0,k0,lb)+facevary(1,i1,j0,k0,lb) & 
     &             +facevary(1,i0,j0,k1,lb)+facevary(1,i1,j0,k1,lb)

            bysum1= facevary(1,i0,j0-2,k0,lb)+facevary(1,i1,j0-2,k0,lb) & 
     &             +facevary(1,i0,j0-2,k1,lb)+facevary(1,i1,j0-2,k1,lb)
            jj = j0-2
            jj1 = jj+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            kk0 = min(k0,k1)
            kk1 = kk0+k3d
            bxsum =facevarx(1,ii0,jj,kk0,lb)+facevarx(1,ii0,jj1,kk0,lb) & 
     &            +facevarx(1,ii0,jj,kk1,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1= & 
     &         facevarx(1,ii0+2,jj,kk0,lb)+facevarx(1,ii0+2,jj1,kk0,lb) & 
     &        +facevarx(1,ii0+2,jj,kk1,lb)+facevarx(1,ii0+2,jj1,kk1,lb)
            if (ndim == 3) then
            bzsum =facevarz(1,ii0,jj,kk0,lb)+facevarz(1,ii0,jj1,kk0,lb) & 
     &            +facevarz(1,ii1,jj,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1= & 
     &         facevarz(1,ii0,jj,kk0+2,lb)+facevarz(1,ii0,jj1,kk0+2,lb) & 
     &        +facevarz(1,ii1,jj,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)
            else
            bzsum = 0.
            bzsum1= 0.
            end if


            divbl  = ( (bxsum1 - bxsum )*delxi & 
     &               + (bysum  - bysum1)*delyi & 
     &               + (bzsum1 - bzsum )*delzi )*.125
            divbnl = ( (bxsum1 - bxsum )*delxi & 
     &               + (4.*recvy(1,i2,j2,k2) -bysum1)*delyi  & 
     &               + (bzsum1 - bzsum )*delzi )*.125


           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then
           iprint = iprint+1
           pattern_buf(iprint) = 5
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bysum
           data_buf(2,iprint) = recvy(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0
           endif

           enddo
           enddo

         elseif(iface.eq.5.and.ndim==3) then
           icoord = 3
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)
            k0 = nguard0+1
            k2 = nzb+nguard0+1

            do j0 = 1+nguard0,nyb+nguard0*k2d,2
            do i0 = 1+nguard0,nxb+nguard0,2

            i2 = (i0-nguard0-1)/2+1+nguard0+ioff
            j2 = (j0-nguard0-1)/2+1+nguard0+joff

            i1 = i0 + 1
            j1 = j0 + k2d

            if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
           endif

            bzsum = facevarz(1,i0,j0,k0,lb)+facevarz(1,i1,j0,k0,lb) & 
     &             +facevarz(1,i0,j1,k0,lb)+facevarz(1,i1,j1,k0,lb)

            bzsum1= facevarz(1,i0,j0,k0+2,lb)+facevarz(1,i1,j0,k0+2,lb) & 
     &             +facevarz(1,i0,j1,k0+2,lb)+facevarz(1,i1,j1,k0+2,lb)
            kk = k0
            kk1 = kk+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            bysum =facevary(1,ii0,jj0,kk,lb)+facevary(1,ii0,jj0,kk1,lb) & 
     &            +facevary(1,ii1,jj0,kk,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1= & 
     &         facevary(1,ii0,jj0+2,kk,lb)+facevary(1,ii0,jj0+2,kk1,lb) & 
     &        +facevary(1,ii1,jj0+2,kk,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            bxsum =facevarx(1,ii0,jj0,kk,lb)+facevarx(1,ii0,jj0,kk1,lb) & 
     &            +facevarx(1,ii0,jj1,kk,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1= & 
     &         facevarx(1,ii0+2,jj0,kk,lb)+facevarx(1,ii0+2,jj0,kk1,lb) & 
     &        +facevarx(1,ii0+2,jj1,kk,lb)+facevarx(1,ii0+2,jj1,kk1,lb)

            divbl  = ( (bxsum1 - bxsum)*delxi & 
     &               + (bysum1 - bysum)*delyi  & 
     &               + (bzsum1 - bzsum)*delzi )*.125
            divbnl = ( (bxsum1 - bxsum)*delxi & 
     &               + (bysum1 - bysum)*delyi  & 
     &               + (bzsum1 - 4.*recvz(1,i2,j2,k2))*delzi )*.125

           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then
           iprint = iprint+1
           pattern_buf(iprint) = 6
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bzsum
           data_buf(2,iprint) = recvz(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0
           endif

           enddo
           enddo

         elseif(iface.eq.6.and.ndim==3) then
           icoord = 3
           idest  = 1
           call mpi_amr_get_remote_block_fvar(mype, & 
     &                                        remote_pe,remote_block, & 
     &                                        icoord,recvx,recvy,recvz, & 
     &                                        idest)

            k0 = nzb+nguard0+1
            k2 = nguard0+1

            do j0 = 1+nguard0,nyb+nguard0*k2d,2
            do i0 = 1+nguard0,nxb+nguard0,2

            i2 = (i0-nguard0-1)/2+1+nguard0+ioff
            j2 = (j0-nguard0-1)/2+1+nguard0+joff

            i1 = i0 + 1
            j1 = j0 + k2d

            if (curvilinear) then
              ndel = nguard - nguard0
              delxi = 1./cell_leng1(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delyi = 1./cell_leng2(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              delzi = 1.
              if(ndim.eq.3) then
              delzi = 1./cell_leng3(i0+ndel,j0+ndel*k2d,k0+ndel*k3d)
              endif
           endif

            bzsum = facevarz(1,i0,j0,k0,lb)+facevarz(1,i1,j0,k0,lb) & 
     &             +facevarz(1,i0,j1,k0,lb)+facevarz(1,i1,j1,k0,lb)

            bzsum1= facevarz(1,i0,j0,k0-2,lb)+facevarz(1,i1,j0,k0-2,lb) & 
     &             +facevarz(1,i0,j1,k0-2,lb)+facevarz(1,i1,j1,k0-2,lb)

            kk = k0-2
            kk1 = kk+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            bysum =facevary(1,ii0,jj0,kk,lb)+facevary(1,ii0,jj0,kk1,lb) & 
     &            +facevary(1,ii1,jj0,kk,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1= & 
     &         facevary(1,ii0,jj0+2,kk,lb)+facevary(1,ii0,jj0+2,kk1,lb) & 
     &        +facevary(1,ii1,jj0+2,kk,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            bxsum =facevarx(1,ii0,jj0,kk,lb)+facevarx(1,ii0,jj0,kk1,lb) & 
     &            +facevarx(1,ii0,jj1,kk,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1= & 
     &         facevarx(1,ii0+2,jj0,kk,lb)+facevarx(1,ii0+2,jj0,kk1,lb) & 
     &        +facevarx(1,ii0+2,jj1,kk,lb)+facevarx(1,ii0+2,jj1,kk1,lb)

            divbl  = ( (bxsum1 - bxsum)*delxi  & 
     &               + (bysum1 - bysum)*delyi  & 
     &               + (bzsum - bzsum1)*delzi )*.125
            divbnl = ( (bxsum1 - bxsum)*delxi  & 
     &               + (bysum1 - bysum)*delyi & 
     &               + (4.*recvz(1,i2,j2,k2) - bzsum1)*delzi )*.125

           divbm = max(abs(divbl),abs(divbnl))

           divbmax = max(divbmax,abs(divbl))
           divbnlmax = max(divbnlmax,abs(divbnl))

           if(test_a.ge.0.5) then
             divbm = max(abs(divbl-test_a),abs(divbnl-test_a))
           endif
           if(divbm.gt.divb_threshold) then
           iprint = iprint+1
           pattern_buf(iprint) = 6
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bzsum
           data_buf(2,iprint) = recvz(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl
           ijk_buf(1,iprint) = i0
           ijk_buf(2,iprint) = j0
           ijk_buf(3,iprint) = k0
           endif

           enddo
           enddo

         endif


      endif                             ! end of neigh type iftest


      enddo                             ! loop over block faces


      endif                             ! end of parent iftest

!-------------


       endif
       enddo                            ! end of loop over blocks
       endif



      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!
! Is there any output to collect?
      jprint = 0
      call comm_int_sum_to_all(jprint,iprint)

      if(mype.eq.0) write(*,*) 'Total of ',jprint, & 
     &              ' divergence errors detected'

      if(mype.eq.0.and.jprint.gt.0) then

        open(unit=22,status='unknown',position='append', & 
     &                                  file='gfacevar.dbg')
        write(22,*) 'Total of ',jprint, & 
     &              ' divergence errors detected', & 
     &              ' test_a = ',test_a

      endif


         l_iprint(1) = iprint
         if(allocated(g_iprint)) deallocate(g_iprint)
         allocate( g_iprint(nprocs) )
         Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

         call mpi_gather(l_iprint,1,MPI_INTEGER,g_iprint,1, & 
     &                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      do iproc=0,nprocs-1
        kprint = 0

            jsrc = iproc
            jdest= 0

            kprint = g_iprint(iproc+1)
            if(mype.eq.0.and.iproc.gt.0) then
              if(kprint.gt.0) then

                itag = jsrc*nprocs*10 + 1
                isize = kprint
                call Mpi_int_recv(lb_buf,isize,MPI_INTEGER, & 
     &                        jsrc,itag,MPI_COMM_WORLD,status,ierr)

                itag = jsrc*nprocs*10 + 2
                isize = kprint
                call Mpi_int_recv(level_buf,isize,MPI_INTEGER, & 
     &                        jsrc,itag,MPI_COMM_WORLD,status,ierr)

                itag = jsrc*nprocs*10 + 3
                isize = kprint
                call Mpi_int_recv(pattern_buf,isize,MPI_INTEGER, & 
     &                        jsrc,itag,MPI_COMM_WORLD,status,ierr)

                itag = jsrc*nprocs*10 + 4
                isize = kprint*4
                call Mpi_real_recv(data_buf,isize,amr_mpi_real, & 
     &                        jsrc,itag,MPI_COMM_WORLD,status,ierr)

                itag = jsrc*nprocs*10 + 5
                isize = kprint*3
                call Mpi_int_recv(ijk_buf,isize,MPI_INTEGER, & 
     &                        jsrc,itag,MPI_COMM_WORLD,status,ierr)

              endif                          ! end of kprint if
            endif                            ! end of mype if

            if(mype.gt.0.and.mype.eq.iproc.and.iprint.gt.0) then

              itag = jsrc*nprocs*10 + 1
              isize = iprint
              call Mpi_int_send(lb_buf,isize,MPI_INTEGER, & 
     &                      jdest,itag,MPI_COMM_WORLD,ierr)

              itag = jsrc*nprocs*10 + 2
              isize = iprint
              call Mpi_int_send(level_buf,isize,MPI_INTEGER, & 
     &                      jdest,itag,MPI_COMM_WORLD,ierr)

              itag = jsrc*nprocs*10 + 3
              isize = iprint
              call Mpi_int_send(pattern_buf,isize,MPI_INTEGER, & 
     &                      jdest,itag,MPI_COMM_WORLD,ierr)

              itag = jsrc*nprocs*10 + 4
              isize = iprint*4
              call Mpi_real_send(data_buf,isize,amr_mpi_real, & 
     &                      jdest,itag,MPI_COMM_WORLD,ierr)

              itag = jsrc*nprocs*10 + 5
              isize = iprint*3
              call Mpi_int_send(ijk_buf,isize,MPI_INTEGER, & 
     &                      jdest,itag,MPI_COMM_WORLD,ierr)

            endif                         ! end of mype+iproc+iprint if

         Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            if(mype.eq.0) then
            do ip = 1,kprint
              if(pattern_buf(ip).eq.1) then
                pattern = 'srl x '
              elseif(pattern_buf(ip).eq.2) then
                pattern = 'srl y '
              elseif(pattern_buf(ip).eq.3) then
                pattern = 'srl z '
              elseif(pattern_buf(ip).eq.4) then
                pattern = 'drl x '
              elseif(pattern_buf(ip).eq.5) then
                pattern = 'drl y '
              elseif(pattern_buf(ip).eq.6) then
                pattern = 'drl z '
              endif
              pattern1 = 'div B '
              write(22,700) pattern,istep,iproc,lb_buf(ip), & 
     &             ijk_buf(1,ip),ijk_buf(2,ip),ijk_buf(3,ip), & 
     &             data_buf(1,ip),data_buf(2,ip),level_buf(ip)
              write(22,700) pattern1,istep,iproc,lb_buf(ip), & 
     &             ijk_buf(1,ip),ijk_buf(2,ip),ijk_buf(3,ip), & 
     &             data_buf(3,ip),data_buf(4,ip),level_buf(ip)
            enddo
            endif
700      format(a6,i5,1x,i4,1x,i4,1x,3(1x,i4),2(1x,1pe12.5),1x,i4)

         Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        enddo                                ! end of loop over iproc

        if(allocated(g_iprint)) deallocate(g_iprint)

      if(mype.eq.0.and.jprint.gt.0) close(unit=22)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      return
      end subroutine gtest_neigh_data
