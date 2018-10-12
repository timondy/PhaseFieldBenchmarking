!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_test_refinement(mype,lrefine_min,lrefine_max)

!--------------------------------------------------------------------------
! 
! This is a template to assist in constructing the routine AMR_TEST_REFINEMENT
! for use in your application. In this illustration we use the workspace
! array WORK to store the data which is used in computing the error measure
! at each grid point. This gives us the freedom to extend the testing
! beyond the normal bounds of individual blocksi, since WORK is declared
! with NGUARD_WORK guard cells at each boundary, which can be set to a
! larger number than NGUARD.

! Arguments:
!      mype integer local processor number
!      lrefine_min integer minimum refinement level to be permitted
!      lrefine_max integer maximum refinement level to be permitted

!--------------------------------------------------------------------------



      use physicaldata
      use tree
      use workspace
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"


      real error(ilw1:iuw1,jlw1:juw1,klw1:kuw1)

      integer pe,lb,icoord
      logical lcc,lfc,l_srl_only,ldiag

!
! Re-initialize the refinement and derefinement flag arrays
      refine(:)   = .false.
      derefine(:) = .false.
      error(:,:,:) = 0.


!
! Set up the workspace array WORK to stroe the variable we wish to examine 
! in order to test the refinement level.

#ifdef NO_PERMANENT_GUARDCELLS
#ifndef ADVANCE_ALL_LEVELS
      iopt = 1
      iempty = 0
      call amr_restrict(mype,iopt,iempty)
#endif
#endif

! Set up the workspace array to store the current solution.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
       do k=1,nzb
       do j=1,nyb
       do i=1,nxb
       work(i,j,k,lb,1) = unk(1,i,j,k,lb)
       end do
       end do
       end do
      endif
      end do
      endif

      call shmem_barrier_all()

!
! Set arguments for the call to amr_1blk_guardcell
      iopt=2
      nlayers=2


!
! Error limits which control the refinement and derefinement requests below.
      ctore = .35
      ctode = .05


!
! Loop over all leaf blocks and all parents of leaf blocks
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then


       pe = mype                   ! pe on which block to be treated is stored
       lcc = .true.                ! fill cell centered guardcell data
       lfc = .false.               ! do not fill cell-face-centered data 
       l_srl_only = .false.        ! do not limit to same refinement level
                                   ! neighbors
       icoord = 0                  ! apply to all block faces
       ldiag = .true.              ! include diagonal cells at edges and corners

       call amr_1blk_guardcell(mype,iopt,nlayers,lb,pe,lcc,lfc, & 
     &                              l_srl_only,icoord,ldiag)


!
! User provided routine which returns an array error, which has some error
! measure computed for each grid cell, based on some computation on the 
! input array WORK.
!      call error_measure(error,work)
       error(:,:,:) = 0.
       do k=klw1+k3d,kuw1-k3d
       do j=jlw1+1,juw1-1
       do i=ilw1+1,iuw1-1
         error1 = abs(work1(i+1,j,k,1)-work1(i,j,k,1))
         error2 = abs(work1(i-1,j,k,1)-work1(i,j,k,1))
         error3 = abs(work1(i,j+k2d,k,1)-work1(i,j,k,1))
         error4 = abs(work1(i,j-k2d,k,1)-work1(i,j,k,1))
         error5 = abs(work1(i,j,k+k3d,1)-work1(i,j,k,1))
         error6 = abs(work1(i,j,k-k3d,1)-work1(i,j,k,1))
         error_num = max( error1,error2,error3,error4, & 
     &                                  error5,error6 )
         error_den = max( work1(i,j,k,1)  ,work1(i+1,j,k,1), & 
     &                    work1(i-1,j,k,1),work1(i,j+k2d,k,1), & 
     &                    work1(i,j-k2d,k,1),work1(i,j,k+k3d,1), & 
     &                    work1(i,j,k-k3d,1), & 
     &                    1.0e-6 )

         error(i,j,k) = error_num/error_den
       enddo
       enddo
       enddo


       error_max = maxval( error )

! Does the error measure on this block anywhere exceed the limit which 
! should trigger refinement?

      if( lrefine(lb).lt.lrefine_max ) then
        if ( error_max . ge. ctore) refine(lb) = .true.
!        if ( error_max . ge. ctore) write(*,*) 
!     .       'refine lb lrefine ',refine(lb),lb,lrefine(lb),
!     .       error_max,maxloc(error)
      endif


! Can we derefine this block?
      if( lrefine(lb).gt.lrefine_min .and. (.not.refine(lb)) ) then
        if ( error_max . lt. ctode) derefine(lb) = .true.
      endif

#ifdef NOTNOW
           xtest = (bnd_box(1,1,lb)+2.05)*(bnd_box(2,1,lb)+2.05)
           ytest = (bnd_box(1,2,lb)+2.05)*(bnd_box(2,2,lb)+2.05)
           ztest = -1.
           if(ndim.eq.3) ztest = (bnd_box(1,3,lb)-.05) & 
     &                          *(bnd_box(2,3,lb)-.05)
           if(nodetype(lb).eq.1.and. & 
     &        xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
          write(*,*) 'block ',lb,nodetype(lb),bnd_box(1:2,1:ndim,lb)
      write(*,*) 'bef flux nsub ',nsub
             k = 1+nguard0*k3d
             do j=1,nyb+2*nguard_work
             write(*,50) j,(work1(i,j,k,1),i=1,nxb+2*nguard_work)
50           format(1x,i3,8(2x,f7.4))
             enddo
           endif
#endif

      endif
      end do                                   ! end of loop over blocks
      endif


      return
      end
