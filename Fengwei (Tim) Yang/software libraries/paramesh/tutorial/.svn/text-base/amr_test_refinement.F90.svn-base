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


      real error(ilw:iuw,jlw:juw,klw:kuw)

!
! Re-initialize the refinement and derefinement flag arrays
      refine(:)   = .false.
      derefine(:) = .false.
      error(:,:,:) = 0.


!
! Set up the workspace array WORK to stroe the variable we wish to examine 
! in order to test the refinement level.

      ndel = nguard_work - nguard

! Set up the workspace array to store the current solution.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
       do k=1+nguard_work*k3d,nzb+nguard_work*k3d
       do j=1+nguard_work*k2d,nyb+nguard_work*k2d
       do i=1+nguard_work,nxb+nguard_work
       work(i,j,k,lb,1) = unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
       end do
       end do
       end do
      endif
      end do
      endif



!
! Fill the guard cell layers of the workspace array.
      iopt=2
      nlayers=2
      call amr_guardcell(mype,iopt,nlayers)


!
! Error limits which control the refinement and derefinement requests below.
      ctore = .35
      ctode = .05


!
! Loop over all leaf blocks and all parents of leaf blocks
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then


!
! User provided routine which returns an array error, which has some error
! measure computed for each grid cell, based on some computation on the 
! input array WORK.
!      call error_measure(error,work)
       error(:,:,:) = 0.
       do k=klw+k3d,kuw-k3d
       do j=jlw+k2d,juw-k2d
       do i=ilw+1  ,iuw-1
         error1 = abs(work(i+1,j,k,lb,1)-work(i,j,k,lb,1))
         error2 = abs(work(i-1,j,k,lb,1)-work(i,j,k,lb,1))
         error3 = abs(work(i,j+k2d,k,lb,1)-work(i,j,k,lb,1))
         error4 = abs(work(i,j-k2d,k,lb,1)-work(i,j,k,lb,1))
         error5 = abs(work(i,j,k+k3d,lb,1)-work(i,j,k,lb,1))
         error6 = abs(work(i,j,k-k3d,lb,1)-work(i,j,k,lb,1))
         error_num = max( error1,error2,error3,error4, & 
     &                                  error5,error6 )
         error_den = max( work(i,j,k,lb,1)  ,work(i+1,j,k,lb,1), & 
     &                    work(i-1,j,k,lb,1),work(i,j+k2d,k,lb,1), & 
     &                    work(i,j-k2d,k,lb,1),work(i,j,k+k3d,lb,1), & 
     &                    work(i,j,k-k3d,lb,1), & 
     &                     1.0e-6 )

         error(i,j,k) = error_num/error_den
       enddo
       enddo
       enddo


       error_max = maxval( error )

! Does the error measure on this block anywhere exceed the limit which 
! should trigger refinement?

      if( lrefine(lb).lt.lrefine_max ) then
        if ( error_max . ge. ctore) refine(lb) = .true.
      endif


! Can we derefine this block?
      if( lrefine(lb).gt.lrefine_min .and. (.not.refine(lb)) ) then
        if ( error_max . lt. ctode) derefine(lb) = .true.
      endif



      endif
      end do                                   ! end of loop over blocks
      endif


      return
      end
