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


      subroutine amr_test_refinement(mype,llrefine_min,llrefine_max)



!--------------------------------------------------------------------------
! 
! This is a template to assist in constructing the routine AMR_TEST_REFINEMENT
! for use in your application. In this illustration we use the workspace
! array WORK to store the data which is used in computing the error measure
! at each grid point. This gives us the freedom to extend the testing
! beyond the normal bounds of individual blocks, since WORK is declared
! with NGUARD_WORK guard cells at each boundary, which can be set to a
! larger number than NGUARD.

! Arguments:
!   mype           integer      local processor number
!   llrefine_min   integer      minimum refinement level to be permitted
!   llrefine_max   integer      maximum refinement level to be permitted

!--------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces, only : amr_guardcell

      integer, intent(in)    ::  mype,llrefine_min,llrefine_max


      real :: error(ilw:iuw,jlw:juw,klw:kuw)

      integer :: ndel

      ndel = (nguard_work - nguard)*npgs

!-----------------------------------------------------------
!
! Re-initialize the refinement and derefinement flag arrays
      refine(:)   = .false.
      derefine(:) = .false.



!
! Set up the workspace array WORK to store the variable we wish to examine 
! in order to test the refinement level.


! Set up the workspace array to store the current solution.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
            do k=kl_bnd+ndel*k3d,ku_bnd+ndel*k3d
            do j=jl_bnd+ndel*k2d,ju_bnd+ndel*k2d
            do i=il_bnd+ndel,iu_bnd+ndel
              work(i,j,k,lb,1) =  & 
     &                 unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
            end do
            end do
            end do
      endif
      end do
      endif



!
! Fill the guard cell layers of the workspace array.
      iopt=2
      nlayers=1
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
! measure computed for each grid cell of the current grid block, based on 
! some computation on the input array WORK(:,:,:,lb,1).
       call error_measure( error, work(1,1,1,lb,1) )

       error_max = maxval( error )

! Does the error measure on this block anywhere exceed the limit which 
! should trigger refinement?

      if( lrefine(lb).lt.llrefine_max ) then
        if ( error_max .ge. ctore) refine(lb) = .true.
      endif


! Can we derefine this block?

      if( lrefine(lb).gt.llrefine_min .and. (.not.refine(lb)) ) then
        if ( error_max .lt. ctode) derefine(lb) = .true.
      endif



      endif
      end do                                   ! end of loop over blocks
      endif


!-----------------------------------------------------------

      return
      end subroutine amr_test_refinement
