!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_restrict
!! NAME
!!
!!   amr_restrict
!! 
!! SYNOPSIS
!!
!!   call amr_restrict (mype, iopt, iempty)
!!   call amr_restrict (mype, iopt, iempty, filling_guardcells)
!!
!!   call amr_restrict (integer, integer, integer, optional logical)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: mype 
!!     Current processor number.
!!
!!   integer, intent(in) :: iopt  
!!     Switch to select which datastructures
!!     are updated. If iopt=1 amr_1blk_restrict
!!     acts on UNK and/or FACEVAR[X][Y][Z]
!!     and/or UNK_E_[X][Y][Z] and/or UNK_N.
!!     If iopt=2 only WORK is updated.
!!
!!   integer, intent(in) iempty 
!!     NO FUNCTION !!!! Dummy argument kept for compatibility with
!!     earlier versions of PARAMESH.
!!
!!   logical, intent(in), optional :: filling_guardcells
!!     An optional argument which can be set to .true. if restriction
!!     is required as part of the guardcell filling step.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!
!! CALLS
!!
!!   mpi_amr_1blk_restrict
!!
!! RETURNS
!!
!!   Upon exit, restriction (i.e. interpolation) from leaf child blocks 
!!   (i.e. those blocks with nodetype = 1) to their parent has been
!!   performed.
!!
!! DESCRIPTION
!!
!!   This routine does the data averaging required when a child block
!!   passes data back to its parent. The parent receives interior data
!!   only, not guard cell data. 
!!
!! AUTHORS
!!
!!   Peter MacNeice (July 1997)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


subroutine amr_restrict(mype,iopt,iempty,filling_guardcells)
      use paramesh_dimensions
      use physicaldata
      use tree
      use timings

      use paramesh_mpi_interfaces, only :  & 
     &                  mpi_amr_1blk_restrict

      implicit none

      include 'mpif.h'
      double precision :: time1

      integer, intent(in)    :: mype,iopt,iempty
      logical, optional, intent(in)  :: filling_guardcells

      logical :: lcc,lfc,lec,lnc,lfulltree,fillingt

!------------------------------------------------------------------------
#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'entered mpi_amr_restrict'
#endif /* DEBUG_FLOW_TRACE */

      if(timing_mpi) then
        time1 = mpi_wtime()
      endif

      lrestrict_in_progress = .true.

      if (present(filling_guardcells)) then
         fillingt = filling_guardcells
      else
         fillingt = .false.
      end if

      if( (iopt.gt.1) .and. (mod(nxb,2).ne.0)  & 
     &     .and. (mype.eq.0) )  write(*,*)  & 
     &      'Restriction Warning ! Applying restriction to ', & 
     &      'WORK may lead to errors at external boundaries! '


      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(nvar.gt.0.or.iopt.ge.2) lcc = .true.
      if(nfacevar.gt.0.and.iopt.eq.1) lfc = .true.
      if(nvaredge.gt.0.and.iopt.eq.1) lec = .true.
      if(nvarcorn.gt.0.and.iopt.eq.1) lnc = .true.

!
! Perform restriction
      lfulltree = .false.
      call mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc, & 
     &                           lfulltree,fillingt)

      lrestrict_in_progress = .false.

      if(timing_mpi) then
        timer_amr_restrict =  timer_amr_restrict  & 
     &               + mpi_wtime() - time1
      endif

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'exiting mpi_amr_restrict'
#endif /* DEBUG_FLOW_TRACE */

      return

end subroutine amr_restrict

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine pf_restrict(mype,iopt,iempty,filling_guardcells, mglevel)

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings

      use paramesh_mpi_interfaces, only :  & 
     &                  mpi_amr_1blk_restrict

      implicit none

      include 'mpif.h'
      double precision :: time1

      integer, intent(in)    :: mype,iopt,iempty, mglevel
      logical, optional, intent(in)  :: filling_guardcells

      logical :: lcc,lfc,lec,lnc,lfulltree,fillingt

!------------------------------------------------------------------------
#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'entered mpi_amr_restrict'
#endif /* DEBUG_FLOW_TRACE */

      if(timing_mpi) then
        time1 = mpi_wtime()
      endif

      lrestrict_in_progress = .true.

      if (present(filling_guardcells)) then
         fillingt = filling_guardcells
      else
         fillingt = .false.
      end if

      if( (iopt.gt.1) .and. (mod(nxb,2).ne.0)  & 
     &     .and. (mype.eq.0) )  write(*,*)  & 
     &      'Restriction Warning ! Applying restriction to ', & 
     &      'WORK may lead to errors at external boundaries! '


      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(nvar.gt.0.or.iopt.ge.2) lcc = .true.
!
! Perform restriction
      lfulltree = .false.

      call pf_mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc, & 
     &                           lfulltree,fillingt, mglevel)

      lrestrict_in_progress = .false.

      if(timing_mpi) then
        timer_amr_restrict =  timer_amr_restrict  & 
     &               + mpi_wtime() - time1
      endif

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'exiting mpi_amr_restrict'
#endif /* DEBUG_FLOW_TRACE */

      return
      end subroutine pf_restrict
