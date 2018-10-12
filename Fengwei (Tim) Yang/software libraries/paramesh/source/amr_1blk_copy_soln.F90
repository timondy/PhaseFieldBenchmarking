!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------


!!****f* source/amr_1blk_copy_soln
!! NAME
!!   
!!   amr_1blk_copy_soln
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_copy_soln(level)
!!   Call amr_1blk_copy_soln(integer)
!!
!! ARGUMENTS
!!
!!    integer, intent(IN) :: level     
!!      If 'level' is -1 then blocks at all refinement
!!      levels are copied, otherwise only blocks
!!      at 'level' are copied.
!!    
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   mpi_amr_global_domain_limits
!!   mpi_morton_bnd
!!
!! RETURNS
!!
!!   Does not return anything.  Upon return data from the gt_* arrays are 
!!   copied to the global arrays (unk, facevarx(y,z), unk_e_x(y,z) and unk_n).
!!
!! DESCRIPTION
!!
!!   This routine copies a global solution update from the time 
!!   synchronized global solution arrays (the gt_* arrays) into (unk, 
!!   facevarx(y,z), unk_e_x(y,z), unk_n) into the arrays used during the 
!!   solution update (gt_* arrays), as is required when 'no_permanent_guardcells 
!!   is set to .TRUE. and the amr_1blk_guardcell routines are called.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          May 1999
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      Subroutine amr_1blk_copy_soln(level)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use paramesh_mpi_interfaces, only :                   & 
                             mpi_amr_global_domain_limits,  & 
                             mpi_morton_bnd

      Implicit None

      include 'mpif.h'

!-----Input/Output Variables
      Integer, Intent(in) :: level

!-----Local Variables
      Integer          :: mype,nprocs,lb,ierr
      Integer          :: tag_offset, ivar
      Integer          :: nguard0
      Double Precision :: time1

!-----Begin Exectuable Code

      If (timing_mpi) Then
         time1 = mpi_wtime()
      End If

      nguard0 = nguard*npgs

      If (no_permanent_guardcells) Then

        Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)


        If (level == -1) Then

        Do lb = 1, lnblocks

!-------Cell Centered Variables
        If (nvar > 0) then
          Do ivar=1,nvar
            If (int_gcell_on_cc(ivar)) Then
              gt_unk(ivar,:,:,:,lb) = unk(ivar,:,:,:,lb)
            End If  ! End If (int_gcell_on_cc(ivar))
          End Do  ! End Do ivar=1,nvar
        End If  ! End If (nvar > 0) 

!-------Face Centered Variables
        If (nfacevar > 0) Then
          Do ivar=1,nfacevar
            If (int_gcell_on_fc(1,ivar)) Then
              gt_facevarx(ivar,:,:,:,lb) =  & 
                 facevarx(ivar,:,:,:,lb)
            End If  ! End If (int_gcell_on_fc)
          End Do  ! End Do ivar=1,nfacevar

          If (ndim >= 2) Then
            Do ivar=1,nfacevar
              If (int_gcell_on_fc(2,ivar)) Then
                gt_facevary(ivar,:,:,:,lb) =  & 
                   facevary(ivar,:,:,:,lb)
              End If  ! End If (int_gcell_on_fc(2,ivar))
            End Do  ! End Do ivar=1,nfacevar
          End If  ! End if (ndim >= 2)

          If (ndim == 3) Then
            Do ivar=1,nfacevar
              If (int_gcell_on_fc(3,ivar)) Then
                gt_facevarz(ivar,:,:,:,lb) =  & 
                   facevarz(ivar,:,:,:,lb)
              End If  ! End If (int_gcell_on_fc(3,ivar))
            End Do  ! End Do ivar=1,nfacevar
          End If  ! End If (ndim == 3)
        End If  ! End If (nfacevar > 0)

!-------Edge Centered Variables
        If (nvaredge > 0) Then
          If (ndim > 1) Then
            Do ivar=1,nvaredge
              If (int_gcell_on_ec(1,ivar)) Then
                gt_unk_e_x(ivar,:,:,:,lb) =  & 
                  unk_e_x(ivar,:,:,:,lb)
              End If  ! End If (int_gcell_on_ec(1,ivar))
            End Do  ! End Do ivar=1,nvaredge
            Do ivar=1,nvaredge
              If (int_gcell_on_ec(2,ivar)) Then
                gt_unk_e_y(ivar,:,:,:,lb) =  & 
                  unk_e_y(ivar,:,:,:,lb)
              End If  ! End If (int_gcell_on_ec(2,ivar))
            End Do  ! End Do ivar=1,nvaredge
            If (ndim == 3) Then
              Do ivar=1,nvaredge
                If (int_gcell_on_ec(3,ivar)) Then
                  gt_unk_e_z(ivar,:,:,:,lb) =  & 
                    unk_e_z(ivar,:,:,:,lb)
                End If  ! End If (int_gcell_on_ec(3,ivar))
              End Do  ! End Do ivar=1,nvaredge
            End If  ! End If (ndim == 2)
          End If  ! End If (ndim > 1)
        End If  ! End If (nvaredge > 0)

!-------Corner Centered Variables
        If (nvarcorn > 0) Then
          Do ivar=1,nvarcorn
          If (int_gcell_on_nc(ivar)) Then
            gt_unk_n(ivar,:,:,:,lb) = unk_n(ivar,:,:,:,lb)
          End If  ! End If (int_gcell_on_nc(ivar))
          End Do  ! End Do ivar=1,nvarcorn
        End If  ! End If (nvarcorn > 0) 

        End Do !  End Do lb = 1, lnblocks

      Else  ! If (level == -1)

        Do lb = 1, lnblocks
        If (lrefine(lb) == level) then

!---------Cell Centered Variables
          If (nvar > 0) Then
          Do ivar=1,nvar
          If (int_gcell_on_cc(ivar)) then
            gt_unk(ivar,:,:,:,lb) = unk(ivar,:,:,:,lb)
          End If  ! End If (int_gcell_on_cc(ivar))
          End Do  ! End Do ivar=1,nvar
          End If  ! End If (nvar > 0)

!---------Face Centered Variables
          If(nfacevar > 0) Then
          Do ivar=1,nfacevar
          If (int_gcell_on_fc(1,ivar)) Then
            gt_facevarx(ivar,:,:,:,lb) = facevarx(ivar,:,:,:,lb)
          End If  ! If (int_gcell_on_fc(1,ivar))
          End Do  ! End Do ivar=1,nfacevar
          If (ndim >= 2) Then
          Do ivar=1,nfacevar
          If (int_gcell_on_fc(2,ivar)) then
            gt_facevary(ivar,:,:,:,lb) = facevary(ivar,:,:,:,lb)
          End If  ! End If (int_gcell_on_fc(2,ivar))
          End Do  ! End Do ivar=1,nfacevar
          End If  ! End If (nvar >= 2)
          If (ndim == 3) Then
          Do ivar=1,nfacevar
          If (int_gcell_on_fc(3,ivar)) Then
            gt_facevarz(ivar,:,:,:,lb) = facevarz(ivar,:,:,:,lb)
          End If  ! If (int_gcell_on_fc(3,ivar))
          End Do  ! End Do ivar=1,nfacevar
          End If  ! End If (ndim == 3)
          End If  ! End If (nfacevar > 0)

!---------Edge Centered Variables
          If(nvaredge > 0) Then
          If (ndim > 1) Then
          Do ivar=1,nvaredge
          If (int_gcell_on_ec(1,ivar)) Then
            gt_unk_e_x(ivar,:,:,:,lb) = unk_e_x(ivar,:,:,:,lb)
          End If  ! End If (int_gcell_on_ec(i,ivar)
          End Do  ! End Do ivar=1,nvaredge
          Do ivar=1,nvaredge
          If (int_gcell_on_ec(2,ivar)) Then
            gt_unk_e_y(ivar,:,:,:,lb) = unk_e_y(ivar,:,:,:,lb)
          End If  ! End If (int_gcell_on_ec(2,ivar))
          End Do  ! End Do ivar=1,nvaredge
          If (ndim == 3) Then
          Do ivar=1,nvaredge
          If (int_gcell_on_ec(3,ivar)) Then
            gt_unk_e_z(ivar,:,:,:,lb) = unk_e_z(ivar,:,:,:,lb)
          End If  ! End If (int_gcell_on_ec(3,ivar)
          End Do  ! End Do ivar=1,nvaredge
          End If  ! End If (ndim == 3)
          End If  ! End If (ndim > 1)
          End If  ! End If (nvaredge >0)

!---------Corner Center Variables
          If (nvarcorn > 0) Then
          Do ivar=1,nvarcorn
          If (int_gcell_on_nc(ivar)) Then
            gt_unk_n(ivar,:,:,:,lb) = unk_n(ivar,:,:,:,lb)
          End If  ! End If (int_gcell_on_nc(ivar))
          End Do  ! End Do ivar=1,nvarcorn
          End If  ! End If (nvarcorn > 0)

        End If  ! End If (lrefine(lb) == level)
        End Do  ! End Do lb = 1, lnblocks

      End If  ! End If (level == -1)

  print *,'CEG'
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!-----Since this routine is a prelude to calling amr_1blk_guardcell
!-----we check here to make sure that surrblks has been computed and
!-----stored for each block.

      If (gsurrblks_set.ne.1) Then

        Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)

!--------Find the computational domain coordinate ranges
        Call mpi_amr_global_domain_limits

!-------Set up surrounding blocks of all local blocks (must not precede
!-------setting of grid_xmin,... etc)
        tag_offset = 100
        If (nprocs.gt.1)                              & 
          Call mpi_morton_bnd(mype,nprocs,tag_offset)

      End If  ! End If (gsurrblks_set.ne.1)

      Else  !  If (no_permanent_guardcells)

!-----Enforce Consistency of face centered variables at block boundaries
      If (force_consistency) Then
      If (nfacevar > 0) Then
        Do lb = 1,lnblocks

          Do ivar=1,nfacevar
          If (int_gcell_on_fc(1,ivar)) Then
            gt_facevarx(ivar,1,:,:,lb) = facevarx(ivar,1+nguard0,:,:,lb)
            gt_facevarx(ivar,2,:,:,lb) =                               & 
                                   facevarx(ivar,nxb+1+nguard0,:,:,lb)
          End if  ! End If (int_gcell_on_fc(1,ivar))
          End Do  ! End Do ivar=1,nfacevar

          If (ndim >= 2) Then
          Do ivar=1,nfacevar
          If (int_gcell_on_fc(2,ivar)) Then
            gt_facevary(ivar,:,1,:,lb) =                               & 
                                     facevary(ivar,:,1+nguard0*k2d,:,lb)
            gt_facevary(ivar,:,1+k2d,:,lb) =                           & 
                               facevary(ivar,:,nyb+(1+nguard0)*k2d,:,lb)
          End If  ! End If (int_gcell_on_fc(2,ivar))
          End Do  ! End Do ivar=1,nfacevar
          End If  ! End If If (ndim >= 2)

          If (ndim == 3) Then
          Do ivar=1,nfacevar
          If (int_gcell_on_fc(3,ivar)) Then
            gt_facevarz(ivar,:,:,1,lb) =                              & 
                                   facevarz(ivar,:,:,1+nguard0*k3d,lb)
            gt_facevarz(ivar,:,:,1+k3d,lb) =                          & 
                             facevarz(ivar,:,:,nzb+(1+nguard0)*k3d,lb)
          End If  ! End If (int_gcell_on_fc(3,ivar))
          End Do  ! End Do ivar=1,nfacevar
          End If  ! If (ndim == 3)

        End Do  ! End Do lb = 1,lnblocks
      End If  ! End If (nfacevar > 0)
      End If  ! End If (force_consistency)


!      print *,'Really in amr_1blk_copy_soln'
      End If ! End If (no_permanent_guardcells)

!  print *,'CEG'
!      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      If (timing_mpi) Then
              timer_amr_1blk_copy_soln =  timer_amr_1blk_copy_soln  & 
                                + mpi_wtime() - time1
      End If

      Return
      End Subroutine amr_1blk_copy_soln
