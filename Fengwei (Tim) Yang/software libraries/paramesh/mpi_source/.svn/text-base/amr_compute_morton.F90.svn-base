!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/amr_compute_morton
!! NAME
!!
!!   amr_compute_morton
!!
!! SYNOPSIS
!!
!!   call amr_compute_morton(mort_no)
!!
!!   call amr_compute_morton(integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(out) ::  mort_no(:,:)
!!     This is a list of the morton numbers computed and returned by this subroutine.  Note
!!     that this is a 2 dimensional array.  For high levels of refinement, more bits are required
!!     to store the morton numbers, we spread these bit patterns over several 32 (or 64) bit
!!     integers.  The first dimension represents the number of these integers which are used and
!!     the 2nd dimension is the length of the block list.
!!
!! INCLUDES
!!
!!   mpif.h
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   io
!!   paramesh_mpi_interfaces
!!   paramesh_comm_data
!!
!! CALLS
!! 
!!   morton_number
!!    
!! RETURNS
!!
!!   Returns the computed morton number for each block.
!!
!! DESCRIPTION
!!
!!   The routine computes the morton numbers for the blocks on a processor.  It
!!   does this by interleaving bits of the integer values of the coordinates of
!!   each individual block.  The subroutine 'morton_number' actuall does the 
!!   computation of the morton number for a particular block and this routine
!!   simply computes these integer coordinates and then loops over the blocks
!!   and calls the subroutine 'morton_number'.
!!
!!   Normally this routine is only called by 'amr_morton_order' and the user will 
!!   need to call it directly.  It is intended to be a support routine for process
!!   reordering blocks during refinement/derefinement.
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2004).
!!   Mike Zingale and Jonathan Dursi (1999).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_compute_morton (mort_no)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use paramesh_mpi_interfaces, only : morton_number
      Use io
      Use paramesh_comm_data

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(out) ::  mort_no(:,:)

!-----Local variables and arrays.
      Real :: xmin,ymin,zmin,xmin_loc,ymin_loc,zmin_loc
      Real :: xyz_loc_vector(3), xyz_min_loc_vector(3)
      Real :: x0,y0,z0
      Integer :: i
      Integer :: ierr
      Integer :: mype

!-----Begin executable code.

      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)

!-----find local minimum + maximum values of x, y, and z
      xmin_loc = 1.e10
      ymin_loc = 1.e10
      zmin_loc = 1.e10
      Do i = 1,lnblocks
         If (nodetype(i) == 1) then
            xmin_loc = min(coord(1,i)-(bsize(1,i)/2.),xmin_loc)
            If (ndim >= 2) then
               ymin_loc = min(coord(2,i)-(bsize(2,i)/2.),ymin_loc)
            End If
            If (ndim == 3) then
               zmin_loc = min(coord(3,i)-(bsize(3,i)/2.),zmin_loc)
            End If
         End If
      End Do

!-----find global min^s across processors
!-----(Changed by M. Zingale and J. Dursi)
!-----pack the 3 allreduces into a single allreduce with a vector of the
!-----minimum in each coordinate
      xyz_loc_vector(1) = xmin_loc
      xyz_loc_vector(2) = ymin_loc
      xyz_loc_vector(3) = zmin_loc

!-----reduce in all ndim dimensions
      Call MPI_ALLREDUCE(xyz_loc_vector, xyz_min_loc_vector, ndim,     & 
           amr_mpi_real,                                               &  
           MPI_MIN, MPI_COMM_WORLD, ierr)

!-----unpack the minimums
      xmin = 0.
      ymin = 0.       ! initialize to avoid NaN in 1D with some compilers
      zmin = 0.       ! initialize to avoid NaN in 1D or 2D with some compilers
      xmin = xyz_min_loc_vector(1)
      If (ndim >= 2) ymin = xyz_min_loc_vector(2)
      If (ndim >= 3) zmin = xyz_min_loc_vector(3)

      Do i = 1,lnblocks
         x0 = coord(1,i)-xmin
         y0 = coord(2,i)-ymin
         z0 = coord(3,i)-zmin
         Call morton_number(x0,y0,z0,bsize(:,i),ndim,                  & 
                            lrefine_max,lrefine(i),mort_no(:,i))
      End Do

      Return
      End Subroutine amr_compute_morton

