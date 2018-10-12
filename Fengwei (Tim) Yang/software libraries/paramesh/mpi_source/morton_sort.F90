!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/morton_sort
!! NAME
!!
!!   morton_sort
!!
!! SYNOPSIS
!!
!!   call morton_sort(mort_no,ix,iend)
!!   call morton_sort(integer array, integer array, integer array)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: iend
!!     size of morton number list.
!!   integer, intent(inout) :: mort_no(6,iend), ix(iend)
!!     mort_no -> List of morton numbers to sort.
!!     ix      -> sorted indeces to be used in permutations after sort is 
!!                called.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_interfaces
!!
!! CALLS
!! 
!!   amr_q_sort
!!    
!! RETURNS
!!
!!   Returns a sorted morton list and a set of keys which can be used 
!!   to permute other lists into morton order.
!!
!! DESCRIPTION
!!
!!   This routine sorts morton numbers stored locally on a single processor.
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2001).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine morton_sort(mort_no,ix,iend)

!-----Use statements.
      Use paramesh_interfaces, Only : amr_q_sort
      Use io

      Implicit None

      Integer, Intent(in) :: iend
      Integer, Intent(inout) :: mort_no(6,iend), ix(iend)

!-----Local arrays and variables.
      Integer :: ix2(iend)
      Integer :: i, is, ie, level
      Integer :: mort_no_old(6,iend)
      Logical :: greater

!-----Begin executeable code.
      If (iend > 0) Then

      Do i = 1,iend
         ix2(i) = i
      End Do

      is = 1
      ie = iend
      mort_no_old(:,is:ie) = mort_no(:,is:ie)

      greater = .False.
      Do i = is, ie-1
         If (mort_no(1,i) > mort_no(1,i+1)) Then
            greater = .True.
            Exit
         End If
      End Do
      If (greater) Then
         Call amr_q_sort(mort_no(1,is:ie),                             & 
                         ie-is+1,                                      & 
                         ix2(is:ie),                                   & 
                         ix(is:ie))
         Do i = is,ie
            mort_no(:,i) = mort_no_old(:,ix2(i))
         End Do
      End If

      Do level = 2,6
         ie = 0
         Do While (ie < iend)
         is = ie + 1
         ie = is
         Do
            If (any(mort_no(1:level-1,is) .ne.                         & 
                    mort_no(1:level-1,ie))) Then
               ie = ie - 1
               Exit
            End If
            ie = ie + 1
            If (ie > iend) Then
               ie = ie - 1
               Exit
            End If
         End Do

         If (ie > is) Then
         mort_no_old(:,is:ie) = mort_no(:,is:ie)

         greater = .False.
         Do i = is, ie-1
            If (mort_no(level,i) > mort_no(level,i+1)) Then
               greater = .True.
               exit
            End If
         End Do
         If (greater) Then
            Do i = is,ie
               ix2(i) = i
            End Do
            Call amr_q_sort(mort_no(level,is:ie),                      & 
                            ie-is+1,                                   & 
                            ix2(is:ie),                                & 
                            ix(is:ie))
            Do i = is,ie
               mort_no(:,i) = mort_no_old(:,ix2(i))
            End Do
         End If
         
         End If

         End Do  ! End Do While (ie < iend)
      End Do  ! End Do level = 2,6

      End If  ! End If (iend > 0)

      Return
      End Subroutine morton_sort
