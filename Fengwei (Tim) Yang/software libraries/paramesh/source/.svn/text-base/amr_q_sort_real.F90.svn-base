!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_q_sort_real
!! NAME
!!
!!   amr_q_sort_real
!!
!! SYNOPSIS
!!
!!   call amr_q_sort_real (ix, n, ia, ib)
!!
!!   call amr_q_sort_real(real, integer, integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   real, dimension(n), intent(inout) :: ix
!!        The keys to be sorted.
!!
!!   integer, intent(in) :: n           
!!        Length of key array, ix.
!!
!!   integer, dimension(n), intent(inout) :: ia, ib
!!        Additional arrays to be permuted into the sorted order determined by
!!        ix.  
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   Uses nothing else.
!!
!! CALLS
!!
!!   Calls no other paramesh subroutines.  Can be used 'stand-alone'.
!!
!! RETURNS
!!
!!   Returns the sorted keys (ix) of length n.  Upon return the arrays ia, ib are 
!!   permuted into the order determined by the sorted keys (ix). 
!!
!! DESCRIPTION
!!
!!   This routine performs a sort using a quik sort algorithm. It takes ix, sorts 
!!   it by calling q_sort_1, which also returns the permutation array iperm, and
!!   then applies iperm to ia and ib also. This routine assumes ix is a 'real' variable.
!!
!! AUTHORS
!!
!!   Kevin Olson (April 2003)
!!   Bug fixes for Intel 8.0 compiler contributed by Breno Imbiriba.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_q_sort_real (ix,n,ia,ib)

      implicit none

      integer, intent(in) :: n
      real, dimension(n),  intent(inout) :: ix
      integer, optional, dimension(n), intent(inout) :: ia, ib
      integer, allocatable, dimension(:) :: ia_temp
      integer, allocatable, dimension(:) :: ib_temp
      integer, dimension(n) :: iperm
      integer :: i

      do i = 1,n
         iperm(i) = i
      end do

      call q_sort_1r(1,n)

      if (present(ia)) then
      allocate(ia_temp(n))
      ia_temp(:) = ia(:)
      do i = 1,n
         ia(i) = ia_temp(iperm(i))
      end do
      deallocate(ia_temp)
      end if

      if (present(ib)) then
      allocate(ib_temp(n))
      ib_temp(:) = ib(:)
      do i = 1,n
         ib(i) = ib_temp(iperm(i))
      end do
      deallocate(ib_temp)
      end if

      contains

      recursive subroutine q_sort_1r (ismall,ibig)

      integer, intent(in) :: ismall, ibig
      integer :: i, j
      real :: pivot, temp
      integer :: itemp
      integer, parameter :: max_qsort_size = 10

      if (ibig < ismall + max_qsort_size) then

         call simple_sortr(ismall, ibig)

      else

         pivot = ix((ismall + ibig)/2)
         i = ismall-1
         j = ibig+1

         do
            do
               i = i + 1
               if (ix(i) >= pivot) exit
            end do

            do
               j = j - 1
               if (ix(j) <= pivot) exit
            end do

            if (i < j) then
               temp = ix(i)
               ix(i) = ix(j)
               ix(j) = temp
               itemp = iperm(i)
               iperm(i) = iperm(j)
               iperm(j) = itemp
            else if (i == j) then
               i = i + 1
               exit
            else
               exit
            end if
         end do

         if (ismall < j) call q_sort_1r(ismall, j)
         if (i < ibig) call q_sort_1r(i, ibig)

      end if

      end subroutine q_sort_1r

      subroutine simple_sortr(ismall, ibig)

      integer, intent(in) :: ismall, ibig
      integer :: i, j
      real :: temp
      integer :: itemp

      do i = ismall, ibig-1
         do j = i + 1, ibig
            if (ix(i) > ix(j)) then
               temp = ix(i)
               ix(i) = ix(j)
               ix(j) = temp
               itemp = iperm(i)
               iperm(i) = iperm(j)
               iperm(j) = itemp
            end if
         end do
      end do

      end subroutine simple_sortr

      end subroutine amr_q_sort_real
