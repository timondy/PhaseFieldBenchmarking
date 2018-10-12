!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/set_f2c_indexes
!! NAME
!!
!!   set_f2c_indexes
!!
!! SYNOPSIS
!!
!!   Call set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)
!!   Call set_f2c_indexes(logical,logical,logical,logical,integer,integer)
!!
!! ARGUMENTS
!!
!!   Logical, Intent(in) :: lcc, lfc, lec, lnc logical switches to indicate which
!!                                             set of variables to work on
!!   Integer, Intent(in) :: nblks_ind no. of input blocks
!!   Integer, Intent(in) :: iopt indicates whether work array is being used
!!
!! INLCUDES
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   workspace
!!
!! CALLS
!!
!! RETURNS
!!
!! DESCRIPTION
!! 
!!   A routine to compute offsets and extents into arrays for fine to course 
!!   interpolcation (i.e. restriction).
!!
!! AUTHORS
!!
!!   Peter MacNeice
!!
!!***

      Subroutine set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)

!-----Use Statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace

      Implicit None

!-----Input/Output variables
      Logical, Intent(in) :: lcc,lfc,lec,lnc
      Integer, Intent(in) :: nblks_ind,iopt

!---------Begin executable code.

          If(lcc) Then
          If(iopt.eq.1) Then
            f2c_ind_unk(1,1,nblks_ind) = iu_bnd1                       & 
                              + nxb + nguard
            f2c_ind_unk(1,2,nblks_ind) = ju_bnd1                       & 
                              + (nyb + nguard)*k2d
            f2c_ind_unk(1,3,nblks_ind) = ku_bnd1                       & 
                              + (nzb + nguard)*k3d
            f2c_ind_unk(2,1,nblks_ind) = il_bnd1
            f2c_ind_unk(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unk(2,3,nblks_ind) = kl_bnd1
          else
            f2c_ind_work(1,1,nblks_ind) = iuw1                         & 
                              + nxb + nguard_work
            f2c_ind_work(1,2,nblks_ind) = juw1                         & 
                              + (nyb + nguard_work)*k2d
            f2c_ind_work(1,3,nblks_ind) = kuw1                         & 
                              + (nzb + nguard_work)*k3d
            f2c_ind_work(2,1,nblks_ind) = ilw1
            f2c_ind_work(2,2,nblks_ind) = jlw1
            f2c_ind_work(2,3,nblks_ind) = klw1
          End If
          End If  ! End If (lcc)
          If(lfc) Then
            f2c_ind_facex(1,1,nblks_ind) = iu_bnd1 + 1                 & 
                              + nxb + nguard
            f2c_ind_facex(1,2,nblks_ind) = ju_bnd1                     & 
                              + (nyb + nguard)*k2d
            f2c_ind_facex(1,3,nblks_ind) = ku_bnd1                     & 
                              + (nzb + nguard)*k3d
            f2c_ind_facex(2,1,nblks_ind) = il_bnd1
            f2c_ind_facex(2,2,nblks_ind) = jl_bnd1
            f2c_ind_facex(2,3,nblks_ind) = kl_bnd1
            f2c_ind_facey(1,1,nblks_ind) = iu_bnd1                     & 
                              + nxb + nguard
            f2c_ind_facey(1,2,nblks_ind) = ju_bnd1 + k2d               & 
                              + (nyb + nguard)*k2d
            f2c_ind_facey(1,3,nblks_ind) = ku_bnd1                     & 
                              + (nzb + nguard)*k3d
            f2c_ind_facey(2,1,nblks_ind) = il_bnd1
            f2c_ind_facey(2,2,nblks_ind) = jl_bnd1
            f2c_ind_facey(2,3,nblks_ind) = kl_bnd1
            f2c_ind_facez(1,1,nblks_ind) = iu_bnd1                     & 
                              + nxb + nguard
            f2c_ind_facez(1,2,nblks_ind) = ju_bnd1                     & 
                              + (nyb + nguard)*k2d
            f2c_ind_facez(1,3,nblks_ind) = ku_bnd1 + k3d               &  
                              + (nzb + nguard)*k3d
            f2c_ind_facez(2,1,nblks_ind) = il_bnd1
            f2c_ind_facez(2,2,nblks_ind) = jl_bnd1
            f2c_ind_facez(2,3,nblks_ind) = kl_bnd1
          End If  ! End If (lfc)
          If(lec) Then
            f2c_ind_unkex(1,1,nblks_ind) = iu_bnd1                     & 
                              + nxb + nguard
            f2c_ind_unkex(1,2,nblks_ind) = ju_bnd1 + k2d               & 
                              + (nyb + nguard)*k2d
            f2c_ind_unkex(1,3,nblks_ind) = ku_bnd1 + k3d               & 
                              + (nzb + nguard)*k3d
            f2c_ind_unkex(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkex(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkex(2,3,nblks_ind) = kl_bnd1
            f2c_ind_unkey(1,1,nblks_ind) = iu_bnd1 + 1                 & 
                              + nxb + nguard
            f2c_ind_unkey(1,2,nblks_ind) = ju_bnd1                     & 
                              + (nyb + nguard)*k2d
            f2c_ind_unkey(1,3,nblks_ind) = ku_bnd1 + k3d               & 
                              + (nzb + nguard)*k3d
            f2c_ind_unkey(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkey(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkey(2,3,nblks_ind) = kl_bnd1
            f2c_ind_unkez(1,1,nblks_ind) = iu_bnd1 + 1                 & 
                              + nxb + nguard
            f2c_ind_unkez(1,2,nblks_ind) = ju_bnd1 + k2d               & 
                              + (nyb + nguard)*k2d
            f2c_ind_unkez(1,3,nblks_ind) = ku_bnd1                     & 
                              + (nzb + nguard)*k3d
            f2c_ind_unkez(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkez(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkez(2,3,nblks_ind) = kl_bnd1
          End If  ! End If (lec)
          If(lnc) Then
            f2c_ind_unkn(1,1,nblks_ind) = iu_bnd1 + 1                  & 
                              + nxb + nguard
            f2c_ind_unkn(1,2,nblks_ind) = ju_bnd1 + k2d                & 
                              + (nyb + nguard)*k2d
            f2c_ind_unkn(1,3,nblks_ind) = ku_bnd1 + k3d                & 
                              + (nzb + nguard)*k3d
            f2c_ind_unkn(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkn(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkn(2,3,nblks_ind) = kl_bnd1
          End If  ! End If (lnc)

          Return

      End Subroutine set_f2c_indexes
