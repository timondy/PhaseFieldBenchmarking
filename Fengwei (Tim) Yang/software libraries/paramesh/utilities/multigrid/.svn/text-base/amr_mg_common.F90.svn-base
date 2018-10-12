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


      module amr_mg_common

      use tree
      use paramesh_dimensions

      public :: nodetype_old
      integer, allocatable, save :: nodetype_old(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! multigrid support
      public :: to_be_data_type
      type to_be_data_type
        integer, pointer, dimension(:,:,:) :: & 
     &      to_be_sent_guard_mg, & 
     &      to_be_sent_prol_mg, & 
     &      to_be_sent_flux_mg, & 
     &      to_be_sent_restrict_mg, & 
     &      to_be_received_guard_mg, & 
     &      to_be_received_prol_mg, & 
     &      to_be_received_flux_mg, & 
     &      to_be_received_restrict_mg
      end type to_be_data_type

      public :: to_be_levels
      type(to_be_data_type), save, dimension(:), allocatable :: & 
     &      to_be_levels

      public :: pe_source_guard_mg
      public :: pe_source_prol_mg
      public :: pe_source_flux_mg
      public :: pe_source_restrict_mg
      integer, save,dimension(:,:),allocatable :: pe_source_guard_mg
      integer, save,dimension(:,:),allocatable :: pe_source_prol_mg
      integer, save,dimension(:,:),allocatable :: pe_source_flux_mg
      integer, save,dimension(:,:),allocatable & 
     &                                   :: pe_source_restrict_mg

      public :: commatrix_guard_mg,commatrix_prol_mg
      public :: commatrix_flux_mg
      public :: commatrix_restrict_mg
      public :: laddress_guard_mg,laddress_prol_mg,laddress_flux_mg
      public :: laddress_restrict_mg
      integer, save,dimension(:,:,:),allocatable :: commatrix_guard_mg
      integer, save,dimension(:,:,:),allocatable :: commatrix_prol_mg
      integer, save,dimension(:,:,:),allocatable :: commatrix_flux_mg
      integer, save,dimension(:,:,:),allocatable :: & 
     &                                         commatrix_restrict_mg

      integer, allocatable, save,dimension(:,:,:) & 
     &                                          :: laddress_guard_mg
      integer, allocatable, save,dimension(:,:,:) & 
     &                                          :: laddress_prol_mg
      integer, allocatable, save,dimension(:,:,:) & 
     &                                          :: laddress_flux_mg
      integer, allocatable, save,dimension(:,:,:) & 
     &                                          :: laddress_restrict_mg

      public :: largest_no_of_blocks_guard_mg
      public :: largest_no_of_blocks_prol_mg
      public :: largest_no_of_blocks_flux_mg
      public :: largest_restrict_mg
      public :: max_no_to_send_guard_mg
      public :: max_no_to_send_prol_mg,max_no_to_send_flux_mg
      public :: max_no_to_send_restrict_mg
      public :: strt_guard_mg,strt_prol_mg,strt_flux_mg
      public :: strt_restrict_mg
      integer,dimension(:),allocatable,save :: & 
     & largest_no_of_blocks_guard_mg
      integer,dimension(:),allocatable,save :: & 
     & largest_no_of_blocks_prol_mg
      integer,dimension(:),allocatable,save :: & 
     & largest_no_of_blocks_flux_mg
      integer,dimension(:),allocatable,save :: & 
     & largest_restrict_mg
      integer,dimension(:),allocatable,save :: max_no_to_send_guard_mg
      integer,dimension(:),allocatable,save :: max_no_to_send_prol_mg
      integer,dimension(:),allocatable,save :: max_no_to_send_flux_mg
      integer,dimension(:),allocatable,save :: & 
     & max_no_to_send_restrict_mg
      integer,dimension(:),allocatable,save :: strt_guard_mg
      integer,dimension(:),allocatable,save :: strt_prol_mg
      integer,dimension(:),allocatable,save :: strt_flux_mg
      integer,dimension(:),allocatable,save :: strt_restrict_mg

    
      end module amr_mg_common
      
