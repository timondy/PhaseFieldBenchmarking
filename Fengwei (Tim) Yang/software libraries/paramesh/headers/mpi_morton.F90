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
!-----------------------------------------------------------------------
! mpi_morton Module



      Module mpi_morton

      Use paramesh_dimensions

      Private

      Public :: npts_neigh
      Integer, Parameter :: npts_neigh = 3000

! 
! variables for storing the morton environment
      Public :: pe_source
      Public :: mortonbnd,r_mortonbnd
      Public :: morton_limits
! CEG removed untyped mortonenv
!      Public :: mortonenv,ir_buf,is_buf
      Public :: ir_buf,is_buf
      Public :: no_of_mortonbnds_received
      Public :: morton_limits_set
      Public :: mpi_tree_set
      Integer, allocatable, Save ::  & 
     &    mortonbnd(:,:,:)
      Integer, Save :: no_of_mortonbnds_received
      Logical, Save :: morton_limits_set
      Logical, Save :: mpi_tree_set


      Integer, Save,dimension(:),allocatable :: pe_remote
      Integer, Save,dimension(:),allocatable :: pe_source
      Integer, Save,dimension(:),allocatable :: pe_destination
      Integer,  & 
     &  Save,dimension(:,:,:,:),allocatable :: r_mortonbnd
      Integer,  & 
     &  Save,dimension(:,:,:,:),allocatable :: morton_limits
      Integer, Save,dimension(:,:),allocatable :: ir_buf
      Integer, Save,dimension(:,:),allocatable :: is_buf

      Public :: commatrix_send, commatrix_recv
      Public :: commatrix_guard,commatrix_prol
      Public :: commatrix_flux
      Public :: commatrix_restrict
      Public :: laddress_guard,laddress_prol,laddress_flux
      Public :: laddress_restrict
      Integer, Save,dimension(:),allocatable :: commatrix_send
      Integer, Save,dimension(:),allocatable :: commatrix_recv
      Integer, Save,dimension(:,:),allocatable :: commatrix_guard
      Integer, Save,dimension(:,:),allocatable :: commatrix_prol
      Integer, Save,dimension(:,:),allocatable :: commatrix_flux
      Integer, Save,dimension(:,:),allocatable ::  & 
     &                                         commatrix_restrict
      Integer, allocatable, Save,dimension(:,:)  & 
     &                                          :: laddress_guard
      Integer, allocatable, Save,dimension(:,:)  & 
     &                                          :: laddress_prol
      Integer, allocatable, Save,dimension(:,:)  & 
     &                                          :: laddress_flux
      Integer, allocatable, Save,dimension(:,:)  & 
     &                                          :: laddress_restrict

! list of block edges which need diagonal info during edge averaging
      Public :: edge_mark,no_of_diagonal_edges
      Integer, Save :: edge_mark(6,4,npts_neigh)
      Integer, Save :: no_of_diagonal_edges

! a list of blocks to be sent from the local processor
      Public :: to_be_sent,to_be_sent_guard,to_be_sent_prol
      Public :: to_be_sent_flux
      Integer,Save,dimension(:,:,:),allocatable :: to_be_sent
      Integer,Save,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_sent_guard
      Integer,Save,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_sent_prol
      Integer,Save,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_sent_flux

      Public :: to_be_sent_restrict
      Integer,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_sent_restrict

! a list of blocks to be received by the local processor
      Public :: to_be_received
      Integer,Save,dimension(:,:,:),allocatable :: to_be_received

! Used to make searching of laddress more efficient
      Public :: ladd_strt,ladd_end
      Integer,Save,dimension(:),allocatable :: ladd_strt,ladd_end

!new code
      Public :: to_be_received_guard
      Public :: to_be_received_prol
      Public :: to_be_received_flux
      Integer,Save,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_received_guard
      Integer,Save,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_received_prol
      Integer,Save,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_received_flux

      Public :: to_be_received_restrict
      Integer,dimension(:,:,:),allocatable  & 
     &                                   :: to_be_received_restrict

      Public :: pe_source_guard
      Public :: pe_source_prol
      Public :: pe_source_flux
      Public :: pe_source_restrict
      Integer, Save,dimension(:),allocatable :: pe_source_guard
      Integer, Save,dimension(:),allocatable :: pe_source_prol
      Integer, Save,dimension(:),allocatable :: pe_source_flux
      Integer, Save,dimension(:),allocatable  & 
     &                                   :: pe_source_restrict

      Public :: message_size_cc
      Public :: message_size_fc
      Public :: message_size_ec
      Public :: message_size_nc
      Public :: message_size_wk
      Public :: mess_segment_loc
      Integer,Save,dimension(2*27) :: message_size_cc,message_size_fc
      Integer,Save,dimension(2*27) :: message_size_ec,message_size_nc
      Integer,Save,dimension(2*27) :: message_size_wk
      Integer,Save,dimension(:),allocatable :: mess_segment_loc

      Public :: temprecv_buf
      Real, Save,dimension(:),allocatable :: temprecv_buf

      Public :: l_datapacked
      Logical,Save,dimension(5) :: l_datapacked

!new code end

      Public :: largest_no_of_blocks,largest_no_of_blocks_guard
      Public :: largest_no_of_blocks_prol,largest_no_of_blocks_flux
      Public :: largest_no_of_blocks_restrict
      Public :: max_no_to_send,max_no_to_send_guard
      Public :: max_no_to_send_prol,max_no_to_send_flux
      Public :: max_no_to_send_restrict
      Public :: strt_guard,strt_prol,strt_flux
      Public :: strt_restrict,no_commatrix_guard
      Public :: no_commatrix_prol,no_commatrix_flux
      Public :: no_commatrix_restrict
      Integer,Save  :: largest_no_of_blocks
      Integer,Save  :: largest_no_of_blocks_guard
      Integer,Save  :: largest_no_of_blocks_prol
      Integer,Save  :: largest_no_of_blocks_flux
      Integer,Save  :: largest_no_of_blocks_restrict
      Integer,Save  :: max_no_to_send
      Integer,Save  :: max_no_to_send_guard
      Integer,Save  :: max_no_to_send_prol
      Integer,Save  :: max_no_to_send_flux
      Integer,Save  :: max_no_to_send_restrict
      Integer,Save  :: strt_guard
      Integer,Save  :: strt_prol
      Integer,Save  :: strt_flux
      Integer,Save  :: strt_restrict
      Integer,Save  :: no_commatrix_guard
      Integer,Save  :: no_commatrix_prol
      Integer,Save  :: no_commatrix_flux
      Integer,Save  :: no_commatrix_restrict

      Public :: lperiodicx,lperiodicy,lperiodicz
      Logical, Save :: lperiodicx
      Logical, Save :: lperiodicy
      Logical, Save :: lperiodicz

      Public :: treeinfo

      Type treeinfo
        Real coord(3)
        Real bsize(3)
        Real bnd_box(2,3)
        Integer parent(2)
        Integer which_child
        Logical newchild
        Integer neigh(2,6)
        Integer lrefine
        Integer nodetype
        Integer empty
      End Type treeinfo


      End Module mpi_morton
!-----------------------------------------------------------------------
