!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      module paramesh_mpi_interfaces


      interface
      subroutine mpi_setup(mype,nprocs)
      integer, intent(in) :: mype,nprocs
      end subroutine mpi_setup
      end interface

      interface
      subroutine pf_mpi_setup(mype,nprocs)
      integer, intent(in) :: mype,nprocs
      end subroutine pf_mpi_setup
      end interface

      interface
      subroutine boundary_locator(x0,y0,z0,lboundary,ibc)
      logical,intent(out) :: lboundary
      integer,intent(out) :: ibc
      real,intent(in)     :: x0,y0,z0
      end subroutine boundary_locator
      end interface

      interface
      subroutine mpi_morton_bnd(mype,nprocs,tag_offset)
      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: tag_offset
      end subroutine mpi_morton_bnd
      end interface

      interface
      subroutine pf_morton_bnd(mype,nprocs,tag_offset)
      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: tag_offset
      end subroutine pf_morton_bnd
      end interface

      interface
      subroutine mpi_amr_get_remote_block(mype,remote_pe,remote_block, & 
     &    idest,iopt,lcc,lfc,lec,lnc, & 
     &    nlayersx,nlayersy,nlayersz)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,iopt
      logical, intent(in) :: lcc,lfc,lec,lnc
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_amr_get_remote_block
      end interface

      interface
      subroutine mpi_amr_get_remote_block_fvar(mype, & 
     &                     remote_pe,remote_block,icoord, & 
     &                     recvx,recvy,recvz,idest)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: icoord,idest
      real, intent(out)   :: recvx(:,:,:,:)
      real, intent(out)   :: recvy(:,:,:,:)
      real, intent(out)   :: recvz(:,:,:,:)
      end subroutine mpi_amr_get_remote_block_fvar
      end interface

      interface
      subroutine mpi_amr_comm_setup(mype,nprocs,lguard,lprolong, & 
     &                              lflux,ledge,lrestrict,lfulltree, & 
     &                              iopt,lcc,lfc,lec,lnc,tag_offset, & 
     &                              nlayersx,nlayersy,nlayersz, & 
     &                              flux_dir)
      integer, intent(in)    :: mype,nprocs,iopt
      integer, intent(inout) :: tag_offset
      logical, intent(in)    :: lcc,lfc,lec,lnc,lfulltree
      logical, intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      integer, intent(in), optional :: flux_dir
      end subroutine mpi_amr_comm_setup
      end interface

      interface
      subroutine pf_mpi_amr_comm_setup(mype,nprocs,lguard,lprolong, & 
     &                              lflux,ledge,lrestrict,lfulltree, & 
     &                              iopt,lcc,lfc,lec,lnc,tag_offset, & 
     &                              mglevel)
      integer, intent(in)    :: mype,nprocs,iopt
      integer, intent(inout) :: tag_offset
      logical, intent(in)    :: lcc,lfc,lec,lnc,lfulltree
      logical, intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict, mglevel
      end subroutine pf_mpi_amr_comm_setup
      end interface

      interface
      subroutine mpi_amr_comm_setup_res(mype,nprocs,lguard,lprolong, & 
     &                              lflux,ledge,lrestrict, & 
     &                              iopt,lcc,lfc,lec,lnc, & 
     &                              tag_offset,level)
      integer, intent(in)    :: mype,nprocs,iopt
      integer, intent(inout) :: tag_offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      logical, intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict
      integer, intent(in)    :: level
      end subroutine mpi_amr_comm_setup_res
      end interface

      interface
      subroutine mpi_amr_tree_setup( & 
     &                         mype,nprocs,tag_offset)
      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: tag_offset
      end subroutine mpi_amr_tree_setup
      end interface

      interface
      subroutine mpi_morton_bnd_prolong(mype,nprocs, & 
     &                                  tag_offset)
      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: tag_offset
      end subroutine mpi_morton_bnd_prolong
      end interface
      interface

      subroutine mpi_morton_bnd_fluxcon(mype,nprocs, & 
     &                                  tag_offset)
      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: tag_offset
      end subroutine mpi_morton_bnd_fluxcon
      end interface

      interface
      subroutine mpi_morton_bnd_restrict(mype,nprocs, & 
     &                                  tag_offset)
      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: tag_offset
      end subroutine mpi_morton_bnd_restrict
      end interface

      interface
      subroutine mpi_amr_local_surr_blks(mype,lb,nprocs, & 
     &                          max_no_of_blocks, & 
     &                          surrblks,l_parent,psurrblks)
      integer, intent(in)    ::  mype,lb,nprocs,max_no_of_blocks
      integer, intent(out)   ::  surrblks(:,:,:,:)
      integer, intent(out)   ::  psurrblks(:,:,:,:)
      logical, intent(in)    ::  l_parent
      end subroutine mpi_amr_local_surr_blks
      end interface

      interface
      subroutine mpi_amr_local_surr_blks_lkup(mype,lb, & 
     &                          surrblks,l_parent,psurrblks)
      integer, intent(in)    ::  mype,lb
      integer, intent(out)   ::  surrblks(:,:,:,:)
      integer, intent(out)   ::  psurrblks(:,:,:,:)
      logical, intent(in)    ::  l_parent
      end subroutine mpi_amr_local_surr_blks_lkup
      end interface

      interface
      subroutine mpi_amr_1blk_guardcell_c_to_f( & 
     &                    mype,lb,pe,iopt,nlayers, & 
     &                    surrblks,lcc,lfc,lec,lnc,icoord,ldiag, & 
     &                    nlayersx,nlayersy,nlayersz,ipolar)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx,nlayersy,nlayersz
      integer, intent(in) :: ipolar(2)
      end subroutine mpi_amr_1blk_guardcell_c_to_f
      end interface

      interface
      subroutine mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc, & 
     &                                 lfulltree,filling_guardcells)
      integer, intent(in)  :: mype,iopt
      logical, intent(in)  :: lcc,lfc,lec,lnc,lfulltree
      logical, intent(in)  :: filling_guardcells
      end subroutine mpi_amr_1blk_restrict
      end interface

      interface
      subroutine pf_mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc, & 
     &                                 lfulltree,filling_guardcells, mglevel)
      integer, intent(in)  :: mype,iopt
      logical, intent(in)  :: lcc,lfc,lec,lnc,lfulltree
      logical, intent(in)  :: filling_guardcells, mglevel
      end subroutine pf_mpi_amr_1blk_restrict
      end interface

      interface
      subroutine mpi_array_allocate(nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_array_allocate
      end interface

      interface
      subroutine mpi_array_deallocate
      end subroutine mpi_array_deallocate
      end interface

      interface
      subroutine mpi_amr_global_domain_limits
      end subroutine mpi_amr_global_domain_limits
      end interface

      interface
      subroutine mpi_set_message_limits(dtype,ia,ib,ja,jb,ka,kb,vtype, & 
     &                                  nlayersx,nlayersy,nlayersz)
      integer, intent(in)  ::  dtype,vtype
      integer, intent(out) ::  ia,ib,ja,jb,ka,kb
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_set_message_limits
      end interface

      interface
      subroutine mpi_pack_blocks(mype,nprocs,iopt, & 
     &                           lcc,lfc,lec,lnc, & 
     &                           buf_dim,S_buffer,offset, & 
     &                           nlayersx,nlayersy,nlayersz)
      integer, intent(in)  ::  mype,nprocs,iopt
      logical, intent(in)  ::  lcc,lfc,lec,lnc
      integer, intent(in)  ::  buf_dim,offset
      real,    intent(out) ::  S_buffer(buf_dim)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_pack_blocks
      end interface

      interface
      subroutine mpi_Sbuffer_size(mype,nprocs,iopt, & 
     &                           lcc,lfc,lec,lnc, & 
     &                           buf_dim,offset, & 
     &                           block_sections, fluxes, edges, &
     &                           flux_dir, &
     &                           nlayersx,nlayersy,nlayersz)
      integer, intent(in)  ::  mype,nprocs,iopt
      logical, intent(in)  ::  lcc,lfc,lec,lnc
      integer, intent(out) ::  buf_dim
      integer, intent(in)  ::  offset
      logical, intent(in)  ::  block_sections, fluxes, edges
      integer, intent(in), optional :: flux_dir
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_Sbuffer_size
      end interface

      interface
      subroutine mpi_pack_edges(mype,nprocs, & 
     &                          buf_dim,S_buffer,offset)
      integer, intent(in)  ::  mype,nprocs
      integer, intent(in)  ::  buf_dim,offset
      real,    intent(out) ::  S_buffer(buf_dim)
      end subroutine mpi_pack_edges
      end interface


      interface
      subroutine mpi_pack_fluxes(mype,nprocs, & 
     &                          buf_dim,S_buffer,offset,flux_dir)
      integer, intent(in)  ::  buf_dim,offset
      real,    intent(out) ::  S_buffer(buf_dim)
      integer, intent(in)  ::  mype,nprocs
      integer, intent(in), optional :: flux_dir
      end subroutine mpi_pack_fluxes
      end interface

      interface
      subroutine mpi_pack_tree_info(mype,nprocs, & 
     &                          buf_dim_bytes,buf_dim,S_buffer)
      integer, intent(in)  ::  mype,nprocs
      integer, intent(in)  ::  buf_dim,buf_dim_bytes
      real,    intent(out) ::  S_buffer(buf_dim)
      end subroutine mpi_pack_tree_info
      end interface


      interface
      subroutine mpi_unpack_blocks(mype,iopt, & 
     &                             lcc,lfc,lec,lnc, & 
     &                             buf_dim,R_buffer, & 
     &                             nlayersx,nlayersy,nlayersz)
      integer, intent(in) :: mype,buf_dim,iopt
      logical, intent(in) :: lcc,lfc,lec,lnc
      real,    intent(inout) ::  R_buffer(buf_dim)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_unpack_blocks
      end interface

      interface
      subroutine mpi_Rbuffer_size(mype,nprocs,iopt, & 
     &                            lcc,lfc,lec,lnc, & 
     &                            buf_dim, & 
     &                            block_sections, fluxes, edges, flux_dir, &
     &                            nlayersx,nlayersy,nlayersz)
      integer, intent(in)  :: mype, nprocs, iopt
      integer, intent(out) :: buf_dim
      logical, intent(in)  :: lcc,lfc,lec,lnc
      logical, intent(in)  :: block_sections, fluxes, edges
      integer, intent(in), optional :: flux_dir
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_Rbuffer_size
      end interface

      interface
      subroutine mpi_unpack_edges(mype, & 
     &                          buf_dim,R_buffer)
      integer, intent(in)  ::  mype
      integer, intent(in)  ::  buf_dim
      real,    intent(inout) ::  R_buffer(buf_dim)
      end subroutine mpi_unpack_edges
      end interface

      interface
      subroutine mpi_unpack_fluxes(mype, & 
     &                          buf_dim,R_buffer,flux_dir)
      integer, intent(in)  ::  mype
      integer, intent(in)  ::  buf_dim
      real,    intent(inout) ::  R_buffer(buf_dim)
      integer, optional, intent(in) :: flux_dir
      end subroutine mpi_unpack_fluxes
      end interface

      interface
      subroutine mpi_unpack_tree_info(mype,nprocs, & 
     &                          buf_dim_bytes,buf_dim,R_buffer)
      integer, intent(in)  ::  mype,nprocs
      integer, intent(in)  ::  buf_dim,buf_dim_bytes
      real,    intent(inout) ::  R_buffer(buf_dim)
      end subroutine mpi_unpack_tree_info
      end interface



      interface
      subroutine mpi_put_buffer(lb,ioptw,offset, & 
     &                          lcc,lfc,lec,lnc, & 
     &                          buffer_size,R_buffer, & 
     &                          nlayersx,nlayersy,nlayersz)
      integer, intent(in)    :: lb,ioptw,buffer_size
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      real,    intent(inout) :: R_buffer(buffer_size)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_put_buffer
      end interface

      interface
      subroutine mpi_get_Rbuffer_size(lb,dtype,ioptw,offset, & 
     &                                lcc,lfc,lec,lnc, & 
     &                                nlayersx,nlayersy,nlayersz)
      integer, intent(in)    :: lb,dtype,ioptw
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_get_Rbuffer_size
      end interface

      interface
      subroutine mpi_put_edge_buffer(lb,offset, & 
     &                          buffer_size,R_buffer)
      integer, intent(in)    :: lb,buffer_size
      integer, intent(inout) :: offset
      real,    intent(inout) :: R_buffer(buffer_size)
      end subroutine mpi_put_edge_buffer
      end interface

      interface
      subroutine mpi_get_Rbuffer_size_edges(lb,dtype,offset)
      integer, intent(in)    :: lb,dtype
      integer, intent(inout) :: offset
      end subroutine mpi_get_Rbuffer_Size_edges
      end interface

      interface
      subroutine mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)
      integer, intent(in)    :: lb,remote_block,remote_pe
      end subroutine mpi_put_edge_buffer_1blk
      end interface

      interface
      subroutine mpi_put_flux_buffer(mype,lb,offset, & 
     &                          buffer_size,R_buffer,flux_dir)
      integer, intent(in)    :: mype,lb,buffer_size
      integer, intent(inout) :: offset
      real,    intent(inout) :: R_buffer(buffer_size)
      integer, optional, intent(in) :: flux_dir
      end subroutine mpi_put_flux_buffer
      end interface
      interface

      subroutine mpi_get_Rbuffer_size_fluxes(lb,dtype,offset, & 
     &                                       flux_dir)
      integer, intent(in)    :: lb,dtype
      integer, intent(inout) :: offset
      integer, optional, intent(in) :: flux_dir
      end subroutine mpi_get_Rbuffer_size_fluxes
      end interface

      interface
      subroutine mpi_get_buffer(mype,lb,dtype,iopt,offset, & 
     &                          lcc,lfc,lec,lnc, & 
     &                          buffer_size,S_buffer, & 
     &                          nlayersx,nlayersy,nlayersz)
      integer, intent(in)    :: dtype
      integer, intent(in)    :: lb,mype,iopt,buffer_size
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      real,    intent(inout) :: S_buffer(buffer_size)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_get_buffer
      end interface

      interface
      subroutine mpi_get_Sbuffer_size(mype,lb,dtype,iopt,offset, & 
     &                                lcc,lfc,lec,lnc, & 
     &                                nlayersx,nlayersy,nlayersz)
      integer, intent(in)    :: dtype
      integer, intent(in)    :: lb,mype,iopt
      integer, intent(inout) :: offset
      logical, intent(in)    :: lcc,lfc,lec,lnc
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_get_Sbuffer_size
      end interface

      interface
      subroutine mpi_get_edge_buffer(mype,lb,dtype,offset, & 
     &                          buffer_size,S_buffer)
      integer, intent(in)    :: dtype
      integer, intent(in)    :: lb,mype,buffer_size
      integer, intent(inout) :: offset
      real,    intent(inout) :: S_buffer(buffer_size)
      end subroutine mpi_get_edge_buffer
      end interface

      interface
      subroutine mpi_get_Sbuffer_size_edges(mype,lb,dtype,offset)
      integer, intent(in)    :: dtype
      integer, intent(in)    :: lb,mype
      integer, intent(inout) :: offset
      end subroutine mpi_get_Sbuffer_size_edges
      end interface

      interface
      subroutine mpi_get_flux_buffer(mype,lb,dtype,offset, & 
     &                          buffer_size,S_buffer,flux_dir)
      integer, intent(in)    :: dtype
      integer, intent(in)    :: lb,mype,buffer_size
      integer, intent(inout) :: offset
      real,    intent(inout) :: S_buffer(buffer_size)
      integer, optional, intent(in) :: flux_dir
      end subroutine mpi_get_flux_buffer
      end interface

      interface
      subroutine mpi_get_Sbuffer_size_fluxes(mype,lb,dtype,offset, & 
     &                                       flux_dir)
      integer, intent(in)    :: dtype
      integer, intent(in)    :: lb,mype
      integer, intent(inout) :: offset
      integer, optional, intent(in) :: flux_dir
      end subroutine mpi_get_Sbuffer_size_fluxes
      end interface

      interface
      subroutine mpi_set_message_sizes(iopt, & 
     &                                 nlayersx,nlayersy,nlayersz)
      integer, intent(in)    :: iopt
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine mpi_set_message_sizes
      end interface

      interface
      subroutine pf_mpi_set_message_sizes(iopt)
      integer, intent(in)    :: iopt
      end subroutine pf_mpi_set_message_sizes
      end interface

      interface
      subroutine mpi_amr_write_guard_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_write_guard_comm
      end interface

      interface
      subroutine mpi_amr_read_guard_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_read_guard_comm
      end interface

      interface
      subroutine mpi_amr_write_prol_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_write_prol_comm
      end interface

      interface
      subroutine mpi_amr_read_prol_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_read_prol_comm
      end interface

      interface
      subroutine mpi_amr_write_flux_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_write_flux_comm
      end interface

      interface
      subroutine mpi_amr_read_flux_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_read_flux_comm
      end interface

      interface
      subroutine mpi_amr_write_restrict_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_write_restrict_comm
      end interface

      interface
      subroutine mpi_amr_read_restrict_comm (nprocs)
      integer, intent(in) :: nprocs
      end subroutine mpi_amr_read_restrict_comm
      end interface

      interface
      subroutine mpi_xchange_blocks(mype,nprocs, tag_offset, & 
     &                              buf_dim_send, S_buffer, &
     &                              buf_dim_recv, R_buffer)
      integer, intent(in)    ::  mype,nprocs,buf_dim_send,buf_dim_recv
      integer, intent(inout) ::  tag_offset
      real,    intent(inout) ::  S_buffer(buf_dim_send)
      real,    intent(inout) ::  R_buffer(buf_dim_recv)
      end subroutine mpi_xchange_blocks
      end interface

      interface
      subroutine mpi_xchange_tree_info(mype,nprocs, tag_offset, & 
     &                              buf_dim, S_buffer, R_buffer)
      integer, intent(in)    ::  mype,nprocs,buf_dim
      integer, intent(inout) ::  tag_offset
      real,    intent(inout) :: S_buffer(buf_dim), R_buffer(buf_dim)
      end subroutine mpi_xchange_tree_info
      end interface

      interface
      subroutine morton_number( & 
     &     x0,y0,z0,bbsize,ndim,lrefine_max,lrefine,mort)

      integer,intent(in ) :: lrefine,lrefine_max,ndim
      real,intent(in)     :: bbsize(3),x0,y0,z0
      integer,intent(out) :: mort(6)
      end subroutine morton_number
      end interface

      interface 
      subroutine mpi_amr_boundary_block_info(mype,nprocs)
      integer,intent(in) :: mype,nprocs
      end subroutine mpi_amr_boundary_block_info
      end interface

      interface
      subroutine mpi_amr_get_bc_settings(blk, & 
     &                           ibxl,ibxr,ibyl,ibyr,ibzl,ibzr)
      integer,intent(in)    :: blk
      integer,intent(inout) :: ibxl,ibxr,ibyl,ibyr,ibzl,ibzr
      end subroutine mpi_amr_get_bc_settings
      end interface      

      interface
      subroutine compress_list(neigh_morts, & 
     &                         istack,no_of_remote_neighs, & 
     &                         mype,nprocs,l_on_pe)

      integer, intent(inout) :: istack
      integer, intent(inout), dimension(:,:,:) :: neigh_morts
      integer, intent(out) :: no_of_remote_neighs
      integer, intent(in)  :: mype, nprocs
      logical, intent(in)  :: l_on_pe

      end subroutine compress_list
      end interface

      interface
      subroutine compress_fetch_list(fetch_list, & 
     &                               istack,no_of_remote_neighs, & 
     &                                mype,nprocs,n_to_left)

      integer, intent(inout), dimension(:,:) :: fetch_list
      integer, intent(out) :: no_of_remote_neighs
      integer, intent(in)  :: istack, mype, nprocs
      integer, intent(in)  :: n_to_left(0:nprocs-1)

      end subroutine compress_fetch_list
      end interface

      interface

      subroutine process_fetch_list(fetch_list,                        & 
                                    istack,mype,nprocs,n_to_left,      &
                                    tag_offset)

      Integer, Intent(inout), dimension(:,:) :: fetch_list
      Integer, Intent(in)  :: istack, mype, nprocs
      Integer, Intent(in)  :: n_to_left(0:nprocs-1)
      Integer, Intent(inout) :: tag_offset

      end subroutine process_fetch_list
      end interface

      interface
      subroutine mpi_amr_singular_line (istrategy,nprocs)

      integer, intent(in) :: istrategy,nprocs

      end subroutine mpi_amr_singular_line
      end interface

      interface

      subroutine morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                            lperiodicx,lperiodicy,lperiodicz, & 
     &                            coord,bbsize,ndim, & 
     &                            lrefine,lrefine_max,mort_neigh, & 
     &                            bndbox)

      real, intent(in)     :: xmin,ymin,zmin,xmax,ymax,zmax
      real, intent(in)     :: coord(3),bbsize(3)
      real,optional, intent(in)     :: bndbox(2,3)
      integer, intent(in)  :: ndim
      integer, intent(in)  :: lrefine,lrefine_max
      logical, intent(in)  :: lperiodicx,lperiodicy,lperiodicz
      integer, intent(out) :: mort_neigh(6,3,3,3)
      end subroutine morton_neighbors
      end interface

      interface

      subroutine morton_child_neighbors(xmin,ymin,zmin,xmax,ymax,zmax, & 
     &                            lperiodicx,lperiodicy,lperiodicz, & 
     &                            coord,bbsize,ndim, & 
     &                            lrefine,lrefine_max,mort_neigh, & 
     &                            bndbox)

      real, intent(in)     :: xmin,ymin,zmin,xmax,ymax,zmax
      real, intent(in)     :: coord(3),bbsize(3)
      real,optional, intent(in)     :: bndbox(2,3)
      integer, intent(in)  :: ndim
      integer, intent(in)  :: lrefine,lrefine_max
      logical, intent(in)  :: lperiodicx,lperiodicy,lperiodicz
      integer, intent(out) :: mort_neigh(6,4,4,4)
      end subroutine morton_child_neighbors
      end interface

      interface
      subroutine rationalize_fetch_list (fetch_list,                   &
                                         istart,                       &
                                         iend,                         &
                                         nptsneigh)
      integer, intent(in) :: istart, iend, nptsneigh
      integer, intent(inout) :: fetch_list(3,nptsneigh)

      end subroutine rationalize_fetch_list
      end interface

      interface
      subroutine amr_morton_process()
      end subroutine amr_morton_process
      end interface


      end module paramesh_mpi_interfaces
