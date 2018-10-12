!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! Modification history:
!
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!
!#ifdef HAVE_CONFIG_H
!#include <config.h>
!#endif

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      module paramesh_interfaces



      interface
      subroutine amr_1blk_bc(mype,cneigh,iopt,nlayers,lb,pe,             & 
     &                       idest,icoord)
      integer, intent(in) :: mype,iopt,nlayers
      integer, intent(in) :: lb,pe,idest,icoord
      integer, dimension(2,6), intent(in) :: cneigh
      end subroutine amr_1blk_bc
      end interface

      interface
      subroutine amr_1blk_bcset(mype,ibc,lb,pe,                          & 
     &    idest,iopt,ibnd,jbnd,kbnd,surrblks)
      integer, intent(in) :: mype,ibc,lb,pe
      integer, intent(in) :: idest,iopt,ibnd,jbnd,kbnd
      integer, intent(in) :: surrblks(:,:,:,:)
      end subroutine amr_1blk_bcset
      end interface

      interface
      subroutine amr_1blk_cc_cp_remote(mype,remote_pe,remote_block,      & 
     &    idest,iopt,id,jd,kd,is,js,ks,ilays,jlays,klays,nblk_ind)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,iopt,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays,nblk_ind
      end subroutine amr_1blk_cc_cp_remote
      end interface

      interface
      subroutine amr_1blk_nc_cp_remote(mype,remote_pe,remote_block,      & 
     &    idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,         & 
     &    ip3,jp3,kp3,nblk_ind)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip3,jp3,kp3
      integer, intent(in) :: nblk_ind
      end subroutine amr_1blk_nc_cp_remote
      end interface

      interface
      subroutine amr_1blk_cc_prol_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb,    & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_cc_prol_gen_unk_fun
      end interface

      interface
      subroutine amr_1blk_cc_prol_inject(recv,ia,ib,ja,jb,ka,kb,         & 
     &       idest,ioff,joff,koff,mype,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_cc_prol_inject
      end interface

      interface
      subroutine amr_1blk_cc_prol_linear(recv,ia,ib,ja,jb,ka,kb,          & 
     &       idest,ioff,joff,koff,mype,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_cc_prol_linear
      end interface

      interface
      subroutine amr_1blk_cc_prol_genorder(recv,ia,ib,ja,jb,ka,kb,        & 
     &       idest,ioff,joff,koff,mype,ivar,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,order
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_cc_prol_genorder
      end interface

      interface
      subroutine amr_1blk_cc_prol_user(recv,ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_cc_prol_user
      end interface

      interface
      subroutine amr_1blk_cc_prol_gen_work_fun(recv,                     & 
     &       ia,ib,ja,jb,ka,kb,                                          & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p,interp
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_gen_work_fun
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_inject(recv,                     & 
     &       ia,ib,ja,jb,ka,kb,                                        & 
     &       idest,ioff,joff,koff,mype)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_work_inject
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_linear(recv,                   & 
     &       ia,ib,ja,jb,ka,kb,                                      & 
     &       idest,ioff,joff,koff,mype)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_work_linear
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_genorder(recv,             & 
     &       ia,ib,ja,jb,ka,kb,                                  & 
     &       idest,ioff,joff,koff,mype,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: order
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_work_genorder
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_user(recv, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_work_user
      end interface

      interface
      subroutine amr_1blk_copy_soln(level)
      integer, intent(in) :: level
      end subroutine amr_1blk_copy_soln
      end interface

      interface
      subroutine amr_1blk_ec_cp_remote(mype,remote_pe,remote_block, & 
     &   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,    & 
     &    ip2,jp2,kp2,ip3,jp3,kp3,iface,nblk_ind)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3,iface
      integer, intent(in) :: nblk_ind
      end subroutine amr_1blk_ec_cp_remote 
      end interface

      interface
      subroutine amr_1blk_fc_cp_remote(mype,remote_pe,remote_block, & 
     &   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,    & 
     &    ip2,jp2,kp2,iface,nblk_ind)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,iface
      integer, intent(in) :: nblk_ind
      end subroutine amr_1blk_fc_cp_remote 
      end interface

      interface
      subroutine amr_1blk_ec_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype,iface)
      integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in)    :: ioff,joff,koff,mype,iface
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_ec_prol_gen_fun
      end interface

      interface
      subroutine amr_1blk_ec_prol_linear & 
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &      mype,ivar,iedge_dir)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_ec_prol_linear
      end interface

      interface
      subroutine amr_1blk_ec_prol_genorder & 
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &      mype,ivar,iedge_dir,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,iedge_dir,order
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_ec_prol_genorder
      end interface

      interface
      subroutine amr_1blk_ec_prol_user()
      end subroutine amr_1blk_ec_prol_user
      end interface
 
      interface
      subroutine amr_1blk_nc_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_nc_prol_gen_fun
      end interface

      interface
      subroutine amr_1blk_nc_prol_linear(recv,ia,ib,ja,jb,ka,kb,idest,   & 
     &       ioff,joff,koff,mype,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_nc_prol_linear
      end interface

      interface
      subroutine amr_1blk_nc_prol_genorder(                              & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,ivar,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,order
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_nc_prol_genorder
      end interface

      interface
      subroutine amr_1blk_nc_prol_user()
      end subroutine amr_1blk_nc_prol_user
      end interface

      interface
      subroutine amr_1blk_fc_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest,  & 
     &       ioff,joff,koff,mype,lb,pe_p,lb_p,iface)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_gen_fun
      end interface

      interface
      subroutine prol_fc_dbz_init(n,i_fc_vars)
      integer, intent(in) :: n, i_fc_vars(:,:)
      end subroutine prol_fc_dbz_init
      end interface

      interface
      function prol_fc_dbz_varp(ivar, iface) result(ldbz)
      integer, intent(in) :: ivar, iface
      logical :: ldbz
      end function prol_fc_dbz_varp
      end interface
      
      interface
      subroutine amr_1blk_fc_prol_dbz(                                   & 
     &        recvfx, recvfy, recvfz,                                    & 
     &        nfacevar, iv1, iv2, iv3                                    & 
     &        ,ia,ib,ja,jb,ka,kb,                                        & 
     &        idest,ioff,joff,koff,                                      & 
     &        mype,lb,parent_pe,parent_blk)
      real, intent(inout), dimension(:,:,:,:) :: recvfx,recvfy,recvfz
      integer, intent(in)    :: nfacevar, iv1, iv2, iv3
      integer, intent(in)    :: ia,ib,ja,jb,ka,kb
      integer, intent(in)    :: idest,ioff,joff,koff
      integer, intent(in)    :: mype,lb,parent_pe,parent_blk
      end subroutine amr_1blk_fc_prol_dbz
      end interface

      interface
      subroutine amr_1blk_fc_prol_inject(                                & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,iface,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_inject
      end interface


      interface
      subroutine amr_1blk_fc_prol_linear(                                & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,iface,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_linear
      end interface

      interface
      subroutine amr_1blk_fc_prol_genorder(                              & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,iface,ivar,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,iface,order
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_genorder
      end interface

      interface
      subroutine amr_1blk_fc_prol_user( & 
     &       recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype,lb,pe_p,lb_p,iface,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)
      integer, intent(in) :: ivar
      end subroutine amr_1blk_fc_prol_user
      end interface

      interface
      subroutine amr_1blk_guardcell(mype,iopt,nlayers,lb,pe,             & 
     &                              lcc,lfc,lec,lnc,                     & 
     &                              l_srl_only,icoord,ldiag,             & 
     &                              nlayersx, nlayersy, nlayersz)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      logical, intent(in) :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      end subroutine amr_1blk_guardcell
      end interface


      interface
      subroutine amr_1blk_guardcell_c_to_f(mype,lb,pe,iopt,nlayers,      & 
     &                         surrblks,lcc,lfc,lec,lnc,icoord,ldiag,    & 
     &                         nlayersx,nlayersy,nlayersz)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx, nlayersy, nlayersz
      end subroutine amr_1blk_guardcell_c_to_f
      end interface


      interface
      subroutine amr_1blk_guardcell_f_to_c(mype,pe,lb,iblock,iopt,       & 
     &                                                    nlayers,       & 
     &                         surrblks,lcc,lfc,lec,lnc,icoord,ldiag,    & 
     &                         nlayersx,nlayersy,nlayersz)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx,nlayersy,nlayersz
      end subroutine amr_1blk_guardcell_f_to_c
      end interface

      interface
      subroutine amr_1blk_guardcell_f_to_c_fil(mype,pe,lb,iblock,iopt,   & 
     &                                                    nlayers,       & 
     &                         surrblks,lcc,lfc,lec,lnc,icoord,ldiag)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      end subroutine amr_1blk_guardcell_f_to_c_fil
      end interface

      interface
      subroutine amr_1blk_guardcell_f_to_c_set(mype,pe,lb,iblock,iopt,   & 
     &                                                    nlayers,       & 
     &                         surrblks,lcc,lfc,lec,lnc,icoord,ldiag,    & 
     &                         nlayersx,nlayersy,nlayersz)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx,nlayersy,nlayersz
      end subroutine amr_1blk_guardcell_f_to_c_set
      end interface

      interface
      subroutine amr_1blk_guardcell_reset
      end subroutine amr_1blk_guardcell_reset
      end interface


      interface
      subroutine amr_1blk_guardcell_srl(mype,pe,lb,iblock,iopt,nlayers,  & 
     &                         surrblks,lcc,lfc,lec,lnc,icoord,ldiag,    & 
     &                         nlayers0x,nlayers0y,nlayers0z)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayers0x, nlayers0y, nlayers0z
      end subroutine amr_1blk_guardcell_srl
      end interface


      interface
      subroutine amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc)
      integer, intent(in)  :: mype,iopt
      logical, intent(in)  :: lcc,lfc,lec,lnc
      end subroutine amr_1blk_restrict
      end interface


      interface
      subroutine amr_1blk_save_soln
      end subroutine amr_1blk_save_soln
      end interface


      interface
      subroutine amr_1blk_t_to_perm( lcc,lfc,lec,lnc,lb,idest)
      integer, intent(in) :: lb,idest
      logical, intent(in) :: lcc,lfc,lec,lnc
      end subroutine amr_1blk_t_to_perm
      end interface

      interface
      subroutine amr_1blk_to_perm(lcc,lfc,lec,lnc,lb,iopt,idest)
      integer, intent(in) :: lb,iopt,idest
      logical, intent(in) :: lcc,lfc,lec,lnc
      end subroutine amr_1blk_to_perm
      end interface


      interface
      subroutine amr_abort()
      end subroutine amr_abort
      end interface

      interface
      subroutine amr_bi_sort(list,gid,npp)
      integer,intent(inout) :: list(:)
      integer,intent(inout) :: gid(:)
      integer,intent(in) :: npp
      end subroutine amr_bi_sort
      end interface

      interface
      subroutine amr_bc_block(jface,ibc,iopt,l,mype)
      integer, intent(in) :: jface,ibc,iopt,l,mype
      end subroutine amr_bc_block
      end interface

      interface
      subroutine amr_bcset_init
      end subroutine amr_bcset_init
      end interface


      interface
      subroutine amr_block_geometry(lb,pe)
      integer, intent(in) :: lb,pe
      end subroutine amr_block_geometry
      end interface

      interface
      subroutine amr_checkpoint_re(iunit1,l_with_guardcells)
      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      end subroutine amr_checkpoint_re
      end interface


      interface
      subroutine amr_checkpoint_wr(iunit1,l_with_guardcells)
      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      end subroutine amr_checkpoint_wr
      end interface

      interface
      subroutine amr_check_derefine(mype)
      integer, intent(in) :: mype
      end subroutine amr_check_derefine
      end interface

      interface
      subroutine amr_check_refine(nprocs,mype,icontinue)
      integer, intent(in) :: nprocs,mype
      integer, intent(in) :: icontinue
      end subroutine amr_check_refine
      end interface


      interface
      subroutine amr_close
      end subroutine amr_close
      end interface

      interface
      subroutine amr_derefine_blocks(lnblocks_old,mype)
      integer, intent(in)    :: mype
      integer, intent(inout) :: lnblocks_old
      end subroutine amr_derefine_blocks
      end interface

      interface
      subroutine amr_compute_morton (mort_no)
      integer, intent(out) ::  mort_no(:,:)
      end subroutine amr_compute_morton 
      end interface


      interface
      subroutine amr_edge_average(mype,lfullblock,nsub)
      integer, intent(in)  ::  mype,nsub
      logical, intent(in)  ::  lfullblock
      end subroutine amr_edge_average
      end interface


      interface
      subroutine amr_edge_average_udt(mype)
      integer, intent(in)  ::  mype
      end subroutine amr_edge_average_udt
      end interface


      interface
      subroutine amr_edge_average_vdt(mype,nsub)
      integer, intent(in)  ::  mype,nsub
      end subroutine amr_edge_average_vdt
      end interface

      interface
      subroutine amr_edge_diagonal_check(mype)
      integer, intent(in)  ::  mype
      end subroutine amr_edge_diagonal_check
      end interface

      interface
      subroutine amr_flush(iunit)
      integer, intent(in)  ::  iunit
      end subroutine amr_flush
      end interface

      interface
      subroutine amr_flux_conserve(mype,nsub,flux_dir)
      integer, optional, intent(in)  ::  flux_dir
      integer, intent(in)  ::  mype,nsub
      end subroutine amr_flux_conserve
      end interface

      interface
      subroutine amr_flux_conserve_udt(mype,flux_dir)
      integer, optional, intent(in) :: flux_dir
      integer, intent(in)  ::  mype
      end subroutine amr_flux_conserve_udt
      end interface

      interface
      subroutine amr_flux_conserve_vdt(mype,nsub)
      integer, intent(in)  ::  mype,nsub
      end subroutine amr_flux_conserve_vdt
      end interface

      interface
      subroutine amr_gsurrounding_blks(mype,ldiag)
      integer, intent(in)    ::  mype
      logical, intent(in)    ::  ldiag
      end subroutine amr_gsurrounding_blks
      end interface


      interface
      subroutine amr_guardcell(mype,iopt,nlayers, & 
     &                         nlayersx,nlayersy,nlayersz)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      integer, intent(in)  ::  mype,iopt,nlayers
      end subroutine amr_guardcell
      end interface

      interface
      subroutine amr_initialize
      end subroutine amr_initialize
      end interface


      interface
      subroutine amr_migrate_tree_data(new_loc,nprocs,mype)
      integer, intent(in)    ::  mype
      integer, intent(in)    ::  nprocs
      integer, intent(inout) ::  new_loc(:,:)
      end subroutine amr_migrate_tree_data
      end interface


      interface
      subroutine amr_morton_order (lnblocks_old,nprocs,mype, & 
     &                             l_move_solution)
      integer, intent(in) ::  mype
      integer, intent(in) ::  nprocs,lnblocks_old
      logical, intent(in) ::  l_move_solution
      end subroutine amr_morton_order
      end interface



      interface
      subroutine amr_perm_to_1blk( lcc,lfc,lec,lnc,lb,pe,iopt,idest)
      integer, intent(in) ::  lb,pe,iopt,idest
      logical, intent(in) ::  lcc,lfc,lec,lnc
      end subroutine amr_perm_to_1blk
      end interface

      interface
      subroutine amr_mpi_find_blk_in_buffer( & 
     &       mype,remote_block,remote_pe,idest,dtype,index,lfound)
      integer, intent(in)  :: mype,remote_pe,remote_block,idest
      integer, intent(out) :: dtype,index
      logical, intent(out) :: lfound
      end subroutine amr_mpi_find_blk_in_buffer
      end interface


      interface
      subroutine amr_prolong(mype,iopt,nlayers)
      integer, intent(in) ::  mype,iopt,nlayers
      end subroutine amr_prolong
      end interface

      interface
      subroutine amr_prolong_cc_fun_init
      end subroutine amr_prolong_cc_fun_init
      end interface

      interface
      subroutine amr_prolong_face_fun_init
      end subroutine amr_prolong_face_fun_init
      end interface

      interface
!      subroutine amr_prolong_fc_divbconsist(mype)
      subroutine amr_prolong_fc_divbconsist(mype,level,nfield)
      integer, intent(in) ::  mype
      integer, intent(in) ::  level
      integer, intent(in) ::  nfield
      end subroutine amr_prolong_fc_divbconsist
      end interface


      interface
      subroutine amr_prolong_fun_init
      end subroutine amr_prolong_fun_init
      end interface

      interface
      subroutine amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)
      integer, intent(in)    ::  nprocs
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in)    ::  lnblocks_old
      integer, intent(in)    ::  mype
      end subroutine amr_redist_blk
      end interface


      interface
      subroutine amr_refine_blocks (nprocs,mype)
      integer, intent(in)    :: nprocs,mype
      end subroutine amr_refine_blocks
      end interface


      interface
      subroutine amr_refine_derefine
      end subroutine amr_refine_derefine
      end interface


      interface
      subroutine amr_restrict(mype,iopt,iempty,filling_guardcells)
      integer, intent(in)    :: mype,iopt,iempty
      logical, optional, intent(in) :: filling_guardcells
      end subroutine amr_restrict
      end interface

      interface
      subroutine amr_restrict_bnd_data(mype,flux_dir)
      integer, intent(in)    :: flux_dir
      integer, intent(in)    :: mype
      end subroutine amr_restrict_bnd_data
      end interface

      interface
      subroutine amr_restrict_bnd_data_vdt(mype)
      integer, intent(in)    :: mype
      end subroutine amr_restrict_bnd_data_vdt
      end interface

      interface
      subroutine amr_restrict_edge(icoord)
      integer, intent(in)    :: icoord
      end subroutine amr_restrict_edge
      end interface

      interface
      subroutine amr_restrict_edge_data(mype)
      integer, intent(in)    :: mype
      end subroutine amr_restrict_edge_data
      end interface

      interface
      subroutine amr_restrict_edge_data_vdt(mype)
      integer, intent(in)    :: mype
      end subroutine amr_restrict_edge_data_vdt
      end interface

      interface
      subroutine amr_restrict_ec_fun(recv,temp,icoord)
      integer, intent(in)    :: icoord
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_ec_fun
      end interface

      interface
      subroutine amr_restrict_ec_genorder(recv,temp,icoord,order,ivar)
      integer, intent(in)    :: icoord, order, ivar
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_ec_genorder
      end interface

      interface
      subroutine amr_restrict_ec_user()
      end subroutine amr_restrict_ec_user
      end interface

      interface
      subroutine amr_restrict_fc_fun(recv,temp,icoord)
      integer, intent(in)    :: icoord
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_fc_fun
      end interface

      interface
      subroutine amr_restrict_fc_genorder(recv,temp,icoord,order,ivar)
      integer, intent(in)    :: icoord, order, ivar
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_fc_genorder
      end interface

      interface
      subroutine amr_restrict_fc_user()
      end subroutine amr_restrict_fc_user
      end interface

      interface
      subroutine amr_restrict_red(icoord)
      integer, intent(in)    :: icoord
      end subroutine amr_restrict_red
      end interface

      interface
      subroutine amr_restrict_unk_fun(datain,dataout,lb)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: lb
      end subroutine amr_restrict_unk_fun
      end interface

      interface
      subroutine amr_restrict_unk_genorder(datain,dataout,order,ivar)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: order, ivar
      end subroutine amr_restrict_unk_genorder
      end interface

      interface
      subroutine amr_restrict_unk_user()
      end subroutine amr_restrict_unk_user
      end interface

      interface
      subroutine amr_restrict_nc_fun(datain,dataout)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      end subroutine amr_restrict_nc_fun
      end interface

      interface
      subroutine amr_restrict_nc_genorder(datain,dataout,ivar)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: ivar
      end subroutine amr_restrict_nc_genorder
      end interface

      interface
      subroutine amr_restrict_nc_user()
      end subroutine amr_restrict_nc_user
      end interface

      interface
      subroutine amr_restrict_work_fun(datain,dataout,iopt)
      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)
      integer, intent(in) :: iopt
      end subroutine amr_restrict_work_fun
      end interface

      interface
      subroutine amr_restrict_work_genorder(datain,dataout,iopt,order)
      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)
      integer, intent(in) :: iopt, order
      end subroutine amr_restrict_work_genorder
      end interface

      interface
      subroutine amr_restrict_work_user()
      end subroutine amr_restrict_work_user
      end interface

      interface
      subroutine amr_restrict_work_fun_recip(datain,dataout)
      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)
      end subroutine amr_restrict_work_fun_recip
      end interface


      interface
      subroutine amr_ser_distribute (nprocs,mype,lnblocks_old)
      integer, intent(in)  ::  nprocs,mype,lnblocks_old
      end subroutine amr_ser_distribute
      end interface


      interface
      subroutine amr_sort_by_work (new_loc,nprocs,mype)
      integer, intent(in)    :: mype
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in)    ::  nprocs
      end subroutine amr_sort_by_work
      end interface

      interface
      subroutine morton_sort(mort_no,ix,iend)
      integer, intent(in) :: iend
      integer, intent(inout) :: mort_no(6,iend), ix(iend)
      end subroutine morton_sort
      end interface

      interface
      subroutine amr_sort_morton (mort_no,new_loc,nprocs)
      integer, intent(inout) ::  mort_no(:,:)
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in)    ::  nprocs
      end subroutine amr_sort_morton
      end interface

      interface
      subroutine amr_surrounding_blks(mype,pe,lb,surrblks,ldiag)
      integer, intent(in)    ::  mype,pe,lb
      integer, intent(inout) ::  surrblks(2,3,3,3)
      logical, intent(in)    ::  ldiag
      end subroutine amr_surrounding_blks
      end interface

      interface
      subroutine amr_test_refinement(mype,lrefine_min,lrefine_max)
      integer, intent(in)    ::  mype,lrefine_min,lrefine_max
      end subroutine amr_test_refinement
      end interface

      interface
      subroutine comm_finish
      end subroutine comm_finish
      end interface

      interface
      subroutine comm_start(MaxProcs,nprocs,mype)
      integer, intent(out) :: nprocs,mype
      integer, intent(in)  :: MaxProcs
      end subroutine comm_start
      end interface

      interface
      subroutine comm_int_sum_to_all(target,source)
      integer, intent(in) :: source
      integer, intent(out)  :: target
      end subroutine comm_int_sum_to_all
      end interface

      interface
      subroutine comm_int_min_to_all(target,source)
      integer, intent(in) :: source
      integer, intent(out)  :: target
      end subroutine comm_int_min_to_all
      end interface

      interface
      subroutine comm_int_max_to_all(target,source)
      integer, intent(in) :: source
      integer, intent(out)  :: target
      end subroutine comm_int_max_to_all
      end interface


      interface
      subroutine comm_real_sum_to_all(target,source)
      real, intent(in) :: source
      real, intent(out)  :: target
      end subroutine comm_real_sum_to_all
      end interface

      interface
      subroutine comm_real_min_to_all(target,source)
      real, intent(in) :: source
      real, intent(out)  :: target
      end subroutine comm_real_min_to_all
      end interface

      interface
      subroutine comm_real_max_to_all(target,source)
      real, intent(in) :: source
      real, intent(out)  :: target
      end subroutine comm_real_max_to_all
      end interface

      interface
      subroutine test_neigh_data(mype,istep)
      integer, intent(in)    ::  mype,istep
      end subroutine test_neigh_data
      end interface

      interface
      subroutine fill_old_loc(new_loc,old_loc,nprocs,mype)
      integer, intent(in)    :: mype,nprocs
      integer, intent(inout) :: new_loc(:,:)
      integer, intent(out)   :: old_loc(:,:)
      end subroutine fill_old_loc
      end interface

      interface
      subroutine gtest_neigh_data(mype,istep,test_a)
      integer, intent(in)    ::  mype,istep
      real,    intent(in)    ::  test_a
      end subroutine gtest_neigh_data
      end interface

      interface
      subroutine mesh_test(mype)
      integer, intent(in)    ::  mype
      end subroutine mesh_test
      end interface

      interface
      subroutine guardcell_test(mype)
      integer, intent(in)    ::  mype
      end subroutine guardcell_test
      end interface

      interface 
      subroutine init_sparse_solver
      end subroutine init_sparse_solver
      end interface

      interface
      subroutine amr_1blk_fc_clean_divb(                                 & 
     &        nfacevar_in,                                               & 
     &        ia,ib,ja,jb,ka,kb,                                         & 
     &        ionea,ioneb,                                               & 
     &        jonea,joneb,                                               & 
     &        konea,koneb,                                               & 
     &        idest,ioff,joff,koff,                                      & 
     &        mype,lb,parent_pe,parent_blk )
      integer, intent(in) :: nfacevar_in
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: ionea,ioneb
      integer, intent(in) :: jonea,joneb
      integer, intent(in) :: konea,koneb
      integer, intent(in) :: idest, ioff, joff, koff
      integer, intent(in) :: mype, lb, parent_pe, parent_blk
      end subroutine amr_1blk_fc_clean_divb 
      end interface

      interface
      subroutine prol_fc_clean_divb_test(flag)
      logical, intent(in) :: flag
      end subroutine prol_fc_clean_divb_test 
      end interface

      interface
      subroutine prol_fc_clean_divb_test_report(nerrors)
      integer, intent(inout) :: nerrors
      end subroutine prol_fc_clean_divb_test_report
      end interface

      interface 
      subroutine amr_q_sort (ix,n,ia,ib)
      integer, intent(in) :: n
      integer, dimension(n),  intent(inout) :: ix
      integer, optional, dimension(n), intent(inout) :: ia, ib
      end subroutine amr_q_sort
      end interface

      end module paramesh_interfaces
