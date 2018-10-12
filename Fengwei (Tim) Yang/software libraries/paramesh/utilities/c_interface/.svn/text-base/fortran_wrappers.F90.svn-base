!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_initialize()
      
      use paramesh_interfaces, only : amr_initialize
      use io
      include 'mpif.h'

      integer :: ierr
      character (len=80) ::stdout

      call amr_initialize()

      stdout = trim(output_dir) // 'stdout'
      open(unit=6,file=stdout,form='formatted',status='unknown')

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end subroutine c_amr_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_close()

      use paramesh_interfaces, only : amr_close

      close(6)
      call amr_close()

      stop
      end subroutine c_amr_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_1blk_copy_soln(level)

      use paramesh_interfaces, only : amr_1blk_copy_soln

      integer ::  level

      call amr_1blk_copy_soln(level)

      return
      end subroutine c_amr_1blk_copy_soln

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_1blk_guardcell (mype,iopt,nlayers,lb,pe, & 
     &                                 lcc,lfc,lec,lnc, & 
     &                                 l_srl_only,icoord,ldiag, & 
     &                                 nlayersx, nlayersy, nlayersz)

      use paramesh_interfaces, only : amr_1blk_guardcell

      integer :: mype,iopt,nlayers,lb,pe,icoord
      logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      integer :: nlayersx,nlayersy,nlayersz
      integer :: lbb

      lbb = lb + 1

      call amr_1blk_guardcell(mype,iopt,nlayers,lbb,pe, & 
     &                        lcc,lfc,lec,lnc, & 
     &                        l_srl_only,icoord,ldiag, & 
     &                        nlayersx, nlayersy, nlayersz)

      return
      end subroutine c_amr_1blk_guardcell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine c_amr_1blk_guardcell_reset()

      use paramesh_interfaces, only : amr_1blk_guardcell_reset

      call amr_1blk_guardcell_reset()

      return
      end subroutine c_amr_1blk_guardcell_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_1blk_save_soln()

      use paramesh_interfaces, only : amr_1blk_save_soln

      call amr_1blk_save_soln()

      return
      end subroutine c_amr_1blk_save_soln

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_1blk_to_perm(lcc,lfc,lec,lnc,lb,iopt,idest)

      use paramesh_interfaces, only : amr_1blk_to_perm

      integer :: lb,iopt,idest
      logical :: lcc,lfc,lec,lnc
      integer :: lbb, iidest

      lbb = lb + 1
      iidest = idest + 1

      call amr_1blk_to_perm(lcc,lfc,lec,lnc,lbb,iopt,iidest)

      return
      end subroutine c_amr_1blk_to_perm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_block_geometry(lb,pe)

      use paramesh_interfaces, only : amr_block_geometry

      integer :: lb, pe
      integer :: lbb

      lbb = lb + 1

      call amr_block_geometry(lbb,pe)

      return
      end subroutine c_amr_block_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_checkpoint_re(iunit1, & 
     &                               l_with_guardcells, & 
     &                               check_format, & 
     &                               user_attr_1, & 
     &                               user_attr_2, & 
     &                               user_attr_3, & 
     &                               user_attr_4, & 
     &                               user_attr_5)

      use paramesh_interfaces, only : amr_checkpoint_re

      integer            :: iunit1
      logical            :: l_with_guardcells
      character (len=80) :: check_format
      real               :: user_attr_1, & 
     &                      user_attr_2, & 
     &                      user_attr_3, & 
     &                      user_attr_4, & 
     &                      user_attr_5

      call amr_checkpoint_re(iunit1, & 
     &                       l_with_guardcells, & 
     &                       check_format, & 
     &                       user_attr_1, & 
     &                       user_attr_2, & 
     &                       user_attr_3, & 
     &                       user_attr_4, & 
     &                       user_attr_5)


      return
      end subroutine c_amr_checkpoint_re

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_checkpoint_wr(iunit1, & 
     &                               l_with_guardcells, & 
     &                               check_format, & 
     &                               user_attr_1, & 
     &                               user_attr_2, & 
     &                               user_attr_3, & 
     &                               user_attr_4, & 
     &                               user_attr_5)

      use paramesh_interfaces, only : amr_checkpoint_wr

      integer            :: iunit1
      logical            :: l_with_guardcells
      character (len=80) :: check_format
      real               :: user_attr_1, & 
     &                      user_attr_2, & 
     &                      user_attr_3, & 
     &                      user_attr_4, & 
     &                      user_attr_5

      call amr_checkpoint_wr(iunit1, & 
     &                       l_with_guardcells, & 
     &                       check_format, & 
     &                       user_attr_1, & 
     &                       user_attr_2, & 
     &                       user_attr_3, & 
     &                       user_attr_4, & 
     &                       user_attr_5)

      return
      end subroutine c_amr_checkpoint_wr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_plotfile_chombo(iunit1)

      use paramesh_interfaces, only : amr_plotfile_chombo

      integer :: iunit1

      call amr_plotfile_chombo(iunit1)

      return
      end subroutine c_amr_plotfile_chombo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_edge_average(mype,lfullblock,nsub)

      use paramesh_interfaces, only : amr_edge_average

      integer  ::  mype,nsub
      logical  ::  lfullblock

      call amr_edge_average(mype,lfullblock,nsub)

      return
      end subroutine c_amr_edge_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_flush(iunit)
      
      use paramesh_interfaces, only : amr_flush

      integer :: iunit

      call amr_flush(iunit)

      return
      end subroutine c_amr_flush

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine c_amr_flux_conserve(mype,nsub,flux_dir)

      use paramesh_interfaces, only : amr_flux_conserve

      integer  ::  flux_dir
      integer  ::  mype,nsub

      call amr_flux_conserve(mype,nsub,flux_dir)

      return
      end subroutine c_amr_flux_conserve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_guardcell(mype,iopt,nlayers, & 
     &                           nlayersx,nlayersy,nlayersz)

      use paramesh_interfaces, only : amr_guardcell

      integer :: nlayersx,nlayersy,nlayersz
      integer :: mype,iopt,nlayers

      call amr_guardcell(mype,iopt,nlayers, & 
     &                   nlayersx,nlayersy,nlayersz)

      return
      end subroutine c_amr_guardcell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_perm_to_1blk(lcc,lfc,lec,lnc,lb,pe,iopt,idest)

      use paramesh_interfaces, only : amr_perm_to_1blk

      integer ::  lb,pe,iopt,idest
      logical ::  lcc,lfc,lec,lnc
      integer ::  lbb, iidest

      lbb = lb + 1
      iidest = idest + 1

      call amr_perm_to_1blk(lcc,lfc,lec,lnc,lbb,pe,iopt,iidest)

      return
      end subroutine c_amr_perm_to_1blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_prolong(mype,iopt,nlayers)

      use paramesh_interfaces, only : amr_prolong

      integer :: mype,iopt,nlayers

      call amr_prolong(mype,iopt,nlayers)

      return
      end subroutine c_amr_prolong

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_refine_derefine()

      use paramesh_interfaces, only : amr_refine_derefine

      call amr_refine_derefine()

      return
      end subroutine c_amr_refine_derefine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_restrict(mype,iopt,iempty,filling_guardcells)

      use paramesh_interfaces, only : amr_restrict

      integer :: mype,iopt,iempty
      logical :: filling_guardcells

      call amr_restrict(mype,iopt,iempty,filling_guardcells)

      return
      end subroutine c_amr_restrict

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_amr_1blk_fc_clean_divb( & 
     &        nfacevar_in, & 
     &        ia,ib,ja,jb,ka,kb, & 
     &        ionea,ioneb, & 
     &        jonea,joneb, & 
     &        konea,koneb, & 
     &        idest,ioff,joff,koff, & 
     &        mype,lb,parent_pe,parent_blk )

      use paramesh_interfaces, only : amr_1blk_fc_clean_divb

      integer :: nfacevar_in
      integer :: ia,ib,ja,jb,ka,kb
      integer :: ionea,ioneb
      integer :: jonea,joneb
      integer :: konea,koneb
      integer :: idest, ioff, joff, koff
      integer :: mype, lb, parent_pe, parent_blk
      integer :: lbb, iidest

      lbb = lb + 1
      iidest = idest + 1

      call amr_1blk_fc_clean_divb( & 
     &        nfacevar_in, & 
     &        ia,ib,ja,jb,ka,kb, & 
     &        ionea,ioneb, & 
     &        jonea,joneb, & 
     &        konea,koneb, & 
     &        iidest,ioff,joff,koff, & 
     &        mype,lbb,parent_pe,parent_blk )

      return
      end subroutine c_amr_1blk_fc_clean_divb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine c_mpi_amr_comm_setup(mype,nprocs,lguard,lprolong, & 
     &                                lflux,ledge,lrestrict,lfulltree, & 
     &                                iopt,lcc,lfc,lec,lnc,tag_offset, & 
     &                                nlayersx,nlayersy,nlayersz, & 
     &                                flux_dir)

      use paramesh_mpi_interfaces, only : mpi_amr_comm_setup

      integer :: mype,nprocs,iopt
      integer :: tag_offset
      logical :: lcc,lfc,lec,lnc,lfulltree
      logical :: lguard,lprolong,lflux,ledge,lrestrict
      integer :: nlayersx,nlayersy,nlayersz
      integer :: flux_dir

      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong, & 
     &                        lflux,ledge,lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset, & 
     &                        nlayersx,nlayersy,nlayersz, & 
     &                        flux_dir)

      return
      end subroutine c_mpi_amr_comm_setup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine c_mpi_morton_bnd_prolong(mype, nprocs, tag_offset)

      use paramesh_mpi_interfaces, only : mpi_morton_bnd_prolong

      integer :: mype, nprocs, tag_offset

      call mpi_morton_bnd_prolong(mype, nprocs, tag_offset)

      return
      end subroutine c_mpi_morton_bnd_prolong




