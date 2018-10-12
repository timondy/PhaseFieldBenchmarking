!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_prolong
!! NAME
!!
!!   amr_prolong
!! 
!! SYNOPSIS
!!
!!   call amr_prolong (mype, iopt, nlayers)
!!
!!   call amr_prolong (integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: mype     
!!      Current processor number
!!
!!   integer, intent(in) :: iopt     
!!      Switch to select which datastructures are updated. If iopt=1 
!!      Then this routine acts on 'unk', 'facevarx(y,z)', 'unk_e_x(y,z)', 
!!      and 'unk_n'.
!!      If iopt=2 only 'work' is updated.
!!
!!   integer, intent(in) :: nlayers 
!!      Number of layers of guard cells at a block boundary.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   workspace
!!   Mpiv_morton
!!   paramesh_mpi_interfaces
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_1blk_cc_prol_gen_work_fun
!!   amr_1blk_fc_prol_gen_fun
!!   amr_1blk_ec_prol_gen_fun
!!   amr_1blk_nc_prol_gen_fun
!!   amr_1blk_copy_soln
!!   amr_1blk_guardcell_reset
!!   amr_1blk_guardcell
!!   amr_1blk_cc_prol_gen_unk_fun
!!   comm_int_max_to_all
!!   comm_int_min_to_all
!!   amr_1blk_fc_prol_dbz
!!   mpi_amr_comm_setup
!!
!! RETURNS
!!
!!   Upon exit, prolongation (i.e. interpolation) from parent blocks to
!!   their newly created child (marked by the 'newchild' flag) has been
!!   performed.
!!
!! DESCRIPTION
!!
!!   This routine interpolates data from parent blocks to any newly created 
!!   child blocks which are marked with the 'newchild' flag set to true.
!!
!! AUTHORS
!!
!!   Peter MacNeice (July 1997)
!!
!!   Modified by Michael L. Rilee, November 2002, *dbz*
!!        Initial support for divergenceless prolongation
!!   Modified by Michael L. Rilee, December 2002, *clean_divb*
!!        Support for projecting field onto divergenceless field
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_prolong(mype,iopt,nlayers)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use mpi_morton
      Use prolong_arrays, Only :                                       & 
           prol_fc_dbz,                                                & 
           prol_fc_dbz_ivar,                                           & 
           prol_fc_dbz_n,                                              & 
           prol_fc_clean_divb
      Use paramesh_interfaces, Only :                                  & 
                        amr_1blk_cc_prol_gen_work_fun,                 & 
                        amr_1blk_fc_prol_gen_fun,                      & 
                        amr_1blk_ec_prol_gen_fun,                      & 
                        amr_1blk_nc_prol_gen_fun,                      & 
                        amr_1blk_copy_soln,                            & 
                        amr_1blk_guardcell_reset,                      & 
                        amr_1blk_guardcell,                            & 
                        amr_1blk_cc_prol_gen_unk_fun,                  & 
                        comm_int_max_to_all,                           & 
                        comm_int_min_to_all,                           & 
                        amr_1blk_fc_prol_dbz,                          &
                        amr_q_sort
      Use paramesh_mpi_interfaces, Only :                              & 
                        mpi_amr_comm_setup

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in) ::  mype,iopt,nlayers

!-----Local arrays and variables.
      Real :: recvf(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,     & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real :: recve(nbndvare,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,    & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real :: recvn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,    & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real :: recvfx(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,        & 
                                            kl_bnd1:ku_bnd1)
      Real :: recvfy(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,      & 
                                            kl_bnd1:ku_bnd1)
      Real :: recvfz(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,          & 
                                            kl_bnd1:ku_bnd1+k3d)
      Integer :: p_cache_addr(2),idest,nguard0,nguard_npgs
      Integer :: parent_blk,parent_pe
!      Integer :: parent_list(2,maxblocks)
!      Integer :: indx(maxblocks),index(maxblocks)
!      Integer :: parent_id(maxblocks),icoord
      Integer,allocatable :: parent_list(:,:)
      Integer,allocatable :: indx(:),index(:)
      Integer,allocatable :: parent_id(:)
      Integer :: icoord
      Integer,Save :: lref_min,lref_max
      Integer,Save :: lref_mint,lref_maxt
      Integer,Save :: nnewchildg,nnewchild
      Integer :: ia,ib,ja,jb,ka,kb,iblock,lb,lreflevel,lbi
      Integer :: ioff,joff,koff ,k,i
      Integer :: tag_offset
      Integer :: p_blk,p_pe,iblk
      Integer :: nprocs, ierr
      Integer :: level,j
      Integer :: nfield
      Integer :: imlr, jmlr, kmlr,idim_mlr
      Integer :: iprol, iv1, iv2, iv3
      Logical :: lnewchild, lfound
      Logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
!-----Begin executable code.
! CEG allocate memory
      allocate(parent_list(2,maxblocks))
      allocate(indx(maxblocks))
      allocate(index(maxblocks))
      allocate(parent_id(maxblocks))

!------set state flag
      lprolong_in_progress = .True.

!------Are there any new children?
      nnewchild = 0
      lnewchild = any(newchild)
      If (lnewchild) nnewchild = 1
      Call comm_int_max_to_all (nnewchildg,nnewchild)

      If (nnewchildg == 0) Return

      Call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!------reset cache addresses
      Call amr_1blk_guardcell_reset

!------Identify variables to be prolonged
      lcc = .False.
      lfc = .False.
      lec = .False.
      lnc = .False.
      If (iopt == 1) Then
         If (nvar > 0)     lcc = .True.
         If (nfacevar > 0) lfc = .True.
         If (nvaredge > 0) lec = .True.
         If (nvarcorn > 0) lnc = .True.
         nguard0 = nguard
      ElseIf (iopt >= 2) Then
         lcc = .True.
         nguard0 = nguard_work
      End If

      nguard_npgs = nguard*npgs

      ia = 1+nguard0       
      ib = nxb+nguard0       
      ja = 1+nguard0*k2d      
      jb = nyb+nguard0*k2d
      ka = 1+nguard0*k3d      
      kb = nzb+nguard0*k3d

!------construct a list of parent blocks for new children
      parent_list(:,:) = -1
      iblock=0
      lref_min = 100
      lref_max = 1
      If (lnblocks > 0) Then
         Do lb = 1,lnblocks
            If (newchild(lb)) Then
               iblock = iblock+1
               index(iblock) = lb
               indx(iblock) = iblock
               parent_list(:,iblock) = parent(:,lb)
               parent_id(iblock) = parent(1,lb)+maxblocks*parent(2,lb)
               lref_min = min(lrefine(lb),lref_min)
               lref_max = max(lrefine(lb),lref_max)
            End If
         End Do
      End If

      lref_maxt = lref_max
      lref_mint = lref_min
      Call comm_int_max_to_all (lref_max,lref_maxt)
      Call comm_int_min_to_all (lref_min,lref_mint)

      If (lref_min > lref_max) lref_min = lref_max

      If (iblock > nchild) Then
!------sort the list of newchildren according to their parents ids
!------This will enable us to avoid costly extra guardcell filling operations
!------on parent blocks
         Call amr_q_sort (parent_id,iblock,ia = indx)
      End If

      p_cache_addr(:) = -1

!------prolongation must be applied in descending order of refinement level
!------so that the case where neighboring blocks at different refinement
!------level are required to be refined at the same time is handled correctly.
      Do lreflevel = lref_min,lref_max

         Call amr_1blk_guardcell_reset
         If (no_permanent_guardcells) Then
!-------Store a copy of the current solution in gt_unk
            level = -1                      ! copy all refinement levels
            Call amr_1blk_copy_soln(level)
         End If
        
         If (.Not.no_permanent_guardcells) Then
            If (force_consistency) Then
               If (lfc) Then  
!-----If using facevars and permanent guardcells Then we will need to ensure
!-----that gt_facevarx, etc are filled before prolongation begins
                  If (nfacevar > 0) Then
                     Do lb = 1,lnblocks
                        gt_facevarx(:,1,:,:,lb) =                                    & 
                            facevarx(:,1+nguard_npgs,:,:,lb)
                        gt_facevarx(:,2,:,:,lb) =                                    & 
                            facevarx(:,nxb+1+nguard_npgs,:,:,lb)
                        If (ndim >= 2) Then
                           gt_facevary(:,:,1,:,lb) =                                    & 
                             facevary(:,:,1+nguard_npgs*k2d,:,lb)
                           gt_facevary(:,:,1+k2d,:,lb) =                                & 
                              facevary(:,:,nyb+(1+nguard_npgs)*k2d,:,lb)
                           If (ndim == 3) Then
                              gt_facevarz(:,:,:,1,lb) =                                    & 
                                 facevarz(:,:,:,1+nguard_npgs*k3d,lb)
                              gt_facevarz(:,:,:,1+k3d,lb) =                                & 
                                 facevarz(:,:,:,nzb+(1+nguard_npgs)*k3d,lb)
                           End If
                        End If
                     End Do
                  End If  ! End If (nfacevar > 0)
               End If  ! End If (lfc)
            End If  ! End If (force_consistency)
         End If  ! End If (.Not.no_permanent_guardcells)

!------moved next block inside loop over refinement levels to ensure
!------that when a block refines while its neighbor is also refining, the
!------more refined of the pre-existing blocks will get good guardcell data.
         tag_offset = 100

         lguard    = .False.
         lprolong  = .True.
         lflux     = .False.
         ledge     = .False.
         lrestrict = .False.
         lfulltree = .False.
         Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,            & 
                               lflux,ledge,lrestrict,lfulltree,        & 
                               iopt,lcc,lfc,lec,lnc,tag_offset)

         If (iblock > 0) Then
            Do lbi = 1,iblock
               lb = index(indx(lbi))

               If (lrefine(lb) == lreflevel) Then

!------compute offset for child cell inside parent
                  ioff = mod(which_child(lb)-1,2)*nxb/2
                  joff = mod((which_child(lb)-1)/2,2)*nyb/2
                  koff = mod((which_child(lb)-1)/4,2)*nzb/2

!------get address of parent block
                  parent_blk = parent(1,lb)
                  parent_pe  = parent(2,lb)

!------Is parent data currently cached?
                  If ( parent_blk.ne.p_cache_addr(1) .Or.                         & 
                      parent_pe.ne.p_cache_addr(2) ) Then

!--------Fetch parent data block into layer 1 of 1blk data structure and fill
!--------its guardcells 
                     ldiag = diagonals
                     l_srl_only = .True.           ! seems to be all that is required ??
                     icoord = 0

!--------If (parent_blk,parent_pe) is not a local block Then it must have a 
!--------local copy available in the buffer space at the end of the local
!--------block list.
                     p_blk = parent_blk
                     p_pe  = parent_pe
                     If (parent_pe.ne.mype) Then

                        lfound = .False.
                        iblk = ladd_strt(parent_pe)
                        Do while(.Not.lfound.And.                                   & 
                           iblk <= ladd_end(parent_pe))
                           If (parent_blk == laddress(1,iblk).And.                   & 
                              parent_pe  == laddress(2,iblk) ) Then
                              p_blk = iblk
                              p_pe  = mype
                              lfound = .True.
                           Else
                              iblk = iblk+1
                           End If
                        End Do
                     End If
!print *, mype,'prolong lb=',lb,p_blk, p_pe
                     Call amr_1blk_guardcell(mype,iopt,nlayers,p_blk,              & 
                                 p_pe,lcc,lfc,lec,lnc,l_srl_only,      & 
                                 icoord,ldiag)

!--------update address of cached parent
                     p_cache_addr(1) = parent_blk
                     p_cache_addr(2) = parent_pe

                  End If  ! End If If ( parent_blk.ne.p_cache_addr(1)


!-----Prolong data from working block to the new child block
      idest = 2

!-----cell-centered data
      If (lcc) Then
         If (iopt == 1) Then
           Call amr_1blk_cc_prol_gen_unk_fun(                          & 
                    unk1(:,:,:,:,1),                                   & 
                    ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,mype,       & 
                    lb,parent_pe,parent_blk)
         Else If (iopt >= 2) Then
          Call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,1),          & 
                    ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,mype,       & 
                    lb,parent_pe,parent_blk,interp_mask_work(iopt-1))
         End If                      
      End If 

!------face-centered data 
      If (lfc) Then 
          If ( prol_fc_dbz) Then
!------------prolong divergenceless b on fc
             recvfx(1:nfacevar, il_bnd1:iu_bnd1+1,                     & 
                  jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)                     & 
               = facevarx1(1:nfacevar, il_bnd1:iu_bnd1+1,              & 
                  jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1, 1)

             recvfy(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,   & 
                  kl_bnd1:ku_bnd1)                                     & 
               = facevary1(1:nfacevar, il_bnd1:iu_bnd1,                & 
                  jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)

             recvfz(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,       & 
                  kl_bnd1:ku_bnd1+k3d)                                 & 
               = facevarz1(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,& 
                  kl_bnd1:ku_bnd1+k3d, 1)
             Do iprol = 1, prol_fc_dbz_n
                iv1 = prol_fc_dbz_ivar(1,iprol)
                iv2 = prol_fc_dbz_ivar(2,iprol)
                iv3 = prol_fc_dbz_ivar(3,iprol)
                Call amr_1blk_fc_prol_dbz(                             & 
                  recvfx, recvfy, recvfz,                              & 
                  nfacevar, iv1, iv2, iv3,                             & 
                  ia,ib,ja,jb,ka,kb,                                   &  
                  idest,ioff,joff,koff,                                & 
                  mype,lb,parent_pe,parent_blk                         & 
                  )
             End Do
          
         Else
!--------prolong using other methods
!--------x-face
         recvf(1:nfacevar, il_bnd1:iu_bnd1+1,                          & 
                     jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)                  & 
            = facevarx1(1:nfacevar, il_bnd1:iu_bnd1+1,                 & 
                     jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1, 1)
         Call amr_1blk_fc_prol_gen_fun(recvf, & 
                    ia,ib+1,ja,jb,ka,kb,idest,ioff,joff,koff,          & 
                    mype,lb,parent_pe,parent_blk,1)

         If (ndim > 1) Then
!--------y-face
         recvf(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,        & 
                    kl_bnd1:ku_bnd1)                                   & 
            = facevary1(1:nfacevar, il_bnd1:iu_bnd1,                   & 
                   jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)
         Call amr_1blk_fc_prol_gen_fun(recvf,                          & 
                    ia,ib,ja,jb+1,ka,kb,idest,ioff,joff,koff,          & 
                    mype,lb,parent_pe,parent_blk,2)
         End If  ! End If (ndim > 1)

         If (ndim == 3) Then
!--------z-face
         recvf(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,            & 
                    kl_bnd1:ku_bnd1+k3d)                               & 
            = facevarz1(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,   & 
                    kl_bnd1:ku_bnd1+k3d, 1)
         Call amr_1blk_fc_prol_gen_fun(recvf,                          & 
                    ia,ib,ja,jb,ka,kb+1,idest,ioff,joff,koff,          & 
                    mype,lb,parent_pe,parent_blk,3)

         End If  ! End If (ndim == 3)

!------clean_divb
       If (prol_fc_clean_divb)Then
          Call amr_1blk_fc_clean_divb(                                 & 
               nfacevar,                                               & 
               ia,ib,ja,jb,ka,kb,                                      & 
               0, 0, 0, 0, 0, 0,                                       & 
               idest,ioff,joff,koff,                                   & 
               mype,lb,parent_pe,parent_blk                            & 
               )
       End If

      End If  ! End If (prol_fc_dbz)

      End If  ! End If (lfc)

       If (ndim > 1) Then
!------edge-centered data
       If (lec) Then
!--------x-edge
         recve(1:nvaredge, il_bnd1:iu_bnd1,                            & 
                     jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d)          & 
            = unk_e_x1(1:nvaredge, il_bnd1:iu_bnd1,                    & 
                     jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d, 1)
        
         Call amr_1blk_ec_prol_gen_fun(recve,                          & 
                                       ia,ib,ja,jb+k2d,ka,kb+k3d,      & 
                                       idest,ioff,joff,koff,mype,1)

!--------y-edge
         recve(1:nvaredge, il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,          & 
                    kl_bnd1:ku_bnd1+k3d)                               & 
            = unk_e_y1(1:nvaredge, il_bnd1:iu_bnd1+1,                  & 
                   jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1+k3d, 1)

         Call amr_1blk_ec_prol_gen_fun(recve,                          & 
                                       ia,ib+1,ja,jb,ka,kb+k3d,        & 
                                       idest,ioff,joff,koff,mype,2)

         If (ndim == 3) Then
!--------z-edge
         recve(1:nvaredge, il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,      & 
                    kl_bnd1:ku_bnd1)                                   & 
            = unk_e_z1(1:nvaredge, il_bnd1:iu_bnd1+1,                  & 
                   jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)
         Call amr_1blk_ec_prol_gen_fun(recve,                          & 
                                       ia,ib+1,ja,jb+k2d,ka,kb,        & 
                                       idest,ioff,joff,koff,mype,3)
         End If

        End If  ! End If (lec)

        End If  ! End If (ndim > 1)

!------cell corner data
       If (lnc) Then
           recvn(:,:,:,:) = unk_n1(:,:,:,:,1)

           Call amr_1blk_nc_prol_gen_fun(recvn,                        & 
                                         ia,ib+1,ja,jb+k2d,ka,kb+k3d,  & 
                                         idest,ioff,joff,koff,mype)
       End If  ! End If (lnc)

!-----copy data back to permanent storage arrays

       Call amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest )

      End If  ! End If (lrefine(lb) == lreflevel)

      End Do  ! End Do lbi = 1,iblock
      End If  ! End If (iblock > 0)

!-----Ensure new blocks inherit data on block face shared with an old
!-----existing neighbor, instead of filling from parent by interpolation.
      If (divergence_free) Then
        If (lfc) Then
        Do nfield = 1,nfield_divf
          Call amr_prolong_fc_divbconsist(mype,lreflevel,              & 
                                          nfield)
        End Do
        End If
      End If

      End Do  ! End Do lreflevel = lref_min,lref_max

!-----reset cache addresses
      Call amr_1blk_guardcell_reset

      newchild(:) = .False.

!-----unset state flag
      lprolong_in_progress = .False.

! CEG deallocate memory
  deallocate(parent_list)
  deallocate(indx)
  deallocate(index)
  deallocate(parent_id)

      Return
      End Subroutine amr_prolong


!-----------------------------------------------------------------------------------------------------------------

      Subroutine pf_prolong(mype,iopt,nlayers,mglevel)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use mpi_morton
      Use prolong_arrays, Only :                                       & 
           prol_fc_dbz,                                                & 
           prol_fc_dbz_ivar,                                           & 
           prol_fc_dbz_n,                                              & 
           prol_fc_clean_divb
      Use paramesh_interfaces, Only :                                  & 
                        amr_1blk_cc_prol_gen_work_fun,                 & 
                        amr_1blk_fc_prol_gen_fun,                      & 
                        amr_1blk_ec_prol_gen_fun,                      & 
                        amr_1blk_nc_prol_gen_fun,                      & 
                        amr_1blk_copy_soln,                            & 
                        amr_1blk_guardcell_reset,                      & 
                        amr_1blk_guardcell,                            & 
                        amr_1blk_cc_prol_gen_unk_fun,                  & 
                        comm_int_max_to_all,                           & 
                        comm_int_min_to_all,                           & 
                        amr_1blk_fc_prol_dbz,                          &
                        amr_q_sort
      Use paramesh_mpi_interfaces, Only :                              & 
                        mpi_amr_comm_setup

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in) ::  mype,iopt,nlayers, mglevel

!-----Local arrays and variables.
      Real :: recvf(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,     & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real :: recve(nbndvare,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,    & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real :: recvn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,    & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real :: recvfx(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,        & 
                                            kl_bnd1:ku_bnd1)
      Real :: recvfy(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,      & 
                                            kl_bnd1:ku_bnd1)
      Real :: recvfz(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,          & 
                                            kl_bnd1:ku_bnd1+k3d)
      Integer :: p_cache_addr(2),idest,nguard0,nguard_npgs
      Integer :: parent_blk,parent_pe
!      Integer :: parent_list(2,maxblocks)
!      Integer :: indx(maxblocks),index(maxblocks)
!      Integer :: parent_id(maxblocks),icoord
      Integer,allocatable :: parent_list(:,:)
      Integer,allocatable :: indx(:),index(:)
      Integer,allocatable :: parent_id(:)
      Integer :: icoord
      Integer,Save :: lref_min,lref_max
      Integer,Save :: lref_mint,lref_maxt
      Integer,Save :: nnewchildg,nnewchild
      Integer :: ia,ib,ja,jb,ka,kb,iblock,lb,lreflevel,lbi
      Integer :: ioff,joff,koff ,k,i
      Integer :: tag_offset
      Integer :: p_blk,p_pe,iblk
      Integer :: nprocs, ierr
      Integer :: level,j
      Integer :: nfield
      Integer :: imlr, jmlr, kmlr,idim_mlr
      Integer :: iprol, iv1, iv2, iv3
      Logical :: lnewchild, lfound
      Logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree

!-----Begin executable code.
! CEG allocate memory
      allocate(parent_list(2,maxblocks))
      allocate(indx(maxblocks))
      allocate(index(maxblocks))
      allocate(parent_id(maxblocks))

!------set state flag
      lprolong_in_progress = .True.

!------Are there any new children?
      nnewchild = 0
      lnewchild = any(newchild)
      If (lnewchild) nnewchild = 1
      Call comm_int_max_to_all (nnewchildg,nnewchild)
      If (nnewchildg == 0) Return

      Call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!------reset cache addresses
      Call amr_1blk_guardcell_reset

!------Identify variables to be prolonged
      lcc = .False.
      lfc = .False.
      lec = .False.
      lnc = .False.
      If (iopt == 1) Then
         If (nvar > 0)     lcc = .True.
         If (nfacevar > 0) lfc = .True.
         If (nvaredge > 0) lec = .True.
         If (nvarcorn > 0) lnc = .True.
         nguard0 = nguard
      ElseIf (iopt >= 2) Then
         lcc = .True.
         nguard0 = nguard_work
      End If

      nguard_npgs = nguard*npgs

      ia = 1+nguard0       
      ib = nxb+nguard0       
      ja = 1+nguard0*k2d      
      jb = nyb+nguard0*k2d
      ka = 1+nguard0*k3d      
      kb = nzb+nguard0*k3d

!------construct a list of parent blocks for new children
      parent_list(:,:) = -1
      iblock=0
      lref_min = 100
      lref_max = 1
      If (lnblocks > 0) Then
         Do lb = 1,lnblocks
            If (newchild(lb)) Then
               iblock = iblock+1
               index(iblock) = lb
               indx(iblock) = iblock
               parent_list(:,iblock) = parent(:,lb)
               parent_id(iblock) = parent(1,lb)+maxblocks*parent(2,lb)
               lref_min = min(lrefine(lb),lref_min)
               lref_max = max(lrefine(lb),lref_max)
            End If
         End Do
      End If

      lref_maxt = lref_max
      lref_mint = lref_min
      Call comm_int_max_to_all (lref_max,lref_maxt)
      Call comm_int_min_to_all (lref_min,lref_mint)

      If (lref_min > lref_max) lref_min = lref_max

      If (iblock > nchild) Then
!------sort the list of newchildren according to their parents ids
!------This will enable us to avoid costly extra guardcell filling operations
!------on parent blocks
         Call amr_q_sort (parent_id,iblock,ia = indx)
      End If

      p_cache_addr(:) = -1

!------prolongation must be applied in descending order of refinement level
!------so that the case where neighboring blocks at different refinement
!------level are required to be refined at the same time is handled correctly.
      Do lreflevel = lref_min,lref_max

         Call amr_1blk_guardcell_reset
         If (no_permanent_guardcells) Then
!-------Store a copy of the current solution in gt_unk
            level = -1                      ! copy all refinement levels
            Call amr_1blk_copy_soln(level)
         End If

        
!------moved next block inside loop over refinement levels to ensure
!------that when a block refines while its neighbor is also refining, the
!------more refined of the pre-existing blocks will get good guardcell data.
         tag_offset = 100

         lguard    = .False.
         lprolong  = .True.
         lflux     = .False.
         ledge     = .False.
         lrestrict = .False.
         lfulltree = .False.
         Call pf_mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,            & 
                               lflux,ledge,lrestrict,lfulltree,        & 
                               iopt,lcc,lfc,lec,lnc,tag_offset, mglevel)

         If (iblock > 0) Then
            Do lbi = 1,iblock
               lb = index(indx(lbi))

               If (lrefine(lb) == lreflevel) Then

!------compute offset for child cell inside parent
                  ioff = mod(which_child(lb)-1,2)*nxb/2
                  joff = mod((which_child(lb)-1)/2,2)*nyb/2
                  koff = mod((which_child(lb)-1)/4,2)*nzb/2

!------get address of parent block
                  parent_blk = parent(1,lb)
                  parent_pe  = parent(2,lb)

!------Is parent data currently cached?
                  If ( parent_blk.ne.p_cache_addr(1) .Or.                         & 
                      parent_pe.ne.p_cache_addr(2) ) Then

!--------Fetch parent data block into layer 1 of 1blk data structure and fill
!--------its guardcells 
                     ldiag = diagonals
                     l_srl_only = .True.           ! seems to be all that is required ??
                     icoord = 0

!--------If (parent_blk,parent_pe) is not a local block Then it must have a 
!--------local copy available in the buffer space at the end of the local
!--------block list.
                     p_blk = parent_blk
                     p_pe  = parent_pe
                     If (parent_pe.ne.mype) Then

                        lfound = .False.
                        iblk = ladd_strt(parent_pe)
                        Do while(.Not.lfound.And.                                   & 
                           iblk <= ladd_end(parent_pe))
                           If (parent_blk == laddress(1,iblk).And.                   & 
                              parent_pe  == laddress(2,iblk) ) Then
                              p_blk = iblk
                              p_pe  = mype
                              lfound = .True.
                           Else
                              iblk = iblk+1
                           End If
                        End Do
                     End If
!print *, mype,'prolong lb=',lb,p_blk, p_pe
                     Call amr_1blk_guardcell(mype,iopt,nlayers,p_blk,              & 
                                 p_pe,lcc,lfc,lec,lnc,l_srl_only,      & 
                                 icoord,ldiag)

!--------update address of cached parent
                     p_cache_addr(1) = parent_blk
                     p_cache_addr(2) = parent_pe

                  End If  ! End If If ( parent_blk.ne.p_cache_addr(1)


!-----Prolong data from working block to the new child block
      idest = 2

!-----cell-centered data
      If (lcc) Then
         If (iopt == 1) Then
           Call amr_1blk_cc_prol_gen_unk_fun(                          & 
                    unk1(:,:,:,:,1),                                   & 
                    ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,mype,       & 
                    lb,parent_pe,parent_blk)
         ElseIf (iopt >= 2) Then
          Call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,1),          & 
                    ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,mype,       & 
                    lb,parent_pe,parent_blk,interp_mask_work(iopt-1))
         End If                      
      End If 

!-----copy data back to permanent storage arrays

       Call amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest )

      End If  ! End If (lrefine(lb) == lreflevel)

      End Do  ! End Do lbi = 1,iblock
      End If  ! End If (iblock > 0)

      End Do  ! End Do lreflevel = lref_min,lref_max

!-----reset cache addresses
      Call amr_1blk_guardcell_reset

      newchild(:) = .False.

!-----unset state flag
      lprolong_in_progress = .False.

! CEG deallocate memory
  deallocate(parent_list)
  deallocate(indx)
  deallocate(index)
  deallocate(parent_id)

      Return
      End Subroutine pf_prolong
