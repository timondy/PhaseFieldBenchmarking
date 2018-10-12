!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_1blk_guardcell
!! NAME
!!
!!   mpi_amr_1blk_restrict
!!
!! SYNOPSIS
!!
!!   call mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc,
!!                              lfulltree,filling_guardcells)
!!
!!   call mpi_amr_1blk_restrict(integer, integer, logical, logical, logical, logical 
!!                              logical, logical)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype           
!!        The local processor number.
!!
!!   integer, intent(in) :: iopt           
!!        A switch to control which data source is to be used:
!!         iopt=1 will use 'unk', 'facevarx', 'facevary', 'facevarz', 
!!                'unk_e_x', 'unk_e_y', 'unk_e_z', and 'unk_n'
!!         iopt>=2 will use 'work'
!!
!!   logical, intent(in) :: lcc, lfc, lec, lnc
!!        Logical switches which indicate which data is to be restricted.
!!         lcc -> cell centered
!!         lfc -> face centered
!!         lec -> edge centered
!!         lnc -> node centered
!!
!!   logical, intent(in) :: lfulltree
!!        A switch to indicate if the entire tree is to be restricted.  If true the
!!        data is restricted from the leaves of the tree all the way to the root at
!!        level = 1.  Otherwise, only the data from the leaves to their parents are
!!        restricted.
!!
!!   logical, intent(in) :: filling_guardcells
!!        A logical switch.  If true then this routine has been called as part of 
!!        the guardcell filling step so that only those blocks at jumps in refinement
!!        need to restrict data to their parents.  This leads to a performance
!!        savings.
!!        
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   workspace
!!   mpi_morton
!!   timings
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_1blk_copy_soln,
!!   amr_1blk_guardcell_reset,
!!   amr_perm_to_1blk,
!!   amr_1blk_guardcell,
!!   amr_restrict_unk_fun,
!!   amr_restrict_nc_fun,
!!   amr_restrict_fc_fun,
!!   amr_restrict_ec_fun,
!!   amr_restrict_work_fun,
!!   amr_restrict_work_fun_recip,
!!   amr_1blk_nc_cp_remote,
!!   comm_int_max_to_all,
!!   comm_int_min_to_all
!!   amr_block_geometry
!!   mpi_amr_comm_setup
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit, data has been restricted from the leaves
!!   of the data tree to their parents.
!!
!! DESCRIPTION
!!
!!   This routine does the data averaging required when a child block
!!   passes data back to its parent. The parent receives interior data
!!   only, not guard cell data. 
!!   The parent gets data for each child and then applies
!!   the restriction operator to it. Guardcell data may be needed for the
!!   child blocks, depending on the particular restriction operator being used.
!!   Thus amr_1blk_guardcell is called below.
!!   This routine calls a user provided routine called restrict_fun
!!   which defines the pattern of restriction which the user wishes to
!!   apply.
!!
!! AUTHORS
!!
!!   Peter MacNeice (February 1999).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc,      & 
                                       lfulltree,filling_guardcells)


!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      use mpi_morton
      Use timings
      Use paramesh_interfaces, Only : amr_1blk_copy_soln,              & 
                                      amr_1blk_guardcell_reset,        & 
                                      amr_perm_to_1blk,                & 
                                      amr_1blk_guardcell,              & 
                                      amr_restrict_unk_fun,            & 
                                      amr_restrict_nc_fun,             & 
                                      amr_restrict_fc_fun,             & 
                                      amr_restrict_ec_fun,             & 
                                      amr_restrict_work_fun,           & 
                                      amr_restrict_work_fun_recip,     & 
                                      amr_1blk_nc_cp_remote,           & 
                                      comm_int_max_to_all,             & 
                                      comm_int_min_to_all,             & 
                                      amr_block_geometry
      Use paramesh_mpi_interfaces, Only :                              & 
                                      mpi_amr_comm_setup

      Implicit None

!-----Input/Output Arguments
      Integer, Intent(in)  :: mype,iopt
      Logical, Intent(in)  :: lcc,lfc,lec,lnc,lfulltree
      Logical, Intent(in)  :: filling_guardcells

!-----Local Variables and Arrays
      Real temp(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

      Real recvn0(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,          & 
                                            kl_bnd:ku_bnd+k3d)
      Real tempn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,       & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real sendn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,       & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real,Allocatable :: tempf(:,:,:,:)
      Real,Allocatable :: sendf(:,:,:,:)
      Real,Parameter :: eps = 1.e-30

      Integer nguard0,nguard_work0,nguard1,nguard_work1
      Integer ::  maxbnd
      Integer :: iproc
      Integer :: remote_pe0,remote_block0
      integer :: remote_pe,remote_block,icoord,nprocs,ierr
      Integer :: lb,level,ich,jchild,ioff,joff,koff,nlayers
      Integer :: idest,i,j,k,ii,jj,kk,ivar,iopt0,jface,ng0
      Integer :: ia,ja,ka,ib,jb,kb,isa,isb,jsa,jsb,ksa,ksb
      Integer :: nlayersx, nlayersy, nlayersz
      Integer :: tag_offset
      Integer :: iblk, ic, jc, kc, iblock, ilays, jlays, klays
      Integer :: id, jd, kd, is, js, ks
      Integer :: ip1, jp1, kp1, ip3, jp3, kp3
      Integer :: ng1, ndel
      Integer :: i1,i2,j1,j2,k1,k2
      Integer,Save ::  anodetype(1),aempty(1)
      Integer,Save :: cnodetype,cempty
      Integer,Save :: llrefine_min,llrefine_max
      Integer,Save :: llrefine_mint,llrefine_maxt

      Logical :: l_srl_only,ldiag
      Logical :: lguard,lprolong,lflux,ledge,lrestrict
      Logical :: lfound

!-----Include Statements
      include 'mpif.h'

!-----Begin Executable Code
      nguard0 = nguard*npgs
      nguard1 = nguard - nguard0
      nguard_work0 = nguard_work*npgs
      nguard_work1 = nguard_work - nguard_work0

      maxbnd = max(1,nbndvare,nbndvar)
      Allocate(                                                        & 
       tempf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,             & 
             kl_bnd1:ku_bnd1+k3d))
      Allocate(                                                        & 
       sendf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,             & 
             kl_bnd1:ku_bnd1+k3d))

      If (.Not.diagonals) Then
         Write(*,*) 'amr_1blk_restrict:  diagonals off'
      End If

!-----For cell-corner data, during a restriction operation, the
!-----data on a block boundary shared with a neighbor at the same
!-----refinement level, needs to be acquired during the operation
!-----of amr_1blk_guardcell_srl for the parent during the call to
!-----amr_1blk_guardcell_c_to_f. This flag tells amr_1blk_guardcell_srl
!-----to get this data. Ordinarily amr_1blk_guardcell_srl does not 
!-----get this data.
      lrestrict_in_progress = .True.

      Call amr_1blk_guardcell_reset()

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!-----Make sure the gt_unk, gt_facevarx, etc copy of the solution exists
!-----if no_permanent_guardcells is .True.. This may be needed to fill
!-----guardcells if the restriction operator needs guardcell data.
      level = -1
      If (iopt == 1) Call amr_1blk_copy_soln(level)

!-----Now parents of leaf nodes get data
!-----from their children and then perform restriction on it.

!-----Cycle through parents in decreasing order of refinement

      If (filling_guardcells) Then

       llrefine_max = 0
       llrefine_min = -1

      Else

       llrefine_max = maxval(lrefine)
       llrefine_maxt = llrefine_max
       Call comm_int_max_to_all (llrefine_max,llrefine_maxt)
       llrefine_min = llrefine_max
       Do lb = 1,lnblocks
          llrefine_min = min(lrefine(lb),llrefine_min)
       End Do
       llrefine_mint = llrefine_min
       Call comm_int_min_to_all (llrefine_min,llrefine_mint)

      End If  ! End If (filling_guardcells)

      If (llrefine_max > llrefine_min .or. filling_guardcells) Then
        Do level = llrefine_max-1,llrefine_min,-1

          tag_offset = 100
          lguard    = .True.
          lprolong  = .False.
          lflux     = .False.
          ledge     = .False.
          lrestrict = .True.
          Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,             & 
                              lflux,ledge,lrestrict,.False.,           & 
                              iopt,lcc,lfc,lec,lnc,                    & 
                              tag_offset,                              & 
                              1,1,1)

          If (lnblocks > 0) Then
            Do lb = 1,lnblocks

!-----Is this a parent block of at least one leaf node?
              If ((nodetype(lb) == 2 .and. lrefine(lb) == level) .or.          & 
                 (nodetype(lb) == 2 .and. filling_guardcells)) Then

!-----If yes then cycle through its children.
                Do ich=1,nchild

                  jchild = ich

!-------Is this child a leaf block? If it is then fetch its data.
                  remote_pe     = child(2,ich,lb)
                  remote_block  = child(1,ich,lb)

!-------if (remote_block,remote_pe) is not a local block then it must have a 
!-------local copy available in the buffer space at the end of the local
!-------block list.
                  If (remote_pe .ne. mype) Then
                    lfound = .False.
                    iblk = ladd_strt(remote_pe)
                    Do While(.not.lfound .and. iblk <= ladd_end(remote_pe))
                      If(remote_block == laddress(1,iblk).and.                   & 
                           remote_pe  == laddress(2,iblk) ) then
                        remote_block = iblk
                        remote_pe    = mype
                        lfound = .True.
                      Else
                        iblk = iblk+1
                      End If
                    End Do 
                  End If

                  cnodetype = nodetype(remote_block)
                  cempty=1
                  cempty = empty(remote_block)

                  If (cnodetype <= 2 .and. cempty == 0 ) Then

!--------compute the offset in the parent block appropriate for this child
                    ioff = mod(jchild-1,2)*nxb/2
                    joff = mod((jchild-1)/2,2)*nyb/2
                    koff = mod((jchild-1)/4,2)*nzb/2
!--------Get the child blocks data and fill its guardcells, putting the result
!--------into the current working block
                    If (iopt == 1) Then
                        nlayers = nguard
                    Else if(iopt >= 2) Then
                        nlayers = nguard_work
                    End If
!--------Put child blocks data into the data_1blk.fh datastructures, with the
!--------appropriate guardcell padding. Note, for even grid sizes the guardcells
!--------do not need to be filled.
                    idest = 1
                    Call amr_perm_to_1blk(lcc,lfc,lec,lnc,                       & 
                                   remote_block,remote_pe,             & 
                                   iopt,idest)

                    If (curvilinear) Then
!--------compute geometry variables for the child block (remote_block,remote_pe)
                      Call amr_block_geometry(remote_block,remote_pe)

                      If (curvilinear_conserve) Then

                        interp_mask_unk_res(:) = 1
                        interp_mask_work_res(:) = 1
                        interp_mask_facex_res(:) = 1
                        interp_mask_facey_res(:) = 1
                        interp_mask_facez_res(:) = 1
                        interp_mask_ec_res(:) = 1
                        interp_mask_nc_res(:) = 1
 
                        If (iopt == 1) then
 
                          i1 = 1 + nguard
                          i2 = nxb + nguard
                          j1 = 1 + k2d*nguard
                          j2 = nyb + k2d*nguard
                          k1 = 1 + k3d*nguard
                          k2 = nzb + k3d*nguard
 
!----------Compute volume weighted cell center data for conservative restriction
                          Do ivar = 1,nvar
                            If (int_gcell_on_cc(ivar))                                & 
                                unk1(ivar,i1:i2,j1:j2,k1:k2,1) =                  & 
                                   unk1(ivar,i1:i2,j1:j2,k1:k2,1)                 & 
                                  *cell_vol(i1:i2,j1:j2,k1:k2)
                          End Do
!----------Compute area weighted cell face-center data for conservative 
!----------restriction
                          Do ivar = 1,nfacevar
                           If (int_gcell_on_fc(1,ivar))                              & 
                              facevarx1(ivar,i1:i2+1,j1:j2,k1:k2,1) =               & 
                              facevarx1(ivar,i1:i2+1,j1:j2,k1:k2,1)              & 
                               *cell_area1(i1:i2+1,j1:j2,k1:k2)
                           If (int_gcell_on_fc(2,ivar))                              & 
                              facevary1(ivar,i1:i2,j1:j2+k2d,k1:k2,1) =             & 
                              facevary1(ivar,i1:i2,j1:j2+k2d,k1:k2,1)            & 
                               *cell_area2(i1:i2,j1:j2+k2d,k1:k2)
                           If (int_gcell_on_fc(2,ivar))                              & 
                              facevarz1(ivar,i1:i2,j1:j2,k1:k2+k3d,1) =             & 
                              facevarz1(ivar,i1:i2,j1:j2,k1:k2+k3d,1)            & 
                               *cell_area3(i1:i2,j1:j2,k1:k2+k3d)
                          End Do
!----------Compute distance weighted cell edge-center data for conservative 
!----------restriction
                          Do ivar = 1,nvaredge
                           If (int_gcell_on_ec(1,ivar))                              & 
                             unk_e_x1(ivar,i1:i2,j1:j2+k2d,k1:k2+k3d,1) =        & 
                              unk_e_x1(ivar,i1:i2,j1:j2+k2d,k1:k2+k3d,1)         & 
                               *cell_leng1(i1:i2,j1:j2+k2d,k1:k2+k3d)
                           If (int_gcell_on_ec(2,ivar))                              & 
                             unk_e_y1(ivar,i1:i2+1,j1:j2,k1:k2+k3d,1) =          & 
                              unk_e_y1(ivar,i1:i2+1,j1:j2,k1:k2+k3d,1)           & 
                               *cell_leng2(i1:i2+1,j1:j2,k1:k2+k3d)
                           If (int_gcell_on_ec(3,ivar))                              & 
                             unk_e_z1(ivar,i1:i2+1,j1:j2+k2d,k1:k2,1) =          & 
                              unk_e_z1(ivar,i1:i2+1,j1:j2+k2d,k1:k2,1)           & 
                               *cell_leng3(i1:i2+1,j1:j2+k2d,k1:k2)
                          End Do

                        Else

!----------Compute volume weighted cell center data for conservative restriction
!----------of work1.
                          ndel = nguard_work - nguard
                          Do k=kl_bnd1+nguard*k3d,ku_bnd1-nguard*k3d
                            Do j=jl_bnd1+nguard*k2d,ju_bnd1-nguard*k2d
                              Do i=il_bnd1+nguard    ,iu_bnd1-nguard
                                work1(i+ndel,j+ndel*k2d,k+ndel*k3d,1) =                   & 
                                   work1(i+ndel,j+ndel*k2d,k+ndel*k3d,1)                & 
                                              *cell_vol(i,j,k)
                              End Do
                            End Do
                          End Do
 
                        End If  ! End If (iopt == 1)
                      End If  ! End If (curvilinear_conserve)

!--------Now reset geometry factors to appropriate values for the 
!--------current block lb
                      Call amr_block_geometry(lb,mype)

                    End If  ! End If (curvilinear)

                    If (iopt == 1) Then

!----------Compute restricted cell-centered data from the data in the buffer
                      If (lcc) Then

                        Call amr_restrict_unk_fun(unk1(:,:,:,:,1),temp)
                        kc = koff + nguard0*k3d
                        jc = joff + nguard0*k2d
                        ic = ioff + nguard0
                        Do k=1+nguard*k3d,nzb+nguard*k3d,2
                          kk = (k-nguard*k3d)/2+1
                          kk = kk + kc
                          Do j=1+nguard*k2d,nyb+nguard*k2d,2
                            jj = (j-nguard*k2d)/2+1
                            jj = jj + jc
                            Do i=1+nguard,nxb+nguard,2
                              ii = (i-nguard)/2+1
                              ii = ii + ic
                              Do ivar=1,nvar
                                If (int_gcell_on_cc(ivar)) Then
                                  If (curvilinear_conserve) Then
                                    unk(ivar,ii,jj,kk,lb) =                         & 
                                    temp(ivar,i,j,k)                                & 
                                    / cell_vol(ii+nguard1,jj+nguard1*k2d,           & 
                                         kk+nguard1*k3d)
                                  Else
                                    unk(ivar,ii,jj,kk,lb) =                           & 
                                    temp(ivar,i,j,k)
                                  End If
                                End If
                              End Do
                            End Do
                          End Do
                        End Do

                      End If  ! End If (lcc)

!----------Compute restricted cell corner data from the data in the buffer
                      If (lnc) Then

                        Call amr_restrict_nc_fun( unk_n1(:,:,:,:,1),                & 
                                     tempn )
 
                        Do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
                          kk = (k-nguard*k3d)/2+1+nguard*k3d
                          Do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
                            jj = (j-nguard*k2d)/2+1+nguard*k2d
                            Do i=1+nguard,nxb+nguard+1,2
                              ii = (i-nguard)/2+1+nguard
                              Do ivar=1,nvarcorn
                                If (int_gcell_on_nc(ivar)) Then
                                  sendn(ivar,ii,jj,kk) = tempn(ivar,i,j,k)
                                End If
                              End Do
                            End Do
                          End Do
                        End Do 

!----------update the parent block
                        Do k=1,nzb+(-nzb/2+1)*k3d
                          Do j=1,nyb+(-nyb/2+1)*k2d
                            do i=1,nxb-nxb/2+1
                              Do ivar=1,nvarcorn
                                If (int_gcell_on_nc(ivar)) Then
                                  unk_n(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,      & 
                                      k+nguard0*k3d+koff,lb)=               & 
                                  sendn(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                End If
                              End Do
                            End Do
                          End Do
                        End Do
 
                      End If  ! End If (lnc)
                    End If  ! End If (iopt == 1)
 
!--------Compute restricted cell-face-centered data from the data in the buffer
                    If (lfc) Then
!--------Compute restricted data from the data in the buffer
                      Call amr_restrict_fc_fun(facevarx1(:,:,:,:,1),tempf,1)

                      Do k=1+nguard*k3d,nzb+nguard*k3d,2
                        kk = (k-nguard*k3d)/2+1+nguard*k3d
                        Do j=1+nguard*k2d,nyb+nguard*k2d,2
                          jj = (j-nguard*k2d)/2+1+nguard*k2d
                          Do i=1+nguard,nxb+nguard+1,2
                            ii = (i-nguard)/2+1+nguard
                            Do ivar=1,nfacevar
                               If(int_gcell_on_fc(1,ivar)) Then
                                sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
                              End If
                            End Do
                          End Do
                        End Do
                      End Do

!--------update the parent block
                      Do k=1,nzb+(-nzb/2)*k3d
                        Do j=1,nyb+(-nyb/2)*k2d
                          Do i=1,nxb-nxb/2+1
                            Do ivar=1,nfacevar
                              If (int_gcell_on_fc(1,ivar)) Then
                                If (curvilinear_conserve) Then
                                  facevarx(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,  & 
                                  k+nguard0*k3d+koff,lb)=             & 
                                  sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)    & 
                                  / cell_area1(i+nguard+ioff,j+nguard*k2d+joff,     & 
                                     k+nguard*k3d+koff)
                                Else
                                  facevarx(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,  & 
                                  k+nguard0*k3d+koff,lb)=             & 
                                  sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                 End If
                               End If
                             End Do
                           End Do
                         End Do
                       End Do

!--------y face next
                       If (ndim >= 2) then

!--------Compute restricted data from the data in the buffer

                         Call amr_restrict_fc_fun(facevary1(:,:,:,:,1),tempf,2)

                         Do k=1+nguard*k3d,nzb+nguard*k3d,2
                           kk = (k-nguard*k3d)/2+1+nguard*k3d
                           Do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
                             jj = (j-nguard*k2d)/2+1+nguard*k2d
                             Do i=1+nguard,nxb+nguard,2
                               ii = (i-nguard)/2+1+nguard
                               Do ivar=1,nfacevar
                                 If (int_gcell_on_fc(2,ivar)) Then
                                   sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
                                 End If
                               End Do
                             End Do
                           End Do
                         End Do

!--------update the parent block
                         Do k=1,nzb+(-nzb/2)*k3d
                           Do j=1,nyb+(-nyb/2+1)*k2d
                             Do i=1,nxb-nxb/2
                               Do ivar=1,nfacevar
                                 If (int_gcell_on_fc(2,ivar)) Then
                                   If (curvilinear_conserve) Then
                                     facevary(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,  & 
                                     k+nguard0*k3d+koff,lb)=                 & 
                                      sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)   & 
                                      / max(cell_area2(i+nguard+ioff,j+nguard*k2d+joff, & 
                                     k+nguard*k3d+koff),eps)
                                    Else
                                      facevary(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,    & 
                                       k+nguard0*k3d+koff,lb)=               & 
                                      sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                   End If
                                 End If
                               End Do
                             End Do
                           End Do
                         End Do

                       End If  ! End If (ndim >= 2)

!--------z face last
                       If (ndim == 3) Then

!--------Compute restricted data from the data in the buffer

                       Call amr_restrict_fc_fun(facevarz1(:,:,:,:,1),tempf,3)
  
                       Do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
                         kk = (k-nguard*k3d)/2+1+nguard*k3d
                         Do j=1+nguard*k2d,nyb+nguard*k2d,2
                           jj = (j-nguard*k2d)/2+1+nguard*k2d
                           Do i=1+nguard,nxb+nguard,2
                             ii = (i-nguard)/2+1+nguard
                             Do ivar=1,nfacevar
                               If (int_gcell_on_fc(3,ivar)) Then
                                 sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
                               End If
                             End Do
                           End Do
                         End Do
                       End Do

!--------update the parent block
                       Do k=1,nzb+(-nzb/2+1)*k3d
                         Do j=1,nyb+(-nyb/2)*k2d
                           Do i=1,nxb-nxb/2
                             Do ivar=1,nfacevar
                               If (int_gcell_on_fc(3,ivar)) Then
                                 If (curvilinear_conserve) Then
                                   facevarz(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,  & 
                                     k+nguard0*k3d+koff,lb)=                 & 
                                    sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)   & 
                                      / cell_area3(i+nguard+ioff,j+nguard*k2d+joff,   & 
                                      k+nguard*k3d+koff)
                                 Else
                                   facevarz(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,  & 
                                       k+nguard0*k3d+koff,lb)=               & 
                                    sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                 End If
                               End If
                             End Do
                           End Do
                         End Do
                       End Do

                     End If  ! End If (ndim == 3)

                   End If  ! End If (lfc)

                   If (ndim > 1) Then
!-----Compute restricted cell-edge-centered data from the data in the buffer
                     If (lec) Then
!------Compute restricted data from the data in the buffer
                       Call amr_restrict_ec_fun(unk_e_x1(:,:,:,:,1),tempf,1)

                       sendf = 0.
                       Do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
                         kk = (k-nguard*k3d)/2+1+nguard*k3d
                         Do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
                           jj = (j-nguard*k2d)/2+1+nguard*k2d
                           Do i=1+nguard,nxb+nguard,2
                             ii = (i-nguard)/2+1+nguard
                             Do ivar=1,nvaredge
                               If (int_gcell_on_ec(1,ivar)) Then
                                 sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
                               End If
                             End Do
                           End Do
                         End Do
                       End Do

!------update the parent block
                       Do k=1,nzb+(-nzb/2+1)*k3d
                         Do j=1,nyb+(-nyb/2+1)*k2d
                           Do i=1,nxb-nxb/2
                             Do ivar=1,nvaredge
                               If (int_gcell_on_ec(1,ivar)) Then
                                 If (curvilinear_conserve) Then
                                    unk_e_x(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,        & 
                                    k+nguard0*k3d+koff,lb)=                  & 
                                    sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)        & 
                                    / cell_leng1(i+nguard+ioff,j+nguard*k2d+joff,        & 
                                         k+nguard*k3d+koff)
                                 Else
                                   unk_e_x(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,        & 
                                          k+nguard0*k3d+koff,lb)=                  & 
                                   sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                 End If
                               End If            
                             End Do
                           End Do
                         End Do
                       End Do

!------y edge next
!------Compute restricted data from the data in the buffer

                       Call amr_restrict_ec_fun(unk_e_y1(:,:,:,:,1),tempf,2)

                       Do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
                         kk = (k-nguard*k3d)/2+1+nguard*k3d
                         Do j=1+nguard*k2d,nyb+nguard*k2d,2
                           jj = (j-nguard*k2d)/2+1+nguard*k2d
                           Do i=1+nguard,nxb+nguard+1,2
                             ii = (i-nguard)/2+1+nguard
                             Do ivar=1,nvaredge
                               If (int_gcell_on_ec(2,ivar)) Then
                                 sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
                               End If
                             End Do
                           End Do
                         End Do
                       End Do

!------update the parent block
                       Do k=1,nzb+(-nzb/2+1)*k3d
                         Do j=1,nyb+(-nyb/2)*k2d
                           Do i=1,nxb-nxb/2+1
                             Do ivar=1,nvaredge
                               If (int_gcell_on_ec(2,ivar)) Then
                                 If (curvilinear_conserve) then
                                   unk_e_y(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,        & 
                                    k+nguard0*k3d+koff,lb)=                  & 
                                    sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)        & 
                                    / cell_leng2(i+nguard+ioff,j+nguard*k2d+joff,        & 
                                           k+nguard*k3d+koff)
                                 Else
                                   unk_e_y(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,         & 
                                         k+nguard0*k3d+koff,lb)=                   & 
                                   sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                 End If
                               End If
                             End Do
                           End Do
                         End Do
                       End Do

                       If (ndim == 3) Then
!------z edge last
!------Compute restricted data from the data in the buffer

                         Call amr_restrict_ec_fun(unk_e_z1(:,:,:,:,1),tempf,3)
 
                         Do k=1+nguard*k3d,nzb+nguard*k3d,2
                           kk = (k-nguard*k3d)/2+1+nguard*k3d
                           Do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
                             jj = (j-nguard*k2d)/2+1+nguard*k2d
                             Do i=1+nguard,nxb+nguard+1,2
                               ii = (i-nguard)/2+1+nguard
                               Do ivar=1,nvaredge
                                 If(int_gcell_on_ec(3,ivar)) Then
                                   sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
                                 End If
                               End Do
                             End Do
                           End Do
                         End Do

!------update the parent block
                         Do k=1,nzb+(-nzb/2)*k3d
                           Do j=1,nyb+(-nyb/2+1)*k2d
                             Do i=1,nxb-nxb/2+1
                               Do ivar=1,nvaredge
                                 If (int_gcell_on_ec(3,ivar)) Then
                                   If (curvilinear_conserve) Then
                                      unk_e_z(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,        & 
                                            k+nguard0*k3d+koff,lb)=                  & 
                                      sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)        & 
                                      / cell_leng3(i+nguard+ioff,j+nguard*k2d+joff,        & 
                                           k+nguard*k3d+koff)
                                   Else
                                     unk_e_z(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,        & 
                                                k+nguard0*k3d+koff,lb)=                  & 
                                     sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                   End If
                                 End If
                               End Do
                             End Do
                           End Do
                         End Do

                       End If ! End If (ndim == 3)

                     End If ! End If (lec)
                   End If ! End If (ndim > 1)

                   If (iopt >= 2) Then
                     iopt0 = iopt-1

!---------Compute restricted cell-centered workspace data from the data 
!---------in the buffer
                     If (mod(nxb,2) == 0) Then
                       Call amr_restrict_work_fun(work1(:,:,:,1),tempw1,iopt)
                     Else
                       Call amr_restrict_work_fun_recip(work1(:,:,:,1),tempw1)
                     End If

                     kc = koff + nguard_work0*k3d
                     jc = joff + nguard_work0*k2d
                     ic = ioff + nguard_work0
                     Do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
                       kk = (k-nguard_work*k3d)/2+1
                       kk = kk + kc
                       Do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
                         jj = (j-nguard_work*k2d)/2+1
                         jj = jj + jc
                         Do i=1+nguard_work,nxb+nguard_work,2
                           ii = (i-nguard_work)/2+1
                           ii = ii + ic
                           If (curvilinear_conserve) Then
                             work(ii,jj,kk,lb,iopt0) =                           & 
                              tempw1(i,j,k)                                  & 
                             / cell_vol(ii+nguard_work1,jj+nguard_work1*k2d,     & 
                                 kk+nguard_work1*k3d)
                           Else
                             work(ii,jj,kk,lb,iopt0) =                           & 
                              tempw1(i,j,k)
                           End If
                         End Do
                       End Do
                     End Do

                   End If  ! End If (iopt >= 2)

                 End If  ! End If (cnodetype <= 2 .and. cempty == 0 )
               End Do  ! End Do ich=1,nchild
           
!-----If using odd sized grid blocks then parent copies any face bounding
!-----a leaf block
               If (iopt == 1) then
                 If (lnc) then

!-------cycle through parents neighbors
                   Do jface = 1,nfaces

!-------get this neighbors nodetype
                     remote_pe     = neigh(2,jface,lb)
                     remote_block  = neigh(1,jface,lb)

                     If (remote_block > 0) then

                       cnodetype = -1
                       If (remote_pe == mype) cnodetype = nodetype(remote_block)

!-------if (remote_block,remote_pe) is not a local block then it must have a
!-------local copy available in the buffer space at the end of the local
!-------block list.
                       remote_pe0     = remote_pe
                       remote_block0  = remote_block

                       If (remote_pe0.ne.mype) Then
                         lfound = .False.
                         iblk = ladd_strt(remote_pe0)
                         Do While (.not.lfound.and.iblk <= ladd_end(remote_pe0))
                           If(remote_block0 == laddress(1,iblk).and.                  & 
                                remote_pe0  == laddress(2,iblk) ) Then
                             remote_block0 = iblk
                             remote_pe0    = mype
                             lfound = .True.
                           Else
                             iblk = iblk+1
                           End If
                         End Do

                         If (lfound) Then
                           cnodetype = nodetype(remote_block0)
                         End If

                       End If  ! End If (remote_pe0.ne.mype)

                       If (cnodetype == 1) Then
                 
                         If(iopt == 1) Then
                            ng0 = nguard0
               	         Else if(iopt >= 2) Then
               	            ng0 = nguard_work0
               	         End If
               	         ia = 1+ng0
               	         ib = nxb+ng0+1
               	         ja = 1+ng0*k2d
               	         jb = nyb+(ng0+1)*k2d
               	         ka = 1+ng0*k3d
               	         kb = nzb+(ng0+1)*k3d
               	         isa = 1+ng0
               	         isb = nxb+ng0+1
               	         jsa = 1+ng0*k2d
               	         jsb = nyb+(ng0+1)*k2d
               	         ksa = 1+ng0*k3d
               	         ksb = nzb+(ng0+1)*k3d

               	          If (jface == 1) Then
               	            ib  = ia
               	            isa = isb
               	          Elseif (jface == 2) Then
               	            ia  = ib
               	            isb = isa
               	          Elseif (jface == 3) Then
               	            jb  = ja
               	            jsa = jsb
               	          Elseif (jface == 4) Then
               	            ja  = jb
               	            jsb = jsa
               	          Elseif (jface == 5) Then
               	            kb  = ka
               	            ksa = ksb
               	          Elseif (jface == 6) Then
               	            ka  = kb
               	            ksb = ksa
               	          End If

!--------Copy neighbor face into this block
               	          If(iopt == 1) Then
               	            If (remote_block <= lnblocks                                & 
                                	      .and.remote_pe == mype) Then
               	              unk_n(:,ia:ib,ja:jb,ka:kb,lb) = &
               	               unk_n(:,isa:isb,jsa:jsb,ksa:ksb,remote_block)
               	            Else
!------------The next section is largely borrowed from amr_1blk_guardcell_srl.
!------------It sets the index ranges for copying the common face unk_n data.
!------------Only 1 layer is required here
               	              iblock = 1
!------------Range - source indeces are initially computed as though there
!------------are no permanent guardcells.
               	              ilays = nxb
               	              jlays = nyb*k2d
               	              klays = nzb*k3d
!------------Starting indeces on destination working block
               	              id = 1 + nguard
               	              jd = 1 + nguard*k2d
               	              kd = 1 + nguard*k3d
!------------Starting indeces on source block
               	              is = 1 + nguard0
               	              js = 1 + nguard0*k2d
               	              ks = 1 + nguard0*k3d

               	              ip1 = 0
               	              jp1 = 0
               	              kp1 = 0
               	              ip3 = 0
               	              jp3 = 0
               	              kp3 = 0

               	              If (jface == 1) Then
               	                is = 1 + nxb + ng0
               	                ilays = 0
               	              Elseif (jface == 2) Then
               	                id = 1 + nxb + nguard
               	                is = 1 + ng0
               	                ilays = 0
               	              Elseif (jface == 3) Then
               	                js = 1 + (nyb+ng0)*k2d
               	                jlays = 0
               	              Elseif (jface == 4) Then
               	                jd = 1 + (nyb + nguard)*k2d
               	                js = 1 + ng0*k2d
               	                jlays = 0
               	              Elseif (jface == 5) Then
               	                ks = 1 + (nzb+ng0)*k3d
               	                klays = 0
               	              Elseif(jface == 6) Then
               	                kd = 1 + (nzb + nguard)*k3d
               	                ks = 1 + ng0*k3d
               	                klays = 0
               	              End If

!------------now we use this index info to extract remote data directly from the
!------------receive buffer. amr_1blk_nc_cp_remote write to unk_n1 so after 
!------------it exits we need to copy the result to unk_n.
               	              Call amr_1blk_nc_cp_remote(                               & 
                                	       mype,remote_pe,remote_block,iblock,             & 
               	                        id,jd,kd,is,js,ks,                              & 
                                	       ilays,jlays,klays,                              & 
               	                        ip1,jp1,kp1,ip3,jp3,kp3,0)

               	              ng1 = nguard*(1-npgs)
               	              unk_n(:,id-ng1:id+ilays-ng1,                              & 
                                	     jd-ng1*k2d:jd+(jlays-ng1)*k2d,                    & 
               	                      kd-ng1*k3d:kd+(klays-ng1)*k3d,lb) =               & 
               	              unk_n1(:,id:id+ilays,                                    & 
                                	       jd:jd+jlays*k2d,                                & 
               	                        kd:kd+klays*k3d,iblock)


               	            End If  ! End If (remote_block <= lnblocks ...

               	          End If  ! End If (iopt == 1)

               	        End If  ! End If (cnodetype == 1)
               	      End If  ! End If (remote_block > 0)

               	    End Do  ! End Do jface = 1,nfaces

                  End If  ! End If (lnc)
                End If  ! End (iopt == 1)

              End If  ! End If ((nodetype(lb) == 2 .and. lrefine(lb) == level) .or. ...
            End Do  ! End Do lb = 1,lnblocks

          End If  ! End If (lnblocks > 0)

!-----Make sure that the global copy of the newly restricted data is
!-----up to date.
          If (.not.filling_guardcells .and. iopt == 1) Then
            Call amr_1blk_copy_soln(level)
          End If

        End Do  ! End Do level = llrefine_max-1,llrefine_min,-1
      End If  ! End If (llrefine_max > llrefine_min .or. filling_guardcells)

      lrestrict_in_progress = .False.

      Deallocate(tempf)
      Deallocate(sendf)

      Return
      End Subroutine mpi_amr_1blk_restrict


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine pf_mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc,      & 
                                       lfulltree,filling_guardcells, mglevel)


!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      use mpi_morton
      Use timings
      Use paramesh_interfaces, Only : amr_1blk_copy_soln,              & 
                                      amr_1blk_guardcell_reset,        & 
                                      amr_perm_to_1blk,                & 
                                      amr_1blk_guardcell,              & 
                                      amr_restrict_unk_fun,            & 
                                      amr_restrict_nc_fun,             & 
                                      amr_restrict_fc_fun,             & 
                                      amr_restrict_ec_fun,             & 
                                      amr_restrict_work_fun,           & 
                                      amr_restrict_work_fun_recip,     & 
                                      amr_1blk_nc_cp_remote,           & 
                                      comm_int_max_to_all,             & 
                                      comm_int_min_to_all,             & 
                                      amr_block_geometry
      Use paramesh_mpi_interfaces, Only :                              & 
                                      mpi_amr_comm_setup

      Implicit None

!-----Input/Output Arguments
      Integer, Intent(in)  :: mype,iopt, mglevel
      Logical, Intent(in)  :: lcc,lfc,lec,lnc,lfulltree
      Logical, Intent(in)  :: filling_guardcells

!-----Local Variables and Arrays
      Real temp(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

      Real recvn0(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,          & 
                                            kl_bnd:ku_bnd+k3d)
      Real tempn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,       & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real sendn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,       & 
                                            kl_bnd1:ku_bnd1+k3d)
      Real,Allocatable :: tempf(:,:,:,:)
      Real,Allocatable :: sendf(:,:,:,:)
      Real,Parameter :: eps = 1.e-30

      Integer nguard0,nguard_work0,nguard1,nguard_work1
      Integer ::  maxbnd
      Integer :: iproc
      Integer :: remote_pe0,remote_block0
      integer :: remote_pe,remote_block,icoord,nprocs,ierr
      Integer :: lb,level,ich,jchild,ioff,joff,koff,nlayers
      Integer :: idest,i,j,k,ii,jj,kk,ivar,iopt0,jface,ng0
      Integer :: ia,ja,ka,ib,jb,kb,isa,isb,jsa,jsb,ksa,ksb
      Integer :: nlayersx, nlayersy, nlayersz
      Integer :: tag_offset
      Integer :: iblk, ic, jc, kc, iblock, ilays, jlays, klays
      Integer :: id, jd, kd, is, js, ks
      Integer :: ip1, jp1, kp1, ip3, jp3, kp3
      Integer :: ng1, ndel
      Integer :: i1,i2,j1,j2,k1,k2
      Integer,Save ::  anodetype(1),aempty(1)
      Integer,Save :: cnodetype,cempty
      Integer,Save :: llrefine_min,llrefine_max
      Integer,Save :: llrefine_mint,llrefine_maxt

      Logical :: l_srl_only,ldiag
      Logical :: lguard,lprolong,lflux,ledge,lrestrict
      Logical :: lfound

!-----Include Statements
      include 'mpif.h'

!-----Begin Executable Code
      nguard0 = nguard*npgs
      nguard1 = nguard - nguard0
      nguard_work0 = nguard_work*npgs
      nguard_work1 = nguard_work - nguard_work0

      maxbnd = max(1,nbndvare,nbndvar)
      Allocate(                                                        & 
       tempf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,             & 
             kl_bnd1:ku_bnd1+k3d))
      Allocate(                                                        & 
       sendf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,             & 
             kl_bnd1:ku_bnd1+k3d))

      If (.Not.diagonals) Then
         Write(*,*) 'amr_1blk_restrict:  diagonals off'
      End If

!-----For cell-corner data, during a restriction operation, the
!-----data on a block boundary shared with a neighbor at the same
!-----refinement level, needs to be acquired during the operation
!-----of amr_1blk_guardcell_srl for the parent during the call to
!-----amr_1blk_guardcell_c_to_f. This flag tells amr_1blk_guardcell_srl
!-----to get this data. Ordinarily amr_1blk_guardcell_srl does not 
!-----get this data.
      lrestrict_in_progress = .True.

      Call amr_1blk_guardcell_reset()

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!-----Make sure the gt_unk, gt_facevarx, etc copy of the solution exists
!-----if no_permanent_guardcells is .True.. This may be needed to fill
!-----guardcells if the restriction operator needs guardcell data.
      level = -1
      If (iopt == 1) Call amr_1blk_copy_soln(level)

!-----Now parents of leaf nodes get data
!-----from their children and then perform restriction on it.

!-----Cycle through parents in decreasing order of refinement

      If (filling_guardcells) Then

         llrefine_min = 1
         llrefine_max = MAX(mglevel+1,lrefine_max)

      Else

         llrefine_min = mglevel-1
         llrefine_max = mglevel

      End If  ! End If (filling_guardcells)

!      If (llrefine_max > llrefine_min .or. filling_guardcells) Then
      If (llrefine_max > llrefine_min) Then
!      If (mglevel > llrefine_min) Then
        Do level = llrefine_max-1,llrefine_min,-1
!        level = mglevel-1

          tag_offset = 100
          lguard    = .True.
          lprolong  = .False.
          lflux     = .False.
          ledge     = .False.
          lrestrict = .True.
          Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,             & 
                              lflux,ledge,lrestrict,.False.,           & 
                              iopt,lcc,lfc,lec,lnc,                    & 
                              tag_offset,                              & 
                              1,1,1)

          If (lnblocks > 0) Then
            Do lb = 1,lnblocks
!print *,'levels',mglevel-1, mglevel, llrefine_min
!print *,block_starts(mglevel-1), block_starts(mglevel)-1
!            Do lb = block_starts(llrefine_min),block_starts(llrefine_max)-1

!-----Is this a parent block of at least one leaf node?
              If ((nodetype(lb) == 2 .and. lrefine(lb) == level) .or.          & 
                 (nodetype(lb) == 2 .and. filling_guardcells)) Then


!-----If yes then cycle through its children.
                Do ich=1,nchild

                  jchild = ich

!-------Is this child a leaf block? If it is then fetch its data.
                  remote_pe     = child(2,ich,lb)
                  remote_block  = child(1,ich,lb)

!-------if (remote_block,remote_pe) is not a local block then it must have a 
!-------local copy available in the buffer space at the end of the local
!-------block list.
                  If (remote_pe .ne. mype) Then
                    lfound = .False.
                    iblk = ladd_strt(remote_pe)
                    Do While(.not.lfound .and. iblk <= ladd_end(remote_pe))
                      If(remote_block == laddress(1,iblk).and.                   & 
                           remote_pe  == laddress(2,iblk) ) then
                        remote_block = iblk
                        remote_pe    = mype
                        lfound = .True.
                      Else
                        iblk = iblk+1
                      End If
                    End Do 
                  End If

                  cnodetype = nodetype(remote_block)
                  cempty=1
                  cempty = empty(remote_block)

                  If (cnodetype <= 2 .and. cempty == 0 ) Then

!--------compute the offset in the parent block appropriate for this child
                    ioff = mod(jchild-1,2)*nxb/2
                    joff = mod((jchild-1)/2,2)*nyb/2
                    koff = mod((jchild-1)/4,2)*nzb/2
!--------Get the child blocks data and fill its guardcells, putting the result
!--------into the current working block
                    If (iopt == 1) Then
                        nlayers = nguard
                    Else if(iopt >= 2) Then
                        nlayers = nguard_work
                    End If
!--------Put child blocks data into the data_1blk.fh datastructures, with the
!--------appropriate guardcell padding. Note, for even grid sizes the guardcells
!--------do not need to be filled.
                    idest = 1
!                    Call amr_perm_to_1blk(lcc,lfc,lec,lnc,                       & 
                    Call pf_perm_to_1blk(lcc,lfc,lec,lnc,                       & 
                                   remote_block,remote_pe,             & 
                                   iopt,idest)

                    If (iopt == 1) Then

!----------Compute restricted cell-centered data from the data in the buffer
                      If (lcc) Then

                        Call amr_restrict_unk_fun(unk1(:,:,:,:,1),temp)
                        kc = koff + nguard0*k3d
                        jc = joff + nguard0*k2d
                        ic = ioff + nguard0
                        Do k=1+nguard*k3d,nzb+nguard*k3d,2
                          kk = (k-nguard*k3d)/2+1
                          kk = kk + kc
                          Do j=1+nguard*k2d,nyb+nguard*k2d,2
                            jj = (j-nguard*k2d)/2+1
                            jj = jj + jc
                            Do i=1+nguard,nxb+nguard,2
                              ii = (i-nguard)/2+1
                              ii = ii + ic
                              Do ivar=1,nvar
                                If (int_gcell_on_cc(ivar)) Then
                                    unk(ivar,ii,jj,kk,lb) = temp(ivar,i,j,k)
                                End If
                              End Do
                            End Do
                          End Do
                        End Do

                      End If  ! End If (lcc)

                   End If  ! End If (iopt == 1)
 
                   If (iopt >= 2) Then
                     iopt0 = iopt-1


!---------Compute restricted cell-centered workspace data from the data 
!---------in the buffer
                     Call amr_restrict_work_fun(work1(:,:,:,1),tempw1,iopt)

                     kc = koff + nguard_work0*k3d
                     jc = joff + nguard_work0*k2d
                     ic = ioff + nguard_work0
                     Do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
                       kk = (k-nguard_work*k3d)/2+1
                       kk = kk + kc
                       Do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
                         jj = (j-nguard_work*k2d)/2+1
                         jj = jj + jc
                         Do i=1+nguard_work,nxb+nguard_work,2
                           ii = (i-nguard_work)/2+1
                           ii = ii + ic
                           work(ii,jj,kk,lb,iopt0) = tempw1(i,j,k)
                         End Do
                       End Do
                     End Do


                   End If  ! End If (iopt >= 2)

                 End If  ! End If (cnodetype <= 2 .and. cempty == 0 )
               End Do  ! End Do ich=1,nchild
           
!-----If using odd sized grid blocks then parent copies any face bounding
!-----a leaf block

              End If  ! End If ((nodetype(lb) == 2 .and. lrefine(lb) == level) .or. ...
            End Do  ! End Do lb = 1,lnblocks

          End If  ! End If (lnblocks > 0)

!-----Make sure that the global copy of the newly restricted data is
!-----up to date.
          If (.not.filling_guardcells .and. iopt == 1) Then
            Call amr_1blk_copy_soln(level)
          End If

        End Do  ! End Do level = llrefine_max-1,llrefine_min,-1
      End If  ! End If (llrefine_max > llrefine_min .or. filling_guardcells)

      lrestrict_in_progress = .False.

      Deallocate(tempf)
      Deallocate(sendf)


      Return
      End Subroutine pf_mpi_amr_1blk_restrict
