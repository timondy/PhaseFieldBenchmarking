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
!!   amr_1blk_guardcell
!!
!! SYNOPSIS
!!
!!   call amr_1blk_guarcell(mype, iopt, nlayers, lb, pe, 
!!                          lcc, lfc, lec, lnc, l_srl_only,
!!                          icoord, ldiag)
!!   call amr_1blk_guarcell(mype, iopt, nlayers, lb, pe, 
!!                          lcc, lfc, lec, lnc, l_srl_only,
!!                          icoord, ldiag, 
!!                          nlayersx, nlayersy, nlayersz)
!!
!!   call amr_1blk_guarcell(integer, integer, integer, integer, integer, 
!!                          logical, logical, logical, logical, logical,
!!                          integer, logical, 
!!                          optional integer, optional integer, optional integer)
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
!!   integer, intent(in) :: nlayers        
!!        The number of guard cell layers at each boundary 
!!         (Note: this argument is deprecated not longer functions !).
!!
!!   integer, intent(in) :: lb             
!!          The block number on processor pe selected for guardcell filling.
!!
!!   integer, intent(in) :: pe             
!!        The processor storing the block selected for guarcell filling.
!!
!!   logical, intent(in) :: lcc            
!!        A logical switch controlling whether unk or work data is filled.
!!
!!   logical, intent(in) :: lfc            
!!        A logical switch controlling whether facevarx(y)(z) data is filled.
!!
!!   logical, intent(in) :: lec            
!!        A logical switch controlling whether unk_e_x(y)(z) data is filled.
!!
!!   logical, intent(in) :: lnc            
!!        A logical switch controlling whether unk_n data is filled.
!!
!!   logical, intent(in) :: l_srl_only     
!!        A logical switch which, if true, switches off the filling from 
!!        coarse neighbors. This is used during restriction when odd  
!!        block sizes are in use.
!!
!!   integer, intent(in) :: icoord         
!!        An integer switch used to select which faces of the block are to be 
!!        considered. If icoord=0 all faces are considered. If icoord=1 only 
!!        faces perp. to the x-axis are considered, if icoord=2 only faces
!!        perp. to the y-axis are considered, and if icoord=3 only faces perp. 
!!        to the z-axis are considered.
!!
!!   logical, intent(in) :: ldiag          
!!        A logical switch which controls whether guardcells corresponding to 
!!        neighbor blocks diagonally opposite block edges and corners are filled.
!!
!!   optional, integer, intent(in) :: nlayersx, nlayersy, nlayersz
!!        Optional arguments to specify the number of guardcells to fill in each 
!!        coordninate direction.  If not specified, then the default is to fill all
!!        available guardcells.
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
!!   timings
!!   workspace
!!   mpi_morton
!!   paramesh_mpi_interfaces
!!   paramesh_interfaces
!!
!! CALLS
!!
!!     mpi_amr_local_surr_blks_lkup
!!     mpi_amr_1blk_guardcell_c_to_f
!!     mpi_amr_get_remote_block
!!     amr_perm_to_1blk
!!     amr_1blk_guardcell_srl
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guardcells are filled with data for 
!!   the block secified in the call to this routine.
!!
!! DESCRIPTION
!!
!!   This routine manages the transfer of guard cell data for a
!!   specific single block. It uses the morton numbering scheme for the
!!   determination of the neighbor relations and differs in that from
!!   the older routine amr_1blk_guardcell.  The user can choose how many layers 
!!   of guardcells are filled on each face of a block by specifying the optional 
!!   arguments 'nlayersx', 'nlayersy', and 'nlayersz'.
!!
!!   The sequence of operations required to fill the guardcells of an
!!   individual block is :
!!
!!   Step 1:
!!   Construct a list of blocks surrounding block lb
!!
!!   Step 2:
!!   Check for coarse neighbors
!!
!!   Step 3:
!!   Put leaf block data into the data_1blk.fh datastructures, with the
!!   appropriate guardcell padding.
!!
!!   Step 4:
!!   Put data from leaf blocks parent into the data_1blk.fh datastructures, 
!!   with the appropriate guardcell padding. Check to see if this data is 
!!   currently cached.
!!
!!   Step 5:
!!   Construct a list of blocks surrounding block lb's parent
!!
!!   Step 6:
!!   Do guardcell filling for lb's parent from any surrounding blocks at 
!!   the same refinement level as this parent.
!!
!!   Step 7:
!!   Do guardcell filling from coarse neigbors into the current block
!!
!!   Step 8:
!!   Do guardcell filling from any surrounding blocks at the same refinement
!!   level as the leaf block.
!!
!!   Step 9:
!!   Apply boundary conditions.
!!
!! NOTES
!!
!!  This routine must be used with care !
!!  This routine was written to be used in a code as illustrated
!!  in the following snippet of pseudo-code
!!
!!              .
!!              .
!!              .
!!        synchronization pt
!!
!!        loop over grid blocks
!!          (copy current block into working block and fill its guardcells)
!!          (perform some set of operations on block)
!!          (store result from working block)
!!        end loop
!!
!!        synchronization pt
!!              .
!!              .
!!              .
!!
!!  Caveat 1:
!!   If you are using this routine, you must remember to save the solution
!!   at the first synchronization point (ie call amr_1blk_copy_soln), so 
!!   that each block uses the same time synchronized solution during its 
!!   guardcell filling.
!!
!!  Caveat 2:
!!   It is implicitly assumed that the parent blocks on all leaf nodes
!!   have valid data. (This is necessary to ensure that a general restriction
!!   operator can be supported in the neighborhood of a jump in refinement.)
!!   If ADVANCE_ALL_LEVELS is defined, then this will generally be true. 
!!   However, if the solution is being time-advanced on leaf blocks only 
!!   this may not be true. In this case you should call amr_restrict  
!!   before the first synchronization pt in the example above.
!!   If you are using blocks with an even number of grid points and the
!!   default restriction operators this is not necessary.
!!
!! AUTHORS
!!
!!   Michael Gehmeyr & Peter MacNeice (November 1999) with modifications by 
!!   Kevin Olson for layered guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_guardcell(                                   & 
                                    mype,iopt,nlayers,lb,pe,           & 
                                    lcc,lfc,lec,lnc,                   & 
                                    l_srl_only,icoord,ldiag,           & 
                                    nlayersx,nlayersy,nlayersz)
!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use workspace
      Use mpi_morton
      Use constants
      Use paramesh_mpi_interfaces, Only :                              & 
                                      mpi_amr_local_surr_blks_lkup,    & 
                                      mpi_amr_1blk_guardcell_c_to_f,   & 
                                      mpi_amr_get_remote_block
      Use paramesh_interfaces, Only : amr_perm_to_1blk,                & 
                                      amr_1blk_guardcell_srl

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguements.
      Logical, intent(in) :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      Integer, intent(in) :: mype,lb,pe,iopt,nlayers,icoord
      Integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!-----Local arrays and variables
      Integer :: nprocs
      Integer :: parent_lb,parent_pe
      Integer :: i,j,k, ll, idest, iblock
      Integer :: ij,ijk
      Integer :: surrblks(3,3,3,3), tsurrblks(3,3,3,3)
      Integer :: psurrblks(3,3,3,3)
      Integer :: pcache_pe,pcache_blk
      Integer :: ierrorcode,ierr
      Integer :: ipolar(2)
      Integer :: ippolar(2),pbnd_box(2,3)
      Integer :: icoord_loc
      Integer :: nlayers0x, nlayers0y, nlayers0z

      Real    :: eps,accuracy
      Double Precision :: time1
      Double Precision :: time2

      Logical :: lcoarse,l_parent
      Logical :: loc_lcc,loc_lfc,loc_lec,loc_lnc
      Logical :: lfound
      Logical :: ldiag_loc
      Logical, Save :: first_cc = .True.
      Logical, Save :: first_nc = .True.
      Logical, Save :: first_ec = .True.
      Logical, Save :: first_fc = .True.

!------Begin Executable code.

! print*, lb,'In amr_1blk_guardcell lb=',lb

! CEG some declarations
      ipolar(1) = -100
      ipolar(2) = -100
      ippolar(1) = -100
      ippolar(2) = -100

      accuracy = 10./10.**precision(accuracy)
      eps = accuracy

      If (timing_mpi) Then
         time1 = mpi_wtime()
         time2 = mpi_wtime()
      End If

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      If (iopt == 1) Then
         If (.Not.present(nlayersx)) Then
            nlayers0x = nguard
         Else
            nlayers0x = nlayersx
         End If
         If (.Not.present(nlayersy)) Then
            nlayers0y = nguard
         Else
            nlayers0y = nlayersy
         End If
         If (.Not.present(nlayersz)) Then
            nlayers0z = nguard
         Else
            nlayers0z = nlayersz
         End If
      Else
         If (.Not.present(nlayersx)) Then
            nlayers0x = nguard_work
         Else
            nlayers0x = nlayersx
         End If
         If (.Not.present(nlayersy)) Then
            nlayers0y = nguard_work
         Else
            nlayers0y = nlayersy
         End If
         If (.Not.present(nlayersz)) Then
            nlayers0z = nguard_work
         Else
            nlayers0z = nlayersz
         End If
      End If  ! End If (iopt == 1)

      If (nxb/nguard < 2) nlayers0x = min(nlayers0x+1,   nguard)
      If (nyb/nguard < 2) nlayers0y = min(nlayers0y+k2d, nguard)
      If (nzb/nguard < 2) nlayers0z = min(nlayers0z+k3d, nguard)
      
      If (iopt == 1 .And.                                              & 
          All(interp_mask_unk(:) < 20)   .And.                         & 
          All(interp_mask_facex(:) < 20) .And.                         & 
          All(interp_mask_facey(:) < 20) .And.                         & 
          All(interp_mask_facez(:) < 20) .And.                         & 
          All(interp_mask_ec(:) < 20)    .And.                         & 
          All(interp_mask_nc(:) < 20) ) Then
         If (Any(nguard < interp_mask_unk(:)-2)   .Or.                 & 
             Any(nguard < interp_mask_facex(:)-2) .Or.                 & 
             Any(nguard < interp_mask_facey(:)-2) .Or.                 & 
             Any(nguard < interp_mask_facez(:)-2) .Or.                 & 
             Any(nguard < interp_mask_ec(:)-2)    .Or.                 & 
             Any(nguard < interp_mask_nc(:)-2)) Then
            Print *,' PARAMESH ERROR: nguard is too small for the' 
            Print *,' chosen interpolation order !!! '
            Print *,' nguard = ',nguard
            Print *,' maxval(interp_mask_unk) = ',                     & 
                      maxval(interp_mask_unk(:))
            Print *,' maxval(interp_mask_facex) = ',                   & 
                      maxval(interp_mask_facex(:))
            Print *,' maxval(interp_mask_facey) = ',                   & 
                      maxval(interp_mask_facey(:))
            Print *,' maxval(interp_mask_facez) = ',                   & 
                      maxval(interp_mask_facez(:))
            Print *,' maxval(interp_mask_ec) = ',                      & 
                      maxval(interp_mask_ec(:))
            Print *,' maxval(interp_mask_nc) = ',                      & 
                      maxval(interp_mask_nc(:))
            Call amr_abort()
         End If
         If ( (Any(nlayers0x < interp_mask_unk(:)-2) .And.             & 
                   nvar > 0) .Or.                                      & 
              (Any(nlayers0x < interp_mask_facex(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0x < interp_mask_facey(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0x < interp_mask_facez(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0x < interp_mask_ec(:)-2) .And.              & 
                   nvaredge > 0) .Or.                                  & 
              (Any(nlayers0x < interp_mask_nc(:)-2) .And.              & 
                   nvarcorn > 0) ) Then
            Print *,' PARAMESH ERROR: nlayersx is too small for the'
            Print *,' chosen interpolation order !!! '
            Print *,' nlayersx = ',nlayers0x
            Print *,' maxval(interp_mask_unk) = ',                     & 
                      maxval(interp_mask_unk(:))
            Print *,' maxval(interp_mask_facex) = ',                   & 
                      maxval(interp_mask_facex(:))
            Print *,' maxval(interp_mask_facey) = ',                   & 
                      maxval(interp_mask_facey(:))
            Print *,' maxval(interp_mask_facez) = ',                   & 
                      maxval(interp_mask_facez(:))
            Print *,' maxval(interp_mask_ec) = ',                      & 
                      maxval(interp_mask_ec(:))
            Print *,' maxval(interp_mask_nc) = ',                      & 
                      maxval(interp_mask_nc(:))
            Call amr_abort()
         End If
         If ( (Any(nlayers0y < interp_mask_unk(:)-2) .And.             & 
                   nvar > 0) .Or.                                      & 
              (Any(nlayers0y < interp_mask_facex(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0y < interp_mask_facey(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0y < interp_mask_facez(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0y < interp_mask_ec(:)-2) .And.              & 
                   nvaredge > 0) .Or.                                  & 
              (Any(nlayers0y < interp_mask_nc(:)-2) .And.              & 
                   nvarcorn > 0) ) Then
            Print *,' PARAMESH ERROR: nlayersy is too small for the' 
            Print *,' chosen interpolation order !!! '
            Print *,' nlayersy = ',nlayers0y
            Print *,' maxval(interp_mask_unk) = ',                     & 
                      maxval(interp_mask_unk(:))
            Print *,' maxval(interp_mask_facex) = ',                   & 
                      maxval(interp_mask_facex(:))
            Print *,' maxval(interp_mask_facey) = ',                   & 
                      maxval(interp_mask_facey(:))
            Print *,' maxval(interp_mask_facez) = ',                   & 
                      maxval(interp_mask_facez(:))
            Print *,' maxval(interp_mask_ec) = ',                      & 
                      maxval(interp_mask_ec(:))
            Print *,' maxval(interp_mask_nc) = ',                      & 
                      maxval(interp_mask_nc(:))
            Call amr_abort()
         End If
         If ( (Any(nlayers0z < interp_mask_unk(:)-2) .And.             & 
                   nvar > 0) .Or.                                      & 
              (Any(nlayers0z < interp_mask_facex(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0z < interp_mask_facey(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0z < interp_mask_facez(:)-2) .And.           & 
                   nfacevar > 0) .Or.                                  & 
              (Any(nlayers0z < interp_mask_ec(:)-2) .And.              & 
                   nvaredge > 0) .Or.                                  & 
              (Any(nlayers0z < interp_mask_nc(:)-2) .And.              & 
                   nvarcorn > 0) ) Then
            Print *,' PARAMESH ERROR: nlayersz is too small for the'
            Print *,' chosen interpolation order !!! '
            Print *,' nlayersz = ',nlayers0z
            Print *,' maxval(interp_mask_unk) = ',                     & 
                      maxval(interp_mask_unk(:))
            Print *,' maxval(interp_mask_facex) = ',                   & 
                      maxval(interp_mask_facex(:))
            Print *,' maxval(interp_mask_facey) = ',                   & 
                      maxval(interp_mask_facey(:))
            Print *,' maxval(interp_mask_facez) = ',                   & 
                      maxval(interp_mask_facez(:))
            Print *,' maxval(interp_mask_ec) = ',                      & 
                      maxval(interp_mask_ec(:))
            Print *,' maxval(interp_mask_nc) = ',                      & 
                      maxval(interp_mask_nc(:))
            Call amr_abort()
         End If
      ElseIf (iopt >= 2 .And. All(interp_mask_work(:) < 20) ) Then
         If (Any(nguard_work < interp_mask_work(:)-2)) Then
           Print *,' PARAMESH ERROR: nguard_work is too small for the'
           Print *,' chosen interpolation order !!! '
           Print *,' nguard_work = ',nguard_work
           Print *,' maxval(interp_mask_work) = ',                     & 
                     maxval(interp_mask_work(:))
            Call amr_abort()
         End If
         If (Any(nlayers0x < interp_mask_work(:)-2) .And.              & 
                 nvar_work > 0) Then
            Print *,' PARAMESH ERROR: nlayersx is too small for the'
            Print *,' chosen interpolation order !!! '
            Print *,' nlayersx = ',nlayers0x
            Print *,' maxval(interp_mask_work) = ',                    & 
                      maxval(interp_mask_work(:))
            Call amr_abort()
         End If
         If (Any(nlayers0y < interp_mask_work(:)-2) .And.              & 
                 nvar_work > 0) Then
            Print *,' PARAMESH ERROR: nlayersy is too small for the'
            Print *,' chosen interpolation order !!! '
            Print *,' nlayersy = ',nlayers0y
            Print *,' maxval(interp_mask_work) = ',                    & 
                      maxval(interp_mask_work(:))
            Call amr_abort()
         End If
         If (Any(nlayers0z < interp_mask_work(:)-2) .And.              & 
                 nvar_work > 0) Then
            Print *,' PARAMESH ERROR: nlayers0z is too small for the'
            Print *,' chosen interpolation order !!! '
            Print *,' nlayersz = ',nlayers0z
            Print *,' maxval(interp_mask_work) = ',                    & 
                      maxval(interp_mask_work(:))
            Call amr_abort()
         End If
      End If  ! End If (iopt == 1 .And. ...)

      If (pe.Ne.mype) Then
          Write(*,*) 'Error : trying to fill guardcells for a ',       & 
                     'remote block - not supported with the mpi ',     & 
                     'version. '
          Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      If (                                                             & 
          mpi_pattern_id == 10  .Or.                                   & 
          (mpi_pattern_id == 20.And.lprolong_in_progress) .Or.         & 
          (mpi_pattern_id == 40.And.lrestrict_in_progress)             & 
                                  ) Then
      ElseIf (nprocs > 1) Then
        Write(*,*) 'Paramesh error : amr_1blk_guardcell : ',           & 
       ' wrong pattern being',                                         & 
       ' used for pre-communication for guardcell fill : Fix',         & 
       ' - insert appropriate Call to mpi_amr_comm_setup '             & 
       ,' mpi_pattern_id ',mpi_pattern_id,                             & 
       ' lprolong_in_progress ',lprolong_in_progress,                  & 
       ' lrestrict_in_progress ',lrestrict_in_progress
        Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      End If

      If (amr_error_checking) Then

      If (iopt == 1) Then   

        loc_lcc = .False.
        loc_lfc = .False.
        loc_lnc = .False.
        loc_lec = .False.
        If (nvar > 0)     loc_lcc = .True.
        If (nfacevar > 0) loc_lfc = .True.
        If (nvaredge > 0) loc_lec = .True.
        If (nvarcorn > 0) loc_lnc = .True.
        If ( (lcc.And.(.Not.loc_lcc)) .Or. (lfc.And.(.Not.loc_lfc))    & 
        .Or.(lec.And.(.Not.loc_lec)) .Or. (lnc.And.(.Not.loc_lnc))     & 
          ) Then
          If (mype == 0) Then
            Write(*,*) 'Paramesh Error: for a Call to',                & 
             '  mpi_amr_1blk_guardcell',                               & 
             ' one of more of the arguments lcc/lfc/lec/lnc are not',  & 
             ' consistent with nvar/nfacevar/nvaredge/nvarcorn.'
            Print *,' nvar = ',nvar,lcc,loc_lcc
          End If
          Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        End If

      Else

        loc_lcc = .False.
        If (nvar_work > 0)     loc_lcc = .True.
        If (lcc.And.(.Not.loc_lcc)) Then
          If (mype == 0) Then
            Write(*,*) 'Paramesh Error: for a Call to',                & 
             '  mpi_amr_1blk_guardcell',                               & 
             ' one of more of the arguments lcc/lfc/lec/lnc are not',  & 
             ' consistent with nvar/nfacevar/nvaredge/nvarcorn.'
            Print *,' nvar = ',nvar,lcc,loc_lcc
          End If
          Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        End If

      End If  ! End If (iopt == 1)

      End If  ! End If (amr_error_checking)

!-----Intialization
      surrblks  = -1         
      tsurrblks = -1  

!-----construct a list of blocks surrounding local block lb
      l_parent = .True.
     If (l_srl_only) l_parent = .False.
      Call mpi_amr_local_surr_blks_lkup(mype,lb,                       & 
                                surrblks,l_parent,psurrblks)

!-----relate surrblks with the guard block indices stored implicitly 
!-----in laddress and update tsurrblks
      tsurrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =                         &
          surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb)

      If (timing_mpi) Then
         time2 = mpi_wtime()
      End If

!------guard block indeces start at strt_buffer after lnblocks, and end at
!------last_buffer as determined in subroutine mpi_commatrix.
      Do k = 2-k3d,2+k3d     ! loop over all its surrounding blocks
       Do j = 2-k2d,2+k2d
       Do i = 1,3
         If ( surrblks(2,i,j,k).Ne.mype .And.                          & 
              surrblks(2,i,j,k) >= 0 ) Then

         lfound = .False.
         ll = ladd_strt(surrblks(2,i,j,k))
         Do While(.Not.lfound.And.ll <= ladd_end(surrblks(2,i,j,k)))
           If ( (surrblks(2,i,j,k) == laddress(2,ll)) .And.            & 
                (surrblks(1,i,j,k) == laddress(1,ll)) ) Then
!------------found the corresponding block id ll
             tsurrblks(1,i,j,k) = ll
             tsurrblks(2,i,j,k) = mype
             tsurrblks(3,i,j,k) = nodetype(ll) 
             lfound = .True.
           Else
             ll = ll+1  
           End If  ! End If ( (surrblks(2,i,j,k) == laddress(2,ll)) .And. ...
         End Do  ! End Do While(.Not.lfound.And.ll <= ladd_end(surrblks(2,i,j,k)))

         End If  ! End If ( surrblks(2,i,j,k).Ne.mype .And. ...

         If ( (tsurrblks(2,i,j,k).Ne.mype) .And.                       & 
              (tsurrblks(2,i,j,k).Ne.-1)   .And.                       & 
              (tsurrblks(1,i,j,k) > -20) ) Then
             write(*,*)tsurrblks(2,i,j,k),mype
             Write(*,*) 'ERROR in mpi_amr_1blk_guardcell : pe ',mype,  & 
               ' working on lb ',lb,' neigh ',i,j,k,                   & 
               ' cannot find surrblk ',                                & 
               surrblks(:,i,j,k),' on this proc ',                     & 
               ' laddress ',laddress(:,strt_buffer:last_buffer),       & 
               ' strt_buffer,last_buffer ',strt_buffer,last_buffer,    & 
               ' tsurrblks ',tsurrblks(1:2,i,j,k)                      & 
             ,' ladd_strt ',ladd_strt,' ladd_end ',ladd_end            & 
             ,' comm pattern id ',mpi_pattern_id
             write(*,*) 'CEG thinks you are probably out of blocks!'

             Call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
         End If

        End Do  ! End Do i = 1, 3
        End Do  ! End Do j = 2-k2d,2+k2d
       End Do  ! End Do k = 2-k3d,2+k3d

!-----update surrblks with the local guard block info
      surrblks = tsurrblks

      If (spherical_pm) Then
       ipolar = 0
       If (lsingular_line) Then
       If (abs(bnd_box(1,2,lb)) < eps) ipolar(1) = -1
       If (abs(bnd_box(2,2,lb)-pi) < eps) ipolar(2) = 1
       End If
      End If

      If (timing_mpi) Then
      timer_amr_1blk_guardcell(1) = timer_amr_1blk_guardcell(1)        & 
                                 +  mpi_wtime() - time2
      time2 = mpi_wtime()
      End If

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      If (.Not.l_srl_only) Then

!-------Are there any coarse neighbors?
        lcoarse = .False.

        If (parent(1,lb) > 0) Then
        Do k = 2-k3d,2+k3d
        Do j = 2-k2d,2+k2d
        Do i = 1,3
          If (surrblks(1,i,j,k) > -20.And.surrblks(1,i,j,k) < 0)       & 
                         lcoarse = .True.
        End Do
        End Do
        End Do
        End If 

      End If  ! End If (.Not.l_srl_only)

!-----Put leaf block lb's data into the '1blk' datastructures, 
!-----with the appropriate guardcell padding.
      idest = 1
      Call amr_perm_to_1blk(lcc,lfc,lec,lnc,lb,pe,iopt,idest)

      If (iopt == 1) Then
          pcache_pe  = pcache_pe_u
          pcache_blk = pcache_blk_u
      ElseIf (iopt >= 2) Then
          pcache_pe  = pcache_pe_w
          pcache_blk = pcache_blk_w
      End If

      If (.Not.l_srl_only) Then
        If (lcoarse) Then

!---------Put data from lb's parent into the data_1blk.fh datastructures, with the
!---------appropriate guardcell padding. Check to see if data is currently cached.
          parent_lb = parent(1,lb)
          parent_pe = parent(2,lb)

        If ( (parent_lb > 0) .And.                                     & 
            ((parent_lb.Ne.pcache_blk).Or.(parent_pe.Ne.pcache_pe) )   & 
            ) Then

!---------record id of new parent block placed in cache
          lnew_parent = .True.
          pcache_blk = parent_lb
          pcache_pe  = parent_pe

        If (lcc) Then
           If (first_cc) Then
              unk1(:,:,:,:,2) = 0.
              if (nvar_work > 0) work1(:,:,:,2) = 0.
              first_cc = .False.
           End If
        End If
        If (lnc) Then
           If (first_nc) Then
              unk_n1(:,:,:,:,2) = 0.
              first_nc = .False.
           End If
        End If
        If (lec) Then
           If (first_ec) Then
             unk_e_x1(:,:,:,:,2) = 0.
             unk_e_y1(:,:,:,:,2) = 0.
             unk_e_z1(:,:,:,:,2) = 0.
             first_ec = .False.
           End If
        End If
        If (lfc) Then
           If (first_fc) Then
             facevarx1(:,:,:,:,2) = 0.
             facevary1(:,:,:,:,2) = 0.
             facevarz1(:,:,:,:,2) = 0.
             first_fc = .False.
           End If
        End If

        idest = 2
        Call mpi_amr_get_remote_block(mype,parent_pe,parent_lb,        & 
                                      idest,iopt,lcc,lfc,lec,lnc,      & 
                                      nlayers0x,nlayers0y,nlayers0z)

!-------Do guardcell filling for lb's parent from any surrounding blocks at 
!-------the same refinement level as this parent.
!-------Diagonal elements are required to ensure that all cells are filled
!-------correctly when icoord is non-zero.
        iblock=2
        icoord_loc = 0
        ldiag_loc = .True.

        If (spherical_pm) Then
        If (parent_pe == mype) Then
          pbnd_box(:,2) = bnd_box(:,2,parent_lb)
        Else
          ijk = 1
          Do ij=strt_buffer,last_buffer
            If (laddress(1,ij) == parent_lb                            & 
            .And.laddress(2,ij) == parent_pe) Then
               ijk = ij
               pbnd_box(:,2) = bnd_box(:,2,ijk)
            End If
          End Do
          If (ijk == -1) Then
            Write(*,*)                                                 & 
             'mpi_amr_1blk_guardcell : ijk still -1 : search failed'
            Call amr_abort()
          End If
        End If ! End If (parent_pe == mype_
        ippolar = 0
        If (lsingular_line) Then
         If (abs(pbnd_box(1,2)) < eps) ippolar(1) = -1
         If (abs(pbnd_box(2,2)-pi) < eps) ippolar(2) = 1
        End If  
        End If  ! End If (spherical_pm)

        Call amr_1blk_guardcell_srl(mype,parent_pe,parent_lb,          & 
                                    iblock,iopt,nlayers,psurrblks,     & 
                                    lcc,lfc,lec,lnc,                   & 
                                    icoord_loc,ldiag_loc,              & 
                                    nlayers0x,                         &
                                    nlayers0y,                         &
                                    nlayers0z,ippolar)
        End If  ! End If ( (parent_lb > 0) .And. ...)

!-------Do guardcell filling from coarse neigbors into the current block
        Call mpi_amr_1blk_guardcell_c_to_f( mype,lb,pe,iopt,nlayers,   & 
                                            surrblks,                  & 
                                            lcc,lfc,lec,lnc,           & 
                                            icoord,ldiag,              & 
                                            nlayers0x,                 & 
                                            nlayers0y,                 & 
                                            nlayers0z,ipolar)

        End If  ! End If (lcoarse)

      End If  ! End If (.Not.l_srl_only)

      If (timing_mpi) Then
      timer_amr_1blk_guardcell(2) = timer_amr_1blk_guardcell(2)        &  
                                 +  mpi_wtime() - time2
      time2 = mpi_wtime()
      End If

!-----Do guardcell filling from any surrounding blocks at the same refinement
!-----level as block lb.
      iblock = 1

      Call amr_1blk_guardcell_srl(mype,mype,lb,                        & 
                                  iblock,iopt,nlayers,surrblks,        & 
                                  lcc,lfc,lec,lnc,                     & 
                                  icoord,ldiag,                        & 
                                  nlayers0x,nlayers0y,nlayers0z,       & 
                                  ipolar)


      If (timing_mpi) Then
      timer_amr_1blk_guardcell(3) = timer_amr_1blk_guardcell(3)        & 
                                 +  mpi_wtime() - time2
      End If

      If (iopt == 1) Then
        pcache_pe_u  = pcache_pe
        pcache_blk_u = pcache_blk
      ElseIf (iopt >= 2) Then
        pcache_pe_w  = pcache_pe
        pcache_blk_w = pcache_blk
      End If

      If (timing_mpi) Then
      timer_amr_1blk_guardcell(0) = timer_amr_1blk_guardcell(0)        & 
                                 +  mpi_wtime() - time1
      End If

      Return
      End Subroutine amr_1blk_guardcell



