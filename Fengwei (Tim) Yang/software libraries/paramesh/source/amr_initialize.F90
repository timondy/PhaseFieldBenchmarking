!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_initialize
!! NAME
!!
!!   amr_initialize
!!
!! SYNOPSIS
!!
!!   Call amr_initialize()
!!
!! ARGUMENTS
!!
!!   No arguments.
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
!!   prolong_arrays
!!   timings
!!   paramesh_comm_data
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_1blk_guardcell_reset
!!   amr_prolong_fun_init
!!   amr_bcset_init
!!   amr_abort
!!   amr_set_runtime_parameters
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit all arrays necessary for PARAMESH to function
!!   are allocated according to the settings in the 'amr_runtime_parameters'.  MPI is also
!!   started.
!!
!! DESCRIPTION
!!
!!   This subroutine initializes the PARAMESH amr package. It performs any
!!   initialization required by the package which is application
!!   independent.  The arrays necessary for PARAMESH to function are allocated
!!   according to the settings in the 'amr_runtime_parameters' file.  MPI is also
!!   started.
!!
!!   NOTE : This routine MUST BE the first executed code in your application!!!!
!!
!! AUTHORS
!!
!!   Peter MacNeice and Kevin Olson
!!
!!***     

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


!----------------------------------------------------------------

      Subroutine amr_initialize 

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace
      Use mpi_morton
!      Use timings
      Use prolong_arrays
      Use timings
      Use paramesh_comm_data

      Use paramesh_interfaces, only : amr_1blk_guardcell_reset,        & 
                                      amr_prolong_fun_init,            & 
                                      amr_bcset_init
      use extvbles

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Integer variables
      Integer :: nfield, nprocs, mype, maxprocs
      Integer :: i
      Integer :: ierr

!-----Real variables
      Real    :: test1,test2
      Real    :: real_test_in(2),real_test_out(2)

!-----Double Precision variables
      Double Precision :: time1
!      Double Precision, allocatable,save :: tmpCEG(:,:,:,:,:)
      character*80 CERRMSG

!----------------------------------------------------------------
!-----Begin Executable code section
!----------------------------------------------------------------

!-----START UP MPI

      Call comm_start(maxprocs,nprocs,mype)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!-----Set the type for passing reals using MPI
#ifdef REAL8
      amr_mpi_real = MPI_DOUBLE_PRECISION
#else
      amr_mpi_real = MPI_REAL
#endif

!-----test mpi real communitcation
      real_test_in(1) = 100. + real(mype)
      real_test_in(2) = 1000. + 2.*real(mype)

      Call mpi_real_allreduce(real_test_in, real_test_out, 2,          & 
           amr_mpi_real,                                               & 
           MPI_MAX, MPI_COMM_WORLD, ierr)
      test1 = 100. + real(nprocs-1)
      test2 = 1000. + 2.*real(nprocs-1)
      If (real_test_out(1).ne.test1.or.real_test_out(2).ne.test2) Then
        If (mype == 0) Then
        Write(*,*) 'PARAMESH ERROR: A test of mpi communication of ',  & 
                   'REAL type failed. Possible cause is using -r8 ',   & 
                   'without defining the preprocessor variable REAL8'  & 
            ,' real_test_in ',real_test_in,' real_test_out ',          & 
              real_test_out,' test1  test2 ', test1, test2
        End If  ! End If (mype == 0)
        Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call amr_abort()
      End If  ! End If(real_test_out(1).ne.test1.or.real_test_out(2).ne.test2)

      ! CEG : This is pointless as coded as timing_mpi not set yet.  Removing 'if'
!      If (timing_mpi) Then
         time1 = mpi_wtime()
!      End If  ! End If (timing_mpi)

!-----Read in 'amr_runtime_parameters' file
      Call amr_set_runtime_parameters()

!-----If nvar, nfacevar, nvaredge, and nvarcorn are all zero, then abort with
!-----a warning message.
      If (nvar <= 0 .and.       &
          nfacevar <= 0 .and.   &
          nvaredge <= 0 .and.   &
          nvarcorn <= 0) Then
      If (mype == 0) then
        print *,' PARAMESH ERROR: nvar = 0, nfacevar = 0, nvaredge = 0, '  
        print *,' and nvarcorn = 0.  Reset one of these values to be >= 1 '
        print *,' rerun '
      End If  ! End If (mype == 0)
      Call amr_abort()
      End If  ! End If (nvar ...

!-----Set other variables based on what was just read in from 'amr_runtime_parameters'.
      If (no_permanent_guardcells) Then
         npgs = 0
      Else
         npgs = 1
      End If  ! End If (no_pemanent_guardcells)

      If (var_dt) Then
         advance_all_levels = .true.
      End If  ! End if (var_dt)

      If (curvilinear .and. .not.cartesian_pm) Then
         consv_fluxes = .true.
         consv_flux_densities = .false.
         edge_value_integ = .true.
         edge_value = .false.
      Else If (.not.curvilinear) Then  ! End If (curvilinear .and. .not. cartesian_pm)
         curvilinear_conserve = .false.
         cartesian_pm = .false.
         cylindrical_pm = .false.
         spherical_pm = .false.
         polar_pm = .false.
      End If  ! End If (.not.curvilinear)

!-----Perform other computation necessary for setting up arrays.
      If (ndim == 1) then
         nyb = 1
      End If  ! End if (ndim == 1)
      If (ndim <= 2) Then
         nzb = 1
      End If  ! End If (ndim == 2)

      nbedges = ndim*2**(ndim-1)
      k3d = (ndim-1)/2
      k2d = ndim/2
      k1d = 1
      red_f = 0.25

      If (consv_fluxes) Then
         If (ndim == 3) Then
            red_f = 1.0
         Else If (ndim == 2) Then
            red_f = 0.5
         Else If (ndim == 1) Then
            red_f = 0.25
         End If
      End If  ! End If (consv_fluxes)

      nchild = 2**ndim
      nfaces = 2*ndim

      If (nboundaries < 2*ndim) Then
         nboundaries=2*ndim
      End If  ! End If (nboundaries < 2*ndim)

      nbndvar  = Max(1,nfacevar)
      nbndvare = Max(1,nvaredge)
      nbndvarc = Max(1,nvarcorn)

      nfluxes  = Max(1,nfluxvar)
      nbndmax  = Max(nbndvar,nfluxes)

      nedgevar = Max(nedgevar1,nvaredge)
      nedges   = Max(1,nedgevar)

      maxdim       = Max(nxb,nyb,nzb)
      gc_off_x     = Mod(nxb,2)
      gc_off_y     = Mod(nyb,2)
      gc_off_z     = Mod(nzb,2)
      il_bnd       = 1
      jl_bnd       = 1
      kl_bnd       = 1
      iu_bnd       = nxb+2*nguard*npgs
      ju_bnd       = nyb+2*nguard*npgs*k2d
      ku_bnd       = nzb+2*nguard*npgs*k3d
      il_bndi      = nguard*npgs+1
      iu_bndi      = nguard*npgs+nxb
      jl_bndi      = nguard*npgs*k2d+1
      ju_bndi      = nguard*npgs*k2d+nyb
      kl_bndi      = nguard*npgs*k3d+1
      ku_bndi      = nguard*npgs*k3d+nzb
      len_block    = iu_bnd*ju_bnd*ku_bnd*nvar
      len_blockfx  = (iu_bnd+1)*ju_bnd*ku_bnd
      len_blockfy  = iu_bnd*(ju_bnd+k2d)*ku_bnd
      len_blockfz  = iu_bnd*ju_bnd*(ku_bnd+k3d)
      len_blockex  = iu_bnd*(ju_bnd+k2d)*(ku_bnd+k3d)
      len_blockey  = (iu_bnd+1)*ju_bnd*(ku_bnd+k3d)
      len_blockez  = (iu_bnd+1)*(ju_bnd+k2d)*ku_bnd
      len_blockn   = (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d)
      len_blockfxf = 2*ju_bnd*ku_bnd
      len_blockfyf = iu_bnd*2*ku_bnd
      len_blockfzf = iu_bnd*ju_bnd*2
      il_bnd1      = 1
      jl_bnd1      = 1
      kl_bnd1      = 1
      iu_bnd1      = nxb+2*nguard
      ju_bnd1      = nyb+2*nguard*k2d
      ku_bnd1      = nzb+2*nguard*k3d
      len_block1   = iu_bnd1*ju_bnd1*ku_bnd1*nvar
      len_blockfx1 = (iu_bnd1+1)*ju_bnd1*ku_bnd1 
      len_blockfy1 = iu_bnd1*(ju_bnd1+k2d)*ku_bnd1
      len_blockfz1 = iu_bnd1*ju_bnd1*(ku_bnd1+k3d)
      len_blockex1 = iu_bnd1*(ju_bnd1+k2d)*(ku_bnd1+k3d)
      len_blockey1 = (iu_bnd1+1)*ju_bnd1*(ku_bnd1+k3d)
      len_blockez1 = (iu_bnd1+1)*(ju_bnd1+1)*ku_bnd1 
      len_blockn1  = (iu_bnd1+1)*(ju_bnd1+1)*(ku_bnd1+k3d)
      ilw          = 1
      jlw          = 1
      klw          = 1
      ngw2         = 2*nguard_work
      iuw          = nxb+ngw2*npgs
      juw          = nyb+ngw2*npgs*k2d
      kuw          = nzb+ngw2*npgs*k3d
      len_wblock   = iuw*juw*kuw
      ilw1         = 1
      jlw1         = 1
      klw1         = 1
      iuw1         = nxb+ngw2
      juw1         = nyb+ngw2*k2d
      kuw1         = nzb+ngw2*k3d
      len_wblock1  = iuw1*juw1*kuw1

      If (ndim == 1) Then
         nmax_lays = nxb/2
      End If ! End if (ndim == 1)

      If (ndim == 2) Then
         nmax_lays = Min(nxb/2,nyb/2)
      End If  ! End If (ndim == 2)

      If (ndim == 3) Then
         nmax_lays = Min(nxb/2,nyb/2,nzb/2)
      End If ! End If (ndim == 3)

      mxblks_buf      = Max(1,maxblocks/100)
      maxblocks_alloc = maxblocks * 4

      maxblocksf      = 1+(maxblocks-1)*Min(1,nfacevar)
      maxblocksue     = 1+(maxblocks-1)*Min(1,nvaredge)
      maxblocksn      = 1+(maxblocks-1)*Min(1,nvarcorn)
      maxblocks_gt    = (maxblocks-1)*(1-npgs)+1
      maxblocksf_gt   = (maxblocksf-1)*(1-npgs)+1
      maxblocksue_gt  = (maxblocksue-1)*(1-npgs)+1
      maxblocksn_gt   = (maxblocksn-1)*(1-npgs)+1
      maxblocksfl     = 1+(maxblocks-1)*Min(1,nfluxvar)
      maxblockse      = 1+(maxblocks-1)*Min(1,nedgevar)

!-----Allocate storage for PARAMESH arrays according to settings read in 
!-----from 'amr_runtime_parameters'

!-----Allocate storage for cell-centered data and their support arrays

      If (nvar <= 0) Then

       Allocate(unk(1,1,1,1,1))
       Allocate(interp_mask_unk(1))
       Allocate(interp_mask_unk_res(1))
       Allocate(gcell_on_cc_pointer(1))
       Allocate(gcell_on_cc(1))
       Allocate(int_gcell_on_cc(1))
       Allocate(checkp_on_cc(1))

      Else
!       print*,'About to allocate tmpCG', mype
!       allocate(tmpCG(75), stat=ierr, ERRMSG=cerrmsg)
!       allocate(tmpCEG(nvar, il_bnd:iu_bnd, jl_bnd:ju_bnd, kl_bnd:ku_bnd, maxblocks), stat=ierr)
!       if (ierr.ne.0) then
!          print*, 'tmpCG = ', ierr, mype,cerrmsg
!       else
!          print*, 'tmpCG successful', mype
!       endif
!
!       print*,'About to allocate tmpCEG', mype
!       allocate(tmpCEG(75), stat=ierr, ERRMSG=cerrmsg)
!!       allocate(tmpCEG(nvar, il_bnd:iu_bnd, jl_bnd:ju_bnd, kl_bnd:ku_bnd, maxblocks), stat=ierr)
!       if (ierr.ne.0) then
!          print*, 'tmpCEG = ', ierr, mype,cerrmsg
!       else
!          print*, 'tmpCEG successful', mype
!       endif

       Allocate(           & 
        unk(nvar,          & 
            il_bnd:iu_bnd, & 
            jl_bnd:ju_bnd, & 
            kl_bnd:ku_bnd, & 
            maxblocks), stat=ierr)
       if (ierr.ne.0) print*, 'unk=', ierr,mype,nvar,il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd,maxblocks
       Allocate(unk1(nvar,            &
                     il_bnd1:iu_bnd1, &
                     jl_bnd1:ju_bnd1, & 
                     kl_bnd1:ku_bnd1, &
                     npblks), stat=ierr)
       if (ierr.ne.0) print *,'unk1=', ierr,mype
       Allocate(               & 
         gt_unk(nvar,          & 
                il_bnd:iu_bnd, &
                jl_bnd:ju_bnd, & 
                kl_bnd:ku_bnd, &
                maxblocks_gt), stat=ierr)
       if (ierr.ne.0) print *,'gt_unk=', ierr

       If (var_dt .or. pred_corr) Then
        Allocate(             & 
         t_unk(nvar,          & 
               il_bnd:iu_bnd, & 
               jl_bnd:ju_bnd, & 
               kl_bnd:ku_bnd, & 
               maxblocks), stat=ierr)
       End If ! End If (var_dt .or. pred_corr)

       Allocate(interp_mask_unk(nvar), stat=ierr)
       Allocate(interp_mask_unk_res(nvar), stat=ierr)
       Allocate(gcell_on_cc_pointer(nvar), stat=ierr)
       Allocate(gcell_on_cc(nvar), stat=ierr)
       Allocate(int_gcell_on_cc(nvar), stat=ierr)
       Allocate(checkp_on_cc(nvar), stat=ierr)

      End If  ! End If (nvar <= 0)
 
!-----Allocate and initialize arrays for face variables

      If (nfacevar <= 0) Then

         Allocate(facevarx(1,1,1,1,1), stat=ierr)
         Allocate(facevary(1,1,1,1,1), stat=ierr)
         Allocate(facevarz(1,1,1,1,1), stat=ierr)

      Else

       Allocate(                  & 
        facevarx(nbndvar,         & 
                 il_bnd:iu_bnd+1, & 
                 jl_bnd:ju_bnd,   & 
                 kl_bnd:ku_bnd,   & 
                 maxblocksf))
       Allocate(                    & 
        facevary(nbndvar,           & 
                 il_bnd:iu_bnd,     & 
                 jl_bnd:ju_bnd+k2d, & 
                 kl_bnd:ku_bnd,     & 
                 maxblocksf))
       Allocate(                    & 
        facevarz(nbndvar,           &  
                 il_bnd:iu_bnd,     & 
                 jl_bnd:ju_bnd,     & 
                 kl_bnd:ku_bnd+k3d, & 
                 maxblocksf))
       facevarx(:,:,:,:,:) = 0.
       facevary(:,:,:,:,:) = 0.
       facevarz(:,:,:,:,:) = 0.

       Allocate(facevarx1(nbndvar,           &
                          il_bnd1:iu_bnd1+1, & 
                          jl_bnd1:ju_bnd1,   & 
                          kl_bnd1:ku_bnd1,   &
                          npblks))
       Allocate(facevary1(nbndvar,             &
                          il_bnd1:iu_bnd1,     & 
                          jl_bnd1:ju_bnd1+k2d, & 
                          kl_bnd1:ku_bnd1,     &
                          npblks))
       Allocate(facevarz1(nbndvar,             &
                          il_bnd1:iu_bnd1,     & 
                          jl_bnd1:ju_bnd1,     & 
                          kl_bnd1:ku_bnd1+k3d, &
                          npblks))
       If (no_permanent_guardcells) Then
        Allocate(                     & 
         gt_facevarx(nbndvar,         &
                     il_bnd:iu_bnd+1, &
                     jl_bnd:ju_bnd,   & 
                     kl_bnd:ku_bnd,   &
                     maxblocksf_gt))
        Allocate(                       & 
         gt_facevary(nbndvar,           &
                     il_bnd:iu_bnd,     &
                     jl_bnd:ju_bnd+k2d, & 
                     kl_bnd:ku_bnd,     &
                     maxblocksf_gt))
        Allocate(                       & 
         gt_facevarz(nbndvar,           &
                     il_bnd:iu_bnd,     &
                     jl_bnd:ju_bnd,     & 
                     kl_bnd:ku_bnd+k3d, &
                     maxblocksf_gt))
       Else ! Else for If (no_permanent_guardcells)
        Allocate(                   & 
         gt_facevarx(nbndvar,       &
                     1:2,           &
                     jl_bnd:ju_bnd, & 
                     kl_bnd:ku_bnd, &
                     maxblocksf))
        Allocate(                   & 
         gt_facevary(nbndvar,       &
                     il_bnd:iu_bnd, &
                     1:1+k2d,       & 
                     kl_bnd:ku_bnd, &
                     maxblocksf))
        Allocate(                   & 
         gt_facevarz(nbndvar,       &
                     il_bnd:iu_bnd, &
                     jl_bnd:ju_bnd, & 
                     1:1+k3d,       &
                     maxblocksf))
       End If  ! End If (no_permanent_guardcells)

       If (var_dt .or. pred_corr) Then
        Allocate(                  & 
        tfacevarx(nbndvar,         & 
                  il_bnd:iu_bnd+1, & 
                  jl_bnd:ju_bnd,   & 
                  kl_bnd:ku_bnd,   & 
                  maxblocksf))
        Allocate(                     & 
         tfacevary(nbndvar,           & 
                   il_bnd:iu_bnd,     & 
                   jl_bnd:ju_bnd+k2d, & 
                   kl_bnd:ku_bnd,     & 
                   maxblocksf))
        Allocate(                     & 
         tfacevarz(nbndvar,           & 
                   il_bnd:iu_bnd,     & 
                   jl_bnd:ju_bnd,     & 
                   kl_bnd:ku_bnd+k3d, & 
                   maxblocksf))

       End If  ! End if (var_dt .or. pred_corr)

      End If ! End If (nfacevar > 0)

      Allocate(gcell_on_fc_pointer(3,nbndvar), stat=ierr)
      Allocate(gcell_on_fc(3,nbndvar), stat=ierr)
      Allocate(int_gcell_on_fc(3,nbndvar), stat=ierr)
      Allocate(interp_mask_facex(nbndvar), stat=ierr)
      Allocate(interp_mask_facey(nbndvar), stat=ierr)
      Allocate(interp_mask_facez(nbndvar), stat=ierr)
      Allocate(interp_mask_facex_res(nbndvar), stat=ierr)
      Allocate(interp_mask_facey_res(nbndvar), stat=ierr)
      Allocate(interp_mask_facez_res(nbndvar), stat=ierr)
      Allocate(checkp_on_fc(3,nbndvar), stat=ierr)

!-----Allocate and intialize edge variables

      If (nvaredge <= 0) Then

       Allocate(unk_e_x(1,1,1,1,1), stat=ierr)
       Allocate(unk_e_y(1,1,1,1,1), stat=ierr)
       Allocate(unk_e_z(1,1,1,1,1), stat=ierr)

      Else

       Allocate(                   & 
        unk_e_x(nbndvare,          & 
                il_bnd:iu_bnd,     & 
                jl_bnd:ju_bnd+k2d, & 
                kl_bnd:ku_bnd+k3d, & 
                maxblocksue))
       Allocate(                   & 
        unk_e_y(nbndvare,          & 
                il_bnd:iu_bnd+1,   & 
                jl_bnd:ju_bnd,     & 
                kl_bnd:ku_bnd+k3d, & 
                maxblocksue))
       Allocate(                   & 
        unk_e_z(nbndvare,          & 
                il_bnd:iu_bnd+1,   & 
                jl_bnd:ju_bnd+k2d, & 
                kl_bnd:ku_bnd,     & 
                maxblocksue))
       unk_e_x(:,:,:,:,:)  = 0.
       unk_e_y(:,:,:,:,:)  = 0.
       unk_e_z(:,:,:,:,:)  = 0.

       Allocate(unk_e_x1(nbndvare,            & 
                         il_bnd1:iu_bnd1,     &
                         jl_bnd1:ju_bnd1+k2d, & 
                         kl_bnd1:ku_bnd1+k3d, & 
                         npblks))
       Allocate(unk_e_y1(nbndvare,            & 
                         il_bnd1:iu_bnd1+1,   &
                         jl_bnd1:ju_bnd1,     & 
                         kl_bnd1:ku_bnd1+k3d, & 
                         npblks))
       Allocate(unk_e_z1(nbndvare,            & 
                         il_bnd1:iu_bnd1+1,   &
                         jl_bnd1:ju_bnd1+k2d, & 
                         kl_bnd1:ku_bnd1,     & 
                         npblks))
       Allocate(                      & 
        gt_unk_e_x(nbndvare,          &
                   il_bnd:iu_bnd,     &
                   jl_bnd:ju_bnd+k2d, & 
                   kl_bnd:ku_bnd+k3d, &
                   maxblocksue_gt))
       Allocate(                      & 
        gt_unk_e_y(nbndvare,          &
                   il_bnd:iu_bnd+1,   &
                   jl_bnd:ju_bnd,     & 
                   kl_bnd:ku_bnd+k3d, &
                   maxblocksue_gt))
       Allocate(                      & 
        gt_unk_e_z(nbndvare,          &
                   il_bnd:iu_bnd+1,   & 
                   jl_bnd:ju_bnd+k2d, & 
                   kl_bnd:ku_bnd,     &
                   maxblocksue_gt))

       If (var_dt .or. pred_corr) Then

        Allocate(                    & 
         t_unk_e_x(nbndvare,         & 
                  il_bnd:iu_bnd,     & 
                  jl_bnd:ju_bnd+k2d, & 
                  kl_bnd:ku_bnd+k3d, & 
                  maxblocksue))
        Allocate(                     & 
         t_unk_e_y(nbndvare,          & 
                   il_bnd:iu_bnd+1,   & 
                   jl_bnd:ju_bnd,     & 
                   kl_bnd:ku_bnd+k3d, & 
                   maxblocksue))
        Allocate(                     & 
         t_unk_e_z(nbndvare,          & 
                   il_bnd:iu_bnd+1,   & 
                   jl_bnd:ju_bnd+k2d, & 
                   kl_bnd:ku_bnd,     & 
                   maxblocksue))

       End If ! End If (var_dt .or. pred_corr)

      End If  ! End If (nvaredge > 0)

      Allocate(interp_mask_ec(nbndvare), stat=ierr)
      Allocate(interp_mask_ec_res(nbndvare), stat=ierr)
      Allocate(gcell_on_ec(3,nbndvare), stat=ierr)
      Allocate(gcell_on_ec_pointer(3,nbndvare), stat=ierr)
      Allocate(int_gcell_on_ec(3,nbndvare), stat=ierr)
      Allocate(checkp_on_ec(3,nbndvare), stat=ierr)

!-----Allocate and initialize corner data

      If (nvarcorn <= 0) Then

       Allocate(unk_n(1,1,1,1,1), stat=ierr)

      Else

       Allocate(unk_n(nbndvarc,          & 
                      il_bnd:iu_bnd+1,   &
                      jl_bnd:ju_bnd+k2d, & 
                      kl_bnd:ku_bnd+k3d, & 
                      maxblocksn))
       unk_n(:,:,:,:,:)    = 0.
       Allocate(unk_n1(nbndvarc,            & 
                       il_bnd1:iu_bnd1+1,   &
                       jl_bnd1:ju_bnd1+k2d, & 
                       kl_bnd1:ku_bnd1+k3d, & 
                       npblks))
       Allocate(                    & 
        gt_unk_n(nbndvarc,          &
                 il_bnd:iu_bnd+1,   &
                 jl_bnd:ju_bnd+k2d, & 
                 kl_bnd:ku_bnd+k3d, &
                 maxblocksn_gt))

       If (var_dt .or. pred_corr) Then
        Allocate(                   & 
         t_unk_n(nbndvarc,          & 
                 il_bnd:iu_bnd+1,   &
                 jl_bnd:ju_bnd+k2d, & 
                 kl_bnd:ku_bnd+k3d, & 
                 maxblocksn))
       End If ! End If (var_dt .or. pred_corr)

      End If  ! End If (nvarcorn <= 0)

      Allocate(interp_mask_nc(nbndvarc), stat=ierr)
      Allocate(interp_mask_nc_res(nbndvarc), stat=ierr)
      Allocate(gcell_on_nc(nbndvarc), stat=ierr)
      Allocate(gcell_on_nc_pointer(nbndvarc), stat=ierr)
      Allocate(int_gcell_on_nc(nbndvarc), stat=ierr)
      Allocate(checkp_on_nc(nbndvarc), stat=ierr)

!-----Allocate arrays for variable time stepping support

      Allocate(time_loc(maxblocks_alloc), stat=ierr)
      Allocate(ldtcomplete(maxblocks_alloc), stat=ierr)

!-----Allocate arrays for flux fix-up at refinement jumps

      Allocate(                & 
       flux_x(nfluxes,         &
              1:2,             & 
              jl_bndi:ju_bndi, &
              kl_bndi:ku_bndi, &
              maxblocksfl), stat=ierr)
      Allocate(                & 
       flux_y(nfluxes,         &
              il_bndi:iu_bndi, & 
              1:2,             &
              kl_bndi:ku_bndi, &
              maxblocksfl), stat=ierr)
      Allocate(                & 
       flux_z(nfluxes,         &
              il_bndi:iu_bndi, & 
              jl_bndi:ju_bndi, &
              1:2,             &
              maxblocksfl), stat=ierr)
      Allocate(                & 
       tflux_x(nfluxes,        &
              1:2,             & 
              jl_bndi:ju_bndi, &
              kl_bndi:ku_bndi, &
              maxblocksfl), stat=ierr)
      Allocate(                 & 
       tflux_y(nfluxes,         &
               il_bndi:iu_bndi, & 
               1:2,             &
               kl_bndi:ku_bndi, &
               maxblocksfl), stat=ierr)
      Allocate(                 & 
       tflux_z(nfluxes,         &
               il_bndi:iu_bndi, & 
               jl_bndi:ju_bndi, &
               1:2,             &
               maxblocksfl), stat=ierr)

!-----Allocate arrrays for edge data fix-up at refinement jumps

      Allocate(                       & 
       bedge_facex_y(nedges,          &
                     1:2,             &
                     jl_bnd:ju_bnd+1, & 
                     kl_bnd:ku_bnd+1, &
                     maxblockse), stat=ierr)
      Allocate(                       & 
       bedge_facex_z(nedges,          &
                     1:2,             &
                     jl_bnd:ju_bnd+1, & 
                     kl_bnd:ku_bnd+1, &
                     maxblockse), stat=ierr)
      Allocate(                       & 
       bedge_facey_x(nedges,          &
                     il_bnd:iu_bnd+1, &
                     1:2,             & 
                     kl_bnd:ku_bnd+1, &
                     maxblockse), stat=ierr)
      Allocate(                       & 
       bedge_facey_z(nedges,          &
                     il_bnd:iu_bnd+1, &
                     1:2,             & 
                     kl_bnd:ku_bnd+1, &
                     maxblockse), stat=ierr)
      Allocate(                       & 
       bedge_facez_x(nedges,          &
                     il_bnd:iu_bnd+1, & 
                     jl_bnd:ju_bnd+1, &
                     1:2,             &
                     maxblockse), stat=ierr)
      Allocate(                       & 
       bedge_facez_y(nedges,          &
                     il_bnd:iu_bnd+1, & 
                     jl_bnd:ju_bnd+1, & 
                     1:2,             &
                     maxblockse), stat=ierr)
      Allocate(                   & 
       recvarx1e(nedges,          &
                 1:2,             &  
                 jl_bnd:ju_bnd+1, & 
                 kl_bnd:ku_bnd+1), stat=ierr)
      Allocate(                   & 
       recvary1e(nedges,          &
                 il_bnd:iu_bnd+1, &
                 1:2,             & 
                 kl_bnd:ku_bnd+1), stat=ierr)
      Allocate(                   & 
       recvarz1e(nedges,          &
                 il_bnd:iu_bnd+1, &
                 jl_bnd:ju_bnd+1, & 
                 1:2), stat=ierr)
      Allocate(                   & 
       recvarx2e(nedges,          &
                 1:2,             &
                 jl_bnd:ju_bnd+1, & 
                 kl_bnd:ku_bnd+1), stat=ierr)
      Allocate(                   & 
       recvary2e(nedges,          &
                 il_bnd:iu_bnd+1, &
                 1:2,             & 
                 kl_bnd:ku_bnd+1), stat=ierr)
      Allocate(                   & 
       recvarz2e(nedges,          &
                 il_bnd:iu_bnd+1, &
                 jl_bnd:ju_bnd+1, & 
                 1:2), stat=ierr)
      Allocate(                        & 
       tbedge_facex_y(nedges,          &
                      1:2,             &
                      jl_bnd:ju_bnd+1, & 
                      kl_bnd:ku_bnd+1, &
                      maxblockse), stat=ierr)
      Allocate(                        & 
       tbedge_facex_z(nedges,          &
                      1:2,             &
                      jl_bnd:ju_bnd+1, & 
                      kl_bnd:ku_bnd+1, &
                      maxblockse), stat=ierr)
      Allocate(                        & 
       tbedge_facey_x(nedges,          &
                      il_bnd:iu_bnd+1, &
                      1:2,             & 
                      kl_bnd:ku_bnd+1, &
                      maxblockse), stat=ierr)
      Allocate(                        & 
       tbedge_facey_z(nedges,          &
                      il_bnd:iu_bnd+1, &
                      1:2,             & 
                      kl_bnd:ku_bnd+1, &
                      maxblockse), stat=ierr)
      Allocate(                        & 
       tbedge_facez_x(nedges,          &
                      il_bnd:iu_bnd+1, & 
                      jl_bnd:ju_bnd+1, &
                      1:2,             &
                      maxblockse), stat=ierr)
      Allocate(                        & 
       tbedge_facez_y(nedges,          &
                      il_bnd:iu_bnd+1, & 
                      jl_bnd:ju_bnd+1, & 
                      1:2,             &
                      maxblockse), stat=ierr)

      Allocate(recvarxf(nfluxes,1:2,jl_bndi:ju_bndi,kl_bndi:ku_bndi), stat=ierr)
      Allocate(recvaryf(nfluxes,il_bndi:iu_bndi,1:2,kl_bndi:ku_bndi), stat=ierr)
      Allocate(recvarzf(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi,1:2), stat=ierr)
      Allocate(bndtempx1(nfluxes,1:2,jl_bndi:ju_bndi,kl_bndi:ku_bndi), stat=ierr)
      Allocate(bndtempy1(nfluxes,il_bndi:iu_bndi,1:2,kl_bndi:ku_bndi), stat=ierr)
      Allocate(bndtempz1(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi,1:2), stat=ierr)

      len_block_bndx = 2*(ju_bndi-jl_bndi+1)*(ku_bndi-kl_bndi+1)
      len_block_bndy = 2*(iu_bndi-il_bndi+1)*(ku_bndi-kl_bndi+1)
      len_block_bndz = 2*(iu_bndi-il_bndi+1)*(ju_bndi-jl_bndi+1)
      len_block_ex   = 2*(ju_bnd+k2d)*(ku_bnd+k3d)
      len_block_ey   = 2*(iu_bnd+1  )*(ku_bnd+k3d)
      len_block_ez   = 2*(iu_bnd+1  )*(ju_bnd+k2d)

!-----Allocate tree data

      maxblocks_tr = 4*maxblocks

!      write (6,*) mchild,maxblocks_tr
      Allocate(neigh(2,mfaces,maxblocks_tr), stat=ierr)
      if (ierr.ne.0) print *,'neigh ierr=',ierr
      Allocate(child(2,mchild,maxblocks_tr), stat=ierr)
      if (ierr.ne.0) print *,'child ierr=',ierr
      child(:,:,:)=0
      Allocate(which_child(maxblocks_tr), stat=ierr)
      if (ierr.ne.0) print *,'which_child ierr=',ierr
      Allocate(parent(2,maxblocks_tr), stat=ierr)
      Allocate(lrefine(maxblocks_tr), stat=ierr)
      Allocate(nodetype(maxblocks_tr), stat=ierr)
      Allocate(empty(maxblocks_tr), stat=ierr)
      Allocate(bflags(mflags,maxblocks_tr), stat=ierr)
      Allocate(newchild(maxblocks_tr), stat=ierr)
      Allocate(derefine(maxblocks_tr), stat=ierr)
      Allocate(refine(maxblocks_tr), stat=ierr)
      Allocate(stay(maxblocks_tr), stat=ierr)
      Allocate(work_block(maxblocks_tr), stat=ierr)
      Allocate(coord(mdim,maxblocks_tr), stat=ierr)
      Allocate(bsize(mdim,maxblocks_tr), stat=ierr)
      Allocate(bnd_box(2,mdim,maxblocks_tr), stat=ierr)
      Allocate(level_cell_sizes(mdim,maxlevels), stat=ierr)
      Allocate(laddress(1:2,1:maxblocks_alloc), stat=ierr)
      Allocate(surr_blks(3,3,1+2*k2d,1+2*k3d,maxblocks_alloc), stat=ierr)
#ifdef SAVE_MORTS
      Allocate(surr_morts(6,3,1+2*k2d,1+2*k3d,maxblocks_alloc), stat=ierr)
#endif
      Allocate(boundary_box(2,mdim,mboundaries), stat=ierr)
      Allocate(boundary_index(mboundaries), stat=ierr)

!-----Allocate workspace data

      if (nvar_work <= 0) then

       Allocate(work(1,1,1,1,1), stat=ierr)
       Allocate(interp_mask_work(1), stat=ierr)
       Allocate(interp_mask_work_res(1), stat=ierr)

      Else

       Allocate(work(ilw:iuw,   &
                     jlw:juw,   &
                     klw:kuw,   &
                     maxblocks, & 
                     nvar_work), stat=ierr)
       Allocate(interp_mask_work(nvar_work), stat=ierr)
       Allocate(interp_mask_work_res(nvar_work), stat=ierr)
       Allocate(recvw(ilw:iuw,jlw:juw,klw:kuw), stat=ierr)
       Allocate(sendw(ilw:iuw,jlw:juw,klw:kuw), stat=ierr)
       Allocate(tempw(ilw:iuw,jlw:juw,klw:kuw), stat=ierr)
       Allocate(work1(ilw1:iuw1,jlw1:juw1,klw1:kuw1,npblks), stat=ierr)
       Allocate(recvw1(ilw1:iuw1,jlw1:juw1,klw1:kuw1,npblks), stat=ierr)
       Allocate(tempw1(ilw1:iuw1,jlw1:juw1,klw1:kuw1), stat=ierr)

      End If  ! End If (nvar_work <= 0)

!-----Allocate morton data

      Allocate(mortonbnd(6,1:3,1:maxblocks), stat=ierr)
      Allocate(laddress_guard(1:2,1:maxblocks_alloc), stat=ierr)
      Allocate(laddress_prol(1:2,1:maxblocks_alloc), stat=ierr)
      Allocate(laddress_flux(1:2,1:maxblocks_alloc), stat=ierr)
      Allocate(laddress_restrict(1:2,1:maxblocks_alloc), stat=ierr)

!-----Allocate prolong_arrays data

      Allocate(prol_dx(il_bnd1:iu_bnd1), stat=ierr)
      Allocate(prol_dy(jl_bnd1:ju_bnd1), stat=ierr)
      Allocate(prol_dz(kl_bnd1:ku_bnd1), stat=ierr)
      Allocate(prol_indexx(2,il_bnd1:iu_bnd1,2), stat=ierr)
      Allocate(prol_indexy(2,jl_bnd1:ju_bnd1,2), stat=ierr)
      Allocate(prol_indexz(2,kl_bnd1:ku_bnd1,2), stat=ierr)
      Allocate(prol_f_dx(il_bnd1:iu_bnd1+1), stat=ierr)
      Allocate(prol_f_dy(jl_bnd1:ju_bnd1+k2d), stat=ierr)
      Allocate(prol_f_dz(kl_bnd1:ku_bnd1+k3d), stat=ierr)
      Allocate(prol_f_indexx(2,il_bnd1:iu_bnd1+1,2), stat=ierr)
      Allocate(prol_f_indexy(2,jl_bnd1:ju_bnd1+k2d,2), stat=ierr)
      Allocate(prol_f_indexz(2,kl_bnd1:ku_bnd1+k3d,2), stat=ierr)
      Allocate(prolw_dx(ilw1:iuw1), stat=ierr)
      Allocate(prolw_dy(jlw1:juw1), stat=ierr)
      Allocate(prolw_dz(klw1:kuw1), stat=ierr)
      Allocate(prolw_indexx(2,ilw1:iuw1,2), stat=ierr)
      Allocate(prolw_indexy(2,jlw1:juw1,2), stat=ierr)
      Allocate(prolw_indexz(2,klw1:kuw1,2), stat=ierr)

!-----Allocate an array for timings

      Allocate(timer_amr_1blk_to_perm(0:1+nvar_work), stat=ierr)

      If (var_dt) Then
       Allocate(                & 
        ttflux_x(nfluxes,       &
                 1:2,           &
                 jl_bnd:ju_bnd, & 
                 kl_bnd:ku_bnd, &
                 maxblocksfl))
       Allocate(                & 
        ttflux_y(nfluxes,       &
                 il_bnd:iu_bnd, & 
                 1:2,           &
                 kl_bnd:ku_bnd, &
                 maxblocksfl))
       Allocate(                & 
        ttflux_z(nfluxes,       &
                 il_bnd:iu_bnd, & 
                 jl_bnd:ju_bnd, &
                 1:2,           &
                 maxblocksfl))
       Allocate(                         & 
        ttbedge_facex_y(nedges,          &
                        1:2,             &
                        jl_bnd:ju_bnd+1, & 
                        kl_bnd:ku_bnd+1, &
                        maxblockse))
       Allocate(                         & 
        ttbedge_facex_z(nedges,          &
                        1:2,             &
                        jl_bnd:ju_bnd+1, & 
                        kl_bnd:ku_bnd+1, &
                        maxblockse))
       Allocate(                         & 
        ttbedge_facey_x(nedges,          &
                        il_bnd:iu_bnd+1, &
                        1:2,             & 
                        kl_bnd:ku_bnd+1, &
                        maxblockse))
       Allocate(                         & 
        ttbedge_facey_z(nedges,          &
                        il_bnd:iu_bnd+1, &
                        1:2,             & 
                        kl_bnd:ku_bnd+1, &
                        maxblockse))
       Allocate(                         & 
        ttbedge_facez_x(nedges,          &
                        il_bnd:iu_bnd+1, & 
                        jl_bnd:ju_bnd+1, &
                        1:2,             &
                        maxblockse))
       Allocate(                         & 
        ttbedge_facez_y(nedges,          &
                        il_bnd:iu_bnd+1, & 
                        jl_bnd:ju_bnd+1, & 
                        1:2,             &
                        maxblockse))
      End If  ! End If (var_dt)


      If (curvilinear) Then
       Allocate(cell_vol(il_bnd1:iu_bnd1, &
                         jl_bnd1:ju_bnd1, &
                         kl_bnd1:ku_bnd1))
       Allocate(cell_area1(il_bnd1:iu_bnd1+1, &
                           jl_bnd1:ju_bnd1,   & 
                           kl_bnd1:ku_bnd1))
       Allocate(cell_area2(il_bnd1:iu_bnd1,     &
                           jl_bnd1:ju_bnd1+k2d, & 
                           kl_bnd1:ku_bnd1))
       Allocate(cell_area3(il_bnd1:iu_bnd1,     &
                           jl_bnd1:ju_bnd1,     & 
                           kl_bnd1:ku_bnd1+k3d))
       Allocate(cell_leng1(il_bnd1:iu_bnd1,     &
                           jl_bnd1:ju_bnd1+k2d, & 
                           kl_bnd1:ku_bnd1+k3d))
       Allocate(cell_leng2(il_bnd1:iu_bnd1+1,   &
                           jl_bnd1:ju_bnd1,     & 
                           kl_bnd1:ku_bnd1+k3d))
       Allocate(cell_leng3(il_bnd1:iu_bnd1+1,   &
                           jl_bnd1:ju_bnd1+k2d, & 
                           kl_bnd1:ku_bnd1))
       Allocate(cell_face_coord1(il_bnd1:iu_bnd1+1))
       Allocate(cell_face_coord2(jl_bnd1:ju_bnd1+k2d))
       Allocate(cell_face_coord3(kl_bnd1:ku_bnd1+k3d))
       Allocate(cell_vol_w(ilw1:iuw1,jlw1:juw1,klw1:kuw1))
      end if  ! End If (curvilinear)

      Allocate(ladd_strt(0:nprocs-1), stat=ierr)
      Allocate(ladd_end(0:nprocs-1), stat=ierr)

! initialize tree data structure
        bsize(:,:)            = -1.
        lrefine(:)            = -1
        nodetype(:)           = -1
        stay(:)               = .TRUE.
        refine(:)             = .FALSE.
        derefine(:)           = .FALSE.
        parent(:,:)           = -1
        child(:,:,:)          = -1
        which_child(:)        = -1
        coord(:,:)            = -1.
        bnd_box(:,:,:)        = -1.
        neigh(:,:,:)          = -1
        empty(:)              = 0
        bflags(:,:)           = -1
        work_block(:)         = 0.
        surr_blks(:,:,:,:,:)  = -1
#ifdef SAVE_MORTS
        surr_morts(:,:,:,:,:) = -1
#endif

!-------initialize solution arrays
        unk(:,:,:,:,:)      = 0.

!-------initialize boundary location arrays for mpi use.
        boundary_box(:,:,nboundaries) = 0.
        boundary_index(nboundaries)   = -1

!-----Initialization required for prolongation routines
      Call amr_prolong_fun_init

!-----Set default values for gcell logical control arrays
      Do i = 1, nvar
        gcell_on_cc_pointer(i) = i
      End Do  ! End Do i = 1,nvar
      gcell_on_cc(:)     = .true.
      int_gcell_on_cc(:) = .true.

      Do i = 1, nfacevar
        gcell_on_fc_pointer(:,i) = i
      End Do  ! End Do i = 1,nfacevar
      gcell_on_fc(:,:)     = .true.
      int_gcell_on_fc(:,:) = .true.

      Do i = 1, nedgevar
        gcell_on_ec_pointer(:,i) = i
      End Do  ! End Do i = 1,nedgevar
      gcell_on_ec(:,:)     = .true.
      int_gcell_on_ec(:,:) = .true.

      Do i = 1, nvarcorn
        gcell_on_nc_pointer(i) = i
      End Do  ! End Do i = 1, nvarcorn
      gcell_on_nc(:)     = .true.
      int_gcell_on_nc(:) = .true.

!-----Set default values for checkpointing 
      checkp_on_cc = .true.
      checkp_on_fc = .true.
      checkp_on_ec = .true.
      checkp_on_nc = .true.

!-----Set default values for interp_masks for prolongation
      interp_mask_unk(:)   = 1
      interp_mask_work(:)  = 1
      interp_mask_facex(:) = 1
      interp_mask_facey(:) = 1
      interp_mask_facez(:) = 1
      interp_mask_ec(:)    = 1
      interp_mask_nc(:)    = 1

!-----Set default values for interp_masks for restriction
      interp_mask_unk_res(:)   = 1
      interp_mask_work_res(:)  = 1
      interp_mask_facex_res(:) = 1
      interp_mask_facey_res(:) = 1
      interp_mask_facez_res(:) = 1
      interp_mask_ec_res(:)    = 1
      interp_mask_nc_res(:)    = 1

!-----Initialize index array defining the variables which
!-----constitute any divergence free fields
      Allocate(i_divf_fc_vars(3,nfield_divf), stat=ierr)
      Do nfield = 1, nfield_divf
        i_divf_fc_vars(:,nfield) = nfield
      End Do  ! End Do nfield = 1, nfield_divf

!-----Initialization required for boundary condition routine
      Call amr_bcset_init
      Call amr_1blk_guardcell_reset

!-----Mark amr_gsurrounding_blks uncalled. This will be set to +1 if
!-----and when amr_gsurrounding_blks is called.
      gsurrblks_set = -1

!-----Initialize grid change marker flag
!-----Initial value = 1 reflects a changed grid.
      grid_changed = 1

!-----Initialize flag to detect if amr_checkpoint_re or amr_refine_derefine
!-----have been called. grid_analysed_mpi = +1 means tha one of them has, so
!-----the mpi version will have control info for any communication dependent
!-----routines.
      grid_analysed_mpi = -1

!-----Initialize mpi communication pattern id
      mpi_pattern_id = 0

!-----set state flags used by mpi communications
      lrestrict_in_progress = .false.
      lprolong_in_progress  = .false.
      lguard_in_progress    = .false.

!-----set state flags used by mpi block boundary info list
!---- -1 is unset. When set, it will have value 100.
      bc_block_neighs_status = -1

!-----Initialize some arrays used in controling mpi communications
      commatrix_recv(:) = 0
      commatrix_send(:) = 0

!-----initialize the instance counter
      instance = 0

! CEG - initialise counters for amr_mg_prolong/restict loops over work variables
      wvar_rest_loop_begin = 1
      wvar_rest_loop_end = nvar_work
      wvar_prol_loop_begin = 1
      wvar_prol_loop_end = nvar_work

      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      If (timing_mpi) Then
         ! CEG added this starting timing call
         call amr_timing_init
         ! Entirely pointless bit
         timer_amr_initialize =  timer_amr_initialize     & 
                                 + mpi_wtime() - time1
      End If  ! End If (timing_mpi)

      Return
      End Subroutine amr_initialize
