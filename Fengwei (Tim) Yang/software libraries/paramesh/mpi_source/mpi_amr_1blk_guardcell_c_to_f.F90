!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_1blk_guardcell_c_to_f
!! NAME
!!
!!   mpi_amr_1blk_guardcell_c_to_f
!!
!! SYNOPSIS
!!
!!   call mpi_amr_1blk_guarcell_c_to_f(mype, lb, pe, iopt, nlayers, surrblks,
!!                                     lcc,lfc,lec,lnc,icoord,ldiag,
!!                                     nlayersx,nlayersy,nlayersz) 
!!
!!   call mpi_amr_1blk_guarcell_c_to_f(integer, integer, integer, integer, integer, 
!!                                     integer,
!!                                     logical, logical, logical, logical, integer, 
!!                                     logical,
!!                                     integer, integer, integer) 
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype           
!!        The local processor number.
!!
!!   integer, intent(in) :: lb
!!        The selected block
!!
!!   integer, intent(in) :: pe
!!        Processor storing the selected block
!!
!!   integer, intent(in) :: iopt     
!!        A switch to control which data source is to be used, iopt = 1 will use 'unk',
!!        'facevarx,y,z', 'unk_e_x,y,z', and 'unk_n'.  iopt >= 2 will use 'work'.
!!
!!   integer, intent(in) :: nlayers
!!        The number of guard cell layers at each boundary.  NOTE: This no longer has
!!        any effect.  Use nlayersx, y and z (see below).
!!
!!   integer, intent(in) :: surrblks(:,:,:,:)
!!        The list of addresses of blocks surrounding block 'lb',
!!
!!   logical, intent(in) :: lcc
!!        A logical switch controlling whether 'unk1' or 'work1' data is filled.
!!        I.e. if true then only 'unk' OR 'work' will have guardcells filled 
!!        depending on the value of 'iopt'.
!!  
!!   logical, intent(in) :: lfc
!!        A logical switch controlling whether 'facevarx1, y1, and z' have 
!!        their guardcells filled.
!!
!!   logical, intent(in) :: lec
!!        A logical switch controlling whether 'unk_e_x1, y1, and z1' have 
!!        their guardcells filled.
!!
!!   logical, intent(in) :: lnc
!!        A logical switch controlling whether 'unk_n' has its guardcells filled.
!!
!!   integer, intent(in) :: icoord
!!        An integer switch used to select which faces of the block are to be considered.
!!        If icoord = 0 all faces are considered.  If icoord = 1 only faces perpendicular
!!        to the x-axis are considered. If icoord = 2 only faces perpendicular to the
!!        y-axis are considered.  If icoord = 3 only faces perpendicular to the
!!        z-axis are considered.
!!
!!   logical, intent(in) :: ldiag
!!        A logical switch which controls whether guardcells corresponding to 
!!        neighbor blocks diagonally opposite block edges and corners are filled.
!!
!!   integer, intent(in) :: nlayersx, nlayersy, nlayersz
!!        Arguments to specify the number of guardcells to fill in each 
!!        coordninate direction.
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
!!   prolong_arrays
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_1blk_cc_prol_gen_work_fun
!!   amr_1blk_cc_prol_gen_unk_fun
!!   amr_1blk_fc_prol_gen_fun
!!   amr_1blk_ec_prol_gen_fun
!!   amr_1blk_nc_prol_gen_fun
!!   amr_1blk_fc_prol_dbz
!!   amr_1blk_fc_clean_divb
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guarcells of the block being operated on 
!!   (which is at a jump in refinement) has its guardcells filled.
!!
!! DESCRIPTION
!!
!!   This routine manages the exchange of guard cell information between
!!   blocks required to fill guard cells on block lb, assuming that 
!!   exchange is only required from surrounding blocks at a coarser 
!!   refinement level.
!!
!!   If you are using an odd number of grid cells, then the interface
!!   condition implicit in this routine is that the finer block uses
!!   data prolonged from its parent on the faces which border coarser
!!   grid blocks.
!!
!! AUTHORS
!!
!!   Peter MacNeice (July 1998) with modifications by 
!!   Kevin Olson (2003) for layered guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_amr_1blk_guardcell_c_to_f(                        & 
                          mype,lb,pe,iopt,nlayers,                     & 
                          surrblks,lcc,lfc,lec,lnc,icoord,ldiag,       & 
                          nlayersx,nlayersy,nlayersz,ipolar)


!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use workspace
      Use prolong_arrays, Only : prol_fc_dbz,                          & 
                                 prol_fc_dbz_ivar,                     & 
                                 prol_fc_dbz_n,                        & 
                                 prol_fc_clean_divb
      Use paramesh_interfaces, Only : amr_1blk_cc_prol_gen_work_fun,   & 
                                      amr_1blk_cc_prol_gen_unk_fun,    & 
                                      amr_1blk_fc_prol_gen_fun,        & 
                                      amr_1blk_ec_prol_gen_fun,        & 
                                      amr_1blk_nc_prol_gen_fun,        & 
                                      amr_1blk_fc_prol_dbz,            & 
                                      amr_1blk_fc_clean_divb

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      Integer, intent(in) :: surrblks(:,:,:,:)
      Logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      Integer, intent(in) :: nlayersx,nlayersy,nlayersz
      Integer, intent(in) :: ipolar(2)

!-----Local arrays and variables
      Integer,Save :: cparent(2)
      Integer,Save :: awhich_child(1)
      Integer,Save :: cwhich_child

      Integer      :: parent_blk,parent_pe
      Integer,Save :: remote_blk,remote_pe
      Integer      :: i, j, k, ii, jj, kk
      Integer      :: igc_off_x, igc_off_y, igc_off_z
      Integer      :: ioff, joff, koff
      Integer      :: nguard0, nlayers0
      Integer      :: jfl, jfu, jface
      Integer      :: ia, ib, ja, jb, ka, kb
      Integer      :: ionea, ioneb, jonea, joneb, konea, koneb
      Integer      :: ia1, ib1, ja1, jb1, ka1, kb1
      Integer      :: imod, jmod, k_mod, idest
      Integer      :: iv1, iv2, iv3, iprol

      Logical      :: lcoarse_neigh

!-----Begin executable code

      igc_off_x = gc_off_x
      igc_off_y = gc_off_y
      igc_off_z = gc_off_z

!-----Does this block have any neighbors at lower refinement level?
      lcoarse_neigh = .False.
      Do k = 2-k3d,2+k3d
      Do j = 2-k2d,2+k2d
      Do i = 1,3
         If(surrblks(1,i,j,k) < 1 .and. surrblks(1,i,j,k) > -20)       & 
              lcoarse_neigh = .True.
      End Do
      End Do
      End Do
      If (.Not.lcoarse_neigh) Return

!-----Does current block have a parent?
      If (lnew_parent) Then
        cparent(:) = parent(:,lb)
      End If
      
      If (cparent(1) > -1) Then

!-----Get parent's address
       parent_blk = cparent(1)
       parent_pe  = cparent(2)
       cwhich_child = which_child(lb)

!------compute the offset in the parent block appropriate for this child
       ioff = mod(cwhich_child-1,2)*nxb/2
       joff = mod((cwhich_child-1)/2,2)*nyb/2
       koff = mod((cwhich_child-1)/4,2)*nzb/2

      End If  ! End If (cparent(1) > -1)

      nguard0 = nguard
      nlayers0 = nguard

      if (iopt >= 2) Then 
           nguard0 = nguard_work
           nlayers0 = nlayers
      End If

!-----First deal with block's regular faces

      jfl = 1
      jfu = nfaces
      If (icoord > 0) Then
        jfl = 1 + 2*(icoord-1)
        jfu = jfl + 1
      End If

!-----cycle through block faces
      Do jface = jfl,jfu

        ia = 1 + nguard0
        ib = nxb + nguard0
        ja = 1 + nguard0*k2d
        jb = nyb + nguard0*k2d
        ka = 1 + nguard0*k3d
        kb = nzb + nguard0*k3d
        ionea = 0
        ioneb = 1
        jonea = 0
        joneb = k2d
        konea = 0
        koneb = k3d

        If (jface == 1) Then
          remote_blk = surrblks(1,1,2,2)
          remote_pe  = surrblks(2,1,2,2)
          ia   = 1 + nguard0 - nlayersx
          ib   = nguard0 + igc_off_x
        Elseif (jface == 2) Then
          remote_blk = surrblks(1,3,2,2)
          remote_pe  = surrblks(2,3,2,2)
          ia   = 1 + nxb + nguard0 - igc_off_x
          ib   = nxb + 2*nguard0
        Elseif (jface == 3) Then
          remote_blk = surrblks(1,2,1,2)
          remote_pe  = surrblks(2,2,1,2)
          ja   = 1 + nguard0 - nlayersy
          jb   = nguard0 + igc_off_y
        Elseif (jface == 4) Then
          remote_blk = surrblks(1,2,3,2)
          remote_pe  = surrblks(2,2,3,2)
          ja   = 1 + nyb + nguard0 - igc_off_y
          jb   = nyb + 2*nguard0
        Elseif (jface == 5) Then
          remote_blk = surrblks(1,2,2,1)
          remote_pe  = surrblks(2,2,2,1)
          ka   = 1 + nguard0 - nlayersz
          kb   = nguard0 + igc_off_z
        Elseif (jface == 6) Then
          remote_blk = surrblks(1,2,2,3)
          remote_pe  = surrblks(2,2,2,3)
          ka   = 1 + nzb + nguard0 - igc_off_z
          kb   = nzb + 2*nguard0
        End If  ! End If (jface == 1)

!-------If a neighbor exists at this blocks refinement level then 
!-------fill guardcells from its data.
        If (remote_blk > -20 .and. remote_blk < 0) Then

!---------Interpolate(prolongate) data from the parent to the child
          If (iopt == 1 .and. lcc) then
             Call amr_1blk_cc_prol_gen_unk_fun(                        & 
                            unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb,         & 
                                        1,ioff,joff,koff,mype,         & 
                                        lb,parent_pe,parent_blk)
          Elseif (iopt >= 2) then
             Call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2),        & 
                                        ia,ib,ja,jb,ka,kb,             & 
                                        1,ioff,joff,koff,mype,         & 
                                        lb,parent_pe,parent_blk,       & 
                                        interp_mask_work(iopt-1))
          End If

          If (lfc.and.iopt == 1) then

            If (prol_fc_dbz) Then
             Do iprol = 1, prol_fc_dbz_n
                iv1 = prol_fc_dbz_ivar(1,iprol)
                iv2 = prol_fc_dbz_ivar(2,iprol)
                iv3 = prol_fc_dbz_ivar(3,iprol)

                call amr_1blk_fc_prol_dbz(                             & 
                     facevarx1(:,:,:,:,2),                             & 
                     facevary1(:,:,:,:,2),                             & 
                     facevarz1(:,:,:,:,2),                             & 
                     nfacevar,                                         & 
                     iv1,iv2,iv3,                                      & 
                     ia+ionea,ib+ioneb,                                & 
                     ja+jonea,jb+joneb,                                & 
                     ka+konea,kb+koneb,1,                              & 
                     ioff,joff,koff,                                   & 
                     mype,lb,parent_pe,parent_blk                      & 
                     )

             End Do
            End If

            If (jface == 2) ionea = 1 
            If (jface == 1) ioneb = 0
            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarx1(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,    & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia+ionea,ib+ioneb,ja,jb,ka,kb,1,             & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,1)

            If (ndim >= 2) Then

            If (jface == 4) jonea = 1
            If (jface == 3) joneb = 0
            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevary1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,  & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia,ib,ja+jonea,jb+joneb,ka,kb,1,             & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,2)
            End If  ! End If (ndim >= 2)

            If (ndim == 3) Then

            If (jface == 6) konea = 1
            If (jface == 5) koneb = 0
            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarz1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,      & 
                       kl_bnd1:ku_bnd1+k3d,2),                         & 
                          ia,ib,ja,jb,ka+konea,kb+koneb,1,             & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,3)
!-----------clean_divb 1
            If (prol_fc_clean_divb) Then
               idest = 2
               Call amr_1blk_fc_clean_divb(                            & 
                  nfacevar,                                            & 
                  ia,ib,ja,jb,ka,kb,                                   & 
                  0, 0, 0, 0, 0, 0,                                    & 
                  idest,ioff,joff,koff,                                & 
                  mype,lb,parent_pe,parent_blk                         & 
                  )
            End If

            End If  ! End If (ndim == 3)

          endif

          If (ndim > 1) Then
          If (lec .and. iopt == 1) Then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d
            If (jface == 2) ionea = 1
            If (jface == 1) ioneb = 0
            If (jface == 4) jonea = 1
            If (jface == 3) joneb = 0
            If (jface == 6) konea = 1
            If (jface == 5) koneb = 0

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,   &  
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,1)

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,     & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,2)

            If (ndim == 3) Then

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
                      kl_bnd1:ku_bnd1,2),                              & 
                          ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
                          ioff,joff,koff,mype,3)
            End If  ! End If (ndim == 3)

          End If  ! End If (lec .and. iopt == 1)

          End If ! End If (ndim > 1)

          If (lnc .and. iopt == 1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            If (jface == 2) ionea = 1
            If (jface == 1) ioneb = 0
            If (jface == 4) jonea = 1
            If (jface == 3) joneb = 0
            If (jface == 6) konea = 1
            If (jface == 5) koneb = 0
            Call amr_1blk_nc_prol_gen_fun(                             & 
             unk_n1(:nvarcorn,il_bnd1:iu_bnd1+1,                       & 
                    jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d,2),        & 
                          ia+ionea,ib+ioneb,ja+jonea,jb+joneb,         & 
                          ka+konea,kb+koneb,1,                         & 
                          ioff,joff,koff,mype)

          End If  ! End If (lnc .and. iopt == 1)

        End If  ! End If (remote_blk > -20 .and. remote_blk < 0)

      End Do  ! End Do jface = jfl,jfu

!-----Now handle diagonal blocks
      If (ldiag) Then

      If (ndim >= 2) Then
      If (icoord .ne. 3) Then

!-----Now fill from edges along the z axis.

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

!-----Loop over the 4 corners
      Do jj = 1,3,2
      Do ii = 1,3,2
    
        remote_blk = surrblks(1,ii,jj,2)
        remote_pe  = surrblks(2,ii,jj,2)

        imod = ii/2
        jmod = jj/2*k2d

        ia = 1 + (nguard0 - nlayersx) +                                & 
                 (nxb + nlayersx - igc_off_x )*imod
        ib = nguard0 + igc_off_x +                                     & 
                 (nxb + nlayersx - igc_off_x )*imod
        ja = 1 + (nguard0 - nlayersy) +                                & 
                 (nyb + nlayersy - igc_off_y )*jmod
        jb = nguard0 + igc_off_y +                                     & 
                (nyb + nlayersy - igc_off_y )*jmod

        ia1 = imod
        ib1 = imod
        ja1 = jmod
        jb1 = jmod

        If (remote_blk > -20 .and. remote_blk < 0) then

!---------interpolate(prolongate) data from the parent to the child
          If (iopt == 1 .and. lcc) Then
             Call amr_1blk_cc_prol_gen_unk_fun(                        & 
                            unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb,         & 
                                          1,ioff,joff,koff,mype,       & 
                                          lb,parent_pe,parent_blk)
          Elseif (iopt >= 2) Then
             Call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2),        & 
                                          ia,ib,ja,jb,ka,kb,           & 
                                          1,ioff,joff,koff,mype,       & 
                                          lb,parent_pe,parent_blk,     & 
                                          interp_mask_work(iopt-1))
          End if

          If (lfc .and. iopt == 1) Then

            If (prol_fc_dbz) Then
            Do iprol = 1, prol_fc_dbz_n
               iv1 = prol_fc_dbz_ivar(1,iprol)
               iv2 = prol_fc_dbz_ivar(2,iprol)
               iv3 = prol_fc_dbz_ivar(3,iprol)

               Call amr_1blk_fc_prol_dbz(                              & 
                     facevarx1(:,:,:,:,2),                             & 
                     facevary1(:,:,:,:,2),                             & 
                     facevarz1(:,:,:,:,2),                             & 
                     nfacevar,                                         & 
                     iv1,iv2,iv3,                                      & 
                     ia+ia1,ib+ib1,                                    & 
                     ja+ja1,jb+jb1,                                    & 
                     ka,kb+1,1,                                        & 
                     ioff,joff,koff,                                   & 
                     mype,lb,parent_pe,parent_blk                      & 
                     )

            End Do  ! End Do iprol = 1, prol_fc_dbz_n
            End If  ! End If (prol_fc_dbz)

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarx1(1:nfacevar,il_bnd1:iu_bnd1+1,                   & 
                       jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1,2),             & 
                          ia+ia1,ib+ib1,ja,jb,ka,kb,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,1)

            If (ndim >= 2) Then

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevary1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia,ib,ja+ja1,jb+jb1,ka,kb,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,2)
            End If

            If (ndim == 3) Then

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarz1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,     & 
                       kl_bnd1:ku_bnd1+k3d,2),                         & 
                          ia,ib,ja,jb,ka,kb+1,1,                       & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,3)

!------clean_divb
            If (prol_fc_clean_divb) Then
               idest = 2
               Call amr_1blk_fc_clean_divb(                            & 
                  nfacevar,                                            & 
                  ia,ib,ja,jb,ka,kb,                                   & 
                  0, 0, 0, 0, 0, 0,                                    & 
                  idest,ioff,joff,koff,                                & 
                  mype,lb,parent_pe,parent_blk                         & 
                  )
            End If

            End If  ! End If (ndim == 3)

          End If  ! End If (lfc .and. iopt == 1)

          If (ndim > 1) Then
          If (lec .and. iopt == 1) Then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            If (ii == 1) ioneb = 0
            If (ii == 3) ionea = 1
            If (jj == 1) joneb = 0
            If (jj == 3) jonea = 1


            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,   & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,1)

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,     & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,2)

            If (ndim == 3) Then
            call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
                      kl_bnd1:ku_bnd1,2),                              & 
                          ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
                          ioff,joff,koff,mype,3)
            End If  ! End If (ndim == 3)

          End If  ! End If (lec .and. iopt == 1)

          End If  ! End If (ndim > 1)

          If (lnc .and. iopt == 1) Then

            Call amr_1blk_nc_prol_gen_fun(                             & 
             unk_n1(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,  & 
                    kl_bnd1:ku_bnd1+k3d,2),                            & 
                          ia+imod,ib+imod,ja+jmod,jb+jmod,             & 
                          ka,kb+k3d,1,                                 & 
                          ioff,joff,koff,mype)

          End If  ! End If (lnc .and. iopt == 1)

        End If  ! End If (remote_blk > -20 .and. remote_blk < 0)

      End Do  ! End Do ii = 1,3,2
      End Do  ! End Do jj = 1,3,2

      End If  ! End If (icoord .ne. 3)
      End If  ! End If (ndim >= 2)

      If (ndim == 3) Then
      If (icoord.ne.2) Then

!-----Now fill from edges along the y axis.

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

!-----Loop over the 4 corners
      Do kk = 1,3,2
      Do ii = 1,3,2
    
        remote_blk = surrblks(1,ii,2,kk)
        remote_pe  = surrblks(2,ii,2,kk)

        imod = ii/2
        k_mod = kk/2*k3d

        ia = 1 + (nguard0 - nlayersx) +                                & 
                 (nxb + nlayersx - igc_off_x )*imod
        ib = nguard0 + igc_off_x +                                     & 
                 (nxb + nlayersx - igc_off_x )*imod
        ka = 1 + (nguard0 - nlayersz) +                                & 
                 (nzb + nlayersz - igc_off_z )*k_mod
        kb = nguard0 + igc_off_z +                                     & 
                 (nzb + nlayersz - igc_off_z )*k_mod

        ia1 = imod
        ib1 = imod
        ka1 = k_mod
        kb1 = k_mod

        If (remote_blk > -20 .and. remote_blk < 0) then

!---------interpolate(prolongate) data from the parent to the child
          If (iopt == 1 .and. lcc) Then
             Call amr_1blk_cc_prol_gen_unk_fun(                        & 
                            unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb,         & 
                                         1,ioff,joff,koff,mype,        & 
                                         lb,parent_pe,parent_blk)
          Elseif (iopt >= 2) Then
             Call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2),        & 
                                         ia,ib,ja,jb,ka,kb,            & 
                                         1,ioff,joff,koff,mype,        & 
                                         lb,parent_pe,parent_blk,      & 
                                         interp_mask_work(iopt-1))
          End If  ! End If (iopt == 1 .and. lcc)

          If (lfc.and.iopt == 1) Then

            If (prol_fc_dbz) Then
              Do iprol = 1, prol_fc_dbz_n
               iv1 = prol_fc_dbz_ivar(1,iprol)
               iv2 = prol_fc_dbz_ivar(2,iprol)
               iv3 = prol_fc_dbz_ivar(3,iprol)

               Call amr_1blk_fc_prol_dbz(                              & 
                  facevarx1(:,:,:,:,2),                                & 
                  facevary1(:,:,:,:,2),                                & 
                  facevarz1(:,:,:,:,2),                                & 
                  nfacevar,                                            & 
                  iv1,iv2,iv3,                                         & 
                  ia+ia1,ib+ib1,                                       & 
                  ja,jb+1,                                             & 
                  ka+ka1,kb+kb1,1,                                     & 
                  ioff,joff,koff,                                      & 
                  mype,lb,parent_pe,parent_blk                         & 
                  )
              End Do  ! End Do iprol = 1, prol_fc_dbz_n
            End If  ! End If (prol_fc_dbz

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarx1(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,   & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia+ia1,ib+ib1,ja,jb,ka,kb,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,1)

            If (ndim >= 2) Then

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevary1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia,ib,ja,jb+1,ka,kb,1,                       & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,2)
            End If

            If (ndim == 3) Then

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarz1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,     & 
                       kl_bnd1:ku_bnd1+k3d,2),                         & 
                          ia,ib,ja,jb,ka+ka1,kb+kb1,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,3)
!-----------clean_divb
            If (prol_fc_clean_divb) Then
             idest = 2
             Call amr_1blk_fc_clean_divb(  & 
                  nfacevar, & 
                  ia,ib,ja,jb,ka,kb, & 
                  0, 0, 0, 0, 0, 0, & 
                  idest,ioff,joff,koff, & 
                  mype,lb,parent_pe,parent_blk & 
                  )
            End If  ! End If (prol_fc_clean_divb

            End If  ! End If (ndim == 3)

          End If  ! End If (lfc.and.iopt == 1)

          If (ndim > 1) Then
          If (lec .and. iopt == 1) Then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            If (ii == 1) ioneb = 0
            If (ii == 3) ionea = 1
            If (kk == 1) koneb = 0
            If (kk == 3) konea = 1

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_x1(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,  & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,1)

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_y1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,    & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,2)

            If (ndim == 3) Then

            Call amr_1blk_ec_prol_gen_fun(                             & 
            unk_e_z1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
                      kl_bnd1:ku_bnd1,2),                              & 
                          ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
                          ioff,joff,koff,mype,3)
            End If 

          End If  ! End If (lec .and. iopt == 1)

          End If  ! End If (ndim > 1)

          If (lnc .and. iopt == 1) Then
            Call amr_1blk_nc_prol_gen_fun(                             & 
             unk_n1(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,  & 
                    kl_bnd1:ku_bnd1+k3d,2),                            & 
                          ia+imod,ib+imod,ja,jb+k2d,                   & 
                          ka+k_mod,kb+k_mod,1,                           & 
                          ioff,joff,koff,mype)
          End If

        End If  ! End If (remote_blk > -20 .and. remote_blk < 0)

      End Do  ! End Do ii = 1,3,2
      End Do  ! End Do kk = 1,3,2

      End If  ! End If (icoord.ne.2)

!-----Now fill from edges along the x axis.
      If (icoord.ne.1) Then

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

! Loop over the 4 corners
      Do kk = 1,3,2
      Do jj = 1,3,2
    
        remote_blk = surrblks(1,2,jj,kk)
        remote_pe  = surrblks(2,2,jj,kk)

        jmod = jj/2*k2d
        k_mod = kk/2*k3d

        ja = 1 + (nguard0 - nlayersy) +                                & 
                 (nyb + nlayersy - igc_off_y )*jmod
        jb = nguard0 + igc_off_y +                                     & 
                 (nyb + nlayersy - igc_off_y )*jmod
        ka = 1 + (nguard0 - nlayersz) +                                & 
                 (nzb + nlayersz - igc_off_z )*k_mod
        kb = nguard0 + igc_off_z +                                     & 
                 (nzb + nlayersz - igc_off_z )*k_mod

        ja1 = jmod
        jb1 = jmod
        ka1 = k_mod
        kb1 = k_mod

        If (remote_blk > -20.and.remote_blk < 0) Then

!---------interpolate(prolongate) data from the parent to the child
          If(iopt == 1 .and. lcc) Then
             Call amr_1blk_cc_prol_gen_unk_fun(                        & 
                            unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb,         & 
                                          1,ioff,joff,koff,mype,       & 
                                          lb,parent_pe,parent_blk)
          Elseif (iopt >= 2) Then
             Call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2),        & 
                                           ia,ib,ja,jb,ka,kb,          & 
                                           1,ioff,joff,koff,mype,      & 
                                           lb,parent_pe,parent_blk,    & 
                                           interp_mask_work(iopt-1))
          End If

          If (lfc.and.iopt == 1) Then

            If (prol_fc_dbz) Then
            Do iprol = 1, prol_fc_dbz_n
             iv1 = prol_fc_dbz_ivar(1,iprol)
             iv2 = prol_fc_dbz_ivar(2,iprol)
             iv3 = prol_fc_dbz_ivar(3,iprol)
             Call amr_1blk_fc_prol_dbz(                                & 
                  facevarx1(:,:,:,:,2),                                & 
                  facevary1(:,:,:,:,2),                                & 
                  facevarz1(:,:,:,:,2),                                & 
                  nfacevar,                                            & 
                  iv1,iv2,iv3,                                         & 
                  ia,ib+1,                                             & 
                  ja+ja1,jb+jb1,                                       & 
                  ka,kb,1,                                             & 
                  ioff,joff,koff,                                      & 
                  mype,lb,parent_pe,parent_blk                         &
                  )
            End Do
            End If

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarx1(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,   & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia,ib+1,ja,jb,ka,kb,1,                       & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,1)

            If (ndim >= 2) Then
            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevary1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia,ib,ja+ja1,jb+jb1,ka,kb,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,2)
            End If

            If (ndim == 3) Then
            Call amr_1blk_fc_prol_gen_fun(                             &
             facevarz1(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,     & 
                       kl_bnd1:ku_bnd1+k3d,2),                         & 
                          ia,ib,ja,jb,ka+ka1,kb+kb1,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,3)
!-----------clean_divb
            If (prol_fc_clean_divb) Then
               idest = 2
               Call amr_1blk_fc_clean_divb(                            & 
               nfacevar,                                               & 
               ia,ib,ja,jb,ka,kb,                                      & 
               0, 0, 0, 0, 0, 0,                                       & 
               idest,ioff,joff,koff,                                   & 
               mype,lb,parent_pe,parent_blk                            & 
               )
            End If
            End If  ! End If (ndim == 3)

          End If  ! End (lfc.and.iopt == 1)

          If (ndim > 1) then
          If (lec.and.iopt == 1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            If (jj == 1) joneb = 0
            If (jj == 3) jonea = 1
            If (kk == 1) koneb = 0
            If (kk == 3) konea = 1

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_x1(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,  & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,1)

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_y1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,    & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,2)

            If (ndim == 3) Then
            Call amr_1blk_ec_prol_gen_fun(                             & 
            unk_e_z1(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
                     kl_bnd1:ku_bnd1,2),                               & 
                          ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
                          ioff,joff,koff,mype,3)
            End If 

          End If ! End If (lec.and.iopt == 1)
          End If ! End If (ndim > 1)

          If (lnc.and.iopt == 1) Then
            call amr_1blk_nc_prol_gen_fun(                             & 
             unk_n1(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,  & 
                    kl_bnd1:ku_bnd1+k3d,2),                            & 
                          ia,ib+1,ja+jmod,jb+jmod,                     & 
                          ka+k_mod,kb+k_mod,1,                           & 
                          ioff,joff,koff,mype)
          End If

        End If  ! End If (remote_blk > -20.and.remote_blk < 0)

      End Do  ! End Do kk = 1,3,2
      End Do  ! End Do jj = 1,3,2

      End If  ! End If (icoord.ne.1)
      End If  ! End If (ndim == 3)

!-----Finally fill corners in 3D.
      If (ndim == 3) then

!-----Loop over the 4 corners
      Do kk = 1,3,2
      Do jj = 1,3,2
      Do ii = 1,3,2
    
        remote_blk = surrblks(1,ii,jj,kk)
        remote_pe  = surrblks(2,ii,jj,kk)

        imod = ii/2
        jmod = jj/2
        k_mod = kk/2

        ia = 1 + (nguard0 - nlayersx) +                                & 
                 (nxb + nlayersx - igc_off_x )*imod
        ib = nguard0 + igc_off_x +                                     & 
                 (nxb + nlayersx - igc_off_x )*imod
        ja = 1 + (nguard0 - nlayersy) +                                & 
                 (nyb + nlayersy - igc_off_y )*jmod
        jb = nguard0 + igc_off_y +                                     & 
                 (nyb + nlayersy - igc_off_y )*jmod
        ka = 1 + (nguard0 - nlayersz) +                                & 
                 (nzb + nlayersz - igc_off_z )*k_mod
        kb = nguard0 + igc_off_z +                                     & 
                 (nzb + nlayersz - igc_off_z )*k_mod

        ia1 = imod
        ib1 = imod
        ja1 = jmod
        jb1 = jmod
        ka1 = k_mod
        kb1 = k_mod

        If (remote_blk > -20 .and. remote_blk < 0) Then

! interpolate(prolongate) data from the parent to the child
          If (iopt == 1.and.lcc) Then
             Call amr_1blk_cc_prol_gen_unk_fun(                        & 
                            unk1(:,:,:,:,2),ia,ib,ja,jb,ka,kb,         & 
                                          1,ioff,joff,koff,mype,       & 
                                          lb,parent_pe,parent_blk)

          Elseif (iopt >= 2) Then
             call amr_1blk_cc_prol_gen_work_fun(work1(:,:,:,2),        & 
                                          ia,ib,ja,jb,ka,kb,           & 
                                          1,ioff,joff,koff,mype,       & 
                                          lb,parent_pe,parent_blk,     & 
                                          interp_mask_work(iopt-1))
          End If

          If (lfc.and.iopt == 1) Then

            If (prol_fc_dbz) Then
            Do iprol = 1, prol_fc_dbz_n
               iv1 = prol_fc_dbz_ivar(1,iprol)
               iv2 = prol_fc_dbz_ivar(2,iprol)
               iv3 = prol_fc_dbz_ivar(3,iprol)
               call amr_1blk_fc_prol_dbz(                              & 
                     facevarx1(:,:,:,:,2),                             & 
                     facevary1(:,:,:,:,2),                             & 
                     facevarz1(:,:,:,:,2),                             & 
                     nfacevar,                                         & 
                     iv1,iv2,iv3,                                      & 
                     ia+ia1,ib+ib1,ja+ja1,jb+jb1,ka+ka1,kb+kb1,1,      & 
                     ioff,joff,koff,                                   & 
                     mype,lb,parent_pe,parent_blk                      & 
                     )
            End Do
            End If

            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarx1(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,    & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia+ia1,ib+ib1,ja,jb,ka,kb,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,1)

            If (ndim >= 2) Then
            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevary1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,  & 
                       kl_bnd1:ku_bnd1,2),                             & 
                          ia,ib,ja+ja1,jb+jb1,ka,kb,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,2)
            End If

            If (ndim == 3) Then
            Call amr_1blk_fc_prol_gen_fun(                             & 
             facevarz1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,      & 
                       kl_bnd1:ku_bnd1+k3d,2),                         & 
                          ia,ib,ja,jb,ka+ka1,kb+kb1,1,                 & 
                          ioff,joff,koff,                              & 
                          mype,lb,parent_pe,parent_blk,3)
!-----------clean_divb
            If (prol_fc_clean_divb) Then
              idest = 2
              call amr_1blk_fc_clean_divb(                             & 
                 nfacevar,                                             & 
                 ia,ib,ja,jb,ka,kb,                                    & 
                 0, 0, 0, 0, 0, 0,                                     & 
                 idest,ioff,joff,koff,                                 & 
                 mype,lb,parent_pe,parent_blk                          & 
                 )
            End If
            End If  ! End If (ndim == 3)

          End If  ! End If (lfc.and.iopt == 1)

          If (ndim > 1) Then
          If (lec.and.iopt == 1) Then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            If (ii == 1) ioneb = 0
            If (ii == 3) ionea = 1
            If (jj == 1) joneb = 0
            If (jj == 3) jonea = 1
            If (kk == 1) koneb = 0
            If (kk == 3) konea = 1

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,   & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia,ib,ja+jonea,jb+joneb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,1)

            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,     & 
                      kl_bnd1:ku_bnd1+k3d,2),                          & 
                          ia+ionea,ib+ioneb,ja,jb,ka+konea,kb+koneb,1, & 
                          ioff,joff,koff,mype,2)

            If (ndim == 3) Then
            Call amr_1blk_ec_prol_gen_fun(                             & 
             unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
                      kl_bnd1:ku_bnd1,2),                              & 
                          ia+ionea,ib+ioneb,ja+jonea,jb+joneb,ka,kb,1, & 
                          ioff,joff,koff,mype,3)
            End If 

          End If  ! End If (lec.and.iopt == 1)
          End If  ! End If (ndim > 1)

          if(lnc.and.iopt == 1) Then
            Call amr_1blk_nc_prol_gen_fun(                             & 
             unk_n1(:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,   & 
                    kl_bnd1:ku_bnd1+k3d,2),                            & 
                          ia+imod,ib+imod,ja+jmod,jb+jmod,             & 
                          ka+k_mod,kb+k_mod,1,                           & 
                          ioff,joff,koff,mype)
          End If

        End If  ! End If (remote_blk > -20 .and. remote_blk < 0)

      End Do  ! End Do ii = 1,3,2
      End Do  ! End Do jj = 1,3,2
      End Do  ! End Do kk = 1,3,2

      End If  ! End If (ndim == 3)

      End If  ! End If (ldiag)

      Return
      End Subroutine mpi_amr_1blk_guardcell_c_to_f


