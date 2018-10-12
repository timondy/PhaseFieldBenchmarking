!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_guardcell_srl
!! NAME
!!
!!   amr_1blk_guardcell_srl
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_guardcell_srl(mype,pe,lb,iblock,iopt,nlayers, 
!!                               surrblks,lcc,lfc,lec,lnc,
!!                               icoord,ldiag,nlayers0x, 
!!                               nlayers0y,nlayers0z, 
!!                               ipolar)
!!   Call amr_1blk_guardcell_srl(integer,integer,integer,integer,
!!                                integer,integer, 
!!                               integer array,logical,logical,
!!                                 logical,logical,
!!                               integer,logical,integer, 
!!                               integer,integer, 
!!                               integer array)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in) :: mype local processor number
!!   Integer, Intent(in) :: pe processor address of the selected block
!!   Integer, Intent(in) :: lb local address on proc. pe of the selected block
!!   Integer, Intent(in) :: iblock selects the storage space in data_1blk.fh which is to
!!                                 be used in this Call. If the leaf node is having its
!!                                 guardcells filled then set this to 1, if its parent
!!                                 is being filled set it to 2.
!!   Integer, Intent(in) :: iopt a switch to control which data source is to be used
!!                               iopt=1 will use 'unk'
!!                               iopt=2 will use 'work'
!!   Integer, Intent(in) :: nlayers the number of guard cell layers at each boundary
!!   Integer, Intent(in) :: surrblks(:,:,:,:) the list of addresses of blocks 
!!                                            surrounding block lb
!!   Logical, Intent(in) :: lcc a logical switch controlling whether unk or work data
!!                              is filled
!!   Logical, Intent(in) :: lfc a logical switch controlling whether facevar data
!!                              is filled
!!   Logical, Intent(in) :: lec a logical switch controlling whether unk_e_x(y)(z) data
!!                              is filled
!!   Logical, Intent(in) :: lnc a logical switch controlling whether unk_n data
!!                              is filled
!!   Integer, Intent(in) :: icoord an integer switch used to select which faces of
!!                                 the block are to be considered. If icoord=0 all
!!                                 faces are considered. If icoord=1 only faces perp.
!!                                 to the y-axis are considered, if icoord=2 only faces
!!                                 perp. to the x-axis are considered, and if icoord=3
!!                                 only faces perp. to the z-axis are considered.
!!   Logical, Intent(in) :: ldiag a logical switch which controls whether guardcells
!!                                corresponding to neighbor blocks diagonally opposite
!!                                block edges and corners are filled.
!!   Integer, Intent(in) :: nlayers0x(yz) The number of guardcell layers to fill in the
!!                                        x, y, and z directions
!!   Integer, Intent(in) :: ipolar(2) this is used to give info about whether this block
!!                                    touches a singular line. For example, in spherical
!!                                    coordinates, ipolar=-1 means the block touches the
!!                                    north polar axis, ipolar=+1 the south polar axis,
!!                                    and ipolar=0 that it is not next to the polar axis.
!!
!! INCLUDES
!!
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   workspace
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_1blk_cc_cp_remote
!!   amr_1blk_fc_cp_remote
!!   amr_1blk_ec_cp_remote
!!   amr_1blk_nc_cp_remote
!!   amr_1blk_bcset
!!
!! RETURNS
!!
!!   Upon exit guardcells of block 'lb' are filled with data from surrounding blocks
!!   at the same refinement level.
!!
!! DESCRIPTION
!!
!!
!!   This routine manages the exchange of guard cell information between
!!   blocks required to fill guard cells on block (pe,lb), assuming that 
!!   exchange is only required between blocks at the same refinement level.
!!   The actual exchanges are performed with Calls to the routines 
!!   amr_1blk_cc_cp_remote and amr_1blk_fc_cp_remote.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          July 1998
!!   Modified:     Rick DeVore             February 2001
!!   Modified:     Peter MacNeice          February 2001
!!   Modified:     Kevin Olson for directional guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f

      Subroutine amr_1blk_guardcell_srl(mype,pe,lb,iblock,iopt,        & 
                                        nlayers,surrblks,              &
                                        lcc,lfc,lec,lnc,               & 
                                        icoord,ldiag,nlayers0x,        & 
                                        nlayers0y,nlayers0z,           & 
                                        ipolar)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use workspace

      Use paramesh_interfaces, Only : amr_1blk_cc_cp_remote,           & 
                                      amr_1blk_fc_cp_remote,           & 
                                      amr_1blk_ec_cp_remote,           & 
                                      amr_1blk_nc_cp_remote,           & 
                                      amr_1blk_bcset

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, Intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      Integer, Intent(in) :: surrblks(:,:,:,:)
      Logical, Intent(in) :: lcc,lfc,lec,lnc,ldiag
      Integer, Intent(in) :: nlayers0x,nlayers0y,nlayers0z
      Integer, Intent(in) :: ipolar(2)

!-----Local arrays and variables

      Integer,Save :: remote_blk,remote_pe,remote_type
      Integer :: nguard0,nlayers0,ng0,jfl,jfu,jface
      Integer :: ilays,jlays,klays,id,jd,kd,is,js,ks
      Integer :: ip1,jp1,kp1,ip2,jp2,kp2,ibc,ii,jj,kk
      Integer :: local_blk_type
      Integer :: nblk_ind
      Integer :: ip3, jp3, kp3, ip4, jp4, kp4
      Integer :: ip5, jp5, kp5, ip6, jp6, kp6
      Integer :: ibnd, jbnd, kbnd
      Integer :: jpolar(2)
      Double Precision :: time1
      Double Precision :: time2

!-----Begin Exectuable Code

      If (timing_mpi) Then
         time1 = mpi_wtime()
      End If  ! End If (timing_mpi)

!-----nblk_ind is not needed for any of the cc_XX type Calls in this routine
      nblk_ind = -1

      local_blk_type = surrblks(3,2,2,2)

      If (iopt == 1) Then
         nguard0 = nguard
         nlayers0 = nguard
      Else If(iopt >= 2) Then 
         nguard0 = nguard_work
!----------nlayers is obsolete !!!
         nlayers0 = max(nlayers0x,nlayers0y,nlayers0z)  
      End If

      ng0 = nguard0*npgs

!-----error trapping!
      If (nlayers0 > nguard0) Then
         Write(*,*) ' nguard = ',nguard
         Write(*,*) ' nguard_work = ',nguard_work
         Write(*,*) ' nlayers0 ',nlayers0,' nguard0 ',nguard0,         & 
                    ' iopt ',iopt
         Write(*,*) 'amr_1blk_guardcell_srl : Too many guardcell ',    & 
                    'layers requested to be filled'
         Call amr_abort()
      End If

      jfl = 1
      jfu = nfaces
      If (icoord > 0) Then
         jfl = 1 + 2*(icoord-1)
         jfu = jfl + 1
      End If

!-----cycle through block faces
      Do jface = jfl,jfu

         jpolar = 0

!-------Default array index limits

!-------Range - source indices are initially computed as though there
!-------are no permanent guardcells.
         ilays = nxb
         jlays = nyb*k2d
         klays = nzb*k3d
!-------Starting indices on destination working block
         id = 1 + nguard0 
         jd = 1 + nguard0*k2d
         kd = 1 + nguard0*k3d
!-------Starting indices on source block
         is = nxb
         js = nyb
         ks = nzb

         ip1 = 0
         jp1 = 0
         kp1 = 0
         ip2 = 0
         jp2 = 0
         kp2 = 0
         ip3 = 0
         jp3 = 0
         kp3 = 0
         ip4 = 0
         jp4 = 0
         kp4 = 0
         ip5 = 0
         jp5 = 0
         kp5 = 0
         ip6 = 0
         jp6 = 0
         kp6 = 0

         ibnd = 0
         jbnd = 0
         kbnd = 0
         If (jface == 1) ibnd = -1
         If (jface == 2) ibnd =  1
         If (jface == 3) jbnd = -1
         If (jface == 4) jbnd =  1
         If (jface == 5) kbnd = -1
         If (jface == 6) kbnd =  1

         If (jface == 1) Then
            remote_blk = surrblks(1,1,2,2)
            remote_pe  = surrblks(2,1,2,2)
            remote_type = surrblks(3,1,2,2)
            id   = 1 + nguard0 - nlayers0x
            is   = 1 + nxb - nlayers0x - gc_off_x
            js   = 1
            ks   = 1
            ip3 = 1
            jp2 = 1
            kp2 = 1
            ilays = nlayers0x
            If (lrestrict_in_progress) ip5 = 1
         Else If (jface == 2) Then
            remote_blk = surrblks(1,3,2,2)
            remote_pe  = surrblks(2,3,2,2)
            remote_type  = surrblks(3,3,2,2)
            id   = 1 + nxb + nguard0 
            is   = 1 + gc_off_x
            js   = 1
            ks   = 1
            ip3 = 1
            ip4 = 1
            ip1 = 1
            jp2 = 1
            kp2 = 1
            ilays = nlayers0x
            If (lrestrict_in_progress) ip5 = 1
            If (lrestrict_in_progress) ip6 = 1
         Else If (jface == 3) Then
            remote_blk = surrblks(1,2,1,2)
            remote_pe  = surrblks(2,2,1,2)
            remote_type  = surrblks(3,2,1,2)
            jd   = 1 + nguard0 - nlayers0y
            js   = 1 + nyb - nlayers0y - gc_off_y
            is   = 1
            ks   = 1
            jp3 = 1
          ip2 = 1
          kp2 = 1
          jlays = nlayers0y
          If (lrestrict_in_progress) jp5 = 1
         Else If (jface == 4) Then
          remote_blk = surrblks(1,2,3,2)
          remote_pe  = surrblks(2,2,3,2)
          remote_type  = surrblks(3,2,3,2)
          jd   = 1 + nyb + nguard0 
          js   = 1 + gc_off_y
          is   = 1
          ks   = 1
          jp3 = 1
          jp4 = 1
          jp1 = 1
          ip2 = 1
          kp2 = 1
          jlays = nlayers0y
          If (lrestrict_in_progress) jp5 = 1
          If (lrestrict_in_progress) jp6 = 1
         Else If (jface == 5) Then
          remote_blk = surrblks(1,2,2,1)
          remote_pe  = surrblks(2,2,2,1)
          remote_type  = surrblks(3,2,2,1)
          kd   = 1 + nguard0 - nlayers0z
          ks   = 1 + nzb - nlayers0z - gc_off_z
          is   = 1
          js   = 1
          kp3 = 1
          ip2 = 1
          jp2 = 1
          klays = nlayers0z
          If (lrestrict_in_progress) kp5 = 1
         Else If (jface == 6) Then
          remote_blk = surrblks(1,2,2,3)
          remote_pe  = surrblks(2,2,2,3)
          remote_type  = surrblks(3,2,2,3)
          kd   = 1 + nzb + nguard0
          ks   = 1 + gc_off_z
          is   = 1
          js   = 1
          kp3 = 1
          kp4 = 1
          kp1 = 1
          ip2 = 1
          jp2 = 1
          klays = nlayers0z
          If (lrestrict_in_progress) kp5 = 1
          If (lrestrict_in_progress) kp6 = 1

         End If  ! End If (jface == 1)

!-------Offset source indices by the no. of permanent guardcells
         is = is + ng0
         js = js + ng0*k2d
         ks = ks + ng0*k3d

!-------If a neighbor exists at this blocks refinement level then fill guardcells
!-------from its data.

         If (remote_blk > 0) Then

            If (timing_mpix) Then
               time2 = mpi_wtime()
            End If

            If (jface == 3.and.ipolar(1) == -1) jpolar(1) = -1
            If (jface == 4.and.ipolar(2) == +1) jpolar(2) = +1
            If (lcc) Call amr_1blk_cc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,iopt,          & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            & 
                       nblk_ind,jpolar)

            If (timing_mpix) Then
               timer_amr_1blk_cc_cp_remote(0) =                           & 
               timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
            Else
               timer_amr_1blk_cc_cp_remote(0) = -999.
            End If
   
            If (lfc) Call amr_1blk_fc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            & 
                       ip1,jp1,kp1,                                    & 
                       ip2,jp2,kp2,jface,                              & 
                       nblk_ind,jpolar)

            If (lec) Call amr_1blk_ec_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,                                    & 
                       ip3,jp3,kp3,ip3,jp3,kp3,jface,                  & 
                       nblk_ind)
            If (lnc) Call amr_1blk_nc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,ip3,jp3,kp3,                        & 
                       nblk_ind)


         Else If (remote_blk <= -20) Then
            ibc = remote_blk
            Call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,            & 
                               ibnd,jbnd,kbnd,surrblks)

         End If  ! End If (remote_blk > 0)

      End Do  ! End Do jface = jfl,jfu


!------------------------------------

      If (ldiag) Then

!------------------------------------
         If (ndim >= 2) Then
            If (icoord.ne.3) Then

!-----Now fill from edges along the z axis.

!-----Loop over the 4 corners
            Do jj = 1,3,2
               Do ii = 1,3,2

          jpolar = 0

!---------Reset default index ranges
          klays = nzb*k3d
          kd = 1 + nguard0*k3d
          ks = 1

          ip1 = 0
          jp1 = 0
          kp1 = 0

          ip2 = 1
          jp2 = 1
          kp2 = 1

          ip3 = 1
          jp3 = 1
          kp3 = 0
          kp4 = 0

          ip5 = 1
          jp5 = 1
          kp5 = 0
          ip6 = 0
          jp6 = 0
          kp6 = 0
          If (lrestrict_in_progress) Then
            If (ii == 3 )ip6 = 1
            If (jj == 3 )jp6 = 1
          End If  ! End If (lrestrict_in_progress)
    
          remote_blk = surrblks(1,ii,jj,2)
          remote_pe  = surrblks(2,ii,jj,2)
          remote_type  = surrblks(3,ii,jj,2)

          ilays = nlayers0x
          jlays = nlayers0y

          is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x) + (ii-2)*gc_off_x
          id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard0
          js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)*gc_off_y
          jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard0

          ip4 = mod(ii/2,2)
          jp4 = mod(jj/2,2)

!-------Offset source indices by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

        ibnd = 0
        jbnd = 0
        kbnd = 0
        If (ii == 1) ibnd = -1
        If (ii == 3) ibnd =  1
        If (jj == 1) jbnd = -1
        If (jj == 3) jbnd =  1

!-------If a neighbor exists at this blocks refinement level then fill guardcells
!-------from its data.

        If (remote_blk > 0) Then

          if (timing_mpix) Then
             time2 = mpi_wtime()
          End If

          If (jj == 1.and.ipolar(1) == -1) jpolar(1) = -1
          If (jj == 3.and.ipolar(2) == +1) jpolar(2) = +1

          If (lcc) Call amr_1blk_cc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,iopt,          & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            & 
                       nblk_ind,jpolar)

          If (timing_mpix) Then
            timer_amr_1blk_cc_cp_remote(0) =                           & 
             timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
          Else
            timer_amr_1blk_cc_cp_remote(0) = -999.
          End If

          If (lfc) Call amr_1blk_fc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            &
                       ip1,jp1,kp1,                                    & 
                       ip2,jp2,kp2,0,                                  & 
                       nblk_ind,jpolar)

          If (lec) Call amr_1blk_ec_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,                                    & 
                       0,0,0,0,0,0,0,                                  & 
                       nblk_ind)
          If (lnc) Call amr_1blk_nc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,0,0,0,                              & 
                       nblk_ind)

        Else If (remote_blk <= -20) Then
          ibc = remote_blk
          Call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,             & 
                               ibnd,jbnd,kbnd,surrblks)
        End If

      End Do  ! End Do ii = 1,3,2
      End Do  ! End Do jj = 1,3,2

      End If  ! End If (icoord.ne.3)
      End If  ! End If (ndim >= 2)


      If (ndim == 3) Then
      If (icoord.ne.2) Then

!-----Now fill from edges along the y axis.

!-----Loop over the 4 corners
      Do kk = 1,3,2
      Do ii = 1,3,2

        jpolar = 0    

!---------Reset default index ranges
        jlays = nyb*k2d
        jd = 1 + nguard0*k2d
        js = 1

        ip1 = 0
        jp1 = 0
        kp1 = 0

        ip2 = 1
        jp2 = 1
        kp2 = 1

        ip3 = 1
        jp3 = 0
        kp3 = 1
        jp4 = 0
    
        remote_blk = surrblks(1,ii,2,kk)
        remote_pe  = surrblks(2,ii,2,kk)
        remote_type  = surrblks(3,ii,2,kk)
        ilays = nlayers0x
        klays = nlayers0z

        is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x)  + (ii-2)*gc_off_x
        id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard0
        ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z)  + (kk-2)*gc_off_z
        kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard0

        ip4 = mod(ii/2,2)
        kp4 = mod(kk/2,2)

!---------Offset source indices by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

        ibnd = 0
        jbnd = 0
        kbnd = 0
        If (ii == 1) ibnd = -1
        If (ii == 3) ibnd =  1
        If (kk == 1) kbnd = -1
        If (kk == 3) kbnd =  1

!-------If a neighbor exists at this blocks refinement level then fill guardcells
!-------from its data.

        If (remote_blk > 0) Then

          If (timing_mpix) Then
             time2 = mpi_wtime()
          End If

          If (lcc) Call amr_1blk_cc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,iopt,          & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            & 
                       nblk_ind,jpolar)

          If (timing_mpix) Then
            timer_amr_1blk_cc_cp_remote(0) =                           & 
             timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
          Else
            timer_amr_1blk_cc_cp_remote(0) = -999.
          End If

          If (lfc) Call amr_1blk_fc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            &
                       ip1,jp1,kp1,                                    & 
                       ip2,jp2,kp2,0,                                  & 
                       nblk_ind,jpolar)

          If (lec) Call amr_1blk_ec_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,                                    & 
                       0,0,0,0,0,0,0,                                  & 
                       nblk_ind)

          If (lnc) Call amr_1blk_nc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,0,0,0,                              & 
                       nblk_ind)

        Else If (remote_blk <= -20) Then
          ibc = remote_blk
          Call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,             & 
                               ibnd,jbnd,kbnd,surrblks)
        End If  ! End If (remote_blk > 0)

      End Do  ! End Do ii = 1,3,2
      End Do  ! End Do kk = 1,3,2

      End If  ! End If (icoord.ne.2)

!-----Now fill from edges along the x axis.
      If (icoord.ne.1) Then

!-----Loop over the 4 corners
      Do kk = 1,3,2
      Do jj = 1,3,2
   
        jpolar = 0 

!---------Reset default index ranges
        ilays = nxb
        id = 1 + nguard0
        is = 1

        ip1 = 0
        jp1 = 0
        kp1 = 0

        jp2 = 1
        kp2 = 1
        ip2 = 1
    
        jp3 = 1
        kp3 = 1
        ip3 = 0
        ip4 = 0
    

        remote_blk = surrblks(1,2,jj,kk)
        remote_pe  = surrblks(2,2,jj,kk)
        remote_type  = surrblks(3,2,jj,kk)
        jlays = nlayers0y
        klays = nlayers0z
        js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)*gc_off_y
        jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard0
        ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z) + (kk-2)*gc_off_z
        kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard0

        jp4 = mod(jj/2,2)
        kp4 = mod(kk/2,2)

!-------Offset source indices by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

        ibnd = 0
        jbnd = 0
        kbnd = 0
        If (jj == 1) jbnd = -1
        If (jj == 3) jbnd =  1
        If (kk == 1) kbnd = -1
        If (kk == 3) kbnd =  1

!-------If a neighbor exists at this blocks refinement level then fill guardcells
!-------from its data.

        If (remote_blk > 0) Then

          If (timing_mpix) Then
             time2 = mpi_wtime()
          End If

          If (jj == 1.and.ipolar(1) == -1) jpolar(1) = -1
          If (jj == 3.and.ipolar(2) == +1) jpolar(2) = +1
          If (lcc) Call amr_1blk_cc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,iopt,          & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            & 
                       nblk_ind,jpolar)

          If (timing_mpix) Then
            timer_amr_1blk_cc_cp_remote(0) =                           & 
             timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
          Else
            timer_amr_1blk_cc_cp_remote(0) = -999.
          End If

          If (lfc) Call amr_1blk_fc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            &
                       ip1,jp1,kp1,                                    & 
                       ip2,jp2,kp2,0,                                  & 
                       nblk_ind,jpolar)

          If (lec) Call amr_1blk_ec_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,                                    & 
                       0,0,0,0,0,0,0,                                  & 
                       nblk_ind)

          If (lnc) Call amr_1blk_nc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,0,0,0,                              & 
                       nblk_ind)

                     Else If (remote_blk <= -20) Then
                        ibc = remote_blk
                        Call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,             & 
                               ibnd,jbnd,kbnd,surrblks)
                     End If  ! End If (remote_blk > 0)

                  End Do  ! End jj = 1,3,2
               End Do  ! End kk = 1,3,2

            End If  ! End If (icoord.ne.1)

         End If  ! End If (ndim == 3)

!------------------------------------

!-----Finally fill corners in 3D.

!------------------------------------
         If (ndim == 3) Then

!-----Loop over the 8 corners
            Do kk = 1,3,2
               Do jj = 1,3,2
                  Do ii = 1,3,2

        jpolar = 0 

        remote_blk = surrblks(1,ii,jj,kk)
        remote_pe  = surrblks(2,ii,jj,kk)
        remote_type  = surrblks(3,ii,jj,kk)

        ilays = nlayers0x
        jlays = nlayers0y
        klays = nlayers0z

        is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x) + (ii-2)* gc_off_x
        id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard0
        js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)* gc_off_y
        jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard0
        ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z) + (kk-2)* gc_off_z
        kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard0

        ip1 = 0
        jp1 = 0
        kp1 = 0

        ip4 = mod(ii/2,2)
        jp4 = mod(jj/2,2)
        kp4 = mod(kk/2,2)

        ip2 = 1
        jp2 = 1
        kp2 = 1

        ip3 = 1
        jp3 = 1
        kp3 = 1
    
!-------Offset source indices by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

        ibnd = 0
        jbnd = 0
        kbnd = 0
        If (ii == 1) ibnd = -1
        If (ii == 3) ibnd =  1
        If (jj == 1) jbnd = -1
        If (jj == 3) jbnd =  1
        If (kk == 1) kbnd = -1
        If (kk == 3) kbnd =  1

!-------If a neighbor exists at this blocks refinement level then fill guardcells
!-------from its data.

        If (remote_blk > 0) Then

          If (timing_mpix) Then
             time2 = mpi_wtime()
          End If

          If (jj == 1.and.ipolar(1) == -1) jpolar(1) = -1
          If (jj == 3.and.ipolar(2) == +1) jpolar(2) = +1
          If (lcc) Call amr_1blk_cc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,iopt,          & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            & 
                       nblk_ind,jpolar)
          If (timing_mpix) Then
            timer_amr_1blk_cc_cp_remote(0) =                           & 
             timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
          Else
            timer_amr_1blk_cc_cp_remote(0) = -999.
          End If

          If (lfc) Call amr_1blk_fc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,ilays,jlays,klays,            &
                       ip1,jp1,kp1,                                    & 
                       ip2,jp2,kp2,0,                                  & 
                       nblk_ind,jpolar)

          If (lec) Call amr_1blk_ec_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,                                    & 
                       0,0,0,0,0,0,0,                                  & 
                       nblk_ind)
          If (lnc) Call amr_1blk_nc_cp_remote(                         & 
                       mype,remote_pe,remote_blk,iblock,               & 
                       id,jd,kd,is,js,ks,                              & 
                       ilays,jlays,klays,                              & 
                       ip1,jp1,kp1,0,0,0,                              & 
                       nblk_ind)

                     Else If (remote_blk <= -20) Then

                        ibc = remote_blk
                        Call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,             & 
                               ibnd,jbnd,kbnd,surrblks)
                     End If  ! End If (remote_blk > 0)

                  End Do  ! End Do ii = 1,3,2
               End Do  ! End DO jj = 1,3,2
            End Do  ! End Do kk = 1,3,2

         End If  ! End If (ndim == 3)

!------------------------------------

      End If  ! End If (ldiag)

!------------------------------------

      If (timing_mpi) Then
       timer_amr_1blk_guardcell_srl =  timer_amr_1blk_guardcell_srl    & 
                                + mpi_wtime() - time1
      End If


      Return
      End Subroutine amr_1blk_guardcell_srl
