!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_guardcell
!! NAME
!!
!!   amr_guardcell
!!
!! SYNOPSIS
!!
!!   call amr_guarcell(mype, iopt, nlayers)
!!   call amr_guarcell(mype, iopt, nlayers, nlayersx, nlayersy, nlayersz)
!!
!!   call amr_guarcell(integer, integer, integer, 
!!                     optional integer, optional integer, optional integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   integer, intent(in) :: iopt  
!!     Selects whether to fill the guardcells for the arrays unk, 
!!     facevarx, facevary, facevarz, unk_e_x, unk_e_y, unk_e_z, and unk_n 
!!     (if iopt = 1) or work (if iopt 2).
!!
!!   integer, intent(in) :: nlayers  
!!     Dummy variable which does nothing.  Included for consistency with older 
!!     PARAMESH versions.
!!
!!   optional, integer, intent(in) :: nlayersx, nlayersy, nlayersz 
!!     Optional integers which select how many guardcell layers to fill in each
!!     direction.  If these varaibles are not passed in, Then the default is to
!!     fill all allocated guardcell space.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   tree
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!! 
!!   amr_1blk_guardcell_reset
!!   amr_restrict
!!   amr_1blk_guardcell
!!   mpi_amr_comm_setup
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guardcells are filled with data for 
!!   all blocks.
!!
!! DESCRIPTION
!!
!!   Routine to perform guardcell filling in the case that permanent guardcells
!!   are user.  Calling this routine will result in the guardcells being filled
!!   for all child blocks if 'advance_all_levels' is turned OFF or in ALL the blocks
!!   having their guardcells filled if 'advance_all_levels' is turned ON.
!!   The user can choose how many layers of guardcells are filled on each face
!!   of a block by specifying the optional arguments 'nlayersx', 'nlayersy', and
!!   'nlayersz'.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_guardcell(mype,iopt,nlayers,         & 
                               nlayersx,nlayersy,nlayersz)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use tree

      Use paramesh_interfaces, Only : amr_1blk_guardcell_reset, & 
                                      amr_restrict,             & 
                                      amr_1blk_guardcell

      Use paramesh_mpi_interfaces, Only : mpi_amr_comm_setup

      Implicit none

!-----Include statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: mype,iopt,nlayers
      Integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!-----Local variables
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      Logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag,l_force_consist
      Integer :: lb,icoord
      Integer :: id,jd,kd
      Integer :: ilays,jlays,klays
      Integer :: nlayers0x, nlayers0y, nlayers0z, nguard0
      Integer :: nguard_npgs
      Integer :: i,j,k,ivar                                  
      Integer :: ip1,ip2,jp1,jp2,kp1,kp2
      Integer :: ilp,iup,jlp,jup,klp,kup
      Integer :: nprocs, ierr, tag_offset, iempty, iu, ju, ku, iopt0

!------------------------------------
!-----Begin Executable code section
!------------------------------------
!CEG Fixed whitespace indenting


      If (.not.diagonals) Then
         Write(*,*) 'amr_guardcell:  diagonals off'
      End if

      If (iopt == 1) Then

!--------set users selections of guardcell variables
         int_gcell_on_cc = gcell_on_cc
         int_gcell_on_fc = gcell_on_fc
         int_gcell_on_ec = gcell_on_ec
         int_gcell_on_nc = gcell_on_nc

         If (.not.present(nlayersx)) Then
            nlayers0x = nguard
         Else
            nlayers0x = nlayersx
         End if
         If (.not.present(nlayersy)) Then
            nlayers0y = nguard
         Else
            nlayers0y = nlayersy
         End if
         If (.not.present(nlayersz)) Then
            nlayers0z = nguard
         Else
            nlayers0z = nlayersz
         End if
      Else
         If (.not.present(nlayersx)) Then
            nlayers0x = nguard_work
         Else
            nlayers0x = nlayersx
         End if
         If (.not.present(nlayersy)) Then
            nlayers0y = nguard_work
         Else
            nlayers0y = nlayersy
         End if
         If (.not.present(nlayersz)) Then
            nlayers0z = nguard_work
         Else
            nlayers0z = nlayersz
         End if
      End if ! End If (iopt == 1)

      If (iopt == 1) Then
        nguard0 = nguard
      Else
        nguard0 = nguard_work
      End if

      nguard_npgs = nguard*npgs

      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

      If (no_permanent_guardcells) Then

        If (mype == 0) Then
          Write(*,*) 'amr_guardcell call ignored!'
          Write(*,*) 'NO_PERMANENT_GUARDCELLS is defined'
        End if
        Return

      Else  ! no_permanent_guardcells

!------make sure that nlayers and iopt are set consistently.
        If (iopt == 1.and.nlayers.ne.nguard) Then
          If (mype == 0) Then
            Write(*,*) 'PARAMESH ERROR !'
            Write(*,*) 'Error in guardcell - iopt and nlayers'
            Write(*,*) 'are not consistent. For iopt=1 you must'
            Write(*,*) 'set nlayers=nguard.'
          Endif
          Call amr_abort
        Else If (iopt >= 2.and.nlayers > nguard_work) Then
          If (mype == 0) Then
            Write(*,*) 'PARAMESH ERROR !'
            Write(*,*) 'Error in guardcell - iopt and nlayers'
            Write(*,*) 'are not consistent. For iopt>=2 you must'
            Write(*,*) 'set nlayers le nguard_work.'
          endif
          Call amr_abort
        Endif   ! end If (iopt == 1.and.nlayers.ne.nguard)

!-----Reinitialize addresses of cached parent blocks
        Call amr_1blk_guardcell_reset

        lcc = .False.
        lfc = .False.
        lec = .False.
        lnc = .False.
        If (iopt == 1) Then
          If (nvar > 0) lcc = .True.
          If (nfacevar > 0) lfc = .True.
          If (nvaredge > 0) lec = .True.
          If (nvarcorn > 0) lnc = .True.
        Else If (iopt >= 2) Then
          lcc = .True.
        Endif

!-----Restrict solution to parent blocks
        If (.not.advance_all_levels) Then
          iempty = 0
          call amr_restrict(mype,iopt,iempty,.True.)
          call amr_1blk_guardcell_reset
        End if

        l_force_consist = .False.
        If (force_consistency) Then
          l_force_consist = .True.
          If (lfc) Then
            Do lb = 1,lnblocks
              gt_facevarx(:,1,:,:,lb) = facevarx(:,1+nguard_npgs,:,:,lb)
              gt_facevarx(:,2,:,:,lb) =                                   & 
                facevarx(:,nxb+1+nguard_npgs,:,:,lb)
              If (ndim >= 2) Then
                gt_facevary(:,:,1,:,lb) =                                   & 
                   facevary(:,:,1+nguard_npgs*k2d,:,lb)
                gt_facevary(:,:,1+k2d,:,lb) =                               & 
                facevary(:,:,nyb+(1+nguard_npgs)*k2d,:,lb)
              End if
              If (ndim == 3) Then
                gt_facevarz(:,:,:,1,lb) =                                   & 
                   facevarz(:,:,:,1+nguard_npgs*k3d,lb)
                gt_facevarz(:,:,:,1+k3d,lb) =                               & 
                facevarz(:,:,:,nzb+(1+nguard_npgs)*k3d,lb)
              End if
            End do  ! end do lb = 1, lnblocks
          End if   ! end if (lfc)
        End if    ! end if force_consistency)

        tag_offset = 100
        lguard    = .True.
        lprolong  = .False.
        lflux     = .False.
        ledge     = .False.
        lrestrict = .False.
        lfulltree = .False.
        Call mpi_amr_comm_setup(mype,nprocs,                             & 
                              lguard,lprolong,lflux,ledge,lrestrict,   & 
                              lfulltree,                               & 
                              iopt,lcc,lfc,lec,lnc,tag_offset,         & 
                              nlayersx,nlayersy,nlayersz)

        If (lnblocks > 0) Then
          Do lb = 1,lnblocks

            If (nodetype(lb) == 1 .or. nodetype(lb) == 2 .or.                & 
              advance_all_levels) Then

!-------Copy this blocks data into the working block, and fill its guardcells
              ldiag = diagonals
              l_srl_only = .False.                     ! fill srl and coarse
              icoord = 0                               ! fill in all coord directions
        
!print *, 'amr_guardcell', mype
              Call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,             & 
                                lcc,lfc,lec,lnc,                       & 
                                l_srl_only,icoord,ldiag,               & 
                                nlayersx,nlayersy,nlayersz)

              Do k = 1,1+2*k3d
                klp = 0
                kup = 0
                If (k == 1) Then
                  klays = nlayers0z*k3d
                  kd = nguard0*k3d+1 - klays
                  kp1 = 0
                  kp2 = 0
                  If (l_force_consist) kup = k3d
                Else if (k == 2) Then
                  klays = nzb*k3d
                  kd = nguard0*k3d+1
                  kp1 = 0
                  kp2 = k3d
                Else if (k == 3) Then
                  klays = nlayers0z*k3d
                  kd = (nguard0+nzb)*k3d + 1
                  kp1 = k3d
                  kp2 = k3d
                  If (l_force_consist) klp = -k3d
                End if
                ku = kd + klays - k3d
                Do j = 1,1+2*k2d
                  jlp = 0
                  jup = 0
                  If (j == 1) Then
                    jlays = nlayers0y*k2d
                    jd = nguard0*k2d+1 - jlays
                    jp1 = 0
                    jp2 = 0
                    If (l_force_consist) jup = k2d
                  Else if (j == 2) Then
                    jlays = nyb*k2d
                    jd = nguard0*k2d+1
                    jp1 = 0
                    jp2 = k2d
                  Else if (j == 3) Then
                    jlays = nlayers0y*k2d
                    jd = (nguard0+nyb)*k2d + 1
                    jp1 = k2d
                    jp2 = k2d
                    If (l_force_consist) jlp = -k2d
                  End if
                  ju = jd + jlays - k2d
                  Do i = 1,3
                    ilp = 0
                    iup = 0
                    If (i == 1) Then
                      ilays = nlayers0x
                      id = nguard0+1 - ilays
                      ip1 = 0
                      ip2 = 0
                      If (l_force_consist) iup = 1
                    Else if (i == 2) Then
                      ilays = nxb
                      id = nguard0+1
                      ip1 = 0
                      ip2 = 1
                    Else if (i == 3) Then
                      ilays = nlayers0x
                      id = nguard0+nxb + 1
                      ip1 = 1
                      ip2 = 1
                      If (l_force_consist) ilp = -1
                    Else if (i == 3) Then
                    End if
                    iu = id + ilays - 1
  
                    If (i  ==  2 .and. j  ==  1+k2d .and. k  ==  1+k3d) Then 
  
                    Else
  
                      If (lcc) Then
                        If (iopt == 1) Then
                          Do ivar=1,nvar
                            If (int_gcell_on_cc(ivar)) Then
                              unk(ivar,id:iu,jd:ju,kd:ku,lb) =                      & 
                              unk1(ivar,id:iu,jd:ju,kd:ku,1)
                            Endif
                          Enddo
                        Else
                          iopt0 = iopt-1
                          work(id:iu,jd:ju,kd:ku,lb,iopt0) =                         & 
                            work1(id:iu,jd:ju,kd:ku,1)
                        End if   ! end if iopt == 1
                      End if    ! end if (lcc)
  
                      If (lfc) Then
                        Do ivar = 1,nfacevar
                          If (int_gcell_on_fc(1,ivar)) Then
                            facevarx( ivar,id+ip1+ilp:iu+ip2+iup,                     & 
                                    jd:ju,kd:ku,lb) =                          & 
                            facevarx1( ivar,id+ip1+ilp:iu+ip2+iup,                    & 
                                    jd:ju,kd:ku,1)
                          End if
                          If (ndim > 1) Then
                            If (int_gcell_on_fc(2,ivar)) Then
                              facevary( ivar,id:iu,jd+jp1+jlp:ju+jp2+jup,              & 
                                      kd:ku,lb) =                               & 
                              facevary1(ivar,id:iu,jd+jp1+jlp:ju+jp2+jup,              & 
                                       kd:ku,1)
                            End if
                            If (ndim == 3) Then
                              If (int_gcell_on_fc(3,ivar)) Then
                                facevarz( ivar,id:iu,jd:ju,kd+kp1+klp:ku+kp2+kup,lb) =  & 
                                facevarz1(ivar,id:iu,jd:ju,kd+kp1+klp:ku+kp2+kup,1)
                              End if
                            End if  ! end if (ndim == 3)
                          End if   ! end if (ndim > 1)
                        End do    ! end do ivar = 1, nfacevar
                      End if     ! end if (lfc)
  
                      If (lec) Then
                        Do ivar = 1, nvaredge
                          If (ndim > 1) Then
                            If (int_gcell_on_ec(1,ivar)) Then
                              unk_e_x( ivar,id:iu,jd+jp1:ju+jp2,kd+kp1:ku+kp2,lb) =    & 
                              unk_e_x1(ivar,id:iu,jd+jp1:ju+jp2,kd+kp1:ku+kp2,1)
                            End if
                            If (int_gcell_on_ec(2,ivar)) Then
                              unk_e_y( ivar,id+ip1:iu+ip2,jd:ju,kd+kp1:ku+kp2,lb) =    & 
                              unk_e_y1(ivar,id+ip1:iu+ip2,jd:ju,kd+kp1:ku+kp2,1)
                            End If
                            If (ndim == 3) Then
                              If (int_gcell_on_ec(3,ivar)) Then
                                unk_e_z( ivar,id+ip1:iu+ip2,jd+jp1:ju+jp2,kd:ku,lb) =   & 
                                unk_e_z1(ivar,id+ip1:iu+ip2,jd+jp1:ju+jp2,kd:ku,1)
                              End If
                            End If ! End If (ndim == 3)
                          End If  ! End If (ndim > 1)
                        End Do   ! End Do ivar = 1, nvaredge
                      End If    ! End If (lec)
 
                      If (lnc) Then
                        Do ivar = 1, nvarcorn
                          If (int_gcell_on_nc(ivar)) Then
                            unk_n( ivar,                                              & 
                            id+ip1:iu+ip2,                                        & 
                            jd+jp1:ju+jp2,                                        & 
                            kd+kp1:ku+kp2,lb) =                                   & 
                            unk_n1(ivar,                                         & 
                            id+ip1:iu+ip2,                                & 
                            jd+jp1:ju+jp2,                                & 
                            kd+kp1:ku+kp2,1)
                          End If ! End If (int_gcell_on_nc(ivar))
                        End Do  ! End Do ivar = 1, nvarcorn
                      End If   ! Dnd If (lnc)

                    End If     ! end if (i  ==  2 .and. j  ==  1+k2d .and. k  ==  1+k3d)

                  End Do  ! End Do i = 1,3
                End Do  ! End Do j = 1,1+2*k2d
              End Do  ! End Do k = 1,1+2*k3d

            End If  ! End If (nodetype(lb) == 1 .or. nodetype(lb) == 2 .or. 
              !         advance_all_levels) 

          End Do  ! End Do lb = 1, lnblocks
        End if  ! End If (lnblocks >= 1)

!-----reinitialize addresses of cached parent blocks
        Call amr_1blk_guardcell_reset

!-----reset selections of guardcell variables to default
        int_gcell_on_cc(:) = .True.
        int_gcell_on_fc(:,:) = .True.
        int_gcell_on_ec(:,:) = .True.
        int_gcell_on_nc(:) = .True.

      Endif ! If (no_permanent_guardcells)

      Return
      End Subroutine amr_guardcell





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine pf_guardcell(mype,iopt,mglevel,         & 
                               nlayersx,nlayersy,nlayersz)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use tree

      Use paramesh_interfaces, Only : amr_1blk_guardcell_reset, & 
                                      amr_restrict,             & 
                                      amr_1blk_guardcell

      Use paramesh_mpi_interfaces, Only : mpi_amr_comm_setup
      use mpi_morton

      Implicit none

!-----Include statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: mype, iopt, mglevel
      Integer, intent(in), optional :: nlayersx, nlayersy, nlayersz

!-----Local variables
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      Logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag,l_force_consist
      Integer :: lb,icoord, nlayers
      Integer :: id,jd,kd
      Integer :: ilays,jlays,klays
      Integer :: nlayers0x, nlayers0y, nlayers0z, nguard0
      Integer :: nguard_npgs
      Integer :: i,j,k,ivar                                  
      Integer :: ip1,ip2,jp1,jp2,kp1,kp2
      Integer :: ilp,iup,jlp,jup,klp,kup
      Integer :: nprocs, ierr, tag_offset, iempty, iu, ju, ku, iopt0

!------------------------------------
!-----Begin Executable code section
!------------------------------------
      !CEG Set old (almost unused) nlayers to be nguard as no longer an input
      nlayers = nguard

!CEG Fixed whitespace indenting
!      print *, present(nlayersx),present(nlayersy),present(nlayersz)

      If (iopt == 1) Then

!--------set users selections of guardcell variables
         int_gcell_on_cc = gcell_on_cc
         int_gcell_on_fc = gcell_on_fc
         int_gcell_on_ec = gcell_on_ec
         int_gcell_on_nc = gcell_on_nc

         If (.not.present(nlayersx)) Then
            nlayers0x = nguard
         Else
            nlayers0x = nlayersx
         End if
         If (.not.present(nlayersy)) Then
            nlayers0y = nguard
         Else
            nlayers0y = nlayersy
         End if
         If (.not.present(nlayersz)) Then
            nlayers0z = nguard
         Else
            nlayers0z = nlayersz
         End if
      Else
         If (.not.present(nlayersx)) Then
            nlayers0x = nguard_work
         Else
            nlayers0x = nlayersx
         End if
         If (.not.present(nlayersy)) Then
            nlayers0y = nguard_work
         Else
            nlayers0y = nlayersy
         End if
         If (.not.present(nlayersz)) Then
            nlayers0z = nguard_work
         Else
            nlayers0z = nlayersz
         End if
      End if ! End If (iopt == 1)

      If (iopt == 1) Then
        nguard0 = nguard
      Else
        nguard0 = nguard_work
      End if

!      print *, nlayers0x,nlayers0y,nlayers0z

      nguard_npgs = nguard*npgs

      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

!------make sure that nlayers and iopt are set consistently.
        If (iopt == 1.and.nlayers.ne.nguard) Then
          If (mype == 0) Then
            Write(*,*) 'PARAMESH ERROR !'
            Write(*,*) 'Error in guardcell - iopt and nlayers'
            Write(*,*) 'are not consistent. For iopt=1 you must'
            Write(*,*) 'set nlayers=nguard.'
          Endif
          Call amr_abort
        Else If (iopt >= 2.and.nlayers > nguard_work) Then
          If (mype == 0) Then
            Write(*,*) 'PARAMESH ERROR !'
            Write(*,*) 'Error in guardcell - iopt and nlayers'
            Write(*,*) 'are not consistent. For iopt>=2 you must'
            Write(*,*) 'set nlayers le nguard_work.'
          endif
          Call amr_abort
        Endif   ! end If (iopt == 1.and.nlayers.ne.nguard)

!-----Reinitialize addresses of cached parent blocks
        Call amr_1blk_guardcell_reset

        lcc = .True.
        lfc = .false.
        lec = .false.
        lnc = .false.

        if (nfacevar > 0 .or.nvaredge > 0.or.nvarcorn > 0) then
          If (mype == 0) then
            Write(*,*) 'PARAMESH ERROR !'
            print *, 'pf_guardcell seems to have non cellcentred data used which has been streamlined out'
          endif
          Call amr_abort
        Endif  

!-----Restrict solution to parent blocks
        If (.not.advance_all_levels) Then
          iempty = 0
!          if (mype.eq.1) print*,'pf_restrict',mglevel,sum(commatrix_recv(:))
          call pf_restrict(mype,iopt,iempty,.True.,mglevel)
          call amr_1blk_guardcell_reset
        End if

        tag_offset = 100
        lguard    = .True.
        lprolong  = .False.
        lflux     = .False.
        ledge     = .False.
        lrestrict = .False.
        lfulltree = .False.
        Call pf_mpi_amr_comm_setup(mype,nprocs,                             & 
                              lguard,lprolong,lflux,ledge,lrestrict,   & 
                              lfulltree,                               & 
                              iopt,lcc,lfc,lec,lnc,tag_offset,         & 
                              mglevel)

        If (lnblocks > 0) Then
          Do lb = 1,lnblocks

!            If (nodetype(lb) == 1 .or. nodetype(lb) == 2 .or.                & 
!              advance_all_levels) Then
             if (mglevel.eq.lrefine(lb)) then

!-------Copy this blocks data into the working block, and fill its guardcells
              ldiag = diagonals
              l_srl_only = .False.                     ! fill srl and coarse
              icoord = 0                               ! fill in all coord directions
!print *, 'pf_guardcell', mype, lb,nlayers0x,nlayers0y,nlayers0z
        
              Call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,             & 
                                lcc,lfc,lec,lnc,                       & 
                                l_srl_only,icoord,ldiag,               & 
                                nlayers0x,nlayers0y,nlayers0z)

!$OMP PARALLEL PRIVATE(i,j,k,klp,kup,klays,kd,kp1,kp2,ku,jlp,jup,jlays,jd,jp1,jp2,ju,ilp,iup,ilays,id,ip1,ip2,iu)
!$OMP DO SCHEDULE(DYNAMIC,1) 
              Do k = 1,1+2*k3d
                klp = 0
                kup = 0
                If (k == 1) Then
                  klays = nlayers0z*k3d
                  kd = nguard0*k3d+1 - klays
                  kp1 = 0
                  kp2 = 0
                Else if (k == 2) Then
                  klays = nzb*k3d
                  kd = nguard0*k3d+1
                  kp1 = 0
                  kp2 = k3d
                Else if (k == 3) Then
                  klays = nlayers0z*k3d
                  kd = (nguard0+nzb)*k3d + 1
                  kp1 = k3d
                  kp2 = k3d
                End if
                ku = kd + klays - k3d
                Do j = 1,1+2*k2d
                  jlp = 0
                  jup = 0
                  If (j == 1) Then
                    jlays = nlayers0y*k2d
                    jd = nguard0*k2d+1 - jlays
                    jp1 = 0
                    jp2 = 0
                  Else if (j == 2) Then
                    jlays = nyb*k2d
                    jd = nguard0*k2d+1
                    jp1 = 0
                    jp2 = k2d
                  Else if (j == 3) Then
                    jlays = nlayers0y*k2d
                    jd = (nguard0+nyb)*k2d + 1
                    jp1 = k2d
                    jp2 = k2d
                  End if
                  ju = jd + jlays - k2d
                  Do i = 1,3
                    ilp = 0
                    iup = 0
                    If (i == 1) Then
                      ilays = nlayers0x
                      id = nguard0+1 - ilays
                      ip1 = 0
                      ip2 = 0
                    Else if (i == 2) Then
                      ilays = nxb
                      id = nguard0+1
                      ip1 = 0
                      ip2 = 1
                    Else if (i == 3) Then
                      ilays = nlayers0x
                      id = nguard0+nxb + 1
                      ip1 = 1
                      ip2 = 1
                    Else if (i == 3) Then
                    End if
                    iu = id + ilays - 1
  
                    If (i  ==  2 .and. j  ==  1+k2d .and. k  ==  1+k3d) Then 
  
                    Else
  
                      If (lcc) Then
                        If (iopt == 1) Then
                          Do ivar=1,nvar
                            If (int_gcell_on_cc(ivar)) Then
                              unk(ivar,id:iu,jd:ju,kd:ku,lb) =                      & 
                              unk1(ivar,id:iu,jd:ju,kd:ku,1)
                            Endif
                          Enddo
                        Else
                          iopt0 = iopt-1
                          work(id:iu,jd:ju,kd:ku,lb,iopt0) =                         & 
                            work1(id:iu,jd:ju,kd:ku,1)
                        End if   ! end if iopt == 1
                      End if    ! end if (lcc)
  
                    End If     ! end if (i  ==  2 .and. j  ==  1+k2d .and. k  ==  1+k3d)

                  End Do  ! End Do i = 1,3
                End Do  ! End Do j = 1,1+2*k2d
              End Do  ! End Do k = 1,1+2*k3d
!$OMP END DO
!$OMP END PARALLEL

            End If  ! End If (nodetype(lb) == 1 .or. nodetype(lb) == 2 .or. 
              !         advance_all_levels) 

          End Do  ! End Do lb = 1, lnblocks
        End if  ! End If (lnblocks >= 1)

!-----reinitialize addresses of cached parent blocks
        Call amr_1blk_guardcell_reset

!-----reset selections of guardcell variables to default
        int_gcell_on_cc(:) = .True.
        int_gcell_on_fc(:,:) = .True.
        int_gcell_on_ec(:,:) = .True.
        int_gcell_on_nc(:) = .True.

      Return
      End Subroutine pf_guardcell



