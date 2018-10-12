!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_block_geometry
!! NAME
!!   amr_block_geometry
!!
!! SYNOPSIS
!!
!!   Call amr_block_geometry (lb,pe)
!!   Call amr_block_geometry (integer,integer)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in) :: lb local block number
!!   
!!   Integer, Intent(in) :: pe processor on which block lb can be found
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
!!
!! CALLS
!!
!!   amr_abort
!!   user_coord_transfm
!!
!! RETURNS
!!
!!   Nothing Returned
!!
!! DESCRIPTION
!!
!!   This routine computes cell volumes, areas and edge lengths
!!   for various grid geometries, for the specified local block lb.
!!
!!   Cartesian :
!!     coord 1      x
!!     coord 2      y
!!     coord 3      z
!!
!!   Cylindrical :  (NOTE: Coordinates are in this order to support a 2-d coordinate
!!                         system which is axisymmetric about the z axis
!!     coord 1      r
!!     coord 2      z
!!     coord 3      theta
!!
!!   Spherical :
!!     coord 1      r
!!     coord 2      theta
!!     coord 3      phi           (azimuthal)
!!
!!   Polar (2D) :
!!     coord 1      r
!!     coord 2      theta
!!     coord 3      z             (has only 1 grid cell in this direction)
!!
!!
!! AUTHORS
!!
!!   Written : Peter MacNeice      December 2001
!!   Cylindrical axisymmetric added by Sergey Pancheshnyi.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_block_geometry(lb,pe)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use workspace

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in) :: lb,pe

!-----Local arrays and variables

      Real :: cell_face_coord1w(ilw1:iuw1+1)
      Real :: cell_face_coord2w(jlw1:juw1+k2d)
      Real :: cell_face_coord3w(klw1:kuw1+k3d)
      Real :: del
      Real :: cell_vol_1  ,cell_vol_2  ,cell_vol_3
      Real :: cell_area1_1,cell_area1_2,cell_area1_3
      Real :: cell_area2_1,cell_area2_2,cell_area2_3
      Real :: cell_area3_1,cell_area3_2,cell_area3_3
      Real :: cell_leng1_1,cell_leng1_2,cell_leng1_3
      Real :: cell_leng2_1,cell_leng2_2,cell_leng2_3
      Real :: cell_leng3_1,cell_leng3_2,cell_leng3_3
      Real :: dx, dy, dz
      Real :: xleft,yleft,zleft
      Real :: eps
      Real,save :: cbnd_box(2,3)

      Integer :: ierr_trap
      Integer :: mype,ierr
      Integer :: lb0,pe0,iloc
      Integer :: i, j, k

      Logical :: lfound

!--------Begin Executable Code

         eps = Tiny(eps)

         Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)

         ierr_trap = 0
         If (cartesian_pm) Then
            ierr_trap = ierr_trap + 1
         End If
         If (spherical_pm) Then
            ierr_trap = ierr_trap + 1
         End If
         If (cylindrical_pm) Then
            ierr_trap = ierr_trap + 1
         End If
         If (polar_pm) Then
            ierr_trap = ierr_trap + 1
            If (ndim.ne.2) ierr_trap = ierr_trap + 1
         End If
         If (ierr_trap > 1) Then
           Write(*,*) 'Paramesh ERROR : amr_block_geometry. ',         &
                 'Inconsistent choice of curvilinear coord.'
           Call amr_abort()
         End If

!-----------------------

        lfound = .False.
        If (pe == mype) Then
          If (lb <= lnblocks) Then
            lfound = .True.
            lb0 = lb
            pe0 = mype
          ElseIf (lb >= strt_buffer.and.lb <= last_buffer) Then
            lfound = .True.
            lb0 = lb
            pe0 = mype
          End If
        Else
          iloc = strt_buffer
          Do while( (iloc <= last_buffer) .and.                        &
                    (.Not.lfound)  )
            If (laddress(1,iloc) == lb.and.                            &
                laddress(2,iloc) == pe ) lfound = .True.
            If (.Not.lfound) iloc = iloc + 1
          End Do
          If (lfound) Then
            lb0 = iloc
            pe0 = mype
          End If
        End If  ! End If (pe == mype)

        If (.Not.lfound) Then
          Write(*,*) 'amr_block_geometry ERROR : blk ',                &
              lb,pe,' not found on pe ',mype,                          &
              ' strt_buffer:last_buffer ',strt_buffer,last_buffer,     &
              ' laddress ',laddress(:,strt_buffer:last_buffer)
          Call amr_abort()
        End If
        cbnd_box(:,:) = bnd_box(:,:,lb0)

!--------compute coords of cell interfaces
!--------for first coordinate direction
         del = (cbnd_box(2,1)-cbnd_box(1,1))/Real(nxb)
         Do i = il_bnd1,iu_bnd1+1
           cell_face_coord1(i) = cbnd_box(1,1) + del*Real(i-1-nguard)
         End Do
         Do i = ilw1,iuw1+1
           cell_face_coord1w(i) = cbnd_box(1,1)                        &
                                  + del*Real(i-1-nguard_work)
         End Do

         dx = (cbnd_box(2,1)-cbnd_box(1,1))/Real(nxb)

!--------for second coordinate direction
         cell_face_coord2 = 0.
         cell_face_coord2w = 0.
         If (ndim >= 2) Then
         del = (cbnd_box(2,2)-cbnd_box(1,2))/Real(nyb)
         yleft = coord(2,lb0) - bsize(2,lb0)/2.
         Do j = jl_bnd1,ju_bnd1+1
           cell_face_coord2(j) = cbnd_box(1,2) + del*Real(j-1-nguard)
         End Do
         Do j = jlw1,juw1+1
           cell_face_coord2w(j) = cbnd_box(1,2)                        &
                                  + del*Real(j-1-nguard_work)
         End Do
         End If  ! End If (ndim >= 2)

         dy = (cbnd_box(2,2)-cbnd_box(1,2))/Real(nyb)

!--------for third coordinate direction
         cell_face_coord3 = 0.
         cell_face_coord3w = 0.
         If (ndim == 3) Then
         del = (cbnd_box(2,3)-cbnd_box(1,3))/Real(nzb)
         Do k = kl_bnd1,ku_bnd1+1
           cell_face_coord3(k) = cbnd_box(1,3) + del*Real(k-1-nguard)
         End Do
         Do k = klw1,kuw1+1
           cell_face_coord3w(k) = cbnd_box(1,3)                        &
                                  + del*Real(k-1-nguard_work)
         End Do
         End If  ! End If (ndim == 3)

         dz = (cbnd_box(2,3)-cbnd_box(1,3))/Real(nzb)

!--------Apply any user specified coordinate transformation
         call user_coord_transfm(lb0,pe0)


!--------Compute cell volumes
!--------Note the style used here to compute cell_vol. We
!--------specify dependence of cell_vol on coord 1 in cell_vol_1,
!--------specify dependence of cell_vol on coord 2 in cell_vol_2,
!--------specify dependence of cell_vol on coord 3 in cell_vol_3.
!--------This style is used throughout this routine.

!--------First compute cell volumes for use with UNK data structure
         Do k = kl_bnd1,ku_bnd1
         Do j = jl_bnd1,ju_bnd1
         Do i = il_bnd1,iu_bnd1

           cell_vol_1 = 1.
           cell_vol_2 = 1.
           cell_vol_3 = 1.

           If (cartesian_pm) Then
           cell_vol_1 = dx
           If (ndim >= 2)                                              &
             cell_vol_2 = dy
           If (ndim == 3)                                              &
             cell_vol_3 = dz
           If (l2p5d == 1) cell_vol_3 =                                &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cartesian_pm)

           If (spherical_pm) Then
           cell_vol_1 = ( cell_face_coord1(i+1)**3 -                   &
                          cell_face_coord1(i)**3 ) /3.
           If (ndim >= 2)                                              &
             cell_vol_2 = cos( cell_face_coord2(j)   ) -               &
                          cos( cell_face_coord2(j+1) )
           If (ndim == 3)                                              &
              cell_vol_3 = cell_face_coord3(k+k3d) -                   &
                           cell_face_coord3(k)
           If (l2p5d == 1) cell_vol_3 =                                &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (spherical_pm)

           If (cylindrical_pm) Then
           cell_vol_1 = ( cell_face_coord1(i+1)**2 -                   &
                          cell_face_coord1(i)**2 )*.5
           If (ndim >= 2)                                              &
             cell_vol_2 =  cell_face_coord2(j+1) -                     &
                           cell_face_coord2(j)
           If (ndim == 3)                                              &
             cell_vol_3 =  cell_face_coord3(k+k3d) -                   &
                           cell_face_coord3(k)
           If (l2p5d == 1) cell_vol_3 =                                &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cylindrical_pm)

           If (polar_pm) Then
           cell_vol_1 = ( cell_face_coord1(i+1)**2 -                   &
                          cell_face_coord1(i)**2 )*.5
           If (ndim >= 2)                                              &
                cell_vol_2 =  cell_face_coord2(j+1) -                  &
                              cell_face_coord2(j)
           End If  ! End If (polar_pm)

           cell_vol(i,j,k) = max(abs(cell_vol_1 * cell_vol_2           &
                                                * cell_vol_3),         &
                                 eps)

         End Do  ! End Do i = il_bnd1,iu_bnd1
         End Do  ! End Do i = jl_bnd1,ju_bnd1
         End Do  ! End Do i = kl_bnd1,ku_bnd1

!--------now cell volumes for use with WORK data structure
         Do k = klw1,kuw1
         Do j = jlw1,juw1
         Do i = ilw1,iuw1

           cell_vol_1 = 1.
           cell_vol_2 = 1.
           cell_vol_3 = 1.

           If (cartesian_pm) Then
           cell_vol_1 = dx
           If (ndim >= 2)                                              &
             cell_vol_2 = dy
           If (ndim == 3)                                              &
             cell_vol_3 = dz
           If (l2p5d == 1) cell_vol_3 =                                &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cartesian_pm)

           If (spherical_pm) Then
           cell_vol_1 = ( cell_face_coord1w(i+1)**3 -                  &
                          cell_face_coord1w(i)**3 ) /3.
           If (ndim >= 2)                                              &
             cell_vol_2 = cos( cell_face_coord2w(j)   ) -              &
                          cos( cell_face_coord2w(j+1) )
           If (ndim == 3)                                              &
             cell_vol_3 = cell_face_coord3w(k+k3d) -                   &
                          cell_face_coord3w(k)
           If (l2p5d == 1) cell_vol_3 =                                &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (spherical_pm)

           If (polar_pm) Then
           cell_vol_1 = ( cell_face_coord1w(i+1)**2 -                  &
                          cell_face_coord1w(i)**2 )*.5
           If (ndim >= 2)                                              &
             cell_vol_2 =  cell_face_coord2w(j+1) -                    &
                           cell_face_coord2w(j) 
           End If  ! End If (polar_pm)

           If (cylindrical_pm) Then
!----------INT(rdr) * INT(dz) * INT(d theta)
           cell_vol_1 = ( cell_face_coord1w(i+1)**2 -                  &
                          cell_face_coord1w(i)**2 )*.5
           If (ndim >= 2)                                              &
             cell_vol_2 =  cell_face_coord2w(j+1) -                    &
                           cell_face_coord2w(j)
           If (ndim == 3)                                              &
             cell_vol_3 =  cell_face_coord3w(k+k3d) -                  &
                           cell_face_coord3w(k)
           If (l2p5d == 1) cell_vol_3 =                                &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cylindrical_pm)

           cell_vol_w(i,j,k) = max(abs(cell_vol_1 * cell_vol_2         &
                                                  * cell_vol_3),       &
                                   eps)

         End Do  ! End Do i = ilw1,iuw1
         End Do  ! End Do j = jlw1,juw1
         End Do  ! End Do k = klw1,kuw1


!--------Compute cell face areas
!--------compute cell area of faces perpendicular to first coord axis
         Do k = kl_bnd1,ku_bnd1
         Do j = jl_bnd1,ju_bnd1
         Do i = il_bnd1,iu_bnd1+1

           cell_area1_1 = 1.
           cell_area1_2 = 1.
           cell_area1_3 = 1.

           If (cartesian_pm) Then
           cell_area1_1 =  1.
           If (ndim >= 2)                                              &
            cell_area1_2 =  dy
           If (ndim == 3)                                              &
             cell_area1_3 = dz
           If (l2p5d == 1) cell_area1_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cartesian_pm)

           If (spherical_pm) Then
           cell_area1_1 = cell_face_coord1(i)**2
           If (ndim >= 2)                                              &
             cell_area1_2 = cos( cell_face_coord2(j)   ) -             &
                            cos( cell_face_coord2(j+1) )
           If (ndim == 3)                                              &
             cell_area1_3 = cell_face_coord3(k+k3d) -                  &
                            cell_face_coord3(k)
           If (l2p5d == 1) cell_area1_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (spherical_pm)

           If (polar_pm) Then
           cell_area1_1 =  cell_face_coord1(i)
           If (ndim >= 2)                                              &
             cell_area1_2 =  cell_face_coord2(j+1) -                   &
                             cell_face_coord2(j)
           End If  ! End If (polar_pm)

           If (cylindrical_pm) Then ! perp to r
!----------INT(dz) * r*INT(d theta)
           cell_area1_1 =  cell_face_coord1(i)
           If (ndim >= 2)                                              &
             cell_area1_2 =  cell_face_coord2(j+1) -                   &
                             cell_face_coord2(j)
           If (l2p5d == 1) cell_area1_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           If (ndim == 3) Then
                cell_area1_3 =  cell_face_coord3(k+k3d) -              &
                                cell_face_coord3(k)
           End If

           End If  ! End If (cylindrical_pm)

           cell_area1(i,j,k) = max(abs(cell_area1_1 * cell_area1_2     &
                                                    * cell_area1_3),   &
                                   eps)

         End Do  ! End Do i = il_bnd1,iu_bnd1+1
         End Do  ! End Do j = jl_bnd1,ju_bnd1
         End Do  ! End Do k = kl_bnd1,ku_bnd1

!--------compute cell area of faces perpendicular to second coord axis
         Do k = kl_bnd1,ku_bnd1
         Do j = jl_bnd1,ju_bnd1+k2d
         Do i = il_bnd1,iu_bnd1

           cell_area2_1 = 1.
           cell_area2_2 = 1.
           cell_area2_3 = 1.

           If (cartesian_pm) Then
           cell_area2_1 = dx
           If (ndim == 3)                                              &
             cell_area2_3 = dz
           If (l2p5d == 1) cell_area2_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cartesiam_pm)

           If (spherical_pm) Then
           cell_area2_1 = (cell_face_coord1(i+1)-cell_face_coord1(i))  &
                  *(cell_face_coord1(i)+cell_face_coord1(i+1))*.5
           cell_area2_2 = sin( cell_face_coord2(j) )
           If (ndim == 3)                                              &
             cell_area2_3 = cell_face_coord3(k+k3d) -                  &
                            cell_face_coord3(k)
           If (l2p5d == 1) cell_area2_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (spherical_pm)

           If (polar_pm) Then
           cell_area2_1 =  cell_face_coord1(i+1) -                     &
                           cell_face_coord1(i)
           End If  ! End If (polar_pm)

           If (cylindrical_pm) Then ! perp to z
!----------INT(rdr) * INT(d theta)
           cell_area2_1 = ( cell_face_coord1(i+1)**2 -                 &
                            cell_face_coord1(i)**2 )*.5
           If (ndim == 3)                                              &
             cell_area2_3 =  cell_face_coord3(k+k3d) -                 &
                             cell_face_coord3(k)
           If (l2p5d == 1) cell_area2_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cylindrical_pm)

           cell_area2(i,j,k) = max(abs(cell_area2_1 * cell_area2_2     &
                                                    * cell_area2_3),   &
                                   eps)

         End Do  ! End Do i = il_bnd1,iu_bnd1
         End Do  ! End Do j = jl_bnd1,ju_bnd1+k2d
         End Do  ! End Do k = kl_bnd1,ku_bnd1

!--------compute cell area of faces perpendicular to third coord axis
         Do k = kl_bnd1,ku_bnd1+k3d
         Do j = jl_bnd1,ju_bnd1
         Do i = il_bnd1,iu_bnd1

           cell_area3_1 =  1.
           cell_area3_2 =  1.
           cell_area3_3 =  1.

           If (cartesian_pm) Then
           cell_area3_1 = dx
           If (ndim >= 2)                                              &
             cell_area3_2 = dy
           End If  ! End If (cartesiam_pm)

           If (spherical_pm) Then
           cell_area3_1 = (cell_face_coord1(i+1)-cell_face_coord1(i))  &
                  *(cell_face_coord1(i)+cell_face_coord1(i+1))*.5
           If (ndim >= 2)                                              &
             cell_area3_2 = cell_face_coord2(j+1) -                    &
                            cell_face_coord2(j)
           End If  ! End If (spherical_pm)

           If (polar_pm) Then
           cell_area3_1 = ( cell_face_coord1(i+1)**2 -                 &
                            cell_face_coord1(i)**2 )*.5
           If (ndim >= 2)                                              &
             cell_area3_2 =  cell_face_coord2(j+1) -                   &
                             cell_face_coord2(j)
           End If  ! End If (polar_pm)

           If (cylindrical_pm) Then  ! perp to theta
!----------INT(dr) * INT(dz)
           cell_area3_1 =  cell_face_coord1(i+1) -                     &
                           cell_face_coord1(i)
           If (ndim >= 2)                                              &
             cell_area3_2 =  cell_face_coord2(j+k2d) -                 &
                             cell_face_coord2(j)
           End If  ! End If (cylindrical_pm)

           cell_area3(i,j,k) = max(abs(cell_area3_1 * cell_area3_2     &
                                                    * cell_area3_3),   &
                                   eps)

         End Do  ! End Do i = il_bnd1,iu_bnd1
         End Do  ! End Do j = jl_bnd1,ju_bnd1
         End Do  ! End Do k = kl_bnd1,ku_bnd1+k3d


!--------Compute cell edge lengths
!--------compute edge length in direction of first coord axis
         Do k = kl_bnd1,ku_bnd1+k3d
         Do j = jl_bnd1,ju_bnd1+k2d
         Do i = il_bnd1,iu_bnd1

           cell_leng1_1 =  1.
           cell_leng1_2 =  1.
           cell_leng1_3 =  1.

           If (cartesian_pm) Then
           cell_leng1_1 = dx
           Else
           cell_leng1_1 =  cell_face_coord1(i+1) -                     &
                           cell_face_coord1(i)
           End If  ! End If (cartesian_pm)

           cell_leng1(i,j,k) = max(cell_leng1_1, eps)

         End Do  ! End Do i = il_bnd1,iu_bnd1
         End Do  ! End Do j = jl_bnd1,ju_bnd1+k2d
         End Do  ! End Do k = kl_bnd1,ku_bnd1+k3d

!--------compute edge length in direction of second coord axis
         Do k = kl_bnd1,ku_bnd1+k3d
         Do j = jl_bnd1,ju_bnd1
         Do i = il_bnd1,iu_bnd1+1

           cell_leng2_1 =  1.
           cell_leng2_2 =  1.
           cell_leng2_3 =  1.

           If (cartesian_pm) Then
           If (ndim >= 2)                                              &
             cell_leng2_2 = dy
           End If  ! End If (cartesiam_pm)

           If (spherical_pm) Then
           cell_leng2_1 =  cell_face_coord1(i)
           If (ndim >= 2)                                              &
             cell_leng2_2 =  cell_face_coord2(j+1) -                   &
                             cell_face_coord2(j)
           End If  ! End If (spherical_pm)

           If (polar_pm) Then
           cell_leng2_1 =  cell_face_coord1(i)
           If (ndim >= 2)                                              &
             cell_leng2_2 =  cell_face_coord2(j+1) -                   &
                             cell_face_coord2(j)
           End If  ! End If (polar_pm)

           If (cylindrical_pm) Then
           If (ndim == 2)                                              &
             cell_leng2_2 =  cell_face_coord2(j+k3d) -                 &
                             cell_face_coord2(j)
           End If  ! End If (cylindrical_pm)

           cell_leng2(i,j,k) = max(cell_leng2_1 * cell_leng2_2         &
                                                * cell_leng2_3,        &
                                   eps)

         End Do  ! End Do i = il_bnd1,iu_bnd1+1
         End Do  ! End Do j = jl_bnd1,ju_bnd1
         End Do  ! End Do k = kl_bnd1,ku_bnd1+k3d

!--------compute edge length in direction of third coord axis
         Do k = kl_bnd1,ku_bnd1
         Do j = jl_bnd1,ju_bnd1+k2d
         Do i = il_bnd1,iu_bnd1+1

           cell_leng3_1 =  1.
           cell_leng3_2 =  1.
           cell_leng3_3 =  1.

           If (cartesian_pm) Then
           If (ndim == 3)                                              &
             cell_leng3_3 = dz
           If (l2p5d == 1) cell_leng3_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cartesian_pm)

           If (spherical_pm) Then
           cell_leng3_1 =  cell_face_coord1(i)
           cell_leng3_2 =  sin( cell_face_coord2(j) )
           If (ndim == 3)                                              &
             cell_leng3_3 =  cell_face_coord3(k+k3d) -                 &
                             cell_face_coord3(k)
           If (l2p5d == 1) cell_leng3_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (spherical_pm)

           If (cylindrical_pm) Then
           cell_leng3_1 =  cell_face_coord1(i)
           If (ndim == 3)                                              &
             cell_leng3_3 =  cell_face_coord3(k+k3d) -                 &
                             cell_face_coord3(k)
           If (l2p5d == 1) cell_leng3_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           End If  ! End If (cylindrical_pm)

           cell_leng3(i,j,k) = max(cell_leng3_1 * cell_leng3_2         &
                                                * cell_leng3_3,        &
                                   eps)

         End Do  ! End Do i = il_bnd1,iu_bnd1+1
         End Do  ! End Do j = jl_bnd1,ju_bnd1+k2d
         End Do  ! End Do k = kl_bnd1,ku_bnd1

      Return
      End Subroutine amr_block_geometry


