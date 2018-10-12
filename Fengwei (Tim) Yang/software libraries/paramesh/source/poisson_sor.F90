!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/poisson_sor
!! NAME
!!   poisson_sor
!!
!! SYNOPSIS
!!
!!   Call poisson_sor(lhs,rhs,  
!!                    nx,ny,nz, 
!!                    iu,ju,ku, 
!!                    dx,dy,dz, 
!!                    nguard,ndim)
!!   Call poisson_sor(Real array,Real array,  
!!                    Integer,Integer,Integer, 
!!                    Integer,Integer,Integer, 
!!                    Real,Real,Real, 
!!                    Integer,Integer)
!!
!! ARGUMENTS
!!
!!      Real,Intent(inout)  :: rhs(iu,ju,ku)  source terms on the mesh
!!      Real,Intent(out)    :: lhs(iu,ju,ku)  potential on the mesh
!!      Integer, Intent(in) :: nx,ny,nz no. of cells in each direction wo. guardcells
!!      Integer, Intent(in) :: iu,ju,ku no. of cells in each direction including guardcells
!!      Real,Intent(in)     :: dx, dy, dz sizes of cells in each direction
!!      Integer, Intent(in) :: nguard,ndim no. of guardcells and dimensions
!!
!! INCLUDES
!!
!! USES
!!
!! CALLS
!!
!! RETURNS
!!
!!   Potential returned in the array 'lhs'.
!!
!! DESCRIPTION
!!
!!   This file contains routines for solving the poisson equation on a uniform mesh.
!!   The solver used is a simple Gauss-Seidel scheme with red-black
!!   ordering. The data is assumed to be at cell centers.  This solver has been written
!!   specifically to work on a single APRAMESH block and is designed to be used for
!!   divergence cleaning on a single block as part to the prolongation operation.
!!   Therefore the boundary condition for the potential are implicitly assumed tobe
!!   Neumann (i.e. grad(phi) = 0. at the boundarys)
!!
!! AUTHOR
!!
!!   Kevin Olson Feb. 2007.
!!
!!***

      Subroutine poisson_sor(lhs,rhs,                                  &
                             nx,ny,nz,                                 &
                             iu,ju,ku,                                 &
                             dx,dy,dz,                                 &
                             nguard,ndim)

      Implicit None

!-----Input/Output arguments.
      Integer, Intent(in) :: nx,ny,nz,iu,ju,ku,nguard,ndim
      Real,Intent(inout)  :: rhs(iu,ju,ku)
      Real,Intent(out)    :: lhs(iu,ju,ku)
      Real,Intent(in)     :: dx, dy, dz

!-----Local arrays and variables.
      Real :: fxf(iu,ju,ku), fyf(iu,ju,ku), fzf(iu,ju,ku), res
      Real :: rhs2(iu,ju,ku)
      Real :: local_norm, local_norm_rhs
      Real :: fac, laplac, dz2i, dy2i, dx2i
      Real :: tot_src, tot_vol, x, y, z
      Integer :: kkk, i, j, k, is, ie, nit, k2d, k3d
      Real :: eps
      Parameter (eps=1.e-20)

!-----Begin executable code.
      k2d = 0
      k3d = 0
      If (ndim >= 2) k2d = 1
      If (ndim == 3) k3d = 1

!-----compute the total volume and total source term
      tot_vol = (dx*dy*dz)*(nx*ny*nz) ! This needs to be modified for curvilinear
      tot_src = 0.
      Do k = nguard*k3d+1,nguard*k3d+nz
         Do j = nguard*k2d+1,nguard*k2d+ny
            Do i = nguard+1,nguard+nx
               tot_src = tot_src + rhs(i,j,k)*(dx*dy*dz)
            End Do
         End Do
      End Do

      rhs2(:,:,:) = rhs(:,:,:) - tot_src/tot_vol

      local_norm = 1.e10
      local_norm_rhs = 1.
      kkk = 0
      Do while((local_norm/local_norm_rhs > eps .or. kkk <= 2)  &
                .and. kkk < 10000)

      dz2i = 1./dz**2
      dy2i = 1./dy**2
      dx2i = 1./dx**2

      fac = dx2i
      If (ndim >= 2) Then
         fac = fac + dy2i
      End If
      If (ndim == 3) Then
         fac = fac + dz2i
      End If
      fac = -2.*fac
      fac = 1./fac

      local_norm = 0.
      Do nit = 0, 1  ! loop twice to cover entire mesh w. red-black ordering

      Call set_zero_grad_bc (lhs,rhs,nx,ny,nz,iu,ju,ku,nguard,ndim,    &
                             k2d,k3d)

      Call fluxes(lhs,fxf,fyf,fzf,dx,dy,dz,nx,ny,nz,iu,ju,ku,nguard,   &
                  k2d,k3d)

      Do k = 1+nguard*k3d,nz+nguard*k3d
        Do j = 1+nguard*k2d,ny+nguard*k2d
           If (mod(nit,2) == 0) Then
           If (mod(j+k,2) == 0) Then
              is = 1+nguard
              ie = nx+nguard-1
           Else
              is = 1+nguard+1
              ie = nx+nguard
           End If
           Else
           If (mod(j+k,2).eq.0) Then
              is = 1+nguard+1
              ie = nx+nguard
           Else
              is = 1+nguard
              ie = nx+nguard-1
           End If
           End If
           Do i = is,ie,2

! compute laplacian -> laplac

             laplac = (fxf(i+1,j,k)-fxf(i,j,k))/dx
             If (ndim >= 2) Then
                laplac = laplac + (fyf(i,j+1,k)-fyf(i,j,k))/dy
             End If
             If (ndim == 3) Then
                laplac = laplac + (fzf(i,j,k+1)-fzf(i,j,k))/dz
             End If

! compute new value of lhs

            res = laplac - rhs2(i,j,k)
            lhs(i,j,k) = lhs(i,j,k) - fac*res
            local_norm = local_norm + res**2

          End Do  ! End Do i = is,ie,2
        End Do  ! End Do j = 1+nguard*k2d,ny+nguard*k2d
      End Do  ! End Do k = 1+nguard*k3d,nz+nguard*k3d

      End Do  ! End Do nit = 0, 1

      kkk = kkk + 1

      End Do  ! End Do while((local_norm/local_norm_rhs > eps .or. kkk <= 2)  
              !               .and. kkk < 10000)

      If (kkk >= 10000) Then
      Print *,' WARNING: poisson solver did not converge during face ' 
      Print *,' variable divergence cleaning.  '
      Print *,' Suggest your check your results '
      End If

      Call set_zero_grad_bc (lhs,rhs,nx,ny,nz,iu,ju,ku,nguard,ndim,    &
                             k2d,k3d)

!-----remove contribution from constant source term.
      fac = (.5/Real(ndim))*tot_src/tot_vol
      Do k = 1,2*nguard*k3d + nz
      If (ndim == 3) Then
         z = (k-(nguard*k3d+1)+.5)*dz
      Else 
         z = 0.
      End If
      Do j = 1,2*nguard*k2d + ny
      If (ndim >= 2) Then
         y = (j-(nguard*k2d+1)+.5)*dy
      Else
         y = 0.
      End If
      Do i = 1,2*nguard + nx
          x = (i-(nguard+1)+.5)*dx
          lhs(i,j,k) = lhs(i,j,k) + fac * (x**2 +                      &
                                           y**2 +                      &
                                           z**2)
      End Do
      End Do
      End Do

      Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine set_zero_grad_bc (lhs,rhs,nx,ny,nz,iu,ju,ku,nguard,   &
                                   ndim,k2d,k3d)

      Implicit none

!-----Input/Output arguments.
      Integer, Intent(in) :: nx, ny, nz, iu, ju, ku, nguard,           &
                             ndim, k2d, k3d
      Real, Intent(inout) :: lhs(iu,ju,ku)
      Real, Intent(in)    :: rhs(iu,ju,ku)
!-----Local variables.
      Integer :: i, j, k, ic, jc, kc

!-----Being executable code.

!-----x faces of block
      Do k = nguard*k3d+1,nguard*k3d+nz
         Do j = nguard*k2d+1,nguard*k2d+ny
            ! -x face
            Do i = 1,nguard
               ic = 2*nguard - (i-1)
               lhs(i,j,k) = lhs(ic,j,k)
            End Do
            ! +x face
            Do i = nguard+1+nx,2*nguard+nx
               ic = nx+nguard - (i-(nx+nguard+1))
               lhs(i,j,k) = lhs(ic,j,k)
            End Do
         End Do
      End Do

      If (ndim >= 2) Then
!-----y faces of block
      Do k = nguard*k3d+1,nguard*k3d+nz
         Do i = nguard+1,nguard+nx
            ! -y face
            Do j = 1,nguard
               jc = 2*nguard - (j-1)
               lhs(i,j,k) = lhs(i,jc,k)
            End Do
            ! +y face
            Do j = nguard+1+ny,2*nguard+ny
               jc = ny+nguard - (j-(ny+nguard+1))
               lhs(i,j,k) = lhs(i,jc,k)
            End Do
         End Do
      End Do
      End If

      If (ndim == 3) Then
!-----z faces of block
      Do j = nguard+1,nguard+ny
         Do i = nguard+1,nguard+nx
            ! -z face
            Do k = 1,nguard
               kc = 2*nguard - (k-1)
               lhs(i,j,k) = lhs(i,j,kc)
!               lhs(i,j,k) = 0.
            End Do
            ! +z face
            Do k = nguard+1+nz,2*nguard+nz
               kc = nz+nguard - (k-(nz+nguard+1))
               lhs(i,j,k) = lhs(i,j,kc)
            End Do
         End Do
      End Do
      End If

      Return
      End Subroutine set_zero_grad_bc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine fluxes(lhs,fxf,fyf,fzf,dx,dy,dz,nx,ny,nz,             &
                        iu,ju,ku,nguard,k2d,k3d)

      Implicit None

!-----Input/Output arguments
      Integer, Intent(in) :: nx, ny, nz, iu, ju, ku, nguard,           &
                             k2d, k3d
      Real, Intent(in) :: lhs(iu,ju,ku), dx, dy, dz
      Real, Intent(out) :: fxf(iu,ju,ku),fyf(iu,ju,ku),fzf(iu,ju,ku)
!-----Local variables
      Integer :: i, j, k

!-----Begin executable code.

!-----compute first derivatives at cell faces

      Do k = nguard*k3d+1, nguard*k3d+nz
      Do j = nguard*k2d+1, nguard*k2d+ny
      Do i = nguard+1, nguard+nx+1

         fxf(i,j,k) = (lhs(i,j,k)-lhs(i-1,j,k))/dx

      End Do
      End Do
      End Do

      Do k = nguard*k3d+1, nguard*k3d+nz
      Do j = nguard*k2d+1, nguard*k2d+ny+k2d
      Do i = nguard+1, nguard+nx

         fyf(i,j,k) = (lhs(i,j,k)-lhs(i,j-k2d,k))/dy

      End Do
      End Do
      End Do


      Do k = nguard*k3d+1, nguard*k3d+nz+k3d
      Do j = nguard*k2d+1, nguard*k2d+ny
      Do i = nguard+1, nguard+nx

         fzf(i,j,k) = (lhs(i,j,k)-lhs(i,j,k-k3d))/dz

      End Do
      End Do
      End Do

      Return
      End Subroutine fluxes

      End Subroutine poisson_sor


