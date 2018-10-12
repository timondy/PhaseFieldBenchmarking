!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/clean_field
!! NAME
!!
!!   clean_field
!!
!! SYNOPSIS
!!
!!   Call clean_field (facex, facey, facez,                  
!!                     iminx,imaxx,jminx,jmaxx,kminx,kmaxx,  
!!                     iminy,imaxy,jminy,jmaxy,kminy,kmaxy,  
!!                     iminz,imaxz,jminz,jmaxz,kminz,kmaxz,  
!!                     ndim, nxb, nyb, nzb, bsize)
!!   Call clean_field (real array, real array, real array,                  
!!                     integer,integer,integer,integer,integer,integer,  
!!                     integer,integer,integer,integer,integer,integer,  
!!                     integer,integer,integer,integer,integer,integer,  
!!                     integer,integer,integer,integer,real array)
!!
!! ARGUMENTS
!!
!!   The following are dimensions for the arrays face, facey, facez
!!    Integer, Intent(in) :: iminx, imaxx, jminx, jmaxx, kminx, kmaxx
!!    Integer, Intent(in) :: iminy, imaxy, jminy, jmaxy, kminy, kmaxy
!!    Integer, Intent(in) :: iminz, imaxz, jminz, jmaxz, kminz, kmaxz
!!
!!    Integer, Intent(in) :: ndim, nxb, nyb, nzb  ndim is the number of dimensions
!!                                                nxb, nyb, nzb are the number of cells
!!                                                in the x, y, and z directions for
!!                                                a blocl
!!   facx, facey, and facez hold the face variable data for the field which is to be
!!   cleaned
!!    Real, Intent(inout) :: facex(iminx:imaxx,jminx:jmaxx,kminx:kmaxx)
!!    Real, Intent(inout) :: facey(iminy:imaxy,jminy:jmaxy,kminy:kmaxy)
!!    Real, Intent(inout) :: facez(iminz:imaxz,jminz:jmaxz,kminz:kmaxz)
!!   bsize is the the physical size of the block in each coordinate direction
!!    Real, Intent(in)    :: bsize(3)
!!
!! INCLUDES
!!
!! USES
!!
!! CALLS
!!
!!   poisson_sor
!!
!! RETURNS
!!
!!   Returns the corrected (cleaned) fields in 'facex', 'facey', and 'facez'.
!!
!! DESCRIPTION
!!
!!   This routine 'cleans' a vector field stored at the faces of cells.
!!   This routine implements the classis 'Brackbill-Barnes' correction algorithm
!!   where the divergence of the field is first calculated and is used as the
!!   source term in poisson's equation.  This poisson's equation is then solved
!!   by calling 'poisson_sor' which returns corrections.  Finally, the corrections
!!   are added to the orinally input fields.
!!
!! AUTHOR
!!
!!   Kevin Olson, May 2007
!!
!!***
    

      Subroutine clean_field (facex, facey, facez,                  &
                              iminx,imaxx,jminx,jmaxx,kminx,kmaxx,  &
                              iminy,imaxy,jminy,jmaxy,kminy,kmaxy,  &
                              iminz,imaxz,jminz,jmaxz,kminz,kmaxz,  &
                              ndim, nxb, nyb, nzb, bsize)

      Implicit None

!-----Input/Output arguments.
      Integer, Intent(in) :: iminx, imaxx, jminx, jmaxx, kminx, kmaxx
      Integer, Intent(in) :: iminy, imaxy, jminy, jmaxy, kminy, kmaxy
      Integer, Intent(in) :: iminz, imaxz, jminz, jmaxz, kminz, kmaxz
      Integer, Intent(in) :: ndim, nxb, nyb, nzb

      Real, Intent(inout) :: facex(iminx:imaxx,jminx:jmaxx,kminx:kmaxx)
      Real, Intent(inout) :: facey(iminy:imaxy,jminy:jmaxy,kminy:kmaxy)
      Real, Intent(inout) :: facez(iminz:imaxz,jminz:jmaxz,kminz:kmaxz)
      Real, Intent(in)    :: bsize(3)
    
   
!-----Local Variables
      Real, allocatable :: rhs(:,:,:)
      Real, allocatable :: lhs(:,:,:)
      Real :: dx, dy, dz
      Real :: divb, divb_max
      Real :: gradx, grady, gradz, laplac
      Real, Parameter :: eps = 1.e-10

      Integer :: i, j, k, k3d, k2d, nit, maxit
      Integer :: imin, imax, jmin, jmax, kmin, kmax
     
!-----Begin Executable code

      dx = bsize(1)/Real(nxb)
      dy = bsize(2)/Real(nyb)
      dz = bsize(3)/Real(nzb)
      k2d = 0
      k3d = 0
      If (ndim >= 2) k2d = 1
      If (ndim == 3) k3d = 1

      imin = min(iminx, iminy, iminz)
      jmin = min(jminx, jminy, jminz)
      kmin = min(kminx, kminy, kminz)

      imax = min(imaxx, imaxy, imaxz)
      jmax = min(jmaxx, jmaxy, jmaxz)
      kmax = min(kmaxx, kmaxy, kmaxz)

      Allocate(rhs(imin-1:imax+1,jmin-k2d:jmax+k2d,kmin-k3d:kmax+k3d))
      Allocate(lhs(imin-1:imax+1,jmin-k2d:jmax+k2d,kmin-k3d:kmax+k3d))

!-----compute divb and set the source terms

      divb_max = 1.e10
      rhs(:,:,:) = 0.

      do k = kmin, kmax
         Do j = jmin, jmax
            Do i = imin, imax

               divb = (facex(i+1,j,k) - facex(i,j,k))/dx
               If (ndim >= 2) then
                  divb = divb + (facey(i,j+1,k) - facey(i,j,k))/dy
               End If
               If (ndim == 3) then
                  divb = divb + (facez(i,j,k+1) - facez(i,j,k))/dz
               End If
            
               rhs(i,j,k) = divb

            End Do
         End Do
      End Do

      lhs(:,:,:) = 0.
!-----Solve the Poisson equation.
      Call poisson_sor(lhs,rhs,imax-imin+1,jmax-jmin+1,kmax-kmin+1,   &
                       imax-imin+3,                                   &
                       jmax-jmin+1+(2*k2d),                           &
                       kmax-kmin+1+(2*k3d),                           &
                       dx,dy,dz,1,ndim)

!-----Now correct the fields.
!-----compute the gradient of the correcting potential (lhs)
!-----and subtract from the divbs at the cell faces
      Do k = kmin, kmax
         Do j = jmin, jmax
            Do i = imin, imax+1
               gradx = (lhs(i,j,k) - lhs(i-1,j,k))/dx
               facex(i,j,k) = facex(i,j,k) - gradx
            End Do
         End Do
      End Do
      If (ndim >= 2) then
      Do k = kmin, kmax
         Do j = jmin, jmax+k2d
            Do i = imin, imax
               grady = (lhs(i,j,k) - lhs(i,j-1,k))/dy
               facey(i,j,k) = facey(i,j,k) - grady
            End Do
         End Do
      End Do
      End If
      If (ndim == 3) then
      Do k = kmin, kmax+k3d
         Do j = jmin, jmax
            Do i = imin, imax
               gradz = (lhs(i,j,k) - lhs(i,j,k-1))/dz
               facez(i,j,k) = facez(i,j,k) - gradz
            End Do
         End Do
      End Do
      End If

! Test by computing divb again and outputing an error statement if it is above a 
! tolerance

!      divb_max = -1.e30
!      Do k = kmin, kmax
!         Do j = jmin, jmax
!            Do i = imin, imax

!               divb = (facex(i+1,j,k) - facex(i,j,k))/dx
!               If (ndim >= 2) then
!                  divb = divb + (facey(i,j+1,k) - facey(i,j,k))/dy
!               End If
!               If (ndim == 3) then
!                  divb = divb + (facez(i,j,k+1) - facez(i,j,k))/dz
!               End If

!               divb_max = max(abs(divb),divb_max)
            
!            End Do
!         End Do
!      End Do<

      Deallocate(rhs)
      Deallocate(lhs)

      Return
      End Subroutine clean_field
