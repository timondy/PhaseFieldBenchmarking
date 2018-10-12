!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_work_fun_recip
!! NAME
!!
!!   amr_restrict_work_fun_recip
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_work_fun_recip(datainw,dataoutw)
!!   Call amr_restrict_work_fun_recip(real array,real array)
!!
!! ARGUMENTS
!!
!!   Real, Intent(in)    :: datainw(:,:,:)   data to restrict
!!   Real, Intent(inout) :: dataoutw(:,:,:)  restricted, returned data
!! 
!! INCLUDES
!! 
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!
!! CALLS
!!
!! DESCRIPTION  
!!
!!   This routine performs restriction on the array datainw and
!!   returns the result in dataoutw. Note that this does not update
!!   guard cell elements of dataoutw.
!!
!!   This particular version applies the 3D generalization of the
!!   restriction operator in eqn (19.6.17) of the 2nd edition of
!!   Numerical recipes.
!!   The 2D case is
!!                  | 1  2  1 |
!!                  | 2  4  2 | /16.
!!                  | 1  2  1 |
!!
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 1997
!!
!!***

      Subroutine amr_restrict_work_fun_recip(datainw,dataoutw)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use workspace

      Implicit None

!-----Input/Output arguments.
      Real, Intent(in)    :: datainw(:,:,:)
      Real, Intent(inout) :: dataoutw(:,:,:)

!-----local arrays and variables
      integer :: i,j,k

!-----Begin executable code

      Do k=1+nguard_work*k3d,nzb+nguard_work*k3d
      Do j=1+nguard_work*k2d,nyb+nguard_work*k2d
      Do i=1+nguard_work,nxb+nguard_work
       dataoutw(i,j,k) = (                                             & 
         ( datainw(i-1,j-k2d,k-k3d) + 2.*datainw(i,j-k2d,k-k3d) +      & 
                    datainw(i+1,j-k2d,k-k3d)   ) +                     & 
        2.*(   datainw(i-1,j,k-k3d) + 2.*datainw(i,j,k-k3d) +          & 
                    datainw(i+1,j,k-k3d)   ) +                         & 
         (   datainw(i-1,j+k2d,k-k3d) + 2.*datainw(i,j+k2d,k-k3d) +    & 
                    datainw(i+1,j+k2d,k-k3d)   ) +                     & 
        2.*(   datainw(i-1,j-k2d,k) + 2.*datainw(i,j-k2d,k) +          & 
                    datainw(i+1,j-k2d,k)   ) +                         & 
        4.*(   datainw(i-1,j,k) +  2.*datainw(i,j,k) +                 & 
                    datainw(i+1,j,k)   ) +                             & 
        2.*(   datainw(i-1,j+k2d,k) + 2.*datainw(i,j+k2d,k) +          & 
                    datainw(i+1,j+k2d,k)   ) +                         & 
         (   datainw(i-1,j-k2d,k+k3d) + 2.*datainw(i,j-k2d,k+k3d) +    & 
                    datainw(i+1,j-k2d,k+k3d)   ) +                     & 
        2.*(   datainw(i-1,j,k+k3d) + 2.*datainw(i,j,k+k3d) +          & 
                    datainw(i+1,j,k+k3d)   ) +                         & 
         (   datainw(i-1,j+k2d,k+k3d) + 2.*datainw(i,j+k2d,k+k3d) +    & 
                    datainw(i+1,j+k2d,k+k3d)   )    )/64.

      End Do
      End Do
      End Do

      Return
      End Subroutine amr_restrict_work_fun_recip
