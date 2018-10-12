!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_restrict_edge
!! NAME
!!
!!   amr_restrict_edge
!! 
!! SYNOPSIS
!!
!!   Call amr_restrict_edge(icoord)
!!   Call amr_restrict_edge(integer)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in) :: icoord Indicates which edge is to be restricted.
!!                                 icoord = 1 -> x edge
!!                                 icoord = 2 -> y edge
!!                                 icoord = 3 -> z edge
!!
!! INCLUDES
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! CALLS
!!
!! RETURNS
!!
!!   Nothing.
!!
!! DESCRIPTION
!!
!!   This routine performs a restriction operation on the arrays 
!!   recvarx(y)(z)[1,2] and returns the result in the same arrays.
!!   These data arrays are defined on block boundaries only.
!!   The restriction operation coded here is for edge data and operates
!!   as part of the edge averaging proceedure for algorithms which need
!!   to enforce edge data consistency at jumps in refinement.
!!
!!   Note that this does not update guard cell elements of recvarx(y)(z)[1,2].
!!
!!   Also note that we use stride 2 along each dimension when computing
!!   reduced data values on block faces, so not all values of the averaged data
!!   have been updated.
!!
!!   This particular version is only appropriate for 2nd order schemes 
!!   using linear interpolation with even number of mesh points along 
!!   each block axis.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          July 1997
!!
!!***

      Subroutine amr_restrict_edge(icoord)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata

      Implicit None

!-----Input/Output arguments.
      Integer, Intent(in) :: icoord

!-----Local variables and arrays.
      Real :: fact
      Integer ::  nguard0
      Integer :: i,j,k

!-----Begin executable code.

      nguard0 = nguard*npgs

      fact = .5
      If (edge_value_integ) Then
         fact = 1.
      End If

      If (icoord == 1) Then                              ! edges on x-face

        Do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
          Do j=1+nguard0,nyb+nguard0,2
            Do i=1,2
!-------------y pointing edge first
              recvarx1e(:,i,j,k) = (                                   & 
                        recvarx1e(:,i,j,k) +                           & 
                        recvarx1e(:,i,j+1,k) )*fact
            End Do
          End Do
        End Do

        If ((ndim == 3).or.(l2p5d == 1)) Then
          Do k=1+nguard0*k3d,nzb+nguard0*k3d,2
            Do j=1+nguard0,nyb+nguard0+1,2
              Do i=1,2
!---------------z pointing edge 
                recvarx2e(:,i,j,k) = (                                 & 
                          recvarx2e(:,i,j,k) +                         & 
                          recvarx2e(:,i,j,k+k3d) )*fact
              End Do
            End Do
          End Do
        End If  ! End If ((ndim == 3).or.(l2p5d == 1))

      ElseIf (icoord == 2) Then                          ! edges on y-face

        If ((ndim == 3).or.(l2p5d == 1)) Then
          Do k=1+nguard0*k3d,nzb+nguard0*k3d,2
            Do j=1,2
              Do i=1+nguard0,nxb+nguard0+1,2
!---------------z pointing edge first
                recvary2e(:,i,j,k) = (                                 & 
                          recvary2e(:,i,j,k) +                         & 
                          recvary2e(:,i,j,k+k3d) )*fact
              End Do
            End Do
          End Do
        End If  ! End If ((ndim == 3).or.(l2p5d == 1))

        Do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
          Do j=1,2
            Do i=1+nguard0,nxb+nguard0,2
!-------------x pointing edge
              recvary1e(:,i,j,k) = (                                   & 
                        recvary1e(:,i,j,k) +                           & 
                        recvary1e(:,i+1,j,k) )*fact
            End Do
          End Do
        End Do

      ElseIf (icoord == 3) Then                          ! edges on z-face

        Do k=1,2
          Do j=1+nguard0,nyb+nguard0+1,2
            Do i=1+nguard0,nxb+nguard0,2
!-------------x pointing edge first
              recvarz1e(:,i,j,k) = (                                   & 
                        recvarz1e(:,i,j,k) +                           & 
                        recvarz1e(:,i+1,j,k) )*fact
            End Do
          End Do
        End Do

        Do k=1,2
          Do j=1+nguard0,nyb+nguard0,2
            Do i=1+nguard0,nxb+nguard0+1,2
!-------------y pointing edge
              recvarz2e(:,i,j,k) = (                                   & 
                        recvarz2e(:,i,j,k) +                           & 
                        recvarz2e(:,i,j+1,k) )*fact
            End Do
          End Do
        End Do

      End If  ! End If (icoord == 1)

      Return
      End Subroutine amr_restrict_edge
