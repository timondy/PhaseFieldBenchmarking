!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/send_block_data
!! NAME
!!
!!   send_block_data
!! 
!! SYNOPSIS
!!
!!   call send_block_data (lb, new_loc, old_loc, free, moved, sent, 
!!                         lnblocks_old, mype, nmoved, test, point_to,
!!                         reqs, nsend, unk_int_type,
!!                         facex_int_type, facey_int_type,
!!                         facez_int_type, edgex_int_type,
!!                         edgey_int_type, edgez_int_type,
!!                         unkn_int_type)
!!
!!   call send_block_data (integer, integer, integer, logical, logical, logical,
!!                         integer, integer, integer, integer, integer,
!!                         integer, integer, integer,
!!                         integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer :: lb
!!     Block being sent
!!
!!   integer :: new_loc(2,maxblocks_tr)
!!     Array which stores the new locations in memory where blocks are to be
!!     moved.  new_loc(1,:) indicates the on-processor location where the block
!!     is to reside and new_loc(2,:) indicates which processor the block data is
!!     to be moved to.
!!
!!   integer :: old_loc(2,maxblocks_tr)
!!     Similar to new_loc above except stores the location in memory where a 
!!     block orginates from.  Necessary for hand shaking during message passing.
!!
!!   logical :: free(maxblocks)
!!     A logical placeholder which indicates if a position in the list of blocks
!!     is free to receive data.
!!
!!   logical :: moved(maxblocks)
!!     A logical which indicates if data has moved from a location in the list of 
!!     blocks.
!!
!!   logical :: sent(maxblocks)
!!     A logical which indicates if data has been sent from a location in the list 
!!     of blocks.
!!
!!   integer :: lnblocks_old
!!     The 'old' number of blocks in the calling processor before data redistribution.
!!
!!   integer :: mype     
!!     Current processor number
!!
!!   integer :: nmoved
!!     Indicates the number of blocks which have been moved during the algorithm. 
!!
!!   integer :: test(maxblocks)
!!     ???
!!
!!   integer :: point_to(maxblocks)
!!     ???
!!   
!!   integer :: reqs(maxblocks_tr)
!!     List of MPI requests.
!!
!!   integer :: nsend
!!     Number MPI sends issued.
!!
!!   integer :: unk_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (unk).
!!
!!   integer :: facex_int_type, facey_int_type, facez_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (face variables)
!!
!!   integer :: edgex_int_type, edgey_int_type, edgez_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (edge variables)
!!
!!   integer :: unkn_int_type
!!     An MPI derived type defining the interior region of block not
!!     including guardcells (unk_n).
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!!
!!
!! RETURNS
!!
!!
!! DESCRIPTION
!!
!!   This routine is part of the redistribution of block data and manages the 
!!   sending of block data to their new locations.  The block data is sent to
!!   its final new location (new_loc) if that location is unoccupied by another
!!   block's data, it is sent to another free location, or nothing is sent.
!!
!!   The routine is only called by the subroutine amr_redist_blk and should never
!!   need to be called by a user's application.
!!
!! AUTHORS
!!
!!   Kevin Olson (1998)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f, unk_test, face[xyz]_test, edge[xyz]_test, unkn_test
#include "paramesh_preprocessor.fh"

      Subroutine send_block_data (lb, new_loc, old_loc, free,          & 
                                  moved, sent,                         & 
                                  lnblocks_old, mype, nmoved,          & 
                                  test, point_to,                      & 
                                  reqs, nsend, unk_int_type,           &
                                  facex_int_type, facey_int_type,      &
                                  facez_int_type, edgex_int_type,      &
                                  edgey_int_type, edgez_int_type,      &
                                  unkn_int_type)
      
!-----Use statements   
      Use paramesh_dimensions
      Use physicaldata
      Use tree

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output statements.
      Integer :: new_loc(2,maxblocks_tr), old_loc(2,maxblocks_tr)
      Logical :: free(maxblocks), moved(maxblocks), sent(maxblocks)
      Integer :: lb, lnblocks_old, mype
      Integer :: reqs(maxblocks_tr), nsend
      Integer :: nmoved
      Integer :: point_to(maxblocks),test(maxblocks)
      Integer :: unk_int_type
      Integer :: facex_int_type, facey_int_type, facez_int_type
      Integer :: edgex_int_type, edgey_int_type, edgez_int_type
      Integer :: unkn_int_type

!-----Local arrays and variables.
      Integer,Save ::  is_unk,js_unk,ks_unk,ie_unk,je_unk,ke_unk
      Integer,Save ::  is_facex,js_facex,ks_facex,ie_facex,je_facex,ke_facex
      Integer,Save ::  is_facey,js_facey,ks_facey,ie_facey,je_facey,ke_facey
      Integer,Save ::  is_facez,js_facez,ks_facez,ie_facez,je_facez,ke_facez
      Integer,Save ::  is_edgex,js_edgex,ks_edgex,ie_edgex,je_edgex,ke_edgex
      Integer,Save ::  is_edgey,js_edgey,ks_edgey,ie_edgey,je_edgey,ke_edgey
      Integer,Save ::  is_edgez,js_edgez,ks_edgez,ie_edgez,je_edgez,ke_edgez
      Integer,Save ::  is_unkn,js_unkn,ks_unkn,ie_unkn,je_unkn,ke_unkn
      Integer :: status(MPI_STATUS_SIZE)
      Integer :: lb2, ierr
      Logical, Save :: first = .True.
      Logical :: success

!-----Begin executable code.

      If (first) Then

      first = .False.

      is_unk = nguard*npgs+1
      js_unk = nguard*k2d*npgs+1
      ks_unk = nguard*k3d*npgs+1
      ie_unk = nguard*npgs+nxb
      je_unk = nguard*k2d*npgs+nyb
      ke_unk = nguard*k3d*npgs+nzb

      is_facex = nguard*npgs+1
      js_facex = nguard*k2d*npgs+1
      ks_facex = nguard*k3d*npgs+1
      ie_facex = nguard*npgs+nxb + 1
      je_facex = nguard*k2d*npgs+nyb
      ke_facex = nguard*k3d*npgs+nzb

      is_facey = nguard*npgs+1
      js_facey = nguard*k2d*npgs+1
      ks_facey = nguard*k3d*npgs+1
      ie_facey = nguard*npgs+nxb 
      je_facey = nguard*k2d*npgs+nyb + k2d
      ke_facey = nguard*k3d*npgs+nzb

      is_facez = nguard*npgs+1
      js_facez = nguard*k2d*npgs+1
      ks_facez = nguard*k3d*npgs+1
      ie_facez = nguard*npgs+nxb 
      je_facez = nguard*k2d*npgs+nyb
      ke_facez = nguard*k3d*npgs+nzb + k3d

      is_edgex = nguard*npgs+1
      js_edgex = nguard*k2d*npgs+1
      ks_edgex = nguard*k3d*npgs+1
      ie_edgex = nguard*npgs+nxb
      je_edgex = nguard*k2d*npgs+nyb + k2d
      ke_edgex = nguard*k3d*npgs+nzb + k3d

      is_edgey = nguard*npgs+1
      js_edgey = nguard*k2d*npgs+1
      ks_edgey = nguard*k3d*npgs+1
      ie_edgey = nguard*npgs+nxb + 1
      je_edgey = nguard*k2d*npgs+nyb
      ke_edgey = nguard*k3d*npgs+nzb + k3d

      is_edgez = nguard*npgs+1
      js_edgez = nguard*k2d*npgs+1
      ks_edgez = nguard*k3d*npgs+1
      ie_edgez = nguard*npgs+nxb + 1
      je_edgez = nguard*k2d*npgs+nyb + k2d
      ke_edgez = nguard*k3d*npgs+nzb

      is_unkn = nguard*npgs+1
      js_unkn = nguard*k2d*npgs+1
      ks_unkn = nguard*k3d*npgs+1
      ie_unkn = nguard*npgs+nxb + 1
      je_unkn = nguard*k2d*npgs+nyb + k2d
      ke_unkn = nguard*k3d*npgs+nzb + k3d

      End If  ! End If (first)

      If (new_loc(1,lb) == lb.And.new_loc(2,lb) == mype) Then
         If (.Not.moved(lb)) moved(lb) = .True.
         Return
      End If

      If (lb <= max(lnblocks_old,new_lnblocks)) Then

         If (lb <= lnblocks_old) Then
           If (new_loc(2,lb).ne.mype) Then
            success = .False.
            Call MPI_IPROBE (new_loc(2,lb),                            & 
                             maxblocks+new_loc(1,lb),                  & 
                             MPI_COMM_WORLD,                           & 
                             success,                                  & 
                             status,                                   & 
                             ierr)
            If (.Not.moved(lb).And.success) Then
               Call MPI_RECV (success,                                 & 
                              1,                                       & 
                              MPI_LOGICAL,                             & 
                              new_loc(2,lb),                           & 
                              maxblocks+new_loc(1,lb),                 & 
                              MPI_COMM_WORLD,                          & 
                              status,                                  & 
                              ierr)
               If (free(lb)) Then

                If (nvar > 0) Then
                Call MPI_SSEND (                                       & 
         unk(1,is_unk,js_unk,ks_unk,point_to(lb)),                     & 
                                1,                                     & 
                                unk_int_type,                          & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb),                         & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If  ! End If (nvar > 0)

!---------------send facevariables
                If (nfacevar > 0) Then
                Call MPI_SSEND (facevarx(1,is_facex,js_facex,ks_facex, &
                                point_to(lb)),                         & 
                                1,                                     & 
                                facex_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+2*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                If (ndim >= 2) Then
                Call MPI_SSEND (facevary(1,is_facey,js_facey,ks_facey, &
                                point_to(lb)),                         & 
                                1,                                     & 
                                facey_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+3*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If
                If (ndim == 3) Then
                Call MPI_SSEND (facevarz(1,is_facez,js_facez,ks_facez, &
                                point_to(lb)),                         & 
                                1,                                     & 
                                facez_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+4*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If

                End If  ! End If (nfacevar > 0)

!---------------send edge variables
                If (nvaredge > 0) Then
                Call MPI_SSEND (unk_e_x(1,is_edgex,js_edgex,ks_edgex,  &
                                point_to(lb)),                         & 
                                1,                                     & 
                                edgex_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+5*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                If (ndim >= 2) Then
                Call MPI_SSEND (unk_e_y(1,is_edgey,js_edgey,ks_edgey,  &
                                point_to(lb)),                         & 
                                1,                                     & 
                                edgey_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+6*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If
                If (ndim == 3) Then
                Call MPI_SSEND (unk_e_z(1,is_edgez,js_edgez,ks_edgez,  &
                                point_to(lb)),                         & 
                                1,                                     & 
                                edgez_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+7*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If
                End If  ! End If (nvaredge)

!---------------send corner variables
                If (nvarcorn > 0) Then
                Call MPI_SSEND (                                       &
                        unk_n(1,is_unkn,js_unkn,ks_unkn,point_to(lb)), & 
                                1,                                     & 
                                unkn_int_type,                         & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+8*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If  ! End If (nvarcorn > 0)
                test(point_to(lb)) = -1

               Else

                If (nvar > 0) Then
                Call MPI_SSEND (unk(1,is_unk,js_unk,ks_unk,lb),        & 
                                1,                                     & 
                                unk_int_type,                          & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb),                         & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If

!---------------send facevariables
                If (nfacevar > 0) Then
                Call MPI_SSEND (                                       &
                            facevarx(1,is_facex,js_facex,ks_facex,lb), & 
                                1,                                     & 
                                facex_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+2*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                If (ndim >= 2) Then
                Call MPI_SSEND (                                       &
                            facevary(1,is_facey,js_facey,ks_facey,lb), & 
                                1,                                     & 
                                facey_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+3*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If
                If (ndim == 3) Then
                Call MPI_SSEND (                                       &
                            facevarz(1,is_facez,js_facez,ks_facez,lb), & 
                                1,                                     & 
                                facez_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+4*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If
                End If  ! End If (nfacevar)

!---------------send edge variables
                If (nvaredge > 0) Then
                Call MPI_SSEND (                                       &
                             unk_e_x(1,is_edgex,js_edgex,ks_edgex,lb), & 
                                1,                                     & 
                                edgex_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+5*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                If (ndim >= 2) Then
                Call MPI_SSEND (                                       &
                             unk_e_y(1,is_edgey,js_edgey,ks_edgey,lb), & 
                                1,                                     & 
                                edgey_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+6*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If
                If (ndim == 3) Then
                Call MPI_SSEND (                                       &
                             unk_e_z(1,is_edgez,js_edgez,ks_edgez,lb), & 
                                1,                                     & 
                                edgez_int_type,                        & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+7*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If
                End If  ! End If (nvaredge > 0)

!---------------send corner variables
                If (nvarcorn > 0) Then
                Call MPI_SSEND (unk_n(1,is_unkn,js_unkn,ks_unkn,lb),   & 
                                1,                                     & 
                                unkn_int_type,                         & 
                                new_loc(2,lb),                         & 
                                new_loc(1,lb)+8*maxblocks,             & 
                                MPI_COMM_WORLD,                        & 
                                ierr)
                End If

                free(lb) = .True.
               End If  ! End If (free(lb))
               moved(lb) = .True.
            End If  ! End If (.Not.moved(lb).And.success)

           Else

            If (.Not.moved(lb).And.free(new_loc(1,lb))) Then
             If (free(lb)) Then

               If (nvar > 0) Then
      unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,new_loc(1,lb)) = & 
          unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,point_to(lb))
               End If

!--------------move facevars
               If (nfacevar > 0) Then
                  facevarx(:,is_facex:ie_facex,js_facex:je_facex,      &
                             ks_facex:ke_facex,new_loc(1,lb)) =        & 
                       facevarx(:,is_facex:ie_facex,js_facex:je_facex, &
                                  ks_facex:ke_facex,point_to(lb))
                  If (ndim >= 2) Then
                  facevary(:,is_facey:ie_facey,js_facey:je_facey,      &
                             ks_facey:ke_facey,new_loc(1,lb)) =        & 
                       facevary(:,is_facey:ie_facey,js_facey:je_facey, &
                                  ks_facey:ke_facey,point_to(lb))
                  End If
                  If (ndim == 3) Then
                  facevarz(:,is_facez:ie_facez,js_facez:je_facez,      &
                             ks_facez:ke_facez,new_loc(1,lb)) =        & 
                       facevarz(:,is_facez:ie_facez,js_facez:je_facez, &
                                  ks_facez:ke_facez,point_to(lb))
                  End If
               End If  ! End If (nfacevar > 0)

!--------------move edgevars
               If (nvaredge > 0) Then
                  unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,       &
                            ks_edgex:ke_edgex,new_loc(1,lb)) =         & 
                       unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
                                 ks_edgex:ke_edgex,point_to(lb))
                  If (ndim >= 2) Then
                  unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,       &
                            ks_edgey:ke_edgey,new_loc(1,lb)) =         & 
                       unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
                                 ks_edgey:ke_edgey,point_to(lb))
                  End If
                  If (ndim == 3) Then
                  unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,       &
                            ks_edgez:ke_edgez,new_loc(1,lb)) =         & 
                       unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
                                 ks_edgez:ke_edgez,point_to(lb))
                  End If
               End If  ! End If (nvaredge > 0)

               If (nvarcorn > 0) Then
                  unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,             &
                          ks_unkn:ke_unkn,new_loc(1,lb)) =             & 
                       unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,        &
                               ks_unkn:ke_unkn,point_to(lb))
               End If

               test(point_to(lb)) = -1

             Else

               If (nvar > 0) Then
      unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,new_loc(1,lb)) = & 
         unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,lb)
               End If

!--------------move facevars
               If (nfacevar > 0) Then
                  facevarx(:,is_facex:ie_facex,js_facex:je_facex,      &
                             ks_facex:ke_facex,new_loc(1,lb)) =        & 
                       facevarx(:,is_facex:ie_facex,js_facex:je_facex, &
                                  ks_facex:ke_facex,lb)
                  If (ndim >= 2) Then
                  facevary(:,is_facey:ie_facey,js_facey:je_facey,      &
                             ks_facey:ke_facey,new_loc(1,lb)) =        & 
                       facevary(:,is_facey:ie_facey,js_facey:je_facey, &
                                  ks_facey:ke_facey,lb)
                  End If
                  If (ndim == 3) Then
                  facevarz(:,is_facez:ie_facez,js_facez:je_facez,      &
                             ks_facez:ke_facez,new_loc(1,lb)) =        & 
                       facevarz(:,is_facez:ie_facez,js_facez:je_facez, &
                                  ks_facez:ke_facez,lb)
                  End If
               End If

!--------------move edgevars
               If (nvaredge > 0) Then
                  unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,       &
                            ks_edgex:ke_edgex,new_loc(1,lb)) =         & 
                       unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
                                 ks_edgex:ke_edgex,lb)
                  If (ndim >= 2) Then
                  unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,       &
                            ks_edgey:ke_edgey,new_loc(1,lb)) =         & 
                       unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
                                 ks_edgey:ke_edgey,lb)
                  End If
                  If (ndim == 3) Then
                  unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,       &
                            ks_edgez:ke_edgez,new_loc(1,lb)) =         &  
                       unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
                                 ks_edgez:ke_edgez,lb)
                  End If
               End If  ! End If (nvaredge > 0)

               If (nvarcorn > 0) Then
               unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,                &
                        ks_unkn:ke_unkn,new_loc(1,lb)) =               &
                 unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,              &
                         ks_unkn:ke_unkn,lb)
               End If

               free(lb) = .True.
             End If  ! End If (free(lb))

             moved(lb) = .True.
            End If  ! End If (.Not.moved(lb).And.free(new_loc(1,lb)))
           End If  ! End If (new_loc(2,lb).ne.mype)
         End If  ! End If (lb <= lnblocks_old)

         If (lb <= new_lnblocks) Then
            If (free(lb).And..Not.sent(lb)) Then
               sent(lb) = .True.
               If (.Not.newchild(lb)) Then
                  If (old_loc(2,lb).ne.mype) Then
                     nsend = nsend + 1
                     Call MPI_ISEND (free(lb),                         & 
                                     1,                                & 
                                     MPI_LOGICAL,                      & 
                                     old_loc(2,lb),                    & 
                                     maxblocks+lb,                     & 
                                     MPI_COMM_WORLD,                   & 
                                     reqs(nsend),                      & 
                                     ierr)
                  End If
               End If
            End If
         End If

         If (lb <= lnblocks_old.And..Not.free(lb)) Then
            nmoved = nmoved + 1
            point_to(lb) = max(lnblocks_old,new_lnblocks)+nmoved
            If (point_to(lb) > maxblocks) Then
               Do lb2 = max(lnblocks_old,new_lnblocks)+1,maxblocks
                  If (test(lb2) == -1) Then
                     point_to(lb) = lb2
                     Go To 22
                  End If
               End Do
            End If  
 22         If (point_to(lb) <= maxblocks) Then
               test(point_to(lb)) = 1

               If (nvar > 0) Then
       unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,point_to(lb)) = & 
               unk(:,is_unk:ie_unk,js_unk:je_unk,ks_unk:ke_unk,lb)
               End If

!--------------move facevars
               If (nfacevar > 0) Then
                  facevarx(:,is_facex:ie_facex,js_facex:je_facex,      &
                             ks_facex:ke_facex,point_to(lb)) =         & 
                       facevarx(:,is_facex:ie_facex,js_facex:je_facex, &
                                  ks_facex:ke_facex,lb)
                  If (ndim >= 2) Then
                  facevary(:,is_facey:ie_facey,js_facey:je_facey,      &
                             ks_facey:ke_facey,point_to(lb)) =         & 
                       facevary(:,is_facey:ie_facey,js_facey:je_facey, &
                                  ks_facey:ke_facey,lb)
                  End If
                  If (ndim == 3) Then
                  facevarz(:,is_facez:ie_facez,js_facez:je_facez,      &
                             ks_facez:ke_facez,point_to(lb)) =         & 
                       facevarz(:,is_facez:ie_facez,js_facez:je_facez, &
                                  ks_facez:ke_facez,lb)
                  End If
               End If  ! End If (nfacvar > 0)

!--------------move edgevars
               If (nvaredge > 0) Then
                  unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,       &
                            ks_edgex:ke_edgex,point_to(lb)) =          & 
                       unk_e_x(:,is_edgex:ie_edgex,js_edgex:je_edgex,  &
                                 ks_edgex:ke_edgex,lb)
                  If (ndim >= 2) Then
                  unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,       &
                            ks_edgey:ke_edgey,point_to(lb)) =          & 
                       unk_e_y(:,is_edgey:ie_edgey,js_edgey:je_edgey,  &
                                 ks_edgey:ke_edgey,lb)
                  End If
                  If (ndim == 3) Then
                  unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,       &
                            ks_edgez:ke_edgez,point_to(lb)) =          & 
                       unk_e_z(:,is_edgez:ie_edgez,js_edgez:je_edgez,  &
                                 ks_edgez:ke_edgez,lb)
                  End If
               End If  ! End If (nvaredge > 0)

               If (nvarcorn > 0) Then
               unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,                &
                       ks_unkn:ke_unkn,point_to(lb)) =                 &
                 unk_n(:,is_unkn:ie_unkn,js_unkn:je_unkn,              &
                         ks_unkn:ke_unkn,lb)
               End If

               free(lb) = .True.
            End If  ! End If (point_to(lb) <= maxblocks)
         End If  ! End If If (lb <= lnblocks_old.And..Not.free(lb))

         Return

      Else

         Return

      End If  ! End If (lb <= max(lnblocks_old,new_lnblocks))

      Return
      End Subroutine send_block_data

