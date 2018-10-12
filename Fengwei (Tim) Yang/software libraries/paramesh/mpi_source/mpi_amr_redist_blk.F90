!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fI mpi_source/amr_redist_blk
!! NAME
!!
!!   amr_redist_blk
!! 
!! SYNOPSIS
!!
!!   call amr_redist_blk (new_loc, nprocs, mype, lnblocks_old)
!!
!!   call amr_redist_blk (integer, integer, integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(inout) :: new_loc(:,:)
!!     Array which stores the new locations in memory where blocks are to be
!!     moved.  new_loc(1,:) indicates the on-processor location where the block
!!     is to reside and new_loc(2,:) indicates which processor the block data is
!!     to be moved to.
!!
!!   integer, intent(in) :: nprocs
!!     The number of processors being used.
!!   
!!   integer, intent (in) :: lnblocks_old
!!     The 'old' number of blocks in the calling processor before data redistribution.
!!
!!   integer, intent(in) :: mype     
!!      Current processor number
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
!!   paramesh_interfaces
!!   paramesh_comm_data
!!
!! CALLS
!!
!!   fill_old_loc
!!   send_block_data
!!
!! RETURNS
!!
!!   Upon exit block data has been redistributed to processors according to
!!   a space filling morton curve.
!!
!! DESCRIPTION
!!
!!   This routine redistributes blocks to processors as part of the refinement/
!!   derefinement process.  This routine actually moves block data to their new
!!   locations if any refinements or derefinements have occured.
!!
!!   For each block a location in memory is passed in using the array new_loc.
!!   new_loc is the new location in memory where a block's data (e.g. in unk, 
!!   facevarx,y,z, ... etc.) is to be moved.  For instance, new_loc(1,lb) is the
!!   location in memory and new_loc(2,lb) is the processor to move the data to
!!   for block lb.
!!
!!   The routine is called by the subroutine amr_refine_derefine and should never 
!!   need to be called by a user's application.
!!
!! AUTHORS
!!
!!   Kevin Olson (1998)
!!   Bug fix contributed by Paul Ricker and Marcus Gross (2003)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f, unk_test, face[xyz]_test, edge[xyz]_test, unkn_test
#include "paramesh_preprocessor.fh"

      Subroutine amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use paramesh_comm_data
      Use paramesh_interfaces, only : fill_old_loc

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(inout) :: new_loc(:,:)
      Integer, Intent(in)    :: nprocs,mype,lnblocks_old

!-----Local variables. and arrays.
      Integer :: lb,ierr,errorcode
      Integer :: nsend, nrecv
!      Integer :: old_loc(2,maxblocks_tr)
!      Integer :: reqr(maxblocks_tr)
!      Integer :: reqs(maxblocks_tr)
!      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
!      Integer :: stats(MPI_STATUS_SIZE,maxblocks_tr)
!      Integer :: test(maxblocks), point_to(maxblocks)
      Integer,allocatable :: old_loc(:,:)
      Integer,allocatable :: reqr(:)
      Integer,allocatable :: reqs(:)
      Integer,allocatable :: statr(:,:)
      Integer,allocatable :: stats(:,:)
      Integer,allocatable :: test(:), point_to(:)

      Integer :: nmoved, nit
      Integer :: nm, nm2, nm2_old
      Integer :: ireduce_datain(1),ireduce_dataout(1)
      Integer, Save :: unk_int_type
      Integer, Save ::  is_unk,js_unk,ks_unk
      Integer, Save :: facex_int_type, facey_int_type, facez_int_type
      Integer, Save :: is_facex,js_facex,ks_facex
      Integer, Save :: is_facey,js_facey,ks_facey
      Integer, Save :: is_facez,js_facez,ks_facez
      Integer, Save :: edgex_int_type, edgey_int_type, edgez_int_type
      Integer, Save :: is_edgex,js_edgex,ks_edgex
      Integer, Save :: is_edgey,js_edgey,ks_edgey
      Integer, Save :: is_edgez,js_edgez,ks_edgez
      Integer, Save :: unkn_int_type
      Integer, Save ::  is_unkn,js_unkn,ks_unkn
      Integer :: type1, type2, type3
      Integer :: udim(4), udim_tot(4), i
      Integer :: nbytes
      Logical :: lreduce_datain(1),lreduce_dataout(1)
      Logical, Save :: first = .True.
!      Logical :: free(maxblocks), moved(maxblocks), sent(maxblocks)
      Logical,allocatable :: free(:), moved(:), sent(:)
      Logical :: repeat, repeatt
      integer :: madetype(8)

! CEG allocate memory
  allocate(old_loc(2,maxblocks_tr))
  allocate(reqr(maxblocks_tr))
  allocate(reqs(maxblocks_tr))
  allocate(statr(MPI_STATUS_SIZE,maxblocks_tr))
  allocate(stats(MPI_STATUS_SIZE,maxblocks_tr))
  allocate(test(maxblocks))
  allocate(point_to(maxblocks))
  allocate(free(maxblocks))
  allocate(moved(maxblocks))
  allocate(sent(maxblocks))

      If (first) Then

         is_unk = nguard*npgs+1
         js_unk = nguard*k2d*npgs+1
         ks_unk = nguard*k3d*npgs+1

         is_unkn = nguard*npgs+1
         js_unkn = nguard*k2d*npgs+1
         ks_unkn = nguard*k3d*npgs+1

         is_facex = nguard*npgs+1
         js_facex = nguard*k2d*npgs+1
         ks_facex = nguard*k3d*npgs+1

         is_facey = nguard*npgs+1
         js_facey = nguard*k2d*npgs+1
         ks_facey = nguard*k3d*npgs+1

         is_facez = nguard*npgs+1
         js_facez = nguard*k2d*npgs+1
         ks_facez = nguard*k3d*npgs+1

         is_edgex = nguard*npgs+1
         js_edgex = nguard*k2d*npgs+1
         ks_edgex = nguard*k3d*npgs+1

         is_edgey = nguard*npgs+1
         js_edgey = nguard*k2d*npgs+1
         ks_edgey = nguard*k3d*npgs+1

         is_edgez = nguard*npgs+1
         js_edgez = nguard*k2d*npgs+1
         ks_edgez = nguard*k3d*npgs+1

         call amr_redist_blk_mpitypeset(madetype)

      unk_int_type=madetype(1)
      unkn_int_type=madetype(2)
      edgex_int_type=madetype(3)
      edgey_int_type=madetype(4)
      edgez_int_type=madetype(5)
      facex_int_type=madetype(6)
      facey_int_type=madetype(7)
      facez_int_type=madetype(8)

      first=.false.

      End If  ! End If (first)


!-----compute old_loc
      Call fill_old_loc (new_loc,old_loc,nprocs,mype)

      nrecv = 0
      nsend = 0


!-----treat unk
      If (nvar > 0) Then
!--------Post all receives for unk
         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
                  Call MPI_IRECV (unk(1,is_unk,js_unk,ks_unk,lb),      & 
                                  1,                                   & 
                                  unk_int_type,                        & 
                                  old_loc(2,lb),                       & 
                                  lb,                                  & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)

               End If
            End If
         End Do

      End If

!-----Treat Facevariables
      If (nfacevar > 0) Then

!--------Treat facevarx
         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
            Call MPI_IRECV (facevarx(1,is_facex,js_facex,ks_facex,lb), & 
                                  1,                                   & 
                                  facex_int_type,                      & 
                                  old_loc(2,lb),                       & 
                                  lb+2*maxblocks,                      & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)
               End If
            End If
         End Do

!--------Treat facevary
         If (ndim >= 2) Then
         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
            Call MPI_IRECV (facevary(1,is_facey,js_facey,ks_facey,lb), & 
                                  1,                                   & 
                                  facey_int_type,                      & 
                                  old_loc(2,lb),                       & 
                                  lb+3*maxblocks,                      & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)
               End If
            End If
         End Do
         End If

!--------Treat Facevarz
         If (ndim == 3) Then
         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
            Call MPI_IRECV (facevarz(1,is_facez,js_facez,ks_facez,lb), & 
                                  1,                                   & 
                                  facez_int_type,                      & 
                                  old_loc(2,lb),                       & 
                                  lb+4*maxblocks,                      & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)
               End If
            End If
         End Do
         End If

      End If  ! End If (nfacevar > 0)

!-----Treat Edge variables
      If (nvaredge > 0) Then

!--------Treat unk_e_x
         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
             Call MPI_IRECV (unk_e_x(1,is_edgex,js_edgex,ks_edgex,lb), & 
                                  1,                                   & 
                                  edgex_int_type,                      & 
                                  old_loc(2,lb),                       & 
                                  lb+5*maxblocks,                      & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)
               End If
            End If
         End Do

!--------Treat unk_e_y
         If (ndim >= 2) Then
         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
             Call MPI_IRECV (unk_e_y(1,is_edgey,js_edgey,ks_edgey,lb), & 
                                  1,                                   & 
                                  edgey_int_type,                      & 
                                  old_loc(2,lb),                       & 
                                  lb+6*maxblocks,                      & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)
               End If
            End If
         End Do
         End If

!--------Treat unk_e_z
         If (ndim == 3) Then
         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
             Call MPI_IRECV (unk_e_z(1,is_edgez,js_edgez,ks_edgez,lb), & 
                                  1,                                   & 
                                  edgez_int_type,                      & 
                                  old_loc(2,lb),                       & 
                                  lb+7*maxblocks,                      & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)
               End If
            End If
         End Do
         End If

      End If  ! End If (nvaredge > 0)

!-----treat unk_n
      If (nvarcorn > 0) Then

         Do lb = 1,new_lnblocks
            If (.Not.newchild(lb)) Then
               If (old_loc(2,lb).ne.mype) Then
                  nrecv = nrecv + 1
                  Call MPI_IRECV (unk_n(1,is_unkn,js_unkn,ks_unkn,lb), & 
                                  1,                                   & 
                                  unkn_int_type,                       & 
                                  old_loc(2,lb),                       & 
                                  lb+8*maxblocks,                      & 
                                  MPI_COMM_WORLD,                      & 
                                  reqr(nrecv),                         & 
                                  ierr)
               End If
            End If
         End Do

      End If  ! End If (nvarcorn > 0)

      moved(:) = .False.
      moved(lnblocks_old+1:maxblocks) = .True.
      free(:) = .False.
      free(lnblocks_old+1:maxblocks) = .True.
      sent(:) = .False.
      repeat = .True.
      nmoved = 0 
      test(:) = 0
      point_to(:) = 0
      
      nit = 0
      nm2 = 0
      nm2_old = 1
      Do While (repeat.And.nit<=100) 
         
         Do lb = 1,max(lnblocks_old,new_lnblocks)
            Call send_block_data (lb, new_loc, old_loc, free,          & 
                                  moved, sent,                         & 
                                  lnblocks_old, mype, nmoved,          & 
                                  test, point_to,                      & 
                                  reqs, nsend, unk_int_type,           &
                                  facex_int_type, facey_int_type,      &
                                  facez_int_type, edgex_int_type,      &
                                  edgey_int_type, edgez_int_type,      &
                                  unkn_int_type)
         End Do
!CEG - moved to here from below
!      If (nsend > 0) Then
!         Call MPI_WAITALL(nsend,reqs,stats,ierr)
!      End If

         repeat = any(.Not.moved(:))
         lreduce_datain(1) = repeat
         Call MPI_ALLREDUCE(                                           & 
                 lreduce_datain(1),lreduce_dataout(1),                 & 
                 1,MPI_LOGICAL,                                        & 
                 MPI_LOR,MPI_COMM_WORLD,ierr)
         repeatt = lreduce_dataout(1)
         repeat = repeatt
         
         nm2_old = nm2
         nm = count(.Not.moved(:))
         ireduce_datain(1) = nm
         Call MPI_ALLREDUCE(                                           & 
              ireduce_datain(1),ireduce_dataout(1),                    & 
              1,MPI_INTEGER,                                           & 
              MPI_SUM,MPI_COMM_WORLD,ierr)
         nm2 = ireduce_dataout(1)
!         If (mype == 0) Then
!            Print *,' iteration, no. not moved = ',nit,nm2
!            Print *,' iteration, no. not moved local,pid = ',nit,nm,mype
!            if (nm.gt.0) print*,moved(:)
!         End If
         
         nit = nit + 1
         
      End Do  ! End Do While (repeat.And.nit<=100)
      
      If (nm2_old == nm2.And.nm2.ne.0.And.nit>=100) Then
         If (mype == 0) Then
          Print *,' ERROR: could not move all blocks in amr_redist_blk'
          Print *,' Try increasing maxblocks or use more processors'
          Print *,' nm2_old, nm2 = ',nm2_old,nm2
          Print *,' ABORTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         End If
         Call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
      End If
      
      If (nrecv > 0) Then
         Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

!-----NOTE: Bug fix by Paul Ricker and Marcus Gross (5/2003).  Added Waitall 
!-----to isends so that SGI MPI buffers do not overflow.

      If (nsend > 0) Then
         Call MPI_WAITALL(nsend,reqs,stats,ierr)
      End If

! CEG deallocate memory
  deallocate(old_loc)
  deallocate(reqr)
  deallocate(reqs)
  deallocate(statr)
  deallocate(stats)
  deallocate(test)
  deallocate(point_to)
  deallocate(free)
  deallocate(moved)
  deallocate(sent)


      Return
      End Subroutine amr_redist_blk
      


! CEG Move the mpi typesetting to a new routine



  subroutine amr_redist_blk_mpitypeset(madetype)
!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use paramesh_comm_data
      Use paramesh_interfaces, only : fill_old_loc

      Implicit None

!-----Include statements
      Include 'mpif.h'

! Input/Output variable set
      integer, intent(inout) :: madetype(8)

!-----Local variables. and arrays.
      Logical, Save :: first = .True.

! Definitely used here
      Integer :: ierr
      Integer :: nbytes
      Integer, Allocatable :: unkn_test(:,:,:,:)
      Integer :: udim(4), udim_tot(4), i
      Integer :: type1, type2, type3
      Integer, Save :: unkn_int_type
      Integer, Save :: unk_int_type
!Definitely local here
      Integer, Allocatable :: unk_test(:,:,:,:)
      Integer, Save :: facex_int_type, facey_int_type, facez_int_type
      Integer, Allocatable :: facex_test(:,:,:,:), facey_test(:,:,:,:),&
                              facez_test(:,:,:,:)
      Integer, Save :: edgex_int_type, edgey_int_type, edgez_int_type
      Integer, Allocatable :: edgex_test(:,:,:,:), edgey_test(:,:,:,:),&
                              edgez_test(:,:,:,:)

! do we need to deallocate the datatypes?
      if (.not.first) then
         if (nvar > 0) call mpi_type_free(unk_int_type, ierr)
         if (nvarcorn > 0) call mpi_type_free(unkn_int_type, ierr)
         if (nfacevar > 0) call mpi_type_free(facex_int_type, ierr)
         if (nfacevar > 0) call mpi_type_free(facey_int_type, ierr)
         if (nfacevar > 0) call mpi_type_free(facez_int_type, ierr)
         if (nvaredge > 0) call mpi_type_free(edgex_int_type, ierr)
         if (nvaredge > 0) call mpi_type_free(edgey_int_type, ierr)
         if (nvaredge > 0) call mpi_type_free(edgez_int_type, ierr)

         return
      endif

#ifdef REAL8
         nbytes = 8
#else
         nbytes = 4
#endif

      If (nvar > 0) Then

!-----This code added for reordering script
      Allocate(unk_test(nvar,nxb,nyb,nzb))
      Do i = 1,4
         udim_tot(i) = size(unk,dim=i) 
         udim(i) = size(unk_test,dim=i)
      End Do
      Deallocate(unk_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              & 
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  & 
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      unk_int_type = type3

      Call MPI_TYPE_COMMIT(unk_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)

      End If  ! End If (nvar > 0)

      If (nfacevar > 0) Then

!-----the following code added for reordering script
      Allocate(facex_test(nfacevar,nxb+1,nyb,nzb))
      Do i = 1,4
         udim_tot(i) = size(facevarx,dim=i) 
         udim(i) = size(facex_test,dim=i)
      End Do
      Deallocate(facex_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              & 
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  &  
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                            udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      facex_int_type = type3

      Call MPI_TYPE_COMMIT(facex_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)


      Allocate(facey_test(nfacevar,nxb,nyb+k2d,nzb))
      Do i = 1,4
         udim_tot(i) = size(facevary,dim=i) 
         udim(i) = size(facey_test,dim=i)
      End Do
      Deallocate(facey_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              &  
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  & 
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                            udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      facey_int_type = type3

      Call MPI_TYPE_COMMIT(facey_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)


      Allocate(facez_test(nfacevar,nxb,nyb,nzb+k3d))
      Do i = 1,4
         udim_tot(i) = size(facevarz,dim=i) 
         udim(i) = size(facez_test,dim=i)
      End Do
      Deallocate(facez_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              & 
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  & 
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                            udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      facez_int_type = type3

      Call MPI_TYPE_COMMIT(facez_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)

      End If  ! End If (nfacevar > 0)

      If (nvaredge > 0) Then

!-----Code added for reordering script
      Allocate(edgex_test(nvaredge,nxb,nyb+k2d,nzb+k3d))
      Do i = 1,4
         udim_tot(i) = size(unk_e_x,dim=i) 
         udim(i) = size(edgex_test,dim=i)
      End Do
      Deallocate(edgex_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              & 
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  & 
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                            udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      edgex_int_type = type3

      Call MPI_TYPE_COMMIT(edgex_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)


      Allocate(edgey_test(nvaredge,nxb+1,nyb,nzb+k3d))
      Do i = 1,4
         udim_tot(i) = size(unk_e_y,dim=i) 
         udim(i) = size(edgey_test,dim=i)
      End Do
      Deallocate(edgey_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              & 
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  & 
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                            udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      edgey_int_type = type3

      Call MPI_TYPE_COMMIT(edgey_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)


      Allocate(edgez_test(nvaredge,nxb+1,nyb+k2d,nzb))
      Do i = 1,4
         udim_tot(i) = size(unk_e_z,dim=i) 
         udim(i) = size(edgez_test,dim=i)
      End Do
      Deallocate(edgez_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              & 
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  & 
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                            udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      edgez_int_type = type3

      Call MPI_TYPE_COMMIT(edgez_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)

      End If  ! End If (nvaredge > 0)

      If (nvarcorn > 0) Then

!-----Code added for reordering script.
      Allocate(unkn_test(nvarcorn,nxb+1,nyb+k2d,nzb+k3d))
      Do i = 1,4
         udim_tot(i) = size(unk_n,dim=i) 
         udim(i) = size(unkn_test,dim=i)
      End Do
      Deallocate(unkn_test)

!-----DEFINE BLOCK INTERIOR
      Call MPI_TYPE_VECTOR (udim(2),                                   & 
                            udim(1),                                   & 
                            udim_tot(1),                               & 
                            amr_mpi_real,                              & 
                            type1,                                     & 
                            ierr)
      Call MPI_TYPE_HVECTOR (udim(3),                                  & 
                             1,                                        &
                             udim_tot(1)*udim_tot(2)*nbytes,           & 
                             type1,                                    & 
                             type2,                                    & 
                             ierr)
      Call MPI_TYPE_HVECTOR (udim(4),                                  &
                             1,                                        &
                            udim_tot(1)*udim_tot(2)*udim_tot(3)*nbytes,&
                             type2,                                    &
                             type3,                                    &
                             ierr)

      unkn_int_type = type3

      Call MPI_TYPE_COMMIT(unkn_int_type,ierr)
! CEG added frees to clean up unneeded types
      call MPI_TYPE_FREE(type1, ierr)
      call MPI_TYPE_FREE(type2, ierr)

      End If  ! End If (nvarcorn > 0)


! Assign types to input/output array
      madetype(1)=unk_int_type
      madetype(2)=unkn_int_type
      madetype(3)=edgex_int_type
      madetype(4)=edgey_int_type
      madetype(5)=edgez_int_type
      madetype(6)=facex_int_type
      madetype(7)=facey_int_type
      madetype(8)=facez_int_type


  first = .false.

  end subroutine amr_redist_blk_mpitypeset
