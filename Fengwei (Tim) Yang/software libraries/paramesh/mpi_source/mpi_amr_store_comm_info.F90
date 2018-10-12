!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!------------------
! Guardcell filling

      Subroutine mpi_amr_Write_guard_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs

!-----Local Variables
      Integer             :: i,istat,ierrorcode,ierr

!-----Begin Executable code
! TEST
      to_be_received(2,:,:) = -1
      to_be_sent(2,:,:) = -1
      pe_source = -1
!!!

      If (.Not.Allocated(pe_source_guard))                             & 
         Allocate(pe_source_guard(1:nprocs),stat = istat) 

      If (Allocated(to_be_received_guard))                             & 
                       Deallocate(to_be_received_guard)
      If (Allocated(to_be_received)) Then
        Allocate(to_be_received_guard(3,Size(to_be_received,2),        & 
               Size(to_be_received,3)),stat = istat) 
      Else
        Write(*,*)'Paramesh error: Write_guard ',                      & 
                  'to_be_received is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If

      to_be_received_guard = to_be_received

      If (Allocated(to_be_sent_guard))                                 & 
                       Deallocate(to_be_sent_guard)
      If (Allocated(to_be_sent)) Then
        Allocate(to_be_sent_guard(3,Size(to_be_sent,2),                & 
               Size(to_be_sent,3)),stat = istat) 
      Else
        Write(*,*)'Paramesh error: Write_guard ',                      & 
                  'to_be_sent is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If

      to_be_sent_guard = to_be_sent

      If (Allocated(commatrix_guard)) Deallocate(commatrix_guard)
      Allocate(commatrix_guard(1:nprocs,2),                            & 
                 stat = istat) 
       If (istat.ne.0) Then
         Write(*,*) 'mpi_amr_write_guard_comm: allocation error'
         Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
       End If

      pe_source_guard = pe_source
      commatrix_guard = 0

      Do i = 1,nprocs
         commatrix_guard(i,1) = commatrix_recv(i)
         commatrix_guard(i,2) = commatrix_send(i)
      End Do

      strt_guard = strt_buffer
      max_no_to_send_guard = max_no_to_send

      laddress_guard = laddress

      Return
      End Subroutine mpi_amr_write_guard_comm





      Subroutine mpi_amr_read_guard_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs

!-----Local Variables
      Integer             :: i,istat,ierrorcode,ierr

!-----Begin Executable Code
      If (grid_analysed_mpi.ne.1) Then
        Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
        Call amr_morton_process()
      End If

      pe_source = pe_source_guard
!      commatrix_recv(:) = 0
!      commatrix_send(:) = 0

      Do i = 1,nprocs
        commatrix_recv(i) = commatrix_guard(i,1)
        commatrix_send(i) = commatrix_guard(i,2)
      End Do

      strt_buffer = strt_guard
      max_no_to_send = max_no_to_send_guard

      laddress = laddress_guard


      If (max_no_to_send > 0) Then

      If (Allocated(to_be_received))                                   & 
                       Deallocate(to_be_received)
      Allocate(to_be_received(3,Size(to_be_received_guard,2),          & 
               Size(to_be_received_guard,3)),stat = istat) 
      to_be_received = to_be_received_guard

      If (Allocated(to_be_sent)) Deallocate(to_be_sent)
      Allocate(                                                        & 
        to_be_sent(3,                                                  &
                   Size(to_be_sent_guard,2),                           &
                   Size(to_be_sent_guard,3)),                          & 
                 stat = istat) 
       If (istat.ne.0) Then
         Write(*,*) 'mpi_amr_read_guard_comm: allocation error'
         Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
       End If

      to_be_sent(:,:,:) = to_be_sent_guard(:,:,:)

      End If  ! End If (max_no_to_send > 0)


      mpi_pattern_id = 10

      Call amr_1blk_guardcell_reset

      Return
      End Subroutine mpi_amr_read_guard_comm


!-----------------
! Prolongation
  
      Subroutine mpi_amr_write_prol_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs

!-----Local Variables
      Integer             :: i,istat,ierrorcode,ierr


!-----Begin Executable code
! TEST
      to_be_received(2,:,:) = -1
      to_be_sent(2,:,:) = -1
      pe_source = -1
!!!

      If (.not.Allocated(pe_source_prol))                              & 
            Allocate(pe_source_prol(1:nprocs),stat = istat)

      If (Allocated(to_be_received_prol))                              & 
                       Deallocate(to_be_received_prol)
      If (Allocated(to_be_received)) Then
       Allocate(to_be_received_prol(3,Size(to_be_received,2),          & 
               Size(to_be_received,3)),stat = istat)
      Else
        Write(*,*)'Paramesh error: write_prol ',                       & 
                  'to_be_received is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      to_be_received_prol = to_be_received

      If (Allocated(to_be_sent_prol))                                  & 
                       Deallocate(to_be_sent_prol)
      If (Allocated(to_be_sent)) Then
       Allocate(to_be_sent_prol(3,Size(to_be_sent,2),                  & 
               Size(to_be_sent,3)),stat = istat)
      Else
        Write(*,*)'Paramesh error: write_prol ',                       & 
                  'to_be_sent is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      to_be_sent_prol = to_be_sent

      If (Allocated(commatrix_prol)) Deallocate(commatrix_prol)
      Allocate(commatrix_prol(1:nprocs,2),                             & 
                 stat = istat)
       If (istat.ne.0) Then
         Write(*,*) 'mpi_amr_write_prol_comm: allocation error'
         Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
       End If

      pe_source_prol = pe_source

      commatrix_prol = 0

      Do i = 1 ,nprocs
         commatrix_prol(i,1) = commatrix_recv(i)
         commatrix_prol(i,2) = commatrix_send(i)
      End Do

      strt_prol = strt_buffer
      max_no_to_send_prol = max_no_to_send

      laddress_prol = laddress

      Return
      End Subroutine mpi_amr_write_prol_comm





      Subroutine mpi_amr_read_prol_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs

!-----Local Variables
      Integer             :: i,istat,ierrorcode,ierr

!-----Begin Executable code
      If (grid_analysed_mpi.ne.1) Then
        Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
        Call amr_morton_process()
      End If

      pe_source = pe_source_prol

!      commatrix_recv(:) = 0
!      commatrix_send(:) = 0

      Do i = 1,nprocs
        commatrix_recv(i) = commatrix_prol(i,1)
        commatrix_send(i) = commatrix_prol(i,2)
      End Do

      strt_buffer = strt_prol
      max_no_to_send = max_no_to_send_prol

      laddress = laddress_prol

      If (max_no_to_send > 0) Then

      If (Allocated(to_be_received))                                   & 
                       Deallocate(to_be_received)
      If (Allocated(to_be_received_prol)) Then
        Allocate(to_be_received(3,Size(to_be_received_prol,2),         & 
               Size(to_be_received_prol,3)),stat = istat)
      Else
        Write(*,*)'Paramesh error: to_be_received_prol ',              & 
                  'is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      to_be_received = to_be_received_prol

      If (Allocated(to_be_sent)) Deallocate(to_be_sent)
      Allocate(                                                        & 
        to_be_sent(3,Size(to_be_sent_prol,2),                          &
                     Size(to_be_sent_prol,3)),                         & 
                     stat = istat) 
       If (istat.ne.0) Then
         Write(*,*) 'mpi_amr_read_prol_comm: allocation error'
         Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
       End If

      to_be_sent(:,:,:) = to_be_sent_prol(:,:,:)

      End If

      mpi_pattern_id = 20

      Return
      End Subroutine mpi_amr_read_prol_comm
 

 
!------------------
! Flux Conservation
  
      Subroutine mpi_amr_write_flux_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs

!-----Local variables
      Integer             :: i,istat,ierrorcode,ierr


!-----Begin Executable code
! TEST
      to_be_received(2,:,:) = -1
      to_be_sent(2,:,:) = -1
      pe_source = -1
!!!

      If (.not.Allocated(pe_source_flux))                              & 
            Allocate(pe_source_flux(1:nprocs),stat = istat)

      If (Allocated(to_be_received_flux))                              & 
                       Deallocate(to_be_received_flux)
      If (Allocated(to_be_received)) Then
       Allocate(to_be_received_flux(3,Size(to_be_received,2),          & 
               Size(to_be_received,3)),stat = istat)
      Else
        Write(*,*)'Paramesh error: write_flux ',                       & 
                  'to_be_received is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      to_be_received_flux = to_be_received

      If (Allocated(to_be_sent_flux))                                  & 
                       Deallocate(to_be_sent_flux)
      If (Allocated(to_be_sent)) Then
       Allocate(to_be_sent_flux(3,Size(to_be_sent,2),                  & 
               Size(to_be_sent,3)),stat = istat)
      Else
        Write(*,*)'Paramesh error: write_flux ',                       & 
                  'to_be_sent is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      to_be_sent_flux = to_be_sent

      If (Allocated(commatrix_flux)) Deallocate(commatrix_flux)
      Allocate(commatrix_flux(1:nprocs,2),                             & 
                 stat = istat) 
       If (istat.ne.0) Then
         Write(*,*) 'mpi_amr_write_flux_comm: allocation error'
         Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
       End If

      pe_source_flux = pe_source

      commatrix_flux = 0

      Do i = 1 ,nprocs
         commatrix_flux(i,1) = commatrix_recv(i)
         commatrix_flux(i,2) = commatrix_send(i)
      End Do

      strt_flux = strt_buffer
      max_no_to_send_flux = max_no_to_send

      laddress_flux = laddress


      Return
      End Subroutine mpi_amr_write_flux_comm





      Subroutine mpi_amr_read_flux_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs

!-----Local variables
      Integer             :: i,istat,ierrorcode,ierr

!-----Begin Executable code
      If (grid_analysed_mpi.ne.1) Then
        Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
        Call amr_morton_process()
      End If

      pe_source = pe_source_flux

!      commatrix_recv(:) = 0
!      commatrix_send(:) = 0

      Do i = 1,nprocs
        commatrix_recv(i) = commatrix_flux(i,1)
        commatrix_send(i) = commatrix_flux(i,2)
      End Do

      strt_buffer = strt_flux
      max_no_to_send = max_no_to_send_flux

      laddress = laddress_flux

      If (max_no_to_send > 0) Then

      If (Allocated(to_be_received))                                   & 
                       Deallocate(to_be_received)
      If (Allocated(to_be_received_flux)) Then
        Allocate(to_be_received(3,Size(to_be_received_flux,2),         & 
               Size(to_be_received_flux,3)),stat = istat)
      Else
        Write(*,*)'Paramesh error: to_be_received_flux ',              & 
                  'is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If
      to_be_received = to_be_received_flux

      If (Allocated(to_be_sent)) Deallocate(to_be_sent)
      Allocate(                                                        & 
        to_be_sent(3,Size(to_be_sent_flux,2),Size(to_be_sent_flux,3)), & 
                 stat = istat) 
       If (istat.ne.0) Then
         Write(*,*) 'mpi_amr_read_flux_comm: allocation error'
         Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
       End If

       to_be_sent(:,:,:) = to_be_sent_flux(:,:,:)

      End If

      mpi_pattern_id = 30

      Return
      End Subroutine mpi_amr_read_flux_comm
 


!------------------
! Restriction

 
      Subroutine mpi_amr_write_restrict_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Statements
      Integer, intent(in) :: nprocs

!-----Local Variables
      Integer             :: i,istat,ierrorcode,ierr

!-----Begin Executable Code
! TEST
      to_be_received(2,:,:) = -1
      to_be_sent(2,:,:) = -1
      pe_source = -1
!!!

      If (.Not.Allocated(pe_source_restrict))                          & 
            Allocate(pe_source_restrict(1:nprocs),                     & 
                     stat = istat)
      pe_source_restrict(:) = pe_source(:)


      If (Allocated(to_be_received)) Then
         if (Allocated(to_be_received_restrict)) Then
            Deallocate(to_be_received_restrict)
         End If
         Allocate(to_be_received_restrict(3,                           & 
                                          Size(to_be_received,2),      & 
                                          Size(to_be_received,3)),     & 
                                          stat = istat)
         to_be_received_restrict(:,:,:) = to_be_received(:,:,:)
      Else
        Write(*,*)'Paramesh error: write_restrict ',                   & 
                  'to_be_received is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If

      If (Allocated(to_be_sent)) Then
         if (Allocated(to_be_sent_restrict)) Then
             Deallocate(to_be_sent_restrict)
         End If
         Allocate(to_be_sent_restrict(3,                               & 
                                      Size(to_be_sent,2),              & 
                                      Size(to_be_sent,3)),             & 
                                      stat = istat)
         to_be_sent_restrict(:,:,:) = to_be_sent(:,:,:)
      Else
        Write(*,*)'Paramesh error: write_restrict ',                   & 
                  'to_be_sent is not Allocated.'
        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      End If


      If (.not.Allocated(commatrix_restrict)) Then
         Allocate(commatrix_restrict(1:nprocs,2),                      & 
                  stat = istat)
         If (istat.ne.0) Then
            Write(*,*) 'mpi_amr_write_restrict_comm: allocation error'
            Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
         End If
      End If
      Do i = 1 ,nprocs
!         commatrix_restrict(i,:) = 0
         commatrix_restrict(i,1) = commatrix_recv(i)
         commatrix_restrict(i,2) = commatrix_send(i)
      End Do


      strt_restrict = strt_buffer
      max_no_to_send_restrict = max_no_to_send

      laddress_restrict(:,:) = laddress(:,:)

      Return
      End Subroutine mpi_amr_write_restrict_comm





      Subroutine mpi_amr_read_restrict_comm(nprocs)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton

      Implicit None
!-----Include Statements
      Include 'mpif.h'

!-----Input/Output Statements
      Integer, intent(in) :: nprocs

!-----Local Variables
      Integer             :: i,istat,ierrorcode,ierr

!-----Begin Executable Code
      If (grid_analysed_mpi.ne.1) Then
        Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
        Call amr_morton_process()
      End If

      pe_source(:) = pe_source_restrict(:)

!      commatrix_recv(:) = 0
!      commatrix_send(:) = 0

      Do i = 1,nprocs
        commatrix_recv(i) = commatrix_restrict(i,1)
        commatrix_send(i) = commatrix_restrict(i,2)
      End Do

      strt_buffer = strt_restrict
      max_no_to_send = max_no_to_send_restrict


      laddress(:,:) = laddress_restrict(:,:)

      If (max_no_to_send > 0) Then

      If (Allocated(to_be_received))                                   & 
                       Deallocate(to_be_received)         
      If (Allocated(to_be_received_restrict)) Then
         Allocate(to_be_received(3,                                    & 
                  Size(to_be_received_restrict,2),                     & 
                  Size(to_be_received_restrict,3)),                    & 
                  stat = istat)
         to_be_received(:,:,:) = to_be_received_restrict(:,:,:)
      End If

      If (Allocated(to_be_sent))                                       & 
                       Deallocate(to_be_sent)
      If (Allocated(to_be_sent_restrict)) Then
         Allocate(to_be_sent(3,                                        &
                             Size(to_be_sent_restrict,2),              &
                             Size(to_be_sent_restrict,3)),             &
                             stat = istat)
         to_be_sent(:,:,:) = to_be_sent_restrict(:,:,:)
      End If

      Call amr_1blk_guardcell_reset

      End If

      mpi_pattern_id = 40

      Return
      End Subroutine mpi_amr_read_restrict_comm
