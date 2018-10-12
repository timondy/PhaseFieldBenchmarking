!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_checkpoint_re_default
!! NAME
!!
!!   amr_checkpoint_re_default
!!
!! SYNOPSIS
!!
!!   call amr_checkpoint_re_default(file_num)
!!   call amr_checkpoint_re_default(file_num, l_with_guardcells, 
!!                                  user_attr1, user_attr2, user_attr3,
!!                                  user_attr4, user_attr5)
!!
!!   call amr_checkpoint_re_default(integer, optional logical, 
!!                                  optional real, optional real, optional real,
!!                                  optional real, optional real)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: file_num
!!     An integer number which be appended to the end of the file name.
!!
!!   optional, logical, intent(in) :: l_with_guardcells
!!     If true, then guardcells are included in the checkpoint file.  Otherwise
!!     (the default) they are not included.
!!
!!   optional, real, intent(in) :: user_attr1(2,3,4,5)
!!     Arguments which allow the user to add some extra information to the file.
!!     Currently only 5 real numbers can be added.
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
!!   timings
!!   io
!!   mpi_morton
!!   paramesh_comm_data
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_morton_order
!!   amr_guardcell
!!   mpi_amr_global_domain_limits
!!   mpi_amr_boundary_block_info
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a checkpoint file in fortran binary 
!!   has been read in.
!!
!!   optional, real, intent(in) :: user_attr1(2,3,4,5)
!!     Arguments which allow the user to add some extra information to the file.
!!     Currently only 5 real numbers can be present.
!!
!! DESCRIPTION
!!
!!  Subroutine to read checkpoint files produced by calling the default 
!!  checkpoint writing routine (amr_checkpoint_wr_default) using PARAMESH.
!!  It read in the tree data structure and data stored in PARAMESH blocks.
!!  Optionally, a user may add a small amout of attribute data to the files written.
!!  This routine excutes the default bevaviour and reads are done serially 
!!  by processor 0.  The data is then sent to other processors.
!!  This routine USES UNFORMATTED DIRECT I/O.  The file which is read in must be in 
!!  fortran binary for the architecture that it is executed on.
!!
!!  The files which are read in must have names of the form 'paramesh_chk_######'.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson (2004).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_checkpoint_re_default (file_num,                  &
                                            l_with_guardcells,         & 
                                            user_attr_1,               &
                                            user_attr_2,               &
                                            user_attr_3,               &
                                            user_attr_4,               &
                                            user_attr_5)
!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_comm_data
      Use paramesh_interfaces, Only : amr_morton_order,                & 
                                      amr_guardcell
      Use paramesh_mpi_interfaces, Only : mpi_amr_global_domain_limits,& 
                                          mpi_amr_boundary_block_info

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output aruments
      Integer, Intent(in) :: file_num
      Logical, Optional, intent(in)  :: l_with_guardcells
      Real, Optional, intent(out) ::                                   & 
       user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5

      Integer :: nguard0, iunit1
      Integer :: block_no
      Integer :: jproc,i,j,ivar,ix,iy,iz,nprocs,mype
      Integer :: lnblockst,ngid
      Integer :: alnblocks,alnblockst
      Integer :: ierr_read,ierr_readt
      Integer :: ierr
      Integer :: tot_blocks

!      Integer :: gid(nfaces+1+nchild,maxblocks_tr)
!      Integer :: lrefinet(maxblocks_tr),nodetypet(maxblocks_tr)
!      Integer :: which_childt(maxblocks_tr)
!      Integer :: gidt(nfaces+1+nchild,maxblocks_tr)
!      Integer :: bflagst(mflags,maxblocks_tr)
      Integer,allocatable :: gid(:,:)
      Integer,allocatable :: lrefinet(:),nodetypet(:)
      Integer,allocatable :: which_childt(:)
      Integer,allocatable :: gidt(:,:)
      Integer,allocatable :: bflagst(:,:)

      Integer :: il0,iu0,jl0,ju0,kl0,ku0
      Integer :: ion_c,ion_f,ion_e,ion_n,iv_c,iv_f,iv_e,iv_n
      Integer,Dimension (:),Allocatable :: glnblocks
      Integer :: isrc,idest,itag,isize,position,position2
      Integer :: ierrorcode
      Integer :: lnblocks_old
      Integer :: status(MPI_STATUS_SIZE)
      Integer :: buf_dim_int 
      Integer :: buf_dim_real
      Integer :: buf_dim1,buf_dim2
      Integer :: no_of_bytes_per_real,no_of_bytes_per_Integer
      Integer :: buf_dim_bytes1,buf_dim_bytes2
      Integer :: nvar_chk_cc,nvar_chk_fc,nvar_chk_ec,nvar_chk_nc
      Integer,Dimension (1) :: ibcast_data

!      Real    :: coordt(mdim,maxblocks_tr)
!      Real    :: work_blockt(maxblocks_tr)
!      Real    :: bnd_boxt(2,mdim,maxblocks_tr)
      Real,allocatable :: coordt(:,:)
      Real,allocatable :: work_blockt(:)
      Real,allocatable :: bnd_boxt(:,:,:)

      Real    :: unkt(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      Real    :: facevarxt(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,      & 
                           kl_bnd:ku_bnd)
      Real    :: facevaryt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,    & 
                           kl_bnd:ku_bnd)
      Real    :: facevarzt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,        & 
                           kl_bnd:ku_bnd+k3d)
      Real    :: unk_nt(nbndvarc,                                      & 
                        il_bnd:iu_bnd+1,                               & 
                        jl_bnd:ju_bnd+k2d,                             & 
                        kl_bnd:ku_bnd+k3d)
      Real    :: unk_e_xt(nbndvare,                                    & 
                          il_bnd:iu_bnd,                               & 
                          jl_bnd:ju_bnd+k2d,                           & 
                          kl_bnd:ku_bnd+k3d)
      Real    :: unk_e_yt(nbndvare,                                    & 
                          il_bnd:iu_bnd+1,                             & 
                          jl_bnd:ju_bnd,                               & 
                          kl_bnd:ku_bnd+k3d)
      Real    :: unk_e_zt(nbndvare,                                    & 
                          il_bnd:iu_bnd+1,                             & 
                          jl_bnd:ju_bnd+k2d,                           & 
                          kl_bnd:ku_bnd)
      Real,Allocatable :: CS_buffer1(:), CR_buffer1(:)
      Real,Allocatable :: CS_buffer2(:), CR_buffer2(:)
      Logical :: l_move_solution, l_with_guardcells2
      Character (len=80) :: filename
      Character (len=6)  :: fnum_string

!-----Begin executable code.
! CEG allocate memory
  allocate(gid(nfaces+1+nchild,maxblocks_tr))
  allocate(lrefinet(maxblocks_tr))
  allocate(nodetypet(maxblocks_tr))
  allocate(which_childt(maxblocks_tr))
  allocate(gidt(nfaces+1+nchild,maxblocks_tr))
  allocate(bflagst(mflags,maxblocks_tr))
  allocate(coordt(mdim,maxblocks_tr))
  allocate(work_blockt(maxblocks_tr))
  allocate(bnd_boxt(2,mdim,maxblocks_tr))

      If (present(l_with_guardcells)) Then
         l_with_guardcells2 = l_with_guardcells
      Else
         l_with_guardcells2 = .False.
      End If

      nguard0 = nguard*npgs
      buf_dim_int = maxblocks*( 3+mflags+(nfaces+1+nchild))
      buf_dim_real =  maxblocks*( 3*mdim + 1 )
      buf_dim1 = buf_dim_real + buf_dim_int
      allocate(CS_buffer1(buf_dim1))
      allocate(CR_buffer1(buf_dim1))

      nvar_chk_cc =  0
      Do i=1,nvar
        If (checkp_on_cc(i)) nvar_chk_cc = nvar_chk_cc + 1
      End Do
      nvar_chk_fc =  0
      Do i=1,nfacevar
        If (checkp_on_fc(1,i)) nvar_chk_fc = nvar_chk_fc + 1
      End Do
      nvar_chk_ec =  0
      Do i=1,nvaredge
        If (checkp_on_ec(1,i)) nvar_chk_ec = nvar_chk_ec + 1
      End Do
      nvar_chk_nc =  0
      Do i=1,nvarcorn
        If (checkp_on_nc(i)) nvar_chk_nc = nvar_chk_nc + 1
      End Do

      buf_dim2 = len_block                                             & 
             + nbndvar*(len_blockfx + len_blockfy + len_blockfz)       & 
             + nbndvarc*len_blockn                                     & 
             + nbndvare*(len_blockex + len_blockey + len_blockez)

      Allocate(CS_buffer2(buf_dim2))
      Allocate(CR_buffer2(buf_dim2))

      Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
      ierr_read = 0
      ierr_readt = 0

      If (allocated(glnblocks)) deallocate( glnblocks )
      Allocate ( glnblocks(0:nprocs-1) )

      Call mpi_type_size(MPI_INTEGER,no_of_bytes_per_Integer,ierr)
      Call mpi_type_size(amr_mpi_real,                                 & 
                         no_of_bytes_per_real   ,ierr)
      buf_dim_bytes1 = buf_dim_real*no_of_bytes_per_real +             & 
                       buf_dim_int *no_of_bytes_per_Integer
      buf_dim_bytes2 = buf_dim2*no_of_bytes_per_real 

!------set limits on data arrays
       il0 = nguard0
       iu0 = nxb-1+nguard0
       jl0 = nguard0*k2d
       ju0 = (nyb-1+nguard0)*k2d
       kl0 = nguard0*k3d
       ku0 = (nzb-1+nguard0)*k3d
       If (.Not.no_permanent_guardcells) Then
       If (l_with_guardcells2) Then
         il0 = 0
         iu0 = nxb+2*nguard0-1
         jl0 = 0
         ju0 = (nyb+2*nguard0-1)*k2d
         kl0 = 0
         ku0 = (nzb+2*nguard0-1)*k3d
       End If
       End If

!-----cell centered data
      iv_c = max(1,nvar)
      ion_c = min(nvar,1)
!-----face-centered data
      iv_f = max(1,nfacevar)
      ion_f = min(nfacevar,1)
!-----edge-centered data
      iv_e = max(1,nvaredge)
      ion_e = min(nvaredge,1)
!-----cell corner data
      iv_n = max(1,nvarcorn)
      ion_n = min(nvarcorn,1)

      iv_c = max(1,nvar_chk_cc)
      ion_c = min(nvar_chk_cc,1)
      iv_f = max(1,nvar_chk_fc)
      ion_f = min(nvar_chk_fc,1)
      iv_e = max(1,nvar_chk_ec)
      ion_e = min(nvar_chk_ec,1)
      iv_n = max(1,nvar_chk_nc)
      ion_n = min(nvar_chk_nc,1)

      If (mype == 0) Then

         Write (fnum_string, '(i6.6)') file_num
         iunit1 = 20
         filename = trim(output_dir) // 'paramesh_chk_' // fnum_string
         Open(unit=iunit1, & 
     &        file=filename, & 
     &        form='unformatted', & 
     &        status='unknown' & 
     &        )

         Read (iunit1) tot_blocks
         Write(*,*) 'blocks to be input ',tot_blocks
         If (present(user_attr_1)) Then
            Read(iunit1) user_attr_1
         End If
         If (present(user_attr_2)) Then
            Read(iunit1) user_attr_2
         End If
         If (present(user_attr_3)) Then
            Read(iunit1) user_attr_3
         End If
         If (present(user_attr_4)) Then
            Read(iunit1) user_attr_4
         End If
         If (present(user_attr_5)) Then
            Read(iunit1) user_attr_5
         End If

!--------compute approximate lnblocks (this will be the number of blocks 
!--------stored on processors 0 -> nprocs-2, nprocs-1 gets tot_blocks - 
!--------the total number on the rest of the blocks)
         alnblocks = int(tot_blocks/nprocs)
         If (tot_blocks < nprocs) Then
            alnblocks = 1
         End If

!--------check for error
         If (tot_blocks-(alnblocks*(nprocs-1)) > maxblocks) Then
           Print *,' ******** ERROR in checkpoint_re: ********'
           Print *,' No. of blocks per processor exceeds maxblocks.'
           Print *,' Suggest you reset maxblocks to a larger number or '
           Print *,' run on a larger no. of processors. '
           ierr_read = 1
           Go To 2
         End If

         Do jproc = 0,nprocs-1
            If (jproc < nprocs-1) Then
               lnblockst = alnblocks
            Else
               lnblockst = tot_blocks - (alnblocks*(nprocs-1))
            End If
!-----------Take care of the case when the no of total blocks is less than
!-----------the number of processors
            If (tot_blocks < nprocs) Then
               If (jproc <= tot_blocks-1) Then
                  lnblockst = 1
               Else
                  lnblockst = 0
               End If
            End If
            glnblocks(jproc) = lnblockst
         End Do

      End If  ! End If (mype == 0) 

      Call MPI_BCAST(glnblocks,nprocs,MPI_INTEGER,0,                   & 
                     MPI_COMM_WORLD,ierr)
      lnblocks = glnblocks(mype)

      If (mype  ==  0) Then
            Do block_no = 1,lnblocks

!--------------Read in data for this block
               Read (iunit1)                                           & 
                    lrefine(block_no),                                 & 
                    nodetype(block_no),                                & 
                    which_child(block_no),                             & 
                    (gid(j,block_no),j=1,nfaces+1+nchild),             & 
                    (bflags(j,block_no),j=1,mflags),                   & 
                    (coord(j,block_no),j=1,ndim),                      & 
                    (bnd_box(1,j,block_no),j=1,ndim),                  & 
                    (bnd_box(2,j,block_no),j=1,ndim),                  & 
                    work_block(block_no)

               ! Read in unk data
               If (nvar_chk_cc > 0) Then
                  Do ivar=1,nvar
                     If (checkp_on_cc(ivar)) Then
                        Do iz=1+kl0*ion_c,1+ku0*ion_c
                           Do iy=1+jl0*ion_c,1+ju0*ion_c
                              Do ix=1+il0*ion_c,1+iu0*ion_c
                                 Read (iunit1) unk(ivar,ix,iy,iz,block_no)
                              End Do
                           End Do
                        End Do
                     Else
                        unk(ivar,:,:,:,block_no) = 0.
                     End If
                  End Do
               End If

               If (nvar_chk_fc > 0) Then
                 Do ivar=1,nfacevar
                 If (checkp_on_fc(1,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   Read (iunit1) facevarx(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 Else
                   facevarx(ivar,:,:,:,block_no) = 0.
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(2,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Read (iunit1) facevary(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 Else
                   facevary(ivar,:,:,:,block_no) = 0.
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(3,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Read (iunit1) facevarz(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 Else
                   facevarz(ivar,:,:,:,block_no) = 0.
                 End If
                 End Do
               End If  ! End If (nvar_chk_fc > 0) 

               If (nvar_chk_ec > 0) Then
                 Do ivar=1,nvaredge
                 If (checkp_on_ec(1,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+iu0*ion_e
                   Read (iunit1) unk_e_x(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_e_x(ivar,:,:,:,block_no) = 0.
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(2,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+ju0*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Read (iunit1) unk_e_y(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_e_y(ivar,:,:,:,block_no) = 0.
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(3,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+ku0*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Read (iunit1) unk_e_z(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_e_z(ivar,:,:,:,block_no) = 0.
                 End If
                 End Do
               End If  ! End If (nvar_chk_ec > 0)

               If (nvar_chk_nc > 0) Then
                 Do ivar=1,nvarcorn
                 If (checkp_on_nc(ivar)) Then
                 Do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 Do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 Do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   Read (iunit1) unk_n(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_n(ivar,:,:,:,block_no) = 0.
                 End If
                 End Do
               End If

            End Do
            If (ndim < 3) Then
              bnd_box(1,3,:) = 0.
              bnd_box(2,3,:) = 1.
              coord(3,:) = .5*(bnd_box(2,3,:)+bnd_box(1,3,:))
            End If
            If (ndim < 2) Then
              bnd_box(1,2,:) = 0.
              bnd_box(2,2,:) = 1.
              coord(2,:) = .5*(bnd_box(2,2,:)+bnd_box(1,2,:))
            End If

      End If  ! End If (mype == 0)

      If (mype  >  0) Then

!--------Post receives on pe mype for messages from proc 0
         lnblockst = glnblocks(mype)
         Do block_no = 1, lnblockst

            isrc = 0
            itag = block_no
            isize = len_block  & 
              + nbndvar*(len_blockfx + len_blockfy + len_blockfz)      & 
              + nbndvare*(len_blockex + len_blockey + len_blockez)     & 
              + nbndvarc*len_blockn 

            Call MPI_RECV(CR_buffer2,isize,amr_mpi_real,               & 
                          isrc,itag,MPI_COMM_WORLD,status,ierr)  

            position2 = 0

            If (nvar > 0) & 
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   unk(1,1,1,1,block_no),len_block,                    & 
                   amr_mpi_real,                                       & 
                   MPI_COMM_WORLD,ierr)
            If (nfacevar > 0) Then
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   facevarx(1,1,1,1,block_no),nbndvar*len_blockfx,     & 
                   amr_mpi_real,MPI_COMM_WORLD,ierr)
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   facevary(1,1,1,1,block_no),nbndvar*len_blockfy,     & 
                   amr_mpi_real,MPI_COMM_WORLD,ierr)
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   facevarz(1,1,1,1,block_no),nbndvar*len_blockfz,     & 
                   amr_mpi_real,MPI_COMM_WORLD,ierr)
            End If
            If (nvaredge > 0) Then
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   unk_e_x(1,1,1,1,block_no),nbndvare*len_blockex,     & 
                   amr_mpi_real,MPI_COMM_WORLD,ierr)
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   unk_e_y(1,1,1,1,block_no),nbndvare*len_blockey,     & 
                   amr_mpi_real,MPI_COMM_WORLD,ierr)
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   unk_e_z(1,1,1,1,block_no),nbndvare*len_blockez,     & 
                   amr_mpi_real,MPI_COMM_WORLD,ierr)
            End If
            If (nvarcorn > 0) Then
              Call MPI_UNPACK(CR_buffer2,                              & 
                   buf_dim_bytes2,position2,                           & 
                   unk_n(1,1,1,1,block_no),nbndvarc*len_blockn,        & 
                   amr_mpi_real,MPI_COMM_WORLD,ierr)
            End If
  
         End Do  ! End Do block_no = 1, lnblockst

         isrc = 0
         idest= mype
         itag = (idest+1)*(maxblocks+1)
         isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild) + 1 )
         Call MPI_RECV(CR_buffer1,isize,amr_mpi_real,                  & 
                       isrc,itag,MPI_COMM_WORLD,status,ierr)

         position = 0
         Do block_no = 1,lnblockst

!-----------fetch data for this block
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              lrefine(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              nodetype(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              which_child(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,      & 
              ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              bflags(1,block_no),mflags,MPI_INTEGER,                   & 
              MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              coord(1,block_no),mdim,amr_mpi_real,                     & 
              MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              bnd_box(1,1,block_no),2*mdim,                            & 
              amr_mpi_real,                                            & 
              MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              work_block(block_no),1,amr_mpi_real,                     & 
              MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
              gid(1,block_no),(nfaces+1+nchild),MPI_INTEGER,           & 
              MPI_COMM_WORLD,ierr)

         End Do  ! End Do block_no = 1,lnblockst

      End If  ! End If (mype > 0)

      If (mype == 0) Then

!-------Read in data and post sends from pe 0 for messages to all other procs
        Do jproc = 1,nprocs-1

            lnblockst = glnblocks(jproc)
            position = 0
            Do block_no = 1,lnblockst

!--------------Read in data for this block
               Read (iunit1)                                           & 
                    lrefinet(block_no),                                & 
                    nodetypet(block_no),                               & 
                    which_childt(block_no),                            & 
                    (gidt(j,block_no),j=1,nfaces+1+nchild),            & 
                    (bflagst(j,block_no),j=1,mflags),                  & 
                    (coordt(j,block_no),j=1,ndim),                     & 
                    (bnd_boxt(1,j,block_no),j=1,ndim),                 & 
                    (bnd_boxt(2,j,block_no),j=1,ndim),                 & 
                    work_blockt(block_no)

               If (nvar_chk_cc > 0) Then
                 Do ivar=1,nvar
                 If (checkp_on_cc(ivar)) Then
                 Do iz=1+kl0*ion_c,1+ku0*ion_c
                 Do iy=1+jl0*ion_c,1+ju0*ion_c
                 Do ix=1+il0*ion_c,1+iu0*ion_c
                   Read (iunit1) unkt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   unkt(ivar,:,:,:) = 0.
                 End If
                 End Do
               End If

               If (nvar_chk_fc > 0) Then
                 Do ivar=1,nfacevar
                 If (checkp_on_fc(1,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   Read (iunit1) facevarxt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   facevarxt(ivar,:,:,:) = 0.
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(2,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Read (iunit1) facevaryt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   facevaryt(ivar,:,:,:) = 0.
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(3,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Read (iunit1) facevarzt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   facevarzt(ivar,:,:,:) = 0.
                 End If
                 End Do
               End If  ! End If (nvar_chk_fc > 0)

               If (nvar_chk_ec > 0) Then
                 Do ivar=1,nvaredge
                 If (checkp_on_ec(1,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+iu0*ion_e
                   Read (iunit1) unk_e_xt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_e_xt(ivar,:,:,:) = 0.
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(2,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+ju0*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Read (iunit1) unk_e_yt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_e_yt(ivar,:,:,:) = 0.
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(3,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+ku0*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Read (iunit1) unk_e_zt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_e_zt(ivar,:,:,:) = 0.
                 End If
                 End Do
               End If  ! End If (nvar_chk_ec > 0)

               If (nvar_chk_nc > 0) Then
                 Do ivar=1,nvarcorn
                 If (checkp_on_nc(ivar)) Then
                 Do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 Do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 Do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   Read (iunit1) unk_nt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 Else
                   unk_nt(ivar,:,:,:) = 0.
                 End If
                 End Do
               End If

               If (ndim < 3) Then
                  bnd_boxt(1,3,:) = 0.
                  bnd_boxt(2,3,:) = 1.
                  coordt(3,:) = .5*(bnd_boxt(2,3,:)+bnd_boxt(1,3,:))
               End If
               If (ndim < 2) Then
                  bnd_boxt(1,2,:) = 0.
                  bnd_boxt(2,2,:) = 1.
                  coordt(2,:) = .5*(bnd_boxt(2,2,:)+bnd_boxt(1,2,:))
               End If

               Call MPI_PACK(lrefinet(block_no),1,MPI_INTEGER,         & 
                 CS_buffer1,buf_dim_bytes1,position,                   & 
                 MPI_COMM_WORLD,ierr)
               Call MPI_PACK(nodetypet(block_no),1,MPI_INTEGER,        & 
                 CS_buffer1,buf_dim_bytes1,position,                   & 
                 MPI_COMM_WORLD,ierr)
               Call MPI_PACK(which_childt(block_no),1,MPI_INTEGER,     & 
                 CS_buffer1,buf_dim_bytes1,position,                   & 
                 MPI_COMM_WORLD,ierr)
               Call MPI_PACK(bflagst(1,block_no),mflags,MPI_INTEGER,   & 
                 CS_buffer1,buf_dim_bytes1,position,                   & 
                 MPI_COMM_WORLD,ierr)
               Call MPI_PACK(coordt(1,block_no),mdim,                  & 
                 amr_mpi_real,                                         & 
                 CS_buffer1,buf_dim_bytes1,position,                   & 
                 MPI_COMM_WORLD,ierr)
               Call MPI_PACK(bnd_boxt(1,1,block_no),2*mdim,            & 
                 amr_mpi_real,                                         & 
                 CS_buffer1,buf_dim_bytes1,position,                   & 
                 MPI_COMM_WORLD,ierr)
               Call MPI_PACK(work_blockt(block_no),1,                  & 
                 amr_mpi_real,                                         & 
                 CS_buffer1,buf_dim_bytes1,position,                   & 
                 MPI_COMM_WORLD,ierr)
               Call MPI_PACK(gidt(1,block_no),(nfaces+1+nchild),       & 
                 MPI_INTEGER,CS_buffer1,buf_dim_bytes1,position,       & 
                 MPI_COMM_WORLD,ierr)

              position2 = 0

              If (nvar > 0) & 
                Call MPI_PACK(unkt(1,1,1,1),                           & 
                  len_block,amr_mpi_real,CS_buffer2,                   & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
              If (nfacevar > 0) Then
                Call MPI_PACK(facevarxt(1,1,1,1),                      & 
                  len_blockfx*nbndvar,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(facevaryt(1,1,1,1),                      & 
                  len_blockfy*nbndvar,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(facevarzt(1,1,1,1),                      & 
                  len_blockfz*nbndvar,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
               End If
              If (nvaredge > 0) Then
                Call MPI_PACK(unk_e_xt(1,1,1,1),                       & 
                  len_blockex*nbndvare,                                & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(unk_e_yt(1,1,1,1),                       & 
                  len_blockey*nbndvare,                                & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(unk_e_zt(1,1,1,1),                       & 
                  len_blockez*nbndvare,                                & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
               End If
              If (nvarcorn > 0) Then
                Call MPI_PACK(unk_nt(1,1,1,1),                         & 
                  len_blockn*nbndvarc,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
               End If

               idest= jproc
               itag = block_no
               isize = len_block                                       & 
                 + nbndvar*(len_blockfx + len_blockfy + len_blockfz)   & 
                 + nbndvare*(len_blockex + len_blockey + len_blockez)  & 
                 + nbndvarc*len_blockn 

               Call MPI_SEND(CS_buffer2,isize,amr_mpi_real,            & 
                             idest,itag,MPI_COMM_WORLD,ierr)

            End Do   ! End Do block_no = 1,lnblockst

            idest= jproc
            itag = (idest+1)*(maxblocks+1)
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild) + 1 )

            Call MPI_SEND(CS_buffer1,isize,amr_mpi_real,               & 
                          idest,itag,MPI_COMM_WORLD,ierr)

        End Do  ! End Do jproc=1,nprocs-1

        Close(iunit1)

      End If  ! End If (mppe == 0)

 2    continue

      Do block_no = 1,lnblocks
        bsize(:,block_no) = bnd_box(2,:,block_no)-bnd_box(1,:,block_no)
      End Do

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!-----all processors fetch error code from proc 0
      ibcast_data(1) = ierr_read
      Call MPI_BCAST(ibcast_data,1,MPI_INTEGER,0,                      & 
                     MPI_COMM_WORLD,ierr)
      ierr_read = ibcast_data(1)
      ierr_readt = ierr_read
      If (ierr_readt == 1)                                             & 
         Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)

!-----COMPUTE TREE DATA FROM gid
      ibcast_data(1) = alnblocks
      Call MPI_BCAST(ibcast_data,1,MPI_INTEGER,0,                      & 
                     MPI_COMM_WORLD,ierr)
      alnblocks = ibcast_data(1)
      alnblockst = alnblocks

      Do block_no = 1,lnblocks
!--------neighbor data
         ngid = 0
         Do j = 1,nfaces
            ngid = ngid + 1
            If (gid(ngid,block_no) > 0) Then
               neigh(2,j,block_no) =                                   & 
                    int((gid(ngid,block_no)-1)/alnblocks)
               If (neigh(2,j,block_no) > nprocs-1)                     & 
                    neigh(2,j,block_no) = nprocs - 1
               neigh(1,j,block_no) = gid(ngid,block_no) -              &  
                    (alnblocks*neigh(2,j,block_no))
            Else
               neigh(1,j,block_no) = gid(ngid,block_no)
               neigh(2,j,block_no) = gid(ngid,block_no)
            End If
         End Do
         
!--------parent data
         ngid = ngid + 1
         If (gid(ngid,block_no) > 0) Then
            parent(2,block_no) =                                       & 
                 int((gid(ngid,block_no)-1)/alnblocks)
            If (parent(2,block_no) > nprocs-1)                         & 
                 parent(2,block_no) = nprocs - 1
            parent(1,block_no) = gid(ngid,block_no) -                  & 
                 (alnblocks*parent(2,block_no))
         Else
            parent(1,block_no) = gid(ngid,block_no)
            parent(2,block_no) = gid(ngid,block_no)
         End If

!--------children data
         Do j = 1,nchild
            ngid = ngid + 1
            If (gid(ngid,block_no) > 0) Then
               child(2,j,block_no) =                                   & 
                    int((gid(ngid,block_no)-1)/alnblocks)
               If (child(2,j,block_no) > nprocs-1)                     & 
                    child(2,j,block_no) = nprocs - 1
               child(1,j,block_no) = gid(ngid,block_no) -              & 
                    (alnblocks*child(2,j,block_no))
            Else
               child(1,j,block_no) = gid(ngid,block_no)
               child(2,j,block_no) = gid(ngid,block_no)
            End If
         End Do
         
      End Do  ! End Do block_no = 1,lnblocks

      If (Allocated(glnblocks))  Deallocate( glnblocks )
      If (Allocated(CS_buffer1)) Deallocate( CS_buffer1 )
      If (Allocated(CR_buffer1)) Deallocate( CR_buffer1 )
      If (Allocated(CS_buffer2)) Deallocate( CS_buffer2 )
      If (Allocated(CR_buffer2)) Deallocate( CR_buffer2 )

!-----Now reorder blocks such that they are better balanced
!-----NOTE: this assumes that the total number of blocks is > nprocs
!-----NOTE: We cannot Do a morton ordering here if l_with_guardcells is 
!-----defined since amr_redist_blks Do not move the guardcells.
      If (.Not.l_with_guardcells2) Then  

      lnblocks_old = lnblocks
      l_move_solution = .True.
      Call amr_morton_order (lnblocks_old,nprocs,mype,                 & 
                             l_move_solution)
!      call ceg_block_balancing_link(mype, nprocs)
      Write(*,*) 'mpi_amr_checkpoint_re_default : after amr_morton_order : pe lnb ',mype,lnblocks

      End If

!-----compute grid_xmax, etc
      Call mpi_amr_global_domain_limits()

      ! CEG included for growing domain cases
      call set_domain_limits_re_hdf5()

      Call amr_morton_process()
!-----Start up an array of cell sizes for each grid refinement level.
!-----These can be used to minimize variation due to rounDoff, but
!-----should ONLY be used with a uniformly spaced grid.
      level_cell_sizes = 0.
      level_cell_sizes(1,1) = (grid_xmax-grid_xmin)/real(nxb)
      If (ndim > 1)                                                    & 
        level_cell_sizes(2,1) = (grid_ymax-grid_ymin)/real(nyb)
      If (ndim == 3)                                                   & 
        level_cell_sizes(3,1) = (grid_zmax-grid_zmin)/real(nzb)
      Do i=2,lrefine_max
        level_cell_sizes(1:ndim,i) = .5*level_cell_sizes(1:ndim,i-1)
      End Do

!-----mark grid as changed
      grid_changed = 1
      grid_analysed_mpi = 1

!-----Now make sure guardcell information is up to date

      If (.Not.l_with_guardcells2)                                     & 
          Call amr_guardcell(mype,1,nguard)

      Call mpi_amr_boundary_block_info(mype,nprocs)

! CEG deallocate memory
  deallocate(gid)
  deallocate(lrefinet)
  deallocate(nodetypet)
  deallocate(which_childt)
  deallocate(gidt)
  deallocate(bflagst)
  deallocate(coordt)
  deallocate(work_blockt)
  deallocate(bnd_boxt)

      Return
      End Subroutine amr_checkpoint_re_default
