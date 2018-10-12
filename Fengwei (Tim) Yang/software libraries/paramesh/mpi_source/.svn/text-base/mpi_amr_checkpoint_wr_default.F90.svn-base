!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_checkpoint_wr_default
!! NAME
!!
!!   amr_checkpoint_wr_default
!!
!! SYNOPSIS
!!
!!   call amr_checkpoint_wr_default(file_num)
!!   call amr_checkpoint_wr_default(file_num, l_with_guardcells, 
!!                                  user_attr1, user_attr2, user_attr3,
!!                                  user_attr4, user_attr5)
!!
!!   call amr_checkpoint_wr_default(integer, optional logical, 
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
!!   io
!!   paramesh_comm_data
!!
!! CALLS
!!
!!   No other PARAMESH routines are called.
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a checkpoint file in fortran binary 
!!   has been written.
!!
!!   optional, real, intent(in) :: user_attr1(2,3,4,5)
!!     Arguments which allow the user to add some extra information to the file.
!!     Currently only 5 real numbers can be added.
!!
!! DESCRIPTION
!!
!!  Subroutine to checkpoint runs using PARAMESH.
!!  It writes out the tree data structure and data stored in PARAMESH blocks.
!!  Optionally, a user may add a small amout of attribute data to the files written.
!!  This routine excutes the default bevaviour and writes are done serially 
!!  by processor 0.  I.e. data is collected from other processors and then written out.
!!  This routine USES UNFORMATTED DIRECT I/O.  The resulting file will be in fortran
!!  binary for the architecture that it is executed on.
!!
!!  The files produced will have names of the form 'paramesh_chk_######'.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson (2004).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_checkpoint_wr_default (file_num,                  &
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

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, Intent(in) :: file_num
      Logical, Optional, Intent(in) :: l_with_guardcells
      Real,    Optional, Intent(in) ::                                 & 
         user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5

!-----Local arrays and variables.
      Integer :: nguard0, iunit1
      Integer :: block_no
      Integer :: jproc,i,j,ivar,ix,iy,iz,nprocs,mype
      Integer :: lnblockst
      Integer :: ngid
      integer,dimension (:), allocatable :: n_to_left

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

      Integer :: ierr
      Integer :: tot_blocks
      Integer :: il0,iu0,jl0,ju0,kl0,ku0
      Integer :: ion_c,ion_f,ion_e,ion_n,iv_c,iv_f,iv_e,iv_n
      Integer,dimension (:),  allocatable :: glnblocks
      Integer :: isrc,idest,itag,isize,ierror,position
      Integer :: status(MPI_STATUS_SIZE)
      Integer :: buf_dim_int
      Integer :: buf_dim_real
      Integer :: buf_dim1, buf_dim2 
      Integer :: no_of_bytes_per_real,no_of_bytes_per_Integer
      Integer :: buf_dim_bytes1,buf_dim_bytes2 
      Integer :: nvar_chk_cc,nvar_chk_fc,nvar_chk_ec,nvar_chk_nc
!      Real    :: coordt(mdim,maxblocks_tr)
!      Real    :: work_blockt(maxblocks_tr)
!      Real    :: bnd_boxt(2,mdim,maxblocks_tr)
      Real,allocatable    :: coordt(:,:)
      Real,allocatable    :: work_blockt(:)
      Real,allocatable    :: bnd_boxt(:,:,:)
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
      Logical :: l_with_guardcells2
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
      buf_dim_int =  maxblocks*( 3+mflags+(nfaces+1+nchild))
      buf_dim_real =  maxblocks*( 3*mdim + 1)
      buf_dim1 = buf_dim_real + buf_dim_int
      Allocate(CS_buffer1(buf_dim1))
      Allocate(CR_buffer1(buf_dim1))
      nvar_chk_cc =  0

      If (nvar > 0) Then
         Do i=1,nvar
            If (checkp_on_cc(i)) nvar_chk_cc = nvar_chk_cc + 1
         End Do
      End If
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

!-----COMPUTE TOTAL NO. OF BLOCKS STORED TO THE 'LEFT' OF THIS PROCESSOR
      If (Allocated(n_to_left)) Deallocate( n_to_left )
      Allocate ( n_to_left(0:nprocs-1) )
      If (Allocated(glnblocks)) Deallocate( glnblocks )
      Allocate ( glnblocks(0:nprocs-1) )

      Call MPI_TYPE_SIZE(MPI_INTEGER,no_of_bytes_per_Integer,ierr)
      Call MPI_TYPE_SIZE(amr_mpi_real,                                 & 
                         no_of_bytes_per_Real   ,ierr)
      buf_dim_bytes1 = buf_dim_real*no_of_bytes_per_Real +             & 
                       buf_dim_int *no_of_bytes_per_Integer
      buf_dim_bytes2 = buf_dim2*no_of_bytes_per_Real

! CEG modified
!      glnblocks(mype) = lnblocks
      Call MPI_AllGATHER(lnblocks, 1,MPI_INTEGER,               & 
                         glnblocks,1,MPI_INTEGER,                      & 
                         MPI_COMM_WORLD,ierror)
      n_to_left = glnblocks

      tot_blocks = 0
      Do i = 0,nprocs-1
         tot_blocks = tot_blocks + n_to_left(i)
      End Do
            
      Do i = nprocs-1,1,-1
         n_to_left(i) = n_to_left(i-1)
      End Do

      n_to_left(0) = 0
      Do i = 2,nprocs-1
         n_to_left(i) = n_to_left(i) + n_to_left(i-1)
      End Do

!-----COMPUTE GLOBAL INDIRECT ADDRESSES FOR TREE DATA (gid)
      Do block_no = 1,lnblocks

         ngid = 0
         Do j = 1,nfaces
            ngid = ngid + 1
            If (neigh(1,j,block_no) > 0) Then
               gid(ngid,block_no) = neigh(1,j,block_no) +              & 
                    n_to_left(neigh(2,j,block_no))
            Else
               gid(ngid,block_no) = neigh(1,j,block_no)
            End If
         End Do
         ngid = ngid + 1
         If (parent(1,block_no) > 0) Then
            gid(ngid,block_no) = parent(1,block_no) +                  & 
                 n_to_left(parent(2,block_no))
         Else
            gid(ngid,block_no) = parent(1,block_no)
         End If
         Do j = 1,nchild
            ngid = ngid + 1
            If (child(1,j,block_no) > 0) Then
               gid(ngid,block_no) = child(1,j,block_no) +              & 
                    n_to_left(child(2,j,block_no))
            Else
               gid(ngid,block_no) = child(1,j,block_no)
            End If
         End Do

      End Do  ! End Do block_no = 1,lnblocks

!-----NOW WRITE OUT THE DATA FROM PROC 0

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
!-----corner data
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

      CR_buffer1 = 0
      CS_buffer1 = 0
      CR_buffer2 = 0
      CS_buffer2 = 0

      If (mype  ==  0) Then

         Write (fnum_string, '(i6.6)') file_num
         iunit1 = 20
         filename = trim(output_dir) // 'paramesh_chk_' // fnum_string
         Open(unit=iunit1,                                             & 
              file=filename,                                           & 
              form='unformatted',                                      & 
              status='unknown'                                         & 
              )

         Write (iunit1) tot_blocks
         If (present(user_attr_1)) Then
            Write(iunit1) user_attr_1
         End If
         If (present(user_attr_2)) Then
            Write(iunit1) user_attr_2
         End If
         If (present(user_attr_3)) Then
            Write(iunit1) user_attr_3
         End If
         If (present(user_attr_4)) Then
            Write(iunit1) user_attr_4
         End If
         If (present(user_attr_5)) Then
            Write(iunit1) user_attr_5
         End If

!--------Write data from processor 0
         Do block_no = 1,lnblocks

               Write (iunit1)                                          & 
                    lrefine(block_no),                                 & 
                    nodetype(block_no),                                & 
                    which_child(block_no),                             & 
                    (gid(j,block_no),j=1,nfaces+1+nchild),             & 
                    (bflags(j,block_no),j=1,mflags),                   & 
                    (coord(j,block_no),j=1,ndim),                      & 
                    (bnd_box(1,j,block_no),j=1,ndim),                  & 
                    (bnd_box(2,j,block_no),j=1,ndim),                  & 
                    work_block(block_no)

               If (nvar_chk_cc > 0) Then
                 Do ivar=1,nvar
                 If (checkp_on_cc(ivar)) Then
                 Do iz=1+kl0*ion_c,1+ku0*ion_c
                 Do iy=1+jl0*ion_c,1+ju0*ion_c
                 Do ix=1+il0*ion_c,1+iu0*ion_c
                   Write (iunit1) unk(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do
               End If

               If (nvar_chk_fc > 0) Then
                 Do ivar=1,nfacevar
                 If (checkp_on_fc(1,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   Write (iunit1) facevarx(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(2,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Write (iunit1) facevary(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(3,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Write (iunit1) facevarz(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

               End If  ! End If (nvar_chk_fc > 0)

               If (nvar_chk_ec > 0) Then
                 Do ivar=1,nvaredge
                 If (checkp_on_ec(1,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+iu0*ion_e
                   Write (iunit1) unk_e_x(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(2,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+ju0*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Write (iunit1) unk_e_y(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(3,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+ku0*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Write (iunit1) unk_e_z(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

               End If  ! End If (nvar_chk_ec > 0)

               If (nvar_chk_nc > 0) Then
                 Do ivar=1,nvarcorn
                 If (checkp_on_nc(ivar)) Then
                 Do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 Do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 Do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   Write (iunit1) unk_n(ivar,ix,iy,iz,block_no)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do
               End If

         End Do  ! End Do block_no = 1,lnblocks

! Now loop over processors and 1) fetch data to proc. 0 and then 2) write the
! data from proc. 0
         Do jproc = 1,nprocs-1

! Post receives on pe 0 for messages from all other procs
! fetch lnblocks from other processors
            lnblockst = glnblocks(jproc)
            isrc = jproc
            idest= 0
            itag = (isrc+1)*(maxblocks+1)
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)      &   
                      + 1) 

            Call MPI_RECV(CR_buffer1,isize,amr_mpi_real,               & 
                          isrc,itag,MPI_COMM_WORLD,status,ierr)

            position = 0
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
               lrefinet,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
                nodetypet,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
                which_childt,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
                bflagst,lnblockst*mflags,MPI_INTEGER,                  & 
                MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
                coordt(1,1),mdim*lnblockst,amr_mpi_real,               & 
                MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
                bnd_boxt(1,1,1),2*mdim*lnblockst,                      & 
                amr_mpi_real,                                          & 
                MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
                work_blockt(1),lnblockst,                              & 
                amr_mpi_real,                                          & 
                MPI_COMM_WORLD,ierr)
            Call MPI_UNPACK(CR_buffer1,buf_dim_bytes1,position,        & 
                gidt(1,1),lnblockst*(nfaces+1+nchild),MPI_INTEGER,     & 
                MPI_COMM_WORLD,ierr)

            Do block_no = 1,lnblockst

               position = 0
               isize = len_block                                       & 
                 + nbndvar*(len_blockfx + len_blockfy + len_blockfz)   & 
                 + nbndvare*(len_blockex + len_blockey + len_blockez)  &  
                 + nbndvarc*len_blockn 
               itag = block_no

!--------------fetch data for this block (= block_no)
               Call MPI_RECV(CR_buffer2,isize,amr_mpi_real,            & 
                             isrc,itag,MPI_COMM_WORLD,status,ierr)

               If (nvar > 0) Then
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   unkt(1,1,1,1),len_block,amr_mpi_real,               & 
                   MPI_COMM_WORLD,ierr)
               End If
               If (nfacevar > 0) Then
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   facevarxt(1,1,1,1),nbndvar*len_blockfx,             & 
                   amr_mpi_real,                                       & 
                   MPI_COMM_WORLD,ierr)
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   facevaryt(1,1,1,1),nbndvar*len_blockfy,             & 
                   amr_mpi_real,                                       & 
                   MPI_COMM_WORLD,ierr)
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   facevarzt(1,1,1,1),nbndvar*len_blockfz,             & 
                   amr_mpi_real, & 
                   MPI_COMM_WORLD,ierr)
               End If
               If (nvaredge > 0) Then
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   unk_e_xt(1,1,1,1),nbndvare*len_blockex,             & 
                   amr_mpi_real,                                       & 
                   MPI_COMM_WORLD,ierr)
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   unk_e_yt(1,1,1,1),nbndvare*len_blockey,             & 
                   amr_mpi_real,                                       & 
                   MPI_COMM_WORLD,ierr)
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   unk_e_zt(1,1,1,1),nbndvare*len_blockez,             & 
                   amr_mpi_real,                                       & 
                   MPI_COMM_WORLD,ierr)
               End If
               If (nvarcorn > 0) Then
                 Call MPI_UNPACK(CR_buffer2,                           & 
                   buf_dim_bytes2,position,                            & 
                   unk_nt(1,1,1,1),nbndvarc*len_blockn,                & 
                   amr_mpi_real,                                       & 
                   MPI_COMM_WORLD,ierr)
               End If

               Write (iunit1)                                          & 
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
                   Write (iunit1) unkt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do
               End If

               If (nvar_chk_fc > 0) Then
                 Do ivar=1,nfacevar
                 If (checkp_on_fc(1,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   Write (iunit1) facevarxt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(2,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+ku0*ion_f
                 Do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Write (iunit1) facevaryt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nfacevar
                 If (checkp_on_fc(3,ivar)) Then
                 Do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 Do iy = 1+jl0*ion_f,1+ju0*ion_f
                 Do ix = 1+il0*ion_f,1+iu0*ion_f
                   Write (iunit1) facevarzt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

               End If  ! End If (nvar_chk_fc > 0)

               If (nvar_chk_ec > 0) Then
                 Do ivar=1,nvaredge
                 If (checkp_on_ec(1,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+iu0*ion_e
                   Write (iunit1) unk_e_xt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(2,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 Do iy = 1+jl0*ion_e,1+ju0*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Write (iunit1) unk_e_yt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

                 Do ivar=1,nvaredge
                 If (checkp_on_ec(3,ivar)) Then
                 Do iz = 1+kl0*ion_e,1+ku0*ion_e
                 Do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 Do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   Write (iunit1) unk_e_zt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do

               End If  ! End If (nvar_chk_ec > 0)

               If (nvar_chk_nc > 0) Then
                 Do ivar=1,nvarcorn
                 If (checkp_on_nc(ivar)) Then
                 Do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 Do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 Do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   Write (iunit1) unk_nt(ivar,ix,iy,iz)
                 End Do
                 End Do
                 End Do
                 End If
                 End Do
               End If

            End Do  ! End Do block_no = 1,lnblockst
         End Do  ! End Do jproc = 1,nprocs-1

         Close(iunit1)

       End If  ! End If (mype == 0)

       If (mype > 0) Then

!-----------Post sends to pe 0 for messages from all other procs
!-----------fetch lnblocks from other processors
            lnblockst = lnblocks
            isrc = mype
            idest= 0
            itag = (isrc+1)*(maxblocks+1)
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)      & 
                                + 1 )
            position = 0
            Call MPI_PACK(lrefine(1),lnblocks,MPI_INTEGER,             & 
                CS_buffer1, & 
                buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            Call MPI_PACK(nodetype(1),lnblocks,MPI_INTEGER,            & 
                CS_buffer1, & 
                buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            Call MPI_PACK(which_child(1),lnblocks,MPI_INTEGER,         & 
                CS_buffer1, & 
                buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            Call MPI_PACK(bflags(1,1),lnblocks*mflags,MPI_INTEGER,     & 
                CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            Call MPI_PACK(coord(1,1),mdim*lnblocks,                    & 
                amr_mpi_real,                                          & 
                CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            Call MPI_PACK(bnd_box(1,1,1),2*mdim*lnblocks,              & 
                amr_mpi_real,                                          & 
                CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            Call MPI_PACK(work_block(1),lnblocks,                      & 
                amr_mpi_real,                                          & 
                CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            Call MPI_PACK(gid(1,1),lnblocks*(nfaces+1+nchild),         & 
                MPI_INTEGER,CS_buffer1,buf_dim_bytes1,position,        & 
                MPI_COMM_WORLD,ierr)

            Call MPI_SEND(CS_buffer1,isize,amr_mpi_real,               & 
             idest,itag,MPI_COMM_WORLD,ierr)

            Do block_no= 1,lnblocks

              position = 0

              If (nvar > 0) & 
                Call MPI_PACK(unk(1,1,1,1,block_no),                   & 
                  len_block,amr_mpi_real,CS_buffer2,                   & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)

              If (nfacevar > 0) Then
                Call MPI_PACK(facevarx(1,1,1,1,block_no),              & 
                  len_blockfx*nbndvar,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(facevary(1,1,1,1,block_no),              & 
                  len_blockfy*nbndvar,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(facevarz(1,1,1,1,block_no),              & 
                  len_blockfz*nbndvar,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
               End If
              If (nvaredge > 0) Then
                Call MPI_PACK(unk_e_x(1,1,1,1,block_no),               & 
                  len_blockex*nbndvare,                                & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(unk_e_y(1,1,1,1,block_no),               & 
                  len_blockey*nbndvare,                                & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                Call MPI_PACK(unk_e_z(1,1,1,1,block_no),               & 
                  len_blockez*nbndvare,                                & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
               End If
               If (nvarcorn > 0) Then
                Call MPI_PACK(unk_n(1,1,1,1,block_no),                 & 
                  len_blockn*nbndvarc,                                 & 
                  amr_mpi_real,CS_buffer2,                             & 
                  buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
               End If

               itag = block_no
               isize = len_block                                       & 
                 + nbndvar*(len_blockfx + len_blockfy + len_blockfz)   & 
                 + nbndvare*(len_blockex + len_blockey + len_blockez)  & 
                 + nbndvarc*len_blockn 

               Call MPI_SEND(CS_buffer2,isize,amr_mpi_real,            & 
                             idest,itag,MPI_COMM_WORLD,ierr)

            End Do  ! End Do block_no= 1,lnblocks

      End If ! End If (mype > 0)
 
      If (Allocated(n_to_left))  Deallocate( n_to_left )
      If (Allocated(glnblocks))  Deallocate( glnblocks )
      If (Allocated(CS_buffer1)) Deallocate( CS_buffer1 )
      If (Allocated(CR_buffer1)) Deallocate( CR_buffer1 )
      If (Allocated(CS_buffer2)) Deallocate( CS_buffer2 )
      If (Allocated(CR_buffer2)) Deallocate( CR_buffer2 )

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
      End Subroutine amr_checkpoint_wr_default



