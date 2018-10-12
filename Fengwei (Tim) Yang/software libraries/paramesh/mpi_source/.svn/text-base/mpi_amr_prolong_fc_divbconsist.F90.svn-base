!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****fi mpi_source/amr_prolong_fc_divbconsist
!! NAME
!!
!!   amr_prolong_fc_divbconsist
!! 
!! SYNOPSIS
!!
!!   call amr_prolong_fc_divbconsist (mype, level, nfield)
!!
!!   call amr_prolong_fc_divbconsist (integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: mype     
!!      Current processor number
!!
!!   integer, intent(in) :: level
!!      level to prolong from.
!!
!!   integer, intent(in) :: nfield
!!     selects which field to prolong (1 - facevarx, 2 - facevary, 3 - facevarz)
!!
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!    paramesh_dimensions
!!    physicaldata
!!    tree
!!    mpi_morton
!!    paramesh_mpi_interfaces
!!
!! CALLS
!!
!!    mpi_amr_get_remote_block_fvar
!!    compute_evalues
!!
!! RETURNS
!!
!!   Upon exit prolongation of the face variable is performed and is
!!   guaranteed to be diveregenge free.
!!
!! DESCRIPTION
!!
!!   Like the routine amr_prolong_fc_consist this routine checks for
!!   existing neighbor to newly created child blocks, and where found,
!!   uses facevar data from the existing neighbor at the common block
!!   boundary, in place of interpolation from the new childs parent.
!!   In addition this routine makes this change while adjusting values
!!   immediately interior to the new face to guarantee that div B is
!!   kept at zero.
!!
!!   Only 1 field constructed from facevar components can be modified 
!!   during this call.
!!   To select it set the values for nfield and i_divf_fc_vars.
!!   Note, the algorithm which modifies the field is not optimal. In particular
!!   it has a tendency to drive oscillations in the neighborhood of strong
!!   field gradients. It is our intenetion to improve this in the future.
!!
!!   This routine is called from amr_prolong.
!!
!! AUTHORS
!!
!!   Peter MacNeice (April 1998)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_prolong_fc_divbconsist(mype,level,nfield)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use mpi_morton
      Use paramesh_mpi_interfaces, only :                              & 
                             mpi_amr_get_remote_block_fvar

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      Integer, intent(in) ::  mype
      Integer, intent(in) ::  level
      Integer, intent(in) ::  nfield

!-----Local arrays and variables

      Integer :: remote_pe,remote_block
      Integer :: ierrorcode,ierr
      Integer :: isw,idvx,idvy,idvz,isg,jf,i,j,k
      Integer :: idest,i_dest,j_dest,k_dest,i_source,j_source,k_source
      Integer :: nguard0
      Integer :: nd, iblk
      Real :: recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,             & 
                    kl_bnd:ku_bnd)
      Real :: recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           & 
                    kl_bnd:ku_bnd)
      Real :: recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,               & 
                    kl_bnd:ku_bnd+k3d)
      Real :: tempx(jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      Real :: tempy(il_bnd:iu_bnd,kl_bnd:ku_bnd)
      Real :: tempz(il_bnd:iu_bnd,jl_bnd:ju_bnd)
      Real :: e1,e2,ea,eb,efact
      Real :: b11,b12,b21,b22
      Real :: dx,dy,dz
      Real :: area11,area12,area21,area22
      Real :: areax1,areax2,areay1,areay2,areaz1,areaz2
      Real :: bsum,divbmax,divb
      Logical :: cnewchild

!-----Begin executable code.

      nguard0 = nguard*npgs
      nd = nguard - nguard0

      If (ndim > 1) Then

      print *,'CEG mpi_amr_prolong_fc_divbconsist.F90'
! CEG removed
!      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!-----Set components of facevar to which div correction is to be applied.
      idvx = i_divf_fc_vars(1,nfield)
      idvy = i_divf_fc_vars(2,nfield)
      idvz = i_divf_fc_vars(3,nfield)

      isw = 0

      efact = 1.
      If (ndim == 2) efact = .5

!-----cycle through the grid blocks on this processor
      If (lnblocks > 0) Then
      Do isg = 1,lnblocks

!-----Is this a newly created leaf block ?
      If (nodetype(isg) == 1.And.newchild(isg)                         & 
         .And. lrefine(isg) == level) Then

!------get block geometry information
       If (curvilinear) Then
       call amr_block_geometry(isg,mype)
       Else
       dx = (bnd_box(2,1,isg)-bnd_box(1,1,isg))/real(nxb)
       dy = 1.
       dz = 1.
       If (ndim >= 2)                                                  & 
               dy = (bnd_box(2,2,isg)-bnd_box(1,2,isg))/real(nyb)
       If (ndim == 3)                                                  & 
               dz = (bnd_box(2,3,isg)-bnd_box(1,3,isg))/real(nzb)
       End If

!------Cycle over the blocks faces
       Do jf = 1,nfaces

          remote_pe = neigh(2,jf,isg)
          remote_block  = neigh(1,jf,isg)

          If (remote_block > 0.And.remote_pe.ne.mype) Then
            Do iblk = strt_buffer,last_buffer
              If ( (remote_pe == laddress(2,iblk)) .And.               & 
                  (remote_block == laddress(1,iblk)) ) Then
                remote_pe = mype
                remote_block = iblk
              End If
            End Do
          End If
          If (remote_block > 0.And.remote_pe.ne.mype) Then
            write(*,*) 'Error : amr_prolong_fc_divbconsist : ',        & 
             ' pe ',mype,' current blk ',isg,                          & 
             ' Remote block ',remote_block,remote_pe,' not ',          & 
             'found loCally.'
            Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
          End If

!---------Is the neighbor to this face a pre-existing block?
          cnewchild = .True.
          If (remote_block > 0) Then
              cnewchild = newchild(remote_block)
          End If
          
          If (.Not.cnewchild) Then

!-----------If the neighbor block is pre-existing Then get its facevar data 
!-----------on the shared block boundary.
            idest = isg

            If (jf == 1) Then
               i_dest   = nguard0 + 1 + iface_off
               i_source = nxb+nguard0 + 1 -gc_off_x + iface_off
               tempx(:,:) = facevarx(idvx,i_dest,:,:,idest)

               Call mpi_amr_get_remote_block_fvar(mype,                & 
                           remote_pe,remote_block,1,                   & 
                           recvx,recvy,recvz,idest)
               facevarx(idvx,i_dest,:,:,idest) =                       & 
                                   recvx(idvx,i_source,:,:)

               Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               Do j=nguard0+1,nguard0 + nyb-1,2

                If (curvilinear) Then
                area11= cell_area1(1+nguard,j+nd*k2d,k+nd*k3d)
                area21= cell_area1(1+nguard,j+nd*k2d,k+(nd+1)*k3d)
                area12= cell_area1(1+nguard,j+(nd+1)*k2d,k+nd*k3d)
                area22= cell_area1(1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d)
                Else
                area11 = dy*dz
                area21 = dy*dz
                area12 = dy*dz
                area22 = dy*dz
                End If

                b11 = (facevarx(idvx,nguard0+1,j,k,idest)              & 
                             -tempx(j,k))*area11
                b21 = (facevarx(idvx,nguard0+1,j,k+k3d,idest)          & 
                             - tempx(j,k+k3d))*area21
                b12 = (facevarx(idvx,nguard0+1,j+1,k,idest)            & 
                             - tempx(j+1,k))*area12
                b22 = (facevarx(idvx,nguard0+1,j+1,k+k3d,idest)        & 
                             - tempx(j+1,k+k3d))*area22
                bsum = b11+b12+b21+b22

                Call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                If (curvilinear) Then
                areay1= cell_area2(1+nguard,j+(nd+1)*k2d,k+nd*k3d) 
                areay2= cell_area2(1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d) 
                areaz1= cell_area3(1+nguard,j+k2d,k+(nd+1)*k3d)
                areaz2= cell_area3(1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d)
                Else
                areay1= dx*dz
                areay2= dx*dz
                areaz1= dx*dy
                areaz2= dx*dy
                End If

                 facevary(idvy,nguard0+1,j+1,k,idest) =                & 
                   facevary(idvy,nguard0+1,j+1,k,idest)                & 
                   - e1*efact/areay1
                 facevary(idvy,nguard0+1,j+1,k+k3d,idest) =            & 
                   facevary(idvy,nguard0+1,j+1,k+k3d,idest)            & 
                   - e2*efact/areay2
                 facevarz(idvz,nguard0+1,j,k+k3d,idest) =              & 
                   facevarz(idvz,nguard0+1,j,k+k3d,idest)              & 
                   + ea/areaz1
                 facevarz(idvz,nguard0+1,j+1,k+k3d,idest) =            & 
                   facevarz(idvz,nguard0+1,j+1,k+k3d,idest)            & 
                   + eb/areaz2

               End Do  ! End Do j=nguard0+1,nguard0 + nyb-1,2
               End Do  ! End Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2

            ElseIf (jf == 2) Then

               i_dest   = nxb+1+nguard0 + iface_off
               i_source = 1+nguard0+gc_off_x + iface_off
               tempx(:,:) = facevarx(idvx,i_dest,:,:,idest)
               Call mpi_amr_get_remote_block_fvar(mype,                & 
                           remote_pe,remote_block,1,                   & 
                           recvx,recvy,recvz,idest)
               facevarx(idvx,i_dest,:,:,idest) =                       & 
                                   recvx(idvx,i_source,:,:)
 
               Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               Do j=nguard0+1,nguard0 + nyb-1,2

                If (curvilinear) Then
                area11= cell_area1(nxb+1+nguard,j+nd*k2d,k+nd*k3d)
                area21= cell_area1(nxb+1+nguard,j+nd*k2d,k+(nd+1)*k3d)
                area12= cell_area1(nxb+1+nguard,j+(nd+1)*k2d,k+nd*k3d)
                area22=                                                & 
                    cell_area1(nxb+1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d)
                Else
                area11 = dy*dz
                area21 = dy*dz
                area12 = dy*dz
                area22 = dy*dz
                End If

                b11 = (facevarx(idvx,nguard0+nxb+1,j,k,idest)          & 
                               - tempx(j,k))*area11
                b21 = (facevarx(idvx,nguard0+nxb+1,j,k+k3d,idest)      & 
                               - tempx(j,k+k3d))*area21
                b12 = (facevarx(idvx,nguard0+nxb+1,j+1,k,idest)        & 
                               - tempx(j+1,k))*area12
                b22 = (facevarx(idvx,nguard0+nxb+1,j+1,k+k3d,idest)    & 
                               - tempx(j+1,k+k3d))*area22
                bsum = b11+b12+b21+b22
                Call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                If (curvilinear) Then
                areay1= cell_area2(nxb+nguard,j+(nd+1)*k2d,k+nd*k3d)
                areay2= cell_area2(nxb+nguard,j+(nd+1)*k2d,k+nd*k3d)
                areaz1= cell_area3(nxb+nguard,j+k2d,k+(nd+1)*k3d)
                areaz2= cell_area3(nxb+nguard,j+nd*k2d,k+(nd+1)*k3d)
                Else
                areay1= dx*dz
                areay2= dx*dz
                areaz1= dx*dy
                areaz2= dx*dy
                End If

                facevary(idvy,nguard0+nxb,j+1,k,idest) =               & 
                   facevary(idvy,nguard0+nxb,j+1,k,idest)              & 
                   + e1*efact/areay1
                facevary(idvy,nguard0+nxb,j+1,k+k3d,idest) =           & 
                   facevary(idvy,nguard0+nxb,j+1,k+k3d,idest)          & 
                   + e2*efact/areay2
                facevarz(idvz,nguard0+nxb,j,k+k3d,idest) =             & 
                   facevarz(idvz,nguard0+nxb,j,k+k3d,idest)            & 
                   - ea/areaz1
                facevarz(idvz,nguard0+nxb,j+1,k+k3d,idest) =           & 
                   facevarz(idvz,nguard0+nxb,j+1,k+k3d,idest)          & 
                   - eb/areaz2
               End Do  ! End Do j=nguard0+1,nguard0 + nyb-1,2
               End Do  ! End Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2

            ElseIf (jf == 3) Then

               j_dest   = nguard0*k2d + 1 + iface_off*k2d
               j_source = nyb+nguard0 + 1 -gc_off_y + iface_off
               tempy(:,:) = facevary(idvy,:,j_dest,:,idest)
               Call mpi_amr_get_remote_block_fvar(mype,                & 
                           remote_pe,remote_block,2,                   & 
                           recvx,recvy,recvz,idest)
               facevary(idvy,:,j_dest,:,idest) =                       & 
                                      recvy(idvy,:,j_source,:)

               Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               Do i=nguard0+1,nguard0 + nxb-1,2

                If (curvilinear) Then
                area11= cell_area2(i+nd,1+nguard,k+nd*k3d)
                area21= cell_area2(i+nd,1+nguard,k+(nd+1)*k3d)
                area12= cell_area2(i+1+nd,1+nguard,k+nd*k3d)
                area22= cell_area2(i+1+nd,1+nguard,k+(nd+1)*k3d)
                Else
                area11 = dx*dz
                area21 = dx*dz
                area12 = dx*dz
                area22 = dx*dz
                End If

                b11 = (facevary(idvy,i,nguard0*k2d+1,k,idest)          & 
                       - tempy(i,k) )*area11
                b21 = (facevary(idvy,i,nguard0*k2d+1,k+k3d,idest)      & 
                       - tempy(i,k+k3d) )*area21
                b12 = (facevary(idvy,i+1,nguard0*k2d+1,k,idest)        & 
                       - tempy(i+1,k) )*area12
                b22 = (facevary(idvy,i+1,nguard0*k2d+1,k+k3d,idest)    & 
                       - tempy(i+1,k+k3d) )*area22
                bsum = b11+b12+b21+b22
                Call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                If (curvilinear) Then
                areax1= cell_area1(i+1+nd,1+nguard,k+nd*k3d)
                areax2= cell_area1(i+1+nd,1+nguard,k+nd*k3d)
                areaz1= cell_area3(i+nd,1+nguard,k+(nd+1)*k3d)
                areaz2= cell_area3(i+nd,1+nguard,k+(nd+1)*k3d)
                Else
                areax1= dy*dz
                areax2= dy*dz
                areaz1= dx*dy
                areaz2= dx*dy
                End If

                facevarx(idvx,i+1,nguard0*k2d+1,k,idest) =             & 
                   facevarx(idvx,i+1,nguard0*k2d+1,k,idest)            & 
                   - e1*efact/areax1
                facevarx(idvx,i+1,nguard0*k2d+1,k+k3d,idest) =         & 
                   facevarx(idvx,i+1,nguard0*k2d+1,k+k3d,idest)        & 
                   - e2*efact/areax2
                facevarz(idvz,i,nguard0*k2d+1,k+k3d,idest) =           & 
                   facevarz(idvz,i,nguard0*k2d+1,k+k3d,idest)          & 
                   + ea/areaz1
                facevarz(idvz,i+1,nguard0*k2d+1,k+k3d,idest) =         & 
                   facevarz(idvz,i+1,nguard0*k2d+1,k+k3d,idest)        & 
                   + eb/areaz2
               End Do  ! End Do i=nguard0+1,nguard0 + nxb-1,2
               End Do  ! End Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2

            ElseIf (jf == 4) Then

               j_dest   = nyb*k2d + 1 + nguard0*k2d + iface_off*k2d
               j_source = 1+nguard0+gc_off_y + iface_off
               tempy(:,:) = facevary(idvy,:,j_dest,:,idest)
               Call mpi_amr_get_remote_block_fvar(mype,                & 
                           remote_pe,remote_block,2,                   & 
                           recvx,recvy,recvz,idest)
               facevary(idvy,:,j_dest,:,idest) =                       & 
                                 recvy(idvy,:,j_source,:)

               Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               Do i=nguard0+1,nguard0 + nxb-1,2

                If (curvilinear) Then
                area11= cell_area2(i+nd,nyb+1+nguard,k+nd*k3d)
                area21= cell_area2(i+nd,nyb+1+nguard,k+(nd+1)*k3d)
                area12= cell_area2(i+1+nd,nyb+1+nguard,k+nd*k3d)
                area22= cell_area2(i+1+nd,nyb+1+nguard,k+(nd+1)*k3d)
                Else
                area11 = dx*dz
                area21 = dx*dz
                area12 = dx*dz
                area22 = dx*dz
                End If

                b11 = ( facevary(idvy,i,nguard0*k2d+nyb+k2d,k,idest)   & 
                            - tempy(i,k) )*area11
                b21 = (facevary(idvy,i,nguard0*k2d+nyb+k2d,k+k3d,idest)& 
                            - tempy(i,k+k3d) )*area21
                b12 = ( facevary(idvy,i+1,nguard0*k2d+nyb+k2d,k,idest) & 
                            - tempy(i+1,k) )*area12
                b22 = ( facevary(idvy,i+1,nguard0*k2d+nyb+k2d,         & 
                              k+k3d,idest) & 
                            - tempy(i+1,k+k3d) )*area22
                bsum = b11+b12+b21+b22
                Call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                If (curvilinear) Then
                areax1= cell_area1(i+1+nd,nyb+1+nguard,k+nd*k3d)
                areax2= cell_area1(i+1+nd,nyb+1+nguard,k+nd*k3d)
                areaz1= cell_area3(i+nd,nyb+1+nguard,k+(nd+1)*k3d)
                areaz2= cell_area3(i+nd,nyb+1+nguard,k+(nd+1)*k3d)
                Else
                areax1= dy*dz
                areax2= dy*dz
                areaz1= dx*dy
                areaz2= dx*dy
                End If

                facevarx(idvx,i+1,nguard0*k2d+nyb,k,idest) =           & 
                   facevarx(idvx,i+1,nguard0*k2d+nyb,k,idest)          & 
                        + e1*efact/areax1
                facevarx(idvx,i+1,nguard0*k2d+nyb,k+k3d,idest) =       & 
                   facevarx(idvx,i+1,nguard0*k2d+nyb,k+k3d,idest)      & 
                        + e2*efact/areax2
                facevarz(idvz,i,nguard0*k2d+nyb,k+k3d,idest) =         & 
                   facevarz(idvz,i,nguard0*k2d+nyb,k+k3d,idest)        & 
                        - ea/areaz1
                facevarz(idvz,i+1,nguard0*k2d+nyb,k+k3d,idest) =       & 
                   facevarz(idvz,i+1,nguard0*k2d+nyb,k+k3d,idest)      & 
                        - eb/areaz2
               End Do  ! End Do i=nguard0+1,nguard0 + nxb-1,2
               End Do  ! End Do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2

            ElseIf (jf == 5 .And. ndim == 3) Then

               k_dest   = nguard0*k3d + 1 + iface_off*k3d
               k_source = nzb+nguard0 + 1 -gc_off_z + iface_off
               tempz(:,:) = facevarz(idvz,:,:,k_dest,idest)
               Call mpi_amr_get_remote_block_fvar(mype,                & 
                           remote_pe,remote_block,3,                   & 
                           recvx,recvy,recvz,idest)
               facevarz(idvz,:,:,k_dest,idest) =                       & 
                                  recvz(idvz,:,:,k_source)

               Do j=nguard0+1,nguard0 + nyb-1,2
               Do i=nguard0+1,nguard0 + nxb-1,2

                If (curvilinear) Then
                area11= cell_area3(i+nd  ,j+nd  ,1+nguard*k3d)
                area21= cell_area3(i+nd  ,j+1+nd,1+nguard*k3d)
                area12= cell_area3(i+1+nd,j+nd  ,1+nguard*k3d)
                area22= cell_area3(i+1+nd,j+1+nd,1+nguard*k3d)
                Else
                area11 = dx*dy
                area21 = dx*dy
                area12 = dx*dy
                area22 = dx*dy
                End If

                b11 = (facevarz(idvz,i,j,nguard0*k3d+1,idest)          & 
                        - tempz(i,j) )*area11
                b21 = (facevarz(idvz,i,j+1,nguard0*k3d+1,idest)        & 
                        - tempz(i,j+1) )*area21
                b12 = (facevarz(idvz,i+1,j,nguard0*k3d+1,idest)        & 
                        - tempz(i+1,j) )*area12
                b22 = (facevarz(idvz,i+1,j+1,nguard0*k3d+1,idest)      & 
                        - tempz(i+1,j+1) )*area22
                bsum = b11+b12+b21+b22
                Call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                If (curvilinear) Then
                areax1= cell_area1(i+1+nd,j+nd  ,1+nguard*k3d)
                areax2= cell_area1(i+1+nd,j+1+nd,1+nguard*k3d)
                areay1= cell_area2(i+nd  ,j+1+nd,1+nguard*k3d)
                areay2= cell_area2(i+nd+1,j+1+nd,1+nguard*k3d)
                Else
                areax1= dy*dz
                areax2= dy*dz
                areay1= dx*dz
                areay2= dx*dz
                End If

                facevarx(idvx,i+1,j,nguard0*k3d+1,idest) =             & 
                   facevarx(idvx,i+1,j,nguard0*k3d+1,idest)            & 
                   - e1/areax1
                facevarx(idvx,i+1,j+1,nguard0*k3d+1,idest) =           & 
                   facevarx(idvx,i+1,j+1,nguard0*k3d+1,idest)          & 
                   - e2/areax2
                facevary(idvy,i,j+1,nguard0*k3d+1,idest) =             & 
                   facevary(idvy,i,j+1,nguard0*k3d+1,idest)            & 
                   + ea/areay1
                facevary(idvy,i+1,j+1,nguard0*k3d+1,idest) =           & 
                   facevary(idvy,i+1,j+1,nguard0*k3d+1,idest)          & 
                   + eb/areay2
               End Do  ! End Do i=nguard0+1,nguard0 + nxb-1,2
               End Do  ! End Do j=nguard0+1,nguard0 + nyb-1,2

            ElseIf (jf == 6 .And. ndim == 3) Then

               k_dest   = nzb*k3d+1+nguard0*k3d + iface_off*k3d
               k_source = 1+nguard0+gc_off_z + iface_off
               tempz(:,:) = facevarz(idvz,:,:,k_dest,idest)
               Call mpi_amr_get_remote_block_fvar(mype,                & 
                           remote_pe,remote_block,3,                   & 
                           recvx,recvy,recvz,idest)
               facevarz(idvz,:,:,k_dest,idest) =                       & 
                                   recvz(idvz,:,:,k_source)

               Do j=nguard0+1,nguard0 + nyb -1,2
               Do i=nguard0+1,nguard0 + nxb -1,2

                If (curvilinear) Then
                area11= cell_area3(i+nd  ,j+nd  ,nzb+(1+nguard)*k3d)
                area21= cell_area3(i+nd  ,j+1+nd,nzb+(1+nguard)*k3d)
                area12= cell_area3(i+1+nd,j+nd  ,nzb+(1+nguard)*k3d)
                area22= cell_area3(i+1+nd,j+1+nd,nzb+(1+nguard)*k3d)
                Else
                area11 = dx*dy
                area21 = dx*dy
                area12 = dx*dy
                area22 = dx*dy
                End If

                b11 = (facevarz(idvz,i,j,nzb+(nguard0+1)*k3d,idest)    & 
                             - tempz(i,j) )*area11
                b21 = (facevarz(idvz,i,j+1,nzb+(nguard0+1)*k3d,idest)  & 
                             - tempz(i,j+1) )*area21
                b12 = (facevarz(idvz,i+1,j,nzb+(nguard0+1)*k3d,idest)  & 
                             - tempz(i+1,j) )*area12
                b22 = (facevarz(idvz,i+1,j+1,nzb+(nguard0+1)*k3d,idest)& 
                             - tempz(i+1,j+1) )*area22
                bsum = b11+b12+b21+b22
                Call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

                If (curvilinear) Then
                areax1= cell_area1(i+1+nd,j+nd  ,nzb+(1+nguard)*k3d)
                areax2= cell_area1(i+1+nd,j+1+nd,nzb+(1+nguard)*k3d)
                areay1= cell_area2(i+nd  ,j+1+nd,nzb+(1+nguard)*k3d)
                areay2= cell_area2(i+nd+1,j+1+nd,nzb+(1+nguard)*k3d)
                Else
                areax1= dy*dz
                areax2= dy*dz
                areay1= dx*dz
                areay2= dx*dz
                End If

                facevarx(idvx,i+1,j,nguard0*k3d+nzb,idest) =           & 
                   facevarx(idvx,i+1,j,nguard0*k3d+nzb,idest)          & 
                   + e1/areax1
                facevarx(idvx,i+1,j+1,nguard0*k3d+nzb,idest) =         & 
                   facevarx(idvx,i+1,j+1,nguard0*k3d+nzb,idest)        & 
                   + e2/areax2
                facevary(idvy,i,j+1,nguard0*k3d+nzb,idest) =           & 
                   facevary(idvy,i,j+1,nguard0*k3d+nzb,idest)          & 
                   - ea/areay1
                facevary(idvy,i+1,j+1,nguard0*k3d+nzb,idest) =         & 
                   facevary(idvy,i+1,j+1,nguard0*k3d+nzb,idest)        & 
                   - eb/areay2
               End Do  ! End Do i=nguard0+1,nguard0 + nxb -1,2
               End Do  ! End Do j=nguard0+1,nguard0 + nyb -1,2

            End If  ! End If (jf == 1)

          End If  ! End If (.Not.newchild)

       End Do  ! End Do jf = 1,nfaces

      End If  ! End If If (nodetype(isg) == 1.And.newchild(isg) ...
      End Do  ! End Do isg = 1, lnblocks
      End If  ! End If (lnblocks > 0)

!-----div B check
      If (lnblocks > 0) Then
      Do isg = 1,lnblocks

!-----Is this a newly created leaf block ?
      If (nodetype(isg) == 1) Then

               divbmax = 0.
               Do k=nguard0*k3d+1,nguard0*k3d + nzb
               Do j=nguard0+1,nguard0 + nyb
               Do i=nguard0+1,nguard0 + nxb
                 divb = (                                              & 
                  facevarx(1,i+1,j,k,isg) - facevarx(1,i,j,k,isg)      & 
                  + facevary(1,i,j+1,k,isg) - facevary(1,i,j,k,isg)    & 
                  + facevarz(1,i,j,k+k3d,isg) - facevarz(1,i,j,k,isg))
                divbmax = max(divbmax,abs(divb))
               End Do
               End Do
               End Do

      End If  ! End If (nodetype(isg) == 10
      End Do  ! End Do isg = 1, lnblocks
      End If  ! End If lnblock > 0)

      End If  ! End If (ndim > 1)

      Return

      Contains

      Subroutine compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

! Computes virtual electric field values required to produce
! changes in facevar's in a manner which will preserve div = 0
! constraints. The b11, b12 etc input arguments specify the
! required change in the vector component on the chosen face.
! The ea,eb,e1,e2 arguments return the virtual electric field
! values required to achieve this adjustment.
!
!
! The relationship between b's and e's is
!
!               -------------------------------------
!              |                  |                  |
!              |                  |                  |
!              |                  ^                  |
!              |       b21        e2     b22         |
!              |                  |                  |
!              |                  |                  |
!              |------ ea ->------|----- eb ->-------|
!              |                  |                  |
!              |                  ^                  |
!              |       b11        e1     b12         |
!              |                  |                  |
!              |                  |                  |
!              |                  |                  |
!               -------------------------------------
!
! These electric fields must be applied to adjust the appropriate
! components of the vector field in the planes perpendicular to
! the chosen face, immediately inside this face.


      Use paramesh_dimensions
      Use physicaldata

!-----Input/Output arguments.
      Integer :: isw

!-----Local variables and arrays.
      Real :: e1,e2,ea,eb
      Real :: b11,b12,b21,b22

#define RICK_DIVB
#ifdef ORIGINAL_DIVB
      ea = ( 2.*b11 +    b12 -    b21          )*.25
      eb = (    b11 + 2.*b12             - b22 )*.25
      e1 = (-2.*b11 +    b12 -    b21          )*.25
      e2 = (   -b11          - 2.*b21 + b22    )*.25
#endif /* ORIGINAL_DIVB */
#ifdef RICK_DIVB
      ea = (+ 5. * b11 + 3. * b12 - 1. * b21 + 1. * b22)*.125
      eb = (- 1. * b11 + 1. * b12 - 3. * b21 - 5. * b22)*.125
      e1 = (- 1. * b11 + 5. * b12 + 1. * b21 + 3. * b22)*.125
      e2 = (- 3. * b11 - 1. * b12 - 5. * b21 + 1. * b22)*.125
#endif /* RICK_DIVB */

      If (isw == 1) Then
       Write(*,*) 'compute : b ',b11,b12,b21,b22
       Write(*,*) 'compute : e ',ea,eb,e1,e2
       isw = 0
      End If

      Return
      End Subroutine compute_evalues

      End Subroutine amr_prolong_fc_divbconsist


