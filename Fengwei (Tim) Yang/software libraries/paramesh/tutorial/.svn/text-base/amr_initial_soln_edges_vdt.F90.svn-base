!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_initial_soln
!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr

      use physicaldata
      use tree

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "amr_shmem.fh"
        include         'shmem_reduce.fh'
!--------------------------------------------------------------

      integer nguard0
      parameter(nguard0 = nguard*npgs)

!--------------------------------------------------------------
      write(*,*) 'Enter nphase '
!      read(*,*) nphase 
      nphase = 2
      nvarp = nfacevar/(nphase+1)

! loop over leaf grid blocks
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif

      dx = size(1,lb)/real(nxb-gc_off_x)
      dy = size(2,lb)/real(nyb-gc_off_y*k2d)
      dz = size(3,lb)/real(nzb-gc_off_z*k3d)


! set values for unk
        do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd-nguard0

             do np = 1,nphase+1
              ivar1 = (np-1)*nvarp+1              
              ivar2 = (np-1)*nvarp+nvarp              
              do iv=ivar1,ivar2
               unk(iv,i,j,k,lb) = 1.0*real(iv-ivar1+1)
              enddo
              xi =  bnd_box(1,1,lb)  & 
     &           + dx*(real(i-nguard0)-.5*real(1+gc_off_x))
              yi =  bnd_box(1,2,lb)  & 
     &           + dy*(real(j-nguard0*k2d)-.5*real(1+gc_off_y*k2d))
              zi = 0.
              if(ndim.eq.3) zi =  bnd_box(1,3,lb)  & 
     &           + dz*(real(k-nguard0*k3d)-.5*real(1+gc_off_z*k3d))
                if( abs(xi).lt.1.01 .and. abs(yi).lt.1.01 .and. & 
     &            abs(zi).lt.1.01) then
                  do iv=ivar1,ivar2
                          unk(iv,i,j,k,lb) = 10.0*real(iv-ivar1+1)
                  enddo
               endif

             enddo

            enddo
          enddo
        enddo

! set values for facevar
        do np = 1,nphase+1
          ivar1 = (np-1)*nvarp+1              
          ivar2 = (np-1)*nvarp+nvarp              
          do iv=ivar1,ivar2
            facevarx(iv,:,:,:,lb) = .0*real(iv-ivar1+1)
            facevary(iv,:,:,:,lb) = .0*real(iv-ivar1+1)
            facevarz(iv,:,:,:,lb) = .0*real(iv-ivar1+1)
          enddo
        enddo

#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo                              ! end loop over grid blocks
      endif

      return
      end
