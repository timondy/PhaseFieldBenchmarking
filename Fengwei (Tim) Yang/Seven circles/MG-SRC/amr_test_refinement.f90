!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  amr_test_refinement
!!!!  * Default mesh refinement testing routine
!!!!  * Looks at gradients in cell to decide if needs to be refined
!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Chris Goodyer, May 2012                                              !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! This file has been modified for greater efficiency


subroutine amr_test_refinement(pid, llrefine_min, llrefine_max, initialadapt)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  use refinement_parameters
  implicit none
  integer, intent(in) :: pid, llrefine_min, llrefine_max
  integer, optional :: initialadapt
  double precision :: dx,error,e1,e2,e3,error_max(total_vars),error_mod(total_vars), dist(8), dist2(4), errnum, errdenom
  integer :: ndel,nlayers,i,j,k,lb,iopt, inadapter, u, v
  logical :: allsolid
  ndel = (nguard_work - nguard)*npgs
  ! Re-initialize the refinement and derefinement flag arrays
  refine(:)   = .false.
  derefine(:) = .false.
  error_max(:) = 0.0
  error_mod(:)=1.0
  iopt=1
  nlayers=1
  inadapter=0
  if(present(initialadapt)) inadapter = initialadapt

! CEG: if doing uniform adaptation
!  return

!  print *, 'Now amr_guardcell'
!  call amr_guardcell(pid,iopt,nlayers)
!  call pf_guardcell(pid, iopt, mg_max_level, 1, 1, 1)

!  call debug_block(pid)

  if (allsolid_test) then
     allsolid = .true.
  else
     allsolid = .false.
  endif

!  print *,'Refinement loop'
  if(lnblocks.gt.0) then
!!!$OMP PARALLEL SHARED(lnblocks) PRIVATE(lb,i,j,k,u,v,e1,e2,e3,error,errdenom,errnum,error_max,dist,dist2)
!!!$OMP DO SCHEDULE(DYNAMIC,1)
     do lb=1,lnblocks

        if (inadapter.eq.1) then
           if (ndim.eq.2) then
              dist2(1) = abs(bnd_box(1,1,lb))+abs(bnd_box(1,2,lb))
              dist2(2) = abs(bnd_box(1,1,lb))+abs(bnd_box(2,2,lb))
              dist2(3) = abs(bnd_box(2,1,lb))+abs(bnd_box(1,2,lb))
              dist2(4) = abs(bnd_box(2,1,lb))+abs(bnd_box(2,2,lb))

              if (MINVAL(dist2).lt.1.0) refine(lb) = .true.

           else           
              dist(1) = abs(bnd_box(1,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(1,3,lb))
              dist(2) = abs(bnd_box(1,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(2,3,lb))
              dist(3) = abs(bnd_box(1,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(1,3,lb))
              dist(4) = abs(bnd_box(1,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(2,3,lb))
              dist(5) = abs(bnd_box(2,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(1,3,lb))
              dist(6) = abs(bnd_box(2,1,lb))+abs(bnd_box(1,2,lb))+abs(bnd_box(2,3,lb))
              dist(7) = abs(bnd_box(2,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(1,3,lb))
              dist(8) = abs(bnd_box(2,1,lb))+abs(bnd_box(2,2,lb))+abs(bnd_box(2,3,lb))

              if (MINVAL(dist).lt.1.0) refine(lb) = .true.

           endif

!        else if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
        else if(nodetype(lb).eq.1) then

           do v = 1,total_vars
              u = 1 + (v-1)*nunkvbles
              error_max(v) = 0.0
!              dx = bsize(1,lb)/real(nxb)
              do k=kl_bnd+nguard*k3d, ku_bnd-nguard*k3d 
                 do j=jl_bnd+nguard*k2d, ju_bnd-nguard*k2d 
                    do i=il_bnd+nguard, iu_bnd-nguard
                       e1 = abs(unk(u,i+1,j,k,lb)-unk(u,i-1,j,k,lb))
                       e2 = abs(unk(u,i,j+1,k,lb)-unk(u,i,j-1,k,lb))
                       if (ndim.eq.3) then
                          e3 = abs(unk(u,i,j,k+1,lb)-unk(u,i,j,k-1,lb))
                          error = 10.0*(e1+e2+e3)
                       else
                          error = 10.0*(e1+e2)
                       endif
                       if (error.gt.error_max(v)) error_max(v)=error
                       if (unk(1,i,j,k,lb).lt.0.999) allsolid = .false.
                    enddo
                 enddo
              enddo
           end do

           if (allsolid) then
!              print *, 'Derefining block ', lb
              refine(lb) = .false.
              derefine(lb) = .true.
       	   else
              errdenom = 0.0
              errnum = 0.0
              do v = 1, total_vars
                 errdenom = errdenom + error_mod(v)
                 errnum = errnum + error_mod(v)*error_max(v)
              end do
!              error=1.0/(error_mod(1)+error_mod(2)+error_mod(3))*&
!                (error_mod(1)*error_max(1)+error_mod(2)*error_max(2)+error_mod(3)*error_max(3))
              error=errnum/errdenom
              if( lrefine(lb).lt.llrefine_max .and.nodetype(lb).eq.1) then
                 if ( error .ge. 0.1*ctore)then
!                    print *, 'Refining block ', lb, pid, parent(:,lb), error, 0.1*ctore
                    refine(lb) = .true.
                    if (derefine(lb))derefine(lb)=.false.
                 else
!                    print *, 'Not refining   ', lb, pid, parent(:,lb), error, 0.1*ctore
                 end if
              endif
              if( lrefine(lb).gt.llrefine_min .and. (.not.refine(lb)) .and. nodetype(lb).eq.1) then
                 if ( error .lt. ctode) then
                    if (lrefine(lb).gt.mg_min_lvl) derefine(lb) = .true.
                 endif
              endif
           endif
        endif
     end do
!!!$OMP END DO
! NOWAIT
!!!$OMP END PARALLEL 
  endif
end subroutine amr_test_refinement

