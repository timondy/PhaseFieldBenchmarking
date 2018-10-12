!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! This file has been modified for greater efficiency
!#include "paramesh_preprocessor.fh"


subroutine amr_test_refinement(pid,llrefine_min,llrefine_max, initialadapt)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  integer, intent(in) :: pid,llrefine_min,llrefine_max
  integer, optional :: initialadapt
  real :: dx,error,e1,e2,e3,error_max(1:3),error_mod(1:3),ctore,ctode, dist(8), dist2(4)
  integer :: ndel,nlayers,i,j,k,lb,iopt, inadapter
  logical :: allsolid
  ndel = (nguard_work - nguard)*npgs
  ! Re-initialize the refinement and derefinement flag arrays
  refine(:)   = .false.
  derefine(:) = .false.
  error_max(:) = 0.0
  ! Error limits which control the refinement and derefinement requests below.
!  ctore = .005!-1e30
!  ctode = .0025!-1e30
  ctore = .1!-1e30
  ctode = .01!-1e30
  error_mod(1)=1
  error_mod(2)=1
  error_mod(3)=1
  iopt=1
  nlayers=1
  inadapter=0
  if(present(initialadapt)) inadapter = initialadapt

! CEG: if doing uniform adaptation
!  return

!  print *, 'Now amr_guardcell'
!  call amr_guardcell(pid,iopt,nlayers)
  call pf_guardcell(pid, iopt, mg_max_level, 1, 1, 1)

!  call debug_block(pid)

  
!  print *,'Refinement loop'
  if(lnblocks.gt.0) then
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

           do current_var = 1,total_vars
              current_unk = 1 + (current_var-1)*nunkvbles
              error_max(current_var) = 0.0
!              dx = bsize(1,lb)/real(nxb)
              do k=kl_bnd+nguard*k3d, ku_bnd-nguard*k3d 
                 do j=jl_bnd+nguard*k2d, ju_bnd-nguard*k2d 
                    do i=il_bnd+nguard, iu_bnd-nguard
                       e1 = abs(unk(current_unk,i+1,j,k,lb)-unk(current_unk,i-1,j,k,lb))
                       e2 = abs(unk(current_unk,i,j+1,k,lb)-unk(current_unk,i,j-1,k,lb))
                       if (ndim.eq.3) then
                          e3 = abs(unk(current_unk,i,j,k+1,lb)-unk(current_unk,i,j,k-1,lb))
                          error = 10.0*(e1+e2+e3)
                       else
                          error = 10.0*(e1+e2)
                       endif
                       if (error.gt.error_max(current_var)) error_max(current_var)=error
                    enddo
                 enddo
              enddo
           end do
           allsolid = .true.
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                 do i=il_bnd+nguard,iu_bnd-nguard
                    if (unk(1,i,j,k,lb).lt.0.999) allsolid = .false.
                 enddo
              enddo
           enddo
           error=1.0/(error_mod(1)+error_mod(2)+error_mod(3))*&
                (error_mod(1)*error_max(1)+error_mod(2)*error_max(2)+error_mod(3)*error_max(3))
           if (allsolid) then
!              print *, 'Derefining block ', lb
              refine(lb) = .false.
              derefine(lb) = .true.
       	   else
              error=1.0/(error_mod(1)+error_mod(2)+error_mod(3))*&
               	(error_mod(1)*error_max(1)+error_mod(2)*error_max(2)+error_mod(3)*error_max(3))
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
  endif
end subroutine amr_test_refinement

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Growing the domain by testing to see if any top level blocks have far field boundaries, but also have children
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amr_grow_domain(pid, noprocs, step)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  use io

  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs, step
  integer :: b, bmax_p, cnt(2), Gcnt(2), lb, ierr, nghbdry, opp
  double precision :: ix, iy, iz
  double precision :: bkdist
  logical :: l_move_solution
  integer :: lnblocks_old, i, j, k, p, fc, bc, tmpI, sending_len
  integer, allocatable :: bklist(:)

  lnblocks_old = lnblocks

  cnt(:) = 0
!  cnt(2) = block_starts(2)-1

  do lb = 1, lnblocks
     if (lrefine(lb).eq.1) cnt(2) = cnt(2) + 1
  end do

  if (cnt(2).gt.0) then 
!print *,'allocate', cnt(2)
     allocate(bklist(cnt(2)))
     bklist(:) = 0
     cnt(2) = 0
     do lb = 1, lnblocks
        if (lrefine(lb).eq.1) then
           cnt(2) = cnt(2) + 1
!print *,'used[1]', cnt(2)
           bklist(cnt(2)) = lb
        endif
     end do

     do lb = 1, cnt(2)
!print *,'used[2]', lb
        if (child(1,1,bklist(lb)).gt.0) then
           do i = 1, 4+2*k3d
!           do i = 2, 4, 2
              if (neigh(1,i,bklist(lb)).eq.bc_cond_far) then
                 cnt(1) = cnt(1)+1
                 exit
              end if
           end do
        end if
     end do
  endif

  call MPI_Allreduce(cnt, Gcnt, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
!  print*,'grow_domain?', cnt, Gcnt 

! Gcnt(1) is the total number of top level blocks with children AND far field boundaries
! Gcnt(2) is the total number of top level blocks 
! bklist  is the list of top level blocks

!  if (newbkcnt.ne.lastbkcnt) then
!  if (Gcnt(1).eq.Gcnt(2)) then ! i.e. are there all blocks parents, hence necessitating domain growing 
  if (Gcnt(1).gt.0) then ! i.e. are there blocks parent with far-field boundaries, hence necessitating domain growing 
 
     if (pid.eq.0)  write (6,'(A,A)') 27, "[31;1;44m"
     if (pid.eq.0)  print *,'Coarse grid blocks grown - should grow domain'

     do p = 0, noprocs-1
 
!     if (cnt(1).eq.cnt(2).and.cnt(1).gt.0) then     ! i.e. this processor needs to grow its neighbours     

        if (pid.eq.p) then
!           do nghbdry = 2, 4, 2
           do nghbdry = 1, 4+2*k3d
!print*,'in ploop,',pid,p,nghbdry
              b = -1
              if (cnt(2).gt.0) then
                 do lb=1, cnt(2)
                    if (neigh(1,nghbdry,bklist(lb)).eq.bc_cond_far.and.child(1,1,bklist(lb)).gt.0) then
                       b = bklist(lb)

!                       print *, 'Process growing domain is', pid, 'direction', nghbdry
                       call get_neigh_dirs(nghbdry, ix, iy, iz, opp, sending_len)

                       ! Create new block
                       !print *, 'Creating block ', lnblocks+1
                       !print *, 'Next to ', b
                       lnblocks = lnblocks+1

                       parent(1,lnblocks) = -1
                       parent(2,lnblocks) = -1

                       bkdist = bnd_box(2,1,b)-bnd_box(1,1,b)         
                       !print *, 'bkdist = ',bkdist

                       ! Set centre of new block using direction information from ix/iy/iz variables
                       coord(1,lnblocks) = coord(1,b)+ix*bkdist    
                       coord(2,lnblocks) = coord(2,b)+iy*bkdist

                       ! Block dimensions
                       bnd_box(1,1,lnblocks) = bnd_box(1,1,b) + ix*bkdist
                       bnd_box(2,1,lnblocks) = bnd_box(2,1,b)	+ ix*bkdist

                       bnd_box(1,2,lnblocks) = bnd_box(1,2,b) + iy*bkdist
                       bnd_box(2,2,lnblocks) = bnd_box(2,2,b) + iy*bkdist

                       if (ndim.eq.3) then
                          coord(3,lnblocks) = coord(3,b)+iz*bkdist
                          bnd_box(1,3,lnblocks) = bnd_box(1,3,b) + iz*bkdist
                          bnd_box(2,3,lnblocks) = bnd_box(2,3,b) + iz*bkdist
                       endif

                       bsize(:,lnblocks) = bsize(:,b)
                       nodetype(lnblocks) = 1
                       lrefine(lnblocks) = 1
                       refine(lnblocks) = 0

!                       print*,bnd_box(:,:,lnblocks)
    
                       ! Now set initial conditions on this block 
                       call amr_initial_soln_blk(lnblocks)

                       ! Now do all the neighbours of the new block
                       ! Initially set all faces to be far boundary
                       do fc = 1, 4+k3d*2
                          bc = bc_cond_far
                          if (neigh(1,fc,b).eq.bc_cond_sym) bc = bc_cond_sym
                          neigh(1,fc,lnblocks) = bc
                          neigh(2,fc,lnblocks) = pid
                       end do
                       surr_blks(:,:,:,:,lnblocks) = -1

                       ! Now do face adjoining neighbouring block
                       neigh(1,opp,lnblocks) = b
                       neigh(2,opp,lnblocks) = pid

                       ! Do the neighbouring face of the old block
                       neigh(1,nghbdry,b) = lnblocks
                       neigh(2,nghbdry,b) = pid

                       ! Now check for other neighbouring blocks
                       call test_other_surr_faces(b, nghbdry, pid, noprocs, sending_len)

                       !print*,neigh(1,2,child(1,1,bmax)),neigh(1,2,child(1,2,bmax)),neigh(1,2,child(1,3,bmax)),neigh(1,2,child(1,4,bmax))

                       ! These bits need to done on all processors irregardless of where the extra blocks were added
                       ! See if anyone added any extra blocks in this direction
                       !call MPI_Bcast(b, 1, MPI_INTEGER, pid, MPI_COMM_WORLD, ierr)
                       !print *, pid, 'I know domain has grown'
                       call mpi_amr_global_domain_limits()
                       call set_domain_limits()
!     call debug_block(pid)
!                       print*,'pf_morton_process - new for ', pid
                       call pf_morton_process
!                       print*,'pf_morton_process - new done ', pid
!     call debug_block(pid)

                    endif !if (b) then

                 end do ! End of loop over blocks needing this direction adding
              else ! i.e no blocks to be grown on this processor
                 !call MPI_Bcast(b, 1, MPI_INTEGER, pid, MPI_COMM_WORLD, ierr) !Sending -1 to all the others
              endif ! end if (cnt(2).gt.0)  ! ie top level blocks

           end do ! End of nbdry=2->6
!print*,'end ploop,',pid,p

           call face_growing_ploop_end(pid, noprocs)
!print*,'face_growing_ploop_end done,',pid

        else ! i.e. pid != p (have we received blocks from another processor?)
!print*,'else ploop,',pid,p
           call get_neigh_dirs(1, ix, iy, iz, opp, sending_len)
           call process_other_surr_face(sending_len, p, pid)
!print*,'else ploop done,',pid,p
        endif 

     end do


!stop
!     if (pid.eq.0) print*,'At end of block adding.'
!     call debug_block(pid)
     refine(:)   = .false.
     derefine(:) = .false.
!     if (pid.eq.0) print*,'Now for rearrangement...'
     call amr_refine_derefine
!     call amr_multigrid_block_types(pid, noprocs)


!        call amr_mg_init() !NOTE FOR CHRIS, THIS IS NEW THING WHICH STOPPED SEG FAULT
!        nodetype_copy(1:lnblocks)=nodetype(1:lnblocks)
!        call amr_mg_get_max_refinement(pid, noprocs)
!        call amr_multigrid_child_set()
!        call reset_nodetype(pid, noprocs,mg_max_level)
        !call amr_mg_init() !NOTE FOR CHRIS, THIS IS OLD THING WHICH NEEDED TO BE MOVED UP
!        call amr_multigrid_block_types(pid, noprocs)


!     if (pid.eq.0) print*,'Order for this step is now...'
     call debug_block(pid)
     call graphviz_output(pid, noprocs, step)
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
     if (pid.eq.0) write (6,'(A,A)') 27, "[0m"

  else
     ! Not growing domain
     !call debug_block(pid)
     !print*,pid,'not grown'
  endif ! end of if for domain needing to grow 

  if (allocated(bklist)) deallocate(bklist)

!     call debug_block(pid)
end subroutine amr_grow_domain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_other_surr_faces(b, nghbdry, pid, noprocs, sending_len)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  use io
  implicit none
  include 'mpif.h'
  integer, intent(in) :: b, nghbdry, pid, noprocs, sending_len
  integer, allocatable :: bk(:), bkp(:), dir(:)
  integer ::  scat_cnt, ierr
  integer, allocatable :: scat_vec(:), sending_data(:)

  allocate(scat_vec(0:noprocs-1))
  scat_vec(:) = 0  
  allocate(sending_data(sending_len))  
  sending_data(1:sending_len) = pid
  if (ndim.eq.3) then 
     allocate(bk(4))
     allocate(bkp(4))
     allocate(dir(4))
  else
     allocate(bk(2))
     allocate(bkp(2))
     allocate(dir(2))
  endif

  sending_data(1) = lnblocks

  if (nghbdry.eq.1) then
     if (ndim.eq.3) then
!        print *, 'Missing 1-3d'
        bk(1)  = surr_blks(1,1,1,2,b)
        bkp(1) = surr_blks(2,1,1,2,b)
        bk(2)  = surr_blks(1,1,3,2,b)
        bkp(2) = surr_blks(2,1,3,2,b)
        dir(1) = 3
        dir(2) = 4
        bk(3)  = surr_blks(1,1,2,1,b)
        bkp(3) = surr_blks(2,1,2,1,b)
        bk(4)  = surr_blks(1,1,2,3,b)
        bkp(4) = surr_blks(2,1,2,3,b)
        dir(3) = 5
        dir(4) = 6

        call set_other_surr_faces_3d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     else
        bk(1)  = surr_blks(1,1,1,1,b)
        bkp(1) = surr_blks(2,1,1,1,b)
        bk(2)  = surr_blks(1,1,3,1,b)
        bkp(2) = surr_blks(2,1,3,1,b)
        dir(1) = 3
        dir(2) = 4

        call set_other_surr_faces_2d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     endif
  else if (nghbdry.eq.2) then         ! Moving right
     if (ndim.eq.3) then
        bk(1)  = surr_blks(1,3,1,2,b)
        bkp(1) = surr_blks(2,3,1,2,b)
        bk(2)  = surr_blks(1,3,3,2,b)
        bkp(2) = surr_blks(2,3,3,2,b)
        dir(1) = 3
        dir(2) = 4
        bk(3)  = surr_blks(1,3,2,1,b)
        bkp(3) = surr_blks(2,3,2,1,b)
        bk(4)  = surr_blks(1,3,2,3,b)
        bkp(4) = surr_blks(2,3,2,3,b)
        dir(3) = 5
        dir(4) = 6
!        print *, 'Missing 2-3d'

        call set_other_surr_faces_3d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     else
        bk(1)  = surr_blks(1,3,1,1,b)
        bkp(1) = surr_blks(2,3,1,1,b)
        bk(2)  = surr_blks(1,3,3,1,b)
        bkp(2) = surr_blks(2,3,3,1,b)
        dir(1) = 3
        dir(2) = 4

        call set_other_surr_faces_2d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     endif
  else if (nghbdry.eq.3) then  
     if (ndim.eq.3) then
        bk(1)  = surr_blks(1,1,1,2,b)
        bkp(1) = surr_blks(2,1,1,2,b)
        bk(2)  = surr_blks(1,3,1,2,b)
        bkp(2) = surr_blks(2,3,1,2,b)
        dir(1) = 1
        dir(2) = 2
        bk(3)  = surr_blks(1,2,1,1,b)
        bkp(3) = surr_blks(2,2,1,1,b)
        bk(4)  = surr_blks(1,2,1,3,b)
        bkp(4) = surr_blks(2,2,1,3,b)
        dir(3) = 5
        dir(4) = 6
!        print *, 'Missing 3-3d'

        call set_other_surr_faces_3d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     else
        bk(1)  = surr_blks(1,1,1,1,b)
        bkp(1) = surr_blks(2,1,1,1,b)
        bk(2)  = surr_blks(1,3,1,1,b)
        bkp(2) = surr_blks(2,3,1,1,b)
        dir(1) = 1
        dir(2) = 2

        call set_other_surr_faces_2d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     endif
  else if (nghbdry.eq.4) then  ! Moving up
     if (ndim.eq.3) then
        bk(1)  = surr_blks(1,1,3,2,b)
        bkp(1) = surr_blks(2,1,3,2,b)
        bk(2)  = surr_blks(1,3,3,2,b)
        bkp(2) = surr_blks(2,3,3,2,b)
        dir(1) = 1
        dir(2) = 2
        bk(3)  = surr_blks(1,2,3,1,b)
        bkp(3) = surr_blks(2,2,3,1,b)
        bk(4)  = surr_blks(1,2,3,3,b)
        bkp(4) = surr_blks(2,2,3,3,b)
        dir(3) = 5
        dir(4) = 6
!        print *, 'Missing 4-3d'

        call set_other_surr_faces_3d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     else
        bk(1)  = surr_blks(1,1,3,1,b)
        bkp(1) = surr_blks(2,1,3,1,b)
        bk(2)  = surr_blks(1,3,3,1,b)
        bkp(2) = surr_blks(2,3,3,1,b)
        dir(1) = 1
        dir(2) = 2

        call set_other_surr_faces_2d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
     endif
  else if (nghbdry.eq.5) then  
        bk(1)  = surr_blks(1,1,2,1,b)
        bkp(1) = surr_blks(2,1,2,1,b)
        bk(2)  = surr_blks(1,3,2,1,b)
        bkp(2) = surr_blks(2,3,2,1,b)
        dir(1) = 1
        dir(2) = 2
        bk(3)  = surr_blks(1,2,1,1,b)
        bkp(3) = surr_blks(2,2,1,1,b)
        bk(4)  = surr_blks(1,2,3,1,b)
        bkp(4) = surr_blks(2,2,3,1,b)
        dir(3) = 3
        dir(4) = 4
!    print *, 'Missing 5'

        call set_other_surr_faces_3d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
  else if (nghbdry.eq.6) then  
        bk(1)  = surr_blks(1,1,2,3,b)
        bkp(1) = surr_blks(2,1,2,3,b)
        bk(2)  = surr_blks(1,3,2,3,b)
        bkp(2) = surr_blks(2,3,2,3,b)
        dir(1) = 1
        dir(2) = 2
        bk(3)  = surr_blks(1,2,1,3,b)
        bkp(3) = surr_blks(2,2,1,3,b)
        bk(4)  = surr_blks(1,2,3,3,b)
        bkp(4) = surr_blks(2,2,3,3,b)
        dir(3) = 3
        dir(4) = 4
!    print *, 'Missing 6'

        call set_other_surr_faces_3d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
  endif

!  print*,'test_other_surr_faces sending', scat_vec(:), pid
  call MPI_Scatter(scat_vec(0), 1, MPI_INTEGER, scat_cnt, 1, MPI_INTEGER, pid, MPI_COMM_WORLD, ierr)

  if (ndim.eq.3) then
     if (sending_data(4).ne.pid)  then
!        print *,pid,'sending to ',sending_data(4)
        call MPI_Send(sending_data, sending_len, MPI_INTEGER, sending_data(4), 136, MPI_COMM_WORLD, ierr)
!        print *,pid,' sent[4] ', sending_data
     endif
     if (sending_data(7).ne.pid.and.sending_data(7).ne.sending_data(4)) then
!        print *,pid,'sending to ',sending_data(7)
        call MPI_Send(sending_data, sending_len, MPI_INTEGER, sending_data(7), 136, MPI_COMM_WORLD, ierr)
!        print *,pid,' sent[7] ', sending_data
     endif
     if (sending_data(10).ne.pid.and. sending_data(10).ne.sending_data(7).and. sending_data(10).ne.sending_data(4)) then
!        print *,pid,'sending to ',sending_data(7)
        call MPI_Send(sending_data, sending_len, MPI_INTEGER, sending_data(10), 136, MPI_COMM_WORLD, ierr)
!        print *,pid,' sent[7] ', sending_data
     endif
     if (sending_data(13).ne.pid.and. sending_data(13).ne.sending_data(10).and. sending_data(13).ne.sending_data(7).and. &
         sending_data(13).ne.sending_data(4)) then
!        print *,pid,'sending to ',sending_data(7)
        call MPI_Send(sending_data, sending_len, MPI_INTEGER, sending_data(13), 136, MPI_COMM_WORLD, ierr)
!        print *,pid,' sent[7] ', sending_data
     endif
!     print *, '3d sender not done'
  else
     if (sending_data(4).ne.pid)  then
!        print *,pid,'sending to ',sending_data(4)
        call MPI_Send(sending_data, sending_len, MPI_INTEGER, sending_data(4), 136, MPI_COMM_WORLD, ierr)
!        print *,pid,' sent[4] ', sending_data
     endif
     if (sending_data(7).ne.pid.and.sending_data(7).ne.sending_data(4)) then
!        print *,pid,'sending to ',sending_data(7)
        call MPI_Send(sending_data, sending_len, MPI_INTEGER, sending_data(7), 136, MPI_COMM_WORLD, ierr)
!        print *,pid,' sent[7] ', sending_data
     endif
  endif

  deallocate(scat_vec)
  deallocate(sending_data)
  deallocate(bk)
  deallocate(bkp)
  deallocate(dir)

end subroutine test_other_surr_faces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Generic routine to process the non-growth-direction adjoining blocks
!!!!  to have correct info on the new one 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_other_surr_faces_3d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  use io
  implicit none
  integer bk(*), bkp(*), dir(*), pid, noprocs, sending_len
  integer :: scat_vec(0:noprocs-1), sending_data(1:sending_len)

  if (bk(1).gt.0) then  
     ! My new block
     neigh(1,dir(1),lnblocks) = bk(1)              
     neigh(2,dir(1),lnblocks) = bkp(1)
     ! This other block
     if (bkp(1).eq.pid) then
        neigh(1,dir(2),bk(1)) = lnblocks
        neigh(2,dir(2),bk(1)) = pid
     else
        scat_vec(bkp(1)) = scat_vec(bkp(1))+1
        sending_data(2) = dir(2)
        sending_data(3) = bk(1)
        sending_data(4) = bkp(1)
     endif
  endif

  if (bk(2).gt.0) then  
     ! My new block
     neigh(1,dir(2),lnblocks) = bk(2)              
     neigh(2,dir(2),lnblocks) = bkp(2)
     ! This other block
     if (bkp(2).eq.pid) then
        neigh(1,dir(1),bk(2)) = lnblocks
        neigh(2,dir(1),bk(2)) = pid
     else
        scat_vec(bkp(2)) = scat_vec(bkp(2))+1
        sending_data(5) = dir(1)
        sending_data(6) = bk(2)
        sending_data(7) = bkp(2)
     endif
  endif

  if (bk(3).gt.0) then  
     ! My new block
     neigh(1,dir(3),lnblocks) = bk(3)              
     neigh(2,dir(3),lnblocks) = bkp(3)
     ! This other block
     if (bkp(3).eq.pid) then
        neigh(1,dir(4),bk(3)) = lnblocks
        neigh(2,dir(4),bk(3)) = pid
     else
        scat_vec(bkp(3)) = scat_vec(bkp(3))+1
        sending_data(8) = dir(4)
        sending_data(9) = bk(3)
        sending_data(10) = bkp(3)
     endif
  endif

  if (bk(4).gt.0) then  
     ! My new block
     neigh(1,dir(4),lnblocks) = bk(4)              
     neigh(2,dir(4),lnblocks) = bkp(4)
     ! This other block
     if (bkp(4).eq.pid) then
        neigh(1,dir(3),bk(4)) = lnblocks
        neigh(2,dir(3),bk(4)) = pid
     else
        scat_vec(bkp(4)) = scat_vec(bkp(4))+1
        sending_data(11) = dir(3)
        sending_data(12) = bk(4)
        sending_data(13) = bkp(4)
     endif
  endif

  return

end subroutine set_other_surr_faces_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Generic routine to process the non-growth-direction adjoining blocks
!!!!  to have correct info on the new one 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_other_surr_faces_2d(bk, bkp, dir, scat_vec, sending_data, pid, noprocs, sending_len)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  use io
  implicit none
  integer bk(*), bkp(*), dir(*), pid, noprocs, sending_len
  integer :: scat_vec(0:noprocs-1), sending_data(1:sending_len)

  if (bk(1).gt.0) then  
     ! My new block
     neigh(1,dir(1),lnblocks) = bk(1)              
     neigh(2,dir(1),lnblocks) = bkp(1)
     ! This other block
     if (bkp(1).eq.pid) then
        neigh(1,dir(2),bk(1)) = lnblocks
        neigh(2,dir(2),bk(1)) = pid
     else
        scat_vec(bkp(1)) = scat_vec(bkp(1))+1
        sending_data(2) = dir(2)
        sending_data(3) = bk(1)
        sending_data(4) = bkp(1)
     endif
  endif

  if (bk(2).gt.0) then  
     ! My new block
     neigh(1,dir(2),lnblocks) = bk(2)              
     neigh(2,dir(2),lnblocks) = bkp(2)
     ! This other block
     if (bkp(2).eq.pid) then
        neigh(1,dir(1),bk(2)) = lnblocks
        neigh(2,dir(1),bk(2)) = pid
     else
        scat_vec(bkp(2)) = scat_vec(bkp(2))+1
        sending_data(5) = dir(1)
        sending_data(6) = bk(2)
        sending_data(7) = bkp(2)
     endif
  endif

  return

end subroutine set_other_surr_faces_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Process the data from every added block.
!!!! Recursive since we don't know how many will be coming from each processor.
!!!! The MPI_Scatter gets how many blocks are needing work on this processor, stores in scat_cnt.
!!!! If scat_cnt=0 then there are new blocks elsewhere, but surr_blks may still need updating.
!!!! If scat_cnt=-1 then we can quit.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine process_other_surr_face(sending_len, p, pid)
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  use io
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, p, sending_len
  integer :: scat_cnt, scat_vec, ierr, tag, stat(MPI_STATUS_SIZE)
  integer, allocatable :: sending_data(:)

!  print *, 'process_other_surr_face pre', pid, p
  call MPI_Scatter(scat_vec, -1, MPI_INTEGER, scat_cnt, 1, MPI_INTEGER, p, MPI_COMM_WORLD, ierr)  !First three arguments dummy

!  print *, 'process_other_surr_face', pid, p, scat_cnt

  if (scat_cnt.gt.0) then
     
     allocate(sending_data(sending_len))  
     tag = 136
!     print *,pid,' receiving from ', p, sending_len
     call MPI_Recv(sending_data, sending_len, MPI_INTEGER, p, tag, MPI_COMM_WORLD, stat, ierr)
!     print *,pid,' received ', sending_data

        if (sending_data(4).eq.pid) then
           neigh(1,sending_data(2),sending_data(3)) = sending_data(1)
           neigh(2,sending_data(2),sending_data(3)) = p
        endif
        if (sending_data(7).eq.pid) then
           neigh(1,sending_data(5),sending_data(6)) = sending_data(1)
           neigh(2,sending_data(5),sending_data(6)) = p
        endif
        if (ndim.eq.3) then
           if (sending_data(10).eq.pid) then
              neigh(1,sending_data(8),sending_data(9)) = sending_data(1)
              neigh(2,sending_data(8),sending_data(9)) = p
           endif
           if (sending_data(13).eq.pid) then
              neigh(1,sending_data(11),sending_data(12)) = sending_data(1)
              neigh(2,sending_data(11),sending_data(12)) = p
           endif
        endif
     deallocate(sending_data)

  else if (scat_cnt.eq.0) then
!     print *, 'New block has no connections to processor ', pid
  endif

  if (scat_cnt.ge.0) then ! New block somewhere, so sort the surr_blks array plus domain boundaries and then recurse

     call mpi_amr_global_domain_limits()
     call set_domain_limits()
!     print *, 'pf_morton_process - from other proc ', pid 
     call pf_morton_process
!     print *, 'pf_morton_process - other proc done', pid 

     call process_other_surr_face(sending_len, p, pid)

  endif

end subroutine process_other_surr_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Sends out a "-1" in the scatter to each other processor to say move on
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine face_growing_ploop_end(pid, noprocs)
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer, allocatable :: scat_vec(:)
  integer :: scat_cnt, ierr

  allocate(scat_vec(0:noprocs-1))
  scat_vec(:) = -1

!  print *, 'face_growing_ploop_end ', pid, noprocs, scat_vec
  call MPI_Scatter(scat_vec(0), 1, MPI_INTEGER, scat_cnt, 1, MPI_INTEGER, pid, MPI_COMM_WORLD, ierr)

end subroutine face_growing_ploop_end  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Does parameter setting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_neigh_dirs(nghbdry, ix, iy, iz, opp, sending_len)
  use paramesh_dimensions
  implicit none
  integer, intent(inout) :: nghbdry, opp, sending_len
  double precision, intent(inout) :: ix, iy, iz

  ix=0.0
  iy=0.0
  iz=0.0
  select case (nghbdry)
     case(1)
        ix = -1.0
        opp = 2
     case(2)
        ix =  1.0
        opp = 1
     case(3)
        iy = -1.0
        opp = 4
     case(4)
        iy =  1.0
        opp = 3
     case(5)
        iz = -1.0
        opp = 6
     case(6)
        iz =  1.0
        opp = 5
  end select

  if (ndim.eq.3) then
     sending_len = 13         ! 4 points to send in 3-d, thus needing to know new block number, appropriate bdry, bk? values, and the remote pid 
  else !2d
     sending_len = 7          ! 2 points to send in 2-d, thus needing to know new block number, appropriate bdry, bk? values, and the remote pid 
  endif

  return

end subroutine get_neigh_dirs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Prints out neigh() and surr_blks() arrays for all top level blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine debug_block(pid)
  use tree
  implicit none
  integer :: lb, i, pid

  do lb=1, lnblocks
     if (lrefine(lb).eq.1) then
!        print *, pid, ' block ', lb, parent(1,lb),parent(2,lb)

        if (1.lt.2) then
           write(6,'(I3,X,I3,X,I3,X,I3,X,I3,X,I3,X,I3,X,I3)')pid, lb,neigh(1,1,lb),neigh(1,2,lb),neigh(1,3,lb),neigh(1,4,lb),neigh(1,5,lb),neigh(1,6,lb)
           write(6,'(I3,X,I3,X,I3,X,I3,X,I3,X,I3,X,I3,X,I3)')pid, lb,neigh(2,1,lb),neigh(2,2,lb),neigh(2,3,lb),neigh(2,4,lb),neigh(2,5,lb),neigh(2,6,lb)

           if (3.lt.2) then
              print *,pid, surr_blks(1,1:3,3,1,lb)
              print *,pid, surr_blks(2,1:3,3,1,lb)

              print *,pid, surr_blks(1,1:3,2,1,lb)
              print *,pid, surr_blks(2,1:3,2,1,lb)

              print *,pid, surr_blks(1,1:3,1,1,lb)
              print *,pid, surr_blks(2,1:3,1,1,lb)
           endif
        endif

        if (3.lt.2) then
           do i=1,4
!           do i=1,2
              if (child(2,i,lb).eq.pid) then
                 print*,pid, child(1,i,lb),neigh(1,1,child(1,i,lb)),neigh(1,2,child(1,i,lb)),neigh(1,3,child(1,i,lb)),neigh(1,4,child(1,i,lb))
                 print*,pid, child(2,i,lb),neigh(2,1,child(1,i,lb)),neigh(2,2,child(1,i,lb)),neigh(2,3,child(1,i,lb)),neigh(2,4,child(1,i,lb))

                 print *,pid, surr_blks(1,1:3,3,1,child(1,i,lb))
                 print *,pid, surr_blks(2,1:3,3,1,child(1,i,lb))

                 print *,pid, surr_blks(1,1:3,2,1,child(1,i,lb))
                 print *,pid, surr_blks(2,1:3,2,1,child(1,i,lb))

                 print *,pid, surr_blks(1,1:3,1,1,child(1,i,lb))
                 print *,pid, surr_blks(2,1:3,1,1,child(1,i,lb))
              else
                 print*,pid, child(1,i,lb), child(2,i,lb),'off proc'
              endif
           end do
        endif
     endif
  end do

end subroutine debug_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Output neigh() information for generation of GraphViz (using neato) pictures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine graphviz_output(pid, noprocs, step)
  use paramesh_dimensions ! for ndim
  use tree                ! for neigh()
  use io                  ! for output_dir
  implicit none
  integer :: lb, i, pid, noprocs, step
  character*80 :: fname
  character (len=2)  :: proc_string
  character (len=2)  :: np_string
  character (len=5)  :: step_string

   Write (proc_string, '(i2.2)') pid
   Write (np_string, '(i2.2)') noprocs
   Write (step_string, '(i5.5)') step

  fname = trim(output_dir)// '/cgblocks_' // np_string // "_" // step_string // "." // proc_string
  open(unit=62,file=fname)

  i=0

  do lb=1, lnblocks
     if (lrefine(lb).eq.1) i = i+1
  end do

  write(62,*) i

  if (ndim.eq.3) then
     do lb=1, lnblocks
        if (lrefine(lb).eq.1) then
           write(62,'(I3,X,I3,X,I3,X,I3,X,I3,X,I3)')neigh(1,1,lb),neigh(1,2,lb),neigh(1,3,lb),neigh(1,4,lb),neigh(1,5,lb),neigh(1,6,lb)
           write(62,'(I3,X,I3,X,I3,X,I3,X,I3,X,I3)')neigh(2,1,lb),neigh(2,2,lb),neigh(2,3,lb),neigh(2,4,lb),neigh(2,5,lb),neigh(2,6,lb)
        endif
     end do
  else
     do lb=1, lnblocks
        if (lrefine(lb).eq.1) then
           write(62,'(I3,X,I3,X,I3,X,I3,X,I3,X,I3)')neigh(1,1,lb),neigh(1,2,lb),neigh(1,3,lb),neigh(1,4,lb),-23,-23
           write(62,'(I3,X,I3,X,I3,X,I3,X,I3,X,I3)')neigh(2,1,lb),neigh(2,2,lb),neigh(2,3,lb),neigh(2,4,lb),-23,-23
        endif
     end do
  endif
  close(62)

end subroutine graphviz_output

