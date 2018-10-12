! This file contains routines for balancing the distribution of block for MLAT versions of PARAMESH

!subroutine ceg_block_balancing(pid, noprocs)
subroutine ceg_block_balancing_link(pid, noprocs)

  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces
  integer, intent(in) :: pid, noprocs
  integer, dimension(2, maxblocks_tr) :: new_loc
!  integer, dimension(2, maxblocks_tr), intent(inout) :: new_loc

interface
  subroutine ceg_block_balancing(pid, noprocs, new_loc)
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces
    integer, intent(in) :: pid, noprocs
    integer, intent(inout) :: new_loc(:,:)
    end subroutine ceg_block_balancing
end interface

  call ceg_block_balancing(pid, noprocs, new_loc)

!  print*,'out of ceg_block_balancing_link'

end subroutine ceg_block_balancing_link


subroutine ceg_block_balancing(pid, noprocs, new_loc)

  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer :: locblocks(lrefine_max), blocks(lrefine_max*noprocs), origblocks(lrefine_max*noprocs)
  integer :: g, b, p, pp, proc, procp, T, ierr, mvd, wanted, total
  integer :: bksend(lrefine_max*noprocs, noprocs), bkrecv(lrefine_max*noprocs, noprocs), tgt(0:noprocs-1)
!  integer, dimension(2, maxblocks_tr) :: new_loc
  integer, dimension(2, maxblocks_tr), intent(inout) :: new_loc
!  integer, intent(inout) :: new_loc(:,:)
  integer :: i, j

interface
    subroutine ceg_setup_newloc(pid, noprocs, bksend, bkrecv, origblocks, new_loc)
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces
    integer, intent(in) :: pid, noprocs
    integer, intent(in) :: origblocks(:)
    integer, intent(in) :: bksend(lrefine_max*noprocs, 0:noprocs-1), bkrecv(lrefine_max*noprocs, 0:noprocs-1)
!    integer, intent(in) :: bksend(:,:), bkrecv(:,:)
    integer, intent(inout) :: new_loc(:,:)
    end subroutine ceg_setup_newloc
end interface


  new_loc(:,:)=-1 

!  print *,pid,' In ceg_mesh_partitioning : maxblocks_tr=',maxblocks_tr
  locblocks(:) = 0

  do b=1,lnblocks
     ! Put all data into blocks on each proc, then gather it on proc 1 into the final data structs
     locblocks(lrefine(b)) = locblocks(lrefine(b))+1
  end do

  call MPI_AllGather(locblocks, lrefine_max, MPI_INTEGER, &
       blocks, lrefine_max, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  do b=1, lrefine_max*noprocs
     origblocks(b) = blocks(b)
  end do

  !print *,pid, ' Gather done'
!  do proc=0, noprocs-1
!     if (pid.eq.proc) then
        !print *, '***** p=', pid
!        call coutput2(lrefine_max, noprocs, blocks)
        !print *, '*****'
!     end if
!     call MPI_Barrier(MPI_COMM_WORLD, ierr)
!  end do

! Let's look at the parents
!  do p = 0, noprocs-1
!    if (pid.eq.p) then
!!     do g = lrefine_max, 2, -1
!        g = lrefine_max
!        do b = 1, lnblocks
!
!           if (lrefine(b).eq.g.and.parent(2,b).ne.pid) then
!              write(6,*)pid,b,g,parent(1,b),parent(2,b), which_child(b)
!              do j = 1, 6
!                 write(6,*) j, neigh(1,j,b), neigh(2,j,b)
!              end do
!           endif
!
!        end do
!!     end do
!    endif
!    call MPI_Barrier(MPI_COMM_WORLD, ierr)
!  end do

! BALANCING

  bksend(:,:) = 0
  bkrecv(:,:) = 0

  do g = 1, lrefine_max

!     if (pid.eq.0) print *,pid, ' On grid ',g
     total = 0
! need to stick right array elements in throughout
     do p = 0, noprocs-1
        total = total + blocks(g+p*lrefine_max)
     end do

! Loop over the processors to balance the blocks
     do proc=0, noprocs-1
        tgt(proc) = total / (noprocs-proc)
        total = total - tgt(proc)
     end do

     do proc=0, noprocs-1
        T = tgt(proc)
!        if (pid.eq.0) then 
!           print *,pid, 'T=', T, proc
!           do p=0, noprocs-1
!              print *,p, 'tgt, act', tgt(p), blocks(g + lrefine_max*p)
!           end do
!	endif
!        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        p = g + lrefine_max*proc
        if (blocks(p).eq.T) then

! do nothing
           wanted = 0

        else if (blocks(p).gt.T) then
           wanted = blocks(p)-T

           procp = proc+1
           pp = p + lrefine_max
           do while (wanted.gt.0)
              do while (blocks(pp).ge.tgt(procp))
                 procp = procp+1
                 pp = pp + lrefine_max
              end do

              if (tgt(procp)-blocks(pp) .lt. wanted) then
                 mvd = tgt(procp)-blocks(pp)
                 bksend(g+procp*lrefine_max, proc+1) = bksend(g+procp*lrefine_max, proc+1) + mvd
                 bkrecv(g+proc*lrefine_max, procp+1) = bkrecv(g+proc*lrefine_max, procp+1) + mvd
                 blocks(pp) = blocks(pp) + mvd
                 blocks(p) = blocks(p) - mvd
              else
                 mvd = blocks(p)-tgt(proc)
                 bksend(g+procp*lrefine_max, proc+1) = bksend(g+procp*lrefine_max, proc+1) + mvd
                 bkrecv(g+proc*lrefine_max, procp+1) = bkrecv(g+proc*lrefine_max, procp+1) + mvd
                 blocks(pp) = blocks(pp) + mvd
                 blocks(p) = blocks(p) - mvd
              endif
              wanted = wanted - mvd

!              if (pid.eq.0) print*,pid,'Moving right ',procp, mvd,blocks(p),blocks(pp)
!              if (procp.ge.noprocs) print *,'Too many procs right'
!              if (procp.ge.noprocs) stop

           end do


        else

           wanted = T - blocks(p)
           !print*,'wanted ',wanted,T,blocks(p)
           procp = proc+1
           pp = p + lrefine_max
           !if (procp.ge.noprocs) print *,'Too many procs left'
           !if (procp.ge.noprocs) stop
           do while (wanted.gt.0)
              if (blocks(pp) .ge. wanted) then
                 mvd = wanted
                 !print*,pid,'Moving left ',mvd,blocks(p),blocks(pp)
                 bksend(g+proc*lrefine_max, procp+1) = bksend(g+proc*lrefine_max, procp+1) + wanted
                 bkrecv(g+procp*lrefine_max, proc+1) = bkrecv(g+procp*lrefine_max, proc+1) + wanted
                 blocks(p) = T
                 blocks(pp) = blocks(pp) - wanted
                 wanted = 0
              else
                 mvd = blocks(pp)
                 !print*,pid,'Shifting left ',mvd,blocks(p),blocks(pp)
                 bksend(g+proc*lrefine_max, procp+1) = bksend(g+proc*lrefine_max, procp+1) + mvd
                 bkrecv(g+procp*lrefine_max, proc+1) = bkrecv(g+procp*lrefine_max, proc+1) + mvd
                 blocks(p) = blocks(p) + mvd
                 blocks(pp) = 0
                 wanted = wanted - mvd
              end if

              procp = procp+1
              pp = pp + lrefine_max
           end do

        end if

     end do

! output loop
!     do proc=0, noprocs-1           ! proc only used here for sequencing the output
!        if (pid.eq.proc) then
!           print *, '***** p,g=', pid, g
!           do p=0, noprocs-1
!              if (pid.ne.p) print *, p, bksend(g+p*lrefine_max, pid+1), bkrecv(g+p*lrefine_max, pid+1)
!           end do        
!           print *, '*****'
!        end if
!        call MPI_Barrier(MPI_COMM_WORLD, ierr)
!     end do
     ! Now move blocks 

  end do

  if (pid.eq.0) then
     print *, '*****'
        call coutput2(lrefine_max, noprocs, blocks)
     print *, '*****'
  end if
 ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

!  print *,'Before new_loc(1)=',new_loc(1,1),new_loc(2,1)
!  do b=1, lrefine_max*noprocs
!     print *, 'orig', origblocks(b), blocks(b)
!  end do

  call ceg_setup_newloc(pid, noprocs, bksend, bkrecv, origblocks, new_loc)

!  do proc=0, noprocs-1           ! proc only used here for sequencing the output
!     if (pid.eq.proc) then
!        print *,'Block distribution on processor ', pid
!!        do b = 1, 10
!        do b = 1, lnblocks
!           print *, b, new_loc(1,b), new_loc(2,b), lrefine(b)
!        end do
!      endif
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!  end do

! These stolen from mpi_amr_refine_derefine.F90 and this should now be called from there so no need
!  call amr_migrate_tree_data (new_loc, noprocs, pid)
!
!  call amr_redist_blk (new_loc, noprocs, pid, lnblocks)
!
!  lnblocks = new_lnblocks

!  call amr_morton_process()

!-----set grid modification flag
!  grid_changed = 1
!  grid_analysed_mpi = 1

!  call mpi_amr_boundary_block_info(pid,noprocs)



  return
end subroutine ceg_block_balancing


subroutine ceg_setup_newloc(pid, noprocs, bksend, bkrecv, origblocks, new_loc)

  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid, noprocs
  integer, intent(in) :: origblocks(:)
!  integer, intent(in) :: bksend(:,:), bkrecv(:,:)
  integer, intent(inout) :: new_loc(:,:)
!  integer, intent(inout) :: origblocks(lrefine_max*noprocs)
  integer, intent(in) :: bksend(lrefine_max*noprocs, 0:noprocs-1), bkrecv(lrefine_max*noprocs, 0:noprocs-1)
!  integer, intent(inout) :: new_loc(2,maxblocks_tr)
  integer :: g, b, c, p, proc, sentdown, sentup, cnt
  integer :: new_cnt(0:noprocs-1)
  integer :: gstart(lrefine_max), gfound, reordered(lnblocks)
  logical :: escape

!  print *,'in ceg_setup_newloc'
!  print *,'pid, noprocs', pid, noprocs

!  do p=0, noprocs-1
!     do b=1, lrefine_max*noprocs
!        print *, 'bksend', b, p, bksend(b, p)
!     end do
!  end do
!  do b=1, lrefine_max*noprocs
!     print *, 'orig', origblocks(b)
!  end do

!  print *,'new_loc(1)=',new_loc(1,1),new_loc(2,1)

! The list of blocks isn't in grid order at the moment, so let's remedy that as we make the list
! Need to find the offsets through the datalist for the next block at that level
  gstart(:) = -1
  gfound = 0
  b=1
  ! escape present to account for cases of zero blocks at that level on this processor
  escape = .false.
  do while (.not.escape.and.gfound.lt.lrefine_max)
     if (gstart(lrefine(b)).eq.-1) then
        gstart(lrefine(b)) = b
        gfound = gfound + 1
     endif
     b = b + 1
     if (b.gt.lnblocks) escape = .true.
  end do

!  do g = 1, lrefine_max
!     print *, 'g, gstart', g, gstart(g)
!  end do 

  ! form a reordered block list in grid order
  cnt = 1
  do g = 1, lrefine_max
     do b = 1, origblocks(g+pid*lrefine_max)
        reordered(cnt) = gstart(g)
        cnt = cnt + 1
        if (b.lt.origblocks(g+pid*lrefine_max)) then
           escape = .false.
           c = gstart(g) + 1
           do while (.not.escape)
              if (lrefine(c).eq.g) then
                 gstart(g) = c
!                 print *, 'g, gstart', g, gstart(g)
                 escape = .true.
              endif
              c = c + 1
           end do
        end if
     end do
  end do

!  do b=1, lnblocks
!     print *,b, reordered(b), lrefine(reordered(b))
!  end do

! Now for the loop making the new_loc array

  cnt=1                    ! cnt is the counter through the existing block list on a processor
  new_cnt(:)=1             ! new_cnt is an array of counters for the new block lists on all processors

  do g=1, lrefine_max

!     do proc = 0, noprocs-1     ! proc is the one doing the checking to see who has sent or received
!        print *,'new_cnt',pid,proc,new_cnt(proc)
!     end do

     ! first establish counter increments from lower processors
     do proc = 0, noprocs-1     ! proc is the one doing the checking to see who has sent or received
        sentdown = 0
        if (proc.gt.0) then
           do p = 0, proc-1      ! p is the one getting or receiving the messages
!           do p = 0, pid-1      ! p is the one getting or receiving the messages
              if (bkrecv(g+p*lrefine_max, proc).gt.0) then
!                 new_cnt(proc) = new_cnt(proc) + bkrecv(g+p*lrefine_max, proc) - bksend(g+p*lrefine_max, proc)
!                 new_cnt(proc) = new_cnt(proc) + bkrecv(g+p*lrefine_max, proc)
              endif
              ! renumber our blocks being moved appropriately
              if (bksend(g+p*lrefine_max, proc).gt.0) then
                 do b = 1, bksend(g+p*lrefine_max, proc)
                    if (proc.eq.pid) then
                       new_loc(1, reordered(cnt)) = new_cnt(p)
                       new_loc(2, reordered(cnt)) = p
                       cnt = cnt + 1
                    endif
                    new_cnt(p) = new_cnt(p) + 1
                    sentdown = sentdown + 1
                 end do
              endif 
           end do
        end if
!     end do        ! end loop over processors

     ! second do the work for processor 'proc'  -  two cases proc==pid or otherwise
!     do proc = 0, noprocs-1
        ! discount the blocks that are being sent elsewhere - down already done 
        ! to do this count how many
        sentup = 0
        if (proc.lt.noprocs-1) then
           do p = proc+1, noprocs-1
              sentup = sentup + bksend(g+p*lrefine_max, proc)
           end do
        endif

        ! renumber the blocks staying here
        do b = 1+sentdown, origblocks(g+proc*lrefine_max)-sentup
!           print*, b, cnt, new_cnt(proc), origblocks(g+proc*lrefine_max)
           if (proc.eq.pid) then
              new_loc(1, reordered(cnt)) = new_cnt(proc)
              new_loc(2, reordered(cnt)) = proc
              cnt = cnt + 1
           endif        
           new_cnt(proc) = new_cnt(proc) + 1
!           if (g.eq.1) print *,pid,proc,'added'
!           if (g.eq.2) print *,pid,proc,'sent',sentdown, sentup, origblocks(g+proc*lrefine_max)
        end do

!     end do        ! end loop over processors

     ! finally establish counter increments from higher processors
!     do proc = 0, noprocs-1
!        do p = pid+1, noprocs-1
        do p = proc+1, noprocs-1
           if (bkrecv(g+p*lrefine_max, proc).gt.0) then
!              new_cnt(proc) = new_cnt(proc) + bkrecv(g+p*lrefine_max, proc) - bksend(g+p*lrefine_max, proc)
!              new_cnt(proc) = new_cnt(proc) + bkrecv(g+p*lrefine_max, proc)
           endif
           ! renumber our blocks being moved appropriately
           if (bksend(g+p*lrefine_max, proc).gt.0) then
              do b = 1, bksend(g+p*lrefine_max, proc)
                 if (proc.eq.pid) then
                    new_loc(1, reordered(cnt)) = new_cnt(p)
                    new_loc(2, reordered(cnt)) = p
                    cnt = cnt + 1
                 endif
                 new_cnt(p) = new_cnt(p) + 1
              end do
           endif 
       end do
     end do        ! end loop over processors
  end do           ! end loop over grids

  new_lnblocks = new_cnt(pid)-1
           

!  print *,'Coming out of ceg_setup_newloc'
           
  return
end subroutine ceg_setup_newloc
