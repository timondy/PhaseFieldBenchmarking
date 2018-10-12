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
!    integer, intent(inout) :: new_loc(:,:)
    integer, dimension(2, maxblocks_tr), intent(inout) :: new_loc
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

  !!!!!!!!!!!
  ! Stage 1 : Discover how many blocks each processor has on each level
  !!!!!!!!!!!
  new_loc(:,:)=-1 

!  print *,pid,' In ceg_mesh_partitioning : maxblocks_tr=',maxblocks_tr
  locblocks(:) = 0

  do b=1,lnblocks
     ! Put all data into blocks on each proc, then gather it from othera
     locblocks(lrefine(b)) = locblocks(lrefine(b))+1
  end do

  call MPI_AllGather(locblocks, lrefine_max, MPI_INTEGER, &
       blocks, lrefine_max, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  do b=1, lrefine_max*noprocs
     origblocks(b) = blocks(b)
  end do

  !Vbles: 
  !       locblocks is maxlevels array of how many blocks per grid on this processor
  !       blocks is maxlevels x np of how many currently on each grid on each processor
  !       origblocks is same as blocks

  !!!!!!!!!!!
  ! Stage 2 : BALANCING
  !!!!!!!!!!!
  bksend(:,:) = 0
  bkrecv(:,:) = 0

  ! Loop over grids  to balance each in turn
  do g = 1, lrefine_max

!     if (pid.eq.0) print *,pid, ' On grid ',g
     total = 0
     ! Count how many on this grid in total
     do p = 0, noprocs-1
        total = total + blocks(g+p*lrefine_max)
     end do

     ! Loop over the processors to establish the target block number for each
     do proc=0, noprocs-1
        tgt(proc) = total / (noprocs-proc)
        total = total - tgt(proc)
     end do

     !!!!!!!!!!!!
     ! Stage 2b : How many do we need to move each way for each procesor
     !!!!!!!!!!!!
     ! Arrays bksend and bkrecv filled here
     !        blocks is now the number after redistribution
     do proc=0, noprocs-1
        T = tgt(proc)

        p = g + lrefine_max*proc
        if (blocks(p).eq.T) then
           ! If at the target number then do nothing
           wanted = 0

        else if (blocks(p).gt.T) then
           ! If I have too many blocks
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
           ! If I have too few blocks
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

!  if (pid.eq.0) then
!     print *, '*****'
!        call coutput2(lrefine_max, noprocs, blocks)
!     print *, '*****'
!  end if
 ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

!  print *,'Before new_loc(1)=',new_loc(1,1),new_loc(2,1)
!  do b=1, lrefine_max*noprocs
!     print *, 'orig', origblocks(b), blocks(b)
!  end do


  !!!!!!!!!!!
  ! Stage 3 : Now move blocks
  !!!!!!!!!!!

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

  !!!!!!!!!!!
  ! Stage 4 : Finally set up block_starts array to be right offsets for new, ordered  distribution
  !!!!!!!!!!!
  block_starts(1) = 1
  do g = 1, lrefine_max
    block_starts(g+1) =  block_starts(g)+blocks(g+pid*lrefine_max)
  end do

!  print*,'out of ceg_block_balancing'

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
  do while (.not.escape.and.gfound.lt.lrefine_max)     ! 
     if (lrefine(b).gt.0) then
        if (gstart(lrefine(b)).eq.-1) then
           gstart(lrefine(b)) = b                         ! gstart() is address of first found block at each level
           gfound = gfound + 1                            ! gfound is just to check if we have found all grid levels
        endif
     endif
     b = b + 1
     if (b.gt.lnblocks) escape = .true.
  end do

!  do g = 1, lrefine_max
!     print *, 'g, gstart', g, gstart(g)
!  end do 

  ! form a reordered block list in grid order
  cnt = 1
  do g = 1, lrefine_max                                ! Loop over grids
     do b = 1, origblocks(g+pid*lrefine_max)           !   Loop over target block numbers
        reordered(cnt) = gstart(g)                     !     reordered id is this levels next position
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

!  print *,'Reordered list'
!  do b=1, lnblocks
!     print *,b, reordered(b), lrefine(reordered(b))
!  end do

! Now for the loop making the new_loc array

  cnt=1                    ! cnt is the counter through the existing block list on a processor
  new_cnt(:)=1             ! new_cnt is an array of counters for the new block lists on all processors

  do g=1, lrefine_max      ! loop over grids

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

!  do b = 1, lnblocks
!     print*, b, new_loc(1, b), new_loc(2, b)
!  end do           

!  print *,'Coming out of ceg_setup_newloc'
           
  return
end subroutine ceg_setup_newloc
