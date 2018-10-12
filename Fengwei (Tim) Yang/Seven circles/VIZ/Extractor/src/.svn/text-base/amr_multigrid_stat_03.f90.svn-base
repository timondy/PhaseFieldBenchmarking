subroutine amr_multigrid_block_types(mype,nprocs, finegridblocks)
!
!    Checks all blocks and sees what ref levels/node types they all are  
!
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: mype,nprocs
  integer :: blocks_type(-1:lrefine_max),blocks_ref(-1:lrefine_max)
  integer :: temp(-1:lrefine_max),recv_buff(1:(lrefine_max+2)*(nprocs))
  integer :: lb,ierr,proc
  integer, intent(out) :: finegridblocks(2)
  do lb=-1,lrefine_max
     temp(lb)=0
     blocks_ref(lb)=0
  end do
  do lb=1,lnblocks
     ! Put all data into temp on each proc, then gather it on proc 1 into the final data structs
     temp(lrefine(lb))=temp(lrefine(lb))+1
  end do
  call MPI_Gather(temp(-1),lrefine_max+2,MPI_INTEGER, &
       recv_buff(1),lrefine_max+2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if(mype.eq.0)then
     call coutput(lrefine_max, nprocs, recv_buff)
     do lb=-1,lrefine_max
        do proc=0,nprocs-1
           blocks_ref(lb)=blocks_ref(lb)+recv_buff(lb+2+proc*(lrefine_max+2))
        end do
     end do    
  end if
  ! Done ref level, now do nodetype
  ! Assuming that nodetype has the same range as lrefine
  do lb=-1,lrefine_max
     temp(lb)=0
     blocks_type(lb)=0
  end do
  do lb=1,lnblocks
     temp(nodetype(lb))=temp(nodetype(lb))+1
  end do
  finegridblocks(2) = 0
  call MPI_Gather(temp(-1),lrefine_max+2,MPI_INTEGER, &
       recv_buff(1),lrefine_max+2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if(mype.eq.0)then
     do lb=-1,lrefine_max
        do proc=0,nprocs-1
           blocks_type(lb)=blocks_type(lb)+recv_buff(lb+2+proc*(lrefine_max+2))
        end do
     end do
     print *,"Ref/Type No. | Blocks by type | Blocks by ref"
     do lb=-1,lrefine_max
        finegridblocks(2) = finegridblocks(2) + blocks_ref(lb)
        print *,lb,blocks_type(lb),blocks_ref(lb)
     end do
     finegridblocks(1) = blocks_ref(lrefine_max)
  end if
end subroutine amr_multigrid_block_types
