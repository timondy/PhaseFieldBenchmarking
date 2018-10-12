!#include "paramesh_preprocessor.fh"


subroutine amr_mg_control(mype,nprocs)
! Wrapper for controlling whole solver process
  integer, intent(in) :: mype, nprocs

  call amr_mg_time_dep_control(mype,nprocs)
end subroutine amr_mg_control

subroutine amr_mg_time_dep_control(mype,nprocs)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use solution_parameters
  use paramesh_interfaces!, only :amr_guardcell, & 
       !&                         amr_prolong, & 
       !&                         amr_restrict
  implicit none
  include 'mpif.h'
  integer, intent(in) :: mype, nprocs
  integer :: back_step_ct,time_step,i_step,i_step_total,ierr,level,iopt, lb
  real defect,end_time,old_defect_1,old_defect_2
  integer output_pe, failer, finegridblocks(2)
  real (KIND=4), dimension(2) :: tarray
  real (KIND=4) :: tresult
  logical chk_gd
  character(len=80) :: chk_checkf

!CEG additions for continuation
  INTEGER::chk_restart=0, chk_chkpt, incset
  DOUBLE PRECISION:: chk_t, chk_dt, chk_dtold

!DEFINE BUFFER HOLDS THE COMMAND LINE ARGUMENT
  CHARACTER *100 BUFFER  

  output_pe=0

  end_time=t0+simulation_time
  incset=0
  if(mype.eq.output_pe) write(6,*) 'end time is ',end_time

!GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT
  if (IArgC().gt.0) then
     if (IArgC().eq.2) then
        call getarg(1,BUFFER)
        read(BUFFER,*) chk_restart
        if (chk_restart.eq.1.or.chk_restart.eq.3) then
           call getarg(2,BUFFER)
           read(BUFFER,*) chk_chkpt
        endif
     else 
        if(mype.eq.0)open (unit=99,file="CHK.out")
           read(99,*) chk_chkpt
        if(mype.eq.0)close(99)
     endif

     if (chk_chkpt.gt.0) then
        if(mype.eq.output_pe) write(6,*) 'Checkpoint on: ', chk_chkpt
        call amr_checkpoint_re(chk_chkpt, user_attr_1=chk_dt, user_attr_2=chk_t, user_attr_3=chk_dtold)
! For Tim
!        call amr_checkpoint_re(chk_chkpt, user_attr_1=chk_dt, user_attr_2=chk_t)
        if(mype.eq.output_pe) write(6,*) 'Checkpoint read at time ', chk_t
        if (chk_restart.eq.3) then
          if(mype.eq.output_pe) write(6,*) 'Produce HDF5 file for Visit'
          call amr_plotfile_chombo(chk_chkpt)
        endif
        if(mype.eq.output_pe) write(6,*) '---------------------------------------'
        call amr_multigrid_block_types(mype,nprocs, finegridblocks)
!        call phase_field_check(1)
        if(mype.eq.output_pe) write(6,*) '---------------------------------------'

     else
        if(mype.eq.output_pe) write(6,*) ' No file found!'
        return
     endif
  else
     if(mype.eq.output_pe) write(6,*) IArgC(), ' command line arguments not found - not 2'
     if(mype.eq.output_pe) write(6,*) ' No file found!'
     return
  endif

  ! Variable setup
  iopt = 1
  print *,mype,iopt,nguard
!  call amr_multigrid_child_set

  allocate(has_children(1:maxblocks),stat=time_step)!Using time_step as dud var
  allocate(has_gchildren(1:maxblocks),stat=time_step)!Using time_step as dud var
  allocate(nodetype_copy(1:maxblocks),stat=time_step)!Using time_step as dud var
  call amr_guardcell(mype,iopt,nguard)
  call amr_mg_init()
  call amr_multigrid_block_types(mype,nprocs, finegridblocks)
 call amr_guardcell(mype, iopt, nguard)

  ! Establish radii, etc
  call phase_field_check(mype, chk_chkpt)
!  call amr_multigrid_block_types(mype, nprocs, finegridblocks)

  ! Now do viz bits
  if (ndim.eq.3) then
     call ceg_loop_over_blocks(chk_chkpt, mype, nprocs)
     call ceg_viz2d(chk_chkpt, mype, nprocs)

     call ceg_viz2dpyr(chk_chkpt, mype, nprocs)
  else
     call ceg_viz2d(chk_chkpt, mype, nprocs)

     call ceg_viz2dpyr(chk_chkpt, mype, nprocs)
  endif

  deallocate(has_children)
  deallocate(has_gchildren)
  deallocate(nodetype_copy)

end subroutine amr_mg_time_dep_control


subroutine ceg_loop_over_blocks(stepno, mype, noprocs)
  ! Not really a user edit but since it's specific to the problem
  ! being simulated it's sat here.
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use paramesh_interfaces
  use solution_parameters
  use time_dep_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: stepno, mype, noprocs
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, result, Nnodes, Ntets, icnt, jcnt, Gnodes, Gtets, Itmp, Itmp2, Itmp3, Itmp4, p
  double precision :: phase,umin,umax,psimin,psimax,u,psi,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dy,dz,p2,p3,p4,phi_val_1,phi_val_2,phi_val_3,phi_val_4,cur_rad=0.0,out_val_1,out_val_2,out_val_3
  double precision :: new_loc, new_rad, temp2, Dtmp, Dtmp2, Dtmp3
  logical :: on_axisx=.false.,on_axisy=.false.,on_axisz=.false.
  character (len=80) :: filename
  integer :: NdProcs(0:noprocs), NdCnts(0:noprocs), NelemProcs(0:noprocs)
  !print *,"Guardcell axis node count:",axis_node_count
  !axis_node_count=0
  phase=0.0
  umin=1e30
  umax=-1e30
  psimin=1e30
  psimax=-1e30

  Nnodes=0
  Ntets=0

  if(mype.eq.0) then
     print *,"Extractor "
  endif

  write(filename,1000) "/tmp/xaxis", stepno,"_",mype
  print *,filename
  open(unit=71,file=filename,status='unknown')
  write(filename,1000) '/tmp/yaxis', stepno,"_",mype
  print *,filename
  open(unit=72,file=filename)
  write(filename,1000) '/tmp/zaxis', stepno,"_",mype
  print *,filename
  open(unit=73,file=filename)

  if (mype.eq.0) then
     write(filename,1003) '/tmp/isosurface', stepno, '.pyr'
     print *,filename
     open(unit=76,file=filename)
     write(filename,1003) '/tmp/isoelem', stepno, '.lat' 
     open(unit=78,file=filename)
     write(filename,1003) '/tmp/isonds', stepno, '.lat'
     open(unit=77,file=filename)
  endif

  write(filename,1001) '/tmp/isoels', stepno, "_",mype, '.con'
  print *,filename
  open(unit=81,file=filename)

  write(filename,1001) '/tmp/isoels', stepno, "_",mype, '.pos'
  print *,filename
  open(unit=82,file=filename)

  write(filename,1001) '/tmp/isonds', stepno, "_",mype, '.pos'
  print *,filename
  open(unit=83,file=filename)

  write(filename,1001) '/tmp/isonds', stepno, "_",mype, '.sol'
  print *,filename
  open(unit=84,file=filename)


! Star by counting nodes and tets that will go into the pyramid

  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then

        call testFaces(lb, result)

        if (result.eq.1) then

            call countCubes(lb, Nnodes, Ntets)

        endif 

     end if
  end do

  call MPI_Reduce(Nnodes, Gnodes, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(Ntets, Gtets, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gather(Nnodes, 1, MPI_INTEGER, NdProcs(0), 1,  MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gather(Ntets, 1, MPI_INTEGER, NelemProcs(0), 1,  MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  if (mype.eq.0) then
     NdCnts(0) = 0

     do p=1, noprocs-1
        NdCnts(p) = NdProcs(p-1) + NdCnts(p-1)
     end do


!     print*,Nnodes, Ntets, Gnodes, Gtets, NdProcs(0)
     call pyrHdrs(stepno, Gnodes, Gtets)

  endif

  call MPI_Scatter(NdCnts(0), 1, MPI_INTEGER, Nnodes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        dz = bsize(3,lb)/real(nzb)
        ! First find out if we're "on axis"
        ! how do we know what the "axis" is?
        ! Simple the nucleate_(x,y,z) parameter holds the co-ords
        ! To be on an axis 2 out of 3 of the block faces must match these co-ords
        ! Only testing lower co-ords so for corner nucleation the nucleus MUST BE at min co-ords for each dimension
        on_axisx=.false.
        on_axisy=.false.
        on_axisz=.false.
        
        if(bnd_box(1,1,lb).eq.nucleate_y)then
           ! lower x bnd matches
           ! yz plane
           on_axisx=.true.
        end if
        if(bnd_box(1,2,lb).eq.nucleate_y)then
           ! lower y bnd matches
           ! xz plane
           on_axisy=.true.
        end if
        if(bnd_box(1,3,lb).eq.nucleate_y)then
           ! lower z bnd matches
           ! xy plane
           on_axisz=.true.
        end if
        j=jl_bnd+nguard
        k=kl_bnd+nguard

        if(on_axisx)then
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
!                 do i=il_bnd+nguard,iu_bnd-nguard
                  i=il_bnd+nguard
                    write(71,1002) bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dy, bnd_box(1,3,lb)+(k-1.5)*dz,  &
                                   unk(1,i,j,k,lb), unk(6,i,j,k,lb), unk(11,i,j,k,lb)
              end do
           end do
        end if

        if(on_axisy)then
           do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
               j=jl_bnd+nguard*k2d
!              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                 do i=il_bnd+nguard,iu_bnd-nguard
                    write(72,1002) bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dy, bnd_box(1,3,lb)+(k-1.5)*dz,  &
                                   unk(1,i,j,k,lb), unk(6,i,j,k,lb), unk(11,i,j,k,lb)
              end do
           end do
        end if

        if(on_axisz)then
           k=kl_bnd+nguard*k3d
!            do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
                 do i=il_bnd+nguard,iu_bnd-nguard
                    write(73,1002) bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dy, bnd_box(1,3,lb)+(k-1.5)*dz,  &
                                   unk(1,i,j,k,lb), unk(6,i,j,k,lb), unk(11,i,j,k,lb)
              end do
           end do
        end if

        call testFaces(lb, result)

        if (result.eq.1) call plotNodes(lb, Nnodes, Ntets)

     end if
  end do

  ! Collate files for .lat and .pyr
  close(81)
  close(82)
  close(83)
  close(84)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (mype.eq.0) then
! Connections into pyramid
     do p=0, noprocs-1
        write(filename,1001) '/tmp/isoels', stepno, "_", p, '.con'
        print *,filename
        open(unit=81,file=filename)
!        print *,p,NelemProcs(p)
        do i = 1, NelemProcs(p)
           read(81,*) Itmp, Itmp2, Itmp3, Itmp4
           write(76,*) Itmp, Itmp2, Itmp3, Itmp4
        end do
        close(81)
     end do
     ! Finish isosurfaceXXXX.pyr
     write(filename,1003) '/tmp/isoelem', stepno, '.lat' 
     write(76,'(A,A)') "include ", filename


! Nodes lattice
     do p=0, noprocs-1
        write(filename,1001) '/tmp/isonds', stepno, "_", p, '.pos'
        print *,filename
        open(unit=83,file=filename)
        do i = 1, NdProcs(p)
           read(83,*) Dtmp, Dtmp2, Dtmp3
           write(77,*) Dtmp, Dtmp2, Dtmp3
        end do
        close(83)
     end do
     do p=0, noprocs-1
        write(filename,1001) '/tmp/isonds', stepno, "_", p, '.sol'
        print *,filename
        open(unit=84,file=filename)
        do i = 1, NdProcs(p)
           read(84,*) Dtmp
           write(77,*) Dtmp
        end do
        close(84)
     end do

! Elements lattice
     do p=0, noprocs-1
        write(filename,1001) '/tmp/isoels', stepno, "_", p, '.pos'
        print *,filename
        open(unit=82,file=filename)
        do i = 1, NelemProcs(p)
           read(82,*) Dtmp, Dtmp2, Dtmp3
           write(78,*) Dtmp, Dtmp2, Dtmp3
        end do
        close(82)
     end do



  endif

  close(71)
  close(72)
  close(73)
  close(76)
  close(77)
  close(78)

1000 format(A15,I0.4,A1,I0.2)    ! axis files
1001 format(A15,I0.4,A1,I0.2,A4) ! per processor part of lat/pyr
1003 format(A15,I0.4,A4)         ! lat and pyr files

1002 format (F10.4,X,F10.4,X,F10.4,X,F10.4,X,F10.4,X,F10.4) 

end subroutine ceg_loop_over_blocks

subroutine testFaces(lb, result)
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use paramesh_interfaces
  use solution_parameters
  use time_dep_parameters
  implicit none
  integer, intent(inout) :: lb, result
  integer :: i,j,k
  double precision :: psimin,psimax

  psimin=1e30
  psimax=-1e30

! Do the min/max from first guard to last real

  do k=kl_bnd, ku_bnd-nguard*k3d+1
     do j=jl_bnd, ju_bnd-nguard*k2d +1
        i=il_bnd
        psimin = MIN(psimin, unk(1,i,j,k,lb))
        psimax = MAX(psimax, unk(1,i,j,k,lb))

        i=iu_bnd-nguard+1
        psimin = MIN(psimin, unk(1,i,j,k,lb))
        psimax = MAX(psimax, unk(1,i,j,k,lb))
     end do
  end do

  do k=kl_bnd, ku_bnd-nguard*k3d +1
     do i=il_bnd, iu_bnd-nguard+1
        j=jl_bnd
        psimin = MIN(psimin, unk(1,i,j,k,lb))
        psimax = MAX(psimax, unk(1,i,j,k,lb))

        j=ju_bnd-nguard*k2d+1
        psimin = MIN(psimin, unk(1,i,j,k,lb))
        psimax = MAX(psimax, unk(1,i,j,k,lb))
     end do
  end do

  do j=jl_bnd, ju_bnd-nguard*k2d +1
     do i=il_bnd, iu_bnd-nguard+1
        k=kl_bnd
        psimin = MIN(psimin, unk(1,i,j,k,lb))
        psimax = MAX(psimax, unk(1,i,j,k,lb))

        k=ku_bnd-nguard*k3d+1
        psimin = MIN(psimin, unk(1,i,j,k,lb))
        psimax = MAX(psimax, unk(1,i,j,k,lb))
     end do
  end do

  if (psimin.le.0.0.and.psimax.ge.0.0) then
     result = 1
  else
     result = 0
  endif

end subroutine testFaces



subroutine countCubes(lb, Nnodes, Ntets)
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use paramesh_interfaces
  use solution_parameters
  use time_dep_parameters
  implicit none
  integer, intent(inout) :: lb, Nnodes, Ntets
  integer :: i,j,k
  double precision :: psimin,psimax

!  do k=kl_bnd, ku_bnd-nguard*k3d-1
!     do j=jl_bnd, ju_bnd-nguard*k2d-1
!        do i=il_bnd, iu_bnd-nguard-1
  do k=kl_bnd, ku_bnd-nguard*k3d
     do j=jl_bnd, ju_bnd-nguard*k2d
        do i=il_bnd, iu_bnd-nguard
           psimin=1.0
           psimax=-1.0

           psimin = MIN(psimin, unk(1,i,j,k,lb))
           psimin = MIN(psimin, unk(1,i+1,j,k,lb))
           psimin = MIN(psimin, unk(1,i+1,j+1,k,lb))
           psimin = MIN(psimin, unk(1,i,j+1,k,lb))
           psimin = MIN(psimin, unk(1,i,j,k+1,lb))
           psimin = MIN(psimin, unk(1,i+1,j,k+1,lb))
           psimin = MIN(psimin, unk(1,i+1,j+1,k+1,lb))
           psimin = MIN(psimin, unk(1,i,j+1,k+1,lb))

           psimax = MAX(psimax, unk(1,i,j,k,lb))
           psimax = MAX(psimax, unk(1,i+1,j,k,lb))
           psimax = MAX(psimax, unk(1,i+1,j+1,k,lb))
           psimax = MAX(psimax, unk(1,i,j+1,k,lb))
           psimax = MAX(psimax, unk(1,i,j,k+1,lb))
           psimax = MAX(psimax, unk(1,i+1,j,k+1,lb))
           psimax = MAX(psimax, unk(1,i+1,j+1,k+1,lb))
           psimax = MAX(psimax, unk(1,i,j+1,k+1,lb))

           if (psimin.le.0.0.and.psimax.ge.0.0) then
              Nnodes = Nnodes + 8
              Ntets = Ntets + 5
           endif
        end do
     end do
  end do

end subroutine countCubes     

subroutine plotNodes(lb, Nnodes, Ntets)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use solution_parameters
  integer, intent(inout) :: lb,  Nnodes, Ntets
  integer :: i, j, k, icnt, jcnt, ii, jj, kk
  double precision :: dx,dy,dz,psimin,psimax

  dx = bsize(1,lb)/real(nxb)
  dy = bsize(2,lb)/real(nyb)
  dz = bsize(3,lb)/real(nzb)


!  do k=kl_bnd, ku_bnd-nguard*k3d-1
!     do j=jl_bnd, ju_bnd-nguard*k2d-1
!        do i=il_bnd, iu_bnd-nguard-1
  do k=kl_bnd, ku_bnd-nguard*k3d
     do j=jl_bnd, ju_bnd-nguard*k2d
        do i=il_bnd, iu_bnd-nguard
           psimin=1d30
           psimax=-1d30
           psimin = MIN(psimin, unk(1,i,j,k,lb))
           psimin = MIN(psimin, unk(1,i+1,j,k,lb))
           psimin = MIN(psimin, unk(1,i+1,j+1,k,lb))
           psimin = MIN(psimin, unk(1,i,j+1,k,lb))
           psimin = MIN(psimin, unk(1,i,j,k+1,lb))
           psimin = MIN(psimin, unk(1,i+1,j,k+1,lb))
           psimin = MIN(psimin, unk(1,i+1,j+1,k+1,lb))
           psimin = MIN(psimin, unk(1,i,j+1,k+1,lb))

           psimax = MAX(psimax, unk(1,i,j,k,lb))
           psimax = MAX(psimax, unk(1,i+1,j,k,lb))
           psimax = MAX(psimax, unk(1,i+1,j+1,k,lb))
           psimax = MAX(psimax, unk(1,i,j+1,k,lb))
           psimax = MAX(psimax, unk(1,i,j,k+1,lb))
           psimax = MAX(psimax, unk(1,i+1,j,k+1,lb))
           psimax = MAX(psimax, unk(1,i+1,j+1,k+1,lb))
           psimax = MAX(psimax, unk(1,i,j+1,k+1,lb))

           if (psimin.le.0.0.and.psimax.ge.0.0) then
              do kk = k, k+1
                 do jj = j, j+1
                    do ii = i, i+1
                       write(83,1003) bnd_box(1,1,lb)+(ii-1.5)*dx, bnd_box(1,2,lb)+(jj-1.5)*dy, bnd_box(1,3,lb)+(kk-1.5)*dz
                       write(84,1004) unk(1,ii,jj,kk,lb)
                    end do
                 end do
              end do
              call tetConnWrite(lb, i, j, k, Nnodes)
              Nnodes = Nnodes + 8
           endif
        end do
     end do
  end do

1003 format (F10.4,X,F10.4,X,F10.4)
1004 format (F10.4)

end subroutine plotNodes

subroutine tetConnWrite(lb, i, j, k, Nnodes)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions
  use time_dep_parameters
  use multigrid_parameters
  use solution_parameters
  integer, intent(in) :: lb, i, j, k
  integer, intent(inout) :: Nnodes
  integer :: n0, n1, n2, n3, n4, n5, n6, n7, icnt, jcnt

  n0 = Nnodes + 0 
  n1 = Nnodes +  1 
  n2 = Nnodes +    3
  n3 = Nnodes +   2 
  n4 = Nnodes +     4
  n5 = Nnodes +      5
  n6 = Nnodes +        7
  n7 = Nnodes +       6

! Tets 0125, 2675, 0475, 0237, 0275
  write(81,1005) n0, n1, n4, n3
  write(81,1005) n1, n6, n5, n4
  write(81,1005) n3, n6, n4, n7
  write(81,1005) n1, n6, n4, n3
  write(81,1005) n1, n2, n6, n3

  write(82,1007) i+0.0, j+0.1, k+0.2
  write(82,1007) i+0.3, j+0.4, k+0.5
  write(82,1007) i+0.6, j+0.7, k+0.8
  write(82,1007) i+0.9, j+0.0, k+0.1
  write(82,1007) i+0.2, j+0.3, k+0.4

1005 format (5I,X,5I,X,5I,X,5I)
1007 format (F10.4,X,F10.4,X,F10.4)

end subroutine tetConnWrite


subroutine pyrHdrs(stepno, Nnodes, Ntets)

  integer, intent(in) :: stepno, Nnodes, Ntets
  integer :: i
  character (len=80) :: fname1, fname2, fname3, fname4, fname5, fname6

  write(fname1,1001) '/tmp/isonds', stepno, '.lat'
  write(fname3,1001) '/tmp/isoelem', stepno, '.lat' 

  write(76,'(A)') "#!/usr/explorer/bin/explorer cxPyramid plain 1.0"
  write(76,'(A,A)') "include ", fname1
  write(76,'(A)') " 3"            !                                     /* NUMBER OF PYRAMID LAYERS */
  write(76,'(A)') " 1"            !                                     /* COMPRESSION TYPE: 1=unique */
  write(76,'(A)') " 4"            !                                     /* COMPRESSION TYPE: 4=tetrahedra */
  write(76,'(A)') " 0 0"          !                                     /* LAYER 1: COMPRESSED */
  write(76,'(A)') " 0 0"          !                                     /* LAYER 2: COMPRESSED */
  write(76,*) Ntets, 4*Ntets          !                                     /* LAYER 2: COMPRESSED */
  do i = 0, Ntets
     write(76,*) 4*i
  end do

  write(77,'(A)') "#!/usr/explorer/bin/explorer cxLattice plain 1.0"
  write(77,'(A)') " 1"    !  /* DIMENSIONS */"
  write(77,*) Nnodes
  write(77,'(A)') " 1"      !/* NUMBER OF DATA VALUES PER ELEMENT */"
  write(77,'(A)') " 4"      !/* DATA FORMAT: 3=float, 4=double */"
  write(77,'(A)') " 2"      !/* COORDINATE TYPE: 2=curvilinear */"
  write(77,'(A)') " 1"      !/* NUMBER OF DATA SETS PRESENT */"
  write(77,'(A)') " 3"      !/* COORDINATE DIMENSIONS */"

  write(78,'(A)') "#!/usr/explorer/bin/explorer cxLattice plain 1.0"
  write(78,'(A)') " 1"      !/* DIMENSIONS */"
! could be Ntets
!  write(78,'(I7)') 0!,       /* NUMBER OF DATA POINTS PER DIMENSION */"
  write(78,'(I7)') Ntets!,       /* NUMBER OF DATA POINTS PER DIMENSION */"
  write(78,'(A)') " 0"!      /* NUMBER OF DATA VALUES PER ELEMENT */"
  write(78,'(A)') " 4"!      /* DATA FORMAT: 3=float, 4=double */"
  write(78,'(A)') " 2"!      /* COORDINATE TYPE: 2=curvilinear */"
  write(78,'(A)') " 1"!      /* NUMBER OF DATA SETS PRESENT */"
  write(78,'(A)') " 3"!      /* COORDINATE DIMENSIONS */"


1001 format(A15,I0.4,A4)


end subroutine pyrHdrs


subroutine ceg_viz2d(stepno, mype, noprocs)
  ! Not really a user edit but since it's specific to the problem
  ! being simulated it's sat here.
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use paramesh_interfaces
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters


  implicit none
  include 'mpif.h'
  integer, intent(in) :: stepno, mype, noprocs
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, result, Nnodes, Ntets, icnt, jcnt, Gnodes, Gtets, Itmp, Itmp2, Itmp3, Itmp4, p
  double precision :: phase,umin,umax,psimin,psimax,u,psi,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dy,dz,p2,p3,p4,phi_val_1,phi_val_2,phi_val_3,phi_val_4,cur_rad=0.0,out_val_1,out_val_2,out_val_3
  double precision :: new_loc, new_rad, temp2, Dtmp, Dtmp2, Dtmp3
  logical :: on_axisx=.false.,on_axisy=.false.,on_axisz=.false.
  character (len=80) :: filename
  integer :: NdProcs(0:noprocs), NdCnts(0:noprocs), NelemProcs(0:noprocs)
  !print *,"Guardcell axis node count:",axis_node_count
  !axis_node_count=0


  write(filename,1000) "/tmp/2dsheet1", stepno,"_",mype
  print *,filename
  open(unit=71,file=filename,status='unknown')
  write(filename,1000) "/tmp/2dsheet2", stepno,"_",mype
  print *,filename
  open(unit=72,file=filename,status='unknown')
	
  if (total_vars.eq.3) then
     write(filename,1000) "/tmp/2dsheet3", stepno,"_",mype
     print *,filename
     open(unit=73,file=filename,status='unknown')
  endif

  icnt=0
  do lb=1,lnblocks
     if(nodetype(lb).eq.1.and.(ndim.eq.2.or.(bnd_box(1,3,lb).lt.1.0e-1.and.bnd_box(2,3,lb).gt.-1.0e-1))) then
!     if(nodetype(lb).eq.1) then
        icnt = icnt+(iu_bnd-nguard - (il_bnd+nguard) +1)*(ju_bnd-nguard*k2d-(jl_bnd+nguard*k2d)+1)
     endif
  end do
!  print *,iu_bnd, nguard, il_bnd
!  print *,(iu_bnd-nguard - (il_bnd+nguard) +1),(ju_bnd-nguard*k2d-(jl_bnd+nguard*k2d)+1)

  write(71,*), icnt
  write(72,*), icnt
  if (total_vars.eq.3) write(73,*), icnt
  do lb=1,lnblocks
     if(nodetype(lb).eq.1.and.(ndim.eq.2.or.(bnd_box(1,3,lb).lt.1.0e-1.and.bnd_box(2,3,lb).gt.-1.0e-1))) then
        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        ! First find out if we're "on axis"
        ! how do we know what the "axis" is?
        ! Simple the nucleate_(x,y,z) parameter holds the co-ords
        ! To be on an axis 2 out of 3 of the block faces must match these co-ords
        ! Only testing lower co-ords so for corner nucleation the nucleus MUST BE at min co-ords for each dimension

        Dtmp = bnd_box(1,1,lb)
        Dtmp2 = bnd_box(1,2,lb)
        k=1
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
           do i=il_bnd+nguard,iu_bnd-nguard
              write(71,1002) Dtmp+(i-1.5)*dx, Dtmp2+(j-1.5)*dy, &
                                   unk(1,i,j,k,lb)
              write(72,1002) Dtmp+(i-1.5)*dx, Dtmp2+(j-1.5)*dy, &
                                   unk(6,i,j,k,lb)
              if (total_vars.eq.3) write(73,1002) Dtmp+(i-1.5)*dx, Dtmp2+(j-1.5)*dy, &
                                   unk(11,i,j,k,lb)
!                    write(71,1003) bnd_box(1,1,lb)+(i-1.5)*dx, bnd_box(1,2,lb)+(j-1.5)*dy, &
!                                   unk(1,i,j,k,lb), unk(6,i,j,k,lb), unk(11,i,j,k,lb)
           end do
           write(71,*) ''
           write(72,*) ''
           if (total_vars.eq.3) write(73,*) ''
        end do
        write(71,*) ''
        write(72,*) ''
        if (total_vars.eq.3) write(73,*) ''
     endif
  end do

  close(71)
  close(72)
  close(73)

1000 format(A15,I0.4,A1,I0.2)    ! axis files
1002 format (F10.4,X,F10.4,X,F10.4) 
1003 format (F10.4,X,F10.4,X,F10.4,X,F10.4,X,F10.4) 

end subroutine ceg_viz2d

subroutine ceg_viz2dpyr(stepno, mype, noprocs)
  ! Not really a user edit but since it's specific to the problem
  ! being simulated it's sat here.
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_comm_data ! Needed for definition of amr_mpi_real
  use paramesh_interfaces
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters


  implicit none
  include 'mpif.h'
  integer, intent(in) :: stepno, mype, noprocs
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, result, Nnodes, Ntets, icnt, jcnt, Gnodes, Gtets, Itmp, Itmp2, Itmp3, Itmp4, p, chan, bcnt
  double precision :: phase,umin,umax,psimin,psimax,u,psi,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dy,dz,p2,p3,p4,phi_val_1,phi_val_2,phi_val_3,phi_val_4,cur_rad=0.0,out_val_1,out_val_2,out_val_3
  double precision :: new_loc, new_rad, temp2, Dtmp, Dtmp2, Dtmp3
  logical :: on_axisx=.false.,on_axisy=.false.,on_axisz=.false.
  character (len=80) :: filename
  integer :: NdProcs(0:noprocs), NdCnts(0:noprocs), NelemProcs(0:noprocs)
  logical :: blockview = .false.
  !print *,"Guardcell axis node count:",axis_node_count
  !axis_node_count=0

! First count how many blocks and points we have
  bcnt=0
  icnt=0
  do lb=1,lnblocks
     if(nodetype(lb).eq.1.and.(ndim.eq.2.or.(bnd_box(1,3,lb).lt.1.0e-1.and.bnd_box(2,3,lb).gt.-1.0e-1))) then
!     if(nodetype(lb).eq.1) then
        bcnt = bcnt+1
        icnt = icnt+(iu_bnd-nguard - (il_bnd+nguard) +1)*(ju_bnd-nguard*k2d-(jl_bnd+nguard*k2d)+1)
     endif
  end do
!  print *,iu_bnd, nguard, il_bnd
!  print *,(iu_bnd-nguard - (il_bnd+nguard) +1),(ju_bnd-nguard*k2d-(jl_bnd+nguard*k2d)+1)



! Open relevant lattice files
  write(filename,1000) "/tmp/2d-data", stepno,"_",mype,'.lat'
  print *,filename
  open(unit=81,file=filename,status='unknown')

  chan=81
  write(chan,'(A)') "#!/usr/explorer/bin/explorer cxLattice plain 1.0"
  write(chan,'(A)') " 1"    !  /* DIMENSIONS */"
  write(chan,*) bcnt*nxb*nyb
  if (total_vars.eq.3) then
     write(chan,'(A)') " 3"      !/* NUMBER OF DATA VALUES PER ELEMENT  */"
  else
     write(chan,'(A)') " 2"      !/* NUMBER OF DATA VALUES PER ELEMENT  */"
  endif
  write(chan,'(A)') " 4"      !/* DATA FORMAT: 3=float, 4=double */"
  write(chan,'(A)') " 2"      !/* COORDINATE TYPE: 2=curvilinear */"
  write(chan,'(A)') " 1"      !/* NUMBER OF DATA SETS PRESENT */"
  write(chan,'(A)') " 2"      !/* COORDINATE DIMENSIONS */"


! Open Element file
  write(filename,1001) '/tmp/2d-elem', stepno, '.lat' 
  open(unit=82,file=filename)

  write(82,'(A)') "#!/usr/explorer/bin/explorer cxLattice plain 1.0"
  write(82,'(A)') " 1"      !/* DIMENSIONS */"
  ! could be Ntets
  !  write(78,'(I7)') 0!,       /* NUMBER OF DATA POINTS PER DIMENSION */"
  write(82,'(I7)') bcnt*(nxb-1)*(nyb-1)  !,       /* NUMBER OF DATA POINTS PER DIMENSION */"
  write(82,'(A)') " 1"!      /* NUMBER OF DATA VALUES PER ELEMENT */"
  write(82,'(A)') " 4"!      /* DATA FORMAT: 3=float, 4=double */"
  write(82,'(A)') " 2"!      /* COORDINATE TYPE: 2=curvilinear */"
  write(82,'(A)') " 1"!      /* NUMBER OF DATA SETS PRESENT */"
  write(82,'(A)') " 2"!      /* COORDINATE DIMENSIONS */"

! Open and fill pyramid file
  write(filename,1001) '/tmp/planeview', stepno, '.pyr'
  print *,filename
  open(unit=86,file=filename)

  write(86,'(A)') "#!/usr/explorer/bin/explorer cxPyramid plain 1.0"
  write(filename,1000) "/tmp/2d-data", stepno,"_",mype,'.lat'
  write(86,'(A,A)') "include ", filename
  write(86,'(A)') " 2"            !                                     /* NUMBER OF PYRAMID LAYERS */
  write(86,'(A)') " 1"            !                                     /* COMPRESSION TYPE: 1=unique */
  write(86,'(A)') " 3"            !                                     /* COMPRESSION TYPE: 4=tetrahedra */
  write(86,'(A)') " 0 0"          !                                     /* LAYER 1: COMPRESSED */
  write(86,*) bcnt*(nxb-1)*(nyb-1), 4*bcnt*(nxb-1)*(nyb-1)  !           /* LAYER 2: COMPRESSED */
  do i = 0, bcnt*(nxb-1)*(nyb-1)
     write(86,*) 4*i
  end do

  bcnt=0
  icnt=0
  do lb=1,lnblocks
     if(nodetype(lb).eq.1.and.(ndim.eq.2.or.(bnd_box(1,3,lb).lt.1.0e-1.and.bnd_box(2,3,lb).gt.-1.0e-1))) then
        if (blockview) then
           dx = bsize(1,lb)/real(nxb)
           dy = bsize(2,lb)/real(nyb)
        else
           dx = bsize(1,lb)/real(nxb-1)
           dy = bsize(2,lb)/real(nyb-1)
        endif

        Dtmp = bnd_box(1,1,lb)
        Dtmp2 = bnd_box(1,2,lb)
        k=1
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
           do i=il_bnd+nguard,iu_bnd-nguard
              if (blockview) then
                 write(81,1005) Dtmp+(i-1.5)*dx, Dtmp2+(j-1.5)*dy
              else
                 write(81,1005) Dtmp+(i-2)*dx, Dtmp2+(j-2)*dy
              endif
           end do
        end do
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d-1
           do i=il_bnd+nguard,iu_bnd-nguard-1
              if (blockview) then
                 write(82,1005) Dtmp+(i-1.0)*dx, Dtmp2+(j-1.0)*dy
              else
                 write(82,1005) Dtmp+(i-1.0)*dx, Dtmp2+(j-1.0)*dy
              endif
           end do
        end do

        bcnt = bcnt+1
        icnt = icnt+(iu_bnd-nguard - (il_bnd+nguard) +1)*(ju_bnd-nguard*k2d-(jl_bnd+nguard*k2d)+1)
     endif
  end do

  bcnt=0
  icnt=0
  do lb=1,lnblocks
     if(nodetype(lb).eq.1.and.(ndim.eq.2.or.(bnd_box(1,3,lb).lt.1.0e-1.and.bnd_box(2,3,lb).gt.-1.0e-1))) then
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
           do i=il_bnd+nguard,iu_bnd-nguard
              if (total_vars.eq.3)  then
                 write(81,1002) unk(1,i,j,k,lb), unk(6,i,j,k,lb), unk(11,i,j,k,lb)
              else
                 write(81,1005) unk(1,i,j,k,lb), unk(6,i,j,k,lb)
              endif
!              write(81,*) unk(1,i,j,k,lb)
              write(82,*) bcnt
           end do
        end do
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d-1
           do i=il_bnd+nguard,iu_bnd-nguard-1
              write(86,1004) icnt+(i-nguard)+(j-nguard-1)*nxb-1,       &
                             icnt+(i+1-nguard)+(j-nguard-1)*nxb-1,     &
                             icnt+(i+1-nguard)+(j+1-nguard-1)*nxb-1,   &
                             icnt+(i-nguard)+(j+1-nguard-1)*nxb-1
           end do
        end do

        bcnt = bcnt+1
        icnt = icnt+(iu_bnd-nguard - (il_bnd+nguard) +1)*(ju_bnd-nguard*k2d-(jl_bnd+nguard*k2d)+1)
     endif
  end do

  write(filename,1001) '/tmp/2d-elem', stepno, '.lat' 
  write(86,'(A,A)') "include ", filename

  close(81)
  close(82)
  close(86)



1000 format(A15,I0.4,A1,I0.2,A4)    ! axis files
1001 format(A15,I0.4,A4)
1002 format (F10.4,X,F10.4,X,F10.4) 
1003 format (F10.4,X,F10.4,X,F10.4,X,F10.4,X,F10.4) 
1004 format(I8,2X,I8,2X,I8,2X,I8)
1005 format (F10.4,X,F10.4) 

end subroutine ceg_viz2dpyr
