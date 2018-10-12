subroutine tim_two_checkpoints(pid, noprocs)
  use workspace
  use tree
  use physicaldata
  use paramesh_dimensions  
  use time_dep_parameters  
  use multigrid_parameters
  use paramesh_interfaces
  use generic_parameters
  use checkpoint_parameters
  implicit none
  include 'mpif.h'
  integer,intent(in) :: pid, noprocs
  integer :: ierr, error
  integer :: i,j,k,lb
  double precision, allocatable :: backup_array (:,:,:,:)
  double precision :: diff, diff_sum, found_count

  allocate(backup_array(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:lnblocks),stat=error)

  do lb=1,lnblocks
    do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
      do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
        do i=il_bnd+nguard,iu_bnd-nguard
          backup_array(i,j,k,lb) = unk(1,i,j,k,lb)
        enddo
      enddo
    enddo
  enddo

  chk_chkpt = 1280
  chk_checkf = "hdf5"
  call app_read_checkpoint(pid, noprocs)
  time = chk_t
  dt = chk_dt
  dtold = chk_dtold

  call mpi_amr_global_domain_limits()
  call set_domain_limits()

  refine(1:lnblocks) = .false.
  call amr_refine_derefine

  diff = 0.0
  diff_sum = 0.0
  found_count = 0.0
  do lb=1,lnblocks
    if(lrefine(lb).eq.lrefine_max)then
      do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
        do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
          do i=il_bnd+nguard,iu_bnd-nguard
            diff = abs(abs(backup_array(i,j,k,lb)) - abs(unk(1,i,j,k,lb)))
            diff_sum = diff_sum+diff*diff
            found_count = found_count + 1.0
          enddo
        enddo
      enddo
    endif
  enddo

  if(pid.eq.0) write(*,*) "-------------------L2 norm----------------------"
  if(pid.eq.0) write(*,*) sqrt((diff_sum)/found_count)
  stop


end subroutine tim_two_checkpoints
