!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_initial_soln_blk                                         REQUIRED
!!!!  * Sets the initial solution on a block
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_test_refinement                                          REQUIRED
!!!!  * Application specific refinement tests - unused in phase field code
!!!!  * Only called if local_adaptation=2
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine app_initial_soln_blk(lb)
  
!--------------------------------------------------------------
! include files 
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters

  use workspace
  use paramesh_interfaces
  use paramesh_mpi_interfaces

  implicit none
  include 'mpif.h'
  
  integer, intent(in) :: lb
  integer :: nguard0, i, j, k, mype, iopt, nlayers
  double precision :: pi, dx, xi, yi, temp, pikk=1.0
  double precision :: r, xx, yy, rr
!--------------------------------------------------------------

  pi=4.0*datan(pikk)
  nguard0 = nguard*npgs

! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
    do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
      do i=il_bnd+nguard0,iu_bnd-nguard0

        xi = (bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5))
        yi = (bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5))

        unk(1,i,j,k,lb) = -1.0

        xx = (1.0/4.0)*2.0*pi
        yy = (1.0/4.0)*2.0*pi
        rr = (1.0/10.0)*2.0*pi
        r = sqrt( (xi-xx)**2.0 + (yi-yy)**2.0 )
        if(r-rr.lt.0.0) unk(1,i,j,k,lb) = -1.0+2.0*exp(-ep*ep/((r-rr)*(r-rr)))

        xx = (1.0/8.0)*2.0*pi
        yy = (3.0/8.0)*2.0*pi
        rr = (1.0/15.0)*2.0*pi
        r = sqrt( (xi-xx)**2.0 + (yi-yy)**2.0 )
        if(r-rr.lt.0.0) unk(1,i,j,k,lb) = -1.0+2.0*exp(-ep*ep/((r-rr)*(r-rr)))

        xx = (1.0/4.0)*2.0*pi
        yy = (5.0/8.0)*2.0*pi
        rr = (1.0/15.0)*2.0*pi
        r = sqrt( (xi-xx)**2.0 + (yi-yy)**2.0 )
        if(r-rr.lt.0.0) unk(1,i,j,k,lb) = -1.0+2.0*exp(-ep*ep/((r-rr)*(r-rr)))

        xx = (1.0/2.0)*2.0*pi
        yy = (1.0/8.0)*2.0*pi
        rr = (1.0/20.0)*2.0*pi
        r = sqrt( (xi-xx)**2.0 + (yi-yy)**2.0 )
        if(r-rr.lt.0.0) unk(1,i,j,k,lb) = -1.0+2.0*exp(-ep*ep/((r-rr)*(r-rr)))

        xx = (3.0/4.0)*2.0*pi
        yy = (1.0/8.0)*2.0*pi
        rr = (1.0/20.0)*2.0*pi
        r = sqrt( (xi-xx)**2.0 + (yi-yy)**2.0 )
        if(r-rr.lt.0.0) unk(1,i,j,k,lb) = -1.0+2.0*exp(-ep*ep/((r-rr)*(r-rr)))

        xx = (1.0/2.0)*2.0*pi
        yy = (1.0/2.0)*2.0*pi
        rr = (1.0/8.0)*2.0*pi
        r = sqrt( (xi-xx)**2.0 + (yi-yy)**2.0 )
        if(r-rr.lt.0.0) unk(1,i,j,k,lb) = -1.0+2.0*exp(-ep*ep/((r-rr)*(r-rr)))

        xx = (3.0/4.0)*2.0*pi
        yy = (3.0/4.0)*2.0*pi
        rr = (1.0/8.0)*2.0*pi
        r = sqrt( (xi-xx)**2.0 + (yi-yy)**2.0 )
        if(r-rr.lt.0.0) unk(1,i,j,k,lb) = -1.0+2.0*exp(-ep*ep/((r-rr)*(r-rr)))

        unk(4,i,j,k,lb) = unk(1,i,j,k,lb)
        unk(5,i,j,k,lb) = unk(1,i,j,k,lb)

        unk(6,i,j,k,lb) = 0.0
        unk(9,i,j,k,lb) = 0.0
        unk(10,i,j,k,lb) = 0.0

      enddo ! i
    enddo ! j
  enddo ! k

  mype = 0
  iopt = 1
  nlayers = nguard
  call amr_guardcell(mype,iopt,nlayers)

  return

end subroutine app_initial_soln_blk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_initial_soln_end
!This subroutine is used to set the inital soln for variable 2 and 3 in the 
!Cahn-Hilliard equations. -- Tim
!--------------------------------------------------------------
! include files 
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters

  use workspace
  use paramesh_interfaces
  use paramesh_mpi_interfaces

  implicit none
  include 'mpif.h'
  
  integer :: lb, l, sweep
  integer :: nguard0, i, j, k, mype, iopt, nlayers
  double precision :: xi, yi, pi, dx, pikk=1.0
  pi=4.0*datan(pikk)

  nguard0 = nguard*npgs
  do lb=1,lnblocks
    if(nodetype(lb).eq.1) then
      dx = bsize(1,lb)/real(nxb)
      do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
        do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
          do i=il_bnd+nguard0,iu_bnd-nguard0

          enddo ! i
        enddo ! j
      enddo ! k
    endif ! nodetype(lb).eq.1
  enddo ! end loop over grid blocks

  mype = 0
  iopt = 1
  nlayers = nguard
  call amr_guardcell(mype,iopt,nlayers)

  return
end subroutine app_initial_soln_end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_test_refinement
  return
end subroutine app_test_refinement
