
!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
!#include "paramesh_preprocessor.fh"

subroutine amr_initial_soln
  
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters
  implicit none
  include 'mpif.h'
  
  integer :: lb

!--------------------------------------------------------------

! loop over leaf grid blocks
  if(lnblocks.gt.0) then
     do lb=1,lnblocks
        
        if(nodetype(lb).eq.1 .or. advance_all_levels) then
           
              call amr_initial_soln_blk(lb)
           
        endif

     enddo ! end loop over grid blocks
  endif

  return
end subroutine amr_initial_soln


subroutine amr_initial_soln_blk(lb)
  
!--------------------------------------------------------------
! include files 
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use solution_parameters
  implicit none
  include 'mpif.h'
  
  integer, intent(in) :: lb
  integer :: nguard0, i, j, k
  double precision :: radius, pi, dx, dy, dz, xi, yi, zi
!--------------------------------------------------------------

  pi=4.0*datan(1.0)
  nguard0 = nguard*npgs

! set values for unk
  dx = bsize(1,lb)/real(nxb)
  do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
     zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5) - nucleate_z
     if (ndim.eq.2) zi = 0.0
     do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5) - nucleate_y
        do i=il_bnd+nguard0,iu_bnd-nguard0
           unk(1,i,j,k,lb) = 0.0
           xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5) - nucleate_x
           
           ! Phase
           unk(1,i,j,k,lb) = -tanh(0.05*(xi*xi+yi*yi+zi*zi-nuc_radius*nuc_radius))
           unk(4,i,j,k,lb) = -tanh(0.05*(xi*xi+yi*yi+zi*zi-nuc_radius*nuc_radius))
           unk(5,i,j,k,lb) = -tanh(0.05*(xi*xi+yi*yi+zi*zi-nuc_radius*nuc_radius))

           ! Temperature
           unk(1+nunkvbles,i,j,k,lb) = delta+abs(delta)*0.5*(-tanh(0.05*(xi*xi+yi*yi+zi*zi-nuc_radius*nuc_radius))+1.0)
           unk(4+nunkvbles,i,j,k,lb) = delta+abs(delta)*0.5*(-tanh(0.05*(xi*xi+yi*yi+zi*zi-nuc_radius*nuc_radius))+1.0)
           unk(5+nunkvbles,i,j,k,lb) = delta+abs(delta)*0.5*(-tanh(0.05*(xi*xi+yi*yi+zi*zi-nuc_radius*nuc_radius))+1.0)

           ! Solute concentration
           if(nvar.eq.3*nunkvbles) then
              ! Solute and thermal
              unk(1+2*nunkvbles,i,j,k,lb) = 0.0
              unk(4+2*nunkvbles,i,j,k,lb) = 0.0
              unk(5+2*nunkvbles,i,j,k,lb) = 0.0
           else if (solute) then
              ! Solute only - no thermal
              unk(1+nunkvbles,i,j,k,lb) = 0.0
              unk(4+nunkvbles,i,j,k,lb) = 0.0
              unk(5+nunkvbles,i,j,k,lb) = 0.0
           endif
        enddo
     enddo
  enddo

  return
end subroutine amr_initial_soln_blk
