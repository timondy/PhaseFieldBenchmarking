!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_mg_smooth                                                 REQUIRED
!!!!  * Main smoothing call
!!!!  * Loops over real variables calling pf_mg_smooth_var
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  pf_mg_smooth_var
!!!!  * Smoothing of each variable
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  
!!!!  app_mg_get_defect_var                                         REQUIRED
!!!!  * Calculates FAS MG defect, max defect and RMS residual
!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Latest version: Chris Goodyer, May 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_mg_smooth(pid,level)
  use paramesh_dimensions
  use tree
  use multigrid_parameters
  use physicaldata
  use paramesh_interfaces
  implicit none
  include 'mpif.h'
  integer, intent(in) :: pid,level
  integer lb, iopt, i, j, k, nguard0

  ! Cycle over all blocks
  do lb=block_starts(level), block_starts(level+1)-1
    call pf_mg_smooth_var(lb,pid,iopt,level,0)
  end do
  iopt=1
  call pf_guardcell(pid,iopt,level)
!!  do lb=block_starts(level), block_starts(level+1)-1
!!    call pf_mg_smooth_var(lb,pid,iopt,level,1)
!!  end do
!!  iopt=1
!!  call pf_guardcell(pid,iopt,level)
end subroutine app_mg_smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pf_mg_smooth_var(lb,pid,iopt,level,red_black)
  ! Smoother for node centered variables
  ! lb is the block number
  ! n is the variable number being smoothed
  use time_dep_parameters
  use solution_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use paramesh_interfaces

  implicit none
  integer,intent(in) :: lb, pid, iopt, level, red_black
  integer :: i,j,k,order
  double precision :: dx,Fi,dFi_dvi, xi, yi, pi, pikk=1.0
  double precision :: rfactor, rf1, rf2, rf
  double precision :: a11, a12, a21, a22, b1, b2, determinant
  !!!-----------------------------------------------------------------------------------------
  !!!-----------------------------------------------------------------------------------------
  rfactor = dt/dtold
  pi=4.0*datan(pikk)

  dx = bsize(1,lb)/real(nxb)

  if(red_black.eq.0)then
    order = 1
  else if(red_black.eq.1)then
    order = 0
  endif

  do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
    do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
      do i=il_bnd+nguard,iu_bnd-nguard
!!      do i=il_bnd+nguard+order,iu_bnd-nguard,2
        xi = bnd_box(1,1,lb) + dx*(real(i-nguard)-.5)
        yi = bnd_box(1,2,lb) + dx*(real(j-nguard)-.5)
        !!!-----------------------------------------------------------------------------------------
        if(mode_1.eq.1)then
          a11 = 1.0
          a12 = (4.0*dt)/(dx*dx)
          a21 = -3.0*unk(1,i,j,k,lb)*unk(1,i,j,k,lb)&
                + 1.0 -((4.0*ep*ep)/(dx*dx))
          a22 = 1.0

          b1 = unk(1,i,j,k,lb) - unk(2,i,j,k,lb)&
              -(dt/(dx*dx))&
              *(unk(6,i+1,j,k,lb)+unk(6,i-1,j,k,lb)&
               +unk(6,i,j+1,k,lb)+unk(6,i,j-1,k,lb)&
               -4.0*unk(6,i,j,k,lb))

          b2 = unk(6,i,j,k,lb)&
             - unk(1,i,j,k,lb)*unk(1,i,j,k,lb)&
              *unk(1,i,j,k,lb) + unk(1,i,j,k,lb)&
              +((ep*ep)/(dx*dx))&
              *(unk(1,i+1,j,k,lb)+unk(1,i-1,j,k,lb)&
              + unk(1,i,j+1,k,lb)+unk(1,i,j-1,k,lb)&
              -4.0*unk(1,i,j,k,lb))&
              -unk(7,i,j,k,lb)

          determinant = a11*a22-a12*a21

          unk(1,i,j,k,lb) = unk(1,i,j,k,lb)&
          - (b1*a22 - b2*a12)/determinant
          unk(6,i,j,k,lb) = unk(6,i,j,k,lb)&
          - (a11*b2 - a21*b1)/determinant
        !!!-----------------------------------------------------------------------------------------
        elseif(mode_1.eq.2)then
          a11 = 1.0
          a12 = (2.0/3.0)*(4.0*dt)/(dx*dx)
          a21 = -3.0*unk(1,i,j,k,lb)*unk(1,i,j,k,lb)&
                + 1.0 -((4.0*ep*ep)/(dx*dx))
          a22 = 1.0

          b1 = unk(1,i,j,k,lb) - unk(2,i,j,k,lb)&
              -(2.0/3.0)*(dt/(dx*dx))&
              *(unk(6,i+1,j,k,lb)+unk(6,i-1,j,k,lb)&
               +unk(6,i,j+1,k,lb)+unk(6,i,j-1,k,lb)&
               -4.0*unk(6,i,j,k,lb))

          b2 = unk(6,i,j,k,lb) - unk(1,i,j,k,lb)&
              *unk(1,i,j,k,lb)*unk(1,i,j,k,lb)&
              + unk(1,i,j,k,lb)&
              +((ep*ep)/(dx*dx))&
              *(unk(1,i+1,j,k,lb)+unk(1,i-1,j,k,lb)&
               +unk(1,i,j+1,k,lb)+unk(1,i,j-1,k,lb)&
               -4.0*unk(1,i,j,k,lb))&
              -unk(7,i,j,k,lb)

          determinant = a11*a22-a12*a21

          unk(1,i,j,k,lb) = unk(1,i,j,k,lb)&
          - (b1*a22 - b2*a12)/determinant
          unk(6,i,j,k,lb) = unk(6,i,j,k,lb)&
          - (a11*b2 - a21*b1)/determinant
        !!!-----------------------------------------------------------------------------------------
        endif
      end do
      if(order.eq.1)then
        order = 0
      else if(order.eq.0)then
        order = 1
      endif
    end do
  end do
end subroutine pf_mg_smooth_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pf_mg_smooth_var_none_time(lb)
! For non time-dependent variables. --Tim
  write(*,*) "------shouldn't be calling this subroutine 'pf_mg_smooth_var_none_time'"
end subroutine pf_mg_smooth_var_none_time


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine app_mg_get_defect_var(lb, max_defect, rms, npts)
  ! Defect for node centered variables
  ! lb is the block number
  ! n is the variable number whose defect is being calculated
  ! defect will be put into work(i,j,k,lb,n+1)
  use multigrid_parameters
  use paramesh_dimensions
  use time_dep_parameters
  use solution_parameters
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: lb
  integer, intent(inout) :: npts
  double precision, intent(inout) :: max_defect, rms
  integer :: i,j,k,level
  double precision :: dx,Fi
  double precision :: rfactor, rf1, rf2, rf3, xi, yi

  dx = bsize(1,lb)/real(nxb)

  if(current_var.eq.1)then
    do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
      do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
          xi = bnd_box(1,1,lb) + dx*(real(i-nguard)-.5)
          yi = bnd_box(1,2,lb) + dx*(real(j-nguard)-.5)
          !-------------------------------------------------------------
          if(mode_1.eq.1)then
            Fi = unk(1,i,j,k,lb) - unk(2,i,j,k,lb)&
                -(dt/(dx*dx))&
                *(unk(6,i+1,j,k,lb)+unk(6,i-1,j,k,lb)&
                 +unk(6,i,j+1,k,lb)+unk(6,i,j-1,k,lb)&
                 -4.0*unk(6,i,j,k,lb))
          elseif(mode_1.eq.2)then
            Fi = unk(1,i,j,k,lb) - unk(2,i,j,k,lb)&
                -(2.0/3.0)*(dt/(dx*dx))&
                *(unk(6,i+1,j,k,lb)+unk(6,i-1,j,k,lb)&
                 +unk(6,i,j+1,k,lb)+unk(6,i,j-1,k,lb)&
                 -4.0*unk(6,i,j,k,lb))
          endif
          !-------------------------------------------------------------

          ! Add new residual to work(2).  
          !   Will be zero at start of coarsening.  
          !   This then sets the value that gets coarsened and then updates with defect
          work(i,j,k,lb,current_work+1) = -work(i,j,k,lb,current_work+1)-Fi

          if(Fi.lt.0.0)then
            Fi=-Fi
          end if
          if(Fi.gt.max_defect)then
            max_defect=Fi
          end if
          rms = rms + Fi*Fi
          npts = npts + 1
        end do
      end do
    end do
  else if(current_var.eq.2)then
    do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 
      do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
        do i=il_bnd+nguard,iu_bnd-nguard
          !-------------------------------------------------------------
            Fi = unk(6,i,j,k,lb)&
            - unk(1,i,j,k,lb)&
            *unk(1,i,j,k,lb)*unk(1,i,j,k,lb)&
            + unk(1,i,j,k,lb)&
            +((ep*ep)/(dx*dx))&
            *(unk(1,i+1,j,k,lb)+unk(1,i-1,j,k,lb)&
            +unk(1,i,j+1,k,lb)+unk(1,i,j-1,k,lb)&
            -4.0*unk(1,i,j,k,lb))&
            -unk(7,i,j,k,lb)
          !-------------------------------------------------------------

          work(i,j,k,lb,current_work+1) = -work(i,j,k,lb,current_work+1)-Fi

          if(Fi.lt.0.0)then
            Fi=-Fi
          end if
          if(Fi.gt.max_defect)then
            max_defect=Fi
          end if
          rms = rms + Fi*Fi
          npts = npts + 1
        end do
      end do
    end do
  endif

  return
end subroutine app_mg_get_defect_var

