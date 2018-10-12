!#include "paramesh_preprocessor.fh"

!!!!!!!  Note that get_stencil_diff and get_stencil are now 3-d only
!!!!!!!    2-d versions follow

subroutine get_stencil_diff(i,j,k,dx,lb,m,dFi_dvi,LapMult)
  use multigrid_parameters
  use solution_parameters
  use time_dep_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  use multigrid_parameters
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,k,lb,m
  real, intent(in) :: dx,LapMult
  real, intent(out) :: dFi_dvi
  double precision :: U,phi,c, divisor, A_func,alewis, phix, phiy, phiz, gradphisq
  
  double precision :: rfactor, rf4, rf5, rf6

  if (mode_1.eq.2) then 
     rfactor = dt/dtold
     rf4 = rfactor*rfactor/(rfactor+1.0)/dt
     rf5 = (rfactor+1.0)/dt
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)/dt
  else
     rf4=0.0
     rf5=1.0/dt
     rf6=1.0/dt
  endif

  ! dFi_dvi=d(N(v))/dvi
  ! N(v) is "whatever equals f" for node i,j,k
  ! Ie solving N(v)=f. N() can be non-linear
  ! d(N(v))/dvi is the result of differentiating N(v) wrt the node
  !    currently being solved (not the same as d(N(v))/dv)
  ! Eg if N(v)=v(x+1)-2v(x)+v(x-1) (1D d2v/dx2)
  ! d(N(v))/dvi=-2
  ! n is the variable number
  ! m is the variable type (cell/face/edge/vertex centered)
  ! m=1 -> cell centered
  !   2 -> vertex
  !   3 -> face
  !   4 -> edge
  ! so get_dFi_dvi(i,j,k,dx,lb,2,1,dFi_dvi)
  !    corresponds to unk(2,i,j,k,lb)
  ! USER EDIT CODE
  U=unk(1+nunkvbles,i,j,k,lb)
  phi=unk(1,i,j,k,lb)
  if(total_vars.eq.3)then
     c=unk(1+2*nunkvbles,i,j,k,lb)
  ! Handler code for when running phase/thermal or phase/solute only (ie not fully coupled)
  else
     if(solute)then
        U=delta
        c=unk(1+nunkvbles,i,j,k,lb)
        if(current_var.eq.2)then
           current_var=3
        end if
     else
        c=0.0
     end if
  end if
  if(m.eq.1)then
     if(current_var.eq.1)then
        if((total_vars.lt.3))then
           ! For thermal only model mcinf=0, alewis=1;
           ! For solute only model alewis=0 and c is actually U_0 which boils term down to 1
           divisor=1.0
        else
           divisor=alewis+mcinf*(1.0+(1.0-ke)*c)
        end if
        call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
        call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
        call get_dunk_n_dz(1,i,j,k,dx,lb,phiz,.false.)
        gradphisq = (phix*phix+phiy*phiy+phiz*phiz)
        call get_w_aniso_2(phix,phiy,phiz,gradphisq,A_func)
! 9pt
!        dFi_dvi=A_0*A_0/divisor*(-2.0/(dx*dx)-2.0/(dy*dy)-2.0/(dz*dz))-1.0/divisor*(3.0*phi*phi-1.0+lambda*U*4.0*phi*(phi*phi-1))
! 19 pt
        dFi_dvi=-LapMult*4.0/(dx*dx*divisor)-(3.0*phi*phi-1.0+lambda*(U+mcinf*c)*4.0*phi*(phi*phi-1))/(A_func*A_func*divisor)
     else if(current_var.eq.2)then
        dFi_dvi=D_therm*(-4.0/(dx*dx))
     else if(current_var.eq.3)then
        if((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb).eq.0)then 
           dFi_dvi=0
        else if (mode_1.eq.1) then
           dFi_dvi=D_solute*(1-unk(1,i,j,k,lb))/2.0*(-2.0/(dx*dx)-2.0/(dx*dx)-2.0/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
                0.5*(1-ke)*(phi-unk(2,i,j,k,lb))/dt/((1+ke)/2-(1-ke)/2*phi)
        else           
           dFi_dvi=D_solute*(1-unk(1,i,j,k,lb))/2.0*(-2.0/(dx*dx)-2.0/(dx*dx)-2.0/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
                0.5*(1-ke)*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/((1+ke)/2-(1-ke)/2*phi)
        end if
     else
        print *,"Unhandled m/cvar combination in subroutine get_dFi_dvi m=",m,"n=",current_var
     end if
  else if(m.eq.2)then
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  else if(m.eq.3)then
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  else if(m.eq.4)then
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  else
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  end if
  if((total_vars.lt.3).and.(solute).and.(current_var.eq.3))then
     current_var=2
  end if
  ! END USER EDIT
end subroutine get_stencil_diff

subroutine get_stencil(i,j,k,dx,lb,m,Fi,LapMult)
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,k,lb,m
  real, intent(in) :: dx
  real, intent(out) :: Fi,LapMult
!  real :: sin_theta,cos_theta,sin_phi,cos_phi,phix,phiy,phiz
  real :: phix,phiy,phiz
!  real :: thetax,thetay,thetaz,phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz
  real :: phixx,phixy,phixz,phiyy,phiyz,phizz
!  real :: A_func,dw_dw_dphi,t1,t2,t3,t4,t5,t6
  real :: A_func,t1,t2,t3, laplacian
  real :: dw_dphix,dw_dphiy,dw_dphiz,dw_dx,dw_dy,dw_dz,gradphisq,cx,cy,cz,cxx,cyy,czz,phi_dot_x,phi_dot_y,phi_dot_z
  double precision :: dh, phi, U, c, alewis, divisor
  double precision :: rfactor, rf4, rf5, rf6
  double precision :: A_func3d,Ai_func3d,Aij_func3d,Aii_func3d,xx,yy,zz
  double precision :: Ai(3),Aij(3,3),Gij(3,3),trace_Gij
  if (mode_1.eq.2) then 
     rfactor = dt/dtold
     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
  else
     rf4=0.0
     rf5=1.0
     rf6=1.0
  endif
  
 

  ! Fi=N(v)
  ! f is the right hand side for node i,j,k
  ! N(v) is "whatever equals f" for node i,j,k
  ! Ie solving N(v)=f. N() can be non-linear
  ! n is the variable number
  ! m is the variable type (cell/face/edge/vertex centered)
  ! m=1 -> cell centered
  !   2 -> vertex
  !   3 -> face
  !   4 -> edge
  ! so get_Fi(i,j,k,dx,lb,2,1,Fi)
  !    corresponds to unk(2,i,j,k,lb)
  ! Ignore above, had such trouble getting mg working for non cell centered node I've not
  ! coded any of the others.
  ! USER EDIT CODE
  U=unk(1+nunkvbles,i,j,k,lb)
  phi=unk(1,i,j,k,lb)
  if(total_vars.eq.3)then
     c=unk(1+2*nunkvbles,i,j,k,lb)
     alewis = 1.0/le    
  else
  ! Handler code for when running phase/thermal or phase/solute only (ie not fully coupled)
     if(solute)then
        alewis = 0.0
        U = delta
        if(current_var.eq.2)then
           current_var=3
        end if
        c=unk(1+nunkvbles,i,j,k,lb)
     else
        alewis = 1.0/le    
        c=0.0
     end if
  end if

  if(m.eq.1)then
     if(current_var.eq.1)then
        call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
        call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
        call get_dunk_n_dz(1,i,j,k,dx,lb,phiz,.false.)
        gradphisq=(phix*phix+phiy*phiy+phiz*phiz)
        call get_w_aniso_2(phix,phiy,phiz,gradphisq,A_func)
        if((total_vars.lt.3))then
           ! For thermal only model mcinf=0, alewis=1;
           ! For solute only model alewis=0 and c is actually U_0 which boils term down to 1
!           divisor=A_func*A_func
           divisor=1.0
        else
!           divisor=A_func*A_func*(alewis+mcinf*(1.0+(1.0-ke)*c))
           divisor=alewis+mcinf*(1.0+(1.0-ke)*c)
        end if
         
!       9 pt stencil
!        laplacian=1.0*((unk(1,i+1,j,k,lb)-2.0*unk(1,i,j,k,lb)+unk(1,i-1,j,k,lb))/(dx*dx)&
!             +(unk(1,i,j+1,k,lb)-2.0*unk(1,i,j,k,lb)+unk(1,i,j-1,k,lb))/(dy*dy)&
!             +(unk(1,i,j,k+1,lb)-2.0*unk(1,i,j,k,lb)+unk(1,i,j,k-1,lb))/(dz*dz))
        ! New 19-point del-squared stencil
!        laplacian=A_func*A_func/(6.0*dx*dx)*(&
        laplacian=1.0/(6.0*dx*dx)*(&
             1.0*(unk(1,i+1,j,k-1,lb)+unk(1,i-1,j,k-1,lb)+unk(1,i,j+1,k-1,lb)+unk(1,i,j-1,k-1,lb)+&
             unk(1,i+1,j+1,k,lb)+unk(1,i+1,j-1,k,lb)+unk(1,i-1,j+1,k,lb)+unk(1,i-1,j-1,k,lb)+&
             unk(1,i+1,j,k+1,lb)+unk(1,i-1,j,k+1,lb)+unk(1,i,j+1,k+1,lb)+unk(1,i,j-1,k+1,lb))+&
             2.0*(unk(1,i,j,k-1,lb)+&!Should be 2*
             unk(1,i+1,j,k,lb)+unk(1,i-1,j,k,lb)+unk(1,i,j+1,k,lb)+unk(1,i,j-1,k,lb)+&
             unk(1,i,j,k+1,lb))&
             -24.0*unk(1,i,j,k,lb))!should be 24*
!        27-point stencil - from page 64 of Jan's thesis  
!             CEG - still a bug as of 24/5/10, probably in amr_1blk_bcset.f90
!        laplacian=1.0/(30.0*dx*dx)*(&
!             1.0*(unk(1,i+1,j+1,k+1,lb)+unk(1,i+1,j+1,k-1,lb)+&
!             unk(1,i+1,j-1,k+1,lb)+unk(1,i+1,j-1,k+1,lb)+&
!             unk(1,i-1,j+1,k+1,lb)+unk(1,i-1,j+1,k-1,lb)+&
!             unk(1,i-1,j-1,k+1,lb)+unk(1,i-1,j-1,k-1,lb))+&
!             3.0*(unk(1,i+1,j,k-1,lb)+unk(1,i-1,j,k-1,lb)+unk(1,i,j+1,k-1,lb)+unk(1,i,j-1,k-1,lb)+&
!             unk(1,i+1,j+1,k,lb)+unk(1,i+1,j-1,k,lb)+unk(1,i-1,j+1,k,lb)+unk(1,i-1,j-1,k,lb)+&
!             unk(1,i+1,j,k+1,lb)+unk(1,i-1,j,k+1,lb)+unk(1,i,j+1,k+1,lb)+unk(1,i,j-1,k+1,lb))+&
!             14.0*(unk(1,i,j,k-1,lb)+&!Should be 14*
!             unk(1,i+1,j,k,lb)+unk(1,i-1,j,k,lb)+unk(1,i,j+1,k,lb)+unk(1,i,j-1,k,lb)+&
!             unk(1,i,j,k+1,lb))&
!             -128.0*unk(1,i,j,k,lb))!should be 128*

        Fi=-(phi*(phi*phi-1.0)+lambda*(U+mcinf*c)*(1-phi*phi)**2) 
 
        LapMult=1.0
        if(gradphisq.gt.0.0.and.epsilon_tilde.gt.0.0) then
           call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
           call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
           call get_d2unk_n_dxdz(1,i,j,k,dx,lb,phixz,.false.)
           call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
           call get_d2unk_n_dydz(1,i,j,k,dx,lb,phiyz,.false.)
           call get_d2unk_n_dz2(1,i,j,k,dx,lb,phizz,.false.)
           
           xx=phix*phix/gradphisq
           yy=phiy*phiy/gradphisq
           zz=phiz*phiz/gradphisq
           
           Ai(1)=Ai_func3d(phix,xx,yy,zz)
           Ai(2)=Ai_func3d(phiy,yy,zz,xx)
           Ai(3)=Ai_func3d(phiz,zz,xx,yy)
           
           Aij(1,1)=Aii_func3d(xx,yy,zz)
           Aij(2,2)=Aii_func3d(yy,zz,xx)
           Aij(3,3)=Aii_func3d(zz,xx,yy)
           
           Aij(1,2)=Aij_func3d(phix,phiy,xx,yy,zz)/gradphisq
           Aij(2,3)=Aij_func3d(phiy,phiz,yy,zz,xx)/gradphisq
           Aij(3,1)=Aij_func3d(phiz,phix,zz,xx,yy)/gradphisq
           
           Gij(1,1)=(Ai(1)*Ai(1)+2.0*A_func*(Ai(1)*phix+Ai(1)*phix))/gradphisq+A_func*Aij(1,1)
           Gij(2,2)=(Ai(2)*Ai(2)+2.0*A_func*(Ai(2)*phiy+Ai(2)*phiy))/gradphisq+A_func*Aij(2,2)
           Gij(3,3)=(Ai(3)*Ai(3)+2.0*A_func*(Ai(3)*phiz+Ai(3)*phiz))/gradphisq+A_func*Aij(3,3)
           Gij(1,2)=(Ai(1)*Ai(2)+2.0*A_func*(Ai(1)*phiy+Ai(2)*phix))/gradphisq+A_func*Aij(1,2)
           Gij(2,3)=(Ai(2)*Ai(3)+2.0*A_func*(Ai(2)*phiz+Ai(3)*phiy))/gradphisq+A_func*Aij(2,3)
           Gij(3,1)=(Ai(3)*Ai(1)+2.0*A_func*(Ai(3)*phix+Ai(1)*phiz))/gradphisq+A_func*Aij(3,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reproduce old version by removing these 6 lines        !
           trace_Gij=(Gij(1,1)+Gij(2,2)+Gij(3,3))/3.0    !
                                                         !
           Gij(1,1)=Gij(1,1)-trace_Gij                   !
           Gij(2,2)=Gij(2,2)-trace_Gij                   !
           Gij(3,3)=Gij(3,3)-trace_Gij                   !
                                                         !
           LapMult = LapMult+trace_Gij/(A_func*A_func)   !
                                                         !
           Fi = Fi+trace_Gij*laplacian                   !
                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
           Fi = Fi +2.0*(Gij(1,2)*phixy+Gij(2,3)*phiyz+Gij(3,1)*phixz)
           Fi = Fi +Gij(1,1)*phixx+Gij(2,2)*phiyy+Gij(3,3)*phizz
           
        end if
        

        Fi=Fi/(divisor*A_func*A_func)
        Fi = Fi + laplacian/divisor
     else if(current_var.eq.2)then
        call get_dh_dpsi(phi,dh)
        if (mode_1.eq.1) then
!        if (zmode_1.ge.1) then
           Fi=D_therm*(((unk(1+nunkvbles,i+1,j,k,lb)-2.0*unk(1+nunkvbles,i,j,k,lb)+unk(1+nunkvbles,i-1,j,k,lb))/(dx*dx)&
             +(unk(1+nunkvbles,i,j+1,k,lb)-2.0*unk(1+nunkvbles,i,j,k,lb)+unk(1+nunkvbles,i,j-1,k,lb))/(dx*dx)&
             +(unk(1+nunkvbles,i,j,k+1,lb)-2.0*unk(1+nunkvbles,i,j,k,lb)+unk(1+nunkvbles,i,j,k-1,lb))/(dx*dx)  &
             ))  +0.5*dh*(phi-unk(2,i,j,k,lb))/dt
        else
           ! "phi" is new phi, i.e. unk1
!9 point
!           Fi=D_therm*(((unk(1+nunkvbles,i+1,j,k,lb)-2.0*unk(1+nunkvbles,i,j,k,lb)+unk(1+nunkvbles,i-1,j,k,lb))/(dx*dx)&
!             +(unk(1+nunkvbles,i,j+1,k,lb)-2.0*unk(1+nunkvbles,i,j,k,lb)+unk(1+nunkvbles,i,j-1,k,lb))/(dx*dx)&
!             +(unk(1+nunkvbles,i,j,k+1,lb)-2.0*unk(1+nunkvbles,i,j,k,lb)+unk(1+nunkvbles,i,j,k-1,lb))/(dx*dx)& 
!             )) +0.5*dh*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
!19 point
        Fi=D_therm*(&
              1.0*(unk(1+nunkvbles,i+1,j,k-1,lb)+unk(1+nunkvbles,i-1,j,k-1,lb)+  &
                    unk(1+nunkvbles,i,j+1,k-1,lb)+unk(1+nunkvbles,i,j-1,k-1,lb)+ &
                  unk(1+nunkvbles,i+1,j+1,k,lb)+unk(1+nunkvbles,i+1,j-1,k,lb)+   &
                    unk(1+nunkvbles,i-1,j+1,k,lb)+unk(1+nunkvbles,i-1,j-1,k,lb)+ &
                  unk(1+nunkvbles,i+1,j,k+1,lb)+unk(1+nunkvbles,i-1,j,k+1,lb)+   &
                    unk(1+nunkvbles,i,j+1,k+1,lb)+unk(1+nunkvbles,i,j-1,k+1,lb))+&
              2.0*(unk(1+nunkvbles,i,j,k-1,lb)+unk(1+nunkvbles,i+1,j,k,lb)+      &
                    unk(1+nunkvbles,i-1,j,k,lb)+unk(1+nunkvbles,i,j+1,k,lb)+     &
                    unk(1+nunkvbles,i,j-1,k,lb)+unk(1+nunkvbles,i,j,k+1,lb))     &
              -24.0*unk(1+nunkvbles,i,j,k,lb)                                    &
              )/(6.0*dx*dx)                                                           &              ! Note assuming dx=dy=dz
             +0.5*dh*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
        endif
     else if(current_var.eq.3)then
        call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
        call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
        call get_dunk_n_dz(1,i,j,k,dx,lb,phiz,.false.)
        call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
        call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
        call get_d2unk_n_dz2(1,i,j,k,dx,lb,phizz,.false.)
        call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
        call get_d2unk_n_dxdz(1,i,j,k,dx,lb,phixz,.false.)
        call get_d2unk_n_dydz(1,i,j,k,dx,lb,phiyz,.false.)
        call get_dunk_n_dx(current_unk,i,j,k,dx,lb,cx,.false.)
        call get_dunk_n_dy(current_unk,i,j,k,dx,lb,cy,.false.)
        call get_dunk_n_dz(current_unk,i,j,k,dx,lb,cz,.false.)
        call get_d2unk_n_dx2(current_unk,i,j,k,dx,lb,cxx,.false.)
        call get_d2unk_n_dy2(current_unk,i,j,k,dx,lb,cyy,.false.)
        call get_d2unk_n_dz2(current_unk,i,j,k,dx,lb,czz,.false.)
! INCORRECTLY USING BDF1 here to make 3rd antitrapping work
!        if (mode_1.eq.1) then
        if (mode_1.ge.1) then
           phi_dot_x=1.0/(2.0*dt*dx)*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)-(unk(4,i+1,j,k,lb)-unk(4,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)-(unk(4,i,j+1,k,lb)-unk(4,i,j-1,k,lb)))
           phi_dot_z=1.0/(2.0*dt*dx)*(unk(1,i,j,k+1,lb)-unk(1,i,j,k-1,lb)-(unk(4,i,j,k+1,lb)-unk(4,i,j,k-1,lb)))
        else 
           phi_dot_x=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)) &
                                      - rf5*(unk(4,i+1,j,k,lb)+unk(4,i-1,j,k,lb))&
                                      + rf4*(unk(5,i+1,j,k,lb)-unk(5,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)) &
                                      - rf5*(unk(4,i,j+1,k,lb)+unk(4,i,j-1,k,lb))&
                                      + rf4*(unk(5,i,j+1,k,lb)-unk(5,i,j-1,k,lb)))
           phi_dot_z=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i,j,k+1,lb)-unk(1,i,j,k-1,lb)) &
                                      - rf5*(unk(4,i,j,k+1,lb)+unk(4,i,j,k-1,lb))&
                                      + rf4*(unk(5,i,j,k+1,lb)-unk(5,i,j,k-1,lb)))
        endif
        Fi=0
! 1st antitrapping term
        t1=phix*cx+phiy*cy+phiz*cz
        if(phix*phix+phiy*phiy+phiz*phiz.eq.0)then
           t2=0
        else
           if (mode_1.eq.1) then
!           if (mode_1.ge.1) then
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy+phiz*phiz)
           else
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (rf6*unk(1,i,j,k,lb)- rf5*unk(4,i,j,k,lb)+ &
                    rf4*unk(5,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy+phiz*phiz)
           endif
        end if
        Fi=t1*t2

! 2nd antitrapping term
        t1=cxx+cyy+czz
        t2=D_solute/2*(1-unk(1,i,j,k,lb))
        Fi=Fi+t1*t2
        if(phix*phix+phiy*phiy+phiz*phiz.eq.0)then
           t1=0
        else
           t1=(phixx+phiyy+phizz)/sqrt(phix*phix+phiy*phiy+phiz*phiz)
           t1=t1-(phix*(phix*phixx+phiy*phixy+phiz*phixz)&
                +phiy*(phix*phixy+phiy*phiyy+phiz*phiyz)&
                +phiz*(phix*phixz+phiy*phiyz+phiz*phizz))/sqrt(phix*phix+phiy*phiy+phiz*phiz)**3
        end if
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/(dt*2.0*sqrt(2.0))
        else
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(rf6*unk(1,i,j,k,lb) -        &
                    rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb))/(dt*2.0*sqrt(2.0))
!           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/(dt*2.0*sqrt(2.0))
!           if (rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb)- (unk(1,i,j,k,lb)-unk(2,i,j,k,lb)).gt.1.0e-12)&
!           print *, mype, rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb), unk(1,i,j,k,lb)-unk(2,i,j,k,lb),&
!                 rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb)- (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))
        endif
        Fi=Fi+anti_trapping_mod*t1*t2


! 3rd antitrapping term
        t1=phix*phi_dot_x+phiy*phi_dot_y+phiz*phi_dot_z
        if(phix*phix+phiy*phiy+phiz*phiz.eq.0)then
           t2=0
        else
           t2=1.0/(2.0*sqrt(2.0))*(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))/sqrt(phix*phix+phiy*phiy+phiz*phiz)
        end if
        Fi=Fi+anti_trapping_mod*t1*t2
        t3=(1+ke)/2-(1-ke)/2*phi
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           Fi=Fi+0.5*(1+(1-ke)*c)*(phi-unk(2,i,j,k,lb))/dt
        else
           Fi=Fi+0.5*(1+(1-ke)*c)*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
        endif
        if(t3.eq.0.0)then
           Fi=0.0
        else
           Fi=Fi/t3
        end if
     else
        print *,"Unhandled m/c_var combination in subroutine get_Fi m=",m,"cvar=",current_var
     end if
  else if(m.eq.2)then
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  else if(m.eq.3)then
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  else if(m.eq.4)then
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  else
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  end if
  if((total_vars.lt.3).and.(solute).and.(current_var.eq.3))then
     current_var=2
  end if
  ! END USER EDIT
end subroutine get_stencil

!!!!!!!!!!!   Let's put the 2-d ones in here

subroutine get_stencil_diff_2d(i,j,dx,lb,m,dFi_dvi,LapMult,tau_val)
  use multigrid_parameters
  use solution_parameters
  use time_dep_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,lb,m
  real, intent(in) :: dx,LapMult,tau_val
  real, intent(out) :: dFi_dvi
  double precision :: U,phi,c, divisor,  alewis, phix, phiy, gradphisq
  double precision :: rfactor, rf4, rf5, rf6,Lap2d
  integer k
  k=1

  if (mode_1.eq.2) then 
     rfactor = dt/dtold
     rf4 = rfactor*rfactor/(rfactor+1.0)/dt
     rf5 = (rfactor+1.0)/dt
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)/dt
  else
     rf4=0.0
     rf5=1.0/dt
     rf6=1.0/dt
  endif

  ! dFi_dvi=d(N(v))/dvi
  ! N(v) is "whatever equals f" for node i,j,k
  ! Ie solving N(v)=f. N() can be non-linear
  ! d(N(v))/dvi is the result of differentiating N(v) wrt the node
  !    currently being solved (not the same as d(N(v))/dv)
  ! Eg if N(v)=v(x+1)-2v(x)+v(x-1) (1D d2v/dx2)
  ! d(N(v))/dvi=-2
  ! n is the variable number
  ! m is the variable type (cell/face/edge/vertex centered)
  ! m=1 -> cell centered
  !   2 -> vertex
  !   3 -> face
  !   4 -> edge
  ! so get_dFi_dvi(i,j,k,dx,dy,lb,2,1,dFi_dvi)
  !    corresponds to unk(2,i,j,k,lb)
  ! USER EDIT CODE
  U=unk(1+nunkvbles,i,j,k,lb)
  phi=unk(1,i,j,k,lb)
  if(total_vars.eq.3)then
     c=unk(1+2*nunkvbles,i,j,k,lb)
     alewis = 1.0/le    
  ! Handler code for when running phase/thermal or phase/solute only (ie not fully coupled)
  else
     if(solute)then
        U=delta
        c=unk(1+nunkvbles,i,j,k,lb)
        if(current_var.eq.2)then
           current_var=3
        end if
     else
        c=0.0
     end if
  end if
  if(m.eq.1)then
     if(current_var.eq.1)then
       dFi_dvi=LapMult*Lap2d(dx,.true.)-(3.0*phi*phi-1.0+lambda*(U+mcinf*c)*4.0*phi*(phi*phi-1))/tau_val
     else if(current_var.eq.2)then
        dFi_dvi=D_therm*Lap2d(dx,.true.)
     else if(current_var.eq.3)then
        if((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb).eq.0)then 
           dFi_dvi=0
        else if (mode_1.eq.1) then
           dFi_dvi=D_solute*(1-unk(1,i,j,k,lb))/2.0*(-2.0/(dx*dx)-2.0/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
                0.5*(1-ke)*(phi-unk(2,i,j,k,lb))/dt/((1+ke)/2-(1-ke)/2*phi)
        else           
           dFi_dvi=D_solute*(1-unk(1,i,j,k,lb))/2.0*(-2.0/(dx*dx)-2.0/(dx*dx))/((1+ke)/2-(1-ke)/2*unk(1,i,j,k,lb))+&
                0.5*(1-ke)*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/((1+ke)/2-(1-ke)/2*phi)
        end if
     else
        print *,"Unhandled m/cvar combination in subroutine get_dFi_dvi m=",m,"n=",current_var
     end if

  else if(m.eq.2)then
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  else if(m.eq.3)then
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  else if(m.eq.4)then
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  else
     print *,"Unhandled m/n combination in subroutine get_dFi_dvi m=",m,"n=",current_var
  end if
  if((total_vars.lt.3).and.(solute).and.(current_var.eq.3))then
     current_var=2
  end if
  ! END USER EDIT
end subroutine get_stencil_diff_2d

subroutine get_stencil_2d(i,j,dx,lb,m,Fi, dx2inv,LapMult,tau_val)
  use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  integer, intent(in) :: i,j,lb,m
  real, intent(in) :: dx, dx2inv
  real, intent(out) :: Fi,LapMult,tau_val
  real :: phix,phiy
  real :: phixx,phixy,phiyy
  real :: A_func,t1,t2,t3
  real :: dw_dphix,dw_dphiy,dw_dx,dw_dy,gradphisq,cx,cy,cxx,cyy,phi_dot_x,phi_dot_y
  double precision :: dh,phi,U,c, divisor, laplacian
  double precision :: rfactor, rf4, rf5, rf6
  double precision :: alewis
  double precision :: A_func2d
  double precision :: M11,M22,M12,A1,A2,A11,A22,A12, xx,yy,Lap2d,Lap_val
  integer k
  k=1

  if (mode_1.eq.2) then 
     rfactor = dt/dtold
     rf4 = rfactor*rfactor/(rfactor+1.0)
     rf5 = (rfactor+1.0)
     rf6 = (2.0*rfactor+1.0)/(rfactor+1.0)
  else
     rf4=0.0
     rf5=1.0
     rf6=1.0
  endif

  ! Fi=N(v)
  ! f is the right hand side for node i,j,k
  ! N(v) is "whatever equals f" for node i,j,k
  ! Ie solving N(v)=f. N() can be non-linear
  ! n is the variable number
  ! m is the variable type (cell/face/edge/vertex centered)
  ! m=1 -> cell centered
  !   2 -> vertex
  !   3 -> face
  !   4 -> edge
  ! so get_Fi(i,j,k,dx,dy,lb,2,1,Fi)
  !    corresponds to unk(2,i,j,k,lb)
  ! Ignore above, had such trouble getting mg working for non cell centered node I've not
  ! coded any of the others.
  ! USER EDIT CODE
  U=unk(1+nunkvbles,i,j,k,lb)
  phi=unk(1,i,j,k,lb)
  if(total_vars.eq.3)then
     c=unk(1+2*nunkvbles,i,j,k,lb)
     alewis = 1.0/le    
  else
  ! Handler code for when running phase/thermal or phase/solute only (ie not fully coupled)
     if(solute)then
        alewis = 0.0
        U = delta
        if(current_var.eq.2)then
           current_var=3
        end if
        c=unk(1+nunkvbles,i,j,k,lb)
     else
        alewis = 1.0/le    
        c=0.0
     end if
  end if
  if(m.eq.1)then
    if(current_var.eq.1)then
      call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
      call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
      gradphisq=phix*phix+phiy*phiy       
      call get_w_aniso_2_2d(phix,phiy,gradphisq,A_func)      
      if((total_vars.lt.3))then
       divisor=1.0
      else
       divisor=(alewis+mcinf*(1.0+(1.0-ke)*c))
      end if
      tau_val=A_func*A_func*divisor
      LapMult=1.0/divisor
      Lap_val=Lap2d(dx,.false.,i,j,k,lb,0)
      Fi=-(phi*(phi*phi-1.0)+lambda*(U+mcinf*c)*(1.0-phi*phi)**2)
! see Docs/AnisotropyPCB.tex for explanation of following      
      if(gradphisq.gt.0.0.and.epsilon_tilde.gt.0.0)then 
	call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
	call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
	call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
	xx=phix*phix/gradphisq
	yy=phiy*phiy/gradphisq
	A1=4.*A_0*epsilon_tilde*phix*yy*(xx-yy)       
	A2=4.*A_0*epsilon_tilde*phiy*xx*(yy-xx)
	A11=-4.*A_0*epsilon_tilde*yy*(3.*xx*xx-8.*xx*yy+yy*yy)
	A22=-4.*A_0*epsilon_tilde*xx*(xx*xx-8.*xx*yy+3.*yy*yy)
	A12=8.*A_0*epsilon_tilde*phix*phiy*(xx*xx-4.*xx*yy+yy*yy)/gradphisq
	M11=(A1*A1+4.*A_func*A1*phix)/gradphisq+A_func*A11
	M22=(A2*A2+4.*A_func*A2*phiy)/gradphisq+A_func*A22
	M12=(A1*A2+2.*A_func*(A1*phiy+A2*phix))/gradphisq+A_func*A12
 	LapMult=LapMult+(M11+M22)*0.5/tau_val	
        Fi=Fi+(M11+M22)*0.5*Lap_val
	Fi=Fi+(M11-M22)*0.5*(phixx-phiyy)+2.*M12*phixy
      endif
      Fi=Fi/tau_val+Lap_val/divisor
    else if(current_var.eq.2)then
      call get_dh_dpsi(phi,dh) 
      if (mode_1.eq.1) then
       Fi=D_therm*Lap2d(dx,.false.,i,j,k,lb,nunkvbles)  +0.5*dh*(phi-unk(2,i,j,k,lb))/dt
      else
       Fi=D_therm*Lap2d(dx,.false.,i,j,k,lb,nunkvbles) +0.5*dh*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
      endif
    else if(current_var.eq.3)then
        call get_dunk_n_dx(1,i,j,k,dx,lb,phix,.false.)
        call get_dunk_n_dy(1,i,j,k,dx,lb,phiy,.false.)
        call get_d2unk_n_dx2(1,i,j,k,dx,lb,phixx,.false.)
        call get_d2unk_n_dy2(1,i,j,k,dx,lb,phiyy,.false.)
        call get_d2unk_n_dxdy(1,i,j,k,dx,lb,phixy,.false.)
        call get_dunk_n_dx(current_unk,i,j,k,dx,lb,cx,.false.)
        call get_dunk_n_dy(current_unk,i,j,k,dx,lb,cy,.false.)
        call get_d2unk_n_dx2(current_unk,i,j,k,dx,lb,cxx,.false.)
        call get_d2unk_n_dy2(current_unk,i,j,k,dx,lb,cyy,.false.)
        if (mode_1.ge.1) then
           phi_dot_x=1.0/(2.0*dt*dx)*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)-(unk(4,i+1,j,k,lb)-unk(4,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)-(unk(4,i,j+1,k,lb)-unk(4,i,j-1,k,lb)))
        else 
           phi_dot_x=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i+1,j,k,lb)-unk(1,i-1,j,k,lb)) &
                                      - rf5*(unk(4,i+1,j,k,lb)+unk(4,i-1,j,k,lb))&
                                      + rf4*(unk(5,i+1,j,k,lb)-unk(5,i-1,j,k,lb)))
           phi_dot_y=1.0/(2.0*dt*dx)*(  rf6*(unk(1,i,j+1,k,lb)-unk(1,i,j-1,k,lb)) &
                                      - rf5*(unk(4,i,j+1,k,lb)+unk(4,i,j-1,k,lb))&
                                      + rf4*(unk(5,i,j+1,k,lb)-unk(5,i,j-1,k,lb)))
        endif
        Fi=0
! 1st antitrapping term
        t1=phix*cx+phiy*cy
        if(phix*phix+phiy*phiy.lt.1.0e-15)then
           t2=0
        else
           if (mode_1.eq.1) then
!           if (mode_1.ge.1) then
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy)
           else
              t2=-D_solute/2.0+1.0*anti_trapping_mod/(2.0*sqrt(2.0))*(1-ke)*&
                   (rf6*unk(1,i,j,k,lb)- rf5*unk(4,i,j,k,lb)+ &
                    rf4*unk(5,i,j,k,lb))/dt/sqrt(phix*phix+phiy*phiy)
           endif
        end if
        Fi=t1*t2

! 2nd antitrapping term
        t1=cxx+cyy
        if (ABS(t1).lt.1.0e-15) t1=0.0
        t2=D_solute/2*(1-unk(1,i,j,k,lb))
        Fi=Fi+t1*t2
!        print *,i,j,lb,Fi

        if(phix*phix+phiy*phiy.lt.1.0e-15)then
           t1=0
        else
           t1=(phixx+phiyy)/sqrt(phix*phix+phiy*phiy)
           t1=t1-(phix*(phix*phixx+phiy*phixy)&
                +phiy*(phix*phixy+phiy*phiyy)&
                )/sqrt(phix*phix+phiy*phiy)**3
        end if
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/(dt*2.0*sqrt(2.0))
        else
           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(rf6*unk(1,i,j,k,lb) -        &
                    rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb))/(dt*2.0*sqrt(2.0))
!           t2=(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))*(unk(1,i,j,k,lb)-unk(2,i,j,k,lb))/(dt*2.0*sqrt(2.0))
!           if (rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb)- (unk(1,i,j,k,lb)-unk(2,i,j,k,lb)).gt.1.0e-12)&
!           print *, mype, rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb), unk(1,i,j,k,lb)-unk(2,i,j,k,lb),&
!                 rf6*unk(1,i,j,k,lb) -rf5*unk(4,i,j,k,lb)  +  rf4*unk(5,i,j,k,lb)- (unk(1,i,j,k,lb)-unk(2,i,j,k,lb))
        endif
        Fi=Fi+anti_trapping_mod*t1*t2


! 3rd antitrapping term
        t1=phix*phi_dot_x+phiy*phi_dot_y
        if(phix*phix+phiy*phiy.lt.1.0e-15)then
           t2=0
        else
           t2=1.0/(2.0*sqrt(2.0))*(1.0+(1.0-ke)*unk(current_unk,i,j,k,lb))/sqrt(phix*phix+phiy*phiy)
        end if
        Fi=Fi+anti_trapping_mod*t1*t2
        t3=(1+ke)/2-(1-ke)/2*phi
        if (mode_1.eq.1) then
!        if (mode_1.ge.1) then
           Fi=Fi+0.5*(1+(1-ke)*c)*(phi-unk(2,i,j,k,lb))/dt
        else
           Fi=Fi+0.5*(1+(1-ke)*c)*(rf4*unk(5,i,j,k,lb)-rf5*unk(4,i,j,k,lb)+rf6*phi)/dt
        endif
        if(t3.eq.0.0)then
           Fi=0.0
        else
           Fi=Fi/t3
        end if
     else
        print *,"Unhandled m/c_var combination in subroutine get_Fi m=",m,"cvar=",current_var
     end if   
   
   

  else if(m.eq.2)then
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  else if(m.eq.3)then
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  else if(m.eq.4)then
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  else
     print *,"Unhandled m/n combination in subroutine get_Fi m=",m,"n=",current_var
  end if
  if((total_vars.lt.3).and.(solute).and.(current_var.eq.3))then
     current_var=2
  end if
  ! END USER EDIT
end subroutine get_stencil_2d




!!!!!!!!!!!   Now back to the main code

subroutine phase_field_check(mype,inumber)
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
  integer, intent(in) :: inumber,mype
  integer :: ierr,lb,i,j,k,mid_x,mid_y,mid_z, found
  double precision :: phase,umin,umax,phimin,phimax,u,phi,temp,gradphi,gradphi2,theta1,theta2
  double precision :: dx,dist,phi_val_1,phi_val_2,phi_val_3,phi_val_4
  double precision :: new_loc, new_rad, temp2
  double precision :: xl, xm, xr, phil, phim, phir, thetal, thetam, thetar, sigmastartemp
  logical :: on_axisx=.false.
!  character (len=19) :: filename
  !print *,"Guardcell axis node count:",axis_node_count
  !axis_node_count=0
  phase=0.0
  umin=1e30
  umax=-1e30
  phimin=1e30
  phimax=-1e30
  mid_x=il_bnd+nxb/2+nguard
  mid_y=jl_bnd+nyb/2+nguard
  if (ndim.eq.3) mid_z=kl_bnd+nzb/2+nguard
  found=0
  if(mype.eq.0)print *,"Phase field summary for time step: ",inumber
! CEG removed profiles
!  if(mype.gt.0)print *,"Warning, solute profiler not yet coded for parallelism"
!  if(mod(inumber,100).eq.0)then
!     write (filename,1001) "solute_profile_",inumber
!     if(mype.eq.0)open (unit=136,file=filename)
!     write (filename,1001) "phase_profile_",inumber
!     if(mype.eq.0)open (unit=137,file=filename)
!     write (filename,1001) "thermal_profile_",inumber
!     if(mype.eq.0)open (unit=138,file=filename)
!  end if
  

  dist = 0.0
  sigmastartemp = 0.0
  if (nvar.eq.2.and.solute) found=1
  do lb=1,lnblocks
     if(nodetype(lb).eq.1) then
        dx = bsize(1,lb)/real(nxb)
        ! First find out if we're "on axis"
        ! how do we know what the "axis" is?
        ! Simple the nucleate_(x,y,z) parameter holds the co-ords
        ! To be on an axis 2 out of 3 of the block faces must match these co-ords
        ! Only testing lower co-ords so for corner nucleation the nucleus MUST BE at min co-ords for each dimension
        on_axisx=.false.
        
        if(bnd_box(1,2,lb).eq.nucleate_y)then
           ! lower y bnd matches
           if(ndim.eq.2.or.bnd_box(1,3,lb).eq.nucleate_z)then
              ! lower z bnd matches
              ! y and z co-ords match so we're on the x axis
              on_axisx=.true.
              !axis_node_count=axis_node_count+1
           end if
        end if
        j=jl_bnd+nguard
        k=kl_bnd+nguard*k3d

        if(on_axisx)then
           ! On x axis now test whether the interface is within this block
           ! Note testing from il_bnd to iu_bnd
           ! This ensure that if interface is between 2 blocks we still find it
           ! Basically we're checking guardcells as well as computational domain
           if(unk(1,il_bnd+nguard,j,k,lb).gt.0.0.and.&
              unk(1,iu_bnd-nguard+1,j,k,lb).le.0.0)then
!              print *, lb, unk(1,il_bnd+nguard,j,k,lb),unk(1,iu_bnd-nguard+1,j,k,lb)
              call phase_loc_rad(lb, new_loc, new_rad)
              found=found+1

!              print *,'fnd', unk(1,iu_bnd,j,k,lb)

           end if

           if (unk(1,il_bnd+nguard,j,k,lb).gt.-0.9.and.&
               unk(1,iu_bnd-nguard+1,j,k,lb).le.-0.9)then
              if (nvar.eq.3.or.(.not.solute)) then
                 do i = il_bnd+nguard, iu_bnd-nguard
                    if (unk(1,i+1,j,k,lb).le.-0.9) exit
                 end do
                 ! Set up local variables to do this double variable interpolation, first in X for xm then in theta for thetam
                 xl = bnd_box(1,1,lb)+(i-1.5)*dx
                 xr = bnd_box(1,1,lb)+(i+1-1.5)*dx
                 phil = unk(1,i,j,k,lb)
                 phim = -0.9
                 phir = unk(1,i+1,j,k,lb)
                 thetal = unk(1+nunkvbles,i,j,k,lb)
                 thetar = unk(1+nunkvbles,i+1,j,k,lb)

                 xm = xl + (phim-phil)*(xr-xl)/(phir-phil)
                 thetam = thetal + (xm-xl)*(thetar-thetal)/(xr-xl)
                 found=found+1
              endif
           end if
        end if
     end if
!     if (found.eq.2) exit
  end do

!  print *,"Phase_field_check axis node count:",axis_node_count
  call MPI_REDUCE(new_rad,temp2,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  if(mype.eq.0)write (6,'(A,A)') 27, "[36;1;47m"
  if(mype.eq.0)print *,"Tip Radius",temp2
  if(mype.eq.0)write(unit=125,fmt=*) time, temp2, dt

  call MPI_REDUCE(dist,temp,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(new_loc,temp2,1,amr_mpi_real,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  if(mype.eq.0)print *,"Tip position",temp2
  if(mype.eq.0)write(unit=126,fmt=*) time, temp2
!  call MPI_REDUCE(phase,temp,1,amr_mpi_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!  if(mype.eq.0)print *,"Amount of solid",temp

  call MPI_REDUCE(thetam,temp,1,amr_mpi_real,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  if ((nvar.eq.3.or.(.not.solute)).and.mype.eq.0)print *,"Sigma star temperature", temp

  if(mype.eq.0)write (6,'(A,A)') 27, "[0m"
  if(mype.eq.0)call flush(125)
  if(mype.eq.0)call flush(126)
! CEG removed profiles and p2, p3, p4, out1, out2, out3
1001 format(A15,I4) 
end subroutine phase_field_check

! CEG routine for calculating tip location and radius using a series of quadratics

subroutine  phase_loc_rad(lb, cur_loc, cur_rad)
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  integer i, j, k, lb
  double precision cur_loc, cur_rad, x1a, x1b, x2a, x2b, a5, b1, b3, b5, dx
  double precision x0, x1, x2, xm1

  dx = bsize(1,lb)/real(nxb)

  j = jl_bnd+nguard
  k = kl_bnd+nguard*k3d
  do i=il_bnd+nguard,iu_bnd-nguard
     if (unk(1,i,j,k,lb).gt.0.0 .and. unk(1,i+1,j,k,lb).le.0.0) then
! Linearly interpolated points
        ! Get x1a - x-position of interface on line y=h/2,z=h/2 using x-y plane
!        x1a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k,lb)/(unk(1,i+1,j,k,lb)-unk(1,i,j,k,lb)))
        ! Get x1b - x-position of interface on line y=3h/2,z=h/2 using x-y plane
!        x1b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k,lb)/(unk(1,i+1,j+1,k,lb)-unk(1,i,j+1,k,lb)))
        ! Get x2a - x-position of interface on line y=h/2,z=3h/2 using x-y plane
!        x2a=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j,k+1,lb)/(unk(1,i+1,j,k+1,lb)-unk(1,i,j,k+1,lb)))
        ! Get x2b - x-position of interface on line y=3h/2,z=3h/2 using x-y plane
!        x2b=bnd_box(1,1,lb)+(i-1.5)*dx+dx*(-unk(1,i,j+1,k+1,lb)/(unk(1,i+1,j+1,k+1,lb)-unk(1,i,j+1,k+1,lb)))

        x1 = bnd_box(1,1,lb)+(i-1.5)*dx
        x0 = x1-dx
        xm1 = x0-dx
        x2 = x1+dx

        call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k,lb), unk(1,i,j,k,lb), unk(1,i+1,j,k,lb), x1a)
        if (unk(1,i,j+1,k,lb).gt.0.0.or.i.le.2) then
           call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), unk(1,i+1,j+1,k,lb), x1b)
        else
           call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k,lb), unk(1,i-1,j+1,k,lb), unk(1,i,j+1,k,lb), x1b)
        endif
        if (ndim.eq.3) then
           if (unk(1,i,j,k+1,lb).gt.0.0.or.i.le.2) then
              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), unk(1,i+1,j,k+1,lb), x2a)
           else
              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j,k+1,lb), unk(1,i-1,j,k+1,lb), unk(1,i,j,k+1,lb), x2a)
           endif
           if (unk(1,i,j+1,k+1,lb).gt.0.0.or.i.le.2) then
              call get_quadratic_interp(x0, x1, x2, unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), unk(1,i+1,j+1,k+1,lb), x2b)
           else
              call get_quadratic_interp(xm1, x0, x1, unk(1,i-2,j+1,k+1,lb), unk(1,i-1,j+1,k+1,lb), unk(1,i,j+1,k+1,lb), x2b)
           endif
        else
           x2a = x1a
           x2b = x1b
        endif
        exit
     endif
  enddo

  b1 = 0.125*(9.0*x1a-x2a)
  b3 = 0.125*(9.0*x1b-x2b)
  b5 = 0.125*(9.0*b1-b3)

  a5 = 0.5*(b3-b1)/(dx*dx)

  cur_loc=b5
! Wolfram-Alpha gives radius of curvature of y=f(x) as R=( (1+(y')^2)^1.5 )/|y''|
  cur_rad = 1.0/abs(2.0*a5)

end subroutine phase_loc_rad

subroutine get_quadratic_interp(x0, x1, x2, y0, y1, y2, xstar)
  double precision x0, x1, x2, y0, y1, y2, xstar, a, b, c

  b = ((y0-y2)*(x1*x1-x2*x2)-(y1-y2)*(x0*x0-x2*x2))/((x0-x2)*(x1*x1-x2*x2)-(x1-x2)*(x0*x0-x2*x2))
  a = (y0-b*(x0-x2)-y2) / (x0*x0-x2*x2)
  c = y2-a*x2*x2-b*x2

! Choose -b MINUS (...)/2a as want left sided root
  xstar = (-b-sqrt(b*b-4.0*a*c))/(2.0*a)

end subroutine get_quadratic_interp


subroutine get_w_aniso_2(phix,phiy,phiz,grad_phi_sq,A_func)
  use solution_parameters
  implicit none
  double precision, intent(out) :: A_func
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq
  if(epsilon_tilde.eq.0.0)then
     A_func=A_0
  else
     if(grad_phi_sq.gt.0)then
        A_func=A_0*(1+epsilon_tilde*(phix**4+phiy**4+phiz**4)/(grad_phi_sq*grad_phi_sq))
     else
        A_func=A_0
     end if
  end if
  if (abs(A_func).lt.1.0e-14) A_func=0.0
end subroutine get_w_aniso_2

subroutine get_dw_dphix(phix,phiy,phiz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*phix*(phix*phix*grad_phi_sq-phix**4-phiy**4-phiz**4)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)

end subroutine get_dw_dphix

subroutine get_dw_dphiy(phix,phiy,phiz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*phiy*(phiy*phiy*grad_phi_sq-phix**4-phiy**4-phiz**4)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphiy

subroutine get_dw_dphiz(phix,phiy,phiz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*phiz*(phiz*phiz*grad_phi_sq-phix**4-phiy**4-phiz**4)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphiz

subroutine get_dw_dx(phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix**3*phixx+phiy**3*phixy+phiz**3*phixz)&
             -(phix**4+phiy**4+phiz**4)*(phix*phixx+phiy*phixy+phiz*phixz))/grad_phi_sq**3
end subroutine get_dw_dx

subroutine get_dw_dy(phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix**3*phixy+phiy**3*phiyy+phiz**3*phiyz)&
             -(phix**4+phiy**4+phiz**4)*(phix*phixy+phiy*phiyy+phiz*phiyz))/grad_phi_sq**3
end subroutine get_dw_dy

subroutine get_dw_dz(phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phiz,phixx,phixy,phixz,phiyy,phiyz,phizz,grad_phi_sq
        val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix**3*phixz+phiy**3*phiyz+phiz**3*phizz)&
             -(phix**4+phiy**4+phiz**4)*(phix*phixz+phiy*phiyz+phiz*phizz))/grad_phi_sq**3
end subroutine get_dw_dz

subroutine get_dw_dpsi(sin_psi,cos_psi,sin_theta,cos_theta,dw_dpsi_val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: dw_dpsi_val
  double precision, intent(in) :: sin_psi,cos_psi,sin_theta,cos_theta
  if(epsilon_tilde.eq.0.0)then
     dw_dpsi_val=0.0
  else
     dw_dpsi_val=4.0*A_0*epsilon_tilde*cos_psi*sin_psi*(sin_psi**2*(1.0-2.0*sin_theta**2*cos_theta**2)-cos_psi**2)
  end if
end subroutine get_dw_dpsi

subroutine get_dunk_n_dx(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf.eqv..true.)then
     val = (work(i+1,j,k,lb,n)-work(i-1,j,k,lb,n))/(2.0*dx)
  else
     val = (unk(n,i+1,j,k,lb)-unk(n,i-1,j,k,lb))/(2.0*dx)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_dunk_n_dx

subroutine get_dunk_n_dy(n,i,j,k,dy,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dy
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j+1,k,lb,n)-work(i,j-1,k,lb,n))/(2.0*dy) 
  else
     val = (unk(n,i,j+1,k,lb)-unk(n,i,j-1,k,lb))/(2.0*dy)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_dunk_n_dy

subroutine get_dunk_n_dz(n,i,j,k,dz,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dz
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j,k+1,lb,n)-work(i,j,k-1,lb,n))/(2.0*dz) 
  else
     val = (unk(n,i,j,k+1,lb)-unk(n,i,j,k-1,lb))/(2.0*dz)
  end if
end subroutine get_dunk_n_dz

subroutine get_d2unk_n_dx2(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i+1,j,k,lb,n)-2.0*work(i,j,k,lb,n)+work(i-1,j,k,lb,n))/(dx*dx) 
  else
     val = (unk(n,i+1,j,k,lb)-2.0*unk(n,i,j,k,lb)+unk(n,i-1,j,k,lb))/(dx*dx)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_d2unk_n_dx2

subroutine get_d2unk_n_dy2(n,i,j,k,dy,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dy
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j+1,k,lb,n)-2.0*work(i,j,k,lb,n)+work(i,j-1,k,lb,n))/(dy*dy)
  else
     val = (unk(n,i,j+1,k,lb)-2.0*unk(n,i,j,k,lb)+unk(n,i,j-1,k,lb))/(dy*dy)
  end if
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_d2unk_n_dy2

subroutine get_d2unk_n_dz2(n,i,j,k,dz,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dz
  logical, intent(in) :: workf
  double precision, intent(out) :: val
  if(workf)then
     val = (work(i,j,k+1,lb,n)-2.0*work(i,j,k,lb,n)+work(i,j,k-1,lb,n))/(dz*dz)
  else
     val = (unk(n,i,j,k+1,lb)-2.0*unk(n,i,j,k,lb)+unk(n,i,j,k-1,lb))/(dz*dz)
  end if
end subroutine get_d2unk_n_dz2

subroutine get_d2unk_n_dxdy(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  double precision, intent(out) :: val
  logical, intent(in) :: workf
  double precision :: temp1,temp2
  call get_dunk_n_dx(n,i,j+1,k,dx,lb,temp1,workf)
  call get_dunk_n_dx(n,i,j-1,k,dx,lb,temp2,workf)
  val = (temp1-temp2)/(2.0*dx)
  if (abs(val).lt.1.0e-14) val=0.0
end subroutine get_d2unk_n_dxdy

subroutine get_d2unk_n_dxdz(n,i,j,k,dx,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dx
  double precision, intent(out) :: val
  logical, intent(in) :: workf
  double precision :: temp1,temp2
  call get_dunk_n_dx(n,i,j,k+1,dx,lb,temp1,workf)
  call get_dunk_n_dx(n,i,j,k-1,dx,lb,temp2,workf)
  val = (temp1-temp2)/(2.0*dx)
end subroutine get_d2unk_n_dxdz

subroutine get_d2unk_n_dydz(n,i,j,k,dy,lb,val,workf)
  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  implicit none
  integer, intent(in) :: n,i,j,k,lb
  double precision, intent(in) :: dy
  double precision, intent(out) :: val
  logical, intent(in) :: workf
  double precision :: temp1,temp2
  call get_dunk_n_dy(n,i,j,k+1,dy,lb,temp1,workf)
  call get_dunk_n_dy(n,i,j,k-1,dy,lb,temp2,workf)
  val = (temp1-temp2)/(2.0*dy)
end subroutine get_d2unk_n_dydz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Alternative versions for 2-d case                                  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_w_aniso_2_2d(phix,phiy,grad_phi_sq,A_func)
  use solution_parameters
  implicit none
  double precision, intent(out) :: A_func
  double precision, intent(in) :: phix,phiy,grad_phi_sq
  if(epsilon_tilde.eq.0.0)then
     A_func=A_0
  else
     if(grad_phi_sq.gt.0)then
        A_func=A_0*(1+epsilon_tilde*(phix*phix*phix*phix+phiy*phiy*phiy*phiy)/(grad_phi_sq*grad_phi_sq))
     else
        A_func=A_0
     end if
  end if
end subroutine get_w_aniso_2_2d

subroutine get_dw_dphix_2d(phix,phiy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,grad_phi_sq
  val=4.0*A_0*epsilon_tilde*phix*(phix*phix*grad_phi_sq-phix**4-phiy**4)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphix_2d

subroutine get_dw_dphiy_2d(phix,phiy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,grad_phi_sq
  val=4.0*A_0*epsilon_tilde*phiy*(phiy*phiy*grad_phi_sq-phix**4-phiy**4)/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dphiy_2d

subroutine get_dw_dx_2d(phix,phiy,phixx,phixy,phiyy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phixx,phixy,phiyy, grad_phi_sq
  val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix*phix*phix*phixx + phiy*phiy*phiy*phixy)&
        -(phix*phix*phix*phix+phiy*phiy*phiy*phiy)*(phix*phixx+phiy*phixy))/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dx_2d

subroutine get_dw_dy_2d(phix,phiy,phixx,phixy,phiyy,grad_phi_sq,val)
  use solution_parameters
  implicit none
  double precision, intent(out) :: val
  double precision, intent(in) :: phix,phiy,phixx,phixy,phiyy,grad_phi_sq
  val=4.0*A_0*epsilon_tilde*(grad_phi_sq*(phix*phix*phix*phixy+phiy*phiy*phiy*phiyy)&
             -(phix*phix*phix*phix+phiy*phiy*phiy*phiy)*(phix*phixy+phiy*phiyy))/(grad_phi_sq*grad_phi_sq*grad_phi_sq)
end subroutine get_dw_dy_2d
double precision function Lap2d(dx,d_flag,i,j,k,lb,n)
 use solution_parameters
  use time_dep_parameters
  use multigrid_parameters
  use paramesh_dimensions
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
integer, intent (in) :: i,j,k,lb,n  
real, intent (in):: dx
logical, intent(in) :: d_flag
real Lap_normal,Lap_diag,dLap,alpha 
!alpha = 0 5_point,alpha=2/3 = 9_point, alpha=1 diagonal Lap  
alpha = 1./3.
if(d_flag) then  
 Lap2d = (2.0*alpha-4.0)/(dx*dx) 
else  
 Lap_normal =   unk(1+n,i+1,j,  k,lb)&
               +unk(1+n,i  ,j+1,k,lb)&
               +unk(1+n,i-1,j  ,k,lb)&
               +unk(1+n,i  ,j-1,k,lb)&
               -4.0* unk(1+n,i ,j,k,lb)
 Lap_diag =     unk(1+n,i+1,j+1,k,lb)&
               +unk(1+n,i-1,j+1,k,lb)&
               +unk(1+n,i-1,j-1,k,lb)&
               +unk(1+n,i+1,j-1,k,lb)&
               -4.0* unk(1+n,i ,j,k,lb)
 Lap_diag = Lap_diag*0.5

 Lap2d = ((1.0-alpha)*Lap_normal+alpha*Lap_diag)/(dx*dx)
end if
end function Lap2d
double precision function A_func2d(X,Y)
  use solution_parameters
  implicit none
  double precision, intent(in) :: X,Y
  double precision :: p,q

  if(epsilon_tilde.eq.0.0)then
     A_func2d=A_0
  else
     if(X+Y.gt.0)then
       p=X/(X+Y)
       q=Y/(X+Y)
       A_func2d=A_0*(1.0+epsilon_tilde*(1.0-2.*p*q))
     else
        A_func2d=A_0
     end if
  end if
  if (abs(A_func2d).lt.1.0e-14) then
  write(*,*)"In A_func2d"
  stop
  endif
end function A_func2d
double precision function A_func3d(X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     A_func3d=A_0
  else
     A_func3d=A_0*(1.0+epsilon_tilde*(X*X+Y*Y+Z*Z))
  end if
end function A_func3d
double precision function Ai_func3d(phix,X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: phix,X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     Ai_func3d=0.0
  else
     Ai_func3d=4.0*A_0*epsilon_tilde*phix*(X*(Y+Z)-Y*Y-Z*Z)
  end if
end function Ai_func3d
double precision function Aij_func3d(phix,phiy,X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: phix,phiy,X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     Aij_func3d=0.0
  else
     Aij_func3d=8.0*A_0*epsilon_tilde*phix*phiy*(X*X+Y*Y-4.0*X*Y-2.0*Z*(X+Y)+3.0*Z*Z)
  end if
end function Aij_func3d
double precision function Aii_func3d(X,Y,Z)
  use solution_parameters
  implicit none
  double precision, intent(in) :: X,Y,Z

  if(epsilon_tilde.eq.0.0)then
     Aii_func3d=0.0
  else
     Aii_func3d=-4.0*A_0*epsilon_tilde*(3.0*X*X*(Y+Z)-8.0*X*(Y*Y+Z*Z)-6.0*X*Y*Z+Y*Y*Y+Y*Y*Z+Z*Z*Y+Z*Z*Z)
  end if
end function Aii_func3d
