      function ftheta(theta,phi) 
      use constants
      implicit none
      real :: ftheta,theta,phi

      ftheta = pi*.5 - theta
      if(theta.le.pi*.25) then
         ftheta = theta
      elseif(theta.ge.pi*.75) then
         ftheta = pi-theta
      endif
      if(phi.gt.pi) ftheta = - ftheta

      end function ftheta

