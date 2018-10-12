      subroutine convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)

! Converts spherical coords (r,theta,phi) into 3D cartesian coords (x,y,z)

      real, intent(out)  :: x,y,z
      real, intent(in) :: r,theta,phi


      sint = sin(theta)
      cost = cos(theta)
      sinp = sin(phi)
      cosp = cos(phi)

      x = r*sint*cosp
      y = r*sint*sinp
      z = r*cost

      return
      end
