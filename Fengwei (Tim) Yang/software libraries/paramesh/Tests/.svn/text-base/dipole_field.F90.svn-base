      Subroutine dipole_field(x, y, z, Ex, Ey, Ez)

      ! computes and returns field of a dipole located at the origin with
      ! dipole moment pointing in the y direction

      Real, Intent(in)  :: x, y, z
      Real, Intent(out) :: Ex, Ey, Ez

      Real :: r
      Real,parameter :: p = .1

!      Real,parameter :: xdi = -0.5, ydi = .5, zdi = .5
      Real,parameter :: xdi = 0., ydi = 0., zdi = 0.

      xx = x - xdi
      yy = y - ydi
      zz = z - zdi

      r = xx**2 + yy**2 + zz**2
      r = sqrt(r)
      if (r <= 1.e-10) then
         print *,' ERROR: r too small in dipole_field'
      end if

!      Ex = p*yy*xx*3./(r**5)
!      Ey = (p*yy*yy*3./(r**5)) - (p/(r**3))
!      Ez = p*yy*zz*3./(r**5)

      n = 3
      Ex = (xx**n)*(yy**(n-1))*(zz**(n-1))
      Ey = -0.5*(xx**(n-1))*(yy**n)*(zz**(n-1))
      Ez = -0.5*(xx**(n-1))*(yy**(n-1))*(zz**n)

      Return
      End Subroutine dipole_field
      
