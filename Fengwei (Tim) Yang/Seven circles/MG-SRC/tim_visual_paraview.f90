! Tim Yang
! 02/04/2013
! This subroutine uses amr_infinity_norm as a template,
! to output coordinates of nodes that are on finest grid.
! It is then used to visualize the mesh and its adaptivity.
#include "paramesh_preprocessor.fh"
      subroutine tim_visual_paraview(output_timestep)
!
!--------------------------------------------------------------
! include files for amr
      use paramesh_dimensions
      use physicaldata
      use tree
      implicit none
      include 'mpif.h'
      integer :: nguard0, lb, k, j, i
      integer, intent(in) :: output_timestep
!--------------------------------------------------------------
      double precision :: xi, yi, dx, zi
      character(len=80) :: file_write, file_write2, file_write3
      nguard0 = nguard*npgs

      if(output_timestep.le.9)then
668     format(A2,I1,A4)
        write(file_write,668) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,668) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,668) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.9.and.output_timestep.le.99)then
669     format(A2,I2,A4)
        write(file_write,669) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,669) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,669) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.99.and.output_timestep.le.999)then
670     format(A2,I3,A4)
        write(file_write,670) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,670) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,670) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.999.and.output_timestep.le.9999)then
671     format(A2,I4,A4)
        write(file_write,671) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,671) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,671) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.9999.and.output_timestep.le.99999)then
672     format(A2,I5,A4)
        write(file_write,672) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,672) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,672) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.99999.and.output_timestep.le.999999)then
673     format(A2,I6,A4)
        write(file_write,673) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,673) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,673) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.999999.and.output_timestep.le.9999999)then
674     format(A2,I7,A4)
        write(file_write,674) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,674) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,674) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.9999999.and.output_timestep.le.99999999)then
675     format(A2,I8,A4)
        write(file_write,675) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,675) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,675) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      elseif(output_timestep.gt.99999999.and.output_timestep.le.999999999)then
676     format(A2,I9,A4)
        write(file_write,676) '1_',output_timestep,'.csv'
        open(unit = 10011, file=file_write)
!!        write(file_write2,676) '2_',output_timestep,'.csv'
!!        open(unit = 10021, file=file_write2)
!!        write(file_write3,676) '3_',output_timestep,'.csv'
!!        open(unit = 10031, file=file_write3)
      else
        write(*,*) "ERROR to form file name."
      endif


667   format(f12.9,',',f12.9,',',f12.9,',',f21.16,',',f21.16,',',f21.16,',',f21.16)
      write(10011,*) 'x coord, y coord, z coord, scalar, pre, prepre, mu'
!!      write(10021,*) 'x coord, y coord, z coord, scalar'
!!      write(10031,*) 'x coord, y coord, z coord, scalar'
      if(lnblocks.gt.0) then
        do lb=1,lnblocks
          dx = bsize(1,lb)/real(nxb)
          !! nodetype = 1 when on finest grid.
          if(nodetype(lb).eq.1) then
            do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
              do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                do i=il_bnd+nguard,iu_bnd-nguard
                  xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
                  yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
                  zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5)
                  write(10011,667) xi,yi,zi,unk(1,i,j,k,lb),unk(4,i,j,k,lb),unk(5,i,j,k,lb),unk(6,i,j,k,lb)
!!                  write(10021,667) xi,yi,zi,unk(6,i,j,k,lb)
!!                  write(10031,667) xi,yi,zi,unk(11,i,j,k,lb)
                enddo ! do i
              enddo ! do j
            enddo ! do k
          endif ! if(nodetype(lb).eq.1 .or. advance_all_levels) then
        enddo ! do lb=1, lnblocks
      endif ! if(lnblocks.gt.0) then
  
      close(10011)
!!      close(10021)
!!      close(10031)
      write(*,*) ""
      write(*,*) "----------visual file output----------"
      write(*,*) ""
      return
      end subroutine tim_visual_paraview
