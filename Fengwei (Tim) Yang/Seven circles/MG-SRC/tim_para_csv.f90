! Tim Yang
! 02/04/2013
! This subroutine uses amr_infinity_norm as a template,
! to output coordinates of nodes that are on finest grid.
! It is then used to visualize the mesh and its adaptivity.
#include "paramesh_preprocessor.fh"
      subroutine tim_para_csv(pid, noprocs, output_timestep)
!
!--------------------------------------------------------------
! include files for amr
      use paramesh_dimensions
      use physicaldata
      use tree
      use multigrid_parameters
      use paramesh_interfaces
      implicit none
      include 'mpif.h'
      integer :: nguard0, lb, k, j, i, x, ierr
      integer, intent(in) :: pid, noprocs, output_timestep
!--------------------------------------------------------------
      double precision :: xi, yi, dx, zi
      character(len=80) :: file_write
      integer :: output_digit=1
      LOGICAL :: file_exists
      character(len=15) :: FMT

      nguard0 = nguard*npgs


      if(pid.eq.0) write(*,*) "parallel output from ", noprocs, "cores"

      output_digit = 1
      if(output_timestep.gt.9) output_digit = 2
      if(output_timestep.gt.99) output_digit = 3
      if(output_timestep.gt.999) output_digit = 4
      if(output_timestep.gt.9999) output_digit = 5
      if(output_timestep.gt.99999) output_digit = 6
      if(output_timestep.gt.999999) output_digit = 7
      write(FMT,'("I", I0, ")")') output_digit
!!112   format(A2,I<output_digit>,A4)
      write(file_write,FMT) '1_',output_timestep,'.csv'

      if(pid.eq.0)then
        inquire(file=file_write, exist=file_exists)
        if(file_exists.eqv..true.)then
          open(unit=1234, file=file_write)
          close(1234, status='delete')
        endif
      endif

667   format(f12.9,',',f12.9,',',f12.9,',',f21.16)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do x=1,noprocs
        if(pid+1.eq.x)then
          write(*,*) "core", pid, "reporting for duty."
          open(unit = 113, file=file_write,&
               position="append", action="write")
          if(pid.eq.0) write(113,*) 'x coord, y coord, z coord, scalar'
          do lb=1,lnblocks
            dx = bsize(1,lb)/real(nxb)
            if(nodetype(lb).eq.1) then
              do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d
                do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d
                  do i=il_bnd+nguard,iu_bnd-nguard
                    xi = bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
                    yi = bnd_box(1,2,lb) + dx*(real(j-nguard0)-.5)
                    zi = bnd_box(1,3,lb) + dx*(real(k-nguard0)-.5)
                    write(113,667) xi,yi,zi,unk(1,i,j,k,lb)
                  enddo ! do i
                enddo ! do j
              enddo ! do k
            endif ! if(nodetype(lb).eq.1 .or. advance_all_levels) then
          enddo ! do lb=1, lnblocks
          if(pid+1.eq.x) close(113)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      enddo ! do x=1,noprocs
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  
      if(pid.eq.noprocs-1) write(*,*) "----------visual file output:----------"
      return
      end subroutine tim_para_csv
