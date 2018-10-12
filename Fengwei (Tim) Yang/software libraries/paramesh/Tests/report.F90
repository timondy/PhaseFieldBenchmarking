      subroutine report(mype,nprocs,ndim,l_no_perm,l_advance,l_mpi, & 
     &                  errors,lrelevant,executable)


      use io

      implicit none

      integer :: mype,nprocs,ndim
      logical :: l_no_perm,l_advance,l_mpi,lrelevant
      integer :: errors
      character*30 :: executable
      character*17 :: crelevant='test not relevant'
      logical :: lexist
      character (len=80) filename


      if(mype.eq.0) then

        filename = trim(output_dir) // 'testn.log'
        inquire(file=filename, exist=lexist)
        if(lexist) then
          open(unit=50,status='old',file=filename,position='append')

        else
          open(unit=50,status='new',file=filename)
 
          write(50,*) 'Test                            ', & 
     &                'ndim  nprocs  No_perm Adv_all MPI    Errors'
          write(50,*) '------------------------------------', & 
     &                '-----------------------------------------'
          write(50,*) ' '

        endif
        if(lrelevant) then
          write(50,100) executable,ndim,nprocs,l_no_perm, & 
     &                  l_advance,l_mpi,errors
100       format(a30,2x,i2,3x,i3,7x,l1,7x,l1,7x,l1,4x,i8)
        else
          write(50,101) executable,ndim,nprocs,l_no_perm, & 
     &                  l_advance,l_mpi,crelevant
101       format(a30,2x,i2,3x,i3,7x,l1,7x,l1,7x,l1,4x,a17)
        endif

      close(unit=50)
      endif

      return
      end subroutine report
