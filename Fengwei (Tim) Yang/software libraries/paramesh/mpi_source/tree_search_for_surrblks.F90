      subroutine tree_search_for_surrblks ()

      use local_tree_common
      use physicaldata
      use tree
      use paramesh_dimensions
      use paramesh_interfaces
      use paramesh_comm_data
      use mpi_morton, only : lperiodicx, lperiodicy, lperiodicz
      use constants

      implicit none

      include 'mpif.h'

      real :: neigh_coord(3), neigh_coord2(3)
      integer :: ierr, nprocs, mype!, convex
      integer :: neigh_lb,neigh_proc, neigh_nodetype
      integer :: i,j,k,ii,jj,kk,iboun,lb

      logical :: found

      real :: time_exe, time_max
      real :: eps,accuracy

      accuracy = 100./10.**precision(accuracy)
      if (accuracy > 1.0e-10) then
         eps = 1.e-10
      else
         eps = accuracy
      end if

      time_exe = MPI_WTIME()

      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

!      print*,boundary_box(:,:,:)

!  convex = 1
! FIND SURROUNDING BLOCKS
      do lb = 1, lnblocks

         kk = -k3d
         do k = 1,1+2*k3d

            neigh_coord(3) = coord(3,lb) + kk*bsize(3,lb)
            if (lperiodicz.and.neigh_coord(3).lt.grid_zmin)               &
                neigh_coord(3) = neigh_coord(3) + (grid_zmax-grid_zmin)
            if (lperiodicz.and.neigh_coord(3).gt.grid_zmax)               &
                neigh_coord(3) = neigh_coord(3) - (grid_zmax-grid_zmin)

            jj = -k2d
            do j = 1,1+2*k2d

               neigh_coord(2) = coord(2,lb) + jj*bsize(2,lb)
               if (lperiodicy.and.neigh_coord(2).lt.grid_ymin)               &
                  neigh_coord(2) = neigh_coord(2) + (grid_ymax-grid_ymin)
               if (lperiodicy.and.neigh_coord(2).gt.grid_ymax)               &
                  neigh_coord(2) = neigh_coord(2) - (grid_ymax-grid_ymin)

               ii = -1
               do i = 1,3

                  neigh_coord(1) = coord(1,lb) + ii*bsize(1,lb)
                  if (lperiodicx.and.neigh_coord(1).lt.grid_xmin)               &
                     neigh_coord(1) = neigh_coord(1) + (grid_xmax-grid_xmin)
                  if (lperiodicx.and.neigh_coord(1).gt.grid_xmax)               &
                     neigh_coord(1) = neigh_coord(1) - (grid_xmax-grid_xmin)

                  neigh_coord2(:) = neigh_coord(:)

!--------Reset coordinates of neighbor in spherical coordinates
                  If (spherical_pm) Then
                     If ( ((jj == -1).and.(abs(bnd_box(1,2,lb)) < eps)) .or.    &  
                         ((jj ==  1).and.(abs(bnd_box(2,2,lb)-pi) < eps))      & 
                         .and. lsingular_line ) Then
                        neigh_coord2(2) = coord(2,lb)
                        If (neigh_coord2(3) < pi) Then
                            neigh_coord2(3) = neigh_coord2(3) + pi
                           Elseif(neigh_coord2(3) > pi)  Then
                           neigh_coord2(3) = neigh_coord2(3) - pi
                        End If
                     End If
                  End If

                  neigh_lb = -1
                  neigh_proc = -1
                  neigh_nodetype = -1
                  found = .false.

!--------First check boundaries
!CEG : Needs modification for non-convex domains
                  ! This bit checks the external bounding box
                  do iboun = 1, nboundaries
                     if (ndim == 1) then
                        if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                            neigh_coord2(1) < boundary_box(2,1,iboun)) then
                           found = .true.
                           neigh_lb = boundary_index(iboun)
                           neigh_proc = boundary_index(iboun)
                        end if
                     elseif (ndim == 2) then
!            if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
!                neigh_coord2(1) < boundary_box(2,1,iboun) .and.        &
!                neigh_coord2(2) > boundary_box(1,1+k2d,iboun) .and.    &
!                neigh_coord2(2) < boundary_box(2,1+k2d,iboun)) then
!CEG rewriting
!               if (neigh_coord2(1).lt.1.0) then
!                   print*,lb, iboun, neigh_coord2(1), boundary_box(1,1,iboun)!, boundary_box(2,1,iboun)
!               endif
                        if (neigh_coord2(1) > boundary_box(1,1,iboun)) then
                           if(neigh_coord2(1) < boundary_box(2,1,iboun)) then
                              if(neigh_coord2(2) > boundary_box(1,1+k2d,iboun)) then
                                 if(neigh_coord2(2) < boundary_box(2,1+k2d,iboun)) then
                                    found = .true.
                                    neigh_lb = boundary_index(iboun)
                                    neigh_proc = boundary_index(iboun)
                                 end if
                              end if
                           end if
                        end if

                     elseif (ndim == 3) then
                        if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                            neigh_coord2(1) < boundary_box(2,1,iboun) .and.        &
                            neigh_coord2(2) > boundary_box(1,1+k2d,iboun) .and.    &
                            neigh_coord2(2) < boundary_box(2,1+k2d,iboun) .and.    &
                            neigh_coord2(3) > boundary_box(1,1+2*k3d,iboun) .and.  &
                            neigh_coord2(3) < boundary_box(2,1+2*k3d,iboun)) then
                           found = .true.
                           neigh_lb = boundary_index(iboun)
                           neigh_proc = boundary_index(iboun)
                        end if
                     end if
                  end do

!if (convex.eq.1)print *,'Did convex'
!  convex = convex + 1
!if (3.gt.2) then
!CEG added this bit for non-convex domains
                  if (lrefine(lb).eq.1) then
                     if (ndim.eq.2.or.k.eq.2) then
                        if (i.eq.1.and.j.eq.1) then!
                           if (neigh(1, 1, lb).le.-20.and.neigh(1, 3, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 1, lb).lt.neigh(1, 3, lb)) then 
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
                              else
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.1) then!
                           if (neigh(1, 3, lb).le.-20) then
                              found = .true.
                              neigh_lb = neigh(1, 3, lb)
                              neigh_proc = neigh(2, 3, lb)
                           endif
                        else if (i.eq.3.and.j.eq.1) then!
                           if (neigh(1, 2, lb).le.-20.and.neigh(1, 3, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 2, lb).lt.neigh(1, 3, lb)) then 
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              endif
                           endif

                        else if (i.eq.1.and.j.eq.2) then!
                           if (neigh(1, 1, lb).le.-20) then
                              found = .true.
                              neigh_lb = neigh(1, 1, lb)
                              neigh_proc = neigh(2, 1, lb)
                           endif
! This one!                     else if (i.eq.2.and.j.eq.2) then!
                        else if (i.eq.3.and.j.eq.2) then!
                           if (neigh(1, 2, lb).le.-20) then
                              found = .true.
                              neigh_lb = neigh(1, 2, lb)
                              neigh_proc = neigh(2, 2, lb)
                           endif

                        else if (i.eq.1.and.j.eq.3) then!
                           if (neigh(1, 1, lb).le.-20.and.neigh(1, 4, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 1, lb).lt.neigh(1, 4, lb)) then 
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
                              else
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.3) then!
                           if (neigh(1, 4, lb).le.-20) then
                              found = .true.
                              neigh_lb = neigh(1, 4, lb)
                              neigh_proc = neigh(2, 4, lb)
                           endif
                        else if (i.eq.3.and.j.eq.3) then!
                           if (neigh(1, 2, lb).le.-20.and.neigh(1, 4, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 2, lb).lt.neigh(1, 4, lb)) then 
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              endif
                           endif
                        endif
                     ! end ndim.eq.2 or k.eq.2
                     !!!!!!!!!!!!!!!!!!!!!!!!!
                     else if (k.eq.1) then!
                        if (i.eq.1.and.j.eq.1) then
                           if (neigh(1, 1, lb).le.-20.and.neigh(1, 3, lb).le.-20.and.neigh(1, 5, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 1, lb).lt.neigh(1, 3, lb).and.neigh(1, 1, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
!!fwy                              else if (neigh(1, 3, lb).and.neigh(1, 1, lb).lt.neigh(1, 5, lb)) then 
                              else if (neigh(1, 3, lb).lt.neigh(1, 5, lb).and.neigh(1, 1, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              else
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.1) then!
                           if (neigh(1, 5, lb).le.-20.and.neigh(1, 3, lb).le.-20) then
                              if (neigh(1, 3, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif
                        else if (i.eq.3.and.j.eq.1) then!
                           if (neigh(1, 2, lb).le.-20.and.neigh(1, 3, lb).le.-20.and.neigh(1, 5, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 2, lb).lt.neigh(1, 3, lb).and.neigh(1, 2, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else if (neigh(1, 3, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              else
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif

                        else if (i.eq.1.and.j.eq.2) then!
                           if (neigh(1, 5, lb).le.-20.and.neigh(1, 1, lb).le.-20) then
                              if (neigh(1, 1, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.2) then!
                           if (neigh(1, 5, lb).le.-20) then
                              found = .true.
                              neigh_lb = neigh(1, 5, lb)
                              neigh_proc = neigh(2, 5, lb)
                           endif
                        else if (i.eq.3.and.j.eq.2) then!
                           if (neigh(1, 5, lb).le.-20.and.neigh(1, 2, lb).le.-20) then
                              if (neigh(1, 2, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif

                        else if (i.eq.1.and.j.eq.3) then!
                           if (neigh(1, 1, lb).le.-20.and.neigh(1, 4, lb).le.-20.and.neigh(1, 5, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 1, lb).lt.neigh(1, 4, lb).and.neigh(1, 1, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
                              else if (neigh(1, 4, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              else
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.3) then!
                           if (neigh(1, 5, lb).le.-20.and.neigh(1, 4, lb).le.-20) then
                              if (neigh(1, 4, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif
                        else if (i.eq.3.and.j.eq.3) then!
                           if (neigh(1, 2, lb).le.-20.and.neigh(1, 4, lb).le.-20.and.neigh(1, 5, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 2, lb).lt.neigh(1, 4, lb).and.neigh(1, 2, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else if (neigh(1, 4, lb).lt.neigh(1, 5, lb)) then 
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              else
                                 neigh_lb = neigh(1, 5, lb)
                                 neigh_proc = neigh(2, 5, lb)
                              endif
                           endif
                        endif
                     ! end k=1
                     !!!!!!!!!!!!!!!!!!!!!!!!!
                     else ! k=3
                        if (i.eq.1.and.j.eq.1) then
                           if (neigh(1, 1, lb).le.-20.and.neigh(1, 3, lb).le.-20.and.neigh(1, 6, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 1, lb).lt.neigh(1, 3, lb).and.neigh(1, 1, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
                              else if (neigh(1, 3, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              else
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.1) then!
                           if (neigh(1, 6, lb).le.-20.and.neigh(1, 3, lb).le.-20) then
                              if (neigh(1, 3, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif
                        else if (i.eq.3.and.j.eq.1) then!
                           if (neigh(1, 2, lb).le.-20.and.neigh(1, 3, lb).le.-20.and.neigh(1, 6, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 2, lb).lt.neigh(1, 3, lb).and.neigh(1, 2, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else if (neigh(1, 3, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 3, lb)
                                 neigh_proc = neigh(2, 3, lb)
                              else
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif

                        else if (i.eq.1.and.j.eq.2) then!
                           if (neigh(1, 6, lb).le.-20.and.neigh(1, 1, lb).le.-20) then
                              if (neigh(1, 1, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.2) then!
                           if (neigh(1, 6, lb).le.-20) then
                              found = .true.
                              neigh_lb = neigh(1, 6, lb)
                              neigh_proc = neigh(2, 6, lb)
                           endif
! To be inserted for side face
                        else if (i.eq.3.and.j.eq.2) then!
                           if (neigh(1, 6, lb).le.-20.and.neigh(1, 2, lb).le.-20) then
                              if (neigh(1, 2, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif

                        else if (i.eq.1.and.j.eq.3) then
                           if (neigh(1, 1, lb).le.-20.and.neigh(1, 4, lb).le.-20.and.neigh(1, 6, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 1, lb).lt.neigh(1, 4, lb).and.neigh(1, 1, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 1, lb)
                                 neigh_proc = neigh(2, 1, lb)
                              else if (neigh(1, 4, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              else
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif
                        else if (i.eq.2.and.j.eq.3) then!
                           if (neigh(1, 6, lb).le.-20.and.neigh(1, 4, lb).le.-20) then
                              if (neigh(1, 4, lb).le.-20) then
                                 found = .true.
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              else
                                 found = .true.
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif
                        else if (i.eq.3.and.j.eq.3) then!
                           if (neigh(1, 2, lb).le.-20.and.neigh(1, 4, lb).le.-20.and.neigh(1, 6, lb).le.-20) then
                              found = .true.
                              if (neigh(1, 2, lb).lt.neigh(1, 4, lb).and.neigh(1, 2, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 2, lb)
                                 neigh_proc = neigh(2, 2, lb)
                              else if (neigh(1, 4, lb).lt.neigh(1, 6, lb)) then 
                                 neigh_lb = neigh(1, 4, lb)
                                 neigh_proc = neigh(2, 4, lb)
                              else
                                 neigh_lb = neigh(1, 6, lb)
                                 neigh_proc = neigh(2, 6, lb)
                              endif
                           endif
                        endif
                     endif
! end k=3
                  endif
!endif
                                    

! CEG If BC is boundary then change from -21 to an appropriate pid
                  if (neigh_lb.le.-20) neigh_proc = mype


!if (lrefine(lb).eq.1) print *,lb,iboun, neigh(:, iboun, lb)
                  call search_sub_tree(local_tree,neigh_coord2,lrefine(lb),     &
                              neigh_lb,neigh_proc,neigh_nodetype,found)
                  surr_blks(1,i,j,k,lb) = neigh_lb
                  surr_blks(2,i,j,k,lb) = neigh_proc
                  surr_blks(3,i,j,k,lb) = neigh_nodetype

!                  if (lrefine(lb).eq.1.and.surr_blks(1,i,j,k,lb).eq.-1) surr_blks(1,i,j,k,lb) = -21

                  ii = ii + 1
               end do
               jj = jj + k2d
            end do
         kk = kk + k3d
         end do

      end do

      end subroutine tree_search_for_surrblks
