!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/local_tree_module
!! NAME
!!
!!   local_tree_module
!!
!! SYNOPSIS
!!
!!   Use local_tree_module
!!
!! DESCRIPTION
!!
!!   This file contains a module which contains the 'local' tree data structure,
!!   the subroutine for adding nodes to the local tree and a subroutine for
!!   search the local tree for surrounding blocks.
!!
!!   Here a 'local' tree means a tree data structured stored on the local 
!!   processor.  Each local tree is unique to that processor and contains 
!!   nodes necessary for the local blocks to find their surrounding blocks.
!!   THIS TREE IS NOT THE SAME AS THAT STORED IN THE module 'tree'.
!!   The tree is constructed using derived data types which utilize fortran
!!   pointers to child and parent nodes.
!!
!! AUTHORS
!!
!!   Kevin Olson
!!
!!***

      Module local_tree_module

!-----The derived data type for a pointer to a tree node.
      Type ptr_to_node
        type(node), pointer :: ptr
      End Type ptr_to_node

!-----The derived data type for a tree node.
      Type node
!-------Pointers to children
!        Type(ptr_to_node), allocatable, dimension(:) :: child
        Type(ptr_to_node), pointer, dimension(:) :: child
!-------Pointer to parent
        Type(ptr_to_node) :: parent
!-------Level in the tree.
        Integer :: level
!-------Block no., processor and nodetype of the block a tree node
        Integer :: lb, proc, nodetype
!-------Number of children
        Integer :: nchild
!-------Physical coodinates are coordinates of tree nodes.
        Real, dimension(3) :: coord  
!-------Physical size of tree nodes.
        Real, dimension(3) :: size   
!-------A logical which indicates if node is local to this processor or not
        Logical :: local
!-------A logical which indicates if node has been searched.
        Logical :: searched
      End Type node

      Contains

!!****f* mpi_source/add_block_to_tree
!! NAME
!!
!!   add_block_to_tree
!!
!! SYNOPSIS
!!
!!   Call Subroutine add_block_to_tree (t, tpar, top,            
!!           coord,bsize,level,                                         
!!           lb, proc, nodetype,                                         
!!           which_child, n_children_at_0, no_blocks_at_0,               
!!           off_proc, searched)
!!   Call Subroutine add_block_to_tree (pointer, pointer, pointer,            
!!           real,real,integer,                                         
!!           integer, integer, integer,                                         
!!           integer, integer, integer,               
!!           logical, logical)
!!
!! ARGUMENTS
!!
!!    t -> A pointer to the (sub)tree 'node' being tested and potentially
!!         which the block is added to.
!!    tpar -> A pointer to the parent node of 't'.
!!    top -> A pointer to the very top of the local tree being searched.
!!    coord -> physical coordinates of block searching the local tree.
!!    bsize -> physical size of block searching the local tree.
!!    level -> level in local tree where the block currently is.
!!    lb -> block id of the searching block.
!!    proc -> processor id of the home processor of the searching block.
!!    nodetype -> nodetype of the searching block.
!!    which_child -> the child number of the current node being searched.
!!    n_children_at_0 -> no. of children at level zero in the tree.
!!                       (level zero has one and only one node, but it can have
!!                        any number of children to accomodate meshes that
!!                        have a number of  level 1 blocks > 1.
!!    no_block_at_0 -> no. of blocks at level zero.
!!    off_proc -> a logical which indicates of the searching block is from a
!!                different processor than where this local tree is stored.
!!    searched -> a logical that indicates if this subtree node has been 
!!                searched of not, helps in pruning.
!!
!!
!! INCLUDES
!!
!!   mpif.h
!!
!! USES
!!
!!   tree
!!   paramesh_dimensions
!!   physicaldata
!!   constants
!!
!! CALLS
!!
!!   self
!!   search_and_prune_local_tree
!!
!! RETURNS
!!
!!   The number of children at level 0.  One return the block has been
!!   added.
!!
!! DESCRIPTION
!!
!!   The routine adds a Paramesh block to a local tree.  It recursively 
!!   searches the tree to find where the block belongs.  Also, the block
!!   is only added if if can possibly be a neighbor block of one of the 
!!   Paramesh blocks stored local to 'this' processor.
!!
!! AUTHORS
!!
!!    Kevin Olson, 2007, 2008.
!!
!!***

      Recursive Subroutine add_block_to_tree (t, tpar, top,            &
           coord,bsize,level,                                          &
           lb, proc, nodetype,                                         &
           which_child, n_children_at_0, no_blocks_at_0,               &
           off_proc, searched)

!-----Use Statements
      Use tree, Only : nchild, grid_xmin, grid_xmax,                   &
                       grid_ymin, grid_ymax, grid_zmin, grid_zmax
      Use paramesh_dimensions, Only : nxb, nyb, nzb, ndim, k2d, k3d
      Use physicaldata, Only : spherical_pm
      Use constants, Only : pi

      Implicit None

!-----Input/Output Arguments
      Type(node), Pointer :: t, tpar
      Type(node), Pointer :: top
      Type(node), Pointer :: tptr
      Real,    Intent(in) :: coord(3),bsize(3)
      Integer, Intent(in) :: lb, proc, level, nodetype
      Logical, Intent(in) :: off_proc
      Integer, Intent(in) :: which_child
      Integer, Intent(in) :: no_blocks_at_0
      Integer, Intent(inout) :: n_children_at_0
      Logical, Intent(inout) :: searched

!-----Local Variables
      Real :: xl, xr, yl, yr, zl, zr, spar(3), xpar, ypar, zpar
      Real :: x_other_part, y_other_part, z_other_part, m_other_part
      Real :: xp1, yp1, zp1, sizep(3), size_parent_cell
      Real :: xpp, ypp, zpp
      Integer :: i, ix, iy, iz
      Integer :: child_no, mype, ierr, curr_level, parent_level, ichild
      Logical :: is_close

!-----Include Statements
      Include 'mpif.h'

!-----Begin Executable code.
      Call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

!-----If the parent exists then store its size and location
      If (Associated(tpar)) Then 
         spar(:) = tpar%size(:)
         spar(2) = spar(2)*k2d
         spar(3) = spar(3)*k3d
         xpar = tpar%coord(1)
         ypar = tpar%coord(2)*k2d
         zpar = tpar%coord(3)*k3d
         parent_level = tpar%level
      Else
!-----If the parent does not exist then compute its size and location
         spar(1) = 2.*(grid_xmax-grid_xmin)
         spar(2) = 2.*(grid_ymax-grid_ymin)*k2d
         spar(3) = 2.*(grid_zmax-grid_zmin)*k3d
         xpar = grid_xmin + spar(1)/2.
         ypar = grid_ymin + spar(2)/2.
         zpar = grid_zmin + spar(3)/2.
         ypar = ypar*k2d
         zpar = zpar*k3d
         parent_level = -1
      End If

!-----If 'this' (sub)tree is empty then add a new subtree for the block to go
!-----to.
      If (.Not.associated(t)) Then

!--------If the searching block originated from off processor then
!--------search the local tree to see if is can possibly be a neighbor of
!--------any of the local blocks (this works because the local blocks are
!--------added first before any off-processor blocks).
         If (off_proc .and. .Not.searched) Then
            curr_level = level
            searched = .true.
            is_close = .false.
            Call search_and_prune_local_tree(top, coord(1),            &
                 coord(2), coord(3),                                   &
                 bsize(:),                                             &
                 curr_level,                                           &
                 is_close)
!Was:
!            If (spherical_pm .and. .Not.is_close) Then
! CEG added extra if
            If (spherical_pm) Then
              If (.Not.is_close) Then
               xpp = coord(1)
               ypp = coord(2)
               zpp = coord(3)
               If (ypp-bsize(2) < 0. .or. ypp+bsize(2) > pi) Then
                 If (zpp < pi) Then
                     zpp = zpp + pi
                  Elseif(zpp > pi)  Then
                     zpp = zpp - pi
                  End If
                  curr_level = level
                  Call search_and_prune_local_tree(top, xpp,           &
                     ypp, zpp,                                         &
                     bsize(:),                                         &
                     curr_level,                                       &
                     is_close)
               End If
              End If  ! End If (spherical_pm .and. .Not.is_close)
            End If  ! End If (spherical_pm .and. .Not.is_close)
         Else
            is_close = .true.
         End If  ! End If (off_proc .and. .Not.searched)

!--------If the searching block can possibly be a neighbor block to any of 
!--------of the local blocks, then create a new tree node
         If (is_close) Then  

            Allocate(t)
            t%local = .False.
            t%searched = .False.
            If (.Not.off_proc) Then
               t%local = .True.
            End If
            t%level = parent_level + 1 
            If (t%level == 0) Then
               t%nchild = no_blocks_at_0
            Else
               t%nchild = 2**ndim
            End If
            Allocate(t%child(t%nchild))
            Do i = 1,t%nchild
               nullify(t%child(i)%ptr)
            End Do
            If (associated(tpar)) Then
               t%parent%ptr => tpar
            End If

            If (t%level == 1) Then
               t%coord(:) = coord(:)
               t%size(:)  = bsize(:)
            Else
               t%size(:) = spar(:)/2.
               If (which_child == 1) Then
                  t%coord(1) = xpar - spar(1)/4.
                  t%coord(2) = ypar - spar(2)/4.
                  t%coord(3) = zpar - spar(3)/4.
               Else If (which_child == 2) Then
                  t%coord(1) = xpar + spar(1)/4.
                  t%coord(2) = ypar - spar(2)/4.
                  t%coord(3) = zpar - spar(3)/4.
               Else If (which_child == 3) Then
                  t%coord(1) = xpar - spar(1)/4.
                  t%coord(2) = ypar + spar(2)/4.
                  t%coord(3) = zpar - spar(3)/4.
               Else If (which_child == 4) Then
                  t%coord(1) = xpar + spar(1)/4.
                  t%coord(2) = ypar + spar(2)/4.
                  t%coord(3) = zpar - spar(3)/4.
               Else If (which_child == 5) Then
                  t%coord(1) = xpar - spar(1)/4.
                  t%coord(2) = ypar - spar(2)/4.
                  t%coord(3) = zpar + spar(3)/4.
               Else If (which_child == 6) Then
                  t%coord(1) = xpar + spar(1)/4.
                  t%coord(2) = ypar - spar(2)/4.
                  t%coord(3) = zpar + spar(3)/4.
               Else if (which_child == 7) Then
                  t%coord(1) = xpar - spar(1)/4.
                  t%coord(2) = ypar + spar(2)/4.
                  t%coord(3) = zpar + spar(3)/4.
               Else If (which_child == 8) Then
                  t%coord(1) = xpar + spar(1)/4.
                  t%coord(2) = ypar + spar(2)/4.
                  t%coord(3) = zpar + spar(3)/4.
               End If  ! End If (which_child == 1)
            End If  ! End If (t%level == 1)

            If (ndim < 2) Then
               t%coord(2) = 0.
               t%size(2)  = 0.
            End If
            If (ndim < 3) Then
               t%coord(3) = 0.
               t%size(3) = 0.
            End If

!--------If the level of the searching block is reached, store the lb, 
!--------proc information and exit
            If (t%level == level) Then
               t%lb = lb
               t%proc = proc
               t%nodetype = nodetype
               Return
            End If

!--------Search down newly created child node (do this until level of block
!--------in main tree is equal to the level in the local tree
            If (coord(1) <= t%coord(1)) Then
               ix = 0
            Else
               ix = 1
            End If
            If (coord(2) <= t%coord(2) .or. ndim < 2) Then
               iy = 0
            Else
               iy = 1
            End If
            If (coord(3) <= t%coord(3) .or. ndim < 3) Then
               iz = 0
            Else
               iz = 1
            End If
            child_no = iz*4 + iy*2 + ix + 1

            If (t%level == 0) Then
               n_children_at_0 = n_children_at_0+1
               child_no = n_children_at_0
            End If

            Call add_block_to_tree (t%child(child_no)%ptr, t,             &
              top,                                                     &
              coord,bsize,level,                                       &
              lb, proc, nodetype,                                      &
              child_no, n_children_at_0, no_blocks_at_0,               &
              off_proc, searched)

         End If  ! End If (is_close)

!-----Otherwise, if the s(sub)tree 't' exists, 
!-----check if block is in this subtree's box
      Else

!--------If the level is reached, store the lb, proc information and exit
         If (t%level == level) Then
            t%lb = lb
            t%proc = proc
            t%nodetype = nodetype
            Return
         End If

         If (coord(1) <= t%coord(1)) then
            ix = 0
         Else
            ix = 1
         End If
         If (coord(2) <= t%coord(2) .or. ndim < 2) then
            iy = 0
         Else
            iy = 1
         End If
         If (coord(3) <= t%coord(3) .or. ndim < 3) then
            iz = 0
         Else
            iz = 1
         End If
         child_no = iz*4 + iy*2 + ix + 1

         If (t%level == 0 .and. level == 1) then
            n_children_at_0 = n_children_at_0 + 1
         End If
         If (t%level == 0) Then
            Do ichild = 1, n_children_at_0
               tptr => t%child(ichild)%ptr
               If (associated(tptr)) then
                If (coord(3) >= (tptr%coord(3) - tptr%size(3)/2.) .and.&
                    coord(3) <= (tptr%coord(3) + tptr%size(3)/2.) .or. &
                    ndim < 3) Then 
                If (coord(2) >= (tptr%coord(2) - tptr%size(2)/2.) .and.&
                    coord(2) <= (tptr%coord(2) + tptr%size(2)/2.) .or. &
                    ndim < 2) Then
                If (coord(1) >= (tptr%coord(1) - tptr%size(1)/2.) .and.&
                    coord(1) <= (tptr%coord(1) + tptr%size(1)/2.)) then
                    child_no = ichild
                End If
                End If
                End If
               Else
                child_no = ichild
               End If  ! End If (associated(tptr))
            end Do  ! End Do ichild = 1, n_children_at_0
         End If  ! End If (t%level == 0) 

         Call add_block_to_tree (t%child(child_no)%ptr, t,             &
              top,                                                     &
              coord,bsize,level,                                       &
              lb, proc, nodetype,                                      &
              child_no, n_children_at_0, no_blocks_at_0,               &
              off_proc, searched)
         
      End If  ! End if (.not.associated(t))

      end subroutine add_block_to_tree

!!****f* mpi_source/search_sub_tree
!! NAME
!!
!!   search_sub_tree
!!
!! SYNOPSIS
!!
!!   Call search_sub_tree (t, neigh_coord,            
!!                         neigh_level,               
!!                         lb, proc, nnodetype, found)
!!
!!   Call search_sub_tree (type(node), real arrray,            
!!                         integer,               
!!                         integer, integer, integer, logical)
!! ARGUMENTS
!!
!!    t -> A pointer to the (sub)tree 'node' being tested and potentially
!!         which the block is added to.
!!    neigh_coord -> physical coordinates of a potential neighbor
!!    neigh_level -> level in tree where potential block is.
!!    lb -> block id of the neighboring block (if found).
!!    proc -> processor id of the neighboring block (if found).
!!    nnodetype -> nodetype of the neighboring block (if found).
!!    found -> a logical that indicates if the neighboring block has been
!!             found or not.
!!
!! INCLUDES
!!
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!
!! CALLS
!!
!!   self
!!
!! RETURNS
!!
!!   The block id and processor of a neighbor (if found).
!!
!! DESCRIPTION
!!
!!   Given the coordinates and level of a potential neighbor of a block,
!!   this routine searches the local tree and determines if a matching block
!!   exists.  If a match is found, the block id and processor where it can 
!!   be found are retruned. 
!!
!! AUTHORS
!!
!!    Kevin Olson, 2007, 2008.
!!
!!***

      Recursive Subroutine search_sub_tree (t, neigh_coord,            &
                                            neigh_level,               &
                                            lb, proc, nnodetype, found)

!-----Use Statements
      Use paramesh_dimensions, Only : ndim

      Implicit None

!-----Input/Output Arguments
      Type(node), Pointer :: t
      Real,    Intent(in) :: neigh_coord(3)
      Integer, Intent(in) :: neigh_level
      Integer, Intent(inout) :: lb, proc  
      Integer, Intent(out)   :: nnodetype
      Logical, Intent(inout) :: found

!-----Local Variables
      Integer :: i
      Real    :: dx, dy, dz

!-----Include Statements
      Include 'mpif.h'

!-----Begin executable code.

!-----If the tree node exists and a neighbor has not been found yet the
!-----check this tree node
      If (Associated(t) .and. .Not.found) Then

!--------If the tree node's level is <= the level of the potential block
!--------then check the node futher.
         If (t%level <= neigh_level) Then

!--------compute the distance between the tree node and the potential neighbor
            dx = abs(t%coord(1) - neigh_coord(1))
            If (ndim < 2) Then
               dy = 0.
            Else
               dy = abs(t%coord(2) - neigh_coord(2))
            End If
            If (ndim < 3) Then
               dz = 0.
            Else
               dz = abs(t%coord(3) - neigh_coord(3))
            End If

!--------If the distances are less than the size of the tree node, then
!--------examine the node further.
            If (dx <= t%size(1)/1.9 .and.                                 &
                dy <= t%size(2)/1.9 .and.                                 &
                dz <= t%size(3)/1.9) Then

!--------If the coordinates of the potential neighbor block are equal to
!--------to the coordintates of the tree node AND their levels are equal,
!--------the neighboring block is found !
               If (all(neigh_coord(1:ndim) <                                 &
                    t%coord(1:ndim)+(t%size(1:ndim)*0.01)) .and.          &
                    all(neigh_coord(1:ndim) >                                 &
                    t%coord(1:ndim)-(t%size(1:ndim)*0.01)) .and.          &
                    neigh_level         == t%level) Then
                  lb = t%lb
                  proc = t%proc
                  nnodetype = t%nodetype
                  found = .True.
!-----------The block is found, so exit.
                  Return
               Else
!-----------Otherwise search down each child subtree.
                  Do i = 1, t%nchild
                     Call search_sub_tree (t%child(i)%ptr,neigh_coord,       &
                                     neigh_level,lb,proc,nnodetype,    &
                                     found)
                  End Do
               End If  ! End If (all(neigh_coord(1:ndim) ...

            End If  ! End If (dx <= t%size(1) ...

         End If  ! End If (neigh_level <= t%level)

      End If  ! End If (associated(t) .and. .Not.found)

      End Subroutine search_sub_tree

!!****f* mpi_source/search_and_prune_local_tree
!! NAME
!!
!!   search_and_prune_local_tree
!!
!! SYNOPSIS
!!
!!    Call search_and_prune_local_tree(loc_t,          
!!                                     xp, yp, zp, sizep,  
!!                                     curr_level, is_close)
!!
!!    Call search_and_prune_local_tree(type(node),          
!!                                     real, real, real, real array,  
!!                                     integer, logical)
!! ARGUMENTS
!!
!!    loc_t -> A pointer to the (sub)tree 'node' being searched.
!!    xp, yp, zp -> coordinates of searching block
!!    sizep -> sizes of searching block
!!    curr_level -> current level in the tree at which the block is searching
!!    is_close -> a logical that indicates if the block should be added to
!!                the local tree or not
!!
!! INCLUDES
!!
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!
!! CALLS
!!
!!   self
!!
!! RETURNS
!!
!!   Returns a local 'is_close' which indicates if the block needs to be
!!   added to the local tree (.True.) or not (.False.).
!!
!! DESCRIPTION
!!
!!   Subroutine to 'search' the tree and determine if 'off-processor' blocks
!!   should trigger the creation of a new node in the local tree during
!!   the execution of the subroutine 'add_particle_to_tree'.
!!
!! AUTHORS
!!
!!    Kevin Olson, 2007, 2008.
!!
!!***

      Recursive Subroutine search_and_prune_local_tree(loc_t,          &
                                                  xp, yp, zp, sizep,   &
                                                  curr_level, is_close)

!-----Use Statements
      Use paramesh_dimensions, Only : ndim

      Implicit None

!-----Include Statements.
      Include 'mpif.h'

!-----Input/Output Arguments
      Type(node), Pointer :: loc_t
      Real, Intent(in) :: xp, yp, zp, sizep(3)
      Integer, Intent(in) :: curr_level
      Logical, Intent(inout) :: is_close

!-----Local Variables
      Real :: xl, yl, zl, xr, yr, zr
      Real :: dx, dy, dz, r, rmin, rcrit, s1, s2
      Integer :: i, j,  mype, ierr, ix, iy, iz, imin
      Logical :: search_it

!-----Begin Executable Code

!-----If the searching block has not been marked to be added, then search
!-----tree further
      If (.not.is_close) Then
!-----If the tree node exists, then search further
      If (Associated(loc_t)) Then
!-----If the tree node is local to the calling processor, then search further
      If (loc_t%local) Then   
!-----If the level of the tree node is <= the level of the searching block is at
      If (loc_t%level <= curr_level) Then

         dx = abs(loc_t%coord(1) - xp)
         If (ndim < 2) Then
            dy = 0.
         Else
            dy = abs(loc_t%coord(2) - yp)
         End If
         If (ndim < 3) Then
            dz = 0.
         Else
            dz = abs(loc_t%coord(3) - zp)
         End If

!--------If the searching block is close enough to be a neighbor to the tree 
!--------node, then search further.
         If (dx <= loc_t%size(1)*1.01 .and.                            &
             dy <= loc_t%size(2)*1.01 .and.                            &
             dz <= loc_t%size(3)*1.01) Then 

!--------If the level the block is at is eqaul to the level in the local tree,
!--------then stop searching, set is_close to .True. and exit
         If (curr_level == loc_t%level) Then
            is_close = .True.
            Return
         Else
!--------Otherwise, search down the children
            Do i = 1, loc_t%nchild
               Call search_and_prune_local_tree(                       &
                    loc_t%child(i)%ptr,                                &
                    xp, yp, zp, sizep, curr_level,                     &
                    is_close)
            End Do
         End If  ! End If (associated(t) .and. .Not.found)

         End If  ! End If (dx <= loc_t%size(1)*1.01 .and. ...

      End If  ! End If (loc_t%level <= curr_level)
      End If  ! End If (loc_t%local)
      End If  ! End If (Associated(loc_t)) 
      End If  ! End If (.not.is_close)

      End Subroutine search_and_prune_local_tree


!!****f* mpi_source/free_local_tree
!! NAME
!!
!!   free_local_tree
!!
!! SYNOPSIS
!!
!!    Call subroutine free_local_tree(t)
!!
!!    Call subroutine free_local_tree(type(node))
!!
!! ARGUMENTS
!!
!!    t -> A pointer to the (sub)tree 'node'
!!
!! INCLUDES
!!
!! USES
!!
!! CALLS
!!
!!   self
!!
!! RETURNS
!!
!!   An empty tree is returned
!!
!! DESCRIPTION
!!
!!   Subroutine to free all memory allocated by a tree 't'.
!!
!! AUTHORS
!!
!!    Kevin Olson, 2007, 2008.
!!
!!***

      Recursive Subroutine free_local_tree(t)

      Implicit None

!-----Input/Output Arguments
      Type(node), Pointer :: t
!-----Local Variables
      Integer :: i

      If (Associated(t)) Then
         Do i = 1,t%nchild
            Call free_local_tree(t%child(i)%ptr)
         End Do
         Deallocate(t%child)
         Deallocate(t)
         Nullify(t)
      End If

      End Subroutine free_local_tree

!!****f* mpi_source/zero_local_tree
!! NAME
!!
!!   zero_local_tree
!!
!! SYNOPSIS
!!
!!    Call subroutine zero_local_tree(t)
!!
!!    Call subroutine zero_local_tree(type(node))
!!
!! ARGUMENTS
!!
!!    t -> A pointer to the (sub)tree 'node'
!!
!! INCLUDES
!!
!! USES
!!
!! CALLS
!!
!!   self
!!
!! RETURNS
!!
!!   An zero'ed-out tree is returned
!!
!! DESCRIPTION
!!
!!   Subroutine to zero-out all variables allocated by a tree 't'.
!!
!! AUTHORS
!!
!!    Kevin Olson, 2007, 2008.
!!
!!***

      Recursive Subroutine zero_local_tree(t)

      Implicit None

!-----Input/Output Variables
      Type(node), Pointer :: t

!-----Local Variables
      Integer :: i

      If (Associated(t)) Then
         Do i = 1,t%nchild
            Call zero_local_tree(t%child(i)%ptr)
         End Do
         t%lb = 0
         t%proc = 0
         t%level = 0
         t%size(:) = 0.
         t%coord(:) = 0.
      End If

      End Subroutine zero_local_tree

      End Module local_tree_module

