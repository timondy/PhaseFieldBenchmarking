!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
!
!!****h* headers/tree
!!
!! NAME
!!
!!   physicaldata
!!
!! SYNOPSIS
!!
!!   Module tree
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!
!! DESCRIPTION
!!
!!   Fortran 90 Module which 'holds' the tree data which PARAMESH uses
!!   for a quad or oct-tree data structure, implemented on a parallel computer.
!!   This data includes information such which blocks are neighbors of other 
!!   blocks, which block is another block's parent, and which blocks are child 
!!   blocks of a particular block.
!!
!!   A convention is established for numbering the neighbors (or faces
!!   of a block. The first neighbor is at lower x coordinate, the 
!!   second at higher x, the third at lower y, fourth at higher y, fifth
!!   at lower z and the sixth at higher z.
!!
!!   The convention by which the children of a block are numbered is the
!!   same as the fortran array ordering, so that the first child is
!!   at lower x, y and z coordinate, the second child is at upper x
!!   but lower y and z, the third is at lower x, upper y and lower z,
!!   and so on.
!!
!!   When a block has a refined neighbor we will need to know which children
!!   of this neighbor are to provide guard cell information. The id's of the
!!   correct children are stored in kchild using the conventions described 
!!   above. For example, if we are working on the 3rd neighbor of the
!!   current block and it is at finer refinement level, then we must access
!!   the children designated by kchild(:,3), in this case children 1, 2, 5
!!   and 6.
!!
!!   The tree organizes a set of up to maxblocks_tr grids on each processor.
!!   All the grids are assumed to be cartesian with a uniform size. Each 
!!   grid has a level of refinement associated with it. The set of level 0
!!   grids cover the computational domain without overlap. Each grid
!!   can be the parent of 2**d offspring grids which completely cover 
!!   the sub-domain of their parents, where d is the physical dimension
!!   of the simulation. The linear resolution varies by a factor of 2 
!!   between successive levels of refinement. At no point do we allow the
!!   resolution to jump by more than one level of refinement.
!!
!!
!!
!! AUTHORS
!!
!!  Peter MacNeice and Kevin Olson
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
!-----------------------------------------------------------------
! tree Module



      Module tree

      Use paramesh_dimensions

      Private

      Public :: maxblocks_tr
      Public :: nchild, nfaces, mchild, mfaces, mdim,  mflags
      Public :: mboundaries,nboundaries

! Block limit to be used in manipulating tree when modifying the grid.
      Integer :: maxblocks_tr

! Number of children of a node
      Integer :: nchild

! Number of faces on a grid block
      Integer :: nfaces

! Parameters used to define array sizes
      Integer, Parameter :: mdim=3,mchild=2**mdim,mfaces=2*mdim

! Parameters used to declare the number of block marker flags needed
      Integer, Save      :: mflags

! Parameters used to declare the number of boundary regions where boundary
! conditions are to be applied.  Typically: 2*ndim
      Integer,Save       :: nboundaries

      Integer, Parameter :: mboundaries=100

      Public :: neigh,child,which_child
      Public :: parent,lrefine,lnblocks,new_lnblocks
      Public :: nodetype,empty,bflags,newchild,derefine,refine
      Public :: stay,work_block,coord,bsize,bnd_box
      Public :: grid_xmin,grid_xmax,grid_ymin,grid_ymax
      Public :: grid_zmin,grid_zmax
      Public :: lrefine_max,lrefine_min
      Public :: level_cell_sizes

! Variables for storing tree datastructure
      Integer, Allocatable, Save :: neigh(:,:,:)
      Integer, Allocatable, Save :: child(:,:,:)
      Integer, Allocatable, Save :: which_child(:)
      Integer, Allocatable, Save :: parent(:,:),lrefine(:)
      Integer, Save :: lnblocks,new_lnblocks
      Integer, Allocatable, Save :: nodetype(:)
      Integer, Allocatable, Save :: empty(:),bflags(:,:)
      Logical, Allocatable, Save :: newchild(:)
      Logical, Allocatable, Save :: derefine(:),refine(:)
      Logical, Allocatable, Save :: stay(:)
      Real, Allocatable, Save :: work_block(:)
      Real, Allocatable, Save :: coord(:,:)
      Real, Allocatable, Save :: bsize(:,:)
      Real, Allocatable, Save :: bnd_box(:,:,:)
      Real,Save :: grid_xmin,grid_xmax
      Real,Save :: grid_ymin,grid_ymax
      Real,Save :: grid_zmin,grid_zmax
      Real, Allocatable, Save :: level_cell_sizes(:,:)
      Integer, Save :: lrefine_max,lrefine_min

! flag to record grid change
      Public :: grid_changed, grid_analysed_mpi
      Integer, Save :: grid_changed, grid_analysed_mpi

! added for surrblks calculation with mpi
      Public :: boundary_box,boundary_index
      Real, Allocatable,Save    :: boundary_box(:,:,:)
      Integer, Allocatable, Save :: boundary_index(:)

! added for use with mpi block buffering
      Public :: strt_buffer,last_buffer
      Public :: strt_buffer_tree,last_buffer_tree
      Public :: laddress,surr_blks
      Public :: surr_morts
      Integer, Save :: strt_buffer,last_buffer
      Integer, Save :: strt_buffer_tree,last_buffer_tree
      Integer, Allocatable, Save :: surr_blks(:,:,:,:,:)
      Integer, Allocatable, Save :: surr_morts(:,:,:,:,:)
      Integer, Allocatable, Save :: laddress(:,:)

! arrays to store info about block neighbors which are boundaries
      Public :: bc_block_neighs,bc_block_neighs_send
      Public :: bc_block_neighs_length
      Public :: bc_block_neighs_status
      Integer,Save,Allocatable :: bc_block_neighs(:,:)
      Integer,Save,Allocatable :: bc_block_neighs_send(:,:)
      Integer,Save             :: bc_block_neighs_length
      Integer,Save             :: bc_block_neighs_status


! DECLARE variables which are Targets
      Target refine, derefine, newchild, empty
      Target lrefine, nodetype, work_block
      Target parent, coord, bsize, neigh
      Target child, bnd_box, stay

!!****v* tree/neigh
!!
!! NAME
!!
!!   neigh
!!
!! SYNOPSIS
!!
!!   Public, Integer :: neigh(2,mfaces,maxblocks_tr)
!!
!!   Public, Integer, Allocatable :: neigh(:,:,:)
!!
!! DESCRIPTION
!!
!!   Local and processor ids of a block's neighbors, at that block's same 
!!   refinement level.  If a neighbor does not exist both values are set to -1,
!!   unless that face is at an external domain boundary where non-periodic 
!!   boundary conditions are to be applied, in which case these are set to 
!!   -20 or less, depending on the boundary conditions to be applied on the 
!!   boundary in question.
!!
!!***

!!****v* tree/child
!!
!! NAME
!!
!!   child
!!
!! SYNOPSIS
!!
!!   Public, Integer :: child(2,mchild,maxblocks_tr)
!!
!!   Public, Integer, Allocatable :: child(:,:,:)
!!
!! DESCRIPTION
!!
!!   Local and processor ids of a block's children.
!!
!!***

!!****v* tree/parent
!!
!! NAME
!! 
!!   parent
!!
!! SYNOPSIS
!!
!!   Public, Integer :: parent(2,maxblocks_tr)
!!
!!   Public, Integer, Allocatable :: parent(:,:)
!!
!! DESCRIPTION
!!  
!!   Local and processor ids of a block's parent.
!!
!!***

!!****v* tree/coord
!!
!! NAME
!!
!!   coord
!!
!! SYNOSIS
!!
!!   Public, Real :: coord(mdim, maxblocks_tr)
!!
!!   Public, Real, Allocatable :: coord(mdim, maxblocks_tr)
!!
!! DESCRIPTION
!!
!!   An array storing x,y and z coordinates of the center of a block.
!!
!!***

!!****v* tree/bnd_box
!!
!! NAME
!!
!!   bnd_box
!!
!! SYNOPSIS
!!
!!   Public, Real :: bnd_box(2,mdim,maxblocks_tr)
!!
!!   Public, Real, Allocatable :: bnd_box(:,:,:)
!!
!! DESCRIPTION
!!
!!   The bounding box information for a block. The lower edge of block i along
!!   along the j-th coordinate axis is at bnd_box(1,j,i) and the upper edge
!!   is at bnd_box(2,j,i).
!!
!!***

!!****v* tree/bsize
!! NAME
!!
!!   bsize
!!
!! SYNOPSIS
!!
!!   Public, Real :: bsize(mdim,maxblocks_tr)
!!
!!   Public, Real, Allocatable :: bsize(:,:)
!!
!! DESCRIPTION
!!
!!   Physical size of a block in the x, y and z directions.
!!
!!***

!!****v* tree/lrefine
!!
!! NAME
!!
!!   lrefine
!!
!! SYNOPSIS
!!
!!   Public, Integer :: lrefine(maxblocks_tr)
!!
!!   Public, Integer, Allocatable :: lrefine(:)
!!
!! DESCRIPTION
!!
!!   The refinement level of a block.
!!
!!***

!!****v* tree/nodetype
!!
!! NAME
!!
!!   nodetype
!!
!! SYNOPSIS
!!
!!   Public, Integer :: nodetype(maxblocks_tr)
!!   
!!   Public, Integer, Allocatable :: nodetype(:)
!!
!! DESCRIPTION
!!
!!   Defines the node type, if 1 then the node is a leaf node, if 2 then the 
!!   node is a leaf node, if 2 then the node is a parent but with at least
!!   1 leaf child, otherwise it is set to 3 and it does not have any up-to-date
!!   data.
!!
!!***

!!****v* tree/empty
!!
!! NAME
!!
!!    empty
!!
!! SYNOPSIS
!!
!!   Public, Integer :: empty(maxblocks_tr)
!!
!!   Public, Integer, Allocatable :: empty(:)
!!
!! DESCRIPTION
!!
!!   Used to designate empty blocks, for example, when an obstacle is inserted 
!!   inside the computational domain. normal blocks have empty=0, empty blocks 
!!   have empty=1.
!!
!!***

!!****v* tree/bflags
!!
!! NAME
!!
!!   bflags
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: bflags(mflags,maxblocks_tr)
!!
!!   Public, Integer, Allocatable :: bflags(:,:)
!!
!! DESCRIPTION
!!
!!   An array of Integer flags which can be used to control computation on the 
!!   grid blocks and which are inherited by children from their parents.
!!
!!***

!!****v* tree/which_child
!!
!! NAME
!!
!!   which_child
!!
!! SYNOPSIS
!!
!!   Public, Integer :: which_child(maxblocks_tr)
!!
!!   Public, Integer, Allocatable :: which_child(:)
!!
!! DESCRIPTION
!!
!!   An Integer identifying which part of the parents volume this child 
!!   corresponds to.
!!
!!***

!!****v* tree/newchild
!!
!! NAME
!!
!! newchild
!!
!! SYNOPSIS
!!
!!   Public, Logical :: new_child(maxblocks_tr)
!!
!!   Public, Logical, Allocatable :: newchild(:)
!!
!! DESCRIPTION
!!   
!!   If true then child has just been produced by a refinement step, otherwise 
!!   false.
!!
!!***

!!****v* tree/lnblocks
!!
!! NAME
!!
!!   lnblocks
!!
!! SYNOPSIS
!!
!!   Public, Integer :: lnblocks
!!
!! DESCRIPTION
!!
!!   The number of blocks on the local processor.
!!
!!***

!!****v* tree/new_lnblocks
!!
!! NAME
!!  
!!   new_lnblocks
!!
!! SYNOPSIS
!! 
!!   Public, Integer :: new_lnblocks
!!
!! DESCRIPTION
!!
!!   The new number of blocks on the local processor after a refinement or 
!!   derefinement step.
!!
!!***

!!****v* tree/refine
!!
!! NAME
!!
!!   refine
!!
!! SYNOPSIS
!!
!!   Public, Logical :: refine(maxblocks_tr)
!!
!!   Public, Logical, Allocatable :: refine(:)
!!    
!! DESCRIPTION
!!
!!   The refinement flag. If set to .true. for a block, then that block will be
!!   refined during the next call to 'amr_refine_derefine'.
!!
!!***

!!****v* tree/derefine
!!
!! NAME
!!
!!   derefine
!!
!! SYNOPSIS
!!
!!   Public, Logical :: derefine(maxblocks_tr)
!!
!!   Public, Logical, Allocatable :: derefine(:)
!!    
!! DESCRIPTION
!!
!!  The derefinement flag. If set to .true. for a block, then PARAMESH will 
!!  attempt to derefine that block at the next call to 'amr_refine_derefine' if 
!!  it comforms to the rules for derefinement (i.e. all siblings marked for 
!!  derefinement and will not produce a refinement jump of more than 1 level).
!!
!!***

!!****v* tree/stay
!!
!! NAME
!!
!!   stay
!!
!! SYNOPSIS
!!
!!   Public, Logical :: stay(maxblocks_tr)
!!
!!   Public, Logical, Allocatable :: stay(:)
!!    
!! DESCRIPTION
!!
!!  If set to .true. for a block, then PARAMESH will neither refine or derefine
!!  the block at the next call to 'amr_refine_derefine'.
!!
!!***

!!****v* tree/work_block
!!
!! NAME
!!
!!   work_block
!!
!! SYNOPSIS
!!
!!   Public, Real :: work_block(maxblocks_tr)
!!
!!   Public, Real, Allocatable :: work_block(:)
!!    
!! DESCRIPTION
!!
!!    The work weighting given to each block from which the load balancing 
!!    across processors is calculated.
!!
!!***

!-----------------------------------------------------------------

      End Module tree
