!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****h* headers/physicaldata
!!
!! NAME
!!
!!   physicaldata
!! 
!! SYNOPSIS
!!
!!   module physicaldata
!!  
!! USES
!!
!!   paramesh_dimensions
!!
!! DESCRIPTION
!!
!!   Fortran 90 Module which 'holds' the main solution data structures of
!!   PARAMESH.  These include 'unk', 'facevarx', 'facevary', 'facevarz', 
!!   'unk_e_x', 'unk_e_y', 'unk_e_z', and 'unk_n'.  This module also holds the
!!   data used in by the flux conservation routines (i.e. flux_x, flux_y, and 
!!   flux_z).  Other data which controls how the solution data arrays are 
!!   manipulated by PARAMESH are also kept in this module such as the control 
!!   arrays for selecting the algorithms used for interpolation on guardcells 
!!   filling (i.e. interp_mask_*) and for selective guardcell filling (i.e. 
!!   gcell_on_*).  Arrays storing geometry information (i.e. cell volumes, 
!!   cell-face areas etc.) are also kept here.
!!
!! AUTHORS
!!
!!  Peter MacNeice and Kevin Olson
!!
!!***


!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      module extvbles

      use paramesh_dimensions

      private

      public :: tmpCG
      DOUBLE PRECISION, allocatable,save ::  tmpCG(:)
      target :: tmpCG

      end module extvbles

!------------------------------------------------------------------------------
! Physicaldata module
!------------------------------------------------------------------------------

      module physicaldata

      use extvbles

      use paramesh_dimensions

      private

!----------------------
! Solution Variables
!----------------------

!---------------------------------------
! Allocate memory for solution variables
!---------------------------------------

!------------------------------------------------------------------------------
! CELL CENTERED CENTERED DATA
!------------------------------------------------------------------------------


! the solution for cell-centered quantities.
      Public :: unk, interp_mask_unk, interp_mask_unk_res, tmpCEG
! CEG
      public :: previt
      Public :: gcell_on_cc,int_gcell_on_cc
      Public :: ngcell_on_cc
      Public :: checkp_on_cc
      Public :: gcell_on_cc_pointer
!
!CEG tried removing "save" to get working on Hector
      Real,Allocatable,Save ::  unk(:,:,:,:,:)
      Real,Allocatable,Save ::  tmpCEG(:)
!      Real,Allocatable,Save ::  tmpCEG(:,:,:,:,:)

!      Real,Allocatable ::  unk(:,:,:,:,:)
      real,allocatable,save ::  previt(:,:,:,:,:)
      Integer,Allocatable,Save :: interp_mask_unk(:)
      Integer,Allocatable,Save :: interp_mask_unk_res(:)
      Integer,Allocatable,Save :: gcell_on_cc_pointer(:)
      Logical,Allocatable,Save :: gcell_on_cc(:)
      Logical,Allocatable,Save :: int_gcell_on_cc(:)
      Logical,Allocatable,Save :: checkp_on_cc(:)
      Integer, Save :: ngcell_on_cc
      target :: unk
      target :: previt


!!****v* physicaldata/unk
!!
!! NAME
!!
!!   unk
!! 
!! SYNOPSIS
!!
!!   Public, Real, target :: 
!!    unk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocks)
!!
!!   Public, Real, Allocatable, target :: unk(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Array for holding cell centered data.  May or may not have permanent 
!!   guardcell storage allocated. 
!!
!!***

!!****v* physicaldata/interp_mask_unk
!!
!! NAME
!!
!!   interp_mask_unk
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_unk(nvar)
!!
!!   Public, Integer, Allocatable :: interp_mask_unk(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the cell 
!!   centered data stored in unk.  The values stored in interp_mask_unk may 
!!   be changed at anytime during execution and they can take on different 
!!   values for each different variable stored in unk.
!!
!!   The default value(s) for interp_mask_unk is 1 for all variables, 
!!   1 through nvar.
!!   The default behaviour for using different values of unk is:
!!     interp_mask_unk =  0  -> direct injection
!!     interp_mask_unk =  1  -> linear interpolation
!!     interp_mask_unk =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_unk =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_unk =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_unk >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in unk, but to use a 2nd order polynomial interpolation 
!!     for variable 4 then set,
!!
!!     interp_mask_unk(4) = 2
!!
!!***

!!****v* physicaldata/interp_mask_unk_res
!!
!! NAME
!!
!!   interp_mask_unk_res
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_unk_res(nvar)
!!
!!   Public, Integer, Allocatable :: interp_mask_unk_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the cell 
!!   centered data stored in unk.  The values stored in interp_mask_unk_res 
!!   may be changed at anytime during execution and they can take on different 
!!   values for each different variable stored in unk.
!!
!!   The default value(s) for interp_mask_unk_res is 1 for all variables, 
!!     1 through nvar.
!!   The default behaviour for using different values of unk is:
!!     interp_mask_unk_res =  0  -> illegal value
!!     interp_mask_unk_res =  1  -> Simple average of child data to parent 
!!                                  (linear interpolation)
!!     interp_mask_unk_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_unk_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_unk_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_unk_res >= 20 -> user defined
!!   
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in unk, but to use a 2nd order polynomial interpolation 
!!     for variable 4 then set,
!!
!!     interp_mask_unk_res(4) = 2
!!
!!***

!!****v* physicaldata/gcell_on_cc
!!
!! NAME
!!
!!   gcell_on_cc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: gcell_on_cc(nvar)
!!
!!   Public, Logical, Allocatable :: gcell_on_cc(:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the unk array have 
!!   their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_cc(ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in unk will be filled, if it is set to 
!!   .false. then variable 'ivar' will be skipped and its guardcells will not 
!!   be filled.  
!!
!!   The values in gcell_on_cc can be changed at any point during a run and 
!!   allow a user to avoid unecessary data movement in cases when only a small
!!   set of the variables in unk need guardcell values at some point during 
!!   their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_cc
!!
!! NAME
!!
!!   int_gcell_on_cc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: int_gcell_on_cc(nvar)
!!
!!   Public, Logical, Allocatable :: int_gcell_on_cc(:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_cc
!!
!! NAME
!!
!!   ngcell_on_cc
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: ngcell_on_cc
!!
!! DESCRIPTION
!!
!!   An Integer variable which holds the current number of variables in unk
!!   selected for guardcell filling by the mask array 'gcell_on_cc'.
!!
!!***

!!****iv physicaldata/gcell_on_cc_pointer
!!
!! NAME
!!
!!   gcell_on_cc_pointer
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: gcell_on_cc_pointer(nvar)
!!
!!   Public, Integer, Allocatable :: gcell_on_cc_pointer(:)
!!
!! DESCRIPTION
!!
!!   I don't know ????
!!
!!***

!!***v* physicaldata/checkp_on_cc
!!
!! NAME
!!
!!   checkp_on_cc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: checkp_on_cc(nvar)
!!
!!   Public, Logical, Allocatable :: checkp_on_cc(:)
!!
!! DESCRIPTION
!!
!!   A Logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_cc(ivar) is set to to true, then the ivar'th variable in the 
!!   unk array is added to the checkpoint file.  The default behaviour has ALL
!!   the values in checkp_on_cc set to .true. so that all variables in unk are
!!   written to the checkpoint file.
!!
!!***

!------------------------------------------------------------------------------
! FACE CENTERED DATA
!------------------------------------------------------------------------------


! the solution for cell-face-centered quantities.
      Public :: facevarx,facevary,facevarz
      Public :: interp_mask_facex,interp_mask_facey,interp_mask_facez
      Public :: interp_mask_facex_res,interp_mask_facey_res, & 
     &          interp_mask_facez_res
      Public :: gcell_on_fc,int_gcell_on_fc
      Public :: ngcell_on_fc
      Public :: gcell_on_fc_pointer
      Public :: checkp_on_fc
      Real,Allocatable,Save ::  facevarx(:,:,:,:,:)
      Real,Allocatable,Save ::  facevary(:,:,:,:,:)
      Real,Allocatable,Save ::  facevarz(:,:,:,:,:)
      Integer,Allocatable,Save :: interp_mask_facex(:)
      Integer,Allocatable,Save :: interp_mask_facey(:)
      Integer,Allocatable,Save :: interp_mask_facez(:)
      Integer,Allocatable,Save :: interp_mask_facex_res(:)
      Integer,Allocatable,Save :: interp_mask_facey_res(:)
      Integer,Allocatable,Save :: interp_mask_facez_res(:)
      Integer,Allocatable,Save :: gcell_on_fc_pointer(:,:)
      Logical,Allocatable,Save :: gcell_on_fc(:,:)
      Logical,Allocatable,Save :: int_gcell_on_fc(:,:)
      Logical,Allocatable,Save :: checkp_on_fc(:,:)
      Integer, Save :: ngcell_on_fc(3)
      target :: facevarx,facevary,facevarz


!!****v* physicaldata/facevars
!!
!! NAME
!!
!!   facevarx, facevary, facevarz
!! 
!! SYNOPSIS
!!
!!   Public,Real ::  facevarx(nbndvar,                                   
!!                            il_bnd:iu_bnd+1,jl_bnd:ju_bnd,             
!!                            kl_bnd:ku_bnd,                             
!!                            maxblocksf)
!!   Public,Real ::  facevary(nbndvar,                                   
!!                            il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           
!!                            kl_bnd:ku_bnd,                             
!!                            maxblocksf)
!!   Public,Real ::  facevarz(nbndvar,                                   
!!                            il_bnd:iu_bnd,jl_bnd:ju_bnd,               
!!                            kl_bnd:ku_bnd+k3d,                         
!!                            maxblocksf)
!!
!!   Public, Real, Allocatable :: facevarx(:,:,:,:,:)
!!   Public, Real, Allocatable :: facevary(:,:,:,:,:)
!!   Public, Real, Allocatable :: facevarz(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Arrays for holding face centered data.  May or may not have permanent 
!!   guardcell storage allocated. 
!!
!!***

!!****v* physicaldata/interp_mask_facex(y,z)
!!
!! NAME
!!
!!   interp_mask_facex, interp_mask_facey, interp_mask_facez
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_facex(nbndvar)
!!   Public, Integer :: interp_mask_facey(nbndvar)
!!   Public, Integer :: interp_mask_facez(nbndvar)
!!
!!   Public, Integer, Allocatable :: interp_mask_facex(:)
!!   Public, Integer, Allocatable :: interp_mask_facey(:)
!!   Public, Integer, Allocatable :: interp_mask_facez(:)
!!
!! DESCRIPTION
!!
!!   Arrays to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the face
!!   centered data stored in facevarx(y,z).  The values stored in 
!!   interp_mask_facex(y,z) may be changed at anytime during execution and 
!!   they can take on different values for each different variable stored in 
!!   facevarx(y,z).
!!
!!   The default value(s) for interp_mask_facex(y,z) is 1 for all variables, 
!!   1 through nbndvar.
!!   The default behaviour for using different values of facevarx(y,z) is:
!!     interp_mask_facex(y,z) =  0  -> direct injection
!!     interp_mask_facex(y,z) =  1  -> linear interpolation
!!     interp_mask_facex(y,z) =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_facex(y,z) =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_facex(y,z) =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_facex(y,z) >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in facevarx(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in facevarx then set,
!!
!!     interp_mask_facex(4) = 2
!!
!!***

!!****v* physicaldata/interp_mask_facex(y,z)_res
!!
!! NAME
!!
!!   interp_mask_facex_res, interp_mask_facey_res, interp_mask_facez_res
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_facex_res(nbndvar)
!!   Public, Integer :: interp_mask_facey_res(nbndvar)
!!   Public, Integer :: interp_mask_facey_res(nbndvar)
!!
!!   Public, Integer, Allocatable :: interp_mask_facex_res(:)
!!   Public, Integer, Allocatable :: interp_mask_facey_res(:)
!!   Public, Integer, Allocatable :: interp_mask_facez_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the face
!!   centered data stored in facevarx(y,z).  The values stored in 
!!   interp_mask_facex(y,z)_res may be changed at anytime during execution and 
!!   they can take on different values for each different variable stored in 
!!   facevarx(y,z).
!!
!!   The default value(s) for interp_mask_facex(y,z)_res is 1 for all 
!!   variables, 1 through nbndvar.
!!   The default behaviour for using different values of facevarx(y,z) is:
!!     interp_mask_facex(y,z)_res =  0  -> illegal value
!!     interp_mask_facex(y,z)_res =  1  -> Simple average of child data to 
!!                                         parent (linear interpolation)
!!     interp_mask_facex(y,z)_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_facex(y,z)_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_facex(y,z)_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_facex(y,z)_res >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in facevarx(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in facevarx then set,
!!
!!     interp_mask_facex(4) = 2
!!
!!***

!!****v* physicaldata/gcell_on_fc
!!
!! NAME
!!
!!   gcell_on_fc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: gcell_on_fc(3,nbndvar)
!!   Public, Logical, Allocatable :: gcell_on_fc(:,:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the facevarx(y,z) 
!!   array have their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_fc(:,ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in facevarx(y,z) will be filled, if it is 
!!   set to .false. then variable 'ivar' will be skipped and its guardcells 
!!   will not be filled.  
!!
!!   The values in gcell_on_fc can be changed at any point during a run and 
!!   allow a user to avoid unecessary data movement in cases when only a small
!!   set of the variables in facevarx(y,z) need guardcell values at some point
!!   during their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_fc
!!
!! NAME
!!
!!   int_gcell_on_fc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: int_gcell_on_fc(3,nbndvar)
!!   Public, Logical, Allocatable :: int_gcell_on_fc(:,:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_fc
!!
!! NAME
!!
!!   ngcell_on_fc
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: ngcell_on_fc(3)
!!
!! DESCRIPTION
!!
!!   An Integer variable which holds the current number of variables in 
!!   facevarx(y,z) selected for guardcell filling by the mask array 
!!   'gcell_on_fc'.
!!
!!***

!!****iv physicaldata/gcell_on_fc_pointer
!!
!! NAME
!!
!!   gcell_on_fc_pointer
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: gcell_on_fc_pointer(3,nbndvar)
!!   Public, Integer, Allocatable :: gcell_on_fc_pointer(:,:)
!!
!! DESCRIPTION
!!
!!   I don't know ???
!!
!!
!!***

!!***v* physicaldata/checkp_on_fc
!!
!! NAME
!!
!!   checkp_on_fc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: checkp_on_fc(3,nbndvar)
!!   Public, Logical, Allocatable :: checkp_on_fc(:,:)
!!
!! DESCRIPTION
!!
!!   A Logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_fc(:,ivar) is set to true, then the ivar'th variable in the
!!   facevarx(y,z) array is added to the checkpoint file.  The default 
!!   behaviour has ALL the values in checkp_on_fc set to .true. so that all 
!!   variables in facevarx(y,z) are written to the checkpoint file.
!!
!!***

!----------------------------------------------------------------------------
! EDGE CENTERED DATA
!----------------------------------------------------------------------------

! the solution for cell-edge-centered quantities.
      Public :: unk_e_x,unk_e_y,unk_e_z,interp_mask_ec, & 
     &     interp_mask_ec_res
      Public :: gcell_on_ec,int_gcell_on_ec
      Public :: ngcell_on_ec
      Public :: gcell_on_ec_pointer
      Public :: checkp_on_ec
      Real,Allocatable,Save ::  unk_e_x(:,:,:,:,:)
      Real,Allocatable,Save ::  unk_e_y(:,:,:,:,:)
      Real,Allocatable,Save ::  unk_e_z(:,:,:,:,:)
      Integer,Allocatable,Save :: interp_mask_ec(:)
      Integer,Allocatable,Save :: interp_mask_ec_res(:)
      Integer,Allocatable,Save :: gcell_on_ec_pointer(:,:)
      Logical,Allocatable,Save :: gcell_on_ec(:,:)
      Logical,Allocatable,Save :: int_gcell_on_ec(:,:)
      Logical,Allocatable,Save :: checkp_on_ec(:,:)
      Integer, Save :: ngcell_on_ec(3)
      target :: unk_e_x,unk_e_y,unk_e_z

!!****v* physicaldata/unk_e_x(y,z)
!!
!! NAME
!!
!!   unk_e_x, unk_e_y, unk_e_z
!! 
!! SYNOPSIS
!!
!!   Public,Real ::  unk_e_x(nbndvare,                                   
!!                           il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,             
!!                           kl_bnd:ku_bnd+k3d,                             
!!                           maxblocksf)
!!   Public,Real ::  unk_e_y(nbndvare,                                   
!!                            il_bnd:iu_bnd+1,jl_bnd:ju_bnd,           
!!                            kl_bnd:ku_bnd+k3d,                             
!!                            maxblocksf)
!!   Public,Real ::  unk_e_z(nbndvare,                                   
!!                           il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,               
!!                           kl_bnd:ku_bnd,                         
!!                           maxblocksf)
!!
!!   Public, Real, Allocatable :: unk_e_x(:,:,:,:,:)
!!   Public, Real, Allocatable :: unk_e_y(:,:,:,:,:)
!!   Public, Real, Allocatable :: unk_e_z(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Arrays for holding edge centered data.  May or may not have permanent 
!!   guardcell storage allocated. 
!!
!!***

!!****v* physicaldata/interp_mask_ec(y,z)
!!
!! NAME
!!
!!   interp_mask_ec
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_ec(nbndvare)
!!
!!   Public, Integer, Allocatable :: interp_mask_ec(:)
!!
!! DESCRIPTION
!!
!!   Arrays to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the edge
!!   centered data stored in unk_e_x, y and z.  The values stored in 
!!   interp_mask_ec may be changed at anytime during execution and they 
!!   can take on different values for each different variable stored in 
!!   unk_e_x, y and z.
!!
!!   The default value(s) for interp_mask_ec is 1 for all variables, 
!!   1 through nbndvare.
!!   The default behaviour for using different values of unk_e_x(y,z) is:
!!     interp_mask_ec =  0  -> direct injection
!!     interp_mask_ec =  1  -> linear interpolation
!!     interp_mask_ec =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_ec =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_ec =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_ec >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in unk_e_x(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_e_x, y, and z then set,
!!
!!     interp_mask_ec(4) = 2
!!
!! NOTE:  Using interp_mask_ec controls all 3 arrays unk_e_x, unk_e_y, and
!!        unk_e_z and you do not have the freedom to control the interpolation
!!        order for these arrays individually as you do for the face centered
!!        variables.
!!
!!***

!!****v* physicaldata/interp_mask_ec_res
!!
!! NAME
!!
!!   interp_mask_ec_res
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_ec_res(nbndvare)
!!
!!   Public, Integer, Allocatable :: interp_mask_ec_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the edge
!!   centered data stored in unk_e_x, y, and z.  The values stored in 
!!   interp_mask_ec_res may be changed at anytime during execution and 
!!   they can take on different values for each different variable stored in 
!!   unk_e_x, y, and z.
!!
!!   The default value(s) for interp_mask_ec_res is 1 for all variables, 
!!     1 through nbndvare.
!!   The default behaviour for using different values of unk_e_x(y,z) is:
!!     interp_mask_ec_res =  0  -> illegal value
!!     interp_mask_ec_res =  1  -> Simple average of child data to parent 
!!                                  (linear interpolation)
!!     interp_mask_ec_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_ec_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_ec_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_ec_res >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in unk_e_x(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_e_x, y, and z then set,
!!
!!     interp_mask_ec_res(4) = 2
!!
!! NOTE:  Using interp_mask_ec_res controls all 3 arrays unk_e_x, unk_e_y, and
!!        unk_e_z and you do not have the freedom to control the interpolation
!!        order for these arrays individually as you do for the face centered
!!        variables.
!!
!!***

!!****v* physicaldata/gcell_on_ec
!!
!! NAME
!!
!!   gcell_on_ec
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: gcell_on_ec(3,nbndvare)
!!   Public, Logical, Allocatable :: gcell_on_ec(:,:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the unk_e_x(y,z) array
!!   have their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_fc(:,ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in unk_e_x(y,z) will be filled, if it 
!!   is set to .false. then variable 'ivar' will be skipped and its guardcells
!!   will not be filled.  
!!
!!   The values in gcell_on_ec can be changed at any point during a run 
!!   and allow a user to avoid unecessary data movement in cases when only 
!!   a small set of the variables in unk_e_x(y,z) need guardcell values at 
!!   some point during their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_ec
!!
!! NAME
!!
!!   int_gcell_on_ec
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: int_gcell_on_ec(3,nbndvare)
!!   Public, Logical, Allocatable :: int_gcell_on_ec(:,:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_ec
!!
!! NAME
!!
!!   ngcell_on_ec
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: ngcell_on_ec(3)
!!
!! DESCRIPTION
!!
!!   An Integer variable which holds the current number of variables in 
!!   unk_e_x(y,z) selected for guardcell filling by the mask array 
!!   'gcell_on_ec'.
!!
!!***

!!****iv physicaldata/gcell_on_ec_pointer
!!
!! NAME
!!
!!   gcell_on_ec_pointer
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: gcell_on_ec_pointer(3,nbndvare)
!!   Public, Integer, Allocatable :: gcell_on_ec_pointer(:,:)
!!
!! DESCRIPTION
!!
!!   I don't know ???
!!
!!
!!***

!!***v* physicaldata/checkp_on_ec
!!
!! NAME
!!
!!   checkp_on_ec
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: checkp_on_ec(3,nbndvare)
!!   Public, Logical, Allocatable :: checkp_on_ec(:,:)
!!
!! DESCRIPTION
!!
!!   A Logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_ec(:,ivar) is set to to true, then the ivar'th variable in 
!!   the unk_e_x(y,z) array is added to the checkpoint file.  The default 
!!   behaviour has ALL the values in checkp_on_ec set to .true. so that all 
!!   variables in unk_e_x(y,z) are written to the checkpoint file.
!!
!!***

!-----------------------------------------------------------------------------
! CELL CORNER DATA
!-----------------------------------------------------------------------------

! the solution for cell-corner based quantities.
      Public :: unk_n, interp_mask_nc, interp_mask_nc_res
      Public :: gcell_on_nc,int_gcell_on_nc
      Public :: ngcell_on_nc
      Public :: gcell_on_nc_pointer
      Public :: checkp_on_nc
      Real,Allocatable,Save ::  unk_n(:,:,:,:,:)
      Integer,Allocatable,Save :: interp_mask_nc(:)
      Integer,Allocatable,Save :: interp_mask_nc_res(:)
      Integer,Allocatable,Save :: gcell_on_nc_pointer(:)
      Logical,Allocatable,Save :: gcell_on_nc(:)
      Logical,Allocatable,Save :: int_gcell_on_nc(:)
      Logical,Allocatable,Save :: checkp_on_nc(:)
      Integer, Save :: ngcell_on_nc
      target :: unk_n

!!****v* physicaldata/unk_n
!!
!! NAME
!!
!!   unk_n
!! 
!! SYNOPSIS
!!
!!   Public,Real ::  unk_n(nbndvarc,                                   
!!                         il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,             
!!                         kl_bnd:ku_bnd+k3d,                             
!!                         maxblocksf)
!!
!!   Public, Real, Allocatable :: unk_n(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Arrays for holding cell corner data.  May or may not have permanent 
!!   guardcell storage allocated. 
!!
!!***

!!****v* physicaldata/interp_mask_nc
!!
!! NAME
!!
!!   interp_mask_nc
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_nc(nbndvarc)
!!
!!   Public, Integer, Allocatable :: interp_mask_nc(:)
!!
!! DESCRIPTION
!!
!!   Arrays to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the cell
!!   corner data stored in unk_n.  The values stored in interp_mask_nc may be 
!!   changed at anytime during execution and they can take on different values
!!   for each different variable stored in unk_n.
!!
!!   The default value(s) for interp_mask_nc is 1 for all variables, 
!!   1 through nbndvarc.
!!   The default behaviour for using different values of unk_n is:
!!     interp_mask_nc =  0  -> direct injection
!!     interp_mask_nc =  1  -> linear interpolation
!!     interp_mask_nc =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_nc =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_nc =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_nc >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in unk_n, but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_n then set,
!!
!!     interp_mask_nc(4) = 2
!!
!!***

!!****v* physicaldata/interp_mask_nc_res
!!
!! NAME
!!
!!   interp_mask_nc_res
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: interp_mask_nc_res(nbndvarc)
!!
!!   Public, Integer, Allocatable :: interp_mask_nc_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the cell
!!   corner data stored in unk_n.  The values stored in interp_mask_nc_res may
!!   be changed at anytime during execution and they can take on different 
!!   values for each different variable stored in unk_n.
!!
!!   The default value(s) for interp_mask_nc_res is 1 for all variables, 
!!     1 through nbndvarc.
!!   The default behaviour for using different values of unk_n is:
!!     interp_mask_nc_res =  0  -> illegal value
!!     interp_mask_nc_res =  1  -> Simple average of child data to parent 
!!                                 (linear interpolation)
!!     interp_mask_nc_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_nc_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_nc_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_nc_res >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in unk_n, but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_n then set,
!!
!!     interp_mask_nc_res(4) = 2
!!
!!
!!***

!!****v* physicaldata/gcell_on_nc
!!
!! NAME
!!
!!   gcell_on_nc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: gcell_on_nc(nbndvarc)
!!   Public, Logical, Allocatable :: gcell_on_nc(:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the unk_e_n array have 
!!   their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_nc(ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in unk_n will be filled, if it is set to 
!!   .false. then variable 'ivar' will be skipped and its guardcells will not 
!!   be filled.  
!!
!!   The values in gcell_on_nc can be changed at any point during a run and 
!!   allow a user to avoid unecessary data movement in cases when only a small
!!   set of the variables in unk_n need guardcell values at some point during 
!!   their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_nc
!!
!! NAME
!!
!!   int_gcell_on_nc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: int_gcell_on_nc(nbndvarc)
!!   Public, Logical, Allocatable :: int_gcell_on_nc(:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_nc
!!
!! NAME
!!
!!   ngcell_on_nc
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: ngcell_on_nc(3)
!!
!! DESCRIPTION
!!
!!   An Integer variable which holds the current number of variables in unk_n
!!   selected for guardcell filling by the mask array 'gcell_on_nc'.
!!
!!***

!!****iv physicaldata/gcell_on_nc_pointer
!!
!! NAME
!!
!!   gcell_on_nc_pointer
!! 
!! SYNOPSIS
!!
!!   Public, Integer :: gcell_on_nc_pointer(3,nbndvarc)
!!   Public, Integer, Allocatable :: gcell_on_nc_pointer(:)
!!
!! DESCRIPTION
!!
!!   I don't know ???
!!
!!
!!***

!!***v* physicaldata/checkp_on_nc
!!
!! NAME
!!
!!   checkp_on_nc
!! 
!! SYNOPSIS
!!
!!   Public, Logical :: checkp_on_nc(nbndvarc)
!!   Public, Logical, Allocatable :: checkp_on_nc(:)
!!
!! DESCRIPTION
!!
!!   A Logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_nc(ivar) is set to true, then the ivar'th variable in the unk_n
!!   array is added to the checkpoint file.  The default behaviour has ALL the
!!   values in checkp_on_nc set to .true. so that all variables in unk_e_n are
!!   written to the checkpoint file.
!!
!!***

!-----------------
! Timestep control
!-----------------

! arrays used for timestep control
      Public :: time_loc,dtlevel,phase_dt,loc_cycle,ncyc_local
      Public :: ldtcomplete
      Real, Save    :: dtlevel(maxlevels)
      Integer, Save :: phase_dt(maxlevels),loc_cycle(maxlevels)
      Integer, Save :: ncyc_local(maxlevels)
      Logical, Allocatable,Save :: ldtcomplete(:)
      Real, Allocatable,Save    :: time_loc(:)

      Public :: t_unk,tfacevarx,tfacevary,tfacevarz
      Public :: t_unk_e_x,t_unk_e_y,t_unk_e_z,t_unk_n
      Real, Allocatable,Save  :: t_unk(:,:,:,:,:)
      Real, Allocatable,Save  :: tfacevarx(:,:,:,:,:)
      Real, Allocatable,Save  :: tfacevary(:,:,:,:,:)
      Real, Allocatable,Save  :: tfacevarz(:,:,:,:,:)
      Real, Allocatable,Save  :: t_unk_e_x(:,:,:,:,:)
      Real, Allocatable,Save  :: t_unk_e_y(:,:,:,:,:)
      Real, Allocatable,Save  :: t_unk_e_z(:,:,:,:,:)
      Real, Allocatable,Save  :: t_unk_n(:,:,:,:,:)

!-----------------------------------------------------------------
! include header file defining 1blk data structure

!------------------------------------------------------------------------------
! data_1blk
!------------------------------------------------------------------------------
!
! This file declares the storage space used to handle the `current
! working block' when the user decides not to reserve permanent
! storage space for guardcells for all blocks, but instead to 
! fill guardcells as needed. This strategy requires 2 working blocks,
! one for the leaf node and one for its parent.

      Public :: unk1,facevarx1,facevary1,facevarz1
      Public :: unk_e_x1,unk_e_y1,unk_e_z1
      Public :: unk_n1
! the solution for cell-centered quantities.
      Real, Save, Allocatable :: unk1(:,:,:,:,:)
! the solution for cell-face-centered quantities.
      Real, Save, Allocatable :: facevarx1(:,:,:,:,:)
      Real, Save, Allocatable :: facevary1(:,:,:,:,:)
      Real, Save, Allocatable :: facevarz1(:,:,:,:,:)
! the solution for cell-edge-centered quantities.
      Real, Save, Allocatable ::  unk_e_x1(:,:,:,:,:)
      Real, Save, Allocatable ::  unk_e_y1(:,:,:,:,:)
      Real, Save, Allocatable ::  unk_e_z1(:,:,:,:,:)
! the solution for cell-corner based quantities.
      Real, Save, Allocatable ::  unk_n1(:,:,:,:,:)
      target :: unk1
      target :: facevarx1, facevary1, facevarz1
      target :: unk_e_x1, unk_e_y1, unk_e_z1
      target :: unk_n1

!------------------------------------------------------------------------------
! temporary copy of solution to be used when storing solution prior
! to use of amr_1blk_guardcell
!------------------------------------------------------------------------------

      Public :: gt_unk,gt_facevarx,gt_facevary,gt_facevarz
      Public :: gt_unk_e_x,gt_unk_e_y,gt_unk_e_z,gt_unk_n
      Real, Allocatable, Save :: gt_unk(:,:,:,:,:)
      Real, Allocatable, Save :: gt_facevarx(:,:,:,:,:)
      Real, Allocatable, Save :: gt_facevary(:,:,:,:,:)
      Real, Allocatable, Save :: gt_facevarz(:,:,:,:,:)
      Real, Allocatable, Save :: gt_unk_e_x(:,:,:,:,:)
      Real, Allocatable, Save :: gt_unk_e_y(:,:,:,:,:)
      Real, Allocatable, Save :: gt_unk_e_z(:,:,:,:,:)
      Real, Allocatable, Save :: gt_unk_n(:,:,:,:,:)
!
! Variables used to control data caching. This helps to avoid unnecessary
! repetition of some guardcell and surrounding-block mapping operations
      Public ::  pcache_blk_u,pcache_pe_u,pcache_blk_w,pcache_pe_w
      Public ::  lnew_parent
      Integer, Save :: pcache_blk_u,pcache_pe_u
      Integer, Save :: pcache_blk_w,pcache_pe_w
      Logical, Save :: lnew_parent

!-----------------------------------------------------------------
! include header file defining data structure on cell faces

!------------------------------------------------------------------------
! block_boundary_data


!
! This file defines a data structure to be used for quantities
! which may need to be defined at grid block interfaces, eg fluxes,
! pressures.
!

! The convention for relating varaibles associated with cell faces to the
! variables defined at cell centers is as follows:

! If iface_off=0 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i-1/2) face of the (i,j,k)-th mesh cell. 
! If iface_off=-1 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i+1/2) face of the (i,j,k)-th mesh cell. 


! storage used for fluxes at block boundaries. This is used when conservation
! constraints need to be imposed.

      Public :: nfluxvar,nfluxes,maxblocksfl
      Public :: flux_x,flux_y,flux_z
      Public :: tflux_x,tflux_y,tflux_z
      Integer,Save :: nfluxvar
      Integer,Save :: nfluxes
      Integer,Save :: maxblocksfl
      Real, Allocatable, Save ::  flux_x(:,:,:,:,:)
      Real, Allocatable, Save ::  flux_y(:,:,:,:,:)
      Real, Allocatable, Save ::  flux_z(:,:,:,:,:)
      Real, Allocatable, Save ::  tflux_x(:,:,:,:,:)
      Real, Allocatable, Save ::  tflux_y(:,:,:,:,:)
      Real, Allocatable, Save ::  tflux_z(:,:,:,:,:)
      target :: flux_x, flux_y, flux_z 


! temporary flux storage needed inside amr_flux_conserve when using
! variable timestep
      Public :: ttflux_x,ttflux_y,ttflux_z
      Real, Allocatable, Save :: ttflux_x(:,:,:,:,:)
      Real, Allocatable, Save :: ttflux_y(:,:,:,:,:)
      Real, Allocatable, Save :: ttflux_z(:,:,:,:,:)

! storage used for cell edges at block boundaries. 
! This is used when quantities located at cell edge centers need to
! be used consistently at the boundaries between blocks at different
! refinement levels.

      Public :: nedgevar1,nedgevar,nedges,maxblockse
      Public :: bedge_facex_y,bedge_facex_z,bedge_facey_x
      Public :: bedge_facey_z,bedge_facez_x,bedge_facez_y
      Public :: recvarx1e,recvary1e,recvarz1e
      Public :: recvarx2e,recvary2e,recvarz2e
      Integer, Save :: nedgevar1
      Integer, Save :: nedgevar
      Integer, Save :: nedges
      Integer, Save :: maxblockse
      Real, Allocatable, Save ::  bedge_facex_y(:,:,:,:,:)
      Real, Allocatable, Save ::  bedge_facex_z(:,:,:,:,:)
      Real, Allocatable, Save ::  bedge_facey_x(:,:,:,:,:)
      Real, Allocatable, Save ::  bedge_facey_z(:,:,:,:,:)
      Real, Allocatable, Save ::  bedge_facez_x(:,:,:,:,:)
      Real, Allocatable, Save ::  bedge_facez_y(:,:,:,:,:)

      Real, Allocatable, Save ::  recvarx1e(:,:,:,:)
      Real, Allocatable, Save ::  recvary1e(:,:,:,:)
      Real, Allocatable, Save ::  recvarz1e(:,:,:,:)
      Real, Allocatable, Save ::  recvarx2e(:,:,:,:)
      Real, Allocatable, Save ::  recvary2e(:,:,:,:)
      Real, Allocatable, Save ::  recvarz2e(:,:,:,:)

      Public :: tbedge_facex_y,tbedge_facex_z,tbedge_facey_x
      Public :: tbedge_facey_z,tbedge_facez_x,tbedge_facez_y      
      Real, Allocatable, Save ::  tbedge_facex_y(:,:,:,:,:)
      Real, Allocatable, Save ::  tbedge_facex_z(:,:,:,:,:)
      Real, Allocatable, Save ::  tbedge_facey_x(:,:,:,:,:)
      Real, Allocatable, Save ::  tbedge_facey_z(:,:,:,:,:)
      Real, Allocatable, Save ::  tbedge_facez_x(:,:,:,:,:)
      Real, Allocatable, Save ::  tbedge_facez_y(:,:,:,:,:)
      Public :: ttbedge_facex_y,ttbedge_facex_z,ttbedge_facey_x
      Public :: ttbedge_facey_z,ttbedge_facez_x,ttbedge_facez_y      
      Real, Allocatable, Save ::  ttbedge_facex_y(:,:,:,:,:)
      Real, Allocatable, Save ::  ttbedge_facex_z(:,:,:,:,:)
      Real, Allocatable, Save ::  ttbedge_facey_x(:,:,:,:,:)
      Real, Allocatable, Save ::  ttbedge_facey_z(:,:,:,:,:)
      Real, Allocatable, Save ::  ttbedge_facez_x(:,:,:,:,:)
      Real, Allocatable, Save ::  ttbedge_facez_y(:,:,:,:,:)




! arrays used to store geometry information for the working block
      Public :: cell_vol
      Public :: cell_area1,cell_area2,cell_area3
      Public :: cell_leng1,cell_leng2,cell_leng3
      Public :: cell_face_coord1,cell_face_coord2,cell_face_coord3
      Real, Allocatable  :: cell_vol(:,:,:)
      Real, Allocatable  :: cell_area1(:,:,:)
      Real, Allocatable  :: cell_area2(:,:,:)
      Real, Allocatable  :: cell_area3(:,:,:)
      Real, Allocatable  :: cell_leng1(:,:,:)
      Real, Allocatable  :: cell_leng2(:,:,:)
      Real, Allocatable  :: cell_leng3(:,:,:)
      Real, Allocatable  :: cell_face_coord1(:)
      Real, Allocatable  :: cell_face_coord2(:)
      Real, Allocatable  :: cell_face_coord3(:)

! workspace arrays used for inter-block communications
      Public :: nbndmax
      Public :: recvarxf,recvaryf,recvarzf
      Public :: bndtempx1,bndtempy1,bndtempz1
      Integer :: nbndmax
      Real, Allocatable, Save :: recvarxf(:,:,:,:)
      Real, Allocatable, Save :: recvaryf(:,:,:,:)
      Real, Allocatable, Save :: recvarzf(:,:,:,:)
      Real, Allocatable, Save :: bndtempx1(:,:,:,:)
      Real, Allocatable, Save :: bndtempy1(:,:,:,:)
      Real, Allocatable, Save :: bndtempz1(:,:,:,:)

! parameters used in communication calls
      Public :: len_block_bndx,len_block_bndy,len_block_bndz
      Public :: len_block_ex,len_block_ey,len_block_ez
      Integer :: len_block_bndx
      Integer :: len_block_bndy
      Integer :: len_block_bndz
      Integer :: len_block_ex
      Integer :: len_block_ey
      Integer :: len_block_ez

!-----------------------------------------------------------------
! Array used to store variables which make up any divergence
! free fields
      Integer,Allocatable,Save,Public :: i_divf_fc_vars(:,:)

!-----------------------------------------------------------------
! Index arrays used in boundary condition routines.
      Integer,Public :: bc_index_i(2,3,5)
      Integer,Public :: bc_index_j(2,3,5)
      Integer,Public :: bc_index_k(2,3,5)


!-----------------------------------------------------------------
! Logical flags required to signal algorithmic states
      Logical, Public :: lrestrict_in_progress
      Logical, Public :: lprolong_in_progress
      Logical, Public :: lguard_in_progress

!-----------------------------------------------------------------
! Logical flag to indicate if restrictionless guardcell filling
! has been selected
      Logical,Public :: l_f_to_c

! Index arrays used to record destination data values for fine layer
! neighbor guardcells
      Integer,Public :: f2c_ind_unk(2,3,27)
      Integer,Public :: f2c_ind_facex(2,3,27)
      Integer,Public :: f2c_ind_facey(2,3,27)
      Integer,Public :: f2c_ind_facez(2,3,27)
      Integer,Public :: f2c_ind_unkex(2,3,27)
      Integer,Public :: f2c_ind_unkey(2,3,27)
      Integer,Public :: f2c_ind_unkez(2,3,27)
      Integer,Public :: f2c_ind_unkn(2,3,27)

!-----------------------------------------------------------------
! Mpi communication pattern identifier
      Integer, Public :: mpi_pattern_id

!-----------------------------------------------------------------
! Error trapping and management

! To record whether amr_gsurrounding_blks has been called.
      Public :: gsurrblks_set
      Integer :: gsurrblks_set

! a counter which can be used to keep track of calls to routines
      Public :: instance
      Integer :: instance

! Diagonals flag
      Public :: diagonals
      Logical, Save :: diagonals

! error checking flag
      Public :: amr_error_checking
      Logical, Save :: amr_error_checking

! no_permanent_guardcells flag
      Public :: no_permanent_guardcells
      Logical, Save :: no_permanent_guardcells

! advance_all_levels flag
      Public :: advance_all_levels
      Logical, Save :: advance_all_levels

! force_consistency flag
      Public :: force_consistency
      Logical, Save :: force_consistency

! consv_fluxes and consv_flux_densities flag
      Public :: consv_fluxes, consv_flux_densities
      Logical, Save :: consv_fluxes, consv_flux_densities

! edge_value and edge_value_integ flag
      Public :: edge_value, edge_value_integ
      Logical, Save :: edge_value, edge_value_integ

! var_dt flag
      Public :: var_dt
      Logical, Save :: var_dt

! pred_corr flag
      Public :: pred_corr
      Logical, Save :: pred_corr

! empty_cells flag
      Public :: empty_cells
      Logical, Save :: empty_cells

! conserve flag
      Public :: conserve
      Logical, Save :: conserve

! divergence_free flag
      Public :: divergence_free
      Logical, Save :: divergence_free

! curvilinear
      Public :: curvilinear, cartesian_pm, cylindrical_pm,  & 
     &          spherical_pm, polar_pm, lsingular_line,     & 
     &          curvilinear_conserve
      Logical, Save :: curvilinear, cartesian_pm, cylindrical_pm,  & 
     &                 spherical_pm, polar_pm, lsingular_line,     & 
     &                 curvilinear_conserve
!-----------------------------------------------------------------

! CEG version of multigrid block list start addresses

      public :: block_starts
      integer, allocatable :: block_starts(:)

      end module physicaldata
