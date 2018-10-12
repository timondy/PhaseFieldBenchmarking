!!
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation

       function minmod(a,b) result(mm)

        implicit none

        real, intent(in) :: a,b
        real :: mm
        if(abs(a)<abs(b))then
         mm=a
        else
         mm=b
        end if
       end function minmod

       subroutine prol_fc_dbz_init(n,i_divf_fc_vars)
        use prolong_arrays, only :  & 
     &   prol_fc_dbz, prol_fc_dbz_ivar, prol_fc_dbz_n
        use physicaldata, only : interp_mask_facex,  & 
     &   interp_mask_facey,interp_mask_facez

        implicit none

        integer, intent(in) :: n, i_divf_fc_vars(:,:)
        integer i,iface
        prol_fc_dbz_n = n ! n should equal nbndvar (or nfacevar?)
        If (.Not.allocated(prol_fc_dbz_ivar))                          &
            allocate(prol_fc_dbz_ivar(3,prol_fc_dbz_n))
        do i = 1,prol_fc_dbz_n
        do iface=1,3
         prol_fc_dbz_ivar(iface,i) = i_divf_fc_vars(iface,i)
        end do
        interp_mask_facex(i_divf_fc_vars(1,i)) = -200 ! Corresponds to fc_dbz.
        interp_mask_facey(i_divf_fc_vars(2,i)) = -200 ! Corresponds to fc_dbz.
        interp_mask_facez(i_divf_fc_vars(3,i)) = -200 ! Corresponds to fc_dbz.
        end do
        prol_fc_dbz = .true.
       end subroutine prol_fc_dbz_init

       function prol_fc_dbz_varp(ivar, iface) result(ldbz)
        use prolong_arrays, only :  & 
     &    prol_fc_dbz, prol_fc_dbz_ivar, prol_fc_dbz_n

        implicit none

        integer, intent(in) :: ivar, iface ! iface from {1 2 3} id'd {x y z}
        logical :: ldbz ! true if this var is prolonged with dbz routines
        integer :: i
        ldbz = .false.
        do i = 1, prol_fc_dbz_n
         ldbz = (prol_fc_dbz_ivar(iface,i) == ivar).or.ldbz
        end do
        return
       end function prol_fc_dbz_varp

       subroutine amr_1blk_fc_prol_dbz(  & 
     &        recvfx, recvfy, recvfz,        & 
     &        nfacevar_in, iv1, iv2, iv3,  & 
     &        ia,ib,ja,jb,ka,kb,    & 
     &        idest,ioff,joff,koff,          & 
     &        mype,lb,parent_pe,parent_blk   & 
     & )

        use paramesh_dimensions 
        use physicaldata
        use tree
        ! use prolong_arrays 

        implicit none

        real, intent(inout), dimension(:,:,:,:) :: recvfx,recvfy,recvfz
        integer, intent(in) :: nfacevar_in
        integer, intent(in) :: iv1,iv2,iv3
        integer, intent(in) :: ia,ib,ja,jb,ka,kb
        integer, intent(in) :: idest, ioff, joff, koff
        integer, intent(in) :: mype, lb, parent_pe, parent_blk

        integer :: nv,n1,n2,n3
        integer :: icl, icu, jcl, jcu, kcl, kcu
        integer :: i, j, k
        real :: x0, y0, z0, x1, y1, z1, x2, y2, z2
        real :: dx, dy, dz, ddx, ddy, ddz
        real :: fx_0, fy_0, fz_0, fx_1, fy_1, fz_1
        real :: tx2, ty2, tz2
        logical :: left, right

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!mlrdbg    
        integer :: ii,jj,kk
        integer :: i1, j1, k1
        logical :: mlrdbg,mlrdbg1,mlrdbg2,mlrdbg3
        logical :: mlrdbg100

        integer :: plb,plp
        real  :: tx1,ty1,tz1,tdx,tdy,tdz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        mlrdbg =.false. ! .true.
        mlrdbg1=.false. ! .true.
        mlrdbg2=.false. ! .true.
        mlrdbg3=.false. ! .true.
        mlrdbg100=.false. ! .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        nv= nfacevar
        n1=iu_bnd1-il_bnd1+1
        n2=ju_bnd1-jl_bnd1+k2d
        n3=ku_bnd1-kl_bnd1+k3d

        plb = parent_blk
        plp = parent_pe

        if(nfacevar.ne.nfacevar_in)then
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'a1fpd:nfacevar inconsistent!'
          call amr_abort
        end if

        !------------------------------------

        ! Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb

    ! Interpolation loop.
    !
    ! Note that the range of indeces used in the facevar plane differs
    ! depending on the value of iface_off. This assumes that the face values
    ! corresponding to index 0 (ie nguard faces to the left of the block
    ! boundary) are never needed, when iface_off=-1. 

    ! Iterate up to upper_bound-1 and flag prolong operator to calculate values
    ! for right hands of intervals.
    !

        kloop: do k=kcl,kcu+k3d
         jloop: do j=jcl,jcu+k2d
          iloop: do i=icl,icu+1

             ! 
             ! Endpoints for the box.  Points to be prolonged to are, e.g. [x0,(y0+y1)/2,(z0+z1)/2]
             ! which are on the faces of the box. Must be reconsidered for non-Cartesian grids.
             !
             ! Choose convenient (local) coordinates.
             ! 

             ii = i - nguard + 1
             i1 = ii/2 + nguard + ioff
             jj = j - nguard + 1
             j1 = jj/2 + nguard + joff
             if (ndim == 3) then
                kk = k - nguard + 1
                k1 = kk/2 + nguard + koff
             else
                kk = k
                k1 = k
             end if

             dx = 2.0*bsize(1,lb)/real(nxb)
             dy = 2.0*bsize(2,lb)/real(nyb)
             dz = 2.0*bsize(3,lb)/real(nzb)

             x0 = -0.5*2.0*bsize(1,lb)
             y0 = -0.5*2.0*bsize(2,lb)
             z0 = -0.5*2.0*bsize(3,lb)
             
             ! Left hand side of intervals, *off included for pedantism
             x1 = x0 + dx*real(i1-1-nguard)
             y1 = y0 + dy*real(j1-1-nguard)
             z1 = z0 + dz*real(k1-1-nguard)

             ! Right hand side of intervals, *off included for pedantism
             x2 = x1 + dx
             y2 = y1 + dy
             z2 = z1 + dz

             ! Distance between coarse points

             ddx = dx
             ddy = dy
             if (ndim == 3) then
                ddz = dz
             else
                ddz = 1.
             end if

             tdx = bsize(1,lb)/real(nxb)
             if (ndim >= 2) then
                tdy = bsize(2,lb)/real(nyb)
             else
                tdy = 0.
             end if
             if (ndim == 3) then
                tdz = bsize(3,lb)/real(nzb)
             else
                tdz = 0.
             end if

             ! Require blocks be cut in half.
             tx1 = -tdx*mod(real(i-nguard),2.0)
             ty1 = -tdy*mod(real(j-nguard*k2d),2.0)
             tz1 = -tdz*mod(real(k-nguard*k3d),2.0)

             tx2 = tx1 + tdx
             ty2 = ty1 + tdy
             tz2 = tz1 + tdz

             ddx = dx
             ddy = dy
             ddz = dz

             left=.true. ! always calculate left sides of box
!!!             right=((i.eq.icu).or.(j.eq.jcu).or.(k.eq.kcu)) ! on the last box, calculate the other side.

             ! compute interpolated values at points referred to above.
             call balpro(  & 
     &                   iv1,iv2,iv3,          & ! Which variables to use in each array.
     &                   i1,j1,k1,   & ! Coords of the cell
     &                   nv, n1, n2, n3,       & ! 
     &                   tx1,ty1,tz1,             & ! 
     &                   tx2,ty2,tz2,             & ! 
     &                   recvfx, recvfy, recvfz,      & ! The field arrays.
     &                   fx_0, fy_0, fz_0,     & ! The output field at p0 points (0-faces).
     &                   fx_1, fy_1, fz_1,     & ! The output field at p1 points (1-faces).
     &                   ddx, ddy, ddz,         & ! Size of cell in each of these directions.
     &                   left, .true.,       & ! Flags to return f*_0 and f*_1 respectively.
     &                   k2d, k3d, mype )

             if (j <= jcu .And. k <= kcu) then
              facevarx1(iv1,i,j,k,idest) = fx_0
             end if
             if (i <= icu .And. k <= kcu) then
              facevary1(iv2,i,j,k,idest) = fy_0
             end if
             if (i <= icu .And. j <= jcu) then
              facevarz1(iv1,i,j,k,idest) = fz_0
             end if

!             if(left)then
!              facevarx1(iv1,i,j,k,idest) = fx_0
!              facevary1(iv2,i,j,k,idest) = fy_0
!              facevarz1(iv3,i,j,k,idest) = fz_0
!             end if
       
!             if(i == icu)then
!              facevarx1(iv1,i+1,j,k,idest)   = fx_1
!             end if
!             if (j == jcu) then
!              facevary1(iv2,i,j+k2d,k,idest) = fy_1
!             end if
!             if (k == kcu) then
!              facevarz1(iv3,i,j,k+k3d,idest) = fz_1
!             end if
       
            enddo iloop
           enddo jloop
          enddo kloop

         end subroutine amr_1blk_fc_prol_dbz

         subroutine balpro(  & 
     & iv1,iv2,iv3,     & ! Which variables to use in each array.
     & i,j,k,           & ! Coords of the cell, shouldn't this be one symbol?
     & nv, n1, n2, n3,      & ! Gibbs & Heavi. popularized the boldface notation.
     & x0,y0,z0,        & ! One hundred years ago...
     & x1,y1,z1,        & ! 
     & fx, fy, fz,      & ! The field arrays.
     & fx_0, fy_0, fz_0,   & 
     & fx_1, fy_1, fz_1,   & 
     & ddx, ddy, ddz,    & ! Size of cell in each of these directions.
     & left, right,       & ! Side(s) of cell to calculate
     & k2d, k3d, mype )

         implicit none

         real, external :: minmod

         real :: x0,y0,z0,x1,y1,z1

         logical :: left, right

         integer :: hp
         integer :: iv1,iv2,iv3,i,j,k
         integer :: nv,n1,n2,n3
         real ::  & 
     &     fx(nv,n1+1,n2,n3),  & 
     &     fy(nv,n1,n2+1,n3),  & 
     &     fz(nv,n1,n2,n3+1)  
         real :: fx_0, fy_0, fz_0
         real :: fx_1, fy_1, fz_1

         real :: ddx, ddy, ddz
         real :: ddx2, ddy2, ddz2

         real :: a0,b0,c0
         real :: ax, ay, az, bx, by, bz, cx, cy, cz

         real :: axx, axy, axz, ayz
         real :: bxy, bxz, byy, byz
         real :: cxy, cxz, cyz, czz

         real :: axyz, bxyz, cxyz
         real :: axxy, axxz, byyx, byyz, cxzz, cyzz

         real :: fxp, fyp, fzp
         real :: fxm, fym, fzm

         real ::  & 
     & dy_fxp, dz_fxp, dx_fyp, dz_fyp, dx_fzp, dy_fzp,  & 
     & dy_fxm, dz_fxm, dx_fym, dz_fym, dx_fzm, dy_fzm 

         real :: dxy_fzp1, dxy_fzp2, dxy_fzp3

         real :: dxy_fzp, dxy_fzm, dyz_fxp,  & 
     &                             dyz_fxm, dxz_fyp, dxz_fym

         integer, intent(in) :: k2d,k3d
         integer :: ihm, jhm, khm, ihp, jhp, khp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!mlrdbg    
         integer :: mype
         logical :: mlrdbg,baldbg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         mlrdbg=.false. ! .true.
         baldbg=.false. ! .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The p-variables can refer to the greater end of the current cell, and the m-variables
    ! refer to the current cell, or the p-variables can refer to the current cell while the
    ! m-variables refer to the previous cell.
         !!! legacy hm=0; hp=1

         ihm= 0;    jhm= 0;      khm= 0
         !! ihp = hp;    jhp = hp;      khp = hp
         ihp= 1;    jhp= k2d;    khp= k3d  ! previously, these were all set to hp=1.

    !!!here!!!

         fxp = fx(iv1,i+ihp,j,k)
         fyp = fy(iv2,i,j+jhp,k)
         fzp = fz(iv3,i,j,k+khp)

         fxm = fx(iv1,i+ihm,j,k)
         fym = fy(iv2,i,j+jhm,k)
         fzm = fz(iv3,i,j,k+khm)

         dy_fxp = minmod(fx(iv1,i+ihp,j+1,k)-fx(iv1,i+ihp,j,k), & 
     &                   fx(iv1,i+ihp,j,k)-fx(iv1,i+ihp,j-1,k)) ! dy(i+ihp,[j,j-1],k)
         dz_fxp = minmod(fx(iv1,i+ihp,j,k+1)-fx(iv1,i+ihp,j,k), & 
     &                   fx(iv1,i+ihp,j,k)-fx(iv1,i+ihp,j,k-1)) ! dz(i+ihp,j,[k,k-1])

         dx_fyp = minmod(fy(iv2,i+1,j+jhp,k)-fy(iv2,i,j+jhp,k), & 
     &                   fy(iv2,i,j+jhp,k)-fy(iv2,i-1,j+jhp,k)) ! dx([i,i-1],j+jhp,k)
         dz_fyp = minmod(fy(iv2,i,j+jhp,k+1)-fy(iv2,i,j+jhp,k), & 
     &                   fy(iv2,i,j+jhp,k)-fy(iv2,i,j+jhp,k-1)) ! dz(i,j+jhp,[k,k-1])

         dx_fzp = minmod(fz(iv3,i+1,j,k+khp)-fz(iv3,i,j,k+khp), & 
     &                   fz(iv3,i,j,k+khp)-fz(iv3,i-1,j,k+khp)) ! dx([i,i-1],j,k+khp)
         dy_fzp = minmod(fz(iv3,i,j+1,k+khp)-fz(iv3,i,j,k+khp), & 
     &                   fz(iv3,i,j,k+khp)-fz(iv3,i,j-1,k+khp)) ! dy(i,[j,j-1],k+khp)

         dy_fxm = minmod(fx(iv1,i+ihm,j+1,k)-fx(iv1,i+ihm,j,k), & 
     &                   fx(iv1,i+ihm,j,k)-fx(iv1,i+ihm,j-1,k)) ! dy(i+ihm,[j,j-1],k)
         dz_fxm = minmod(fx(iv1,i+ihm,j,k+1)-fx(iv1,i+ihm,j,k), & 
     &                   fx(iv1,i+ihm,j,k)-fx(iv1,i+ihm,j,k-1)) ! dz(i+ihm,j,[k,k-1])

         dx_fym = minmod(fy(iv2,i+1,j+jhm,k)-fy(iv2,i,j+jhm,k), & 
     &                   fy(iv2,i,j+jhm,k)-fy(iv2,i-1,j+jhm,k)) ! dx([i,i-1],j+jhm,k)
         dz_fym = minmod(fy(iv2,i,j+jhm,k+1)-fy(iv2,i,j+jhm,k), & 
     &                   fy(iv2,i,j+jhm,k)-fy(iv2,i,j+jhm,k-1)) ! dz(i,j+jhm,[k,k-1])

         dx_fzm = minmod(fz(iv3,i+1,j,k+khm)-fz(iv3,i,j,k+khm), & 
     &                   fz(iv3,i,j,k+khm)-fz(iv3,i-1,j,k+khm)) ! dx([i,i-1],j,k+khm)
         dy_fzm = minmod(fz(iv3,i,j+1,k+khm)-fz(iv3,i,j,k+khm), & 
     &                   fz(iv3,i,j,k+khm)-fz(iv3,i,j-1,k+khm)) ! dy(i,[j,j-1],k+khm)

         ! Is the following correct?
         dxy_fzp =  & 
     & minmod(  & 
     &  fz(iv3,i+1,j+1,k+khp) + fz(iv3,i,j,k+khp) -  & 
     &  fz(iv3,i+1,j,k+khp) - fz(iv3,i,j+1,k+khp),  & 
     &  fz(iv3,i,j,k+khp) + fz(iv3,i-1,j-1,k+khp) -  & 
     &  fz(iv3,i,j-1,k+khp) - fz(iv3,i-1,j,k+khp)  & 
     &  )

         dxy_fzm =  & 
     & minmod(  & 
     &  fz(iv3,i+1,j+1,k+khm) + fz(iv3,i,j,k+khm) -  & 
     &  fz(iv3,i+1,j,k+khm) - fz(iv3,i,j+1,k+khm),  & 
     &  fz(iv3,i,j,k+khm) + fz(iv3,i-1,j-1,k+khm) -  & 
     &  fz(iv3,i,j-1,k+khm) - fz(iv3,i-1,j,k+khm)  & 
     &  )

         dyz_fxp =  & 
     & minmod(  & 
     &  fx(iv1,i+ihp,j+1,k+1) + fx(iv1,i+ihp,j,k) -  & 
     &  fx(iv1,i+ihp,j+1,k) - fx(iv1,i+ihp,j,k+1),  & 
     &  fx(iv1,i+ihp,j,k) + fx(iv1,i+ihp,j-1,k-1) -  & 
     &  fx(iv1,i+ihp,j,k-1) - fx(iv1,i+ihp,j-1,k)  & 
     & )

         dyz_fxm =  & 
     & minmod(  & 
     &  fx(iv1,i+ihm,j+1,k+1) + fx(iv1,i+ihm,j,k) -  & 
     &  fx(iv1,i+ihm,j+1,k) - fx(iv1,i+ihm,j,k+1),  & 
     &  fx(iv1,i+ihm,j,k) + fx(iv1,i+ihm,j-1,k-1) -  & 
     &  fx(iv1,i+ihm,j,k-1) - fx(iv1,i+ihm,j-1,k) & 
     & )

         dxz_fyp =  & 
     & minmod(  & 
     &  fy(iv2,i+1,j+jhp,k+1) + fy(iv2,i,j+jhp,k) -  & 
     &  fy(iv2,i+1,j+jhp,k) - fy(iv2,i,j+jhp,k+1), & 
     &  fy(iv2,i,j+jhp,k) + fy(iv2,i-1,j+jhp,k-1) -  & 
     &  fy(iv2,i,j+jhp,k-1) - fy(iv2,i-1,j+jhp,k)  & 
     & )

         dxz_fym =  & 
     & minmod(  & 
     &  fy(iv2,i+1,j+jhm,k+1) + fy(iv2,i,j+jhm,k) -  & 
     &  fy(iv2,i+1,j+jhm,k) - fy(iv2,i,j+jhm,k+1),  & 
     &  fy(iv2,i,j+jhm,k) + fy(iv2,i-1,j+jhm,k-1) -  & 
     &  fy(iv2,i,j+jhm,k-1) - fy(iv2,i-1,j+jhm,k) & 
     & )

        ! ddx := ddx(xp-xm)
        ! For our case, hm == 0, so ddx(...) := ddx(i1,j1,k1)

        ddx2=ddx*ddx
        ddy2=ddy*ddy
        ddz2=ddz*ddz

        axyz =       ( ( dyz_fxp - dyz_fxm ) / (ddy * ddz) ) / ddx !
        bxyz =       ( ( dxz_fyp - dxz_fym ) / (ddz * ddx) ) / ddy !
        cxyz =       ( ( dxy_fzp - dxy_fzm ) / (ddx * ddy) ) / ddz !

        ayz  = 0.5 * ( dyz_fxp + dyz_fxm ) / (ddy * ddz) !
        bxz  = 0.5 * ( dxz_fyp + dxz_fym ) / (ddz * ddx) !
        cxy  = 0.5 * ( dxy_fzp + dxy_fzm ) / (ddx * ddy) !

        ax   =       ( fxp - fxm ) / ddx   !
        by   =       ( fyp - fym ) / ddy   !
        cz   =       ( fzp - fzm ) / ddz   !

        axy  =       ( ( dy_fxp - dy_fxm ) / ddy ) / ddx  !
        byz  =       ( ( dz_fyp - dz_fym ) / ddz ) / ddy  !
        cxz  =       ( ( dx_fzp - dx_fzm ) / ddx ) / ddz  !

        axz  =       ( ( dz_fxp - dz_fxm ) / ddz ) / ddx  !
        bxy  =       ( ( dx_fyp - dx_fym ) / ddx ) / ddy  !
        cyz  =       ( ( dy_fzp - dy_fzm ) / ddy ) / ddz  !

        ay   = 0.5 * ( dy_fxp + dy_fxm ) / ddy + cxyz * ddx2 / 16.0 !
        bz   = 0.5 * ( dz_fyp + dz_fym ) / ddz + axyz * ddy2 / 16.0 !
        cx   = 0.5 * ( dx_fzp + dx_fzm ) / ddx + bxyz * ddz2 / 16.0 !

        az   = 0.5 * ( dz_fxp + dz_fxm ) / ddz + bxyz * ddx2 / 16.0 !
        bx   = 0.5 * ( dx_fyp + dx_fym ) / ddx + cxyz * ddy2 / 16.0 !
        cy   = 0.5 * ( dy_fzp + dy_fzm ) / ddy + axyz * ddz2 / 16.0 !

        axx  = - 0.5 * ( bxy + cxz )  !
        byy  = - 0.5 * ( cyz + axy )  !
        czz  = - 0.5 * ( axz + byz )  !

        a0   = 0.5 * ( fxp + fxm ) - 0.25 * axx * ddx2 !
        b0   = 0.5 * ( fyp + fym ) - 0.25 * byy * ddy2 !
        c0   = 0.5 * ( fzp + fzm ) - 0.25 * czz * ddz2 !

        cyzz = -axyz * 0.25 !
        byyz = cyzz !

        axxz = -bxyz * 0.25 !
        cxzz = axxz !

        byyx = -cxyz * 0.25 
        axxy = byyx

        if(left)then
          fx_0 = ( a0 + ax * x0 + axx * x0 * x0 )  & 
     &     + ( ay + axy * x0 + axxy * x0 * x0 ) * 0.5 * (y1+y0)  & 
     &     + ( az + axz * x0 + axxz * x0 * x0 ) * 0.5 * (z1 + z0)  & 
     &     + (ayz + axyz * x0) * 0.25 * (y1 + y0) * (z1 + z0)

          fy_0 = ( b0 + by * y0 + byy * y0 * y0 )  & 
     &     + ( bz + byz * y0 + byyz * y0 * y0) * 0.5 * (z1+z0)  & 
     &     + ( bx + bxy * y0 + byyx * y0 * y0) * 0.5 * (x1 + x0) &
     &     + (bxz + bxyz * y0) * 0.25 * (z1 + z0) * (x1 + x0)
       
          fz_0 = ( c0 + cz * z0 + czz * z0 * z0 )  & 
     &     + ( cx + cxz * z0 + cxzz * z0 * z0) * 0.5 * (x1+x0)  & 
     &     + ( cy + cyz * z0 + cyzz * z0 * z0) * 0.5 * (y1 + y0)  & 
     &     + (cxy + cxyz * z0) * 0.25 * (x1 + x0) * (y1 + y0)
         end if

    ! 
    ! These may have to be reconsidered for non-orthogonal cells.
    ! 
         if(right)then
           fx_1 = ( a0 + ax * x1 + axx * x1 * x1 )  & 
     &      + ( ay + axy * x1 + axxy * x1 * x1) * 0.5 * (y1+y0)  & 
     &      + ( az + axz * x1 + axxz * x1 * x1) * 0.5 * (z1 + z0)  & 
     &      + (ayz + axyz * x1) * 0.25 * (y1 + y0) * (z1 + z0)

           fy_1 = ( b0 + by * y1 + byy * y1 * y1 )  &
     &      + ( bz + byz * y1 + byyz * y1 * y1) * 0.5 * (z1+z0)  & 
     &      + ( bx + bxy * y1 + byyx * y1 * y1) * 0.5 * (x1 + x0)  &
     &      + (bxz + bxyz * y1) * 0.25 * (z1 + z0) * (x1 + x0)
       
           fz_1 = ( c0 + cz * z1 + czz * z1 * z1 )  &
     &      + ( cx + cxz * z1 + cxzz * z1 * z1) * 0.5 * (x1+x0)  & 
     &      + ( cy + cyz * z1 + cyzz * z1 * z1) * 0.5 * (y1 + y0)  & 
     &      + (cxy + cxyz * z1) * 0.25 * (x1 + x0) * (y1 + y0)
          end if
    
         end subroutine balpro
