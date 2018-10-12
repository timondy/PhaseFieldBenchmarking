

!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

       subroutine prol_fc_clean_divb_init(n,i_divf_fc_vars)
       use prolong_arrays, only :  & 
     & prol_fc_clean_divb, prol_fc_clean_divb_ivar, prol_fc_clean_divb_n
       implicit none
       integer, intent(in) :: n, i_divf_fc_vars(3,n)
       integer i,iface
       prol_fc_clean_divb_n = n ! n should equal nbndvar (or nfacevar?)
       allocate(prol_fc_clean_divb_ivar(3,prol_fc_clean_divb_n))
       do i = 1,prol_fc_clean_divb_n
        do iface=1,3
         prol_fc_clean_divb_ivar(iface,i) = i_divf_fc_vars(iface,i)
        end do
       end do
       prol_fc_clean_divb = .true.
       end subroutine prol_fc_clean_divb_init

       subroutine amr_1blk_fc_clean_divb(  & 
     &        nfacevar_in,  & 
     &        ia,ib,ja,jb,ka,kb,    & 
     &        ionea,ioneb,  & 
     &        jonea,joneb,  & 
     &        konea,koneb,  & 
     &        idest,ioff,joff,koff,          & 
     &        mype,lb,parent_pe,parent_blk   & 
     & )

       use prolong_arrays, only :  & 
     & prol_fc_clean_divb, prol_fc_clean_divb_ivar, prol_fc_clean_divb_n

       use paramesh_dimensions 
       use physicaldata
       use tree

       implicit none

       integer, intent(in) :: nfacevar_in
       integer, intent(in) :: ia,ib,ja,jb,ka,kb
       integer, intent(in) :: ionea,ioneb
       integer, intent(in) :: jonea,joneb
       integer, intent(in) :: konea,koneb
       integer, intent(in) :: idest, ioff, joff, koff
       integer, intent(in) :: mype, lb, parent_pe, parent_blk

       integer :: iv1,iv2,iv3,iprol
       integer :: nv,n1,n2,n3,i1,i2,i3,j1,k1,i,ii,j,jj,k,kk
       integer :: ik1,ik2,ik3
       integer :: iminx, imaxx, jminx, jmaxx, kminx, kmaxx
       integer :: iminy, imaxy, jminy, jmaxy, kminy, kmaxy
       integer :: iminz, imaxz, jminz, jmaxz, kminz, kmaxz

       integer :: icl,icu,jcl,jcu,kcl,kcu,ifo
       integer :: imin, imax
       integer :: jmin, jmax
       integer :: kmin, kmax
       
       real :: dx,dy,dz,x0,y0,z0

       logical :: status

       real, allocatable, dimension(:) ::  & 
     & x1, y1, z1, xc, yc, zc

       ifo = iface_off ! Do I need this?

       icl=ia; icu=ib; jcl=ja; jcu=jb; kcl=ka; kcu=kb

       nv= nfacevar
       n1=icu+ioneb+ifo-(icl+ionea)+1
       n2=jcu+joneb+ifo-(jcl+jonea)+1
       n3=kcu+koneb+ifo-(kcl+konea)+1

       allocate(xc(n1), yc(n2), zc(n3))
       allocate(x1(n1+1),y1(n2+1),z1(n3+1))

       dx = bsize(1,lb)/real(nxb)
       dy = bsize(2,lb)/real(nyb)
       dz = bsize(3,lb)/real(nzb)

       x0 = -0.5*bsize(1,lb)
       y0 = -0.5*bsize(2,lb)
       z0 = -0.5*bsize(3,lb)

       ik3=0
       kloop: do k=kcl+konea,kcu+koneb+ifo+1
        ik3=ik3+1
        k1=k
        z1(ik3) = z0 + dz*real(k1-1-nguard)
       end do kloop
       do ik3=1,n3
        zc(ik3)=0.5d0*(z1(ik3+1)+z1(ik3))
       end do

       ik2=0
       jloop: do j=jcl+jonea,jcu+joneb+ifo+1
        ik2=ik2+1
        j1=j
        y1(ik2) = y0 + dy*real(j1-1-nguard)
       end do jloop
       do ik2=1,n2
        yc(ik2)=0.5d0*(y1(ik2+1)+y1(ik2))
       end do

       ik1=0
       iloop: do i=icl+ionea,icu+ioneb+ifo+1
        ik1=ik1+1
        i1=i
        x1(ik1) = x0 + dx*real(i1-1-nguard)
       end do iloop
       do ik1=1,n1
        xc(ik1)=0.5d0*(x1(ik1+1)+x1(ik1))
       end do

        iminx = icl+ionea
        imaxx = icu+ioneb+ifo + 1
        jminx = jcl
        jmaxx = jcu
        kminx = kcl
        kmaxx = kcu

        iminy = icl
        imaxy = icu
        jminy = jcl+jonea
        jmaxy = jcu+joneb+ifo + k2d
        kminy = kcl
        kmaxy = kcu

        iminz = icl
        imaxz = icu
        jminz = jcl+jonea
        jmaxz = jcu+joneb+ifo
        kminz = kcl+konea
        kmaxz = kcu+koneb+ifo + k3d

       do iprol = 1, prol_fc_clean_divb_n
        iv1 = prol_fc_clean_divb_ivar(1,iprol)
        iv2 = prol_fc_clean_divb_ivar(2,iprol)
        iv3 = prol_fc_clean_divb_ivar(3,iprol)

        call clean_field(                                          &  
        facevarx1(iv1,iminx:imaxx,jminx:jmaxx,kminx:kmaxx,idest),  & 
        facevary1(iv2,iminy:imaxy,jminy:jmaxy,kminy:kmaxy,idest),  &
        facevarz1(iv3,iminz:imaxz,jminz:jmaxz,kminz:kmaxz,idest),  &
        iminx,imaxx,jminx,jmaxx,kminx,kmaxx,   &
        iminy,imaxy,jminy,jmaxy,kminy,kmaxy,   &
        iminz,imaxz,jminz,jmaxz,kminz,kmaxz,   &
        ndim, nxb, nyb, nzb, bsize(:,lb)  )       
       end do

        deallocate(x1,y1,z1)
        deallocate(xc,yc,zc)

       end subroutine amr_1blk_fc_clean_divb
