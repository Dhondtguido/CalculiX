!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine e_corio(co,nk,konl,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &  ff,nelem,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,
     &  t0,t1,ithermal,vold,iperturb,nelemload,
     &  sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &  matname,mi,ncmat_,ttime,time,istep,iinc,nmethod,ielprop,prop)
!
!     computation of the element matrix and rhs for the element with
!     the topology in konl
!
!     ff: rhs without temperature and eigenstress contribution
!
!     nmethod=0: check for positive Jacobian
!     nmethod=1: stiffness matrix + right hand side
!     nmethod=2: stiffness matrix + mass matrix
!     nmethod=3: static stiffness + buckling stiffness
!     nmethod=4: stiffness matrix + mass matrix
!
      implicit none
!
      logical mass,stiffness,buckling,rhsi
!
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
!
      integer konl(20),ifaceq(8,6),nelemload(2,*),nk,nbody,nelem,
     &  mi(*),mattyp,ithermal(*),iperturb(*),nload,idist,i,j,k,l,i1,i2,
     &  j1,nmethod,k1,l1,ii,jj,ii1,jj1,id,ipointer,ig,m1,m2,m3,m4,kk,
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),null,ielprop(*),
     &  ntmat_,nope,nopes,norien,ihyper,iexpl,kode,imat,mint2d,
     &  mint3d,ifacet(7,4),nopev,iorien,istiff,ncmat_,
     &  ifacew(8,5),intscheme,n,ipointeri,ipointerj,istep,iinc,
     &  layer,kspt,jltyp,iflag,iperm(60),m,iscale,ne0,
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),prop(*),
     &  s(60,60),w(3,3),p1(3),p2(3),bodyf(3),bodyfx(3),ff(60),
     &  bf(3),q(3),shpj(4,20),elcon(0:ncmat_,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),xkl(3,3),eknlsign,reltime,
     &  alcon(0:6,ntmat_,*),alzero(*),orab(7,*),t0(*),t1(*),
     &  anisox(3,3,3,3),voldl(0:mi(2),20),vo(3,3),
     &  xl2(3,8),xsj2(3),shp2(7,8),vold(0:mi(2),*),xload(2,*),
     &  v(3,3,3,3),
     &  om,omx,e,un,al,um,xi,et,ze,tt,const,xsj,xsjj,sm(60,60),
     &  sti(6,mi(1),*),stx(6,mi(1),*),s11,s22,s33,s12,s13,s23,s11b,
     &  s22b,s33b,s12b,s13b,s23b,t0l,t1l,
     &  senergy,senergyb,rho,elas(21),
     &  sume,factorm,factore,alp,elconloc(ncmat_),eth(6),
     &  weight,coords(3),dmass,xl1(3,8),term
!
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),plconloc(802),dtime,ttime,time,
     &  sax(60,60),ffax(60),gs(8,4),a
!
      data iflag /3/
      data iperm /13,14,-15,16,17,-18,19,20,-21,22,23,-24,
     &            1,2,-3,4,5,-6,7,8,-9,10,11,-12,
     &            37,38,-39,40,41,-42,43,44,-45,46,47,-48,
     &            25,26,-27,28,29,-30,31,32,-33,34,35,-36,
     &            49,50,-51,52,53,-54,55,56,-57,58,59,-60/
!
      include "gauss.f"
!
      null=0
!
      imat=ielmat(1,nelem)
      amat=matname(imat)
!
c     Bernhardi start
      if(lakonl(1:5).eq.'C3D8I') then
         nope=11
      elseif(lakonl(4:4).eq.'2') then
c     Bernhardi end
         nope=20
      elseif(lakonl(4:4).eq.'8') then
         nope=8
      elseif(lakonl(4:5).eq.'10') then
         nope=10
      elseif(lakonl(4:4).eq.'4') then
         nope=4
      elseif(lakonl(4:5).eq.'15') then
         nope=15
      elseif(lakonl(4:4).eq.'6') then
         nope=6
      elseif(lakonl(1:2).eq.'ES') then
         read(lakonl(8:8),'(i1)') nope
         nope=nope+1
      endif
!
      if(lakonl(4:5).eq.'8R') then
         mint3d=1
      elseif(lakonl(4:7).eq.'20RB') then
         if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
            mint3d=50
         else
            call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &           null,xi,et,ze,weight)
         endif
      elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
         if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
     &        (lakonl(7:7).eq.'E')) then
            mint3d=4
         else
            mint3d=8
         endif
      elseif(lakonl(4:4).eq.'2') then
         mint3d=27
      elseif(lakonl(4:5).eq.'10') then
         mint3d=4
      elseif(lakonl(4:4).eq.'4') then
         mint3d=1
      elseif(lakonl(4:5).eq.'15') then
         mint3d=9
      elseif(lakonl(4:4).eq.'6') then
         mint3d=2
      else
         mint3d=0
      endif
!     
!     computation of the coordinates of the local nodes
!     
      do i=1,nope
         do j=1,3
            xl(j,i)=co(j,konl(i))
         enddo
      enddo
!     
!     initialisation of s
!
      do i=1,3*nope
        do j=1,3*nope
          s(i,j)=0.d0
        enddo
      enddo
!
!     computation of the matrix: loop over the Gauss points
!
      do kk=1,mint3d
         if(lakonl(4:5).eq.'8R') then
            xi=gauss3d1(1,kk)
            et=gauss3d1(2,kk)
            ze=gauss3d1(3,kk)
            weight=weight3d1(kk)
         elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
               xi=gauss3d13(1,kk)
               et=gauss3d13(2,kk)
               ze=gauss3d13(3,kk)
               weight=weight3d13(kk)
            else
               call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &              kk,xi,et,ze,weight)
            endif
         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &           then
            xi=gauss3d2(1,kk)
            et=gauss3d2(2,kk)
            ze=gauss3d2(3,kk)
            weight=weight3d2(kk)
         elseif(lakonl(4:4).eq.'2') then
            xi=gauss3d3(1,kk)
            et=gauss3d3(2,kk)
            ze=gauss3d3(3,kk)
            weight=weight3d3(kk)
         elseif(lakonl(4:5).eq.'10') then
            xi=gauss3d5(1,kk)
            et=gauss3d5(2,kk)
            ze=gauss3d5(3,kk)
            weight=weight3d5(kk)
         elseif(lakonl(4:4).eq.'4') then
            xi=gauss3d4(1,kk)
            et=gauss3d4(2,kk)
            ze=gauss3d4(3,kk)
            weight=weight3d4(kk)
         elseif(lakonl(4:5).eq.'15') then
            xi=gauss3d8(1,kk)
            et=gauss3d8(2,kk)
            ze=gauss3d8(3,kk)
            weight=weight3d8(kk)
         elseif(lakonl(4:4).eq.'6') then
            xi=gauss3d7(1,kk)
            et=gauss3d7(2,kk)
            ze=gauss3d7(3,kk)
            weight=weight3d7(kk)
         endif
!     
!     calculation of the shape functions and their derivatives
!     in the gauss point
!     
c     Bernhardi start
         if(lakonl(1:5).eq.'C3D8R') then
            call shape8hr(xl,xsj,shp,gs,a)
         elseif(lakonl(1:5).eq.'C3D8I') then
            call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
         else if(nope.eq.20) then
c     Bernhardi end
            if(lakonl(7:7).eq.'A') then
               call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
            elseif((lakonl(7:7).eq.'E').or.(lakonl(7:7).eq.'S')) then
               call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            endif
         elseif(nope.eq.8) then
            call shape8h(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.10) then
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.4) then
            call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.15) then
            call shape15w(xi,et,ze,xl,xsj,shp,iflag)
         else
            call shape6w(xi,et,ze,xl,xsj,shp,iflag)
         endif
!
!           check the jacobian determinant
!
         if(xsj.lt.1.d-20) then
            write(*,*) '*ERROR in e_c3d: nonpositive jacobian'
            write(*,*) '       determinant in element',nelem
            write(*,*)
            xsj=dabs(xsj)
            nmethod=0
         endif
!
!           material data and local stiffness matrix
!
         istiff=1
         call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &        imat,amat,iorien,coords,orab,ntmat_,elas,rho,
     &        nelem,ithermal,alzero,mattyp,t0l,t1l,
     &        ihyper,istiff,elconloc,eth,kode,plicon,
     &        nplicon,plkcon,nplkcon,npmat_,
     &        plconloc,mi(1),dtime,kk,
     &        xstiff,ncmat_)
!
!
!           initialisation for the body forces
!
         om=2.d0*rho*dsqrt(omx)*weight
!
!           incorporating the jacobian determinant in the shape
!           functions
!
         xsjj=dsqrt(xsj)
         do i1=1,nope
            shpj(1,i1)=shp(1,i1)*xsjj
            shpj(2,i1)=shp(2,i1)*xsjj
            shpj(3,i1)=shp(3,i1)*xsjj
            shpj(4,i1)=shp(4,i1)*xsj
         enddo
!
         jj1=1
         do jj=1,nope
!     
            ii1=1
            do ii=1,jj
!     
!     Coriolis matrix
!     
               dmass=
     &              om*shpj(4,ii)*shp(4,jj)
               s(ii1,jj1+1)=s(ii1,jj1+1)-p2(3)*dmass
               s(ii1,jj1+2)=s(ii1,jj1+2)+p2(2)*dmass
               s(ii1+1,jj1)=s(ii1+1,jj1)+p2(3)*dmass
               s(ii1+1,jj1+2)=s(ii1+1,jj1+2)-p2(1)*dmass
               s(ii1+2,jj1)=s(ii1+2,jj1)-p2(2)*dmass
               s(ii1+2,jj1+1)=s(ii1+2,jj1+1)+p2(1)*dmass
!     
               ii1=ii1+3
            enddo
            jj1=jj1+3
         enddo
      enddo
!     
!     for axially symmetric and plane stress/strain elements: 
!     complete s and sm
!     
      if((lakonl(6:7).eq.'RA').or.(lakonl(6:7).eq.'RS').or.
     &     (lakonl(6:7).eq.'RE')) then
         do i=1,60
            do j=i,60
               k=abs(iperm(i))
               l=abs(iperm(j))
               if(k.gt.l) then
                  m=k
                  k=l
                  l=m
               endif
               sax(i,j)=s(k,l)*iperm(i)*iperm(j)/(k*l)
            enddo
         enddo
         do i=1,60
            do j=i,60
               s(i,j)=s(i,j)+sax(i,j)
            enddo
         enddo
!     
      endif
!     
      return
      end
      
