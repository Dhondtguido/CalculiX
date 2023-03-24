!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine e_c3d_u(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &     ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &     alzero,ielmat,ielorien,norien,orab,ntmat_,
     &     t0,t1,ithermal,vold,iperturb,nelemload,
     &     sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &     ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &     ne0,ipkon,thicke,
     &     integerglob,doubleglob,tieset,istartset,iendset,ialset,ntie,
     &     nasym,ielprop,prop)
!     
!     computation of the element matrix and rhs for the user element
!     of type u1
!     
!     U1: This is a beam type element. Reference:
!     Yunhua Luo, An Efficient 3D Timoshenko Beam Element with
!     Consistent Shape Functions, Adv. Theor. Appl. Mech., Vol. 1,
!     2008, no. 3, 95-106
!     
!     special case for which the beam axis goes through the
!     center of gravity of the cross section and the 1-direction
!     corresponds with a principal axis
!     
!     US45: Flat shell 4 nodes 5 dof + (artificial rotz stiff)
!     
      implicit none
!     
      integer mass(*),stiffness,buckling,rhsi,coriolis
!     
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
      character*81 tieset(3,*)
!     
      integer konl(26),nelemload(2,*),nbody,nelem,mi(*),kon(*),
     &     ielprop(*),null,index,mattyp,ithermal(*),iperturb(*),nload,
     &     idist,
     &     i,j,i1,nmethod,kk,nelcon(2,*),nrhcon(*),nalcon(2,*),
     &     ielmat(mi(3),*),ielorien(mi(3),*),ipkon(*),indexe,
     &     ntmat_,nope,norien,ihyper,iexpl,kode,imat,iorien,istiff,
     &     ncmat_,intscheme,istep,iinc,iflag,ipompc(*),nodempc(3,*),
     &     nmpc,ikmpc(*),ilmpc(*),ne0,ndof,istartset(*),iendset(*),
     &     ialset(*),ntie,integerglob(*),nasym,nplicon(0:ntmat_,*),
     &     nplkcon(0:ntmat_,*),npmat_
!     
      real*8 co(3,*),xl(3,26),veold(0:mi(2),*),rho,s(60,60),bodyfx(3),
     &     ff(60),elconloc(ncmat_),coords(3),p1(3),
     &     p2(3),eth(6),rhcon(0:1,ntmat_,*),reltime,prop(*),
     &     alcon(0:6,ntmat_,*),alzero(*),orab(7,*),t0(*),t1(*),
     &     xloadold(2,*),vold(0:mi(2),*),xload(2,*),omx,e,un,um,tt,
     &     sm(60,60),sti(6,mi(1),*),stx(6,mi(1),*),t0l,t1l,coefmpc(*),
     &     elas(21),thicke(mi(3),*),doubleglob(*),dl,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),plconloc(802),dtime,ttime,time,
     &     a,xi11,xi12,xi22,xk,e1(3),offset1,offset2,y1,y2,y3,z1,z2,z3,
     &     elcon(0:ncmat_,ntmat_,*)
!     
      if(lakonl(2:2).eq.'1') then
        call e_c3d_u1(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &       ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &       alzero,ielmat,ielorien,norien,orab,ntmat_,
     &       t0,t1,ithermal,vold,iperturb,nelemload,
     &       sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &       nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &       matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &       ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &       ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &       ne0,ipkon,thicke,integerglob,doubleglob,tieset,istartset,
     &       iendset,ialset,ntie,nasym,ielprop,prop)
!        
      elseif(lakonl(2:4).eq.'S45') then
        call e_c3d_us45(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &       ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &       alzero,ielmat,ielorien,norien,orab,ntmat_,
     &       t0,t1,ithermal,vold,iperturb,nelemload,
     &       sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &       nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &       matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &       ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &       ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &       ne0,ipkon,thicke,integerglob,doubleglob,tieset,istartset,
     &       iendset,ialset,ntie,nasym,ielprop,prop)
!        
      elseif(lakonl(2:3).eq.'S3') then
        call e_c3d_us3(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &       ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &       alzero,ielmat,ielorien,norien,orab,ntmat_,
     &       t0,t1,ithermal,vold,iperturb,nelemload,
     &       sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &       nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &       matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &       ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &       ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &       ne0,ipkon,thicke,integerglob,doubleglob,tieset,istartset,
     &       iendset,ialset,ntie,nasym,ielprop,prop)     
      else
        write(*,*) '*ERROR in e_c3d_u.f: user element'
        write(*,*) '       ',lakonl(1:5),' is not defined'
        call exit(201)
      endif
!     
      return
      end
