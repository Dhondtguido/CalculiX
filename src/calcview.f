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
!
!     center of gravity of the projection of the vertices for
!     visibility purposes
!     exact integration for one triangle: routine cubtri
!     if the surfaces are far enough away, one-point integration
!     is used
! 
      subroutine calcview(sideload,vold,co,pmid,e1,e2,e3,
     &     kontri,nloadtr,adview,auview,dist,idist,area,
     &     ntrit,mi,jqrad,irowrad,nzsrad,sidemean,ntria,ntrib,
     &     covered,ng)
!     
      implicit none
!     
      character*1 covered(ng,ng)
!
      character*20 sideload(*)
!     
      integer ntr,i,j,k,l,mi(*),ntria,ntrib,
     &     ncovered,kontri(4,*),nloadtr(*),
     &     i1,j1,istart,iend,jstart,jend,imin,imid,imax,
     &     k1,kflag,idist(*),ndist,i2,i3,ng,idi,idj,ntri,
     &     ix,iy,ntrit,limev,ier,nw,idata(1),ncalls,
     &     irowrad(*),jqrad(*),nzsrad,i0,nzsradv(3)
!     
      real*8 w(239),vold(0:mi(2),*),co(3,*),xn(3),xxn,xy(ng),
     &     pmid(3,*),e3(4,*),e1(3,*),e2(3,*),p1(3),p2(3),p3(3),
     &     x,y,porigin(3),etmin,etmax,ximin,anglemin,distance,
     &     ximid,ximax,dummy,a(3,3),b(3,3),c(3,3),ddd(3),p31(3),
     &     xi(3),et(3),ftij,adview(*),auview(*),dxi,dir(3),
     &     dirloc(3),dist(*),area(*),dd,p21(3),sidemean,r(3,3),
     &     fform,pl(2,3),epsabs,epsrel,abserr,q(3,3),
     &     rdata(21),factor,argument
!
      external fform
!     
      nzsradv(3)=nzsrad
!     
!     dxi is the grid size in xi- and et-direction     
!     
!     -1.le.xi.le.1 and -1.le.et.le.1  
!     
      dxi=2.d0/ng
      anglemin=dacos((ng/2.d0-1.d0)/(ng/2.d0))
!
!     factor determines when the numerical integration using cubtri
!     is replaced by a simplified formula using only the center
!     of gravity of one of the triangles. The integration over the
!     other triangle is exact (analytical formula, see
!     "Radiosity: a Programmer's Perspective", by Ian Ashdown, Wiley, 1994)
!     If the distance between the center of gravity of the triangles
!     is larger then factor*the projected sqrt(area) of the triangle on the 
!     hemisphere, the simplified formula is taken
!
      factor=0.d0
!     
      do i=ntria,ntrib
         if(area(i).lt.1.d-20) cycle
!     
!        pl contains the coordinates of the vertices of triangle i
!        in the local coordinate system attached to triangle i
!
!        This local coordinate system has its origin in the first node
!        of triangle i (i1 underneath), its local x-axis along the
!        edge from node 1 to node 2 and its local y-axis such that
!        e1 x e2 agrees with the usual corkscrew rule
!     
         i1=kontri(1,i)
         if(i1.eq.0) cycle
         i2=kontri(2,i)
         i3=kontri(3,i)
         do j=1,3
            porigin(j)=co(j,i1)+vold(j,i1)
            p2(j)=co(j,i2)+vold(j,i2)
            p3(j)=co(j,i3)+vold(j,i3)
            p21(j)=p2(j)-porigin(j)
            p31(j)=p3(j)-porigin(j)
         enddo
         pl(1,1)=0.d0
         pl(2,1)=0.d0
         pl(1,2)=dsqrt(p21(1)**2+p21(2)**2+p21(3)**2)
         pl(2,2)=0.d0
         pl(1,3)=p31(1)*e1(1,i)+p31(2)*e1(2,i)+p31(3)*e1(3,i)
         pl(2,3)=p31(1)*e2(1,i)+p31(2)*e2(2,i)+p31(3)*e2(3,i)
!     
!     checking which triangles face triangle i
!     
         ndist=0
         idi=kontri(4,i)
         do j=1,ntrit
            if((kontri(1,j).eq.0).or.(area(j).lt.1.d-20).or.
     &         (i.eq.j)) cycle
!
!           if the angle between the connection of the two triangles
!           and the base plane of the hemisphere is too small (i.e. 
!           triangle j is too close to the horizon) the viewfactor
!           for this triangle is not taken into account
!
            distance=dsqrt((pmid(1,j)-pmid(1,i))**2+
     &           (pmid(2,j)-pmid(2,i))**2+
     &           (pmid(3,j)-pmid(3,i))**2)
            if((pmid(1,j)*e3(1,i)+pmid(2,j)*e3(2,i)+
     &           pmid(3,j)*e3(3,i)+e3(4,i))/distance.le.anglemin) cycle
            if((pmid(1,i)*e3(1,j)+pmid(2,i)*e3(2,j)+
     &           pmid(3,i)*e3(3,j)+e3(4,j))/distance.le.anglemin) cycle
!     
            idj=kontri(4,j)
            if(sideload(nloadtr(idi))(18:20).ne.
     &           sideload(nloadtr(idj))(18:20)) cycle
!     
            ndist=ndist+1
            dist(ndist)=distance
            idist(ndist)=j
         enddo
         if(ndist.eq.0) cycle
!     
!     ordering the triangles
!     
         kflag=2
         call dsort(dist,idist,ndist,kflag)
!     
!     initializing the coverage matrix
!  
         do i1=1,ng
            xy(i1)=((i1-0.5d0)*dxi-1.d0)**2
         enddo
!
         ncovered=0
         do i1=1,ng
            x=xy(i1)
            do j1=1,ng
               if(x+xy(j1).gt.1.d0) then
                  covered(i1,j1)='T'
                  ncovered=ncovered+1
               else
                  covered(i1,j1)='F'
               endif
            enddo
         enddo
!     
         do k1=1,ndist
            j=idist(k1)
!     
!           determining the 2-D projection of the vertices
!           of triangle j
!
!           r is the connection of cg of triangle i with the
!           vertices of triangle j
!
!           xi and et are the local coordinates in the local 
!           coordinate system attached of triangle i of the projection
!           of the vertices of triangle j on the base of the hemisphere
!           at triangle i
!     
            do k=1,3
               r(k,1)=co(k,kontri(1,j))+vold(k,kontri(1,j))-pmid(k,i)
            enddo
            ddd(1)=dsqrt(r(1,1)*r(1,1)+r(2,1)*r(2,1)+r(3,1)*r(3,1))
            do k=1,3
               r(k,1)=r(k,1)/ddd(1)
            enddo
            xi(1)=r(1,1)*e1(1,i)+r(2,1)*e1(2,i)+r(3,1)*e1(3,i)
            et(1)=r(1,1)*e2(1,i)+r(2,1)*e2(2,i)+r(3,1)*e2(3,i)
!     
            do k=1,3
               r(k,2)=co(k,kontri(2,j))+vold(k,kontri(2,j))-pmid(k,i)
            enddo
            ddd(2)=dsqrt(r(1,2)*r(1,2)+r(2,2)*r(2,2)+r(3,2)*r(3,2))
            do k=1,3
               r(k,2)=r(k,2)/ddd(2)
            enddo
            xi(2)=r(1,2)*e1(1,i)+r(2,2)*e1(2,i)+r(3,2)*e1(3,i)
            et(2)=r(1,2)*e2(1,i)+r(2,2)*e2(2,i)+r(3,2)*e2(3,i)
!     
            do k=1,3
               r(k,3)=co(k,kontri(3,j))+vold(k,kontri(3,j))-pmid(k,i)
            enddo
            ddd(3)=dsqrt(r(1,3)*r(1,3)+r(2,3)*r(2,3)+r(3,3)*r(3,3))
            do k=1,3
               r(k,3)=r(k,3)/ddd(3)
            enddo
            xi(3)=r(1,3)*e1(1,i)+r(2,3)*e1(2,i)+r(3,3)*e1(3,i)
            et(3)=r(1,3)*e2(1,i)+r(2,3)*e2(2,i)+r(3,3)*e2(3,i)
!
!     determining the center of gravity of the projected
!     triangle
!     
            dirloc(1)=(xi(1)+xi(2)+xi(3))/3.d0
            dirloc(2)=(et(1)+et(2)+et(3))/3.d0
!     
!     check whether this direction was already covered
!     
            ix=int((dirloc(1)+1.d0)/dxi)+1
            iy=int((dirloc(2)+1.d0)/dxi)+1
            if(covered(ix,iy).eq.'T') then
               cycle
            endif
!     
!     determine the direction vector in global coordinates
!     
            do k=1,3
               dir(k)=(pmid(k,j)-pmid(k,i))/dist(k1)
            enddo
!     
!     direction vector in local coordinates of triangle i
!     
            dirloc(3)=dir(1)*e3(1,i)+dir(2)*e3(2,i)+dir(3)*e3(3,i)
!     
!     if surfaces are close, numerical integration with
!     cubtri is performed
!     
            if(dist(k1).le.factor*dsqrt(area(i))*dirloc(3)) then
!     
!     vertices of triangle j
!     
               do k=1,3
                  do l=1,3
                     q(l,k)=co(l,kontri(k,j))+vold(l,kontri(k,j))
                  enddo
               enddo
!     
!     formfactor contribution
!     
               epsrel=0.01d0
               epsabs=0.d0
               limev=100
               nw=239
               ncalls=0
!
!              storing common data into field rdata
!
               rdata(1)=q(1,1)
               rdata(2)=q(1,2)
               rdata(3)=q(1,3)
               rdata(4)=q(2,1)
               rdata(5)=q(2,2)
               rdata(6)=q(2,3)
               rdata(7)=q(3,1)
               rdata(8)=q(3,2)
               rdata(9)=q(3,3)
               rdata(10)=e1(1,i)
               rdata(11)=e1(2,i)
               rdata(12)=e1(3,i)
               rdata(13)=e2(1,i)
               rdata(14)=e2(2,i)
               rdata(15)=e2(3,i)
               rdata(16)=e3(1,i)
               rdata(17)=e3(2,i)
               rdata(18)=e3(3,i)
               rdata(19)=porigin(1)
               rdata(20)=porigin(2)
               rdata(21)=porigin(3)
!     
!     max 1000 evaluations for nw=239
!     
               call cubtri(fform,pl,epsrel,limev,ftij,abserr,ncalls,
     &              w,nw,idata,rdata,ier)
               ftij=ftij/2.d0
            endif
!
!     ftij is the viewfactor between triangle i and j. In the end
!     the viewfactor between the parent faces is needed. This
!     transformation, and the division by pi is performed in postview.f
!
c            if(dabs(xi(2)-xi(1)).lt.1.d-5) xi(2)=xi(1)+1.d-5
c            if(dabs(xi(2)-xi(1)).lt.1.d-5) xi(2)=xi(1)+1.d-5
!     
!     if the surfaces are far enough away, one-point
!     integration is used
!     
            if(dist(k1).gt.factor*dsqrt(area(i))*dirloc(3)) then
               ftij=0.d0
               do k=1,3
                  l=k-1
                  if(l.lt.1) l=3
                  xn(1)=r(2,k)*r(3,l)-r(2,l)*r(3,k)
                  xn(2)=r(3,k)*r(1,l)-r(3,l)*r(1,k)
                  xn(3)=r(1,k)*r(2,l)-r(1,l)*r(2,k)
                  xxn=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
!     
!     argument of dacos must have an absolute value
!     smaller than or equal to 1.d0; due to
!     round-off the value can slightly exceed one
!     and has to be cut-off.
!     
                  argument=
     &                 r(1,k)*r(1,l)+r(2,k)*r(2,l)+r(3,k)*r(3,l)
                  if(dabs(argument).gt.1.d0) then
                     if(argument.gt.0.d0) then
                        argument=1.d0
                     else
                        argument=-1.d0
                     endif
                  endif
                  ftij=ftij+
     &                 (e3(1,i)*xn(1)
     &                 +e3(2,i)*xn(2)
     &                 +e3(3,i)*xn(3))/xxn
     &                 *dacos(argument)
               enddo
               ftij=ftij*area(i)/2.d0
            endif
!     
            idj=kontri(4,j)
            i0=0
            call add_sm_st_as(auview,adview,jqrad,irowrad,
     &           idi,idj,ftij,i0,i0,nzsradv)          
!
!     ftij is the viewfactor between triangle i and j. In the end
!     the viewfactor between the parent faces is needed. This
!     transformation, and the division by pi is performed in postview.f
!     
!     determining maxima and minima
!     
            ximin=2.d0
            ximax=-2.d0
            do k=1,3
               if(xi(k).lt.ximin) then
                  ximin=xi(k)
                  imin=k
               endif
               if(xi(k).gt.ximax) then
                  ximax=xi(k)
                  imax=k
               endif
            enddo
!     
            if(((imin.eq.1).and.(imax.eq.2)).or.
     &           ((imin.eq.2).and.(imax.eq.1))) then
               imid=3
               ximid=xi(3)
            elseif(((imin.eq.2).and.(imax.eq.3)).or.
     &              ((imin.eq.3).and.(imax.eq.2))) then
               imid=1
               ximid=xi(1)
            else
               imid=2
               ximid=xi(2)
            endif
!     
!     check for equal x-values
!     
            if(ximid-ximin.lt.1.d-5) then
               ximin=ximin-1.d-5
               xi(imin)=ximin
            endif
            if(ximax-ximid.lt.1.d-5) then
               ximax=ximax+1.d-5
               xi(imax)=ximax
            endif
!     
!     equation of the straight lines connecting the
!     triangle vertices in the local x-y plane
!     
            a(1,2)=et(2)-et(1)
            b(1,2)=xi(1)-xi(2)
            c(1,2)=et(1)*xi(2)-xi(1)*et(2)
!     
            a(2,3)=et(3)-et(2)
            b(2,3)=xi(2)-xi(3)
            c(2,3)=et(2)*xi(3)-xi(2)*et(3)
!     
            a(3,1)=et(1)-et(3)
            b(3,1)=xi(3)-xi(1)
            c(3,1)=et(3)*xi(1)-xi(3)*et(1)
!     
            a(2,1)=a(1,2)
            b(2,1)=b(1,2)
            c(2,1)=c(1,2)
            a(3,2)=a(2,3)
            b(3,2)=b(2,3)
            c(3,2)=c(2,3)
            a(1,3)=a(3,1)
            b(1,3)=b(3,1)
            c(1,3)=c(3,1)
!
!           xmin.le.xi.le.xmid
!
            istart=int((ximin+1.d0+dxi/2.d0)/dxi)+1
            iend=int((ximid+1.d0+dxi/2.d0)/dxi)
            do i1=istart,iend
               x=dxi*(i1-0.5d0)-1.d0
               etmin=-(a(imin,imid)*x+c(imin,imid))/b(imin,imid)
               etmax=-(a(imin,imax)*x+c(imin,imax))/b(imin,imax)
               if(etmin.gt.etmax) then
                  dummy=etmin
                  etmin=etmax
                  etmax=dummy
                endif
!
!              et is between lower line and upper line
!              lower and upper line connect ximin with ximid or ximax
!
               jstart=int((etmin+1.d0+dxi/2.d0)/dxi)+1
               jend=int((etmax+1.d0+dxi/2.d0)/dxi)
               do j1=jstart,jend
                  covered(i1,j1)='T'
               enddo
               ncovered=ncovered+jend-jstart+1
            enddo
!
!           xmid.le.xi.le.xmax
!
            istart=int((ximid+1.d0+dxi/2.d0)/dxi)+1
            iend=int((ximax+1.d0+dxi/2.d0)/dxi)
            do i1=istart,iend
               x=dxi*(i1-0.5d0)-1.d0
               etmin=-(a(imid,imax)*x+c(imid,imax))/b(imid,imax)
               etmax=-(a(imin,imax)*x+c(imin,imax))/b(imin,imax)
               if(etmin.gt.etmax) then
                  dummy=etmin
                  etmin=etmax
                  etmax=dummy
               endif
!
!              et is between lower line and upper line
!              lower and upper line connect xmin or xmid with ximax
!
               jstart=int((etmin+1.d0+dxi/2.d0)/dxi)+1
               jend=int((etmax+1.d0+dxi/2.d0)/dxi)
               do j1=jstart,jend
                  covered(i1,j1)='T'
               enddo
               ncovered=ncovered+jend-jstart+1
            enddo
            if(ncovered.eq.ng*ng)exit
!     
         enddo
      enddo
!
      return
      end
!     
!     function to be integrated
!     
      real*8 function fform(x,y,idata,rdata)
!
      implicit none
!     
      integer k,l,idata(1)
!     
      real*8 pint(3),ddd(3),xn(3),q(3,3),
     &   unitvec(3,3),r(3,3),xxn,x,y,porigin(3),rdata(21)
!
!     retrieving common data from field rdata
!
      q(1,1)=rdata(1)
      q(1,2)=rdata(2)
      q(1,3)=rdata(3)
      q(2,1)=rdata(4)
      q(2,2)=rdata(5)
      q(2,3)=rdata(6)
      q(3,1)=rdata(7)
      q(3,2)=rdata(8)
      q(3,3)=rdata(9)
      unitvec(1,1)=rdata(10)
      unitvec(2,1)=rdata(11)
      unitvec(3,1)=rdata(12)
      unitvec(1,2)=rdata(13)
      unitvec(2,2)=rdata(14)
      unitvec(3,2)=rdata(15)
      unitvec(1,3)=rdata(16)
      unitvec(2,3)=rdata(17)
      unitvec(3,3)=rdata(18)
      porigin(1)=rdata(19)
      porigin(2)=rdata(20)
      porigin(3)=rdata(21)
!
      do k=1,3
         pint(k)=porigin(k)+x*unitvec(k,1)+y*unitvec(k,2)
      enddo
!
      do k=1,3
         r(k,1)=q(k,1)-pint(k)
      enddo
      ddd(1)=dsqrt(r(1,1)*r(1,1)+r(2,1)*r(2,1)+r(3,1)*r(3,1))
!
      do k=1,3
         r(k,2)=q(k,2)-pint(k)
      enddo
      ddd(2)=dsqrt(r(1,2)*r(1,2)+r(2,2)*r(2,2)+r(3,2)*r(3,2))
!
      do k=1,3
         r(k,3)=q(k,3)-pint(k)
      enddo
      ddd(3)=dsqrt(r(1,3)*r(1,3)+r(2,3)*r(2,3)+r(3,3)*r(3,3))
!     
!     calculating the contribution
!     
      fform=0.d0
      do k=1,3
         l=k-1
         if(l.lt.1) l=3
         xn(1)=r(2,k)*r(3,l)-r(2,l)*r(3,k)
         xn(2)=r(3,k)*r(1,l)-r(3,l)*r(1,k)
         xn(3)=r(1,k)*r(2,l)-r(1,l)*r(2,k)
         xxn=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
         fform=fform+
     &        (unitvec(1,3)*xn(1)
     &        +unitvec(2,3)*xn(2)
     &        +unitvec(3,3)*xn(3))/xxn
     &        *dacos(
     &        (r(1,k)*r(1,l)+r(2,k)*r(2,l)+r(3,k)*r(3,l))/
     &        (ddd(k)*ddd(l)))
      enddo
!
      return
      end
      
