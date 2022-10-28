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
      subroutine updategeodata(nktet,netet_,h,d,dmin,ipoed,iedg,cotet,
     &     planfa,bc,cg,kontet,ifac,ipofa,doubleglob,integerglob,ipoeln)
!     
!     calculating the desired size of h in all nodes of the actual
!     mesh
!     
      implicit none
!     
      integer nktet,netet_,ipoed(*),iedg(3,*),kontet(4,*),ifac(4,*),
     &     ipofa(*),integerglob(*),ipoeln(*),
     &     nodef(3),nodes(4),iselect(1),nselect,nterms,nktri,nkon,
     &     nfield,nfaces,netet,nelem,ne,n1,n2,i,j,ialset(1),ielement,
     &     iendset(1),iface,imastset,istartset(1),index,konl(20),
     &     loopa
!     
      real*8 h(*),d(*),dmin,cotet(3,*),planfa(4,*),bc(4,*),cg(3,*),
     &     doubleglob(*),x,y,z,r,ratio(20),p1x,p1y,p1z,p2x,p2y,p2z,p3x,
     &     p3y,p3z,p4x,p4y,p4z,a11,a12,a13,a21,a22,a23,a31,a32,a33,d1,
     &     d2,d3,det,coords(3),dist
!     
!     determine the size of all edges
!     
      dmin=1.d30
!     
      loop: do i=1,nktet
        index=ipoed(i)
        do
          if(index.eq.0) cycle loop
!     
          n1=iedg(1,index)
          n2=iedg(2,index)
!     
          d(index)=dsqrt((cotet(1,n1)-cotet(1,n2))**2+
     &         (cotet(2,n1)-cotet(2,n2))**2+
     &         (cotet(3,n1)-cotet(3,n2))**2)
!     
          if(d(index).lt.dmin) dmin=d(index)
!     
          index=iedg(3,index)
        enddo
      enddo loop
!     
!     calculating the desired edge size through interpolation
!     
!     initializing fields
!     
      nktri=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=1
      nselect=1
      iselect(1)=1
      imastset=0
!     
      do i=1,nktet
        if(ipoeln(i).eq.0) cycle
!     
!     perform the interpolation for the internal node
!     
        do j=1,3
          coords(j)=cotet(j,i)
        enddo
        loopa=8
        call basis(doubleglob(1),doubleglob(netet+1),
     &       doubleglob(2*netet+1),
     &       doubleglob(3*netet+1),doubleglob(4*netet+1),
     &       doubleglob(5*netet+1),integerglob(6),
     &       integerglob(netet+6),
     &       integerglob(2*netet+6),doubleglob(6*netet+1),
     &       integerglob(3*netet+6),nktri,netet,
     &       doubleglob(4*nfaces+6*netet+1),nfield,
     &       doubleglob(nktri+4*nfaces+6*netet+1),
     &       integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &       integerglob(2*ne+7*netet+6),
     &       integerglob(nkon+2*ne+7*netet+6),
     &       coords(1),coords(2),coords(3),h(i),ratio,iselect,
     &       nselect,istartset,iendset,ialset,imastset,
     &       integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,
     &       loopa,dist)
      enddo
!     
!     updating the cg and bc fields
!     
      do i=1,netet_
        ielement=i
!     
        if(kontet(1,ielement).eq.0) cycle
!     
        do j=1,4
          nodes(j)=kontet(j,ielement)
        enddo
!     
!     finding the center and radius of the circumscribed sphere
!     
        p1x=cotet(1,nodes(1))
        p1y=cotet(2,nodes(1))
        p1z=cotet(3,nodes(1))
!     
        p2x=cotet(1,nodes(2))
        p2y=cotet(2,nodes(2))
        p2z=cotet(3,nodes(2))
!     
        p3x=cotet(1,nodes(3))
        p3y=cotet(2,nodes(3))
        p3z=cotet(3,nodes(3))
!     
        p4x=cotet(1,nodes(4))
        p4y=cotet(2,nodes(4))
        p4z=cotet(3,nodes(4))
!     
        a11=p1x-p2x
        a12=p1y-p2y
        a13=p1z-p2z
        d1=((p1x+p2x)*a11+(p1y+p2y)*a12+(p1z+p2z)*a13)/2.d0
!     
        a21=p1x-p3x
        a22=p1y-p3y
        a23=p1z-p3z
        d2=((p1x+p3x)*a21+(p1y+p3y)*a22+(p1z+p3z)*a23)/2.d0
!     
        a31=p1x-p4x
        a32=p1y-p4y
        a33=p1z-p4z
        d3=((p1x+p4x)*a31+(p1y+p4y)*a32+(p1z+p4z)*a33)/2.d0
!     
        det=a11*(a22*a33-a32*a23)-a12*(a21*a33-a31*a23)+
     &       a13*(a21*a32-a31*a22)
        x=(d1*(a22*a33-a23*a32)-d2*(a12*a33-a32*a13)+
     &       d3*(a12*a23-a22*a13))/det
        y=(-d1*(a21*a33-a31*a23)+d2*(a11*a33-a31*a13)-
     &       d3*(a11*a23-a21*a13))/det
        z=(d1*(a21*a32-a31*a22)-d2*(a11*a32-a31*a12)+
     &       d3*(a11*a22-a21*a12))/det
!     
        r=dsqrt((x-p1x)**2+(y-p1y)**2+(z-p1z)**2)
!     
        bc(1,ielement)=x
        bc(2,ielement)=y
        bc(3,ielement)=z
        bc(4,ielement)=r
!     
!     determining the center of gravity of the element
!     
        cg(1,ielement)=(p1x+p2x+p3x+p4x)/4.d0
        cg(2,ielement)=(p1y+p2y+p3y+p4y)/4.d0
        cg(3,ielement)=(p1z+p2z+p3z+p4z)/4.d0
!     
      enddo
!     
!     update planfa field
!     
      loop2:do i=1,nktet
        iface=ipofa(i)
        do
          if(iface.eq.0) cycle loop2
          nodef(1)=ifac(1,iface)
          nodef(2)=ifac(2,iface)
          nodef(3)=ifac(3,iface)
!     
          call planeeq(cotet,nodef,planfa(1,iface))
!     
          iface=ifac(4,iface)
        enddo
      enddo loop2
!     
      return
      end

