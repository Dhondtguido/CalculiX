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
      subroutine createlocalsys(nnfront,istartfront,iendfront,ifront,
     &     co,xt,xn,xa,nfront,ifrontrel,stress,iedno,ibounedg,ieled,
     &     kontri,isubsurffront,istartcrackfro,iendcrackfro,ncrack,
     &     angle,nstep,ier)
!     
!     creating a local coordinate system at the crack front
!     
      implicit none
!     
      integer nnfront,i,j,k,j1,j2,istartfront(*),iendfront(*),ifront(*),
     &     node,nodenext,nodelast,noderel,ier,iedge,iedgerel,ielem,
     &     nfront,ifrontrel(*),iedno(2,*),ibounedg(*),ieled(2,*),
     &     kontri(3,*),isubsurffront(*),istartcrackfro(*),ncrack,
     &     iendcrackfro(*),nstep,n1,n2,n3
!     
      real*8 co(3,*),xt(3,*),xn(3,*),xa(3,*),dd,xnor(3),
     &     stress(6,nstep,*),xtj2m1(3),angle(*),
     &     xn1(3),xn2(3),p12(3),p23(3),det,pi
!
      pi=4.d0*datan(1.d0)
!
      do i=1,nnfront
        j1=istartfront(i)
        j2=iendfront(i)
!     
!     calculate the tangent vector t
!     at the end nodes: tangent to the adjacent front edge        
!     in all other nodes: mean of the tangent to the adjacent front edges
!     
        do j=j1,j2
          node=ifront(j)
          if(j.lt.j2) then
            nodenext=ifront(j+1)
            do k=1,3
              xt(k,j)=co(k,nodenext)-co(k,node)
            enddo
          elseif(j.eq.j2) then
            if(isubsurffront(i).eq.1) then
              nodenext=ifront(j1)
              do k=1,3
                xt(k,j)=co(k,nodenext)-co(k,node)
              enddo
            else
              nodelast=ifront(j-1)
              do k=1,3
                xt(k,j)=co(k,node)-co(k,nodelast)
              enddo
            endif
          endif
!     
!     normalizing
!     
          dd=dsqrt(xt(1,j)*xt(1,j)+xt(2,j)*xt(2,j)+xt(3,j)*xt(3,j))
          do k=1,3
            xt(k,j)=xt(k,j)/dd
          enddo
        enddo
!     
!     taking the mean and normalizing (due to normalizing the factor of
!     2 is not important)
!     
        if(isubsurffront(i).eq.1) then
          do k=1,3
            xtj2m1(k)=xt(k,j2-1)
          enddo
        endif
!     
        do j=j2-1,j1+1,-1
          do k=1,3
            xt(k,j)=(xt(k,j-1)+xt(k,j))
          enddo
          dd=dsqrt(xt(1,j)*xt(1,j)+xt(2,j)*xt(2,j)+xt(3,j)*xt(3,j))
          do k=1,3
            xt(k,j)=xt(k,j)/dd
          enddo
        enddo
!     
!     for subsurface cracks: adjust the tangent at the starting and
!     end node
!     
        if(isubsurffront(i).eq.1) then
          do k=1,3
            xt(k,j1)=xt(k,j1)+xt(k,j2)
          enddo
          dd=dsqrt(xt(1,j1)*xt(1,j1)+xt(2,j1)*xt(2,j1)+
     &         xt(3,j1)*xt(3,j1))
          do k=1,3
            xt(k,j1)=xt(k,j1)/dd
          enddo
          do k=1,3
            xt(k,j2)=xt(k,j2)+xtj2m1(k)
          enddo
          dd=dsqrt(xt(1,j2)*xt(1,j2)+xt(2,j2)*xt(2,j2)+
     &         xt(3,j2)*xt(3,j2))
          do k=1,3
            xt(k,j2)=xt(k,j2)/dd
          enddo
        else
!     
!     for surface cracks: calculate the angle between the tangents
!     at the crossing points of the crack fronts with the free
!     surface; can be used to determine an appropriate shape factor
!     in shapefactor.f
!     
          angle(i)=xt(1,j1)*xt(1,j2)+xt(2,j1)*xt(2,j2)+xt(3,j1)*xt(3,j2)
          angle(i)=min(max(angle(i),-1.d0),1.d0)
          angle(i)=dacos(angle(i))
        endif
      enddo
!     
!     creating a local system based on
!     - the tangential vector t
!     - the normal vector on the adjacent triangles of the crack mesh
!     - a = t x n
!     the direction of n corresponds according to the corkscrew rule
!     with the node numbering of the crack elements      
!     
      do i=1,ncrack
        do j=istartcrackfro(i),iendcrackfro(i)
!     
!     normal on the triangle belonging to the one adjacent edge
!     
          noderel=ifrontrel(j)
          iedgerel=iedno(1,noderel)
          iedge=ibounedg(iedgerel)
          ielem=ieled(1,iedge)
!     
          n1=kontri(1,ielem)
          n2=kontri(2,ielem)
          n3=kontri(3,ielem)
          do k=1,3
            p12(k)=co(k,n2)-co(k,n1)
            p23(k)=co(k,n3)-co(k,n2)
          enddo
          xn1(1)=p12(2)*p23(3)-p12(3)*p23(2)
          xn1(2)=p12(3)*p23(1)-p12(1)*p23(3)
          xn1(3)=p12(1)*p23(2)-p12(2)*p23(1)
!     
!     normal on the triangle belonging to the other adjacent edge
!     
          iedgerel=iedno(2,noderel)
          iedge=ibounedg(iedgerel)
          ielem=ieled(1,iedge)
!     
          n1=kontri(1,ielem)
          n2=kontri(2,ielem)
          n3=kontri(3,ielem)
          do k=1,3
            p12(k)=co(k,n2)-co(k,n1)
            p23(k)=co(k,n3)-co(k,n2)
          enddo
          xn2(1)=p12(2)*p23(3)-p12(3)*p23(2)
          xn2(2)=p12(3)*p23(1)-p12(1)*p23(3)
          xn2(3)=p12(1)*p23(2)-p12(2)*p23(1)
!     
!     taking the mean (factor of 2 is not important due to
!     subsequent normalization)
!     
          do k=1,3
            xn(k,j)=xn1(k)+xn2(k)
          enddo
!     
!     projection on a plane orthogonal to the local tangent vector
!     xm.xt=0 must apply
!     
          dd=xn(1,j)*xt(1,j)+xn(2,j)*xt(2,j)+xn(3,j)*xt(3,j)
          do k=1,3
            xn(k,j)=xn(k,j)-dd*xt(k,j)
          enddo
!     
!     normalizing vector xm
!     
          dd=dsqrt(xn(1,j)*xn(1,j)+xn(2,j)*xn(2,j)+xn(3,j)*xn(3,j))
          do k=1,3
            xn(k,j)=xn(k,j)/dd
          enddo
!     
!     propagation direction a=t x n 
!     
          xa(1,j)=xt(2,j)*xn(3,j)-xt(3,j)*xn(2,j)
          xa(2,j)=xt(3,j)*xn(1,j)-xt(1,j)*xn(3,j)
          xa(3,j)=xt(1,j)*xn(2,j)-xt(2,j)*xn(1,j)
        enddo
      enddo
!
!     calculate a fictituous opening angle for each crack front        
!     
      do i=1,nnfront
        j1=istartfront(i)
        j2=iendfront(i)
!     
        if(isubsurffront(i).ne.1) then
!     
!     for surface cracks: calculate the angle between the tangents
!     at the crossing points of the crack fronts with the free
!     surface; can be used to determine an appropriate shape factor
!     in shapefactor.f
!     
          angle(i)=xt(1,j1)*xt(1,j2)+xt(2,j1)*xt(2,j2)+xt(3,j1)*xt(3,j2)
          angle(i)=min(max(angle(i),-1.d0),1.d0)
!     
!     if cos(angle(i))<0, i.e. pi/2 < angle(i) < 3*pi/2, a
!     check is performed whether angle(i) > pi or
!     angle (i) < pi 
!     (the cosine of alpha and (2*pi-alpha) is the same)
!     
          if(angle(i).ge.0.d0) then
            angle(i)=dacos(angle(i))
          else
            angle(i)=dacos(angle(i))
!     
!     mean normal on the front
!     
            do j=1,3
              xnor(j)=(xn(j,j1)+xn(j,j2))
            enddo
!     
!     check whether angle(i) is less than or greater than pi
!     (the cosine of alpha and (2*pi-alpha) is the same
!     
            det=xnor(1)*(xt(2,j1)*xt(3,j2)-xt(3,j1)*xt(2,j2))+
     &           xnor(2)*(xt(3,j1)*xt(1,j2)-xt(1,j1)*xt(3,j2))+
     &           xnor(3)*(xt(1,j1)*xt(2,j2)-xt(2,j1)*xt(1,j2))
            if(det.lt.0.d0) then
              angle(i)=2.d0*pi-angle(i)
            endif
          endif
!          
        endif
      enddo
!     
c     do i=1,nfront
c     write(*,*) 'createlocalsys xt'
c     write(*,*) 'xt ',i,xt(1,i),xt(2,i),xt(3,i)
c     enddo
c     do i=1,nfront
c     write(*,*) 'createlocalsys xn'
c     write(*,*) 'xn ',i,xn(1,i),xn(2,i),xn(3,i)
c     enddo
c     do i=1,nfront
c     write(*,*) 'createlocalsys xa'
c     write(*,*) 'xa ',i,xa(1,i),xa(2,i),xa(3,i)
c     enddo
c     write(*,*)
!     
      return
      end

