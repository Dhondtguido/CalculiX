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
!     cuts a triangle of the master surface with a slave surface
!     inserts new active edges into iactiveline for current triangle
!
      subroutine treatmasterface(
     &  nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &  ipe,ime,iactiveline,nactiveline,
     &  ifreeintersec,nelemm,nintpoint,pslavsurf,
     &  xl2m,nnodelem,xl2m2,nmp,nodem,areaslav)
!     
!    Author: Saskia Sitzmann     
!     
      implicit none
!
      integer nvertex,nopes,ipe(*),ime(4,*),iactiveline(3,*),
     &  nactiveline,ifreeintersec,nmp,i,j,k,nintpoint,
     &  nnodelem,nodem(*),modf,nelemm,k_max,nipold
!
      real*8 pvertex(3,13),slavstraight(36),xn(3),xilm,etlm,xnl(3),
     &  xl2s(3,*),p1(2),p2(2),pslavsurf(3,*),xil,etl,p(3),dist,
     &  area,xl2m(3,8),xl2m2(3,8),al,err,xns(3,8),
     &  xl2sp(3,*),xl2mp(3,8),cgp(3),pm(3),ps(3),xit(3),etat(3),areaslav
!
!
!     
      include "gauss.f"
!     
      err=1.d-6
      nvertex=0
      nipold=nintpoint
!     
!     Project master nodes to meanplane, needed for Sutherland-Hodgman
!     
      do j=1, nmp
         al=-xn(1)*xl2m2(1,j)-xn(2)*
     &        xl2m2(2,j)-xn(3)*
     &        xl2m2(3,j)-slavstraight(nopes*4+4)
         do k=1,3
            xl2mp(k,j)= xl2m2(k,j)+al*xn(k)    
         enddo
      enddo 
!     
!     call Sutherland-Hodgman Algo
!     
      call sutherland_hodgman(nopes,xn,xl2sp,xl2mp,nodem,
     &     ipe,ime,iactiveline,nactiveline,
     &     ifreeintersec,nelemm,nmp,
     &     nvertex,pvertex) 
!     
!     
      do k=1,3
         cgp(k)=0.0d0
      enddo
      if(nvertex.lt.3) return       
!     
      if(nvertex.eq.3)then
         do k=1,3
            cgp(k)=pvertex(k,nvertex)
         enddo
         nvertex=nvertex-1
         k_max=1
      else
         do i=1,nvertex
            do k=1,3
               cgp(k)=cgp(k)+pvertex(k,i)/nvertex
            enddo
         enddo
         k_max=nvertex
      endif 
!     
!     Project center point back on slave face
!     
      call attachline(xl2s,cgp,nopes,xit(3),etat(3),xn,p,dist)
!     
!     generating integration points on the slave surface S
!     
      do k=1,k_max
!     
!     Project back on slave surface
!     
         call attachline(xl2s,pvertex(1:3,modf(nvertex,k)),
     &        nopes,xit(1),etat(1),xn,p,dist)
         call attachline(xl2s,pvertex(1:3,modf(nvertex,k+1)),
     &        nopes,xit(2),etat(2),xn,p,dist)
!
         p1(1)=xit(1)-xit(3)
         p1(2)=etat(1)-etat(3)
!
         p2(1)=xit(2)-xit(3)
         p2(2)=etat(2)-etat(3)
!
         area=dabs(p1(1)*p2(2)-p2(1)*p1(2))/2.d0
!
         if(area.lt.1.d-4) cycle
         if(nopes.eq.4.and.areaslav+area-4.0d0.gt.1.d-3
     &        .and.nactiveline.gt.0)then
           nactiveline=0
           return
         endif
         if(nopes.eq.3.and.areaslav+area-0.5d0.gt.1.d-4
     &       .and.nactiveline.gt.0)then
           nactiveline=0
           return
         endif
         areaslav=areaslav+area
!     
!     7 points scheme
!     
         do i=1,7
            xil= xit(3)*gauss2d6(1,i)+
     &           xit(1)*gauss2d6(2,i)+
     &           xit(2)*(1-gauss2d6(1,i)-gauss2d6(2,i))
            
            etl= etat(3)*gauss2d6(1,i)+
     &           etat(1)*gauss2d6(2,i)+
     &           etat(2)*(1-gauss2d6(1,i)-gauss2d6(2,i))
!     
            call evalshapefunc(xil,etl,xns,nopes,xnl)
            call evalshapefunc(xil,etl,xl2s,nopes,ps)
!     
            nintpoint=nintpoint+1
!     
!     projection of the integration point in the mean
!     slave plane onto the slave surface
!     
!     projection of the master integration point onto the
!     master surface in order to get the local coordinates
!     own xn for every integration point?
!     
            call attachline(xl2m,ps,nnodelem,xilm,etlm,xn,p,dist)
            call evalshapefunc(xilm,etlm,xl2m,nnodelem,pm)   
!
            pslavsurf(1,nintpoint)=xil
            pslavsurf(2,nintpoint)=etl
!
!           weights sum up to 0.5 for triangles in local coordinates
!
            pslavsurf(3,nintpoint)=2.d0*area*weight2d6(i)
         enddo
      enddo
!
      return
      end
