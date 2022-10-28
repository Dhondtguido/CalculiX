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
      subroutine crackprop(ifrontrel,ibounnod,domphi,da,co,costruc,nk,
     &     xa,xn,nnfront,istartfront,iendfront,doubleglob,integerglob,
     &     isubsurffront,dadn,ncyc,ifrontprop,nstep,acrack,acrackglob,
     &     datarget,iincglob,iinc,dnglob,ncyctot,ier)
!     
!     calculate the crack propagation increment
!
!     acrackglob, iincglob and dnglob are the crack length, the
!     increment number and the total number of cycles at the END
!     of the present increment and are therefore attached to the
!     propagated front, therefore they get values assigned in the
!     present routine; notice that the propagated front MAY be
!     changed in eqspacednodes.f
!      
      implicit none
!     
      integer i,j,k,m,node,ifrontrel(*),nk,nnfront,istartfront(*),
     &     iendfront(*),integerglob(*),nktet,netet,ne,nkon,nfaces,
     &     nfield,nselect,imastset,iselect(6),nterms,nelem,
     &     ialset(1),iendset(1),istartset(1),konl(20),loopa,
     &     noderel,ibounnod(*),isubsurffront(*),ifrontprop(*),
     &     ncyc,nstep,nodep,nodeq,iincglob(*),iinc,
     &     ncyctot,ier
!     
      real*8 da(*),domphi(*),co(3,*),costruc(3,*),xa(3,*),xn(3,*),
     &     doubleglob(*),coords(3),ratio(20),dist,pi,theta,ctheta,
     &     stheta,acrack(*),c(3,3),r0(3),r(3),dadn(*),dnglob(*),
     &     value(1),acrackglob(*),datarget,p(3),q(3),dd,al
!     
      nktet=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=13
!     
      nselect=0
      imastset=0
      loopa=8
!     
!     Definition variable theta,ctheta,stheta
!     
      pi=4.d0*datan(1.d0)
      theta=5.d0*pi/180.d0
!     
      do k=1,nnfront
!     
        do i=istartfront(k),iendfront(k)
!     
          noderel=ifrontrel(i)
          node=ibounnod(noderel)
          da(i)=dadn(i)*ncyc
!     
!     if the crack propagation is too small (less than 1 % of
!     the target):
!     1) 1 % of the target is taken
!     2) the deflection angle is set to zero
!     3) equal spacing is deactivated (if the crack propagation
!     is too small, moving the nodes may lead to overlap for
!     large front curvatures
!     
          if(da(i).lt.datarget*0.01d0) then
            da(i)=datarget*0.01d0
            domphi(i)=0.d0
          endif
!     
!     increase da(i) such that the node extends by 1.2d-6 outside
!     the structure when rotated such that the propagation orthogonal
!     to the crack front is da(i)
!     
          da(i)=dsqrt((1.2d-6)**2+da(i)**2)*dsign(1.d0,da(i))
!     
          nk=nk+1
          ifrontprop(i)=nk
!     
!     store the total crack length in the node
!     
          acrackglob(nk)=acrack(i)+da(i)
          iincglob(nk)=iinc+1
          dnglob(nk)=1.d0*ncyctot
!     
          do j=1,3
            co(j,nk)=co(j,node)+(xa(j,i)*dcos(domphi(i))+
     &           xn(j,i)*dsin(domphi(i)))*da(i)
          enddo
!
!     rotating the propagated end nodes of each front about the
!     non-propagated nodes in such a way that they lie outside 
!     the structure
!
          if(isubsurffront(k).eq.1) cycle
!     
          if (i.eq.istartfront(k)) then
!
!     correcting the projection of the node outside the structure
!     onto the structure: this node can lie way outside the local
!     tangent plane to the crack if the crack plane is not
!     orthogonal to the free surface.
!
            nodep=node
            nodeq=ibounnod(ifrontrel(i+1))
            do j=1,3
              p(j)=co(j,nodep)
              q(j)=co(j,nodeq)
              r(j)=costruc(j,noderel)
            enddo
            dd=(r(1)-q(1))**2+(r(2)-q(2))**2+(r(3)-q(3))**2
            al=dd/((r(1)-q(1))*(p(1)-q(1))+
     &           (r(2)-q(2))*(p(2)-q(2))+
     &           (r(3)-q(3))*(p(3)-q(3)))
!            
!           check for the range of al            
!            
            if((al.lt.0.d0).or.(al.gt.1.d0)) then
              al=1.d0
            endif
!
            do j=1,3
              costruc(j,noderel)=q(j)+al*(p(j)-q(j))
            enddo
!
!           r is the propagation vector
!
            do j=1,3
              r(j)=(xa(j,i)*dcos(domphi(i))+
     &             xn(j,i)*dsin(domphi(i)))*da(i)
              co(j,nk)=costruc(j,noderel)+r(j)
              r0(j)=r(j)
            enddo
!     
            ctheta=dcos(theta)
            stheta=-dsin(theta)
!
!           rotating vector in steps of 5 degrees
!
            ier=1
            do m=1,72
!     
              do j=1,3
                coords(j)=co(j,nk)
              enddo
              call basis(doubleglob(1),doubleglob(netet+1),
     &             doubleglob(2*netet+1),doubleglob(3*netet+1),
     &             doubleglob(4*netet+1),doubleglob(5*netet+1),
     &             integerglob(6),integerglob(netet+6),
     &             integerglob(2*netet+6),doubleglob(6*netet+1),
     &             integerglob(3*netet+6),nktet,netet,
     &             doubleglob(4*nfaces+6*netet+1),nfield,
     &             doubleglob(nstep*13*nktet+4*nfaces+6*netet+1),
     &             integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &             integerglob(2*ne+7*netet+6),
     &             integerglob(nkon+2*ne+7*netet+6),coords(1),
     &             coords(2),coords(3),value,ratio,iselect,
     &             nselect,istartset,iendset,ialset,imastset,
     &             integerglob(nkon+2*ne+8*netet+6),nterms,konl,
     &             nelem,loopa,dist)
!     
              if(dist.ge.1.d-6) then
                ier=0
                exit
              endif
!     
              c(1,1)=ctheta+(1-ctheta)*(xn(1,i)**2)
              c(1,2)=-stheta*xn(3,i)+(1-ctheta)*xn(1,i)*xn(2,i)
              c(1,3)=stheta*xn(2,i)+(1-ctheta)*xn(1,i)*xn(3,i)
              c(2,1)=stheta*xn(3,i)+(1-ctheta)*xn(2,i)*xn(1,i)
              c(2,2)=ctheta+(1-ctheta)*(xn(2,i)**2)
              c(2,3)=-stheta*xn(1,i)+(1-ctheta)*xn(2,i)*xn(3,i)
              c(3,1)=-stheta*xn(2,i)+(1-ctheta)*xn(3,i)*xn(1,i)
              c(3,2)=stheta*xn(1,i)+(1-ctheta)*xn(3,i)*xn(2,i)
              c(3,3)=ctheta+(1-ctheta)*(xn(3,i)**2)
!     
              do j=1,3
                r0(j)=r(j)
              enddo
!     r=C*r0
              do j=1,3
                r(j)=c(j,1)*r0(1)+c(j,2)*r0(2)+c(j,3)*r0(3)
              enddo
!     Rotation of node nk 
              do j=1,3
                co(j,nk)=costruc(j,noderel)+r(j)
              enddo
            enddo
!     
!           if no position external to the structure was found: error     
!     
            if(ier.eq.1) then
              return
            endif
!     
          elseif (i.eq.iendfront(k)) then
!
!     correcting the projection of the node outside the structure
!     onto the structure: this node can lie way outside the local
!     tangent plane to the crack if the crack plane is not
!     orthogonal to the free surface.
!
            nodep=ibounnod(ifrontrel(i-1))
            nodeq=node
            do j=1,3
              p(j)=co(j,nodep)
              q(j)=co(j,nodeq)
              r(j)=costruc(j,noderel)
            enddo
            dd=(r(1)-q(1))**2+(r(2)-q(2))**2+(r(3)-q(3))**2
            al=dd/((r(1)-q(1))*(p(1)-q(1))+
     &           (r(2)-q(2))*(p(2)-q(2))+
     &           (r(3)-q(3))*(p(3)-q(3)))
!            
!           check for the range of al            
!            
            if((al.lt.0.d0).or.(al.gt.1.d0)) then
              al=0.d0
            endif
!            
            do j=1,3
              costruc(j,noderel)=q(j)+al*(p(j)-q(j))
            enddo
!
!           r is the propagation vector
!
            do j=1,3
              r(j)=(xa(j,i)*dcos(domphi(i))+
     &             xn(j,i)*dsin(domphi(i)))*da(i)
              co(j,nk)=costruc(j,noderel)+r(j)
              r0(j)=r(j)
            enddo
!     
            ctheta=dcos(theta)
            stheta=dsin(theta)
!
!           rotating vector in steps of 5 degrees
!
            ier=1
            do m=1,72
!
              do j=1,3
                coords(j)=co(j,nk)
              enddo
              call basis(doubleglob(1),doubleglob(netet+1),
     &             doubleglob(2*netet+1),doubleglob(3*netet+1),
     &             doubleglob(4*netet+1),doubleglob(5*netet+1),
     &             integerglob(6),integerglob(netet+6),
     &             integerglob(2*netet+6),doubleglob(6*netet+1),
     &             integerglob(3*netet+6),nktet,netet,
     &             doubleglob(4*nfaces+6*netet+1),nfield,
     &             doubleglob(nstep*13*nktet+4*nfaces+6*netet+1),
     &             integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &             integerglob(2*ne+7*netet+6),
     &             integerglob(nkon+2*ne+7*netet+6),coords(1),
     &             coords(2),coords(3),value,ratio,iselect,
     &             nselect,istartset,iendset,ialset,imastset,
     &             integerglob(nkon+2*ne+8*netet+6),nterms,konl,
     &             nelem,loopa,dist)
!     
              if(dist.ge.1.d-6) then
                ier=0
                exit
              endif
!     
              c(1,1)=ctheta+(1-ctheta)*(xn(1,i)**2)
              c(1,2)=-stheta*xn(3,i)+(1-ctheta)*xn(1,i)*xn(2,i)
              c(1,3)=stheta*xn(2,i)+(1-ctheta)*xn(1,i)*xn(3,i)
              c(2,1)=stheta*xn(3,i)+(1-ctheta)*xn(2,i)*xn(1,i)
              c(2,2)=ctheta+(1-ctheta)*(xn(2,i)**2)
              c(2,3)=-stheta*xn(1,i)+(1-ctheta)*xn(2,i)*xn(3,i)
              c(3,1)=-stheta*xn(2,i)+(1-ctheta)*xn(3,i)*xn(1,i)
              c(3,2)=stheta*xn(1,i)+(1-ctheta)*xn(3,i)*xn(2,i)
              c(3,3)=ctheta+(1-ctheta)*(xn(3,i)**2)
!     
              do j=1,3
                r0(j)=r(j)
              enddo
!     r=C*r0
              do j=1,3
                r(j)=c(j,1)*r0(1)+c(j,2)*r0(2)+c(j,3)*r0(3)
              enddo
!     Rotation of node nk 
              do j=1,3
                co(j,nk)=costruc(j,noderel)+r(j)
              enddo
            enddo
!     
!           if no position external to the structure was found: error     
!     
            if(ier.eq.1) then
              return
            endif
          endif         
!     
        enddo
      enddo
!     
      return
      end

      
