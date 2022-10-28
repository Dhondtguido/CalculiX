!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine distattach_3d(xig,etg,zeg,pneigh,pnode,a,p,ratio,
     &          nterms)
!
!     calculates the distance between the node with coordinates
!     in "pnode" and the node with local coordinates xig and etg
!     in a face described by "nterms" nodes with coordinates
!     in pneigh
!
      implicit none
!
      integer nterms,i,j,n
!
      real*8 ratio(20),pneigh(3,*),pnode(3),a,xi,et,xig,etg,p(3),
     &  dummy,ze,zeg,omg,omh,omr,opg,oph,opr,dx(3),al,
     &  tpgphpr,tmgphpr,tmgmhpr,tpgmhpr,tpgphmr,tmgphmr,tmgmhmr,tpgmhmr,
     &  omgopg,omhoph,omropr,omgmopg,omhmoph,omrmopr
!
      n=3
!
      if(nterms.eq.4) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         ze=(zeg+1.d0)/2.d0
         dx(1)=xi
         dx(2)=et
         dx(3)=ze
         call insertsortd(dx,n)
c         call dsort(dx,iy,n,kflag)
         if(dx(3).gt.1.d-30) then
            al=dx(3)/(xi+et+ze)
            xi=al*xi
            et=al*et
            ze=al*ze
         endif
!
!        shape functions
!
         ratio(1)=1.d0-xi-et-ze
         ratio(2)=xi
         ratio(3)=et
         ratio(4)=ze
      elseif(nterms.eq.6) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
!
         ze=zeg

         a=1.d0-xi-et
!
!     shape functions
!
         ratio(1)=0.5d0*a *(1.d0-ze)
         ratio(2)=0.5d0*xi*(1.d0-ze)
         ratio(3)=0.5d0*et*(1.d0-ze)
         ratio(4)=0.5d0*a *(1.d0+ze)
         ratio(5)=0.5d0*xi*(1.d0+ze) 
         ratio(6)=0.5d0*et*(1.d0+ze) 
      elseif(nterms.eq.8) then
         xi=xig
         et=etg
         ze=zeg
!
         omg=1.d0-xi
         omh=1.d0-et
         omr=1.d0-ze
         opg=1.d0+xi
         oph=1.d0+et
         opr=1.d0+ze
!     
!        shape functions
!
         ratio(1)=omg*omh*omr/8.d0
         ratio(2)=opg*omh*omr/8.d0
         ratio(3)=opg*oph*omr/8.d0
         ratio(4)=omg*oph*omr/8.d0
         ratio(5)=omg*omh*opr/8.d0
         ratio(6)=opg*omh*opr/8.d0
         ratio(7)=opg*oph*opr/8.d0
         ratio(8)=omg*oph*opr/8.d0
      elseif(nterms.eq.10) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         ze=(zeg+1.d0)/2.d0
         dx(1)=xi
         dx(2)=et
         dx(3)=ze
         call insertsortd(dx,n)
c         call dsort(dx,iy,n,kflag)
         if(dx(3).gt.1.d-30) then
            al=dx(3)/(xi+et+ze)
            xi=al*xi
            et=al*et
            ze=al*ze
         endif
c         if(xi+et+ze.gt.1.d0) then
c            dummy=2.d0*(1.d0-xi-et-ze)/3.d0
c            xi=dummy+xi
c            et=dummy+et
c            ze=dummy+ze
c         endif
!
!        shape functions
!
         a=1.d0-xi-et-ze
         ratio( 1)=(2.d0*a-1.d0)*a
         ratio( 2)=xi*(2.d0*xi-1.d0)
         ratio( 3)=et*(2.d0*et-1.d0)
         ratio( 4)=ze*(2.d0*ze-1.d0)
         ratio( 5)=4.d0*xi*a
         ratio( 6)=4.d0*xi*et
         ratio( 7)=4.d0*et*a
         ratio( 8)=4.d0*ze*a
         ratio( 9)=4.d0*xi*ze
         ratio(10)=4.d0*et*ze
      elseif(nterms.eq.15) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
!
         ze=zeg
         a=1.d0-xi-et
!
!     shape functions
!
         ratio(1)=-0.5d0*a*(1.d0-ze)*(2.d0*xi+2.d0*et+ze)
         ratio(2)=0.5d0*xi*(1.d0-ze)*(2.d0*xi-2.d0-ze)
         ratio(3)=0.5d0*et*(1.d0-ze)*(2.d0*et-2.d0-ze)
         ratio(4)=-0.5d0*a*(1.d0+ze)*(2.d0*xi+2.d0*et-ze)
         ratio(5)=0.5d0*xi*(1.d0+ze)*(2.d0*xi-2.d0+ze)
         ratio(6)=0.5d0*et*(1.d0+ze)*(2.d0*et-2.d0+ze)
         ratio(7)=2.d0*xi*a*(1.d0-ze)
         ratio(8)=2.d0*xi*et*(1.d0-ze)
         ratio(9)=2.d0*et*a*(1.d0-ze)
         ratio(10)=2.d0*xi*a*(1.d0+ze)
         ratio(11)=2.d0*xi*et*(1.d0+ze) 
         ratio(12)=2.d0*et*a*(1.d0+ze)
         ratio(13)= a*(1.d0-ze*ze)
         ratio(14)=xi*(1.d0-ze*ze)
         ratio(15)=et*(1.d0-ze*ze)
      elseif(nterms.eq.20) then
         xi=xig
         et=etg
         ze=zeg
!
         omg=1.d0-xi
         omh=1.d0-et
         omr=1.d0-ze
         opg=1.d0+xi
         oph=1.d0+et
         opr=1.d0+ze
         tpgphpr=opg+oph+ze
         tmgphpr=omg+oph+ze
         tmgmhpr=omg+omh+ze
         tpgmhpr=opg+omh+ze
         tpgphmr=opg+oph-ze
         tmgphmr=omg+oph-ze
         tmgmhmr=omg+omh-ze
         tpgmhmr=opg+omh-ze
         omgopg=omg*opg/4.d0
         omhoph=omh*oph/4.d0
         omropr=omr*opr/4.d0
         omgmopg=(omg-opg)/4.d0
         omhmoph=(omh-oph)/4.d0
         omrmopr=(omr-opr)/4.d0
!     
!     shape functions
!
         ratio( 1)=-omg*omh*omr*tpgphpr/8.d0
         ratio( 2)=-opg*omh*omr*tmgphpr/8.d0
         ratio( 3)=-opg*oph*omr*tmgmhpr/8.d0
         ratio( 4)=-omg*oph*omr*tpgmhpr/8.d0
         ratio( 5)=-omg*omh*opr*tpgphmr/8.d0
         ratio( 6)=-opg*omh*opr*tmgphmr/8.d0
         ratio( 7)=-opg*oph*opr*tmgmhmr/8.d0
         ratio( 8)=-omg*oph*opr*tpgmhmr/8.d0
         ratio( 9)=omgopg*omh*omr
         ratio(10)=omhoph*opg*omr
         ratio(11)=omgopg*oph*omr
         ratio(12)=omhoph*omg*omr
         ratio(13)=omgopg*omh*opr
         ratio(14)=omhoph*opg*opr
         ratio(15)=omgopg*oph*opr
         ratio(16)=omhoph*omg*opr
         ratio(17)=omropr*omg*omh
         ratio(18)=omropr*opg*omh
         ratio(19)=omropr*opg*oph
         ratio(20)=omropr*omg*oph
      else
         write(*,*) '*ERROR in distattach: case with ',nterms
         write(*,*) '       terms is not covered'
         call exit(201)
      endif
!
!     calculating the position in the face
!      
      do i=1,3
         p(i)=0.d0
         do j=1,nterms
            p(i)=p(i)+ratio(j)*pneigh(i,j)
         enddo
      enddo
!
!     calculating the distance
!
      a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+(pnode(3)-p(3))**2
!
      return
      end
      
