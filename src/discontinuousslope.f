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
!     Solve the Bresse equation for the turbulent stationary flow
!     in channels with a non-erosive bottom: sluice gate
!     
      subroutine discontinuousslope(nelem,ielprop,prop,nup,nmid,ndo,dg,
     &     mode,xflow,rho,nelup,neldo,istack,nstack,
     &     mi,v,inv,epsilon,co)
!
!     treats the channel element DISCONTINUOUS SLOPE
!
      implicit none
!
      character*1 mode
!
      integer nelem,ielprop(*),index,nup,ndo,nelup,neldo,nstack,
     &     istack(2,*),mi(*),inv,nmid
!
      real*8 prop(*),s01,s02,sqrts01,sqrts02,hdo,hup,xflow,dg,rho,hk,
     &     v(0:mi(2),*),areado,areaup,b,e,theta,
     &     tth,epsilon,co(3,*),sqrts0up,sqrts0do
!
!     determining the properties
!
      index=ielprop(nelem)
!
!     width of the channel at node 1
!     
      b=prop(index+1)
!
!     trapezoidal angle of the channel at node 1
!
      theta=prop(index+2)
!
!     s01: sine of slope in first node of element
!          (the slope is the angle phi between the
!           channel bottom and a plane orthogonal to the gravity vector)
!     sqrts01: cosine of slope in first node of element needed 
!              calculate the normal depth (he)      
!
      s01=prop(index+4)
      if(s01.lt.-1.d0) then
        write(*,*) '*ERROR in contraction: sine of slope'
        write(*,*) '       must be given explicitly'
        write(*,*) '       for a contraction, enlargement,'
        write(*,*) '       step or drop'
        call exit(201)
      endif
      sqrts01=1.d0-s01*s01
      if(sqrts01.lt.0.d0) then
        sqrts01=0.d0
      else
        sqrts01=dsqrt(sqrts01)
      endif
!
!     s02: sine of slope in third node of element
!          (the slope is the angle phi between the
!           channel bottom and a plane orthogonal to the gravity vector)
!     sqrts02: cosine of slope in third node of element needed 
!              calculate the normal depth (he)      
!
      s02=prop(index+5)
      if(s02.lt.-1.d0) then
        write(*,*) '*ERROR in contraction: sine of slope'
        write(*,*) '       must be given explicitly'
        write(*,*) '       for a contraction, enlargement,'
        write(*,*) '       step or drop'
        call exit(202)
      endif
      sqrts02=1.d0-s02*s02
      if(sqrts02.lt.0.d0) then
        sqrts02=0.d0
      else
        sqrts02=dsqrt(sqrts02)
      endif
!
!     check the direction of the flow
!
      if(inv.eq.1) then
        sqrts0up=sqrts01
        sqrts0do=sqrts02
      else
        sqrts0up=sqrts02
        sqrts0do=sqrts01
      endif
!
      v(1,nmid)=inv*xflow
!
      if(mode.eq.'F') then
!
!        frontwater curve
!
        tth=dtan(theta)
        hup=v(2,nup)
        if(hup.le.0.d0)then
!
!         take the critical depth upstream
!
          call hcrit(xflow,rho,b,theta,dg,sqrts0up,hk)
          areaup=(b+hk*tth)*hk
          e=(xflow/(areaup*rho))**2/(2.d0*dg)+hk*sqrts0up
        else
          areaup=(b+hup*tth)*hup
          e=(xflow/(areaup*rho))**2/(2.d0*dg)+hup*sqrts0up
        endif
!
!       calculate the downstream height
!
        call henergy(xflow,rho,b,theta,dg,sqrts0do,e,mode,hdo)
!
        if(hdo.gt.0.d0) then
!
!         first calculate the backwater curve starting in nup
!
          if(hup.le.0.d0) then
            v(2,nup)=hk
            ndo=nup
            nelem=nelup
            mode='B'
            nstack=nstack+1
            istack(1,nstack)=nelup
            istack(2,nstack)=nup
            return
          endif
!
          v(2,ndo)=hdo
!     
!         calculate the critical depth for output purposes
!     
          call hcrit(xflow,rho,b,theta,dg,sqrts0up,hk)
          v(3,nup)=hk
!
          nelup=nelem
          nelem=0
          nup=ndo
        else
!
!         no solution, raise downstream to the critical height corresponding to
!         the fluid flow
!     
          call hcrit(xflow,rho,b,theta,dg,sqrts0do,hk)
          v(3,ndo)=hk
!
          v(2,ndo)=hk
!
!         store the actual element and downstream node as start of a
!         frontwater curve          
!
          nstack=nstack+1
          istack(1,nstack)=nelem
          istack(2,nstack)=ndo
!
!         repeat the calculation of the actual element with the
!         downstream node as the start of a backwater curve
!
          mode='B'
        endif
      else
!
!       backwater curve
!
        hdo=v(2,ndo)
        tth=dtan(theta)
        areado=(b+hdo*tth)*hdo
        e=(xflow/(areado*rho))**2/(2.d0*dg)+hdo*sqrts0do
!
!       calculate the upstream height
!
        call henergy(xflow,rho,b,theta,dg,sqrts0up,e,mode,hup)
!
        if(hup.gt.0.d0) then
          v(2,nup)=hup
!     
!         calculate the critical depth for output purposes
!     
          call hcrit(xflow,rho,b,theta,dg,sqrts0do,hk)
          v(3,ndo)=hk
!          
          ndo=nup
          neldo=nelem
          nelem=0
        else
!
!         no solution, drop upstream to the critical height corresponding to
!         the fluid flow
!     
          call hcrit(xflow,rho,b,theta,dg,sqrts0up,hk)
          v(3,nup)=hk
!          
          v(2,nup)=hk
!          
          nstack=nstack+1
          istack(1,nstack)=nelup
          istack(2,nstack)=nup
!
          ndo=nup
          nelem=nelup
          neldo=nelem
        endif
      endif
!
      return
      end
