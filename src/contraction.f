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
      subroutine contraction(nelem,ielprop,prop,nup,nmid,ndo,dg,
     &     mode,xflow,rho,nelup,neldo,istack,nstack,
     &     mi,v,inv,epsilon,co)
!
!     treats the channel elements CONTRACTION, ENLARGEMENT, STEP
!     and DROP
!
      implicit none
!
      character*1 mode
!
      integer nelem,ielprop(*),index,nup,ndo,nelup,neldo,nstack,
     &     istack(2,*),mi(*),inv,nmid
!
      real*8 prop(*),s0,sqrts0,hdo,hup,xflow,dg,rho,hk,v(0:mi(2),*),
     &     areado,areaup,b1,b2,bdo,bup,d,e,theta1,theta2,thetado,
     &     thetaup,tthdo,tthup,epsilon,co(3,*),pi,dgmod,dl,xk,alpha
!
      pi=4.d0*datan(1.d0)
!
!     determining the properties
!
      index=ielprop(nelem)
!
!     width of the channel at node 1
!     
      b1=prop(index+1)
!
!     trapezoidal angle of the channel at node 1
!
      theta1=prop(index+2)
!
!     width of the channel at the other end node (node 3)
!     
      b2=prop(index+3)
!
!     trapezoidal angle of the channel at the other end node (node 3)
!
      theta2=prop(index+4)
!
!     size of the step going from node 1 to node 3 (a negative value
!     corresponds to a drop of the channel bottom)
!
      d=prop(index+5)
      if((d.ne.0.d0).and.
     &     (((b2.ne.0.d0).and.(b1.ne.b2)).or.
     &     (((theta2.ne.0.d0).and.(theta2.ne.theta1))))) then
        write(*,*) '*ERROR in contraction'
        write(*,*) '       step height is nonzero and'
        write(*,*) '       cross section is changing at the'
        write(*,*) '       same time; this is not allowed'
        call exit(201)
      endif
!
!     step/drop: section does not change
!
      if(d.ne.0.d0) then
        b2=b1
        theta2=theta1
      endif
!
!     if the length of the element is negative, it is determined from
!     the coordinates
!
      dl=prop(index+6)
      if(dl.le.0.d0) then
        dl=dsqrt((co(1,nup)-co(1,ndo))**2+
     &           (co(2,nup)-co(2,ndo))**2+
     &       (co(3,nup)-co(3,ndo))**2)
      endif
!
!     s0: sine of downstream slope (the slope is the angle phi between the
!         channel bottom and a plane orthogonal to the gravity vector
!     sqrts0: cosine of downstream slope
!     needed to calculate the normal depth (he)      
!
      s0=prop(index+7)
      if(s0.lt.-1.d0) then
        write(*,*) '*ERROR in contraction: sine of slope'
        write(*,*) '       must be given explicitly'
        write(*,*) '       for a contraction, enlargement,'
        write(*,*) '       step or drop'
        call exit(201)
      endif
      sqrts0=1.d0-s0*s0
      if(sqrts0.lt.0.d0) then
        sqrts0=0.d0
      else
        sqrts0=dsqrt(sqrts0)
      endif
!
!     check the direction of the flow
!
      if(inv.eq.1) then
        bup=b1
        thetaup=theta1
        bdo=b2
        thetado=theta2
      else
        bup=b2
        thetaup=theta2
        bdo=b1
        thetado=theta1
        d=-d
      endif
!
!     calculating the contraction/expansion angle
!
      if(d.eq.0.d0) then
        if(dl.eq.0.d0) then
          if(bdo.gt.bup) then
            alpha=pi/2.d0
          elseif(bdo.lt.bup) then
            alpha=-pi/2.d0
          else
            alpha=0.d0
          endif
        else
          alpha=datan((bdo-bup)/(2.d0*dl))
        endif
      endif
!
!     modifying g by the head loss coefficient
!
      dgmod=dg
!
!     head loss coefficient: contraction
!
      if(d.eq.0.d0) then
        if(alpha.le.0.d0) then
          if(mode.eq.'F') then
            dgmod=pi*dg/(pi+alpha)
          else
            dgmod=pi*dg/(pi-alpha)
          endif
        endif
!     
!       head loss coefficient: enlargement
!     
        if(alpha.gt.0.d0) then
          if(alpha.ge.0.79d0) then
            xk=0.87d0
          elseif(alpha.ge.0.46d0) then
            xk=0.68d0+0.5757d0*(alpha-0.46d0)
          elseif(alpha.ge.0.32d0) then
            xk=0.41d0+1.9286d0*(alpha-0.32d0)
          elseif(alpha.ge.0.25d0) then
            xk=0.27d0+2.d0*(alpha-0.25d0)
          else
            xk=0.27d0*alpha/0.25d0
          endif
          if(mode.eq.'F') then
            dgmod=dg/(1.d0+xk)
          else
            dgmod=dg/(1.d0-xk)
          endif
        endif
      endif
!***
c      dgmod=dg
!***
!
      v(1,nmid)=inv*xflow
!
      if(mode.eq.'F') then
!
!        frontwater curve
!
        tthup=dtan(thetaup)
        hup=v(2,nup)/sqrts0
        if(hup.le.0.d0)then
!
!         take the critical depth upstream
!
          call hcrit(xflow,rho,bup,thetaup,dg,sqrts0,hk)
          areaup=(bup+hk*tthup)*hk
          e=(xflow/(areaup*rho))**2/(2.d0*dgmod)+(hk-d)*sqrts0
        else
          areaup=(bup+hup*tthup)*hup
          e=(xflow/(areaup*rho))**2/(2.d0*dgmod)+(hup-d)*sqrts0
        endif
!
!       calculate the downstream height
!
        call henergy(xflow,rho,bdo,thetado,dgmod,sqrts0,e,mode,hdo)
!
        if(hdo.gt.0.d0) then
!
!         first calculate the backwater curve starting in nup
!
          if(hup.le.0.d0) then
            v(2,nup)=hk*sqrts0
            ndo=nup
            nelem=nelup
            mode='B'
            nstack=nstack+1
            istack(1,nstack)=nelup
            istack(2,nstack)=nup
            return
          endif
!
          v(2,ndo)=hdo*sqrts0
!     
!         calculate the critical depth for output purposes
!     
          call hcrit(xflow,rho,bup,thetaup,dg,sqrts0,hk)
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
          call hcrit(xflow,rho,bdo,thetado,dg,sqrts0,hk)
          v(3,ndo)=hk
!
          v(2,ndo)=hk*sqrts0
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
        hdo=v(2,ndo)/sqrts0
        tthdo=dtan(thetado)
        areado=(bdo+hdo*tthdo)*hdo
        e=(xflow/(areado*rho))**2/(2.d0*dgmod)+(hdo+d)*sqrts0
!
!       calculate the upstream height
!
        call henergy(xflow,rho,bup,thetaup,dgmod,sqrts0,e,mode,hup)
!
        if(hup.gt.0.d0) then
          v(2,nup)=hup*sqrts0
!     
!         calculate the critical depth for output purposes
!     
          call hcrit(xflow,rho,bdo,thetado,dg,sqrts0,hk)
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
          call hcrit(xflow,rho,bup,thetaup,dg,sqrts0,hk)
          v(3,nup)=hk
!          
          v(2,nup)=hk*sqrts0
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
