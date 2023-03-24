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
      subroutine henergy(xflow,rho,b,theta,dg,sqrts0,e,mode,h)
!
!     determine the depth corresponding to a given value of the
!     specific energy and mass flow; either there is no solution
!     (xflow > xflowmax), one solution (xflow = xflowmax) or
!     two solutions (xflow < xflowmax), one for a frontwater curve
!     and one for a backwater curve.
!
!     INPUT:
!
!     xflow:     fluid mass flow (may be negative)
!     rho:       fluid density
!     b:         width at bottom
!     theta:     angle describing the width increase with fluid depth
!     dg:        earth acceleration
!     sqrts0:    sqrt(1.d0-s0*s0), where s0 is the change of bottom
!                height with channel length
!     e:         specific energy 
!     mode:      if 'F': frontwater curve
!                if 'B': backwater curve      
!     
!     OUTPUT:
!     
!     h:        depth: if -1 there is no solution, else it is the
!               solution for the corresponding mode
!
      implicit none
!      
      character*1 mode
!
      integer iit
!      
      real*8 xflow,rho,b,dg,sqrts0,h,theta,tth,e,cc,area,
     &     hkmax,xflowmax,discharge,c1,c2,dadh,v,f,df,dh,aa,bb
!
      tth=dtan(theta)
!
!     determine the critical depth
!
      if(theta.lt.1.d-10) then
        hkmax=2.d0*e/3.d0
      else
        aa=5.d0*tth*sqrts0
        bb=-4.d0*e*tth+3.d0*b*sqrts0
        cc=-2.d0*b*e
        hkmax=(-bb+dsqrt(bb*bb-4.d0*aa*cc))/(2.d0*aa)
      endif
!
      area=(b+hkmax*tth)*hkmax
!
      c1=2.d0*dg
      c2=c1*sqrts0
      c1=c1*e
!
      xflowmax=rho*area*dsqrt(c1-hkmax*c2)
!
!     if xflow > xflowmax: no solution
!
      if(xflow.gt.xflowmax) then
        h=-1.d0
        return
      endif
!
!     starting value for the iterations: in the middle between
!     hkmax and the points for zero flow (0 or e/sqrts0)
!
      if(mode.eq.'F') then
        h=hkmax/2.d0
      else
        h=(hkmax+e/sqrts0)/2.d0
      endif
!
      discharge=xflow/rho
      iit=0
!
!     solving area*dsqrt(2.d0*dg*(e-h*sqrts0))-discharge=0
!     with the Newton-Raphson method
!
      do
        area=h*(b+h*tth)
        dadh=b+2.d0*h*tth
        v=dsqrt(c1-h*c2)
        f=area*v-discharge
        df=dadh*v-c2*area/(2.d0*v)
        dh=-f/df
        if((dabs(dh).lt.1.d-10).or.(dabs(dh).lt.1.d-5*h)) exit
        h=h+dh
        h=max(0.d0,h)
        if(h.lt.0.d0) h=(h-dh)/2.d0
        if(h.ge.e/sqrts0) h=(h-dh+e/sqrts0)/2.d0
        iit=iit+1
        if(iit.gt.100) then
          write(*,*) '*ERROR in henergy: too many iterations'
          call exit(201)
        endif
      enddo
!
      return
      end
