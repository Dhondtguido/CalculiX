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
!     in channels with a non-erosive bottom
!     
      subroutine channeljointback(neldo,ndo,iponoel,inoel,ipkon,kon,
     &     mi,v,istackb,nstackb,co,ielprop,prop,g,dg,xflow,rho)
!
!     treats a channel joint for a backwater calculation
!
      implicit none
!
      logical stacked
!
      integer nel1,nel2,nup1,nup2,iponoel(*),inoel(2,*),ipkon(*),kon(*),
     &     index,inv,mi(*),i,istackb(2,*),nstackb,nel1sav,nup1sav,
     &     ndoo,indexp,ielprop(*),niter,neldo,ndo
!
      real*8 v(0:mi(2),*),xflow1,xflow2,xflow1sav,vec(3),vec1(3),dl,dl1,
     &     vec2(3),co(3,*),alpha1,alpha2,pi,prop(*),g(3),dg,h,u,dl2,
     &     up2,area,u1,u1p2,area1,u2,u2p2,area2,h1new,h2new,theta,h1,h2,
     &     theta1,theta2,xflow,zdo,zup,sqrts0,sqrts01,sqrts02,s0,rho,b,
     &     b1,b2
!
      nel1=0
      nel2=0
!
      index=iponoel(ndo)
!
!     find the two upstream elements nel1 and nel2 of neldo and
!     the upstream nodes nup1 and nup2 of these elements
!
      do
        if(inoel(1,index).ne.neldo) then
          if(nel1.eq.0) then
            nel1=inoel(1,index)
            if(kon(ipkon(nel1)+1).eq.ndo) then
              nup1=kon(ipkon(nel1)+3)
              inv=-1
            else
              nup1=kon(ipkon(nel1)+1)
              inv=1
            endif
            xflow1=inv*v(1,kon(ipkon(nel1)+2))
          else
            nel2=inoel(1,index)
            if(kon(ipkon(nel2)+1).eq.ndo) then
              nup2=kon(ipkon(nel2)+3)
              inv=-1
            else
              nup2=kon(ipkon(nel2)+1)
              inv=1
            endif
            xflow2=inv*v(1,kon(ipkon(nel2)+2))
          endif
        endif
        index=inoel(2,index)
        if(index.eq.0) exit
      enddo
!
!     check whether nel1 or nel2 is already on istackb;
!     1) if nel1 is on istackb, do nothing
!     2) if nel2 is on istackb, switch nel1 with nel2
!     3) if none is on istack, put nel1 on istackb
!
      stacked=.false.
      do i=1,nstackb
        if(istackb(1,i).eq.nel1) then
          stacked=.true.
        elseif(istackb(1,i).eq.nel2) then
!
!         switch nel1 and nel2
!
          nel1sav=nel1
          nup1sav=nup1
          xflow1sav=xflow1
!
          nel1=nel2
          nup1=nup2
          xflow1=xflow2
!
          nel2=nel1sav
          nup2=nup1sav
          xflow2=xflow1sav
!
          stacked=.true.
!
        endif
      enddo
!
      if(.not.stacked) then
        nstackb=nstackb+1
        istackb(1,nstackb)=nel1
        istackb(2,nstackb)=nup1
      endif
!
!     loop for the downstream node of element neldo
!
      if(kon(ipkon(neldo)+1).eq.ndo) then
        ndoo=kon(ipkon(neldo)+3)
      else
        ndoo=kon(ipkon(neldo)+1)
      endif
!
!     determine the properties of neldo
!
      indexp=ielprop(neldo)
!     
!     width of the channel at zero depth
!     
      b=prop(indexp+1)
!     
!     trapezoidal angle of the channel cross section
!     
      theta=prop(indexp+2)
!     
!     if the length of the element is negative, it is determined from
!     the coordinates
!     
      dl=prop(indexp+3)
      if(dl.le.0.d0) then
        dl=dsqrt((co(1,ndo)-co(1,ndoo))**2+
     &       (co(2,ndo)-co(2,ndoo))**2+
     &       (co(3,ndo)-co(3,ndoo))**2)
      endif
!     
!     s0: sine of slope (the slope is the angle phi between the channel
!     bottom and a plane orthogonal to the gravity vector
!     sqrts0: cosine of slope
!     
!     determining the sine of the slope; if the sine is less than
!     -1.d0 it is calculated from the coordinates
!     
      s0=prop(indexp+4)
      if(s0.lt.-1.d0) then
        zup=(-g(1)*co(1,ndo)-g(2)*co(2,ndo)-g(3)*co(3,ndo))/dg
        zdo=(-g(1)*co(1,ndoo)-g(2)*co(2,ndoo)-g(3)*co(3,ndoo))/dg
        s0=(zup-zdo)/dl
      endif
      sqrts0=1.d0-s0*s0
      if(sqrts0.lt.0.d0) then
        sqrts0=0.d0
      else
        sqrts0=dsqrt(sqrts0)
      endif
!
!     determine the properties of nel1
!
      indexp=ielprop(nel1)
!     
!     width of the channel at zero depth
!     
      b1=prop(indexp+1)
!     
!     trapezoidal angle of the channel cross section
!     
      theta1=prop(indexp+2)
!     
!     if the length of the element is negative, it is determined from
!     the coordinates
!     
      dl1=prop(indexp+3)
      if(dl1.le.0.d0) then
        dl1=dsqrt((co(1,nup1)-co(1,ndo))**2+
     &       (co(2,nup1)-co(2,ndo))**2+
     &       (co(3,nup1)-co(3,ndo))**2)
      endif
!     
!     s0: sine of slope (the slope is the angle phi between the channel
!     bottom and a plane orthogonal to the gravity vector
!     sqrts0: cosine of slope
!     
!     determining the sine of the slope; if the sine is less than
!     -1.d0 it is calculated from the coordinates
!     
      s0=prop(indexp+4)
      if(s0.lt.-1.d0) then
        zup=(-g(1)*co(1,nup1)-g(2)*co(2,nup1)-g(3)*co(3,nup1))/dg
        zdo=(-g(1)*co(1,ndo)-g(2)*co(2,ndo)-g(3)*co(3,ndo))/dg
        s0=(zup-zdo)/dl1
      endif
      sqrts01=1.d0-s0*s0
      if(sqrts01.lt.0.d0) then
        sqrts01=0.d0
      else
        sqrts01=dsqrt(sqrts01)
      endif
!
!     angle between nel1 and neldo
!
      alpha1=prop(indexp+6)
!
!     determine the properties of nel2
!
      indexp=ielprop(nel2)
!     
!     width of the channel at zero depth
!     
      b2=prop(indexp+1)
!     
!     trapezoidal angle of the channel cross section
!     
      theta2=prop(indexp+2)
!     
!     if the length of the element is negative, it is determined from
!     the coordinates
!     
      dl2=prop(indexp+3)
      if(dl2.le.0.d0) then
        dl2=dsqrt((co(1,nup2)-co(1,ndo))**2+
     &       (co(2,nup2)-co(2,ndo))**2+
     &       (co(3,nup2)-co(3,ndo))**2)
      endif
!     
!     s0: sine of slope (the slope is the angle phi between the channel
!     bottom and a plane orthogonal to the gravity vector
!     sqrts0: cosine of slope
!     
!     determining the sine of the slope; if the sine is less than
!     -1.d0 it is calculated from the coordinates
!     
      s0=prop(indexp+4)
      if(s0.lt.-1.d0) then
        zup=(-g(1)*co(1,nup2)-g(2)*co(2,nup2)-g(3)*co(3,nup2))/dg
        zdo=(-g(1)*co(1,ndo)-g(2)*co(2,ndo)-g(3)*co(3,ndo))/dg
        s0=(zup-zdo)/dl2
      endif
      sqrts02=1.d0-s0*s0
      if(sqrts02.lt.0.d0) then
        sqrts02=0.d0
      else
        sqrts02=dsqrt(sqrts02)
      endif
!
!     angle between nel2 and neldo
!
      alpha2=prop(indexp+6)
!
!     if alpha1<0 or alpha2<0 it means that the angle at stake has
!     to be calculated from the coordinates
!
!     calculating the normalized vectors connecting ndo0 with ndo (vec),
!     ndo with nup1 (vec1) and ndo with nup2 (vec2)
!
      if((alpha1.lt.0.d0).or.(alpha2.lt.0.d0)) then
        do i=1,3
          vec(i)=(co(i,ndo)-co(i,ndoo))/dl
        enddo
        pi=4.d0*datan(1.d0)
      endif
!     
      if(alpha1.lt.0.d0) then
        do i=1,3
          vec1(i)=(co(i,nup1)-co(i,ndo))/dl1
        enddo
!     
!     determine the angle alpha1 between vec and vec1 and the angle
!     alpha2 between vec and vec2
!     
        alpha1=vec(1)*vec1(1)+vec(2)*vec1(2)+vec(3)*vec1(3)
        if(alpha1.gt.1.d0) then
          alpha1=0.d0
        elseif(alpha1.lt.-1.d0) then
          alpha1=pi
        else
          alpha1=dacos(alpha1)
        endif
      endif
!     
      if(alpha2.lt.0.d0) then
        do i=1,3
          vec2(i)=(co(i,nup2)-co(i,ndo))/dl2
        enddo
!     
!     determine the angle alpha1 between vec and vec1 and the angle
!     alpha2 between vec and vec2
!     
        alpha2=vec(1)*vec2(1)+vec(2)*vec2(2)+vec(3)*vec2(3)
        if(alpha2.gt.1.d0) then
          alpha2=0.d0
        elseif(alpha2.lt.-1.d0) then
          alpha2=pi
        else
          alpha2=dacos(alpha2)
        endif
      endif
!
!     determine the loss coefficient from the angle
!
      alpha1=alpha1/pi
      alpha2=alpha2/pi
!
!     depth in ndo and velocity in neldo
!
      h=v(2,ndo)/sqrts0
      area=h*(b+h*dtan(theta))
      u=xflow/(rho*area)
      up2=u*u
!
!     initial guess for u1^2 and u2^2
!
      u1p2=up2*xflow1/xflow
      u2p2=up2*xflow2/xflow
!
!     calculation of h1 and h2 from the energy equation
!
      h1=h+(up2+alpha1*dabs(u1p2-up2)-u1p2)/(2.d0*dg)
      h2=h+(up2+alpha2*dabs(u2p2-up2)-u2p2)/(2.d0*dg)
!
      niter=0
!
!     iterative loop to find h1 and h2
!
      do
        area1=h1*(b1+h1*dtan(theta1))
        u1=xflow1/(rho*area1)
        u1p2=u1*u1
!        
        area2=h2*(b2+h2*dtan(theta2))
        u2=xflow2/(rho*area2)
        u2p2=u2*u2
!
        h1new=h+(up2+alpha1*dabs(u1p2-up2)-u1p2)/(2.d0*dg)
        h2new=h+(up2+alpha2*dabs(u2p2-up2)-u2p2)/(2.d0*dg)
!
!       convergence check
!
        if(((dabs(h1-h1new).lt.1.d-3).or.(dabs(h1-h1new).lt.1.d-3*h1))
     &       .and.
     &     ((dabs(h2-h2new).lt.1.d-3).or.(dabs(h2-h2new).lt.1.d-3*h2)))
     &       exit
!
        h1=h1new
        h2=h2new
        niter=niter+1
        if(niter.gt.100) then
          write(*,*) '*ERROR in channeljointback; more than 100'
          write(*,*) '       iterations: stop'
          call exit(201)
        endif
      enddo
!
      v(2,nup1)=h1*sqrts01
      v(2,nup2)=h2*sqrts02
      neldo=nel2
      ndo=nup2
!
      return
      end
