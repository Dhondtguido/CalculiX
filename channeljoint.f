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
      subroutine channeljoint(nelem,nelup,nup,iponoel,inoel,ielprop,
     &     prop,ipkon,kon,mi,v,g,dg,nstackb,istackb,rho,xflow,co)
!
!     treats a channel joint for a frontwater calculation
!
      implicit none
!
      integer nelem,iponoel(*),inoel(2,*),nelup,nup,index,ielprop(*),
     &     indexp,indexe,ipkon(*),kon(*),nup1,nup2,mi(*),nel1,nel,
     &     nstackb,istackb(2,*)
!
      real*8 xflow1,xflow2,v(0:mi(2),*),prop(*),b1,b2,theta1,theta2,dl,
     &     co(3,*),sqrts01,sqrts02,h1,h2,hk,hnsj,rho,xflow,s0,zup,zdo,
     &     g(3),dg
!
c      nel=0
      nel1=0
!
!     loop over all elements to which node nup belongs
!
      index=iponoel(nup)
      do
!
!       treating element nelup: branch number 2
!
        if(inoel(1,index).eq.nelup) then
!
          indexe=ipkon(nelup)
          if(kon(indexe+1).eq.nup) then
            nup2=kon(indexe+3)
            xflow2=-v(1,kon(indexe+2))
          else
            nup2=kon(indexe+1)
            xflow2=v(1,kon(indexe+2))
          endif
!
          indexp=ielprop(nelup)
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
          dl=prop(indexp+3)
          if(dl.le.0.d0) then
            dl=dsqrt((co(1,nup2)-co(1,nup))**2+
     &           (co(2,nup2)-co(2,nup))**2+
     &           (co(3,nup2)-co(3,nup))**2)
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
            zdo=(-g(1)*co(1,nup)-g(2)*co(2,nup)-g(3)*co(3,nup))/dg
            s0=(zup-zdo)/dl
          endif
          sqrts02=1.d0-s0*s0
          if(sqrts02.lt.0.d0) then
            sqrts02=0.d0
          else
            sqrts02=dsqrt(sqrts02)
          endif
!
!         height at the joint as calculated from branch 2
!
          h2=v(2,nup)/sqrts02
!
          index=inoel(2,index)
          if(index.eq.0) exit
          cycle
        endif
!
!       element to be treated is not nelup
!
        nel=inoel(1,index)
        indexe=ipkon(nel)
!
!       if mass flow in branch is zero: cycle
!
        if(v(1,kon(indexe+2)).eq.0.d0) then
          nelem=nel
          index=inoel(2,index)
          if(index.eq.0) exit
          cycle
        endif
!     
!     mass flow is not zero in new branch: must be branch 1
!
        nel1=nel
        if(kon(indexe+1).eq.nup) then
          nup1=kon(indexe+3)
          xflow1=-v(1,kon(indexe+2))
        else
          nup1=kon(indexe+1)
          xflow1=v(1,kon(indexe+2))
        endif
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
        dl=prop(indexp+3)
        if(dl.le.0.d0) then
          dl=dsqrt((co(1,nup1)-co(1,nup))**2+
     &         (co(2,nup1)-co(2,nup))**2+
     &         (co(3,nup1)-co(3,nup))**2)
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
          zdo=(-g(1)*co(1,nup)-g(2)*co(2,nup)-g(3)*co(3,nup))/dg
          s0=(zup-zdo)/dl
        endif
        sqrts01=1.d0-s0*s0
        if(sqrts01.lt.0.d0) then
          sqrts01=0.d0
        else
          sqrts01=dsqrt(sqrts01)
        endif
!
!     depth at the intersection calculated from branch 1.
!     it is assumed that element nel1 is very short and that 
!     the depth in nup is the same as in nup1
!
        if(kon(indexe+1).eq.nup) then
          h1=v(2,kon(indexe+3))/sqrts01
        else
          h1=v(2,kon(indexe+1))/sqrts01
        endif
        index=inoel(2,index)
        if(index.eq.0) exit
      enddo
!
!     if no element for branch 1 was found it means that of the
!     3 intersecting branched the mass flow is known in only one of
!     them. Return and look for another branch to treat
!
      if(nel1.eq.0) then
        nelem=0
        return
      endif
!
!     the mass flow in two of the three branches in known
!
      if((h1.gt.0.d0).or.(h2.gt.0.d0)) then
!
!     if h1>h2: branch 1 is continued in forward direction
!
        if(h1.gt.h2) then
!
!       if h2<0: store info in istackb for later treatment as 
!       backwater curve with starting point: critical depth at
!       the intersection calculated within branch 2
!
          if(h2.lt.0.d0) then
            call hcrit(xflow2,rho,b2,theta2,dg,sqrts02,hk)
            v(2,nup2)=hk*sqrts02
!     
!     stackb stores the node from which to start the backwater
!     curve (in istackb(2,*) and the element downstream of this
!     node (in istackb(1,*)
!     
            nstackb=nstackb+1
            istackb(1,nstackb)=nelup
            istackb(2,nstackb)=nup2
!
!        0<=h2<h1: check for depth after jump
!
          else
            call hns(xflow2,rho,b2,theta2,dg,sqrts02,h2,hnsj)
!
!           if hns(h2)<=h1: backwater curve starting at h1 in branch2
!
            if(hnsj.le.h1) then
              v(2,nup2)=h1*sqrts02
              nstackb=nstackb+1
              istackb(1,nstackb)=nelup
              istackb(2,nstackb)=nup2
!
!           else: do nothing (frontwater curve in branch 1 and branch 2)
!
            endif
          endif
!
!         starting point of forward curve in downstream branch of joint
!
          v(2,nup)=h1*sqrts01
        else
!
!     if h1<=h2: branch 2 is continued in forward direction
!
!     if h1<0: store info in istackb for later treatment as 
!     backwater curve with starting point: critical depth at
!     the intersection calculated within branch 1
!
          if(h1.lt.0.d0) then
            call hcrit(xflow1,rho,b1,theta1,dg,sqrts01,hk)
            v(2,nup1)=hk*sqrts01
!     
!     stackb stores the node from which to start the backwater
!     curve (in istackb(2,*) and the element downstream of this
!     node (in istackb(1,*)
!     
            nstackb=nstackb+1
            istackb(1,nstackb)=nel1
            istackb(2,nstackb)=nup1
!
!        0<=h1<=h2: check for depth after jump
!
          else
            call hns(xflow1,rho,b1,theta1,dg,sqrts01,h1,hnsj)
!
!           if hns(h1)<=h2: backwater curve starting at h2 in branch 1
!
            if(hnsj.le.h2) then
              v(2,nup1)=h2*sqrts01
              nstackb=nstackb+1
              istackb(1,nstackb)=nel1
              istackb(2,nstackb)=nup1
!
!           else: do nothing (frontwater curve in branch 1 and branch 2)
!
            endif
          endif
!
!         starting point of forward curve in downstream branch of joint
!         is already correct (depth of branch 2)
!
        endif
      endif
!
!     downstream mass flow is the sum of the upstream mass flows
!
      xflow=xflow1+xflow2
!     
      return
      end
