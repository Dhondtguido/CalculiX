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
!     in channels with a non-erosive bottom: reservoir
!     
      subroutine reservoir(nelem,ielprop,prop,nup,nmid,ndo,co,g,
     &     dg,mode,xflow,rho,dvi,nelup,mi,v,inv,epsilon,istack,
     &     nstack)
!
      implicit none
!
      character*1 mode
!
      integer nelem,ielprop(*),index,nup,ndo,nmid,nelup,inv,
     &     mi(*),istack(2,*),nstack,ialreadystacked,i
!
      real*8 prop(*),co(3,*),g(3),b,theta,dl,s0,sqrts0,xks,hdo,
     &     hup,xflow,dg,rho,hk,tth,reynolds,dvi,friction,hd,
     &     form_fact,he,v(0:mi(2),*),zdo,zup,epsilon,hr
!
!     determining the properties
!
      index=ielprop(nelem)
!
!     width of the channel at zero depth
!
      b=prop(index+1)
!
!     trapezoidal angle of the channel cross section
!
      theta=prop(index+2)
      tth=dtan(theta)
!
!     length of the element (actually not needed)
!
      dl=prop(index+3)
!
!     s0: sine of slope (the slope is the angle phi between the channel
!         bottom and a plane orthogonal to the gravity vector
!     sqrts0: cosine of slope
!      
      s0=prop(index+4)
      if(s0.lt.-1.d0) then
        write(*,*) '*ERROR in reservoir: sine of slope'
        write(*,*) '       must be given explicitly'
        write(*,*) '       for sluice gate or wear'
        call exit(201)
      endif
      sqrts0=1.d0-s0*s0
      if(sqrts0.lt.0.d0) then
        sqrts0=0.d0
      else
        sqrts0=dsqrt(sqrts0)
      endif
!
!     grain size (if positive: White-Colebrook friction, if negative:
!     Manning)
!
      xks=prop(index+5)
!
!     depth of the reservoir
!
      hr=v(2,ndo)
!
      hup=v(2,nup)/sqrts0
      v(1,nmid)=inv*xflow
!
      if(hup.le.0.d0) then
!
!       no forward solution: force backward solution
!
c        hdo=0.d0
        hdo=-1.d0
      else
        call hns(xflow,rho,b,theta,dg,sqrts0,hup,hdo)
      endif
!
      if(hdo.gt.hr) then
!
!        fall or jump in reservoir
!
        v(2,ndo)=hr
!     
!       calculate the critical depth for output purposes
!     
        call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
        v(3,nup)=hk
!     
        nelup=nelem
        nelem=0
        nup=ndo
        return
      else
!     
!       calculate the critical depth
!     
        call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
        v(3,nup)=hk
!     
!       calculate the normal depth
!     
        if(xks.gt.0.d0) then
          reynolds=xflow/(b*dvi)
          form_fact=1.d0
          hd=4.d0*hr
          call friction_coefficient(dl,hd,xks,reynolds,form_fact,
     &         friction)
        endif
        call hnorm(xflow,rho,b,theta,dg,s0,friction,xks,he)
!     
        if(hr.gt.hk) then
!
!       backwater curve with raccordation
!
          v(2,ndo)=hr
          v(2,nup)=hr*sqrts0
!          
        elseif(hk.lt.he) then
!
!       backwater curve starting at hk (weak slope)
!
          v(2,ndo)=hr
c          v(2,nup)=(hr+epsilon)*sqrts0
c          v(2,nup)=hr*sqrts0
          v(2,nup)=hk*sqrts0
!
        else
          write(*,*) '*ERROR in reservoir: strong slope at element',
     &         nelem
          write(*,*) '       backwater curve starting at the critical'
          write(*,*) '       is not feasible'
          call exit(201)
        endif
!
!       putting the actual element and downstream node on stack
!
        ialreadystacked=0
        do i=1,nstack
          if(istack(1,i).eq.nelem) then
            ialreadystacked=1
            exit
          endif
        enddo
        if(ialreadystacked.eq.0) then  
          nstack=nstack+1
          istack(1,nstack)=nelem
          istack(2,nstack)=ndo
        endif
!
!       starting the calculation of a backwater curve
!
        nelem=nelup
        ndo=nup
        mode='B'
      endif
!
      return
      end
      
