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
      subroutine sluicegate(nelem,ielprop,prop,nup,nmid,ndo,co,g,dg,
     &     mode,xflow,rho,dvi,nelup,neldo,istack,nstack,ikboun,nboun,
     &     mi,v,ipkon,kon,inv,epsilon,lakon)
!
!     treats the channel elements SLUICE GATE and WEAR
!
      implicit none
!
      character*1 mode
      character*8 lakon(*)
!
      integer nelem,ielprop(*),index,nup,ndo,nmid,nelup,neldo,nstack,
     &     istack(2,*),ikboun(*),nboun,mi(*),id,idof,ipkon(*),kon(*),
     &     inv
!
      real*8 prop(*),co(3,*),g(3),b,theta,dl,s0,sqrts0,xks,hdo,hw,h2,
     &     ha,hup,xflow,dg,rho,hk,tth,area,reynolds,dvi,friction,hkmax,
     &     aa,bb,cc,form_fact,he,v(0:mi(2),*),zdo,zup,epsilon,xflowcor,
     &     cd,hd
!
!     determining the properties
!
      index=ielprop(nelem)
!
!     width of the downstream channel at zero depth
!     needed to calculate the critical depth (hk)
!     
      b=prop(index+1)
!
!     trapezoidal angle of the downstream channel cross section
!     needed to calculate the critical depth (hk)
!
      theta=prop(index+2)
      tth=dtan(theta)
!
!     length of the downstream element (actually not needed for this
!     element)
!
      dl=prop(index+3)
!
!     s0: sine of downstream slope (the slope is the angle phi between the
!         channel bottom and a plane orthogonal to the gravity vector
!     sqrts0: cosine of downstream slope
!     needed to calculate the normal depth (he)      
!
      s0=prop(index+4)
      if(s0.lt.-1.d0) then
        write(*,*) '*ERROR in sluicegatewear: sine of slope'
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
!     grain size of downstream channel (if positive: White-Colebrook friction,
!     if negative: Manning); needed to calculate the normal depth (he)
!
      xks=prop(index+5)
!
      if(lakon(nelem)(6:7).eq.'SG') then
!     
!       depth underneath the sluice gate
!     
        ha=prop(index+6)
        cd=1.d0
        hw=0.d0
      elseif(lakon(nelem)(6:7).eq.'WE') then
!     
!       height of the wear crest
!     
        hw=prop(index+6)
        cd=prop(index+7)*(1.5d0)**(1.5d0)/dsqrt(dg)
        ha=1.d30
      endif
!
      hup=v(2,nup)
!
!     forward mode (frontwater curve)
!
      if(mode.eq.'F') then
!
!       frontwater curve
!
        if(hup.eq.0.d0) then
!
!         upstream height not known: flow must be known
!        
          v(1,nmid)=inv*xflow
!
!         determine the critical depth
!
          call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
          v(3,ndo)=hk
!
          if(ha.lt.hk) then
!
!           A3 or B2 or B3
!
            area=(b+ha*tth)*ha
            v(2,ndo)=ha*sqrts0
            v(2,nup)=(xflow/(area*rho))**2/(2.d0*dg)+ha
            nelup=nelem
            nelem=0
            nup=ndo
          else
!
!           depth underneath sluice gate exceeds critical depth
!
!           calculate the normal depth
!
            if(xks.gt.0.d0) then
              reynolds=xflow/(b*dvi)
              form_fact=1.d0
              hd=4.d0*hk
               call friction_coefficient(dl,hd,xks,reynolds,form_fact,
     &               friction)
            endif
            call hnorm(xflow,rho,b,theta,dg,s0,friction,xks,he)
!
            if(he.lt.hk) then
!
!             B2
!
              area=(b+hk*tth)*hk
c              v(2,ndo)=(hk-epsilon)*sqrts0
              v(2,ndo)=hk*sqrts0
              v(2,nup)=(xflow/(cd*area*rho))**2/(2.d0*dg)+hk*sqrts0+hw
              nelup=nelem
              nelem=0
              nup=ndo
            else
!
!             no frontwater solution
!
              v(2,ndo)=-1.d0
              nelup=nelem
              nelem=0
              nup=ndo
            endif
          endif
        elseif(((hup-hw).gt.0.d0).and.(xflow.eq.0.d0)) then
!
!         upstream depth known: flow not known (start of a new branch)
!
!         calculate hk(Qmax) and compare with ha
!
          if(theta.lt.1.d-10) then
!
!           rectangular cross section
!
            hkmax=2.d0*(hup-hw)/3.d0
          else
!
!           trapezoidal cross section
!
            aa=5.d0*tth*sqrts0
            bb=-4.d0*(hup-hw)*tth+3.d0*b*sqrts0
            cc=-2.d0*b*(hup-hw)
            hkmax=(-bb+dsqrt(bb*bb-4.d0*aa*cc))/(2.d0*aa)
          endif
!
          if(ha.lt.hkmax) then
!
!           A3 or B2 or B3
!
            area=(b+ha*tth)*ha
            xflow=area*dsqrt(2.d0*dg*(hup-ha*sqrts0))*rho
            v(1,nmid)=inv*xflow
            if(kon(ipkon(nelup)+1).eq.0) then
              v(1,kon(ipkon(nelup)+2))=xflow
            else
              v(1,kon(ipkon(nelup)+2))=-xflow
            endif
            v(2,ndo)=ha*sqrts0
            nelup=nelem
            nelem=0
            nup=ndo
          else
!
!           calculate maximal flow
!
            area=(b+hkmax*tth)*hkmax
            xflow=rho*cd*area*dsqrt(2.d0*dg*(hup-hw-hkmax))
!
!           calculate the normal depth
!
            if(xks.gt.0.d0) then
               reynolds=xflow/(b*dvi)
               form_fact=1.d0
               hd=4.d0*hkmax
               call friction_coefficient(dl,hd,xks,reynolds,
     &               form_fact,friction)
            endif
            call hnorm(xflow,rho,b,theta,dg,s0,friction,xks,he)
!
            if(he.lt.hkmax) then
!
!             B2
!
c              v(2,ndo)=(hkmax-epsilon)*sqrts0
              v(2,ndo)=hkmax*sqrts0
              v(1,nmid)=inv*xflow
              if(kon(ipkon(nelup)+1).eq.0) then
                v(1,kon(ipkon(nelup)+2))=xflow
              else
                v(1,kon(ipkon(nelup)+2))=-xflow
              endif
              nelup=nelem
              nelem=0
              nup=ndo
            else
!
!             no frontwater solution
!
              v(2,ndo)=-1.d0
              nelup=nelem
              nelem=0
              nup=ndo
            endif
          endif
        else
!
!         upstream depth known and flow known: sluice gate in the
!         middle of a channel: fluid depth upstream of gate is
!         recalculated; a wear in the middle of a channel is
!         simulated by a step and straight channel
!
          if(lakon(nelem)(6:7).eq.'WE') then
            write(*,*) '*ERROR in sluicegatewear: WEAR element'
            write(*,*) '       must on one side be connected'
            write(*,*) '       to exactly one CHANNEL INOUT'
            write(*,*) '       element; faulty element:',nelem
            call exit(201)
          endif
!          
          v(1,nmid)=inv*xflow
!
c          if((hup.lt.0.d0).or.(hup.gt.ha)) then
          if(hup.gt.ha) then
!     
!           determine the critical depth
!     
            call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
            v(3,ndo)=hk
!     
c            if(ha.lt.hk) then
!     
!     A3 or B2 or B3
!     
            area=(b+ha*tth)*ha
            v(2,ndo)=ha*sqrts0
            mode='B'
            nstack=nstack+1
            istack(1,nstack)=nelem
            istack(2,nstack)=ndo
c            else
c!     
c!             no disturbance of the flow
c!     
c              v(2,ndo)=v(2,nup)
c              v(1,nmid)=inv*xflow
c              nelup=nelem
c              nelem=0
c              nup=ndo
c!     
c!     depth underneath sluice gate exceeds critical depth
c!     
c!     calculate the normal depth
c!     
c              if(xks.gt.0.d0) then
c                reynolds=xflow/(b*dvi)
c                form_fact=1.d0
c                hd=4.d0*hk
c                call friction_coefficient(dl,hd,xks,reynolds,form_fact,
c     &               friction)
c              endif
c              call hnorm(xflow,rho,b,theta,dg,s0,friction,xks,he)
c!     
c              if(he.lt.hk) then
c!     
c!     B2
c!     
c                area=(b+hk*tth)*hk
c                v(2,ndo)=(hk-epsilon)*sqrts0
c                mode='B'
c                nstack=nstack+1
c                istack(1,nstack)=nelem
c                istack(2,nstack)=ndo
c              else
c!     
c!     no frontwater solution
c!     
c                v(2,ndo)=-1.d0
c                nelup=nelem
c                nelem=0
c                nup=ndo
c              endif
c!
c!             calculate the depth in the upstream node;
c!             for a gate element or wear in between other elements
c!             the upstream velocity is not assumed to be zero,              
c!             i.e. no big reservoir since this makes no sense for
c!             steady state calculations
c!
c              if(mode.eq.'B') then
c                if(v(2,ndo).ge.ha*sqrts0) then
c                  call hns(xflow,rho,b,theta,dg,sqrts0,ha,h2)
c                  v(2,nup)=h2/sqrts0
cc                else
cc                  v(2,nup)=v(2,ndo)+hw
c                endif
c                neldo=nelem
c                ndo=nup
c                nelem=0
c              endif
c            endif
          else
!
!           no disturbance of the flow
!
            v(2,ndo)=v(2,nup)
            v(1,nmid)=inv*xflow
            nelup=nelem
            nelem=0
            nup=ndo
          endif
        endif
      else
!
!       mode = 'B': backwater curve
!
        hdo=v(2,ndo)/sqrts0
!
        idof=8*(nup-1)+2
        call nident(ikboun,idof,nboun,id)
        if(id.gt.0) then
          if(ikboun(id).eq.idof) then
!
!         hup is a boundary condition
!
            if(hdo.le.ha) then
              area=(b+hdo*tth)*hdo
            else
              area=(b+ha*tth)*ha
            endif
            xflowcor=rho*area*dsqrt(2.d0*dg*(v(2,nup)-hdo*sqrts0))
            if(dabs(xflow-xflowcor).le.1.d-3*xflow) then
!
!     corrected flow is sufficiently close to assumed flow:
!     solution found
!
              if(nstack.gt.0) then
                v(1,nmid)=xflow
                ndo=nup
                neldo=nelem
                nelem=0
                return
              endif
            else
!
!             new guess: mean of xflow and xflowcorr
!
              xflow=(xflowcor+xflow)/2.d0
              nelem=istack(1,nstack)
              ndo=istack(2,nstack)
              return
            endif
          endif
        endif
!
!       xflow is given
!
        if(hdo.le.ha) then
          area=(b+hdo*tth)*hdo
        else
          area=(b+ha*tth)*ha
        endif
        v(2,nup)=hdo*sqrts0+(xflow/(rho*area))**2/(2.d0*dg)
        ndo=nup
        neldo=nelem
        nelem=0
      endif
!      
      return
      end
      
