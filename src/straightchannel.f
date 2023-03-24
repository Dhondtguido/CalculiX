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
      subroutine straightchannel(nelem,ielprop,prop,nup,nmid,ndo,co,g,
     &     dg,mode,xflow,rho,dvi,nelup,neldo,istack,nstack,ikboun,nboun,
     &     mi,v,ipkon,kon,ndata,nel,sfr,hfr,sba,hba,jumpup,jumpdo,inv,
     &     epsilon)
!
      implicit none
!
      logical jump
!
      character*1 mode
!
      integer nelem,ielprop(*),index,nup,ndo,nmid,nelup,neldo,nstack,
     &     istack(2,*),ikboun(*),nboun,mi(*),id,idof,ipkon(*),kon(*),
     &     ndata,nel,jumpdo(*),jumpup(*),i,j,k,inv
!
      real*8 prop(*),co(3,*),g(3),b,theta,dl,s0,sqrts0,xks,hdo,
     &     ha,hup,xflow,dg,rho,hk,tth,area,reynolds,dvi,friction,hkmax,
     &     aa,bb,cc,form_fact,he,v(0:mi(2),*),zdo,zup,epsilon,
     &     sfr(*),hfr(*),sba(*),hba(*),sg,s,hback,h1,h2,dhds1,dhds2,
     &     dhdsm,dh,cthi,dhds,delta,hd
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
      cthi=1.d0/dcos(theta)
!
!     if the length of the element is negative, it is determined from
!     the coordinates
!
      dl=prop(index+3)
      if(dl.le.0.d0) then
        dl=dsqrt((co(1,nup)-co(1,ndo))**2+
     &           (co(2,nup)-co(2,ndo))**2+
     &       (co(3,nup)-co(3,ndo))**2)
      endif
!
!     s0: sine of slope (the slope is the angle phi between the channel
!         bottom and a plane orthogonal to the gravity vector
!     sqrts0: cosine of slope
!      
!     determining the sine of the slope; if the sine is less than
!     -1.d0 it is calculated from the coordinates
!
      s0=prop(index+4)
      if(s0.lt.-1.d0) then
        zup=(-g(1)*co(1,nup)-g(2)*co(2,nup)-g(3)*co(3,nup))/dg
        zdo=(-g(1)*co(1,ndo)-g(2)*co(2,ndo)-g(3)*co(3,ndo))/dg
        s0=(zup-zdo)/dl
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
      v(1,nmid)=inv*xflow
!
!     calculate the critical depth
!
      call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
!
!     forward mode (frontwater curve)
!
      if(mode.eq.'F') then
!
        v(3,nup)=hk
        hup=v(2,nup)/sqrts0
!     
!     calculate the normal depth
!     
        if(xks.gt.0.d0) then
          reynolds=xflow/(b*dvi)
          form_fact=1.d0
          if(hup.lt.0.d0) then
            hd=4.d0*hk
            call friction_coefficient(dl,hd,xks,reynolds,form_fact,
     &           friction)
          else
            hd=4.d0*hup
            call friction_coefficient(dl,hd,xks,reynolds,form_fact,
     &           friction)
          endif
        endif
        call hnorm(xflow,rho,b,theta,dg,s0,friction,xks,he)
!     
        if((hup.lt.0.d0).or.(hup.gt.hk-epsilon)) then
!
          if(hk.lt.he) then
!
!           weak slope: no solution
!
            v(2,ndo)=-1.d0
            nelup=nelem
            nelem=0
            nup=ndo
            return
          else
!
!           strong slope: B2
!
            if(hup.lt.0.d0) then
!
!             first calculate the backwater curve starting in nup
!
c              v(2,nup)=(hk+epsilon)*sqrts0
              v(2,nup)=(hk)*sqrts0
              ndo=nup
              nelem=nelup
              mode='B'
              nstack=nstack+1
              istack(1,nstack)=nelup
              istack(2,nstack)=nup
              return
            endif
          endif
        endif
!     
        if(hup.gt.hk-epsilon) hup=hk-epsilon
!     
!     determine the h-increment
!     
        if(hk.lt.he) then
!     
!         A3
!     
          dh=(hk-epsilon-hup)/(ndata-1)
        else
!     
!         B2 or B3
!     
          dh=(he-hup)/(ndata-1)
        endif
!     
        sfr(1)=0.d0
        hfr(1)=hup
!     
!     calculate dhds for hfr(1)
!     
        call calcdhds(xflow,b,tth,cthi,s0,sqrts0,friction,xks,hup,
     &       dg,rho,dhds1)
!        
        do k=2,ndata
          hfr(k)=hfr(k-1)+dh
!     
!     calculate dhds for hfr(k)
!     
          call calcdhds(xflow,b,tth,cthi,s0,sqrts0,friction,xks,hfr(k),
     &         dg,rho,dhds2)
!     
!     taking the mean
!     
          dhdsm=(dhds1+dhds2)/2.d0
          dhds1=dhds2
!     
          sfr(k)=sfr(k-1)+dh/dhdsm
!     
!     check whether length coordinate extends beyond end of element
!     
          if(sfr(k).gt.dl) then
            hfr(k)=hfr(k-1)+(dl-sfr(k-1))/(sfr(k)-sfr(k-1))*
     &           (hfr(k)-hfr(k-1))
            sfr(k)=dl
            jumpup(nel)=k
!     
!           initialization jumpdo(nel)
!     
            jumpdo(nel)=ndata+1
            v(2,ndo)=hfr(k)*sqrts0
            exit
          endif
        enddo
!     
        if(k.le.ndata) then
!     
!     complete frontwater curve
!
          nelup=nelem
          nelem=0
          nup=ndo
          return
        else
!     
!     partially backwater curve
!     
          v(2,ndo)=-1.d0
          jumpup(nel)=ndata
!     
!         initialization jumpdo(nel)
!     
          jumpdo(nel)=ndata+1
          nelup=nelem
          nelem=0
          nup=ndo
          return
        endif
      else
!     
!     mode='B'
!     
        v(3,ndo)=hk
        hdo=v(2,ndo)/sqrts0
        if(hdo.lt.hk+epsilon) hdo=hk+epsilon
!     
!     calculate the normal depth
!     
        if(xks.gt.0.d0) then
          reynolds=xflow/(b*dvi)
          form_fact=1.d0
          hd=4.d0*hdo
          call friction_coefficient(dl,hd,xks,reynolds,form_fact,
     &         friction)
        endif
        call hnorm(xflow,rho,b,theta,dg,s0,friction,xks,he)
!     
!     determine the h-increment
!     
        if(hk.gt.he) then
! 
!         B1
!
          dh=(hdo-hk-epsilon)/(ndata-1)
        else
! 
!         A1 or A2
!
          dh=(hdo-he)/(ndata-1)
        endif
!
        hba(ndata)=hdo
        sba(ndata)=dl
!     
!       calculate dhds for hfr(ndata)
!     
        call calcdhds(xflow,b,tth,cthi,s0,sqrts0,friction,xks,hdo,
     &       dg,rho,dhds2)
!
        jumpdo(nel)=1
        do k=ndata-1,1,-1
          hba(k)=hba(k+1)-dh
          call calcdhds(xflow,b,tth,cthi,s0,sqrts0,friction,xks,hba(k),
     &         dg,rho,dhds1)
          dhdsm=(dhds1+dhds2)/2.d0
c          dhdsm=dhds2
          dhds2=dhds1
!
          sba(k)=sba(k+1)-dh/(dhdsm)
          if(sba(k).lt.0.d0) then
            hba(k)=hba(k+1)-sba(k+1)/(sba(k)-sba(k+1))*
     &           (hba(k)-hba(k+1))
            sba(k)=0.d0
            jumpdo(nel)=k
            exit
          endif
        enddo
!
!       check for the location of a jump
!
        sg=0.d0
        jump=.false.
!  
!       search for a jump only if there is a frontwater curve
!
        if(jumpup(nel).gt.0) then
          do k=ndata,jumpdo(nel),-1
            if(k.eq.ndata) then
              s=sba(ndata)-epsilon
            else
              s=sba(k)
            endif
            hback=hba(k)
            call ident(sfr,s,jumpup(nel),id)
            if(id.lt.jumpup(nel)) then
              h1=hfr(id)+(s-sfr(id))/(sfr(id+1)-sfr(id))*
     &             (hfr(id+1)-hfr(id))
              call hns(xflow,rho,b,theta,dg,sqrts0,h1,h2)
c              write(*,*) 'B jump ',s,h1,h2,hback
              delta=h2-hback
              if(delta*sg.eq.0.d0) then
                sg=delta
              elseif(delta*sg.lt.0.d0) then
!     
!     jump between point id in the frontwater curve and
!     point k+1 in the backwater curve
!
                write(*,*) '*INFO in straightchannel'
                write(*,*) '      jump detected between ',sfr(id)
                write(*,*) '      and ',sba(k+1),' length units'
                write(*,*) '      from upstream node ',nup
                write(*,*) '      in element ',nelem
                write(*,*)
                jumpup(nel)=id
                jumpdo(nel)=k+1
                jump=.true.
                exit
              endif
            endif
          enddo
        endif
!  
!       complete backwater curve in this element
!
        if(.not.jump) then
          jumpup(nel)=0
          v(2,nup)=hba(jumpdo(nel))*sqrts0
          ndo=nup
          neldo=nelem
          nelem=0
        else
          mode='F'
          nelup=istack(1,nstack)
          nelem=0
          nup=istack(2,nstack)
          nstack=nstack-1
        endif
      endif
!
      return
      end
      
