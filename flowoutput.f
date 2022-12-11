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
!     author: Yannick Muller
!     
      subroutine flowoutput(itg,ieg,ntg,nteq,bc,lakon,ntmat_,
     &     v,shcon,nshcon,ipkon,kon,co,nflow, dtime,ttime,time,
     &     ielmat,prop,ielprop,nactdog,nacteq,iin,physcon,
     &     camt,camf,camp,rhcon,nrhcon,vold,jobnamef,set,istartset,
     &     iendset,ialset,nset,mi,iaxial,istep,iit,ipobody,ibody,
     &     xbodyact,nbody,ndata,sfr,sba,jumpup,jumpdo)
!     
      implicit none
!     
      logical identity
!
      character*8 lakon(*)
      character*81 set(*)
      character*132 jobnamef(*),fnnet
!     
      integer mi(*),itg(*),ieg(*),ntg,nflow,ielmat(mi(3),*),i,j,
     &     nrhcon(*),node,iaxial,ider,idirf(8),ieq,imat,kflag,
     &     ntmat_,nteq,nshcon(*),nelem,index,ipkon(*),kon(*),iin,
     &     nactdog(0:3,*),nacteq(0:3,*),ielprop(*),node1,nodem,node2,
     &     istartset(*),iendset(*),ialset(*),nset,nodef(8),numf,
     &     istep,iit,iplausi,nup,ndo,inv,ipobody(2,*),ibody(3,*),
     &     nbody,nel,ndata,jumpup(*),jumpdo(*)
!     
      real*8 physcon(*),v(0:mi(2),*),shcon(0:3,ntmat_,*),co(3,*),
     &     prop(*),dtime,ttime,time,xflow,camp(*),camt(*),camf(*),
     &     rhcon(0:1,ntmat_,*),vold(0:mi(2),*),eta,xks,zup,zdo,
     &     bc(*),cp,dvi,df(8),gastemp,f,g(3),r,rho,ts1,ts2,thup,thdo,
     &     hup,hdo,thetaup,thetado,seup,sedo,theta1,theta2,temp,
     &     theta,sqrts0,s0,reynoldsup,reynoldsdo,hkup,hkdo,heup,hedo,
     &     xbodyact(7,*),hw,hd,ha,frictionup,frictiondo,form_fact,d,
     &     dg,dl,cd,bup,bdo,b,b1,b2,areaup,areado,uup,udo,sfr(*),sba(*)
!
!     output element per element in a dedicated file (jobname.net)      
!     
      do i=1,132
         if(jobnamef(1)(i:i).eq.' ') exit
      enddo
      i=i-1
      fnnet=jobnamef(1)(1:i)//'.net'
      open(1,file=fnnet,status='unknown')
!     
      do i=1,nflow
        nelem=ieg(i)
!     
!     output for channel networks
!     
        if((lakon(nelem)(2:5).ne.'LICH').or.
     &     (lakon(nelem)(6:7).eq.'IO')) cycle
!     
        index=ipkon(nelem)
        node1=kon(index+1)
        nodem=kon(index+2)
        node2=kon(index+3)
!
!       mass flow
!
        xflow=v(1,nodem)
!
!       check the flow direction; the temperature for the
!       material properties is the temperature of the upstream
!       node
!
        if(xflow.ge.0.d0) then
          inv=1
          nup=node1
          ndo=node2
        else
          inv=-1
          nup=node2
          ndo=node1
        endif
!
!       material of the element
!
        temp=v(0,nup)
        imat=ielmat(1,nelem)
!     
        call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
     &         dvi,rhcon,nrhcon,rho)
!     
!       determine the gravity vector
!     
        do j=1,3
          g(j)=0.d0
        enddo
        if(nbody.gt.0) then
          index=nelem
          do
            j=ipobody(1,index)
            if(j.eq.0) exit
            if(ibody(1,j).eq.2) then
              g(1)=g(1)+xbodyact(1,j)*xbodyact(2,j)
              g(2)=g(2)+xbodyact(1,j)*xbodyact(3,j)
              g(3)=g(3)+xbodyact(1,j)*xbodyact(4,j)
            endif
            index=ipobody(2,index)
            if(index.eq.0) exit
          enddo
        endif
        dg=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
!     
!       elevation of upstream and downstream floor
!     
        zup=(-g(1)*co(1,nup)-g(2)*co(2,nup)-g(3)*co(3,nup))/dg
        zdo=(-g(1)*co(1,ndo)-g(2)*co(2,ndo)-g(3)*co(3,ndo))/dg
!
!       determining the properties
!
        index=ielprop(nelem)
!
        if((lakon(nelem)(6:7).eq.'SG').or.
     &       (lakon(nelem)(6:7).eq.'WE').or.
     &       (lakon(nelem)(6:7).eq.'  ').or.
     &       (lakon(nelem)(6:7).eq.'RE')) then
          b=prop(index+1)
          theta=prop(index+2)
          dl=prop(index+3)
          s0=prop(index+4)
          xks=prop(index+5)
          if(lakon(nelem)(6:7).eq.'SG') then
            ha=prop(index+6)
          elseif(lakon(nelem)(6:7).eq.'WE') then
            hw=prop(index+6)
            cd=prop(index+7)*(1.5d0)**(1.5d0)/dsqrt(dg)
          elseif(lakon(nelem)(6:7).eq.'  ') then
            if(dl.le.0.d0) then
              dl=dsqrt((co(1,nup)-co(1,ndo))**2+
     &             (co(2,nup)-co(2,ndo))**2+
     &             (co(3,nup)-co(3,ndo))**2)
            endif
            if(s0.lt.-1.d0) then
              s0=(zup-zdo)/dl
            endif
          endif
          bup=b
          bdo=b
          thetaup=theta
          thetado=theta
        elseif((lakon(nelem)(6:7).eq.'CO').or.
     &         (lakon(nelem)(6:7).eq.'EL').or.
     &         (lakon(nelem)(6:7).eq.'ST')) then
c     &         (lakon(nelem)(6:7).eq.'DR')) then
          b1=prop(index+1)
          theta1=prop(index+2)
          b2=prop(index+3)
          theta2=prop(index+4)
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
          endif
          d=prop(index+5)
          dl=prop(index+6)
          if(dl.le.0.d0) then
            dl=dsqrt((co(1,nup)-co(1,ndo))**2+
     &           (co(2,nup)-co(2,ndo))**2+
     &           (co(3,nup)-co(3,ndo))**2)
          endif
          s0=prop(index+7)
        endif
!        
        sqrts0=1.d0-s0*s0
        if(sqrts0.lt.0.d0) then
          sqrts0=0.d0
        else
          sqrts0=dsqrt(sqrts0)
        endif
!
!       upstream and downstream depth
!
        if(lakon(nelem)(6:7).ne.'SG') then
          hup=v(2,nup)/sqrts0
        else
          hup=v(2,nup)
        endif
        if(lakon(nelem)(6:7).ne.'WE') then
          hdo=v(2,ndo)/sqrts0
        else
          hdo=v(2,ndo)
        endif
!
!       calculate the critical depth
!
        if((lakon(nelem)(6:7).eq.'SG').or.
     &     (lakon(nelem)(6:7).eq.'WE')) then
          hkup=0.d0
        else
          call hcrit(xflow,rho,bup,thetaup,dg,sqrts0,hkup)
        endif
        call hcrit(xflow,rho,bdo,thetado,dg,sqrts0,hkdo)
!     
!       calculate the normal depth
!     
        form_fact=1.d0
        reynoldsup=xflow/(bup*dvi)
        if(xks.gt.0.d0) then
          hd=4.d0*hup
          call friction_coefficient(dl,hd,xks,reynoldsup,form_fact,
     &         frictionup)
        endif
        call hnorm(xflow,rho,bup,thetaup,dg,s0,frictionup,xks,heup)
        reynoldsdo=xflow/(bdo*dvi)
!     
        if(lakon(nelem)(6:7).eq.'RE') then
          hedo=0.d0
          if(xks.gt.0.d0) frictiondo=0.d0
        else
          if(xks.gt.0.d0) then
            hd=4.d0*hdo
            call friction_coefficient(dl,hd,xks,reynoldsdo,form_fact,
     &           frictiondo)
          endif
          call hnorm(xflow,rho,bdo,thetado,dg,s0,frictiondo,xks,hedo)
        endif
!     
!       calculate the velocity
!     
        if((lakon(nelem)(6:7).eq.'SG').or.
     &     (lakon(nelem)(6:7).eq.'WE')) then
          uup=0.d0
        else
          areaup=(bup+hup*dtan(thetaup))*hup
          uup=xflow/(rho*areaup)
        endif
        if(lakon(nelem)(6:7).eq.'RE') then
          udo=0.d0
        else
          areado=(bdo+hdo*dtan(thetado))*hdo
          udo=xflow/(rho*areado)
        endif
!     
!       calculate the specific energy
!     
        if((lakon(nelem)(6:7).eq.'SG').or.
     &     (lakon(nelem)(6:7).eq.'WE')) then
          seup=hup
        else
          seup=hup*sqrts0+uup**2/(2.d0*dg)
        endif
        if(lakon(nelem)(6:7).eq.'RE') then
          sedo=hdo
        else
          sedo=hdo*sqrts0+udo**2/(2.d0*dg)
        endif
!     
!       calculate the total head
!     
        thup=seup+zup
        thdo=sedo+zdo
!     
!       output     
!     
        write(1,*) ''
        write(1,55) ' from node ',node1,
     &       ' to node ', node2,': massflow rate = ',xflow
!     
        write(1,*)'             Element ',nelem,lakon(nelem)
        if(xks.lt.0.d0) then
          write(1,57)'             dvi = ',dvi,', Re = '
     &         ,reynoldsup,', Manning = ',-xks
        else
          write(1,57)'             dvi = ',dvi,', Re = '
     &         ,reynoldsup
        endif
        write(1,57)'             S0 = ',s0,', Length = '
     &       ,dl
        call nident(ieg,nelem,nflow,nel)
        if((jumpup(nel).gt.0).and.(jumpdo(nel).lt.ndata+1)) then
          write(1,57)'             jump between ',
     &         sfr((nel-1)*ndata+jumpup(nel)),
     &         ' and ',sba((nel-1)*ndata+jumpdo(nel)),
     &         ' from upstream node '
        endif
!     
        write(1,*) ''
        write(1,53)' Inlet node ',nup,':  T = ',v(0,nup),
     &       ', fluid depth = ',hup*sqrts0
        write(1,54)'                     b = ',bup,
     &       ', theta = ',thetaup
        write(1,54)'                     velocity = ',uup,
     &       ', Froude = ',uup/dsqrt(dg*hup*sqrts0)
        write(1,54)'                     critical depth = ',hkup,
     &       ', normal depth = ',heup
        if(xks.gt.0) then
          write(1,54)'                     total head = ',thup,
     &         ', specific energy = ',seup,', f = ',frictionup
        else
          write(1,54)'                     total head = ',thup,
     &         ', specific energy = ',seup
        endif
!     
        write(1,*) ''
        write(1,53)' Outlet node ',ndo,': T = ',v(0,ndo),
     &       ', fluid depth  = ',hdo*sqrts0
        write(1,54)'                     b = ',bdo,
     &       ', theta = ',thetado
        write(1,54)'                     velocity = ',udo,
     &       ', Froude = ',udo/dsqrt(dg*hdo*sqrts0)
        write(1,54)'                     critical depth = ',hkdo,
     &       ', normal depth = ',hedo
        if(xks.gt.0) then
          write(1,54)'                     total head = ',thdo,
     &         ', specific energy = ',sedo,', f = ',frictiondo
        else
          write(1,54)'                     total head = ',thdo,
     &         ', specific energy = ',sedo
        endif
        write(1,54) '***************************************************
     &***************************'
      enddo
 54   format(1X,a,e11.4,a,e11.4)
 53   format(1X,a,i6,a,e11.4,a,e11.4,a,e11.4,a,e11.4)
 55   format(1X,a,i6,a,i6,a,e11.4,a,e11.4)
 57   format(1X,a,e11.4,a,e11.4,a,e11.4)
!     
       kflag=3
!     
       do i=1,nflow
         nelem=ieg(i)
!     
!     output for gas networks
!     
         if(lakon(nelem)(2:5).ne.'LICH') then
!     
           index=ipkon(nelem)
           node1=kon(index+1)
           nodem=kon(index+2)
           node2=kon(index+3)
!     
           xflow=v(1,nodem)
!     
           if((lakon(nelem)(2:3).ne.'LP').and.
     &          (lakon(nelem)(2:3).ne.'LI')) then
!     
!     compressible
!     
             if(node1.eq.0) then
               ts1=v(3,node2)
               ts2=ts1
             elseif(node2.eq.0) then
               ts1=v(3,node1)
               ts2=ts1
             else
               ts1=v(3,node1)
               ts2=v(3,node2)
             endif
             gastemp=(ts1+ts2)/2.d0
           else
!     
!     incompressible
!     
             if((xflow.gt.0).and.(node1.ne.0)) then
               gastemp=v(3,node1)
             else
               gastemp=v(3,node2)
             endif
           endif
!     
           imat=ielmat(1,nelem)
!     
           call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,
     &          dvi,rhcon,nrhcon,rho)
!     
           if(nacteq(2,nodem).ne.0) then
             ieq=nacteq(2,nodem)
!     
!     dummy set number
!     
             numf=1
!     
             call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &            nactdog,identity,
     &            ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &            cp,r,rho,physcon,g,co,dvi,numf,vold,set,shcon,
     &            nshcon,rhcon,nrhcon,ntmat_,mi,ider,ttime,time,
     &            iaxial,iplausi)
           endif
         endif
       enddo
!     
       close(1)
!     
       return
       end
