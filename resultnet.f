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
!     construction of the b matrix      
!     
      subroutine resultnet(itg,ieg,ntg,
     &     bc,nload,sideload,nelemload,xloadact,lakon,ntmat_,
     &     v,shcon,nshcon,ipkon,kon,co,nflow,
     &     iinc,istep,dtime,ttime,time,
     &     ikforc,ilforc,xforcact,nforc,ielmat,nteq,prop,ielprop,
     &     nactdog,nacteq,iin,physcon,camt,camf,camp,rhcon,nrhcon,
     &     ipobody,ibody,xbodyact,nbody,dtheta,vold,xloadold,
     &     reltime,nmethod,set,mi,ineighe,cama,vamt,vamf,vamp,vama,
     &     nmpc,nodempc,ipompc,coefmpc,labmpc,iaxial,qat,qaf,ramt,
     &     ramf,ramp,cocon,ncocon,iponoel,inoel,iplausi)
!     
      implicit none
!     
      logical identity
      character*8 lakonl,lakon(*)
      character*20 sideload(*),labmpc(*)
      character*81 set(*)
!     
      integer mi(*),itg(*),ieg(*),ntg,nteq,nflow,nload,
     &     ielmat(mi(3),*),iflag,ider,iaxial,nalt,nalf,
     &     nelemload(2,*),nope,nopes,mint2d,i,j,k,m,nrhcon(*),
     &     node,imat,ntmat_,id,ifaceq(8,6),ifacet(6,4),numf,
     &     ifacew(8,5),node1,node2,nshcon(*),nelem,ig,index,konl(20),
     &     ipkon(*),kon(*),idof,ineighe(*),idir,ncocon(2,*),
     &     iinc,istep,jltyp,nfield,ikforc(*),ipobody(2,*),
     &     ilforc(*),nforc,nodem,idirf(8),ieq,nactdog(0:3,*),nbody,
     &     nacteq(0:3,*),ielprop(*),nodef(8),iin,kflag,ibody(3,*),icase,
     &     inv, index2,nmethod,nelem0,nodem0,nelem1,nodem1,nelem2,
     &     nodem2,nelemswirl,nmpc,nodempc(3,*),ipompc(*),
     &     iponoel(*),inoel(2,*),iplausi,indexe
!     
      real*8 bc(nteq),xloadact(2,*),cp,h(2),physcon(*),r,dvi,rho,
     &     xl2(3,8),coords(3),dxsj2,temp,xi,et,weight,xsj2(3),
     &     gastemp,v(0:mi(2),*),shcon(0:3,ntmat_,*),co(3,*),shp2(7,8),
     &     field(1),prop(*),tg1,tg2,dtime,ttime,time,g(3),eta,
     &     xforcact(*),areaj,xflow,tvar(2),f,df(8),camt(*),camf(*),
     &     camp(*),tl2(8),cama(*),vamt,vamf,vamp,vama,term,
     &     rhcon(0:1,ntmat_,*),xbodyact(7,*),sinktemp,kappa,A,T,Tt,pt,
     &     dtheta,ts1,ts2,xs2(3,7),pt1,pt2,Qred_crit,xflow_crit,
     &     xflow0,xflow1,reltime,coefmpc(*),qat,qaf,ramt,ramf,ramp,
     &     xflow2,heat,cocon(0:6,ntmat_,*),
     &     Cp_cor,U,Ct,vold(0:mi(2),*),xloadold(2,*),bnac,
     &     xflow360,Tt1,Tt2,T1,T2,heatnod,heatfac,A2,d,l,s
!     
      include "gauss.f"
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
      data iflag /2/
!     
      kflag=2
      ider=0
!     
      tvar(1)=time
      tvar(2)=ttime+time
!     
      heatnod=0.d0
!     
!     calculating the maximum correction to the solution in the
!     present network iteration
!     
      camt(1)=0.d0
      camf(1)=0.d0
      camp(1)=0.d0
      cama(1)=0.d0
!     
      camt(2)=0.5d0
      camf(2)=0.5d0
      camp(2)=0.5d0
      cama(2)=0.5d0
!     
!     maximum change in the solution since the start of the
!     network iterations (= entering radflowload)
!     
      vamt=0.d0
      vamf=0.d0
      vamp=0.d0
      vama=0.d0
!     
      do i=1,ntg
        node=itg(i)
        do j=0,3
          if(nactdog(j,node).eq.0) cycle
          idof=nactdog(j,node)
!     
          if(dabs(v(j,node)).gt.1.d-30) then
            if(dabs(bc(idof)/v(j,node)).lt.1.d-13) bc(idof)=0.d0
          endif
          bnac=bc(idof)
!     
          if(j.eq.0) then
            if(dabs(bnac).gt.dabs(camt(1))) then
              camt(1)=bnac
              camt(2)=node+0.5d0
            endif
          elseif(j.eq.1) then
            if(dabs(bnac).gt.dabs(camf(1))) then
              camf(1)=bnac
              camf(2)=node+0.5d0
            endif
          elseif(j.eq.2) then
            if(dabs(bnac).gt.dabs(camp(1))) then
              camp(1)=bnac
              camp(2)=node+0.5d0
            endif
          else
            if(dabs(bnac).gt.dabs(cama(1))) then
              cama(1)=bnac
              cama(2)=node+0.5d0
            endif
          endif
        enddo
      enddo
!     
!     updating v
!     vold is the solution at the entry of radflowload (= at
!     the start of the network iterations)
!     
      do i=1,ntg
        node=itg(i)
        do j=0,2
          if(nactdog(j,node).eq.0) cycle
          v(j,node)=v(j,node)+bc(nactdog(j,node))*dtheta
!     
          if((j.eq.0).and.(dabs(v(j,node)-vold(j,node)).gt.vamt)) then
            vamt=dabs(v(j,node)-vold(j,node))
          elseif((j.eq.1).and.
     &           (dabs(v(j,node)-vold(j,node)).gt.vamf)) then
            vamf=dabs(v(j,node)-vold(j,node))
          elseif((j.eq.2).and.
     &           (dabs(v(j,node)-vold(j,node)).gt.vamp)) then
            vamp=dabs(v(j,node)-vold(j,node))
          endif
!     
        enddo
      enddo
c     write(*,*) 'resultnet'
c     write(*,*) 'mass flow in node 3: ',v(1,3)
c     write(*,*) 'pressure in node 4: ',v(2,4)
c     write(*,*) 'mass flow in node 5: ',v(1,5)
!     
!     update geometry changes
!     
      do i=1,nflow
        if(lakon(ieg(i))(6:7).eq.'GV') then
          index=ipkon(ieg(i))
          node=kon(index+2)
          if(nactdog(3,node).eq.0) cycle
          index=ielprop(ieg(i))
          v(3,node)=v(3,node)+bc(nactdog(3,node))*dtheta
          if(v(3,node).gt.0.99999d0) then
            v(3,node)=0.99999d0
          elseif(v(3,node).lt.0.12501d0) then
            v(3,node)=0.12501d0
          endif
          if(dabs(v(3,node)).gt.vama) vama=dabs(v(3,node))
!     
!     Geometry factor of an ACC tube
!     for all versions of the ACC tube
!     
        elseif(lakon(ieg(i))(2:7).eq.'ACCTUB') then
          index=ipkon(ieg(i))
          node=kon(index+2)
          if(nactdog(3,node).eq.0) cycle
          index=ielprop(ieg(i))
!     Interval factor unknown
          if(nint(prop(index+1)).eq.2) then
!     Using a smaller step gets better convergence(factor 0.5)
            v(3,node)=v(3,node)+bc(nactdog(3,node))*
     &           dtheta
            if(v(3,node).lt.0.1d0)then
              v(3,node) = 0.1d0
            endif
!     Hole diameter factor unknown               
          elseif(nint(prop(index+1)).eq.3) then
            v(3,node)=v(3,node)+bc(nactdog(3,node))*dtheta
            if(v(3,node).lt.0.1d0)then
              v(3,node) = 0.1d0
            endif
          endif
          if(dabs(v(3,node)).gt.vama) vama=dabs(v(3,node))
        endif
      enddo
!     
c     write(*,*) 'resultnet'
c     do i=1,ntg
c     node=itg(i)
c     write(*,'(i10,3(1x,e11.4))') node,(v(j,node),j=0,2)
c     enddo
!     
!     testing the validity of the pressures
!     
      do i=1,ntg
        node=itg(i)
        if(v(2,node).lt.0) then
          write(*,*) '*WARNING in resultnet:'
          write(*,*) '         negative pressure; node ',node,v(2,node)
          iin=0
          return
        endif
      enddo
!     
!     testing validity of temperatures
!     
      do i=1,ntg
        node=itg(i)
        if(v(0,node).lt.0) then
          write(*,*) '*WARNING in resultnet:'
          write(*,*)
     &         '         negative temperature; node ',node,v(0,node)
          iin=0
          return
        endif
      enddo
!     
!     testing the validity of the solution for branches elements
!     and restrictor. Since the element properties are dependent on
!     a predefined flow direction a change of this will lead to
!     wrong head losses
!     
      do i=1,nflow
        nelem=ieg(i)
        if ((lakon(nelem)(4:5).eq.'ATR').or. 
     &       (lakon(nelem)(4:5).eq.'RTA')) then
          xflow=v(1,kon(ipkon(nelem)+2))
          if(xflow.lt.0d0)then
            Write(*,*) '*WARNING in resultnet.f'
            write(*,*) '         Element',nelem,
     &           ' of TYPE ABSOLUTE TO RELATIVE'
            write(*,*) '         The flow direction is no more conform '
            write(*,*) '         to element definition'
            write(*,*) '         Check the pertinence of the results'
          endif
        elseif(lakon(nelem)(2:3).eq.'RE') then
!     
          if(lakon(nelem)(4:5).ne.'BR') then
            nodem=kon(ipkon(nelem)+2)
            xflow=v(1,nodem)
            if (xflow.lt.0) then
              Write(*,*) '*WARNING in resultnet.f'
              write(*,*) '         Element',nelem,'of TYPE RESTRICTOR'
              write(*,*)
     &             '         The flow direction is no more conform '
              write(*,*) '         to element definition'
              write(*,*) '         Check the pertinence of the results'
            endif
!     
          elseif(lakon(nelem)(4:5).eq.'BR') then
            index=ielprop(nelem)
!     
            nelem0=nint(prop(index+1))
            nodem0=kon(ipkon(nelem0)+2)
            xflow0=v(1,nodem0)
!     
            nelem1=nint(prop(index+2))
            nodem1=kon(ipkon(nelem1)+2)
            xflow1=v(1,nodem1)
!     
            nelem2=nint(prop(index+3))
            nodem2=kon(ipkon(nelem2)+2)
            xflow2=v(1,nodem2)
!     
            if((xflow0.lt.0).or.(xflow1.lt.0).or.(xflow2.lt.0)) then
              Write(*,*) '*WARNING in resultnet.f'
              write(*,*) '         Element',nelem,'of TYPE BRANCH'
              write(*,*)
     &             '         The flow direction is no more conform '
              write(*,*) '         to element definition'
              write(*,*) '         Check the pertinence of the results'
            endif
          endif
        endif
      enddo
!     
!     determining the static temperature and checking for
!     critical conditions (only for non-chambers, i.e. for
!     nodes purely belonging to gas pipes or restrictors except
!     restrictor wall orifice)
!     
      iplausi=1
      do i=1,ntg
        node=itg(i)
        nelem=ineighe(i)
        if(nelem.eq.-1) then
          v(3,node)=v(0,node)
        elseif(nelem.gt.0) then
!     
          nodem=kon(ipkon(nelem)+2)
          T=v(3,node)
          Tt=v(0,node)
          pt=v(2,node)
          xflow=v(1,nodem)
!     
          icase=0
          inv=1
          imat=ielmat(1,nelem)
          call materialdata_tg(imat,ntmat_,v(3,node),
     &         shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
!     
          index=ielprop(nelem)
          indexe=ipkon(nelem)
!     
          kappa=(cp/(cp-R))
!     
          if(lakon(nelem)(2:5).eq.'GAPF') then
!     
!     distinguish between standard pipes and flexible
!     pipes
!     
            if((lakon(nelem)(6:8).eq.'AFR').or.
     &           (lakon(nelem)(6:8).eq.'ARL').or.
     &           (lakon(nelem)(6:8).eq.'ARG').or.
     &           (lakon(nelem)(6:8).eq.'IFR').or.
     &           (lakon(nelem)(6:8).eq.'IRL')) then
              call calcgeomelemnet(vold,co,prop,lakon,nelem,
     &             ttime,time,ielprop,mi,A,A2,d,l,s)
            else
              A=prop(index+1)
            endif
!     
            if(lakon(nelem)(2:6).eq.'GAPFA') then
              icase=0
            elseif(lakon(nelem)(2:6).eq.'GAPFI') then
              icase=1
            endif  
!     
          elseif(lakon(nelem)(2:5).eq.'GAPR') then
!     
!     rotating adiabatic gas pipe with variable cross section:
!     check whether section 1 or 2                  
!     
            if(kon(ipkon(nelem)+1).eq.node) then
              A=prop(index+1)
            else
              A=prop(index+2)
            endif
            icase=0
!     
          elseif(lakon(nelem)(2:3).eq.'RE') then
            node1=kon(indexe+1)
            node2=kon(indexe+3)
!     
            if(lakon(nelem)(4:5).eq.'EX') then
              if(lakon(nint(prop(index+4)))(2:6).eq.'GAPFA') then
                icase=0
              elseif(lakon(nint(prop(index+4)))(2:6).eq.'GAPFI')then
                icase=1
              endif
            else
              icase=0
            endif
!     
!     defining the sections
!     
            if(lakon(nelem)(4:5).eq.'BE') then
              A=prop(index+1)
!     
            elseif(lakon(nelem)(4:5).eq.'BR') then
              if(lakon(nelem)(4:6).eq.'BRJ') then
                if(nelem.eq.nint(prop(index+2)))then
                  A=prop(index+5)
                elseif(nelem.eq.nint(prop(index+3))) then
                  A=prop(index+6)
                endif
              elseif(lakon(nelem)(4:6).eq.'BRS') then
                if(nelem.eq.nint(prop(index+2)))then
                  A=prop(index+5)
                elseif(nelem.eq.nint(prop(index+3))) then
                  A=prop(index+6)
                endif
              endif
!     
            else
!     
              if((lakon(nelem)(2:5).eq.'RBSF')) then
                call calcgeomelemnet(vold,co,prop,lakon,nelem,
     &               ttime,time,ielprop,mi,A,A2,d,l,s)
                if(node.eq.node2) then
                  A=A2
                endif
              else
                if(node.eq.node1) then
                  A=prop(index+1)
                elseif(node.eq.node2) then
                  A=prop(index+2)
                endif
              endif
            endif
          elseif(lakon(nelem)(2:3).eq.'UP') then
!     
!     user elements whose names start with UP are assumed
!     to be Pipe-like (static and total temperatures differ)
!     
            call calcgeomelemnet(vold,co,prop,lakon,nelem,
     &           ttime,time,ielprop,mi,A,A2,d,l,s)
c     A=prop(index+1)
            icase=0
          endif
!     
          xflow360=dabs(xflow*iaxial)
          call ts_calc(xflow360,Tt,pt,kappa,r,A,T,icase)
!     
          v(3,node)=T
!     
          if(lakon(nelem)(2:3).eq.'RE') then
!     
!     check whether the mass flow does not exceed
!     critical conditions (only for restrictor elements)
!     
            if(xflow.lt.0d0) then
              inv=-1
            else
              inv=1
            endif
!     
            if(icase.eq.0) then
              Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*(kappa-1.d0))
     &             **(-0.5d0*(kappa+1.d0)/(kappa-1.d0))
            else
              Qred_crit=dsqrt(1/R)*(1.d0+0.5d0*(kappa-1.d0)/kappa)
     &             **(-0.5d0*(kappa+1.d0)/(kappa-1.d0))
            endif
            xflow_crit=Qred_crit*pt*A/dsqrt(Tt)   
!     
            if(xflow360.ge.xflow_crit) then
!     
!     next 3 lines: change on 20200929 (Thomas Petzke)
!     
              if(dabs((xflow_crit-xflow360)/xflow_crit).gt.1.d-5)
     &             then
                iplausi=0
              endif
!     
              v(1,nodem)=inv*xflow_crit/iaxial
              if(icase.eq.1) then
!     
                if(nactdog(0,node2).ne.0) then
                  node1=kon(indexe+1)
                  node2=kon(indexe+3)
                  v(3,node2)=v(3,node1)
                  v(0,node2)=v(3,node2)
     &                 *(1+0.5d0*(kappa-1)/kappa)
                  
                endif
              endif
            endif
          endif
!     
        endif
!     
      enddo
!     
!     testing if the flow direction is conform to the pressure
!     gradient
!     
      do i=1,nflow
        nelem=ieg(i)
        indexe=ipkon(nelem)
        node1=kon(indexe+1)
        nodem=kon(indexe+2)
        node2=kon(indexe+3)
!     
        if((node1.ne.0).and.(node2.ne.0)) then
!     
!     the mass flow rates must always be oriented in the direction
!     of the decreasing pressure gradient except for ACCTUBE, 
!     CROSS SPLIT, MOEHRING, ROTATING CAVITY, RADIAL OUTFLOW and 
!     INLET, RIMSEAL, S-PUMP and VORTEX
!     
          if((lakon(nelem)(2:8).ne.'ACCTUBE').and.
     &         (lakon(nelem)(2:5).ne.'ATR').and.
     &         (lakon(nelem)(2:5).ne.'RTA').and.   
     &         (lakon(nelem)(2:5).ne.'CROS').and.
     &         (lakon(nelem)(2:3).ne.'LI').and.
     &         (lakon(nelem)(2:3).ne.'LP').and.     
     &         (lakon(nelem)(2:4).ne.'MRG').and.
     &         (lakon(nelem)(2:4).ne.'RCV').and.
     &         (lakon(nelem)(2:3).ne.'RO').and.
     &         (lakon(nelem)(2:5).ne.'RIMS').and.
     &         (lakon(nelem)(2:5).ne.'FDPF').and.
     &         (lakon(nelem)(2:5).ne.'FCVF').and.
     &         (lakon(nelem)(2:5).ne.'MFPC').and.
     &         (lakon(nelem)(2:6).ne.'SPUMP').and.
     &         (lakon(nelem)(2:5).ne.'GAPR').and.
     &         (lakon(nelem)(2:3).ne.'VO')) then
            if(nactdog(1,nodem).ne.0) then
              if(((v(1,nodem).gt.0.d0).and.
     &             (v(2,node1)/v(2,node2).lt.(1.d0-1.d-10))).or.
     &             ((v(1,nodem).lt.0.d0).and.
     &             (v(2,node2)/v(2,node1).lt.(1.d0-1.d-10)))) then
                iplausi=0
                write(*,*) '*WARNING in resultnet: the flow'
                write(*,*) '         direction is not conform'
                write(*,*) '         to the pressure gradient'
                write(*,*) '         element',nelem
              endif
            endif
          endif
        endif
      enddo
!     
!     reinitialisation of the Bc matrix
!     
      do i=1,nteq
        bc(i)=0.d0
      enddo
!     
      qat=0.d0
      qaf=0.d0
      nalt=0
      nalf=0
!     
!     determining the residual
!     
      do i=1,nflow
        nelem=ieg(i)
        index=ipkon(nelem)
        node1=kon(index+1)
        nodem=kon(index+2)
        node2=kon(index+3)
        xflow=v(1,nodem)
!     
!     gas: the property temperature is the static temperature
!     
        if((lakon(nelem)(2:3).ne.'LP').and.
     &       (lakon(nelem)(2:3).ne.'LI')) then
          if(node1.eq.0) then
            tg1=v(0,node2)
            tg2=tg1
            ts1=v(3,node2)
            ts2=ts1
          elseif(node2.eq.0) then
            tg1=v(0,node1)
            tg2=tg1
            ts1=v(3,node1)
            ts2=ts1
          else
            tg1=v(0,node1)
            tg2=v(0,node2)
            ts1=v(3,node1)
            ts2=v(3,node2)
          endif
          gastemp=(ts1+ts2)/2.d0
        else
!     
!     liquid: only one temperature
!     
          if(xflow.gt.0) then
            if(node1.eq.0) then
              gastemp=v(0,node2)
            else
              gastemp=v(0,node1)
            endif
          else
            if(node2.eq.0) then
              gastemp=v(0,node1)
            else
              gastemp=v(0,node2)
            endif
          endif
!     
          if(node1.eq.0) then
            tg2=v(0,node2)
            tg1=tg2
          elseif(node2.eq.0) then
            tg1=v(0,node1)
            tg2=tg1
          else
            tg1=v(0,node1)
            tg2=v(0,node2)
          endif
        endif
!     
        imat=ielmat(1,nelem)
!     
        call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,dvi,
     &       rhcon,nrhcon,rho)
        kappa=Cp/(Cp-R)
!     
!     Definitions of the constant for isothermal flow elements
!     
        if(lakon(nelem)(2:6).eq.'GAPFI') then
          if((node1.ne.0).and.(node2.ne.0)) then
!     
            if((lakon(nelem)(7:8).eq.'FR').or.
     &           (lakon(nelem)(7:8).eq.'RL')) then
              call calcgeomelemnet(vold,co,prop,lakon,nelem,
     &             ttime,time,ielprop,mi,A,A2,d,l,s)
            else
              A=prop(ielprop(nelem)+1)
            endif
!     
            pt1=v(2,node1)
            pt2=v(2,node2)
!     
            xflow360=xflow*iaxial
            icase=1
!     
            if(pt1.ge.pt2)then
              if(dabs(tg2/ts2-(1.d0+0.5d0*(kappa-1)/kappa)).lt.1d-5)
     &             then
                pt2=dabs(xflow360)*dsqrt(Tg2*R)/A
     &               *(1.d0+0.5d0*(kappa-1.d0)/kappa)
     &               **(0.5d0*(kappa+1.d0)/(kappa-1.d0))
!     
              endif
!
!              Guido 17.06.2022 Option 2 (cf. initialnet.f)
!
c              Tt1=v(0,node1)-physcon(1)
c              call nident(itg,node1,ntg,id1)
c              if(ineighe(id1).eq.-1) then
c                T1=Tt1
c              else
c                call ts_calc(xflow360,Tt1,pt1,kappa,r,A,T1,icase)
c              endif
c              v(3,node1)=T1
!     
c              Tt2=v(0,node2)-physcon(1)
c              call nident(itg,node2,ntg,id2)
c              if(ineighe(id2).eq.-1) then
c                T2=Tt2
c              else
c                call ts_calc(xflow360,Tt2,pt2,kappa,r,A,T2,icase)
c              endif
c              v(3,node2)=T2
!                
              Tt1=v(0,node1)-physcon(1)
              Tt2=v(0,node2)-physcon(1)
              T1=v(3,node1)
              call ts_calc(xflow360,Tt1,pt1,kappa,r,A,T1,icase)
              call ts_calc(xflow360,Tt2,pt2,kappa,r,A,T2,icase)
              v(3,node1)=T1
              v(3,node2)=T2
            else
              pt1=v(2,node2)
              pt2=v(2,node1)
              if(v(2,nodem).ge.(pt2/pt1))then
                pt2=v(2,nodem)*pt1
              endif
!     
!             Option 2 ????  T not stored?
!     
!
              Tt1=v(0,node2)-physcon(1)
              call ts_calc(xflow360,Tt1,pt1,kappa,r,A,T1,icase)
              Tt2=v(0,node1)-physcon(1)
              call ts_calc(xflow360,Tt2,pt2,kappa,r,A,T2,icase)
            endif
!     
          endif
        endif
!     
        if(node1.ne.0) then
!     
!     energy equation contribution node1
!     
          if (nacteq(0,node1).ne.0) then
            ieq=nacteq(0,node1)
!     
            if(nacteq(3,node1).eq.0) then
              if (xflow.lt.0d0)then
                term=cp*(tg1-tg2)*xflow
                bc(ieq)=bc(ieq)+term
                qat=qat+dabs(term)
                nalt=nalt+1
              endif
!     
            elseif(lakon(nelem)(2:6).eq.'GAPFI') then
              if((nacteq(3,node1).eq.node2)) then
!     
                bc(ieq)=(T2-T1)
!     
              endif
            endif
          endif
!     
!     mass equation contribution node1
!     
          if (nacteq(1,node1).ne.0) then
            ieq=nacteq(1,node1)
            bc(ieq)=bc(ieq)-xflow
            qaf=qaf+dabs(xflow)
            nalf=nalf+1
          endif
        endif
!     
        if(node2.ne.0) then
!     
!     energy equation contribution node2
!     
          if (nacteq(0,node2).ne.0) then
            ieq=nacteq(0,node2)
!     
            if(nacteq(3,node2).eq.0) then
              if (xflow.gt.0d0)then
                term=cp*(tg2-tg1)*xflow
                bc(ieq)=bc(ieq)-term
                qat=qat+dabs(term)
                nalt=nalt+1
              endif
!     
            elseif(lakon(nelem)(2:6).eq.'GAPFI') then
              if(nacteq(3,node2).eq.node1) then
!     
                bc(ieq)=(T2-T1)
!     
              endif
            endif
          endif
!     
!     mass equation contribution node2
!     only in the case of an ACCTUBE but not in the case of an ACCTUBO
!     
          if(lakon(nelem)(2:8).ne.'ACCTUBE') then
            if (nacteq(1,node2).ne.0) then
              ieq=nacteq(1,node2)
              bc(ieq)=bc(ieq)+xflow
              qaf=qaf+dabs(xflow)
              nalf=nalf+1
            endif
          else
            if (nacteq(1,node2).ne.0) then
              if(nelem.ne.nint(prop(ielprop(nelem)+14))) then
                ieq=nacteq(1,node2)
                bc(ieq)=bc(ieq)+xflow-v(0,nodem)
                qaf=qaf+dabs(xflow)
                nalf=nalf+1
              else
                ieq=nacteq(1,node2)
                bc(ieq)=bc(ieq)+xflow
                qaf=qaf+dabs(xflow)
                nalf=nalf+1
              endif
            endif
          endif
        endif
!     
!     element equation
!     
        if (nacteq(2,nodem).ne.0) then
          ieq=nacteq (2,nodem)
!     
!     for liquids: determine the gravity vector
!     
          if((lakon(nelem)(2:3).eq.'LI').or.
     &         (lakon(nelem)(2:3).eq.'LP')) then
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
          endif
!     
          call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &         nactdog,identity,
     &         ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &         cp,r,rho,physcon,g,co,dvi,numf,vold,set,shcon,
     &         nshcon,rhcon,nrhcon,ntmat_,mi,ider,ttime,time,
     &         iaxial,iplausi)
          bc(ieq)=-f
        endif
      enddo
!     
!     convection with the walls: contribution to the energy equations
!     
      do i=1,nload
        if(sideload(i)(3:4).eq.'FC') then
          nelem=nelemload(1,i)
          index=ipkon(nelem)
          if(index.lt.0) cycle
          lakonl=lakon(nelem)
          node=nelemload(2,i)
          ieq=nacteq(0,node)
          if(ieq.eq.0) then 
            cycle
          endif
!     
c          call nident(itg,node,ntg,id)
!     
!     calculate the area
!     
          read(sideload(i)(2:2),'(i1)') ig
!     
!     number of nodes and integration points in the face
!     
          if(lakonl(4:4).eq.'2') then
            nope=20
            nopes=8
          elseif(lakonl(4:4).eq.'8') then
            nope=8
            nopes=4
          elseif(lakonl(4:5).eq.'10') then
            nope=10
            nopes=6
          elseif(lakonl(4:4).eq.'4') then
            nope=4
            nopes=3
          elseif(lakonl(4:5).eq.'15') then
            nope=15
          else
            nope=6
          endif
!     
          if(lakonl(4:5).eq.'8R') then
            mint2d=1
          elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R'))
     &           then
            if(lakonl(7:7).eq.'A') then
              mint2d=2
            else
              mint2d=4
            endif
          elseif(lakonl(4:4).eq.'2') then
            mint2d=9
          elseif(lakonl(4:5).eq.'10') then
            mint2d=3
          elseif(lakonl(4:4).eq.'4') then
            mint2d=1
          endif
!     
          if(lakonl(4:4).eq.'6') then
            mint2d=1
            if(ig.le.2) then
              nopes=3
            else
              nopes=4
            endif
          endif
          if(lakonl(4:5).eq.'15') then
            if(ig.le.2) then
              mint2d=3
              nopes=6
            else
              mint2d=4
              nopes=8
            endif
          endif
!     
!     connectivity of the element
!     
          do k=1,nope
            konl(k)=kon(index+k)
          enddo
!     
!     coordinates of the nodes belonging to the face
!     
          if((nope.eq.20).or.(nope.eq.8)) then
            do k=1,nopes
              tl2(k)=v(0,konl(ifaceq(k,ig)))
              do j=1,3
                xl2(j,k)=co(j,konl(ifaceq(k,ig)))+
     &               v(j,konl(ifaceq(k,ig)))
              enddo
            enddo
          elseif((nope.eq.10).or.(nope.eq.4)) then
            do k=1,nopes
              tl2(k)=v(0,konl(ifacet(k,ig)))
              do j=1,3
                xl2(j,k)=co(j,konl(ifacet(k,ig)))+
     &               v(j,konl(ifacet(k,ig)))
              enddo
            enddo
          else
            do k=1,nopes
              tl2(k)=v(0,konl(ifacew(k,ig)))
              do j=1,3
                xl2(j,k)=co(j,konl(ifacew(k,ig)))+
     &               v(j,konl(ifacew(k,ig)))
              enddo
            enddo
          endif
!     
!     integration to obtain the area and the mean
!     temperature
!     
          do m=1,mint2d
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
              xi=gauss2d1(1,m)
              et=gauss2d1(2,m)
              weight=weight2d1(m)
            elseif((lakonl(4:4).eq.'8').or.
     &             (lakonl(4:6).eq.'20R').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
              xi=gauss2d2(1,m)
              et=gauss2d2(2,m)
              weight=weight2d2(m)
            elseif(lakonl(4:4).eq.'2') then
              xi=gauss2d3(1,m)
              et=gauss2d3(2,m)
              weight=weight2d3(m)
            elseif((lakonl(4:5).eq.'10').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
              xi=gauss2d5(1,m)
              et=gauss2d5(2,m)
              weight=weight2d5(m)
            elseif((lakonl(4:4).eq.'4').or.
     &             ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
              xi=gauss2d4(1,m)
              et=gauss2d4(2,m)
              weight=weight2d4(m)
            endif
!     
            if(nopes.eq.8) then
              call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4) then
              call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
              call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            else
              call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            endif
!     
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))
            areaj=dxsj2*weight
!     
            temp=0.d0
            do k=1,3
              coords(k)=0.d0
            enddo
            do j=1,nopes
              temp=temp+tl2(j)*shp2(4,j)
              do k=1,3
                coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
              enddo
            enddo
!     
            sinktemp=v(0,node)
            if(sideload(i)(5:6).ne.'NU') then
              h(1)=xloadact(1,i)
            else
              read(sideload(i)(2:2),'(i1)') jltyp
              jltyp=jltyp+10
              call film(h,sinktemp,temp,istep,
     &             iinc,tvar,nelem,m,coords,jltyp,field,nfield,
     &             sideload(i),node,areaj,v,mi,ipkon,kon,lakon,
     &             iponoel,inoel,ielprop,prop,ielmat,shcon,nshcon,
     &             rhcon,nrhcon,ntmat_,cocon,ncocon,
     &             ipobody,xbodyact,ibody,heatnod,heatfac)
            endif
            if(lakonl(5:7).eq.'0RA') then
              term=2.d0*((temp-sinktemp)*h(1)+heatnod)*dxsj2*weight
              bc(ieq)=bc(ieq)+term
              qat=qat+dabs(term)
              nalt=nalt+1
            else
              term=((temp-sinktemp)*h(1)+heatnod)*dxsj2*weight
              bc(ieq)=bc(ieq)+term
              qat=qat+dabs(term)
              nalt=nalt+1
            endif
          enddo
        endif
      enddo
!     
!     prescribed heat generation: contribution to the energy equations        
!     
      do i=1,ntg
        node=itg(i)
        idof=8*(node-1)
        call nident(ikforc,idof,nforc,id)
        if(id.gt.0) then
          if(ikforc(id).eq.idof) then
            ieq=nacteq(0,node)
            if(ieq.ne.0) then
              term=xforcact(ilforc(id))
              bc(ieq)=bc(ieq)+term
              qat=qat+dabs(term)
              nalt=nalt+1
            endif
            cycle
          endif
        endif
      enddo
!     
!     in the case of forced vortices, when temperature change 
!     is required, additional heat input is added in the energy 
!     equation for the downstream node
!     
      do i=1,nflow
        nelem=ieg(i)
!     
!     free vortex and no temperature change
!     
        if(lakon(nelem)(2:3).eq.'VO') then
          if(lakon(nelem)(4:5).eq.'FR') then
            if((nint(prop(ielprop(nelem)+8)).eq.0).or.
     &           (nint(prop(ielprop(nelem)+8)).eq.1)) cycle
          elseif(lakon(nelem)(4:5).eq.'FO') then
            if(nint(prop(ielprop(nelem)+6)).eq.0) cycle
          endif
!     
          call calcheatnet(nelem,lakon,ipkon,kon,v,ielprop,prop,
     &         ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,ipobody,ibody,
     &         xbodyact,mi,nacteq,bc,qat,nalt)
!     
        elseif(lakon(nelem)(2:5).eq.'GAPR') then
!     
          call calcheatnet(nelem,lakon,ipkon,kon,v,ielprop,prop,
     &         ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,ipobody,ibody,
     &         xbodyact,mi,nacteq,bc,qat,nalt)
!     
        elseif(lakon(nelem)(2:2).eq.'U') then
          call calcheatnet(nelem,lakon,ipkon,kon,v,ielprop,prop,
     &         ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,ipobody,ibody,
     &         xbodyact,mi,nacteq,bc,qat,nalt)
        else
          cycle
        endif
      enddo
!     
!     transfer element ABSOLUTE TO RELATIVE / RELATIVE TO ABSOLUTE
!     
      do i=1,nflow
        nelem=ieg(i)
!     
        if((lakon(nelem)(2:4).eq.'ATR').or.
     &       (lakon(nelem)(2:4).eq.'RTA'))  then
!     
          nodem=kon(ipkon(nelem)+2)
          xflow=v(1,nodem)
          if(xflow.gt.0d0) then
            node1=kon(ipkon(nelem)+1)
            node2=kon(ipkon(nelem)+3)
          else
            node1=kon(ipkon(nelem)+1)
            node2=kon(ipkon(nelem)+3)
          endif
!     
!     computing temperature corrected Cp=Cp(T) coefficient 
!     
          Tg1=v(0,node1)
          Tg2=v(0,node2)
          if((lakon(nelem)(2:3).ne.'LP').and.
     &         (lakon(nelem)(2:3).ne.'LI')) then
            gastemp=(tg1+tg2)/2.d0
          else
            if(xflow.gt.0) then
              gastemp=tg1
            else
              gastemp=tg2
            endif
          endif
!     
          imat=ielmat(1,nelem)
          call materialdata_tg(imat,ntmat_,gastemp,
     &         shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
!     
          call cp_corrected(cp,Tg1,Tg2,cp_cor)
!     
          index=ielprop(nelem)
          U=prop(index+1)
          ct=prop(index+2)
!     
!     if a swirl element was given by the user, this takes
!     precendence to a value of ct
!     
          if(nint(prop(index+3)).ne.0) then
            nelemswirl=nint(prop(index+3))
            index2=ielprop(nelemswirl)
!     
!     previous element is a preswirl nozzle
!     
            if(lakon(nelemswirl)(2:5).eq.'ORPN') then
              ct=prop(index2+4)
!     
!     previous element is a forced vortex
!     
            elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
              ct=prop(index2+7)
!     
!     previous element is a free vortex
!     
            elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
              ct=prop(index2+9)
            endif
          endif                
!     
          if(lakon(nelem)(2:4).eq.'ATR') then
            heat=Cp/Cp_cor*(0.5d0*(U**2-2d0*U*Ct)*xflow)
!     
          elseif(lakon(nelem)(2:4).eq.'RTA') then
            heat=Cp/Cp_cor*(-0.5d0*(U**2-2d0*U*Ct)*xflow)
          endif
!     
!     including the resulting additional heat flux in the energy equation
!     
          if(xflow.gt.0d0)then
            ieq=nacteq(0,node2)
            if(ieq.ne.0) then
              bc(ieq)=bc(ieq)+heat
              qat=qat+dabs(heat)
              nalt=nalt+1
            endif
          else
            ieq=nacteq(0,node1)
            if(ieq.ne.0) then
              bc(ieq)=bc(ieq)+heat
              qat=qat+dabs(heat)
              nalt=nalt+1
            endif
          endif
        endif
      enddo 
!     
!     additional multiple point constraints
!     
      j=nteq+1
      do i=nmpc,1,-1
        if(labmpc(i)(1:7).ne.'NETWORK') cycle
        j=j-1
        if(labmpc(i)(8:8).eq.' ') then
!     
!     linear equation
!     
          index=ipompc(i)
!     
          do
            node=nodempc(1,index)
            idir=nodempc(2,index)
            bc(j)=bc(j)-v(idir,node)*coefmpc(index)
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        else
!     
!     user-defined network equation
!     
          call networkmpc_rhs(i,ipompc,nodempc,coefmpc,labmpc,
     &         v,bc,j,mi,ipkon,kon,lakon,iponoel,
     &         inoel,ielprop,prop,ielmat,
     &         shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
        endif
      enddo
!     
!     determining the typical energy flow
!     
      if(nalt.gt.0) qat=qat/nalt
!     
!     determining the typical mass flow
!     
      if(nalf.gt.0) qaf=qaf/nalf
!     
!     max. residual energy flow, residual mass flow
!     and residuals from the element equations
!     
      ramt=0.d0
      ramf=0.d0
      ramp=0.d0
      do i=1,ntg
        node=itg(i)
        if(nacteq(0,node).ne.0) then
          ramt=max(ramt,dabs(bc(nacteq(0,node))))
        endif
        if(nacteq(1,node).ne.0) then
          ramf=max(ramf,dabs(bc(nacteq(1,node))))
        endif
        if(nacteq(2,node).ne.0) then
          ramp=max(ramp,dabs(bc(nacteq(2,node))))
        endif
      enddo
!     
c     write(*,*) 'bc in resultnet'
c     do i=1,3
c     write(*,'(1x,e11.4)') bc(i)
c     enddo
!     
      return
      end
