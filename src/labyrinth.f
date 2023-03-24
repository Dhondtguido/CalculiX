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
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine labyrinth(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,co,dvi,numf,vold,set,
     &     kon,ipkon,mi,ttime,time,iaxial,iplausi)
!     
!     labyrinth element
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),index,kflag,mi(*),
     &     inv,kgas,n,iaxial,nodea,nodeb,ipkon(*),kon(*),i,
     &     itype,iplausi
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),kappa,R,a,d,
     &     p1,p2,T1,Aeff,C1,C2,C3,cd,cp,physcon(*),p2p1,km1,dvi,
     &     kp1,kdkm1,tdkp1,km1dk,x,y,ca1,cb1,ca2,cb2,dT1,alambda,
     &     rad,reynolds,pi,ppkrit,co(3,*),ttime,time,
     &     carry_over,dlc,hst,e,szt,num,denom,t,s,b,h,cdu,
     &     cd_radius,cst,dh,cd_honeycomb,cd_lab,bdh,
     &     pt0zps1,cd_1spike,cdbragg,rzdh,
     &     cd_correction,p1p2,xflow_oil,T2,vold(0:mi(2),*)
!    
!
!
      itype=1
      pi=4.d0*datan(1.d0)
      e=2.718281828459045d0
!
      index=ielprop(nelem)
!    
      if(kflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif(kflag.eq.1)then
         if(v(1,nodem).ne.0.d0) then
            xflow=v(1,nodem)
            return
         endif
!     
          kappa=(cp/(cp-R))
!
!     Usual Labyrinth
!     
          if(lakon(nelem)(2:5).ne.'LABF') then
             t=prop(index+1)
             s=prop(index+2)
             d=prop(index+4)
             n=nint(prop(index+5))
             b=prop(index+6)
             h=prop(index+7)
             dlc=prop(index+8)
             rad=prop(index+9)
             X=prop(index+10)
             Hst=prop(index+11)
!
             A=pi*D*s
!
!    "flexible" labyrinth for thermomechanical coupling
!
          elseif(lakon(nelem)(2:5).eq.'LABF') then
             nodea=nint(prop(index+1))
             nodeb=nint(prop(index+2))
             t=prop(index+4)
             d=prop(index+5)
             n=nint(prop(index+6))
             b=prop(index+7)
             h=prop(index+8)
             dlc=prop(index+9)
             rad=prop(index+10)
             X=prop(index+11)
             Hst=prop(index+12)

!
!     gap definition
             s=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &            co(1,nodea)-vold(1,nodea))**2)
             a=pi*d*s
          endif
!     
         p1=v(2,node1)
         p2=v(2,node2)
         if(p1.ge.p2) then
            inv=1
            T1=v(0,node1)-physcon(1)
         else
            inv=-1
            p1=v(2,node2)
            p2=v(2,node1)
            T1=v(0,node2)-physcon(1)
         endif
!     
         cd=1.d0
         Aeff=A*cd
         p2p1=p2/p1
!     
!************************
!     one fin 
!*************************
         if(n.eq.1) then
!     
            km1=kappa-1.d0
            kp1=kappa+1.d0
            kdkm1=kappa/km1
            tdkp1=2.d0/kp1
            C2=tdkp1**kdkm1
!     
!     subcritical
!     
            if(p2p1.gt.C2) then
               xflow=inv*p1*Aeff*dsqrt(2.d0*kdkm1*p2p1**(2.d0/kappa)
     &              *(1.d0-p2p1**(1.d0/kdkm1))/r)/dsqrt(T1)
!     
!     critical
!     
            else
               xflow=inv*p1*Aeff*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &              dsqrt(T1)
            endif
         endif
!     
!***********************
!     straight labyrinth and stepped labyrinth
!     method found in "Air system Correlations Part1 Labyrinth Seals"
!     H.Zimmermann and K.H. Wolff
!     ASME 98-GT-206
!**********************
!     
         if(n.ge.2) then
!     
            call lab_straight_ppkrit(n,ppkrit)
!     
!     subcritical case
!     
            if(p2p1.gt.ppkrit) then
               xflow=inv*p1*Aeff/dsqrt(T1)*dsqrt((1.d0-p2p1**2.d0)
     &              /(R*(n-log(p2p1)/log(e))))
!     
!     critical case
!     
            else
               xflow=inv*p1*Aeff/dsqrt(T1)*dsqrt(2.d0/R)*ppkrit
            endif
         endif
!
      elseif(kflag.eq.2)then
         numf=4
         alambda=10000.d0
!     
         p1=v(2,node1)
         p2=v(2,node2)
         if(p1.ge.p2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            T1=v(0,node1)-physcon(1)
            T2=v(0,node2)-physcon(1)
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
         else
            inv=-1
            p1=v(2,node2)
            p2=v(2,node1)
            xflow=-v(1,nodem)*iaxial
            T1=v(0,node2)-physcon(1)
            T2=v(0,node1)-physcon(1)
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
         endif
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
!     Usual labyrinth
!
         if(lakon(nelem)(2:5).ne. 'LABF') then
            kappa=(cp/(cp-R))
            t=prop(index+1)
            s=prop(index+2)
            d=prop(index+4)
            n=nint(prop(index+5))
            b=prop(index+6)
            h=prop(index+7)
            dlc=prop(index+8)
            rad=prop(index+9)
            X=prop(index+10)
            Hst=prop(index+11)
            A=pi*D*s
!
!     Flexible labyrinth for coupled calculations
!
         elseif(lakon(nelem)(2:5).eq.'LABF') then
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
c            iaxial=nint(prop(index+3))
            t=prop(index+4)
            d=prop(index+5)
            n=nint(prop(index+6))
            b=prop(index+7)
            h=prop(index+8)
            dlc=prop(index+9)
            rad=prop(index+10)
            X=prop(index+11)
            Hst=prop(index+12)
!     
!     gap definition
             s=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &            co(1,nodea)-vold(1,nodea))**2)
             a=pi*d*s
          endif
!     
         p2p1=p2/p1
         dT1=dsqrt(T1)
!     
         Aeff=A
!     
!     honeycomb stator correction
!     
         cd_honeycomb=1.d0
         if(dlc.ne.0.d0)then
            call cd_lab_honeycomb(s,dlc,cd_honeycomb)
            cd_honeycomb=1+cd_honeycomb/100
         endif
!     
!     inlet radius correction
!     
         cd_radius=1.d0
         if((rad.ne.0.d0).and.(n.ne.1d0)) then
            call cd_lab_radius(rad,s,Hst,cd_radius)
         endif
!     
!     carry over factor (only for straight throught labyrinth)
!     
         if((n.ge.2).and.(hst.eq.0.d0)) then
            cst=n/(n-1.d0)
            szt=s/t
            carry_over=cst/dsqrt(cst-szt/(szt+0.02d0))
            Aeff=Aeff*carry_over
         endif
!     
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in labyrinth: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif     
!     
!     calculation of the number of reynolds for a gap
!     
         reynolds=dabs(xflow)*2.d0*s/(dvi*A*cd_honeycomb/cd_radius)
!     
!**************************************
!     single fin labyrinth 
!     the resolution procedure is the same as for the restrictor
!**************************************
!     
         if(n.eq.1)then
!     
!     single fin labyrinth
!     
!     incompressible basis cd , reynolds correction,and radius correction
!     
!     "Flow Characteristics of long orifices with rotation and corner radiusing"
!     W.F. Mcgreehan and M.J. Schotsch
!     ASME 87-GT-162
!     
            dh=2*s
            bdh=b/dh
            rzdh=rad/dh
!     
            call cd_Mcgreehan_Schotsch(rzdh,bdh,reynolds,cdu) 
!     
!     compressibility correction factor
!     
!     S.L.Bragg
!     "Effect of conpressibility on the discharge coefficient of orifices and convergent nozzles"
!     Journal of Mechanical engineering vol 2 No 1 1960
!     
            call cd_bragg(cdu,p2p1,cdbragg,itype)
            cd=cdbragg
            Aeff=Aeff*cd 
!     
            km1=kappa-1.d0
            kp1=kappa+1.d0
            kdkm1=kappa/km1
            tdkp1=2.d0/kp1
            C2=tdkp1**kdkm1
!     
            if(p2p1.gt.C2) then
               C1=dsqrt(2.d0*kdkm1/r)*Aeff
               km1dk=1.d0/kdkm1
               y=p2p1**km1dk
               x=dsqrt(1.d0-y)
               ca1=-C1*x/(kappa*p1*y)
               cb1=C1*km1dk/(2.d0*p1)
               ca2=-ca1*p2p1-xflow*dT1/(p1*p1)
               cb2=-cb1*p2p1
               f=xflow*dT1/p1-C1*p2p1**(1.d0/kappa)*x
               if(cb2.le.-(alambda+ca2)*x) then
                  df(1)=-alambda
               elseif(cb2.ge.(alambda-ca2)*x) then
                  df(1)=alambda
               else
                  df(1)=ca2+cb2/x
               endif
               df(2)=xflow/(2.d0*p1*dT1)
               df(3)=inv*dT1/p1
               if(cb1.le.-(alambda+ca1)*x) then
                  df(4)=-alambda
               elseif(cb1.ge.(alambda-ca1)*x) then
                  df(4)=alambda
               else
                  df(4)=ca1+cb1/x
               endif
            else
               C3=dsqrt(kappa/r)*(tdkp1)**(kp1/(2.d0*km1))*Aeff
               f=xflow*dT1/p1-C3
               df(1)=-xflow*dT1/(p1)**2
               df(2)=xflow/(2*p1*dT1)
               df(3)=inv*dT1/p1
               df(4)=0.d0
            endif
         endif
!     
!****************************************
!     straight labyrinth & stepped labyrinth
!     method found in "Air system Correlations Part1 Labyrinth Seals"
!     H.Zimmermann and K.H. Wolff
!     ASME 98-GT-206
!****************************************
!     
         if(n.ge.2) then
            num=(1.d0-p2p1**2)
            denom=R*(n-log(p2p1)/log(e))
!     
!     straight labyrinth
!     
            if((hst.eq.0.d0).and.(n.ne.1)) then
               call cd_lab_straight(n,p2p1,s,b,reynolds,cd_lab)
               Aeff=Aeff*cd_lab*cd_honeycomb*cd_radius
!     
!     Stepped Labyrinth
!     
            else 
!     corrective term for the first spike
               p1p2=p1/p2
               pt0zps1=(p1p2)**(1/prop(index+4))
               call cd_lab_1spike(pt0zps1,s,b,cd_1spike)
!     
!     corrective term for cd_lab_1spike
!     
               call cd_lab_correction(p1p2,s,b,cd_correction)
!     
!     calculation of the discharge coefficient of the stepped labyrinth
!     
               cd=cd_1spike*cd_correction
               cd_lab=cd
!     
               Aeff=Aeff*cd_lab*cd_radius*cd_honeycomb
            endif
!     
            call lab_straight_ppkrit(n,ppkrit)
!     
!     subcritical case
!     
            if(p2p1.gt.ppkrit) then
!     
               f=xflow*dT1/p1-dsqrt(num/denom)*Aeff
!     
               df(1)=xflow*dt1/p1**2.d0-Aeff/2.d0
     &              *dsqrt(denom/num)*(2.d0*(p2**2.d0/p1**3.d0)/denom)
     &              +num/denom**2.d0*r/p1
               df(2)=xflow/(2.d0*p1*dT1)
               df(3)=inv*dT1/p1
               df(4)=-Aeff/2.d0*dsqrt(denom/num)*(-2.d0*(p2/p1**2.d0)
     &              /denom)+num/denom**2.d0*r/p2
!     
!     critical case
!     
            else
               C2=dsqrt(2/R)*Aeff*ppkrit
!     
               f=xflow*dT1/p1-C2
               df(1)=-xflow*dT1/(p1**2)
               df(2)=xflow/(2.d0*p1*dT1)
               df(3)=inv*dT1/p1
               df(4)=0.d0
            endif
         endif
!
!     output
!
      elseif(kflag.eq.3)then
!
         p1=v(2,node1)
         p2=v(2,node2)
         if(p1.ge.p2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            T1=v(0,node1)-physcon(1)
            T2=v(0,node2)-physcon(1)
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
         else
            inv=-1
            p1=v(2,node2)
            p2=v(2,node1)
            xflow=-v(1,nodem)*iaxial
            T1=v(0,node2)-physcon(1)
            T2=v(0,node2)-physcon(1)
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
         endif
!     
         kappa=(cp/(cp-R))
         t=prop(index+1)
         s=prop(index+2)
         d=prop(index+3)
         n=nint(prop(index+4))
         b=prop(index+5)
         h=prop(index+6)
         dlc=prop(index+7)
         rad=prop(index+8)
         X=prop(index+9)
         Hst=prop(index+10)
!     
         p2p1=p2/p1
         dT1=dsqrt(T1)
!     
         pi=4.d0*datan(1.d0)
         A=pi*D*s
         Aeff=A
         e=2.718281828459045d0
!     
!     honeycomb stator correction
!     
         if(dlc.ne.0.d0)then
            call cd_lab_honeycomb(s,dlc,cd_honeycomb)
            Aeff=Aeff*(1.d0+cd_honeycomb/100.d0)
         else
            cd_honeycomb=0
         endif
!     
!     inlet radius correction
!     
         if((rad.ne.0.d0).and.(n.ne.1d0)) then
            call cd_lab_radius(rad,s,Hst,cd_radius)
            Aeff=Aeff*cd_radius
         else
            cd_radius=1
         endif
!     
!     carry over factor (only for straight throught labyrinth)
!     
         if((n.gt.1).and.(hst.eq.0.d0)) then
            cst=n/(n-1.d0)
            szt=s/t
            carry_over=cst/dsqrt(cst-szt/(szt+0.02d0))
            Aeff=Aeff*carry_over
         endif
!     
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in labyrinth: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif     
!     
!     calculation of the number of reynolds for a gap
!     
         reynolds=dabs(xflow)*2.d0*s/(dvi*A)
!**************************************
!     single fin labyrinth 
!     the resolution procedure is the same as for the restrictor
!**************************************
!     
         if(n.eq.1)then
!     
!     single fin labyrinth
!     
!     incompressible basis cd , reynolds correction,and radius correction
!     
!     "Flow Characteristics of long orifices with rotation and corner radiusing"
!     W.F. Mcgreehan and M.J. Schotsch
!     ASME 87-GT-162
!     
            dh=2*s
            bdh=b/dh
            rzdh=rad/dh
!     
            call cd_Mcgreehan_Schotsch(rzdh,bdh,reynolds,cdu) 
!     
!     compressibility correction factor
!     
!     S.L.Bragg
!     "Effect of conpressibility on the discharge coefficient of orifices and convergent nozzles"
!     Journal of Mechanical engineering vol 2 No 1 1960
!     
            call cd_bragg(cdu,p2p1,cdbragg,itype)
            cd=cdbragg
            Aeff=Aeff*cd 
         endif
!     
!****************************************
!     straight labyrinth & stepped labyrinth
!     method found in "Air system Correlations Part1 Labyrinth Seals"
!     H.Zimmermann and K.H. Wolff
!     ASME 98-GT-206
!****************************************
!     
         if(n.ge.2) then
            num=(1.d0-p2p1**2)
            denom=R*(n-log(p2p1)/log(e))
!     
!     straight labyrinth
!     
            if((hst.eq.0.d0).and.(n.ne.1)) then
               call cd_lab_straight(n,p2p1,s,b,reynolds,cd_lab)
               Aeff=Aeff*cd_lab*cd_honeycomb*cd_radius
!     
!     Stepped Labyrinth
!     
            else 
!     corrective term for the first spike
               p1p2=p1/p2
               pt0zps1=(p1p2)**(1/prop(index+4))
               call cd_lab_1spike(pt0zps1,s,b,cd_1spike)
!     
!     corrective term for cd_lab_1spike
!     
               call cd_lab_correction(p1p2,s,b,cd_correction)
!     
!     calculation of the discharge coefficient of the stepped labyrinth
!     
               cd=cd_1spike*cd_correction
               cd_lab=cd
!     
               Aeff=Aeff*cd_lab*cd_radius*cd_honeycomb
            endif
!     
            call lab_straight_ppkrit(n,ppkrit)
!
         endif
!
         xflow_oil=0
!
         write(1,*) ''
         write(1,55) ' from node',node1,
     &' to node', node2,':   air massflow rate= ',xflow,
     &', oil massflow rate= ',xflow_oil
 55      FORMAT(1X,A,I6,A,I6,A,e11.4,A,A,e11.4,A)
         
         if(inv.eq.1) then
          write(1,56)'       Inlet node ',node1,':   Tt1=',T1,
     &           ', Ts1=',T1,', Pt1=',p1

            write(1,*)'             Element ',nelem,lakon(nelem)
            write(1,57)'             dyn.visc.= ',dvi,', Re= ' ,
     &           reynolds,
     &', Cd_radius= ',cd_radius,', Cd_honeycomb= ', 1+cd_honeycomb/100
!
!     straight labyrinth
           if((hst.eq.0.d0).and.(n.ne.1)) then
              write(1,58)'             COF= ',carry_over,
     &             ', Cd_lab= ',cd_lab,', Cd= ',carry_over*cd_lab

!     stepped labyrinth
           elseif(hst.ne.0d0) then
              write(1,59)'             Cd_1_fin= ',
     &             cd_1spike, ', Cd= ',cd,', pt0/ps1= ',pt0zps1,
     &             ', p0/pn= ',p1/p2

!     single fin labyrinth
           elseif(n.eq.1) then
              write(1,60) '             Cd_Mcgreehan= ',cdu,
     &             ', Cd= ',cdbragg
           endif
                 
            write(1,56)'      Outlet node ',node2,':   Tt2= ',T2,
     &           ', Ts2= ',T2,', Pt2= ',p2

!     
         else if(inv.eq.-1) then
            write(1,56)'       Inlet node ',node2,':    Tt1= ',T1,
     &           ', Ts1= ',T1,', Pt1= ',p1 
         
            write(1,*)'             element ',nelem,lakon(nelem)
            write(1,57)'             dyn.visc.=',dvi,', Re= '
     &           ,reynolds,
     & ', Cd_radius= ',cd_radius,', Cd_honeycomb= ',1+cd_honeycomb/100
!
!     straight labyrinth
            if((hst.eq.0.d0).and.(n.ne.1)) then
               write(1,58)'                  COF = ',carry_over,
     &              ', Cd_lab= ',cd_lab,', Cd= ',carry_over*cd_lab
!
!     stepped labyrinth
            elseif(hst.ne.0d0) then
               write(1,59)'                 Cd_1_fin= ',
     &              cd_1spike,', Cd= ',cd,', pt0/ps1= ',pt0zps1,
     &             ', p0/pn= ',p1/p2

!     single fin labyrinth
            elseif(n.eq.1) then
               write(1,60) '              Cd_Mcgreehan= ',
     & cdu,' Cd= ',cdbragg
           endif
           write(1,56)'      Outlet node ',node1,':    Tt2= ',T2,
     &          ', Ts2= ',T2,', Pt2= ',p2

        endif
!         
 56      FORMAT(1X,A,I6,A,e11.4,A,e11.4,A,e11.4,A)
 57      FORMAT(1X,A,E11.5,A,e11.4,A,e11.4,A,e11.4)
 58      FORMAT(1X,A,e11.4,A,e11.4,A,e11.4)
 59      FORMAT(1X,A,e11.4,A,e11.4,A,e11.4,A,e11.4)
 60      FORMAT(1X,A,e11.4,A,e11.4)
      endif
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!         
      return
      end
