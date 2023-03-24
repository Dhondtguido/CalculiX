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
      subroutine orifice(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
!     
!     orifice element
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
     &     ielprop(*),nodef(*),idirf(*),index,kflag,
     &     inv,ipkon(*),kon(*),number,kgas,nelemswirl,
     &     nodea,nodeb,iaxial,mi(*),i,itype,iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),kappa,R,a,d,dl,
     &     p1,p2,T1,Aeff,C1,C2,C3,cd,cp,physcon(*),p2p1,km1,dvi,
     &     kp1,kdkm1,tdkp1,km1dk,x,y,ca1,cb1,ca2,cb2,dT1,alambda,
     &     rad,beta,reynolds,theta,k_phi,c2u_new,u,pi,xflow_oil,
     &     ps1pt1,uref,cd_chamf,angle,vid,cdcrit,T2,radius,
     &     initial_radius,co(3,*),vold(0:mi(2),*),offset,ttime,time,
     &     x_tab(100), y_tab(100),x_tab2(100),y_tab2(100),curve,xmach
!
!
!
      pi=4.d0*datan(1.d0)   
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
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         a=prop(index+1)
         d=prop(index+2)
         dl=prop(index+3)
!     
         if(lakon(nelem)(2:5).eq.'ORFL') then
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            offset=prop(index+4)
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)-offset
            initial_radius=dsqrt((co(1,nodeb)-co(1,nodea))**2)-offset
               A=pi*radius**2
            d=2*radius
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
!     
         p2p1=p2/p1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
         Aeff=A*cd
         if(p2p1.gt.C2) then
            xflow=inv*p1*Aeff*dsqrt(2.d0*kdkm1*p2p1**(2.d0/kappa)
     &           *(1.d0-p2p1**(1.d0/kdkm1))/r)/dsqrt(T1)
         else
            xflow=inv*p1*Aeff*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &           dsqrt(T1)
         endif
!     
      elseif(kflag.eq.2)then
!     
         numf=4
         alambda=10000.d0
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         a=prop(index+1)
!     
         p1=v(2,node1)
         p2=v(2,node2)
         if(p1.ge.p2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            T1=v(0,node1)-physcon(1)
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
!     calculation of the dynamic viscosity 
!     
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in orifice: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif 
!     
         if((lakon(nelem)(4:5).ne.'BT').and.
     &        (lakon(nelem)(4:5).ne.'PN').and.
     &        (lakon(nelem)(4:5).ne.'C1').and.
     &        (lakon(nelem)(4:5).ne.'FL') ) then
            d=prop(index+2)
            dl=prop(index+3)
!     circumferential velocity of the rotating hole (same as disc @ given radius)
            u=prop(index+7)
            nelemswirl=nint(prop(index+8))
            if(nelemswirl.eq.0) then
               uref=0.d0
            else
!     swirl generating element
!     
!     preswirl nozzle
               if(lakon(nelemswirl)(2:5).eq.'ORPN') then
                  uref=prop(ielprop(nelemswirl)+5)
!     rotating orifices
               else if((lakon(nelemswirl)(2:5).eq.'ORMM').or.
     &                 (lakon(nelemswirl)(2:5).eq.'ORMA').or.
     &                 (lakon(nelemswirl)(2:5).eq.'ORPM').or.
     &                 (lakon(nelemswirl)(2:5).eq.'ORPA')) then
                  uref=prop(ielprop(nelemswirl)+7)
!     forced vortex
               elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
                  uref=prop(ielprop(nelemswirl)+7)
!     free vortex 
               elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
                  uref=prop(ielprop(nelemswirl)+9)
!     Moehring 
               elseif(lakon(nelemswirl)(2:4).eq.'MRG') then
                  uref=prop(ielprop(nelemswirl)+10)
!     RCAVO 
               elseif((lakon(nelemswirl)(2:4).eq.'ROR').or.
     &                 (lakon(nelemswirl)(2:4).eq.'ROA'))then
                  uref=prop(ielprop(nelemswirl)+6)
!     RCAVI 
               elseif(lakon(nelemswirl)(2:4).eq.'RCV') then
                  uref=prop(ielprop(nelemswirl)+5)
!     
               else
                  write(*,*) '*ERROR in orifice:'
                  write(*,*) ' element',nelemswirl
                  write(*,*) ' refered by element',nelem
                  write(*,*) ' is not a swirl generating element'
               endif
            endif
!     write(*,*) 'nelem',nelem, u, uref
            u=u-uref
            angle=prop(index+5)
!     
         endif
!     
!     calculate the discharge coefficient using Bragg's Method
!     "Effect of Compressibility on the discharge coefficient 
!     of orifices and convergent nozzles"   
!     journal of mechanical Engineering
!     vol2 No 1 1960
!     
         if((lakon(nelem)(2:5).eq.'ORBG')) then
!     
            p2p1=p2/p1
            cdcrit=prop(index+2)
!     
            itype=2
            call cd_bragg(cdcrit,p2p1,cd,itype)
!     
         elseif(lakon(nelem)(2:5).eq.'ORMA') then
!     
!     calculate the discharge coefficient using own table data and 
!     using Dr.Albers method for rotating cavities
!     
            call cd_own_albers(p1,p2,dl,d,cd,u,T1,R,kappa)
!     
!     outlet circumferential velocity of the fluid is equal to the circumferential velocity of the hole
!     as the holes are perpendicular to the rotating surface and rotating with it
!     prop(index+7)
!     
!     chamfer correction
!     
            if(angle.gt.0.d0)then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif(lakon(nelem)(2:5).eq.'ORMM') then
!     
!     calculate the discharge coefficient using McGreehan and Schotsch method
!     
            rad=prop(index+4)
!     
            reynolds=dabs(xflow)*d/(dvi*a)
!     
!     outlet circumferential velocity of the fluid is equal to the circumferential velocity of the hole
!     as the holes are perpendicular to the rotating surface and rotating with it
!     prop(index+7)
!     
            call cd_ms_ms(p1,p2,T1,rad,d,dl,kappa,r,reynolds,u,vid,cd)
!     
            if(cd.ge.1) then
               write(*,*) ''
               write(*,*) '**WARNING**'
               write(*,*) 'in RESTRICTOR ',nelem
               write(*,*) 'Calculation using' 
               write(*,*) ' McGreehan and Schotsch method:'
               write(*,*) ' Cd=',Cd,'>1 !'
               write(*,*) 'Calcultion will proceed will Cd=1'
               write(*,*) 'l/d=',dl/d,'r/d=',rad/d,'u/vid=',u/vid
               cd=1.d0
            endif
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif  (lakon(nelem)(2:5).eq.'ORPA') then
!     
!     calculate the discharge coefficient using Parker and Kercher method
!     and using Dr. Albers method for rotating cavities
!     
            rad=prop(index+4)
!     
            beta=prop(index+6)
!     
            reynolds=dabs(xflow)*d/(dvi*a)
!     
            call cd_pk_albers(rad,d,dl,reynolds,p2,p1,beta,kappa,
     &           cd,u,T1,R)
!     
!     outlet circumferential velocity of the fluid is equal to the circumferential velocity of the hole
!     as the holes are perpendicular to the rotating surface and rotating with it
!     prop(index+7)
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif(lakon(nelem)(2:5).eq.'ORPM') then
!     
!     calculate the discharge coefficient using Parker and Kercher method
!     and using Mac Grehan and Schotsch method for rotating cavities
!     
            rad=prop(index+4)
!     
            beta=prop(index+6)
            reynolds=dabs(xflow)*d/(dvi*a)
!     
            call cd_pk_ms(rad,d,dl,reynolds,p2,p1,beta,kappa,cd,
     &           u,T1,R)
!     
!     outlet circumferential velocity of the fluid is equal to the circumferential velocity of the hole
!     as the holes are perpendicular to the rotating surface and rotating with it
!     prop(index+7)
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif(lakon(nelem)(2:5).eq.'ORC1') then
!     
            d=dsqrt(a*4/Pi)
            reynolds=dabs(xflow)*d/(dvi*a)
            cd=1.d0
!     
         elseif(lakon(nelem)(2:5).eq.'ORBT') then
!     
!     calculate the discharge coefficient of bleed tappings (OWN tables)
!     
            ps1pt1=prop(index+2)
            curve=nint(prop(index+3))
            number=nint(prop(index+4))
!     
            if(number.ne.0.d0)then
               do i=1,number
                  x_tab(i)=prop(index+2*i+3)
                  y_tab(i)=prop(index+2*i+4)
               enddo
            endif
!     
            call cd_bleedtapping(p2,p1,ps1pt1,number,curve,x_tab,y_tab,
     &           cd)
!     
         elseif(lakon(nelem)(2:5).eq.'ORPN') then
!     
!     calculate the discharge coefficient of preswirl nozzle (OWN tables)
!     
            d=dsqrt(4*A/pi)
            reynolds=dabs(xflow)*d/(dvi*a)
            curve=nint(prop(index+4))
            number=nint(prop(index+6))
            if(number.ne.0.d0)then
               do i=1,number
                  x_tab2(i)=prop(index+2*i+5)
                  y_tab2(i)=prop(index+2*i+6)
               enddo
            endif
            call cd_preswirlnozzle(p2,p1,number,curve,x_tab2,y_tab2
     &           ,cd)
!     
            theta=prop(index+2)
            k_phi=prop(index+3)
!     
            if(p2/p1.gt.(2/(kappa+1.d0))**(kappa/(kappa-1.d0))) then
               c2u_new=k_phi*cd*sin(theta*Pi/180.d0)*r*
     &              dsqrt(2.d0*kappa/(r*(kappa-1)))*
     &              dsqrt(T1*(1.d0-(p2/p1)**((kappa-1)/kappa)))
!     
            else
               c2u_new=k_phi*cd*sin(theta*Pi/180.d0)*r*
     &              dsqrt(2.d0*kappa/(r*(kappa-1)))*
     &              dsqrt(T1*(1.d0-2/(kappa+1)))
            endif
            prop(index+5)=c2u_new
!     
         elseif(lakon(nelem)(2:5).eq.'ORFL') then
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
c            iaxial=nint(prop(index+3))
            offset=prop(index+4)
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)-offset
!     
            initial_radius=dsqrt((co(1,nodeb)-co(1,nodea))**2)-offset
!     
c            if(iaxial.ne.0) then
c               A=pi*radius**2/iaxial
c            else
               A=pi*radius**2
c            endif
            d=2*radius
            reynolds=dabs(xflow)*d/(dvi*a)
            cd=1.d0
!     
         endif
!     
         if(cd.gt.1.d0) then
            Write(*,*) '*WARNING:'
            Write(*,*) 'In RESTRICTOR',nelem
            write(*,*) 'Cd greater than 1'
            write (*,*) 'Calculation will proceed using Cd=1'
            cd=1.d0
         endif
!     
         p2p1=p2/p1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
         Aeff=A*cd
         dT1=dsqrt(T1)
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
!     
!     output
!     
      elseif(kflag.eq.3) then
!     
         pi=4.d0*datan(1.d0)
         p1=v(2,node1)
         p2=v(2,node2)
         if(p1.ge.p2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            T1=v(0,node1)-physcon(1)
            T2=v(0,node2)-physcon(1)
         else
            inv=-1
            p1=v(2,node2)
            p2=v(2,node1)
            xflow=-v(1,nodem)*iaxial 
            T1=v(0,node2)-physcon(1)
            T2=v(0,node1)-physcon(1)
         endif
!     
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in orifice: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif 
!     
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         a=prop(index+1)
!     
         if((lakon(nelem)(4:5).ne.'BT').and.
     &        (lakon(nelem)(4:5).ne.'PN').and.
     &        (lakon(nelem)(4:5).ne.'C1')) then
            d=prop(index+2)
            dl=prop(index+3)
            u=prop(index+7)
            nelemswirl=nint(prop(index+8))
            if(nelemswirl.eq.0) then
               uref=0.d0
            else
!     swirl generating element
!     
!     preswirl nozzle
               if(lakon(nelemswirl)(2:5).eq.'ORPN') then
                  uref=prop(ielprop(nelemswirl)+5)
!     rotating orifices
               else if((lakon(nelemswirl)(2:5).eq.'ORMM').or.
     &                 (lakon(nelemswirl)(2:5).eq.'ORMA').or.
     &                 (lakon(nelemswirl)(2:5).eq.'ORPM').or.
     &                 (lakon(nelemswirl)(2:5).eq.'ORPA')) then
                  uref=prop(ielprop(nelemswirl)+7)
!     forced vortex
               elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
                  uref=prop(ielprop(nelemswirl)+7)
!     
!     free vortex 
               elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
                  uref=prop(ielprop(nelemswirl)+9)
!     Moehring 
               elseif(lakon(nelemswirl)(2:4).eq.'MRG') then
                  uref=prop(ielprop(nelemswirl)+10)
!     RCAVO 
               elseif((lakon(nelemswirl)(2:4).eq.'ROR').or.
     &                 (lakon(nelemswirl)(2:4).eq.'ROA'))then
                  uref=prop(ielprop(nelemswirl)+6)
!     RCAVI 
               elseif(lakon(nelemswirl)(2:4).eq.'RCV') then
                  uref=prop(ielprop(nelemswirl)+5)
               else
                  write(*,*) '*ERROR in orifice:'
                  write(*,*) ' element',nelemswirl
                  write(*,*) 'refered by element',nelem
                  write(*,*) 'is not a swirl generating element'
               endif
            endif
!     write(*,*) 'nelem',nelem, u, uref
            u=u-uref
            angle=prop(index+5)
!     
         endif
!     
!     calculate the discharge coefficient using Bragg's Method
!     "Effect of Compressibility on the discharge coefficient 
!     of orifices and convergent nozzles"   
!     journal of mechanical Engineering
!     vol2 No 1 1960
!     
         if((lakon(nelem)(2:5).eq.'ORBG')) then
!     
            p2p1=p2/p1
            d=dsqrt(a*4/Pi)           
            reynolds=dabs(xflow)*d/(dvi*a)
            cdcrit=prop(index+2)
!     
            itype=2
            call cd_bragg(cdcrit,p2p1,cd,itype)
!     
         elseif(lakon(nelem)(2:5).eq.'ORMA') then
!     
!     calculate the discharge coefficient using own table data and 
!     using Dr.Albers method for rotating cavities
!     
            reynolds=dabs(xflow)*d/(dvi*a)
!     
            call cd_own_albers(p1,p2,dl,d,cd,u,T1,R,kappa)
!     
!     chamfer correction
!     
            if(angle.gt.0.d0)then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif(lakon(nelem)(2:5).eq.'ORMM') then
!     
!     calculate the discharge coefficient using McGreehan and Schotsch method
!     
            rad=prop(index+4)
!     
            reynolds=dabs(xflow)*d/(dvi*a)
!     
            call cd_ms_ms(p1,p2,T1,rad,d,dl,kappa,r,reynolds,u,vid,cd)
!     
            if(cd.ge.1) then
               write(*,*) ''
               write(*,*) '**WARNING**'
               write(*,*) 'in RESTRICTOR ',nelem
               write(*,*) 'Calculation using' 
               write(*,*) ' McGreehan and Schotsch method:'
               write(*,*) ' Cd=',Cd,'>1 !'
               write(*,*) 'Calcultion will proceed will Cd=1'
               write(*,*) 'l/d=',dl/d,'r/d=',rad/d,'u/vid=',u/vid
               cd=1.d0
            endif
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif  (lakon(nelem)(2:5).eq.'ORPA') then
!     
!     calculate the discharge coefficient using Parker and Kercher method
!     and using Dr. Albers method for rotating cavities
!     
            rad=prop(index+4)
!     
            beta=prop(index+6)
!     
            reynolds=dabs(xflow)*d/(dvi*a)
!     
            call cd_pk_albers(rad,d,dl,reynolds,p2,p1,beta,kappa,
     &           cd,u,T1,R)
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif(lakon(nelem)(2:5).eq.'ORPM') then
!     
!     calculate the discharge coefficient using Parker and Kercher method
!     and using Mac Grehan and Schotsch method for rotating cavities
!     
            rad=prop(index+4)
!     
            beta=prop(index+6)
            reynolds=dabs(xflow)*d/(dvi*a)
!     
            call cd_pk_ms(rad,d,dl,reynolds,p2,p1,beta,kappa,cd,
     &           u,T1,R)
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(dl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif(lakon(nelem)(2:5).eq.'ORC1') then
!     
            d=dsqrt(a*4/Pi)
            reynolds=dabs(xflow)*d/(dvi*a)
            cd=1.d0
!     
         elseif(lakon(nelem)(2:5).eq.'ORBT') then
!     
!     calculate the discharge coefficient of bleed tappings (OWN tables)
!     
            d=dsqrt(A*Pi/4)
            reynolds=dabs(xflow)*d/(dvi*a)
            ps1pt1=prop(index+2)
            curve=nint(prop(index+3))
            number=nint(prop(index+4))
            reynolds=dabs(xflow)*d/(dvi*a)
            if(number.ne.0.d0)then
               do i=1,number
                  x_tab(i)=prop(index+2*i+3)
                  y_tab(i)=prop(index+2*i+4)
               enddo
            endif
!     
            call cd_bleedtapping(p2,p1,ps1pt1,number,curve,x_tab,y_tab,
     &           cd)
!     
         elseif(lakon(nelem)(2:5).eq.'ORPN') then
!     
!     calculate the discharge coefficient of preswirl nozzle (OWN tables)
!     
            d=dsqrt(4*A/pi)
            reynolds=dabs(xflow)*d/(dvi*a)
            curve=nint(prop(index+4))
            number=nint(prop(index+6))
!     
            if(number.ne.0.d0)then             
               do i=1,number
                  x_tab2(i)=prop(index+2*i+5)
                  y_tab2(i)=prop(index+2*i+6)
               enddo
            endif
!     
            call cd_preswirlnozzle(p2,p1,number,curve,x_tab2,y_tab2,cd)
!     
            theta=prop(index+2)
            k_phi=prop(index+3)
!     
            if(p2/p1.gt.(2/(kappa+1.d0))**(kappa/(kappa-1.d0))) then
               c2u_new=k_phi*cd*sin(theta*Pi/180.d0)*r*
     &              dsqrt(2.d0*kappa/(r*(kappa-1)))*
     &              dsqrt(T1*(1.d0-(p2/p1)**((kappa-1)/kappa)))
!     
            else
               c2u_new=k_phi*cd*sin(theta*Pi/180.d0)*r*
     &              dsqrt(2.d0*kappa/(r*(kappa-1)))*
     &              dsqrt(T1*(1.d0-2/(kappa+1)))
            endif
            prop(index+5)=c2u_new
         endif
!     
         if(cd.gt.1.d0) then
            Write(*,*) '*WARNING:'
            Write(*,*) 'In RESTRICTOR',nelem
            write(*,*) 'Cd greater than 1'
            write(*,*) 'Calculation will proceed using Cd=1'
            cd=1.d0
         endif
         xflow_oil=0
!     
         if(kflag.eq.3)then
!     
          xmach=dsqrt(((p1/p2)**((kappa-1.d0)/kappa)-1.d0)*2.d0/
     &          (kappa-1.d0))
          write(1,*) ''
          write(1,55) ' from node ',node1,
     &   ' to node ', node2,' :   air massflow rate = '
     &,inv*xflow,' ',
     &   ', oil massflow rate = ',xflow_oil,' '
!          
          if(inv.eq.1) then
             write(1,56)'       Inlet node ',node1,' :   Tt1 = ',T1,
     &           '  , Ts1 = ',T1,'  , Pt1 = ',p1, ' '
!             
             write(1,*)'             Element ',nelem,lakon(nelem)
             write(1,57)'             dyn.visc = ',dvi,'  , Re = '
     &           ,reynolds,', M = ',xmach
             if(lakon(nelem)(2:5).eq.'ORC1') then
                write(1,58)'             CD = ',cd
             else if((lakon(nelem)(2:5).eq.'ORMA').or.
     &          (lakon(nelem)(2:5).eq.'ORMM').or.
     &          (lakon(nelem)(2:5).eq.'ORPM').or.
     &          (lakon(nelem)(2:5).eq.'ORPA'))then
                write(1,59)'             CD = ',cd,' , C1u = ',u,
     &  '  , C2u = ', prop(index+7), ' '
             endif
!     special for bleed tappings
             if(lakon(nelem)(2:5).eq.'ORBT') then
                write(1,63) '             P2/P1 = ',p2/p1,
     &' , ps1pt1 = ', ps1pt1, ' , DAB = ',(1-p2/p1)/(1-ps1pt1),
     &' , curve N° = ', curve,' , cd = ',cd
!     special for preswirlnozzles
             elseif(lakon(nelem)(2:5).eq.'ORPN') then
                write(1,62) '             cd = ', cd,
     &'  , C2u = ',c2u_new,' '
!     special for recievers
             endif 
!
             write(1,56)'      Outlet node ',node2,' :   Tt2 = ',T2,
     &           '  , Ts2 = ',T2,'  , Pt2 = ',p2,' '
!     
          else if(inv.eq.-1) then
             write(1,56)'       Inlet node ',node2,':    Tt1 = ',T1,
     &           '  , Ts1 = ',T1,' , Pt1 = ',p1, ' '
     &          
             write(1,*)'             element R    ',set(numf)(1:30)
             write(1,57)'             dyn.visc. = ',dvi,' , Re ='
     &           ,reynolds,', M = ',xmach
             if(lakon(nelem)(2:5).eq.'ORC1') then
                write(1,58)'             CD = ',cd
             else if((lakon(nelem)(2:5).eq.'ORMA').or.
     &          (lakon(nelem)(2:5).eq.'ORMM').or.
     &          (lakon(nelem)(2:5).eq.'ORPM').or.
     &          (lakon(nelem)(2:5).eq.'ORPA'))then
                write(1,59)'             CD = ',cd,' , C1u = ',u,
     &  '  , C2u = ', prop(index+7), ' '
             endif
!     special for bleed tappings
             if(lakon(nelem)(2:5).eq.'ORBT') then
                write(1,63) '             P2/P1 = ',p2/p1,
     &' , ps1pt1 = ', ps1pt1, ' , DAB = ',(1-p2/p1)/(1-ps1pt1),
     &' , curve N° = ', curve,' , cd = ',cd
!     special for preswirlnozzles
             elseif(lakon(nelem)(2:5).eq.'ORPN') then
                write(1,*) ' cd = ', cd,' , C2u = '
     &,c2u_new,' '
             endif
! 
             write(1,56)'      Outlet node ',node1,':    Tt2 = ',T2,
     &           '  , Ts2 = ',T2,'  , Pt2 = ',p2, ' '
         endif
       endif
!
!
      endif
!     
 55   format(1x,a,i6,a,i6,a,e11.4,a,a,e11.4,a)
 56   format(1x,a,i6,a,e11.4,a,e11.4,a,e11.4,a)
 57   format(1x,a,e11.4,a,e11.4,a,e11.4)
 58   format(1x,a,e11.4)
 59   format(1x,a,e11.4,a,e11.4,a,e11.4,a)
 62   format(1x,a,e11.4,a,e11.4,a,e11.4,a)
 63   format(1x,a,e11.4,a,e11.4,a,e11.4,a,i2,a,e11.4)
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!     
      return
      end
