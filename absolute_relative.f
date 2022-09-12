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
      subroutine absolute_relative(node1,node2,nodem,nelem,lakon,
     &     kon,ipkon, nactdog,identity,ielprop,prop,kflag,v,
     &     xflow,f,nodef,idirf,df,cp,R,physcon,numf,set,mi,ttime,time,
     &     iaxial,iplausi)
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
     &     ipkon(*),kon(*),nelemswirl,mi(*),iaxial,iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),kappa,R,
     &     p1,p2,T1,T2,cp,physcon(*),km1,kp1,kdkm1,
     &     kdkp1,u,pi,Qred_crit,pt1,pt2,Tt1,Tt2,ct,fact,
     &     Cp_cor,ttime,time
!
!
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
         xflow=0.d0
!     
      elseif(kflag.eq.2) then
!     
         numf=4
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!     
         index=ielprop(nelem)
!         
         u=prop(index+1)
         ct=prop(index+2)
         nelemswirl=nint(prop(index+3))
!
        if(nelemswirl.ne.0) then
!
!     previous element is a preswirl nozzle
!
            if(lakon(nelemswirl)(2:5).eq.'ORPN') then
               ct=prop(ielprop(nelemswirl)+5)
!
!     previous element is a forced vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
               ct=prop(ielprop(nelemswirl)+7)
!
!     previous element is a free vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
               ct=prop(ielprop(nelemswirl)+9)
            endif
         endif                
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         if(lakon(nelem)(2:4).eq.'ATR') then
!     
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)-physcon(1)
            Tt2=v(0,node2)-physcon(1)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
!     in the case of a negative flow direction
!     
            if(xflow.le.0d0) then
               write(*,*)''
               write(*,*)'*WARNING:'
               write(*,*)'in element',nelem
               write(*,*)'TYPE=ABSOLUTE TO RELATIVE'
               write(*,*)'mass flow negative!'
               write(*,*)'check results and element definition'
            endif
!       
         elseif(lakon(nelem)(2:4).eq.'RTA') then
!     
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)-physcon(1)
            Tt2=v(0,node2)-physcon(1)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
            if(xflow.le.0d0) then
               write(*,*)''
               write(*,*)'*WARNING:'
               write(*,*)'in element',nelem
               write(*,*)'TYPE=RELATIVE TO ABSOLUTE'
               write(*,*)'mass flow negative!'
               write(*,*)'check results and element definition'
            endif
         endif
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
!     computing temperature corrected Cp=Cp(T) coefficient 
         call cp_corrected(cp,Tt1,Tt2,cp_cor)
!     
            if(Tt1.lt.273) then
               Tt1= Tt2
            endif
!
            if(cp_cor.eq.0.d0) then
               cp_cor=cp
            endif
!     
!     transformation from absolute system to relative system
!     
         if(lakon(nelem)(2:4).eq.'ATR') then
!     
            fact=1+(u**2-2*u*ct)/(2*Cp_cor*Tt1)
!     
            f=pt2-pt1*(fact)**kdkm1
!     
!     pressure node 1
!
            df(1)=-fact**kdkm1
!
!     temperature node1
!
            df(2)=-pt1*Kdkm1*(-(u**2-2*u*ct)/(2*Cp_cor*Tt1**2))
     &           *fact**(kdkm1-1)
!
!     mass flow node m
!
            df(3)=0
!
!     pressure node 2
!     
            df(4)=1
!     
!     transformation from relative system to absolute system
!     
         elseif(lakon(nelem)(2:4).eq.'RTA') then
!     
            fact=1-(u**2-2*u*ct)/(2*Cp*Tt1)
!     
            f=pt2-pt1*(fact)**kdkm1
!     
            df(1)=-fact**kdkm1
!     
            df(2)=-pt1*Kdkm1*((u**2-2*u*ct)/(2*Cp*Tt1**2))
     &           *fact**(kdkm1-1)
!     
            df(3)=0
!     
            df(4)=1
!     
         endif

      elseif(kflag.eq.3) then
            
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!     
         index=ielprop(nelem)
!         
         u=prop(index+1)
         ct=prop(index+2)
         nelemswirl=nint(prop(index+3))
!
        if(nelemswirl.ne.0) then
!
!     previous element is a preswirl nozzle
!
            if(lakon(nelemswirl)(2:5).eq.'ORPN') then
               ct=prop(ielprop(nelemswirl)+5)
!
!     previous element is a forced vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
               ct=prop(ielprop(nelemswirl)+7)
!
!     previous element is a free vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
               ct=prop(ielprop(nelemswirl)+9)
            endif
         endif                
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         if(lakon(nelem)(2:4).eq.'ATR') then
!     
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)-physcon(1)
            Tt2=v(0,node2)-physcon(1)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
!     in the case of a negative flow direction
!     
            if(xflow.le.0d0) then
               write(*,*)''
               write(*,*)'*WARNING:'
               write(*,*)'in element',nelem
               write(*,*)'TYPE=ABSOLUTE TO RELATIVE'
               write(*,*)'mass flow negative!'
               write(*,*)'check results and element definition'
            endif
!       
         elseif(lakon(nelem)(2:4).eq.'RTA') then
!     
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)-physcon(1)
            Tt2=v(0,node2)-physcon(1)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
            if(xflow.le.0d0) then
               write(*,*)''
               write(*,*)'*WARNING:'
               write(*,*)'in element',nelem
               write(*,*)'TYPE=RELATIVE TO ABSOLUTE'
               write(*,*)'mass flow negative!'
               write(*,*)'check results and element definition'
            endif
         endif
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
!     computing temperature corrected Cp=Cp(T) coefficient 
         call cp_corrected(cp,Tt1,Tt2,cp_cor)
!     
         if(Tt1.lt.273) then
            Tt1= Tt2
         endif
!     
         if(cp_cor.eq.0.d0) then
            cp_cor=cp
         endif

                  write(1,*) ''
         write(1,55) ' from node',node1,
     &' to node', node2,':   air massflow rate=',xflow,''
 55      FORMAT(1X,A,I6,A,I6,A,e11.4,A,A,e11.4,A)

            write(1,56)'       Inlet node ',node1,':     Tt1= ',Tt1,
     &           ', Ts1= ',Tt1,', Pt1= ',pt1
            write(1,*)'             Element ',nelem,lakon(nelem)
            write(1,57)'             u= ',u,' ,Ct= ',Ct,''
            write(1,56)'      Outlet node ',node2,':     Tt2= ',Tt2,
     &           ', Ts2= ',Tt2,', Pt2= ',pt2
!     
 56      FORMAT(1X,A,I6,A,e11.4,A,e11.4,A,e11.4,A,e11.4)  
 57      FORMAT(1X,A,e11.4,A,e11.4,A)

      endif
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!     
      return 
      end
