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
      subroutine carbon_seal(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,R,physcon,dvi,numf,set,mi,ttime,time,
     &     iaxial,iplausi)
!     
!     carbon seal element calculated with Richter method
!      Richter "Rohrhydraulik", Springer ,1971,p. 175
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
     &     inv,mi(*),iaxial,iplausi
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),R,d,dl,
     &     p1,p2,T1,physcon(*),dvi,pi,s,T2,ttime,time
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
!     
         index=ielprop(nelem)
         d=prop(index+1)
         s=prop(index+2)
         dl=prop(index+3)
         pi=4.d0*datan(1.d0)
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
         if(lakon(nelem)(2:6).eq.'CARBS') then
!     
!     gapflow
!     Richter "Rohrhydraulik", Springer ,1971,p. 175
!     
            xflow=inv*Pi*d*s**3*(p1**2-p2**2)/(24.d0*R*T1*dvi*dl)
            
         elseif(lakon(nelem)(2:6).ne.'CARBS') then
            write(*,*) '*WARNING in Carbon_seal.f'
            write(*,*) 'unable to perform carbon seal calculation'
            write(*,*) 'check input file'
         endif
!     
      elseif(kflag.eq.2)then
!     
         numf=4
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
         index=ielprop(nelem)
         d=prop(index+1)
         s=prop(index+2)
         dl=prop(index+3)
         pi=4.d0*datan(1.d0)
         
!     
         if(lakon(nelem)(2:6).eq.'CARBS') then
!     
            f=xflow*T1-pi*d*s**3*(p1**2-p2**2)/(24.d0*R*dvi*dl)
!     
            df(1)=-(pi*d*s**3*p1)/(12.d0*R*dvi*dl)
            df(2)=xflow
            df(3)=T1
            df(4)=(pi*d*s**3*p2)/(12.d0*R*dvi*dl)
!       
         endif

      elseif(kflag.eq.3) then
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

         write(1,*) ''
         write(1,55) ' from node',node1,
     &' to node', node2,':   air massflow rate=',xflow
 55      FORMAT(1X,A,I6,A,I6,A,e11.4,A,A,e11.4,A)

         if(inv.eq.1) then
            write(1,56)'       Inlet node ',node1,':   Tt1=',T1,
     &           ' , Ts1=',T1,' , Pt1=',p1
         
            write(1,*)'             Element',nelem,lakon(nelem)

            write(1,56)'      Outlet node ',node2,':   Tt2=',T2,
     &           ' , Ts2=',T2,' , Pt2=',p2
!     
         else if(inv.eq.-1) then
            write(1,56)'       Inlet node ',node2,':    Tt1=',T1,
     &           ' , Ts1=',T1,' , Pt1=',p1
     &          
            write(1,*)'             Element',nelem,lakon(nelem)

            write(1,56)'      Outlet node ',node1,':    Tt2=',T2,
     &           ' , Ts2=',T2,' , Pt2=',p2

         endif
      
 56      FORMAT(1X,A,I6,A,e11.4,A,e11.4,A,e11.4,A)
      endif
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!
      return
      end
      

