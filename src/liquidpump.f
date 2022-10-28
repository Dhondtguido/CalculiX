!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine liquidpump(node1,node2,nodem,nelem,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,rho,g,co,numf,mi,ttime,time,
     &     iaxial,iplausi)
!
!     pump for incompressible media
!     
      implicit none
!     
      logical identity
!      
      integer nelem,nactdog(0:3,*),node1,node2,nodem,
     &     ielprop(*),nodef(*),idirf(*),index,kflag,
     &     inv,id,numf,npu,i,mi(*),iaxial,iplausi
!      
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),ttime,time,
     &     p1,p2,rho,g(3),dg,z1,z2,co(3,*),
     &     xpu(10),ypu(10),xxpu(10),yypu(10),dh
!
!
!
      numf=3
!
      if (kflag.eq.0) then
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
      elseif((kflag.eq.1).or.(kflag.eq.2)) then
         if(kflag.eq.1) then
            if(v(1,nodem).ne.0.d0) then
               xflow=v(1,nodem)
               return
            endif
         endif
!     
         index=ielprop(nelem)
!
         npu=nint(prop(index+1))
         do i=1,npu
            xpu(i)=prop(index+2*i)
            ypu(i)=prop(index+2*i+1)
         enddo
!     
         p1=v(2,node1)
         p2=v(2,node2)
!     
         z1=-g(1)*co(1,node1)-g(2)*co(2,node1)-g(3)*co(3,node1)
         z2=-g(1)*co(1,node2)-g(2)*co(2,node2)-g(3)*co(3,node2)
!     
         if(kflag.eq.2) then
            xflow=v(1,nodem)*iaxial
            if(xflow.ge.0.d0) then
               inv=1
            else
               inv=-1
            endif
            nodef(1)=node1
            nodef(2)=nodem
            nodef(3)=node2
            idirf(1)=2
            idirf(2)=1
            idirf(3)=2
         endif
!     
         dg=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
!     
         if(kflag.eq.1) then
            dh=(z2-z1+(p2-p1)/rho)/dg
!
!           reverting the order in xpu and ypu and storing the
!           result in xxpu and yypu
!
            do i=1,npu
               xxpu(i)=xpu(npu+1-i)
               yypu(i)=ypu(npu+1-i)
            enddo
            call ident(yypu,dh,npu,id)
            if(id.eq.0) then
               xflow=xxpu(1)
            elseif(id.eq.npu) then
               xflow=0.d0
            else
               xflow=xxpu(id)+(xxpu(id+1)-xxpu(id))*(dh-yypu(id))/
     &               (yypu(id+1)-yypu(id))
            endif
         else
            df(1)=1.d0/rho
            df(3)=-df(1)
            xflow=xflow/rho
            call ident(xpu,xflow,npu,id)
            if(id.eq.0) then
               if(xflow.ge.0.d0) then
                  f=z1-z2+(p1-p2)/rho+dg*ypu(1)
                  df(2)=0.d0
               else
                  df(2)=-1.d10
                  f=z1-z2+(p1-p2)/rho+dg*(ypu(1)+xflow*df(2))
                  df(2)=df(2)*dg/rho
               endif
            elseif(id.eq.npu) then
               df(2)=-1.d10
               f=z1-z2+(p1-p2)/rho+dg*(ypu(npu)+df(2)*(xflow-xpu(npu)))
               df(2)=df(2)*dg/rho
            else
               df(2)=(ypu(id+1)-ypu(id))/(xpu(id+1)-xpu(id))
               f=z1-z2+(p1-p2)/rho+dg*(ypu(id)+(xflow-xpu(id))*df(2))
               df(2)=df(2)*dg/rho
            endif
         endif
!     
      endif
!     
      xflow=xflow/iaxial
      df(2)=df(2)*iaxial
!     
      return
      end
