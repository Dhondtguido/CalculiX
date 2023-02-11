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
      subroutine transition(feasdir,nobject,nk,nodedesi,ndesi,
     &           objectset,xo,yo,zo,x,y,z,nx,ny,nz,co,ifree,
     &           ndesia,ndesib)                         
!
!     scaling of sensitivitites between the designspace
!     and non-designspace      
!
      implicit none
!
      character*81 objectset(5,*)

      integer nobject,nk,nodedesi(*),ndesi,j,m,neighbor(10),nx(*),
     &   ny(*),nz(*),istat,ndesia,ndesib,ifree,nnodes,irefnode,        
     &   iactnode
!
      real*8 feasdir(2,*),xo(*),yo(*),zo(*),x(*),y(*),z(*),
     &   trans,co(3,*),xdesi,ydesi,zdesi,scale,actdist,dd,xrefnode,      
     &   yrefnode,zrefnode
!
!     Read transition distance.
!   
      read(objectset(1,1)(21:40),'(f20.0)',iostat=istat) trans
!      
      nnodes=1
!
      do j=ndesia,ndesib
!         
         iactnode=nodedesi(j)
         xdesi=co(1,iactnode)
         ydesi=co(2,iactnode)
         zdesi=co(3,iactnode)
!
         call near3d(xo,yo,zo,x,y,z,nx,ny,nz,xdesi,ydesi,zdesi,
     &        ifree,neighbor,nnodes)
!  
!        Calculate scaled sensitivity 
!
         irefnode=neighbor(1)
         xrefnode=xo(irefnode)
         yrefnode=yo(irefnode)
         zrefnode=zo(irefnode)
         dd=(xrefnode-xdesi)**2+(yrefnode-ydesi)**2
     &      +(zrefnode-zdesi)**2
         actdist=dsqrt(dd)
!
         if(actdist.ge.trans) cycle
!
!        Linear scaling
         scale=actdist/trans
!        Exponential scaling
!        scale=1/(1+dexp(-actdist/trans*10+5))
!
         feasdir(2,iactnode)=feasdir(2,iactnode)*scale
      enddo
!
      return        
      end
