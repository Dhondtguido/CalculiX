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
      subroutine packaging(nodedesiboun,ndesiboun,objectset,
     &     xo,yo,zo,x,y,z,nx,ny,nz,co,ifree,ndesia,ndesib,
     &     iobject,ndesi,dgdxglob,nk,extnor,g0,nodenum)                       
!     
!     calcualtion of the actual wall thickness      
!     
      implicit none
!     
      character*81 objectset(5,*)

      integer nodedesiboun(*),ndesi,iactnode,iobject,ndesiboun,j,
     &   neighbor(10),nx(*),ny(*),nz(*),nk,ndesia,ndesib,ifree,
     &   nnodes,irefnode,nodenum(*),irefnodeid
!     
      real*8 xo(*),yo(*),zo(*),x(*),y(*),z(*),extnor(3,*),co(3,*),
     &   xdesi,ydesi,zdesi,actdist,dgdxglob(2,nk,*),g0(*),funcval,     
     &   deltax,deltay,deltaz   
!
!     
!     PACKAGING is related to a GE constraint:
!     actdist>=0.0
!     
!     find minimum distance 
!     
      nnodes=1
      do j=ndesia,ndesib
!     
        iactnode=nodedesiboun(j) 
        xdesi=co(1,iactnode)
        ydesi=co(2,iactnode)
        zdesi=co(3,iactnode)
!     
        call near3d(xo,yo,zo,x,y,z,nx,ny,nz,xdesi,ydesi,zdesi,
     &       ifree,neighbor,nnodes)
!     
!       Calculate the distance between the design variable 
!       and the closest node from the reference node set
!       
        irefnode=neighbor(1)
        deltax=xo(irefnode)-xdesi
        deltay=yo(irefnode)-ydesi
        deltaz=zo(irefnode)-zdesi
!
        actdist=dsqrt(deltax**2+deltay**2+deltaz**2)
        dgdxglob(1,iactnode,iobject)=actdist
!    
!       Calculate the function value of the objective function
!
        irefnodeid=nodenum(irefnode)
        funcval=deltax*extnor(1,irefnodeid)+
     &          deltay*extnor(2,irefnodeid)+
     &          deltaz*extnor(3,irefnodeid)
        dgdxglob(2,iactnode,iobject)=funcval
! 
!       count number of active nodes
!
        if(dgdxglob(2,iactnode,iobject).ge.0.d0) then
           g0(iobject)=g0(iobject)+1.d0
        endif
!
      enddo
!
      return        
      end
