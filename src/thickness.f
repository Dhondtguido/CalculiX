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
      subroutine thickness(nobject,nodedesiboun,ndesiboun,
     &     objectset,xo,yo,zo,x,y,z,nx,ny,nz,co,ifree,ndesia,
     &     ndesib,iobject,ndesi,dgdxglob,nk,normdesi)                       
!     
!     calcualtion of the actual wall thickness      
!     
      implicit none
!     
      character*81 objectset(5,*)

      integer nobject,nodedesiboun(*),ndesi,iactnode,iobject,
     &     ndesiboun,j,neighbor(10),nx(*),ny(*),nz(*),nk,
     &     istat,ndesia,ndesib,ifree,nnodes,irefnode
!     
      real*8 xo(*),yo(*),zo(*),x(*),normdesi(3,*),scalar,
     &     y(*),z(*),refdist,co(3,*),xdesi,ydesi,zdesi,      
     &     actdist,dd,xrefnode,yrefnode,zrefnode,
     &     numericdist,dgdxglob(2,nk,nobject)
!     
      read(objectset(1,iobject)(61:80),'(f20.0)',iostat=istat) refdist
!     
      nnodes=1
      numericdist=0.02d0
!     
!     check if upper or lower constraint is defined
!     
      if(objectset(1,iobject)(19:20).eq.'LE') then
        refdist=refdist*(1-numericdist)
      else
        refdist=refdist*(1+numericdist)
      endif
!     
!     find minimum distance 
!     
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
!     Calculate distance and check if node is active
!     
        irefnode=neighbor(1)
        xrefnode=xo(irefnode)
        yrefnode=yo(irefnode)
        zrefnode=zo(irefnode)
        dd=(xrefnode-xdesi)**2+(yrefnode-ydesi)**2
     &       +(zrefnode-zdesi)**2
        actdist=dsqrt(dd)
        scalar=xrefnode*normdesi(1,j)+
     &         yrefnode*normdesi(2,j)+
     &         zrefnode*normdesi(3,j)
!    
!       distinction whether the check nodeset lies in positive or
!       negative normal direction.
!       positive normal direction: dgdxglob(2,node,iobject)=1
!       negative normal direction: dgdxglob(2,node,iobject)=-1
!       important to know for the correct sign (constraint active or
!       nor active) of the lagrange multiplier
!  
        dgdxglob(1,iactnode,iobject)=actdist
        if(objectset(1,iobject)(19:20).eq.'LE') then
          if(actdist.gt.refdist) then
            if(scalar.gt.0.d0) then  
              dgdxglob(2,iactnode,iobject)=1.0d0
            else
              dgdxglob(2,iactnode,iobject)=-1.0d0
            endif
          endif
        else
          if(actdist.lt.refdist) then
            if(scalar.gt.0.d0) then  
              dgdxglob(2,iactnode,iobject)=1.0d0
            else
             dgdxglob(2,iactnode,iobject)=-1.0d0
            endif
          endif
        endif      
      enddo
!     
      return        
      end
