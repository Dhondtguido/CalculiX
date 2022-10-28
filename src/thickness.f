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
      subroutine thickness(nodedesiboun,ndesiboun,objectset,
     &     xo,yo,zo,x,y,z,nx,ny,nz,co,ifree,ndesia,ndesib,
     &     iobject,dgdxglob,nk,extnor,g0,coini)                       
!!     
!     calcualtion of the actual wall thickness      
!     
      implicit none
!     
      character*81 objectset(5,*)

      integer nodedesiboun(*),iactnode,iobject,ndesiboun,neighbor(10),
     &   j,nx(*),ny(*),nz(*),nk,istat,ndesia,ndesib,ifree,nnodes,
     &   irefnode
!     
      real*8 xo(*),yo(*),zo(*),x(*),y(*),z(*),scalar,bound,co(3,*),
     &   xdesi,ydesi,zdesi,actdist,dgdxglob(2,nk,*),coini(3,*),    
     &   inivector(3),g0(*),deltax,deltay,deltaz,extnor(3,*)
!
!     
      read(objectset(1,iobject)(61:80),'(f20.0)',iostat=istat) bound
!         
!     MAXMEMBERSIZE is related to a LE constraint:
!     actdist<=bound --> actdist-bound<=0
!     MINMEMBERSIZE is related to a GE constraint:
!     actdist>=bound --> bound-actdist<=0
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
!       Calculate the vector between the design variable 
!       and the closest node from the reference node set
!       
        irefnode=neighbor(1)
        deltax=xo(irefnode)-xdesi
        deltay=yo(irefnode)-ydesi
        deltaz=zo(irefnode)-zdesi
!
!       check if node lies in negativ normal direction w.r.t. the
!       design variable
!       --> simplified check that node lies on other side of the volume
!
        scalar=deltax*extnor(1,iactnode)+deltay*extnor(2,iactnode)+
     &         deltaz*extnor(3,iactnode)   
!
        if(scalar.lt.0.d0) then
!
!          calculated distance
!  
           actdist=dsqrt(deltax**2+deltay**2+deltaz**2)
           dgdxglob(1,iactnode,iobject)=actdist
!    
!          calculate the function value of the objective function
!
           if(objectset(1,iobject)(1:13).eq.'MINMEMBERSIZE') then
              dgdxglob(2,iactnode,iobject)=bound-actdist
           else if(objectset(1,iobject)(1:13).eq.'MAXMEMBERSIZE') then
              dgdxglob(2,iactnode,iobject)=actdist-bound
           endif
!
!          Calculate the normalized objective function and check if vectors
!          of actual and inital design variable position still point in the
!          same direction --> measure that both points lie on the same side
!          of the bound verifiying that the bound has not been crossed 
!          (only relevant for minmembersize)
!     
           if(objectset(1,iobject)(1:13).eq.'MINMEMBERSIZE') then
!
              inivector(1)=xo(irefnode)-coini(1,iactnode)
              inivector(2)=yo(irefnode)-coini(2,iactnode)
              inivector(3)=zo(irefnode)-coini(3,iactnode)
              scalar=deltax*inivector(1)+deltay*inivector(2)+
     &               deltaz*inivector(3)
              if((scalar.lt.0.d0).and.
     &           (dgdxglob(2,iactnode,iobject).gt.0.d0)) then
                 dgdxglob(2,iactnode,iobject)=
     &              -dgdxglob(2,iactnode,iobject)
              endif
           endif
        else
           write(*,*) '*WARNING no reference node found in negative'
           write(*,*) '         normal direction for node ',iactnode
           write(*,*) '         node ',iactnode,'ignored for'
           write(*,*) '         MEMBERSIZE constraint'
           dgdxglob(1,iactnode,iobject)=-1.0
        endif      
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
