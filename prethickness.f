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
      subroutine prethickness(co,xo,yo,zo,x,y,z,nx,ny,nz,ifree,
     &      nodedesiinv,ndesiboun,nodedesiboun,set,nset,objectset,
     &      iobject,istartset,iendset,ialset)       
!
      implicit none
!
      character*81 objectset(5,*),set(*)
!
      integer j,k,i,ifree,nx(*),ny(*),nz(*),kflag,ndesinode,id,
     &  nodedesiinv(*),ndesiboun,nodedesiboun(*),nset,
     &  iobject,istartset(*),iendset(*),ialset(*)
!
      real*8 co(3,*),xo(*),yo(*),zo(*),x(*),y(*),z(*)
!
!     determining the set of boundary nodes
!
c      do i=1,nset
c         if(objectset(4,iobject).eq.set(i)) exit
c      enddo
      call cident81(set,objectset(4,iobject),nset,id)
      i=nset+1
      if(id.gt.0) then
        if(objectset(4,iobject).eq.set(id)) then
          i=id
        endif
      endif
!
      if(i.le.nset) then
!
!     all nodes which define the boundary are put in a set
!    
         ifree=0 
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               k=ialset(j)
               ifree=ifree+1
               xo(ifree)=co(1,k)
               x(ifree)=xo(ifree)
               nx(ifree)=ifree
               yo(ifree)=co(2,k)
               y(ifree)=yo(ifree)      
               ny(ifree)=ifree
               zo(ifree)=co(3,k)
               z(ifree)=zo(ifree)
               nz(ifree)=ifree
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  ifree=ifree+1
                  xo(ifree)=co(1,k)
                  x(ifree)=xo(ifree)
                  nx(ifree)=ifree
                  yo(ifree)=co(2,k)
                  y(ifree)=yo(ifree)      
                  ny(ifree)=ifree
                  zo(ifree)=co(3,k)
                  z(ifree)=zo(ifree)
                  nz(ifree)=ifree
               enddo
            endif      
         enddo
!
!     Sorting of nodes set w.r.t. coordinates
!      
         kflag=2
         call dsort(x,nx,ifree,kflag)
         call dsort(y,ny,ifree,kflag)
         call dsort(z,nz,ifree,kflag)
!
      endif
!
!     determining the set of designvariables which have a wall
!     thickness constraint
!
c      do i=1,nset
c         if(objectset(3,iobject).eq.set(i)) exit
c      enddo
      call cident81(set,objectset(3,iobject),nset,id)
      i=nset+1
      if(id.gt.0) then
        if(objectset(3,iobject).eq.set(id)) then
          i=id
        endif
      endif
!
      if(i.le.nset) then
!
!     all designvariables which have the constraint are put in the set
!     
         ndesiboun=0
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               k=ialset(j)
               if(nodedesiinv(k).ne.1) cycle
               ndesiboun=ndesiboun+1
               nodedesiboun(ndesiboun)=k
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if(nodedesiinv(k).ne.1) cycle
                  ndesiboun=ndesiboun+1
                  nodedesiboun(ndesiboun)=k
               enddo
            endif
         enddo
      endif  
!
      return
      end
