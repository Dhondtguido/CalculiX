!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine findslavcfd(nmpc,labmpc,ipompc,nodempc,islav,
     &  nslav,inoslav,inomast,ics,cs,imast,nmast,co,inomat,
     &  nr,nz,rcs,zcs,rcs0,zcs0,ncs)
!
!     find the slave nodes in a cyclic symmetry constraint for
!     CFD calculations
!
      implicit none
!
      character*20 labmpc(*)
!
      integer i,j,nmpc,nodeslav,nodempc(3,*),ipompc(*),nslav,id,
     &  islav(*),inoslav(*),inomast(*),nodemast,ics(*),imast(*),
     &  nmast,inomat(*),nr(*),nz(*),noden(1),nneigh,ncs,kflag,node
!
      real*8 cs(17,*),co(3,*),xa(3),xn(3),rcs(*),zcs(*),rcs0(*),
     &  zcs0(*),rp,zp,xap,yap,zap,dd
!
      nslav=0
      nmast=0
!
!     cyclic symmetry slave nodes are the dependent nodes in a 
!     CYCLIC MPC; The fluid cyclic MPC's are assumed to be one-to-one
!     node-MPC's (slave and master side must match)
!
!     they are stored and ordered in field islav
!
!     the fluid master nodes are stored in imast
!
!     inoslav contains for each fluid node the corresponding node on the
!             slave side, if any
!     inomast contains for each fluid node the corresponding node on the
!             master side, if any
!
!     xn is an normed vector along the axis
!     xa is a point on the axis
!
      do i=1,3
         xa(i)=cs(5+i,1)
         xn(i)=cs(8+i,1)-xa(i)
      enddo
      dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
      do i=1,3
         xn(i)=xn(i)/dd
      enddo
!
!     catalogue all master nodes (fluid and solid) by their
!     radial and axial coordinate
!
      do i=1,ncs
         node=ics(i)
!
         xap=co(1,node)-xa(1)
         yap=co(2,node)-xa(2)
         zap=co(3,node)-xa(3)
!
         zcs(i)=xap*xn(1)+yap*xn(2)+zap*xn(3)
         rcs(i)=dsqrt((xap-zcs(i)*xn(1))**2+
     &                (yap-zcs(i)*xn(2))**2+
     &                (zap-zcs(i)*xn(3))**2)
      enddo
!
!     initialization for near2d
!
      do i=1,ncs
         nr(i)=i
         nz(i)=i
         rcs0(i)=rcs(i)
         zcs0(i)=zcs(i)
      enddo
      kflag=2
      call dsort(rcs,nr,ncs,kflag)
      call dsort(zcs,nz,ncs,kflag)
!
!     only one neighbor is looked for (one-to-one correspondence
!     of slave and master side assumed)
!
      nneigh=1
!
      do i=1,nmpc
         if(labmpc(i)(1:6).ne.'CYCLIC') cycle
         nodeslav=nodempc(1,ipompc(i))
!
!        check whether slave node is cfd-node
!
         if(inomat(nodeslav).eq.0) cycle
!
         call nident(islav,nodeslav,nslav,id)
         if(id.gt.0) then
            if(islav(id).eq.nodeslav) cycle
         endif
!
         xap=co(1,nodeslav)-xa(1)
         yap=co(2,nodeslav)-xa(2)
         zap=co(3,nodeslav)-xa(3)
!     
         zp=xap*xn(1)+yap*xn(2)+zap*xn(3)
         rp=dsqrt((xap-zp*xn(1))**2+(yap-zp*xn(2))**2+(zap-zp*xn(3))**2)
!
         call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncs,noden,nneigh)
!
         nodemast=ics(noden(1))
!
         if(dabs((rcs0(noden(1))-rp)**2+
     &           (zcs0(noden(1))-zp)**2).gt.1.d-10) then
            write(*,*) '*ERROR in findslavcfd: cyclic symmetric fluid'
            write(*,*) '       slave and master nodes do not match'
            stop
         endif
!
         nslav=nslav+1
         do j=nslav,id+2,-1
            islav(j)=islav(j-1)
         enddo
         islav(id+1)=nodeslav
!
         call nident(imast,nodemast,nmast,id)
         nmast=nmast+1
         do j=nmast,id+2,-1
            imast(j)=imast(j-1)
         enddo
         imast(id+1)=nodemast
!
         inomast(nodeslav)=nodemast
         inoslav(nodemast)=nodeslav
      enddo
!
c      do i=1,nslav
c         write(*,*) i,islav(i)
c      enddo
c      do i=1,nmast
c         write(*,*) i,imast(i)
c      enddo
c      do i=1,38
c         write(*,*) i,inoslav(i),inomast(i)
c      enddo
      return
      end

