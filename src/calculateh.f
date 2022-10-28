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
      subroutine calculateh(nk,v,veold,stn,een,emn,epn,enern,qfn,
     &     errn,h,filab,mi,d,nh,dmin,ipoed,iedg,cotet,jfix)
!     
!     calculating the desired size of h in all nodes of the
!     original mesh
!     
      implicit none
!     
      character*4 label
      character*87 filab(*)
!     
      integer mi(*),i,nk,istat,nh(*),ipoed(*),iedg(3,*),index,n1,n2,
     &     jfix(*)
!     
      real*8 v(0:mi(2),*),veold(0:mi(2),*),stn(6,*),een(6,*),emn(6,*),
     &     epn(*),enern(*),qfn(3,*),errn(6,*),h(*),targetsize,size,d(*),
     &     dmin,cotet(3,*)
!     
!     
!     
      label(1:4)=filab(48)(3:6)
!     
      read(filab(48)(7:26),'(f20.0)',iostat=istat) targetsize
      if(istat.gt.0) then
        write(*,*) '*ERROR in calculateh:'
        write(*,*) '       targetsize not readable'
        call exit(201)
      endif
!     
!     determine the size of all edges
!     
      dmin=1.d30
!     
      loop: do i=1,nk
        index=ipoed(i)
        do
          if(index.eq.0) cycle loop
!     
          n1=iedg(1,index)
          n2=iedg(2,index)
!     
          d(index)=dsqrt((cotet(1,n1)-cotet(1,n2))**2+
     &         (cotet(2,n1)-cotet(2,n2))**2+
     &         (cotet(3,n1)-cotet(3,n2))**2)
!     
          h(n1)=h(n1)+d(index)
          h(n2)=h(n2)+d(index)
          nh(n1)=nh(n1)+1
          nh(n2)=nh(n2)+1
!     
          if(d(index).lt.dmin) dmin=d(index)
!     
          index=iedg(3,index)
        enddo
      enddo loop
!     
      do i=1,nk
!     
!     take the mean at each node
!     
        if(nh(i).le.0) cycle
        h(i)=h(i)/nh(i)
!
!     if the node is on the boundary with part of the mesh which
!     is not being refined: h is set to the mean edge length
!
        if(jfix(i).eq.1) cycle
!     
        if(label.eq.'U   ') then
          size=dsqrt(v(1,i)**2+v(2,i)**2+v(3,i)**2)
        elseif(label.eq.'V   ') then
          size=dsqrt(veold(1,i)**2+veold(2,i)**2+veold(3,i)**2)
        elseif(label.eq.'S   ') then
          size=dsqrt(stn(1,i)**2+stn(2,i)**2+stn(3,i)**2+
     &         2.d0*(stn(4,i)**2+stn(5,i)**2+stn(6,i)**2))
        elseif(label.eq.'E   ') then
          size=dsqrt(een(1,i)**2+een(2,i)**2+een(3,i)**2+
     &         2.d0*(een(4,i)**2+een(5,i)**2+een(6,i)**2))
        elseif(label.eq.'ME  ') then
          size=dsqrt(emn(1,i)**2+emn(2,i)**2+emn(3,i)**2+
     &         2.d0*(emn(4,i)**2+emn(5,i)**2+emn(6,i)**2))
        elseif(label.eq.'PEEQ') then
          size=dabs(epn(i))
        elseif(label.eq.'ENER') then
          size=dabs(enern(i))
        elseif(label.eq.'HFL ') then
          size=dsqrt(qfn(1,i)**2+qfn(2,i)**2+qfn(3,i)**2)
        elseif(label.eq.'ERR ') then
          size=dabs(errn(1,i))
        elseif(label.eq.'USER') then
c     call ucalculateh(v,veold,stn,een,emn,epn,enern,qfn,
c     &           errn,size,mi)
        endif
!     
        if(size/targetsize.gt.1.d0) then
          h(i)=h(i)*targetsize/size
        endif
      enddo
!     
      return
      end
