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
      subroutine dashdamp(xl,elas,konl,voldl,s,imat,elcon,nelcon,
     &  ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,plicon,
     &  nplicon,npmat_,iperturb,time,nmethod)
!
!     calculates the damping coefficient of a dashpot
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),i,j,imat,ncmat_,ntmat_,nope,iperturb(*),niso,
     &  kode,npmat_,nelcon(2,*),nplicon(0:ntmat_,*),nmethod,id
!
      real*8 xl(3,9),elas(21),s(60,60),voldl(3,9),xn(3),dd,
     &  elcon(0:ncmat_,ntmat_,*),t0l,t1l,elconloc(*),damp,
     &  plicon(0:2*npmat_,ntmat_,*),plconloc(802),pl(3,9),time,
     &  xiso(200),yiso(200)
!
!     original positions of the nodes belonging to the dashpot
!
      if(iperturb(1).eq.0) then
         do i=1,nope
            do j=1,3
               pl(j,i)=xl(j,i)
            enddo
         enddo
      else
         do i=1,nope
            do j=1,3
               pl(j,i)=xl(j,i)+voldl(j,i)
            enddo
         enddo
      endif
!     
      dd=dsqrt((pl(1,2)-pl(1,1))**2
     &     +(pl(2,2)-pl(2,1))**2
     &     +(pl(3,2)-pl(3,1))**2)
      do i=1,3
         xn(i)=(pl(i,2)-pl(i,1))/dd
      enddo
!
!     interpolating the material data
!
      call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!
!     calculating the damping force and damping coefficient
!
      if(kode.eq.2) then
         damp=elconloc(1)
      else
         if(nmethod.ne.5) then
            write(*,*) '*ERROR in dashdamp: the damping coefficient'
            write(*,*) '       may depend on temperature and frequency'
            write(*,*) '       only; the latter is only allowed for'
            write(*,*) '       steady state dynamics calculations'
            call exit(201)
         endif
         niso=int(plconloc(801))
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
         call ident(xiso,time,niso,id)
         if(id.eq.0) then
            damp=yiso(1)
         elseif(id.eq.niso) then
            damp=yiso(niso)
         else
            damp=yiso(id)+(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))*
     &          (time-xiso(id))
         endif
      endif
c      write(*,*) 'dashdamp ',time,damp
!     
      do i=1,3
         do j=1,3
            s(i,j)=damp*xn(i)*xn(j)
         enddo
      enddo
      do i=1,3
         do j=1,3
            s(i+3,j)=-s(i,j)
            s(i,j+3)=-s(i,j)
            s(i+3,j+3)=s(i,j)
         enddo
      enddo
!     
      return
      end

