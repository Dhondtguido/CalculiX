!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine springforc_f2f_th(xl,vl,imat,elcon,nelcon,
     &  tnl,ncmat_,ntmat_,nope,lakonl,kode,elconloc,
     &  plicon,nplicon,npmat_,mi,
     &  springarea,nmethod,reltime,jfaces,igauss,
     &  pslavsurf,pmastsurf,clearini,timeend,istep,iinc,
     &  plkcon,nplkcon,node,noel,matname)
!
!     calculates the heat flux across a contact area
!
      implicit none
!
      character*8 lakonl
!
      character*80 matname(*),slname,msname
!
      integer i,j,k,imat,ncmat_,ntmat_,nope,iflag,mi(*),noel,
     &  kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),
     &  nmethod,jfaces,istep,iinc,npred,node,
     &  igauss,nopes,nopem,nopep,nplkcon(0:ntmat_,*)
!
      real*8 xl(3,19),al(3),vl(0:mi(2),19),conductance,
     &  pl(3,19),xn(3),alpha,beta,
     &  elcon(0:ncmat_,ntmat_,*),pproj(3),clear,
     &  xi,et,elconloc(*),plconloc(802),xk,xiso(20),yiso(20),
     &  plicon(0:2*npmat_,ntmat_,*),coords(3),
     &  springarea(2),overlap,clearini(3,9,*),
     &  reltime,weight,xsj2m(3),xs2m(3,7),shp2m(7,9),
     &  xsj2s(3),xs2s(3,7),shp2s(7,9),pslavsurf(3,*),pmastsurf(6,*),
     &  t1ls,t1lm,tmean,pressure,temp(2),timeend(2),ak(5),d(2),tnl(19),
     &  constant,dtemp,flowm(2),predef(2),plkcon(0:2*npmat_,ntmat_,*)
!
      include "gauss.f"
!
      iflag=2
!     
!     # of master nodes
!
      nopem=ichar(lakonl(8:8))-48
!
!     # of slave nodes
!
      nopes=nope-nopem
!
!     actual positions of the master nodes belonging to the contact spring
!     
      do i=1,nopem
         do j=1,3
            pl(j,i)=xl(j,i)+vl(j,i)
         enddo
      enddo
!
!     actual positions of the slave nodes belonging to the contact spring
!     
      do i=nopem+1,nope
         do j=1,3
            pl(j,i)=xl(j,i)+clearini(j,i-nopem,jfaces)*reltime
     &           +vl(j,i)
         enddo
      enddo
!
!     location of integration point in slave face
!
      xi=pslavsurf(1,igauss)
      et=pslavsurf(2,igauss)
      weight=pslavsurf(3,igauss)
!
      iflag=1
c      if(nopes.eq.9) then
c          call shape9q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      if(nopes.eq.8) then
          call shape8q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      elseif(nopes.eq.4) then
          call shape4q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      elseif(nopes.eq.6) then
          call shape6tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
c      elseif(nopes.eq.7) then
c          call shape7tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      else
          call shape3tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      endif
!
      nopep=nope+1
!
      do k=1,3
         pl(k,nopep)=0.d0
      enddo
      t1ls=0.d0
      do j=1,nopes
         do k=1,3
            pl(k,nopep)=pl(k,nopep)+shp2s(4,j)*pl(k,nopem+j)
         enddo
         t1ls=t1ls+shp2s(4,j)*vl(0,nopem+j)
      enddo    
!
!     corresponding location in the master face
!
      xi=pmastsurf(1,igauss)
      et=pmastsurf(2,igauss)
!
!     determining the jacobian vector on the surface 
!
      iflag=2
c      if(nopem.eq.9) then
c         call shape9q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      if(nopem.eq.8) then
         call shape8q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      elseif(nopem.eq.4) then
         call shape4q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      elseif(nopem.eq.6) then
         call shape6tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
c      elseif(nopem.eq.7) then
c         call shape7tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      else
         call shape3tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      endif
!
      t1lm=0.d0
      do i=1,3
         pproj(i)=0.d0
      enddo
      do j=1,nopem
         do i=1,3
            pproj(i)=pproj(i)+shp2m(4,j)*pl(i,j)
         enddo
         t1lm=t1lm+shp2m(4,j)*vl(0,j)
      enddo
!     
!     distance vector between both
!
      do i=1,3
         al(i)=pl(i,nopep)-pproj(i)
      enddo
!
!     normal on the master face
!
      xn(1)=pmastsurf(4,igauss)
      xn(2)=pmastsurf(5,igauss)
      xn(3)=pmastsurf(6,igauss)
!
!     distance from surface along normal (= clearance)
!
      clear=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!
!     check for a reduction of the initial penetration, if any
!
      if(nmethod.eq.1) then
         clear=clear-springarea(2)*(1.d0-reltime)
      endif
!
!     pressure-overclosure relationship
!
      if(int(elcon(3,1,imat)).eq.1) then
!
!        exponential overclosure
!
         if(dabs(elcon(2,1,imat)).lt.1.d-30) then
            pressure=0.d0
            beta=1.d0
         else
!     
            alpha=elcon(2,1,imat)
            beta=elcon(1,1,imat)
            if(-beta*clear.gt.23.d0-dlog(alpha)) then
               beta=(dlog(alpha)-23.d0)/clear
            endif
            pressure=dexp(-beta*clear+dlog(alpha))
         endif
      elseif((int(elcon(3,1,imat)).eq.2).or.
     &       (int(elcon(3,1,imat)).eq.4)) then
!     
!        linear overclosure
!     
         pressure=-elcon(2,1,imat)*clear
      elseif(int(elcon(3,1,imat)).eq.3) then
!     
!        tabular overclosure
!
!        interpolating the material data
!
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,tmean,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
         overlap=-clear
         niso=int(plconloc(81))
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
         call ident(xiso,overlap,niso,id)
         if(id.eq.0) then
            pressure=yiso(1)
         elseif(id.eq.niso) then
            pressure=yiso(niso)
         else
            xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            pressure=yiso(id)+xk*(overlap-xiso(id))
         endif
      endif
!
!     calculating the temperature difference across the contact
!     area and the mean temperature for the determination of the
!     conductance
!
      dtemp=t1lm-t1ls
      tmean=(t1lm+t1ls)/2.d0
!
!     interpolating the material data according to temperature
!
      call materialdata_sp(elcon,nelcon,imat,ntmat_,i,tmean,
     &     elconloc,kode,plkcon,nplkcon,npmat_,plconloc,ncmat_)
!
!     interpolating the material data according to pressure
!
      niso=int(plconloc(801))
!
      if(niso.eq.0) then
         d(1)=clear
         d(2)=pressure
         temp(1)=t1ls
         temp(2)=t1lm
         do k=1,3
            coords(k)=0.d0
            do j=1,nopes
               coords(k)=coords(k)+shp2s(4,j)*xl(k,nopem+j)
            enddo
         enddo
         call gapcon(ak,d,flowm,temp,predef,timeend,matname(imat),
     &               slname,msname,coords,noel,node,npred,istep,iinc,
     &               springarea)
         conductance=ak(1)
      else
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
         call ident(xiso,pressure,niso,id)
         if(id.eq.0) then
            xk=0.d0
            conductance=yiso(1)
         elseif(id.eq.niso) then
            xk=0.d0
            conductance=yiso(niso)
         else
            xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            conductance=yiso(id)+xk*(pressure-xiso(id))
         endif
      endif
!
!     calculating the concentrated heat flow
!
      constant=springarea(1)*conductance*dtemp
!
!     master nodes
!
      do j=1,nopem
         tnl(j)=shp2m(4,j)*constant
      enddo
!
!     slave nodes
!
      do j=1,nopes
         tnl(nopem+j)=-shp2s(4,j)*constant
      enddo
!
      return
      end

