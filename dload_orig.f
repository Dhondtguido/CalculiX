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
      subroutine dload(f,kstep,kinc,time,noel,npt,layer,kspt,
     &     coords,jltyp,loadtype,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,veold,
     &     rho,amat,mi)
!
!     user subroutine dload
!
!
!     INPUT:
!
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     layer              currently not used
!     kspt               currently not used
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        21 = face 1 
!                        22 = face 2 
!                        23 = face 3 
!                        24 = face 4 
!                        25 = face 5 
!                        26 = face 6
!     loadtype           load type label
!     vold(0..4,1..nk)   solution field in all nodes
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!     veold(0..3,1..nk)  derivative of the solution field w.r.t.
!                        time in all nodes
!                        0: temperature rate
!                        1: velocity in global x-direction
!                        2: velocity in global y-direction
!                        3: velocity in global z-direction
!     co(3,1..nk)        coordinates of all nodes
!                        1: coordinate in global x-direction
!                        2: coordinate in global y-direction
!                        3: coordinate in global z-direction
!     lakonl             element label
!     konl(1..20)        nodes belonging to the element
!     ipompc(1..nmpc))   ipompc(i) points to the first term of
!                        MPC i in field nodempc
!     nodempc(1,*)       node number of a MPC term
!     nodempc(2,*)       coordinate direction of a MPC term
!     nodempc(3,*)       if not 0: points towards the next term
!                                  of the MPC in field nodempc
!                        if 0: MPC definition is finished
!     coefmpc(*)         coefficient of a MPC term
!     nmpc               number of MPC's
!     ikmpc(1..nmpc)     ordered global degrees of freedom of the MPC's
!                        the global degree of freedom is
!                        8*(node-1)+direction of the dependent term of
!                        the MPC (direction = 0: temperature;
!                        1-3: displacements; 4: static pressure;
!                        5-7: rotations)
!     ilmpc(1..nmpc)     ilmpc(i) is the MPC number corresponding
!                        to the reference number in ikmpc(i)   
!     rho                local density 
!     amat               material name
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!
!     OUTPUT:
!
!     f                  magnitude of the distributed load
!     iscale             determines whether the flux has to be
!                        scaled for increments smaller than the 
!                        step time in static calculations
!                        0: no scaling
!                        1: scaling (default)
!           
      implicit none
!
      character*8 lakonl
      character*20 loadtype
      character*80 amat
!
      integer kstep,kinc,noel,npt,jltyp,layer,kspt,konl(20),iscale,
     &  mi(*)
!
      real*8 f,time(2),coords(3),vold(0:mi(2),*),co(3,*),rho
!
!
!
!     the code starting here up to the end of the file serves as
!     an example for combined mechanical-lubrication problems. 
!     Please replace it by your own code for your concrete application.
!
      integer ifaceq(8,6),ifacet(6,4),ifacew(8,5),ig,nelem,nopes,
     &  iflag,i,j,nope,ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),
     &  node,idof,id
!
      real*8 xl2(3,8),pres(8),xi,et,xsj2(3),xs2(3,7),shp2(7,8),
     &  coefmpc(*),veold(0:mi(2),*)
!
      include "gauss.f"
!
      ifaceq=reshape((/4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/),(/8,6/))
      ifacet=reshape((/1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/),(/6,4/))
      ifacew=reshape((/1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/),(/8,5/))
      iflag=2
!
      nelem=noel
      ig=jltyp-20
!
      if(lakonl(4:4).eq.'2') then
         nope=20
         nopes=8
      elseif(lakonl(4:4).eq.'8') then
         nope=8
         nopes=4
      elseif(lakonl(4:5).eq.'10') then
         nope=10
         nopes=6
      elseif(lakonl(4:4).eq.'4') then
         nope=4
         nopes=3
      elseif(lakonl(4:5).eq.'15') then
         nope=15
      elseif(lakonl(4:4).eq.'6') then
         nope=6
      endif
!     
!     treatment of wedge faces
!     
      if(lakonl(4:4).eq.'6') then
         if(ig.le.2) then
            nopes=3
         else
            nopes=4
         endif
      endif
      if(lakonl(4:5).eq.'15') then
         if(ig.le.2) then
            nopes=6
         else
            nopes=8
         endif
      endif
!     
      do i=1,nopes
         do j=1,3
            xl2(j,i)=0.d0
         enddo
      enddo
!
      if((nope.eq.20).or.(nope.eq.8)) then
         do i=1,nopes
            node=konl(ifaceq(i,ig))
            idof=8*(node-1)
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in dload: node ',node
               write(*,*) '       is not connected to the oil film'
               call exit(201)
            endif
            node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
            pres(i)=vold(0,node)
         enddo
      elseif((nope.eq.10).or.(nope.eq.4)) then
         do i=1,nopes
            node=konl(ifacet(i,ig))
            node=konl(ifaceq(i,ig))
            idof=8*(node-1)
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in dload: node ',node
               write(*,*) '       is not connected to the oil film'
               call exit(201)
            endif
            node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
            pres(i)=vold(0,node)
         enddo
      else
         do i=1,nopes
            node=konl(ifacew(i,ig))
            node=konl(ifaceq(i,ig))
            idof=8*(node-1)
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in dload: node ',node
               write(*,*) '       is not connected to the oil film'
               call exit(201)
            endif
            node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
            pres(i)=vold(0,node)
         enddo
      endif
!
      i=npt
!     
      if((lakonl(4:5).eq.'8R').or.
     &     ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
         xi=gauss2d1(1,i)
         et=gauss2d1(2,i)
      elseif((lakonl(4:4).eq.'8').or.
     &        (lakonl(4:6).eq.'20R').or.
     &        ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
         xi=gauss2d2(1,i)
         et=gauss2d2(2,i)
      elseif(lakonl(4:4).eq.'2') then
         xi=gauss2d3(1,i)
         et=gauss2d3(2,i)
      elseif((lakonl(4:5).eq.'10').or.
     &        ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
         xi=gauss2d5(1,i)
         et=gauss2d5(2,i)
      elseif((lakonl(4:4).eq.'4').or.
     &        ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
         xi=gauss2d4(1,i)
         et=gauss2d4(2,i)
      endif
!     
      if(nopes.eq.8) then
         call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.4) then
         call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.6) then
         call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
      else
         call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
      endif
!
!     determining the pressure
!
      f=0.d0
      do j=1,nopes
         f=f+pres(j)*shp2(4,j)
      enddo
!
      iscale=0
!
      return
      end

