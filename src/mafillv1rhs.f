!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine mafillv1rhs(co,nk,kon,ipkon,lakon,ne,
     &     ipompc,nodempc,coefmpc,nmpc,
     &     nelemload,sideload,xload,nload,xbody,ipobody,
     &     nbody,b2,nactdoh,nmethod,ikmpc,
     &     ilmpc,rhcon,nrhcon,ielmat,ntmat_,ithermal,
     &     vold,vcon,mi,physcon,shcon,nshcon,ttime,
     &     time,istep,ibody,xloadold,iturbulent,
     &     nelemface,sideface,nface,compressible,nea,neb,dtimef,ipvar,
     &     var,ipvarf,varf,ipface,ifreesurface,depth,dgravity,cocon,
     &     ncocon,iinc,theta1,reltimef,b1)
!     
!     filling the rhs b1 for equations of steps 1,2,4 and 5
!     filling the rhs b2 for equations of step 3
!     corrections to step 2 (all fluids) and 3 (only incompressible
!     fluids) are done in mafillprhs and mafillv2rhs
!     
      implicit none
!     
      integer iturbulent,compressible
!     
      character*1 sideface(*)
      character*8 lakon(*)
      character*20 sideload(*)
!     
      integer kon(*),ipompc(*),nodempc(3,*),nelemload(2,*),
     &     ikmpc(*),ilmpc(*),nactdoh(*),
     &     nrhcon(*),mi(*),ielmat(mi(3),*),ipkon(*),nshcon(*),
     &     ipobody(2,*),ipface(*),ifreesurface,ncocon(2,*),
     &     nbody,ibody(3,*),nelemface(*),nface,nea,neb,kstart,
     &     nk,ne,nmpc,nload,nmethod,ithermal(*),i,j,k,id,
     &     ist,index,jdof1,node,ntmat_,indexe,nope,istep,
     &     ipvar(*),ipvarf(*),iinc
!     
      real*8 co(3,*),coefmpc(*),xload(2,*),p1(3),
     &     p2(3),bodyf(3),xloadold(2,*),rhcon(0:1,ntmat_,*),
     &     vold(0:mi(2),*),vcon(nk,0:mi(2)),ff(0:mi(2),8),
     &     physcon(*),shcon(0:3,ntmat_,*),xbody(7,*),var(*),varf(*),
     &     om,dtimef,ttime,time,dgravity,b2(nk,3),
     &     cocon(0:6,ntmat_,*),theta1,bb(3,8),reltimef,b1(nk,0:mi(2)),
     &     depth(*)
!
!     check whether energy equation is needed
!
      if(ithermal(1).gt.1) then
        kstart=0
      else
        kstart=1
      endif
!     
      do i=nea,neb
!     
        indexe=ipkon(i)
        if(lakon(i)(4:4).eq.'8') then
          nope=8
        elseif(lakon(i)(4:4).eq.'4') then
          nope=4
        elseif(lakon(i)(4:4).eq.'6') then
          nope=6
        else
          cycle
        endif
!     
        om=0.d0
!     
        if(nbody.gt.0) then
!     
!     assigning centrifugal forces
!     
          bodyf(1)=0.
          bodyf(2)=0.
          bodyf(3)=0.
!     
          index=i
          do
            j=ipobody(1,index)
            if(j.eq.0) exit
            if(ibody(1,j).eq.1) then
              om=xbody(1,j)
              p1(1)=xbody(2,j)
              p1(2)=xbody(3,j)
              p1(3)=xbody(4,j)
              p2(1)=xbody(5,j)
              p2(2)=xbody(6,j)
              p2(3)=xbody(7,j)
!     
!     assigning gravity forces
!     
            elseif(ibody(1,j).eq.2) then
              bodyf(1)=bodyf(1)+xbody(1,j)*xbody(2,j)
              bodyf(2)=bodyf(2)+xbody(1,j)*xbody(3,j)
              bodyf(3)=bodyf(3)+xbody(1,j)*xbody(4,j)
            endif
            index=ipobody(2,index)
            if(index.eq.0) exit
          enddo
        endif
!     
        call e_c3d_v1rhs(co,nk,kon(indexe+1),lakon(i),p1,p2,om,
     &       bodyf,nbody,ff,i,nmethod,rhcon,nrhcon,ielmat,ntmat_,vold,
     &       vcon,dtimef,mi(1),
     &       ttime,time,istep,shcon,nshcon,
     &       iturbulent,nelemface,sideface,nface,compressible,
     &       ipvar,var,ipvarf,varf,ithermal,ipface,nelemload,
     &       sideload,xload,nload,ifreesurface,depth,dgravity,cocon,
     &       ncocon,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iinc,
     &       theta1,bb,physcon,reltimef,xloadold)
!     
!     distributed forces
!     
        if(compressible.eq.1) then
!
!         compressible fluid: dof = node number
!
          do j=1,nope
            node=kon(indexe+j)
            do k=0,mi(2)
              b1(node,k)=b1(node,k)+ff(k,j)
            enddo
            do k=1,3
              b2(node,k)=b2(node,k)+bb(k,j)
            enddo
          enddo
        else
!
!     incompressible fluid: dof = node number except for the
!     pressure equation, for which the dofs with SPC/MPC are
!     removed
!
          do j=1,nope
            node=kon(indexe+j)
            do k=kstart,3
              b1(node,k)=b1(node,k)+ff(k,j)
            enddo
            do k=1,3
              b2(node,k)=b2(node,k)+bb(k,j)
            enddo
          enddo
!     
!     pressure equation
!     
          do j=1,nope
!     
            node=kon(indexe+j)
            jdof1=nactdoh(node)
            if(jdof1.le.0) then
              if(nmpc.ne.0) then
                if(jdof1.ne.2*(jdof1/2)) then
                  id=(-jdof1+1)/2
                  ist=ipompc(id)
                  index=nodempc(3,ist)
                  if(index.eq.0) cycle
                  do
                    jdof1=nactdoh(nodempc(1,index))
                    if(jdof1.gt.0) then
                      b1(jdof1,4)=b1(jdof1,4)
     &                     -coefmpc(index)*ff(4,j)
     &                     /coefmpc(ist)
                    endif
                    index=nodempc(3,index)
                    if(index.eq.0) exit
                  enddo
                endif
              endif
              cycle
            endif
            b1(jdof1,4)=b1(jdof1,4)+ff(4,j)
          enddo
!
!         turbulent equations
!
          if(iturbulent.gt.0) then
            do j=1,nope
              node=kon(indexe+j)
              do k=5,mi(2)
                b1(node,k)=b1(node,k)+ff(k,j)
              enddo
            enddo
          endif
        endif
      enddo
c      write(*,*) 'mafilltrhs '
c      do i=1,nk
c        write(*,*) i,b1(i,0)
c      enddo
!     
      return
      end
