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
!     cuts a 3- or 4-noded polygon of the master surface with a slave surface
!     inserts new active edges into iactiveline for current triangle
!     
      subroutine treatmasterface_load(nopes,slavstraight,xn,xl2m,xl3sp,
     &     ipe,ime,iactiveline,nactiveline,ifacem,nintpoint,imastload,
     &     pmastload,xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &     xload,nload,nload_,iamload,nam,ifaces)
!     
!     based on treatmasterface.f    
!     
      implicit none
!
      character*20 sideload(*),label
!     
      integer nvertex,nopes,ipe(*),ime(4,*),iactiveline(2,*),nopem,
     &     nactiveline,npg,i,j,k,nintpoint,imastload(2,*),jfacem,ier,
     &     nodepg(*),modf,ifacem,k_max,nelemload(2,*),iamload(2,*),
     &     nload,nload_,nam,nelemm,id,iloadm,ifaces
!     
      real*8 pvertex(3,13),slavstraight(36),xn(3),xl2m(3,*),p1(2),p2(2),
     &     pmastload(3,*),xil,etl,p(3),dist,area,xlpg(3,8),al,err,
     &     xl3sp(3,*),xlpgp(3,8),cgp(3),xit(3),etat(3),areaslav,
     &     xload(2,*)
!     
      include "gauss.f"
!     
      err=1.d-6
      nvertex=0
!     
!     Project master nodes to meanplane, needed for Sutherland-Hodgman
!     
      do j=1, npg
        al=-xn(1)*xlpg(1,j)-xn(2)*
     &       xlpg(2,j)-xn(3)*
     &       xlpg(3,j)-slavstraight(nopes*4+4)
        do k=1,3
          xlpgp(k,j)=xlpg(k,j)+al*xn(k)    
        enddo
      enddo 
!     
!     call Sutherland-Hodgman Algo
!     
      call sutherland_hodgman(nopes,xn,xl3sp,xlpgp,nodepg,
     &     ipe,ime,iactiveline,nactiveline,
     &     ifacem,npg,nvertex,pvertex) 
!     
      do k=1,3
        cgp(k)=0.0d0
      enddo
      if(nvertex.lt.3) return       
!     
      if(nvertex.eq.3)then
        do k=1,3
          cgp(k)=pvertex(k,nvertex)
        enddo
        nvertex=nvertex-1
        k_max=1
      else
        do i=1,nvertex
          do k=1,3
c     cgp(k)=cgp(k)+pvertex(k,i)/nvertex
            cgp(k)=cgp(k)+pvertex(k,i)
          enddo
        enddo
        do k=1,3
          cgp(k)=cgp(k)/nvertex
        enddo
        k_max=nvertex
      endif 
!     
!     Project center point back on slave face
!     
      call attachline(xl2m,cgp,nopem,xit(3),etat(3),xn,p,dist)
!     
!     generating integration points on the slave surface S
!     
      do k=1,k_max
!     
!     Project back on slave surface
!     
        call attachline(xl2m,pvertex(1:3,modf(nvertex,k)),
     &       nopem,xit(1),etat(1),xn,p,dist)
        call attachline(xl2m,pvertex(1:3,modf(nvertex,k+1)),
     &       nopem,xit(2),etat(2),xn,p,dist)
!     
        p1(1)=xit(1)-xit(3)
        p1(2)=etat(1)-etat(3)
!     
        p2(1)=xit(2)-xit(3)
        p2(2)=etat(2)-etat(3)
!     
        area=dabs(p1(1)*p2(2)-p2(1)*p1(2))/2.d0
!     
        if(area.lt.1.d-4) cycle
        if(((nopes.eq.4).or.(nopes.eq.8))
     &       .and.(areaslav+area-4.0d0.gt.1.d-3)
     &       .and.(nactiveline.gt.0))then
          nactiveline=0
          return
        endif
        if(((nopes.eq.3).or.(nopes.eq.6))
     &       .and.(areaslav+area-0.5d0.gt.1.d-4)
     &       .and.(nactiveline.gt.0))then
          nactiveline=0
          return
        endif
        areaslav=areaslav+area
!     
!     7 points scheme
!     
        do i=1,7
          xil=xit(3)*gauss2d6(1,i)+
     &         xit(1)*gauss2d6(2,i)+
     &         xit(2)*(1-gauss2d6(1,i)-gauss2d6(2,i))
          
          etl=etat(3)*gauss2d6(1,i)+
     &         etat(1)*gauss2d6(2,i)+
     &         etat(2)*(1-gauss2d6(1,i)-gauss2d6(2,i))
!     
          nintpoint=nintpoint+1
!     
          pmastload(1,nintpoint)=xil
          pmastload(2,nintpoint)=etl
!     
!     weights sum up to 0.5 for triangles in local coordinates
!     
          pmastload(3,nintpoint)=2.d0*area*weight2d6(i)
!
!         creating load
!
          nelemm=int(ifacem/10.d0)
          jfacem=ifacem-10*nelemm
          label='I                   '
          write(label(2:2),'(i1)') jfacem
          ier=1
          call nident2(nelemload,nelemm,nload,id)
          if(id.gt.0) then
            if(nelemload(1,id).eq.nelemm) then
              do
                if(sideload(id).eq.label) then
                  ier=0
                  exit
                endif  
                id=id-1
                if(id.eq.0) exit
                if(nelemload(1,id).ne.nelemm) exit
              enddo
            endif
          endif
!
!           new load
!
          if(ier.eq.1) then
            nload=nload+1
            if(nload.gt.nload_) then
              write(*,*)
     &             '*ERROR in treatmasterface_load: increase nload_'
              call exit(201)
            endif
!     
!     shifting existing loading
!     
            do j=nload,id+2,-1
              nelemload(1,j)=nelemload(1,j-1)
              nelemload(2,j)=nelemload(2,j-1)
              sideload(j)=sideload(j-1)
              xload(1,j)=xload(1,j-1)
              xload(2,j)=xload(2,j-1)
              if(nam.gt.0) then
                iamload(1,j)=iamload(1,j-1)
                iamload(2,j)=iamload(2,j-1)
              endif
            enddo
!     
!     inserting new loading
!     
            nelemload(1,id+1)=nelemm
            nelemload(2,id+1)=0
            sideload(id+1)=label
            xload(1,id+1)=0.d0
            xload(2,id+1)=0.d0
            if(nam.gt.0) then
              iamload(1,id+1)=0
              iamload(2,id+1)=0
            endif
            iloadm=id+1
          else
            iloadm=id
          endif
!          
          imastload(1,nintpoint)=ifaces
          imastload(2,nintpoint)=nelemload(2,iloadm)
          nelemload(2,iloadm)=nintpoint
!          
        enddo
      enddo
!     
      return
      end
