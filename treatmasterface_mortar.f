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
!     
!     Cuts a convex part of the master surface with a slave surface
!     inserts new active edges into iactiveline for current triangle
!     calling sutherland-hodgeman algorithm
!     
!     [in]       nopes		number of slave nodes for current slave surface
!     [in]       slavstraight	plane equations for mean slave plane	
!     [in]       xn			mean slave normal
!     [in]       xns		slave normals in nodes of slave surface
!     [in]       xl2s		current positions of lsave nodes
!     [in]       xl2sp		projected positions of slave nodes
!     [in]       ipe	 	(i) pointer to ime for node i 
!     [in]       ime     		... cataloging the edges with node i
!     [in,out]   iactiveline 	field storing the active master triange lines for the active line search
!     [in,out]   nactiveline 	number of active lines
!     [in]       nelemm		current master element
!     [in,out]   nintpoint		number of generated integration points
!     [in,out]   pslavsurf		field storing position xil, etal and weight for integration point on slave side
!     [in,out]   imastsurf		pointer into pmastsurf 
!     [in,out]   pmastsurf		field storing position and weights for integration points on master side 
!     [in,out]   xl2m		coordinates of master nodes for current master face
!     [in,out]   nnodelem		number of nodes in curretn master element
!     [in,out]   xl2m2		resorted coordinates of master nodes for current master face
!     [in]       nmp		number of master nodes for current master face
!     [in]       nodem		node numbers for current master face
!     [in,out]   gapmints		gap evaluated at the integration points
!     [in]       issurf		current slave surface index	
!     [in,out]   areaslav		current covering of the slave surface (in the reference element) 
!     [in]       debug		debug output parameter
!     
      subroutine treatmasterface_mortar(
     &     nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &     ipe,ime,iactiveline,nactiveline,
     &     nelemm,nintpoint,pslavsurf,
     &     imastsurf,pmastsurf,xl2m,nnodelem,xl2m2,nmp,
     &     nodem,gapmints,issurf,areaslav,debug,clearance,
     &     cl2s,cl2m,shrink,reltime)
!     
!     Author: Saskia Sitzmann     
!     
      implicit none
!     
      logical debug,shrink
!     
      integer nvertex,nopes,ipe(*),ime(4,*),iactiveline(3,*),
     &     nactiveline,ifreeintersec,nmp,i,j,k,nintpoint,imastsurf(*),
     &     nnodelem,ijk,nodem(*),modf,nelemm,k_max,issurf,ii,jj,nipold
!     
      real*8 pvertex(3,13),slavstraight(36),xn(3),
     &     xilm,etlm,xnl(3),clearance,p(3),dist,
     &     xl2s(3,*),p1(3),p2(3),pslavsurf(3,*),
     &     xil,etl,area,areax,areay,areaz,pmastsurf(2,*),
     &     xl2m(3,8),xl2m2(3,8),al,gapmints(*),err,xns(3,8),
     &     xl2sp(3,*),xl2mp(3,8),cgp(3),spm,spmc,
     &     pm(3),ph(3),reltime,ps(3),xit(3),etat(3),phc(3),
     &     areaslav,cl2s(3,8),cl2m(3,8),psc(3),pmc(3)
!     
!     
!     
      data ijk /0/
      save ijk
!     
      include "gauss.f"
!     
      ifreeintersec=0
      err=1.d-6
      nvertex=0
      nipold=nintpoint
!     
      if(debug) write(20,*) 'TT:slavsurf',issurf,' melem',nelemm
      if(debug) write(20,*) 'TT:xn',(xn(1:3))
!     
!     Project master nodes to meanplane, needed for Sutherland-Hodgman
!     
      do j=1,nmp
        al=-xn(1)*xl2m2(1,j)-xn(2)*
     &       xl2m2(2,j)-xn(3)*
     &       xl2m2(3,j)-slavstraight(nopes*4+4)
        do k=1,3
          xl2mp(k,j)=xl2m2(k,j)+al*xn(k)    
        enddo
      enddo 
!     
      if(debug) then
        write(20,*) 'TT: xl2sp'
        do j=1,nopes
          write(20,*)(xl2sp(k,j),k=1,3)
        enddo       
        write(20,*)'TT:xl2mp'
        do j=1,nmp
          write(20,*)(xl2mp(k,j),k=1,3)
        enddo  
      endif
!     
 111  format(3(e27.20))
!     
!     call Sutherland-Hodgman Algo
!     
      call sutherland_hodgman(nopes,xn,xl2sp,xl2mp,nodem,
     &     ipe,ime,iactiveline,nactiveline,
     &     ifreeintersec,nelemm,nmp,
     &     nvertex,pvertex) 
!     
!     do we have a degenerated triangle?
!     
      if(debug)then
        write(20,*) 'nelemm',nelemm ,'p_new'
        do k=1,nvertex
          write(20,111)(pvertex(i,k),i=1,3)
        enddo
      endif
!     
      do k=1,3
        cgp(k)=0.0
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
            cgp(k)=cgp(k)+pvertex(k,i)/nvertex
          enddo
        enddo
        k_max=nvertex
      endif 
!     
      if(debug)then
        write(20,*) 'TT: nactiveline',nactiveline,'iactiveline'
        write(20,*) (iactiveline(1,k),k=1,nactiveline)
        write(20,*) 'cg' ,(cgp(k),k=1,3)
        write(20,*)'********************************************' 
      endif    
!     
!     Project center point back on slave face
!     
      call attachline(xl2s,cgp,nopes,xit(3),etat(3),xn,p,dist)
!     
!     generating integration points on the slave surface S
!     
      do k=1,k_max
!     
!     storing the triangulation of the slave surfaces
!     
        if(debug) then
          ijk=ijk+1
          write(40,100) ijk,(cgp(i),i=1,3)
          ijk=ijk+1
          write(40,100) ijk,(pvertex(i,modf(nvertex,k)),i=1,3)
          ijk=ijk+1
          write(40,100) ijk,(pvertex(i,modf(nvertex,k+1)),i=1,3)
          write(40,101) ijk-2,ijk-2,ijk-1
          write(40,101) ijk-1,ijk-1,ijk
          write(40,101) ijk,ijk,ijk-2
        endif
 100    format('PNT ',i10,'P',3(1x,e21.14))
 101    format('LINE ',i10,'L',i10,'P ',i10,'P')
!     
!     Project back on slave surface
!     
        call attachline(xl2s,pvertex(1:3,modf(nvertex,k)),
     &       nopes,xit(1),etat(1),xn,p,dist)
        call attachline(xl2s,pvertex(1:3,modf(nvertex,k+1)),
     &       nopes,xit(2),etat(2),xn,p,dist)
        p1(1)=xit(1)-xit(3)
        p1(2)=etat(1)-etat(3)
        p1(3)=0
        p2(1)=xit(2)-xit(3)
        p2(2)=etat(2)-etat(3)
        p2(3)=0
        areax=((p1(2)*p2(3))-(p2(2)*p1(3)))**2
        areay=(-(p1(1)*p2(3))+(p2(1)*p1(3)))**2
        areaz=((p1(1)*p2(2))-(p2(1)*p1(2)))**2
        area=dsqrt(areax+areay+areaz)/2.
        if(area.lt.1.e-4) cycle
        if((nopes.eq.4 .or. nopes.eq.8)
     &       .and.areaslav+area-4.0.gt.1.e-3
     &       .and.nactiveline.gt.0)then
          write(*,*)'TT: face',issurf,'loop in slavintmortar'
          write(*,*)'area',areaslav,'+',area,'.gt.4!',nactiveline
          nactiveline=0
          return
        endif
        if((nopes.eq.3.or.nopes.eq.6)
     &       .and.areaslav+area-0.5.gt.1.e-4
     &       .and.nactiveline.gt.0)then
          write(*,*)'TT: face',issurf,'loop in slavintmortar'
          write(*,*)'area',areaslav,'+',area,'.gt.0.5!'
          write(20,*)'TT: face',issurf,'loop in slavintmortar'
          write(20,*)'area',areaslav,'+',area,'.gt.0.5!'
          nactiveline=0
          return
        endif
        areaslav=areaslav+area
        if(debug)then
          write(20,106) k,area,areaslav
 106      format('tri',i10,' area',e15.8,' atot',e15.8)
          write(20,*) '# fuer itri',k,'werden 7intp gen.'
        endif
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
          call evalshapefunc(xil,etl,xns,nopes,xnl)
          call evalshapefunc(xil,etl,xl2s,nopes,ps)
!     
          nintpoint=nintpoint+1
!     
!     projection of the master integration point onto the
!     master surface in order to get the local coordinates
!     own xn for every integration point?
!     
          call attachline(xl2m,ps,nnodelem,xilm,etlm,xn,p,dist)
          call evalshapefunc(xilm,etlm,xl2m,nnodelem,pm)   
!     
          if(debug)then     
            ijk=ijk+1
            write(40,100) ijk,(ps(jj),jj=1,3)
            ijk=ijk+1
            write(40,100) ijk,(pm(jj),jj=1,3)
            write(40,101) ijk-1,ijk-1,ijk 
          endif
!     
!     Calculation of the gap function at the integration point
!     
          do ii=1,3
            ph(ii)=pm(ii)-ps(ii)
          enddo

          al=sqrt(ph(1)**2+ph(2)**2+ph(3)**2)
          spm=ph(1)*xn(1)
     &         +ph(2)*xn(2)
     &         +ph(3)*xn(3)
!     
          if(.not.((clearance.gt.1.2357111316d0).and.
     &         (clearance.lt.1.2357111318d0))) then
!     
!     assuming zero displacements at the beginning of the calculation
!     only valid for small tangential movement in the contact zone
!     
            call evalshapefunc(xil,etl,cl2s,nopes,psc)
            call evalshapefunc(xilm,etlm,cl2m,nnodelem,pmc)
            do ii=1,3
              phc(ii)=pmc(ii)-psc(ii)
            enddo
            spmc=phc(1)*xn(1)
     &           +phc(2)*xn(2)
     &           +phc(3)*xn(3)
            spm=spm-spmc+clearance
          endif
          if(shrink)then
            spm=reltime*spm
          endif 
          gapmints(nintpoint)=spm
          pslavsurf(1,nintpoint)=xil
          pslavsurf(2,nintpoint)=etl
!     
!     weight add to 0.5 \hat{w}_p=A_triaref*w_p=0.5*w_p
!     division by 0.5 to get to original weights...
!     
          pslavsurf(3,nintpoint)=area*weight2d6(i)/0.5
          pmastsurf(1,nintpoint)=xilm
          pmastsurf(2,nintpoint)=etlm
          imastsurf(nintpoint)=nelemm
          if(debug)then
            write(30,201) xil,etl,pslavsurf(3,nintpoint),spm,nelemm 
          endif 
 201      format('xil ',1x,e15.8,2x,'etl',2x,e15.8,2x,'w ',e15.8,2x,
     &         e15.8,2x,'M',i10)
        enddo
!     
      enddo
c     write(20,*)'********************'
      return
      end
