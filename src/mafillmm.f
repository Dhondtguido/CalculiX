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
      subroutine mafillmm(ndesi,nodedesi,co,nodedesiinv,iregion,au,ad,
     &   aub,adb,irow,icol,jq,ipkonfa,konfa,lakonfa,nodedesipos,
     &   ipkonfadesi,nsurfa,nsurfb,area)
!
!     calculation of the entries of the mass matrix and it's derivative
!
      implicit none
!
      character*8 lakonfa(*)
!
      integer idesvar1,idesvar2,kk,l,ipos1,ipos2,n,ndesi,nodedesi(*),
     &   nopes,ifour,ithree,indexel,mint2d,iflag,nodes1(8),indexs,
     &   i,node1,node2,nopedesi,nnodes,nodedesiinv(*),iregion,irow(*),
     &   icol(*),jq(*),nsurfa,nsurfb,isurf,ipkonfa(*),konfa(*),jj,
     &   konfal(8),nodedesipos(*),id,ipointer,ipkonfadesi(*),ipostemp
!
      real*8 xi,et,xl(3,9),xs(3,2),xsj(3),shp(7,9),co(3,*),au(*),ad(*),
     &   aub(*),adb(*),xsjj,weight,val,area(*),dshp1(3),dshp2(3),
     &   shpprj1(3),shpprj2(3),sclprd1,sclprd2
!
!     flag for shape functions
!
      data iflag /3/
      data indexel /0/
      save indexel
      include "gauss.f"
!
!     flag for shape functions
      ifour=4
      ithree=3
!
!---------------------------------------------------------------------
!     Calculating the entries of the matrices
!---------------------------------------------------------------------
!
!     Loop over all external surfaces
      do i=nsurfa,nsurfb
        isurf=ipkonfadesi(i)
        indexs=ipkonfa(isurf)
!     
!       mint2d: # of integration points on the surface
!       nopes:  # of nodes in the surface
!       nope:   # of nodes in the element   
        if(lakonfa(isurf)(1:2).eq.'S3') then
          mint2d=3
          nopes=3
          nopedesi=3
        elseif(lakonfa(isurf)(1:2).eq.'S4') then
          mint2d=4
          nopes=4
          nopedesi=3
        elseif(lakonfa(isurf)(1:2).eq.'S6') then
          mint2d=7
          nopes=6
          nopedesi=4
        elseif(lakonfa(isurf)(1:2).eq.'S8') then
          mint2d=9
          nopes=8
          nopedesi=5
        else
          exit
        endif
        if(iregion.eq.0) nopedesi=0
!     
!       Check if sufficient designnodes on that surface  
        nnodes=0
        do l=1,nopes
          if(nodedesiinv(konfa(indexs+l)).eq.1) then
            nnodes=nnodes+1
          endif
        enddo
!
        if(nnodes.ge.nopedesi) then
!
!         storing node numbers and coordinates of nodes of the surface   
          do l=1,nopes
            konfal(l)=konfa(indexs+l)
            do n=1,3
              xl(n,l)=co(n,konfal(l))
            enddo
          enddo 
!     
!         Evaluate Shape functions and their derivatives 
          do kk=1,mint2d   
            if(lakonfa(isurf)(1:2).eq.'S3') then
              xi=gauss2d5(1,kk)
              et=gauss2d5(2,kk)
              weight=weight2d5(kk)
              call shape3tri(xi,et,xl,xsj,xs,shp,iflag)
            elseif(lakonfa(isurf)(1:2).eq.'S4') then
              xi=gauss2d2(1,kk)
              et=gauss2d2(2,kk)
              weight=weight2d2(kk)
              call shape4q(xi,et,xl,xsj,xs,shp,iflag)
            elseif(lakonfa(isurf)(1:2).eq.'S6') then
              xi=gauss2d6(1,kk)
              et=gauss2d6(2,kk)
              weight=weight2d6(kk)
              call shape6tri(xi,et,xl,xsj,xs,shp,iflag)
            elseif(lakonfa(isurf)(1:2).eq.'S8') then
              xi=gauss2d3(1,kk)
              et=gauss2d3(2,kk)
              weight=weight2d3(kk)
              call shape8q(xi,et,xl,xsj,xs,shp,iflag)       
            endif
!
!           Loop over all nodes on the surface element
            do idesvar1=1,nopes
              node1=konfal(idesvar1)
              ipos1=nodedesipos(node1)
!
!             Loop over all nodes on the surface element
              do idesvar2=1,nopes
                node2=konfal(idesvar2)
                ipos2=nodedesipos(node2)
!
!               check if both nodes are designvariables
                if((nodedesiinv(node1).eq.1).and.
     &             (nodedesiinv(node2).eq.1)) then      
!
!                 Calculate Jacobian determinant     
                  xsjj=dsqrt(xsj(1)**2+xsj(2)**2+xsj(3)**2)
!
!                 entry on the main diagonal
                  if(node1.eq.node2) then
                    adb(ipos1)=adb(ipos1)
     &                        +shp(4,idesvar1)**2*weight*xsjj
!
                    area(ipos1)=area(ipos1)
     &                         +shp(4,idesvar1)*weight*xsjj    
!
                    sclprd1=(shp(1,idesvar1)*xsj(1)
     &                      +shp(2,idesvar1)*xsj(2)
     &                      +shp(3,idesvar1)*xsj(3))/xsjj
                    shpprj1(1)=shp(1,idesvar1)-xsj(1)*sclprd1/xsjj
                    shpprj1(2)=shp(2,idesvar1)-xsj(2)*sclprd1/xsjj
                    shpprj1(3)=shp(3,idesvar1)-xsj(3)*sclprd1/xsjj
!
                    val=(shpprj1(1)**2
     &              +shpprj1(2)**2
     &              +shpprj1(3)**2)*weight*xsjj
                    ad(ipos1)=ad(ipos1)+val

C                   write(5,*) node1,node2,kk
C                   write(5,*) konfal(1),konfal(2),konfal(3)
C                   write(5,*) xl(1,1),xl(2,1),xl(3,1)
C                   write(5,*) xl(1,2),xl(2,2),xl(3,2)
C                   write(5,*) xl(1,3),xl(2,3),xl(3,3)
C                   write(5,*) sclprd1
C                   write(5,*) shp(1,idesvar1),shp(2,idesvar1),shp(3,idesvar1)
C                   write(5,*) shpprj1(1),shpprj1(2),shpprj1(3)
!
!                 entry on sub diagonal
                  elseif(ipos2.gt.ipos1) then
!
                    call nident(irow(jq(ipos1)),ipos2,
     &                          jq(ipos1+1)-jq(ipos1),id)
!
                    ipointer=jq(ipos1)+id-1
!
                    if(irow(ipointer).ne.ipos2) then
                      write(*,*) '*ERROR in mafillmm: 
     &                            coefficient should be 0'
                      call exit(201)
                    else
                      aub(ipointer)=aub(ipointer)
     &                  +shp(4,idesvar1)*shp(4,idesvar2)*weight*xsjj
! 
                      sclprd1=(shp(1,idesvar1)*xsj(1)
     &                        +shp(2,idesvar1)*xsj(2)
     &                        +shp(3,idesvar1)*xsj(3))/xsjj
                      shpprj1(1)=shp(1,idesvar1)-xsj(1)*sclprd1/xsjj
                      shpprj1(2)=shp(2,idesvar1)-xsj(2)*sclprd1/xsjj
                      shpprj1(3)=shp(3,idesvar1)-xsj(3)*sclprd1/xsjj
!           
                      sclprd2=(shp(1,idesvar2)*xsj(1)
     &                        +shp(2,idesvar2)*xsj(2)
     &                        +shp(3,idesvar2)*xsj(3))/xsjj
                      shpprj2(1)=shp(1,idesvar2)-xsj(1)*sclprd2/xsjj
                      shpprj2(2)=shp(2,idesvar2)-xsj(2)*sclprd2/xsjj
                      shpprj2(3)=shp(3,idesvar2)-xsj(3)*sclprd2/xsjj
!
                      val=(shpprj1(1)*shpprj2(1)
     &                    +shpprj1(2)*shpprj2(2)
     &                    +shpprj1(3)*shpprj2(3))*weight*xsjj
                      au(ipointer)=au(ipointer)+val
!
!                     write(5,*) node1,node2,kk
!                     write(5,*) shp(1,idesvar1),shp(2,idesvar1),shp(3,idesvar1)
!                     write(5,*) shp(1,idesvar2),shp(2,idesvar2),shp(3,idesvar2)
!
                    endif            
                  endif
                endif
              enddo
            enddo
          enddo
        endif 
      enddo   
!              
!      do l=1,ndesi
!       write(5,*) nodedesi(l),nodedesi(l),ad(l)
!       do kk=jq(l),jq(l+1)-1
!       write(5,*) nodedesi(l),nodedesi(irow(kk)),au(l)
!       enddo
!      enddo
!
      return
      end
