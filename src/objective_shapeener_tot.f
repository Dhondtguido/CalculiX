!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine objective_shapeener_tot(ne,kon,ipkon,lakon,
     &   fint,vold,iperturb,mi,nactdof,dgdx,df,ndesi,iobject,
     &   jqs,irows,vec,nod1st)
!
      implicit none
!
      character*8 lakon(*)
!
      integer ndesi,iobject,idesvar,i,j,l,jqs(*),irows(*),idof,
     &   ne,ipkon(*),ielem,iperturb(*),indexe,konl(26),kon(*),mi(*),
     &   nope,nactdof(0:mi(2),*),node,nod1st(*)
!      
      real*8 dgdx(ndesi,*),df(*),vec(*),vold(0:mi(2),*),fint(*)
!
!     ----------------------------------------------------------------
!     Calculation of the total differential:
!     non-linear:  dgdx = dgdx + fint^(T) * ( df )
!     linear:      dgdx = dgdx + vold^(T) * ( df )
!     ----------------------------------------------------------------
!
!     copying the entries of vold (linear) or fint (nonlinear) of the nodes
!     belonging to the active element set in the field vec
!
      do ielem=1,ne
!         
         if(ipkon(ielem).lt.0) cycle
!   
         indexe=ipkon(ielem)
!   
         if(lakon(ielem)(4:4).eq.'8') then
            nope=8
         elseif(lakon(ielem)(4:5).eq.'20') then
            nope=20
         elseif(lakon(ielem)(4:5).eq.'10') then
            nope=10
         elseif(lakon(ielem)(4:4).eq.'4') then
            nope=4
         elseif(lakon(ielem)(4:4).eq.'6') then          
            nope=6
         elseif(lakon(ielem)(4:5).eq.'15') then
            nope=15
         else
            exit
         endif
!   
         do l=1,nope
            konl(l)=kon(indexe+l)
         enddo
!
!        field nod1st points for each expanded 3d-node
!        to the node in the mid surface; this is the only node
!        with degrees of freedom (for plane stress/strain/axi)         
!         
         if(iperturb(2).eq.1) then
            do i=1,nope
c              if((lakon(ielem)(7:7).eq.'A').or.
c     &             (lakon(ielem)(7:7).eq.'S').or.       
c     &             (lakon(ielem)(7:7).eq.'E')) then
c                node=nod1st(konl(i))+1     
c              else
                node=konl(i)
c              endif
               do j=1,3
                  idof=nactdof(j,node)
                  if(idof.gt.0) then
                     vec(idof)=fint(idof)
                  endif               
               enddo
            enddo
         else
            do i=1,nope
c              if((lakon(ielem)(7:7).eq.'A').or.
c     &             (lakon(ielem)(7:7).eq.'S').or.       
c     &             (lakon(ielem)(7:7).eq.'E')) then
c                node=nod1st(konl(i))+1
c              else
                node=konl(i)
c              endif
               do j=1,3
                  idof=nactdof(j,node)
                  if(idof.gt.0) then      
                     vec(idof)=vold(j,node)
                  endif              
               enddo
            enddo
         endif
      enddo
!
!     Calculation of the total differential:    
!
      do idesvar=1,ndesi
         do j=jqs(idesvar),jqs(idesvar+1)-1
            idof=irows(j)
            dgdx(idesvar,iobject)=dgdx(idesvar,iobject) 
     &            +vec(idof)*df(j) 
         enddo
      enddo     
!      
      return
      end
