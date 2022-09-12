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
      subroutine objective_shapeener_tot(ne,kon,ipkon,lakon,
     &   fint,vold,iperturb,mi,nactdof,dgdx,df,ndesi,iobject,
     &   jqs,irows,vec,iponk2dto3d)
!
      implicit none
!
      character*8 lakon(*)
!
      integer ndesi,iobject,idesvar,i,j,l,jqs(*),irows(*),idof,
     &   ne,ipkon(*),ielem,iperturb(*),indexe,konl(26),kon(*),mi(*),
     &   nope,nactdof(0:mi(2),*),node,iponk2dto3d(*)
!      
      real*8 dgdx(ndesi,*),df(*),vec(*),vold(0:mi(2),*),fint(*)
!
!
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
!        field iponk2dto3d points for each expanded 3d-node
!        to the node in the mid surface; this is the only node
!        with degrees of freedom (for plane stress/strain/axi)         
!         
         if(iperturb(2).eq.1) then
            do i=1,nope
               do j=1,3
                  if((lakon(ielem)(7:7).eq.'A').or.
     &               (lakon(ielem)(7:7).eq.'S').or.       
     &               (lakon(ielem)(7:7).eq.'E')) then
                     node=iponk2dto3d(konl(i))     
                  else
                     node=konl(i)
                  endif
                  idof=nactdof(j,node)
                  if(idof.gt.0) then
                     vec(idof)=fint(idof)
                  endif               
               enddo
            enddo
         else
            do i=1,nope
               do j=1,3
                  if((lakon(ielem)(7:7).eq.'A').or.
     &               (lakon(ielem)(7:7).eq.'S').or.       
     &               (lakon(ielem)(7:7).eq.'E')) then
                     node=iponk2dto3d(konl(i))
c                write(*,*) 'objective_shapeener_tot node',ielem,i,node
                  else
                     node=konl(i)
                  endif
c                     write(*,*) 'objective_shapeener_tot',j,node,idof,
c     &vold(j,node),nactdof(j,node)
                  idof=nactdof(j,node)
                  if(idof.gt.0) then      
                     vec(idof)=vold(j,node)
c                     write(*,*) 'objective_shapeener_tot',j,node,idof,
c     &vold(j,node)
                  endif              
               enddo
            enddo
         endif
      enddo
!
!     Calculation of the total differential:    
!
c      do idesvar=1,ndesi
c         dgdx(idesvar,iobject)=0.d0
c      enddo
      do idesvar=1,ndesi
         do j=jqs(idesvar),jqs(idesvar+1)-1
            idof=irows(j)
            dgdx(idesvar,iobject)=dgdx(idesvar,iobject) 
     &            +vec(idof)*df(j) 
c     &           +vec(idof)
c            write(*,*) idesvar,j,idof,vec(idof),dgdx(idesvar,iobject)
         enddo
      enddo     
!      
      return
      end
