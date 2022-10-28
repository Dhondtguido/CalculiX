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
      subroutine createnodeneigh(nk,istartnk,ialnk,
     &   istartnneigh,ialnneigh,ichecknodes,lakon,ipkon,kon,
     &   nkinsetinv,neielemtot)                   
!
      implicit none
!
      character*8 lakon(*)
!
      integer nk,istartnneigh(*),ialnneigh(*),
     &   istartnk(*),ialnk(*),ifree,index,i,j,k,nea,neb,elem,
     &   ipkon(*),kon(*),ipos,nope,ichecknodes(*),node,
     &   inode,nkinsetinv(*),neielemtot,nka,nkb
!
!     determining all the OBJECTIVE nodes (and only those;
!     is verified by use of field nkinsetinv) of the 
!     neighboring elements of node i.
!     They are stored in ialnneigh(istartnneigh(i))..
!     ...up to..... ialnneigh(istartnneigh(i+1)-1)
!
      ifree=1
      do i=1,nk
!        
         istartnneigh(i)=ifree 
c         index=iponoel(i)
c         if(index.eq.0) cycle 
         nea=istartnk(i)
         neb=istartnk(i+1)-1
!   
         do j=nea,neb
!   
            elem=ialnk(j)
            ipos=ipkon(elem)
!
            if(lakon(elem)(4:4).eq.'8') then
               nope=8
            elseif(lakon(elem)(4:5).eq.'20') then
                nope=20
            elseif(lakon(elem)(4:5).eq.'10') then
               nope=10
            elseif(lakon(elem)(4:4).eq.'4') then
               nope=4
            elseif(lakon(elem)(4:4).eq.'6') then
               nope=6
            elseif(lakon(elem)(4:5).eq.'15') then
               nope=15
            endif
!
            do k=1,nope                     
               if(ichecknodes(kon(ipos+k)).eq.i) cycle
               if(nkinsetinv(kon(ipos+k)).eq.1) then
                  inode=kon(ipos+k)
                  ialnneigh(ifree)=kon(ipos+k)
                  ifree=ifree+1  
                  ichecknodes(kon(ipos+k))=i
               endif       
            enddo
         enddo
      enddo
      istartnneigh(nk+1)=ifree
!
!     determining an upper limit of the number of elements
!     to which the [objective nodes belonging to the elements
!     adjacent of node nk] belong
!
!     needed for allocation purposes
!
      neielemtot=0
      do i=1,nk
!
!        loop over all neighboring objective nodes of node i
!         
         nka=istartnneigh(i)
         nkb=istartnneigh(i+1)-1
         do j=nka,nkb
            node=ialnneigh(j)
!
!           neighboring elements
!
            neielemtot=neielemtot+istartnk(node+1)-istartnk(node)
         enddo
      enddo
!
      return
      end
