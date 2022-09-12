!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine totalcontact(tieset,ntie,ne,ipkon,kon,
     &     lakon,islavsurf,itiefac,pmastsurf,ne0,nkon0)
!     
!     generate contact elements for the slave contact nodes
!     
      implicit none
!     
      character*8 lakon(*)
      character*81 tieset(3,*)
!     
      integer ntie,ifree,ipkon(*),kon(*),ne,nel,i,k,jj,ne0,
     &     ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),nelemm,
     &     jfacem,indexe,nface,nope,nodefm(9),ifaces,jfaces,ifacem,
     &     id,islavsurf(2,*),itiefac(2,*),nelems,m,mint2d,nopes,
     &     igauss,nopem,nodefs(9),indexf,nkon0
!     
      real*8 pmastsurf(6,*)
!     
!     nodes per face for hex elements
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
!     
!     nodes per face for tet elements
!     
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
!     
!     nodes per face for linear wedge elements
!     
      data ifacew1 /1,3,2,0,
     &     4,5,6,0,
     &     1,2,5,4,
     &     2,3,6,5,
     &     3,1,4,6/
!     
!     nodes per face for quadratic wedge elements
!     
      data ifacew2 /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!
      ifree=nkon0
      ne=ne0
!      
      igauss=0
!     
!     loop over all active contact ties
!     
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         do jj=itiefac(1,i),itiefac(2,i)
            ifaces=islavsurf(1,jj)
            nelems=int(ifaces/10)
            jfaces=ifaces-nelems*10            
!     
            if(lakon(nelems)(4:5).eq.'8R') then
               nopes=4
               nope=8
            elseif(lakon(nelems)(4:4).eq.'8') then
               nopes=4
               nope=8
            elseif(lakon(nelems)(4:6).eq.'20R') then
               nopes=8
               nope=20
            elseif(lakon(nelems)(4:5).eq.'20') then
               nopes=8
               nope=20
            elseif(lakon(nelems)(4:5).eq.'10') then
               nopes=6
               nope=10
            elseif(lakon(nelems)(4:4).eq.'4') then
               nopes=3
               nope=4
!     
!     treatment of wedge faces
!     
            elseif(lakon(nelems)(4:4).eq.'6') then
               nope=6
               if(jfaces.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
            elseif(lakon(nelems)(4:5).eq.'15') then
               nope=15
               if(jfaces.le.2) then
                  nopes=6
               else
                  nopes=8
               endif
            endif
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
            if((nope.eq.20).or.(nope.eq.8)) then
               do m=1,nopes
                  nodefs(m)=kon(ipkon(nelems)+ifaceq(m,jfaces))
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do m=1,nopes
                  nodefs(m)=kon(ipkon(nelems)+ifacet(m,jfaces))
               enddo
            elseif(nope.eq.15) then
               do m=1,nopes
                  nodefs(m)=kon(ipkon(nelems)+ifacew2(m,jfaces))
               enddo
            else
               do m=1,nopes
                  nodefs(m)=kon(ipkon(nelems)+ifacew1(m,jfaces))
               enddo
            endif
!     
            mint2d=islavsurf(2,jj+1)-islavsurf(2,jj)
            if(mint2d.eq.0) cycle
            indexf=islavsurf(2,jj)
!     
            do m=1,mint2d
               igauss=indexf+m
!     
!     identifying the element face to which the
!     triangle belongs
!     
               ifacem=int(pmastsurf(3,igauss))
!
!              if no opposite master face: cycle
!
               if(ifacem.eq.0) cycle
!               
               nelemm=int(ifacem/10.d0)
               jfacem=ifacem-10*nelemm
!     
               indexe=ipkon(nelemm)
               if(lakon(nelemm)(4:5).eq.'20') then
                  nopem=8
                  nface=6
               elseif(lakon(nelemm)(4:4).eq.'8') then
                  nopem=4
                  nface=6
               elseif(lakon(nelemm)(4:5).eq.'10') then
                  nopem=6
                  nface=4
               elseif(lakon(nelemm)(4:4).eq.'4') then
                  nopem=3
                  nface=4
               elseif(lakon(nelemm)(4:5).eq.'15') then
                  if(jfacem.le.2) then
                     nopem=6
                  else
                     nopem=8
                  endif
                  nface=5
                  nope=15
               elseif(lakon(nelemm)(4:4).eq.'6') then
                  if(jfacem.le.2) then
                     nopem=3
                  else
                     nopem=4
                  endif
                  nface=5
                  nope=6
               else
                  cycle
               endif
!     
!     determining the nodes of the master face
!     
               if(nface.eq.4) then
                  do k=1,nopem
                     nodefm(k)=kon(indexe+ifacet(k,jfacem))
                  enddo
               elseif(nface.eq.5) then
                  if(nope.eq.6) then
                     do k=1,nopem
                        nodefm(k)=kon(indexe+ifacew1(k,jfacem))
                     enddo
                  elseif(nope.eq.15) then
                     do k=1,nopem
                        nodefm(k)=kon(indexe+ifacew2(k,jfacem))
                     enddo
                  endif
               elseif(nface.eq.6) then
                  do k=1,nopem
                     nodefm(k)=kon(indexe+ifaceq(k,jfacem))
                  enddo
               endif
!     
!     generation of a contact spring element
!     
               ne=ne+1
               nel=ne
!     
               id=ifree+1
               ifree=ifree+nopem+nopes+4
!     
               ipkon(nel)=id
               lakon(nel)(1:7)='ESPRNGC'
               lakon(nel)(8:8)=char(nopem+48)
!     
               kon(id)=nopes+nopem
!     
               do k=1,nopem
                  kon(id+k)=nodefm(k) 
               enddo
               id=id+nopem
               do k=1,nopes
                  kon(id+k)=nodefs(k)
               enddo
               id=id+nopes
               kon(id+1)=igauss
               kon(id+2)=jj
               kon(id+3)=indexf+m
!     
            enddo               ! m
         enddo                  ! jj
      enddo                     ! ntie
!     
      return
      end
