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
      subroutine actideactistr(set,nset,istartset,iendset,ialset,
     &           objectset,ipkon,iobject,ne,neinset,iponoel,inoel,
     &           nepar,nkinsetinv,nk)
!
!     deactivates the elements which are not adjacent to the nodes in
!     the STRESS objective function set
!
      implicit none
!
      character*81 objectset(5,*),set(*)
!
      integer i,j,k,nset,istartset(*),iendset(*),ialset(*),ipkon(*),
     &  iobject,ne,index,nelem,iponoel(*),inoel(2,*),neinset(*),
     &  nepar,nkinsetinv(*),nk,id
!
!     determining the nodes set corresponding to the STRESS
!     objective function
!
c      do i=1,nset
c         if(objectset(3,iobject).eq.set(i)) exit
c      enddo
      call cident81(set,objectset(3,iobject),nset,id)
      i=nset+1
      if(id.gt.0) then
        if(objectset(3,iobject).eq.set(id)) then
          i=id
        endif
      endif
!
      nepar=0
!
      if(i.le.nset) then
!
!        deactivate all elements
!
         do j=1,ne
            if(ipkon(j).lt.0) cycle
            ipkon(j)=-2-ipkon(j)
         enddo
!
!        reactivate the elements adjacent to the nodes in the
!        STRESS objective function set (the stress is extrapolated
!        to the nodes, therefore only those elements are needed to
!        which the nodes in the STRESS objective function belong)
!
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               index=iponoel(ialset(j))
               nkinsetinv(ialset(j))=1
               do
                  if(index.eq.0) exit 
                  nelem=inoel(1,index)
                  if(neinset(nelem).eq.0) then
                     ipkon(nelem)=-ipkon(nelem)-2
                     neinset(nelem)=1
                     nepar=nepar+1
                  endif
                  index=inoel(2,index)
               enddo
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  index=iponoel(k)
                  nkinsetinv(k)=1
                  do
                     if(index.eq.0) exit 
                     nelem=inoel(1,index)
                     if(neinset(nelem).eq.0) then
                        ipkon(nelem)=-ipkon(nelem)-2
                        neinset(nelem)=1
                        nepar=nepar+1
                     endif
                     index=inoel(2,index)
                  enddo
               enddo
            endif
         enddo
      else
!
!     all elements are taken into account
!
         do i=1,ne
            if(ipkon(i).lt.0) cycle
            neinset(i)=1
            nepar=nepar+1
         enddo
         do i=1,nk
            nkinsetinv(i)=1
         enddo
      endif
!
!     putting all active elements in ascending order in field
!     neinset
!
      nepar=0
      do i=1,ne
         if(neinset(i).eq.1) then
            nepar=nepar+1
            neinset(nepar)=i
         endif
      enddo
!
      return
      end
