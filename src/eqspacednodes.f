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
      subroutine eqspacednodes(co,istartfront,iendfront,nnfront,
     &     ifrontprop,nk,nfront,ifronteq,charlen,
     &     istartfronteq,iendfronteq,nfronteq,acrackglob,ier,
     &     iendcrackfro,iincglob,iinc,dnglob,ncyctot)
!     
!     determine the mesh characteristic length for each front
!     
!     acrackglob, iincglob and dnglob are the crack length, the
!     increment number and the total number of cycles at the END
!     of the present increment and are therefore attached to the
!     propagated front, therefore they get values assigned in the
!     present routine; 
!      
      implicit none
!     
      integer i,k,m,n1,n2,istartfront(*),iendfront(*),iendcrackfro(*),
     &     nnfront,ifrontprop(*),nodesnum,ier,icrack,id,nk,nfront,
     &     ifronteq(*),kend,istartfronteq(*),iendfronteq(*),nfronteq,
     &     iincglob(*),iinc,ncyctot,nklim
!     
      real*8 co(3,*),dist,charlen(*),x(nfront),px,delta,x1,x2,
     &     acrackglob(*),dnglob(*)
!
      nklim=nk+2*nfront
!     
!     loop over all front(s) 
!     
      icrack=1
      nfronteq=0
      do i=1,nnfront
        istartfronteq(i)=nfronteq+1
!     
!     loop over all nodes belonging to the propagated front
!     
        k=1
        x(1)=0.d0
!     
        if(iendcrackfro(icrack).lt.istartfront(i)) then
          icrack=icrack+1
        endif
!        
        do m=istartfront(i),iendfront(i)-1
!     
!     distance between two adjacent propagated nodes
!     
          n1=ifrontprop(m)
          n2=ifrontprop(m+1)
          dist=dsqrt((co(1,n2)-co(1,n1))**2+
     &         (co(2,n2)-co(2,n1))**2+
     &         (co(3,n2)-co(3,n2))**2)
          k=k+1
          x(k)=x(k-1)+dist
        enddo
        kend=k
!
!     first node of front (position is not changed); done
!     for consistency in the node numbering along the crack front        
!
        nk=nk+1
        if(nk.gt.nklim) then
          write(*,*) '*ERROR in eqspacednodes: nfronteq > 2*nfront'
          ier=1
          return
        endif
        do k=1,3
          co(k,nk)=co(k,ifrontprop(istartfront(i)))
        enddo
        acrackglob(nk)=acrackglob(ifrontprop(istartfront(i)))
        iincglob(nk)=iinc+1
        dnglob(nk)=1.d0*ncyctot
        ifronteq(istartfronteq(i))=nk
!     
!     nodesnum is the new number of new nodes on the propagated front
!     
        nodesnum=nint(x(kend)/charlen(icrack))+1
        if(nk+nodesnum-1.gt.nklim) then
          write(*,*) '*ERROR in eqspacednodes: nfronteq > 2*nfront'
          ier=1
          return
        endif
        delta=x(kend)/(nodesnum-1)
!     
!     treating the nodes in between start and end
!     
        do m=1,nodesnum-2
          px=m*delta
          call ident(x,px,kend,id)
          nk=nk+1
          x1=(x(id+1)-px)/(x(id+1)-x(id))
          x2=1.d0-x1
          n1=ifrontprop(istartfront(i)-1+id)
          n2=ifrontprop(istartfront(i)+id)
          do k=1,3
            co(k,nk)=x1*co(k,n1)+x2*co(k,n2)
          enddo
          acrackglob(nk)=x1*acrackglob(n1)+x2*acrackglob(n2)
          iincglob(nk)=iinc+1
          dnglob(nk)=1.d0*ncyctot
          ifronteq(istartfronteq(i)+m)=nk
        enddo
!     
!     ifronteq contains the equivalent front nodes
!     
        nfronteq=nfronteq+nodesnum
        iendfronteq(i)=nfronteq
!
!     last node of front (position is not changed); done
!     for consistency in the node numbering along the crack front        
!
        nk=nk+1
        do k=1,3
          co(k,nk)=co(k,ifrontprop(iendfront(i)))
        enddo
        acrackglob(nk)=acrackglob(ifrontprop(iendfront(i)))
        iincglob(nk)=iinc+1
        dnglob(nk)=1.d0*ncyctot
        ifronteq(iendfronteq(i))=nk
!        
c        if(nfronteq.gt.(2*nfront)) then
c          write(*,*) '*ERROR in calccharlength: nfronteq.gt.2*nfront'
c          ier=1
c        endif
!
      enddo
      return
      end
      
      
