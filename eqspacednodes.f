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
     &     nnfront,ifrontprop(*),nodesnum,ier,icrack,nk,nfront,
     &     ifronteq(*),istartfronteq(*),iendfronteq(*),nfronteq,
     &     iincglob(*),iinc,ncyctot,nklim,mm,j
!     
      real*8 co(3,*),dist,charlen(*),x1,x2,
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
        if(iendcrackfro(icrack).lt.istartfront(i)) then
          icrack=icrack+1
        endif
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
        m=0
        n1=0
        do mm=istartfront(i)+1,iendfront(i)
c          n1=ifrontprop(mm-1)
          if(n1.eq.0) n1=ifrontprop(mm-1)
          n2=ifrontprop(mm)
          dist=dsqrt((co(1,n2)-co(1,n1))**2+
     &         (co(2,n2)-co(2,n1))**2+
     &         (co(3,n2)-co(3,n1))**2)
          nodesnum=nint(dist/charlen(icrack))
          if(nodesnum.eq.0) cycle
          do j=1,nodesnum
            x2=(1.d0*j)/nodesnum
            x1=1.d0-x2
            nk=nk+1
            if(nk.gt.nklim) then
              write(*,*)
     &             '*ERROR in eqspacednodes: nfronteq>2*nfront'
              ier=1
              return
            endif
            m=m+1
            do k=1,3
              co(k,nk)=x1*co(k,n1)+x2*co(k,n2)
            enddo
            acrackglob(nk)=x1*acrackglob(n1)+x2*acrackglob(n2)
            iincglob(nk)=iinc+1
            dnglob(nk)=1.d0*ncyctot
            ifronteq(istartfronteq(i)+m)=nk
          enddo
          n1=n2
        enddo
!     
!       ifronteq contains the equivalent front nodes
!     
        nfronteq=nfronteq+m+1
        iendfronteq(i)=nfronteq
!
      enddo
      return
      end
      
      
