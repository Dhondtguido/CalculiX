      
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
      subroutine cracklength(ncrack,istartcrackfro,iendcrackfro,
     &     co,istartcrackbou,iendcrackbou,costruc,ibounnod,xt,acrack,
     &     istartfront,iendfront,nnfront,isubsurffront,ifrontrel,
     &     ifront,posfront,doubleglob,integerglob,
     &     nproc,iinc,acrackglob,ier)
!     
!     determine for each crack front node
!     1. a crack length based on the previous crack length augmented
!     by the last crack increment
!     2. a crack length based on the
!     intersection of a plane orthogonal to the local tangent with
!     the boundary line of the crack mesh
!     
      implicit none
!     
      integer i,j,k,istartcrackfro(*),iendcrackfro(*),m,n1,n2,
     &     ncrack,istartcrackbou(*),iendcrackbou(*),ibounnod(*),
     &     istartfront(*),iendfront(*),isubsurffront(*),ifront(*),
     &     nnfront,jrel,ifrontrel(*),integerglob(*),nterms,nselect,
     &     node,nktet,nkon,nfield,nfaces,netet,nelem,ne,loopa,
     &     konl(20),istartset(1),iendset(1),ialset(1),iselect(1),
     &     imastset,nproc,iinc,ier
!     
      real*8 co(3,*),costruc(3,*),xt(3,*),xsec(3),x1,x2,ca,cb,cc,cd,a,
     &     acrack(*),p(3),q(3),r(3),coords(3),posfront(*),
     &     acrackglob(*),dlength,alambda,
     &     doubleglob(*),dummy(1),ratio(20),dist,cosang        
!     
      nktet=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=13
!     
      nselect=0
      imastset=0
      loopa=8 
!     
!     nproc=1: incremental summation of the crack length
!     nproc=2: crack length on basis of the intersection of a plane
!     orthogonal to the tangent with the opposite crack boundary
!     
!     crack length from last increment (incremental summation)
!     (nproc=1 in all but first increment)
!
      if((nproc.eq.1).and.(iinc.gt.0)) then
!     
c        write(*,*)
c        write(*,*) 'cracklength based on last increment'
c        write(*,*) '(nproc=1 in all but first increment) '
c        write(*,*)
!        
        do i=1,nnfront
          do j=istartfront(i),iendfront(i)-1
            node=ifront(j)
            acrack(j)=acrackglob(node)
          enddo
        enddo
      endif
!
!     crack length based on intersection with opposite crack front
!      
      if((nproc.ge.2).or.(iinc.eq.0)) then
!     
c        write(*,*)
c        write(*,*) 'cracklength based on intersection of crack boundary'
c        write(*,*) 'with plane orthogonal to tangent vector'
c        write(*,*) '(nproc=1 in the first increment or '
c        write(*,*) 'nproc=2 or first guess for nproc=3)'
c        write(*,*)
!     
        do i=1,ncrack
!     
!     loop over all nodes belonging to the crack front(s)
!     
          do j=istartcrackfro(i),iendcrackfro(i)          
!     
!     equation of plane through node and orthogonal to the local
!     tangent
!     
            ca=xt(1,j)
            cb=xt(2,j)
            cc=xt(3,j)
!     
!     position of front node j in ibounnod
!     
            jrel=ifrontrel(j)
            cd=-ca*costruc(1,jrel)-cb*costruc(2,jrel)-cc*costruc(3,jrel)
!     
!     loop over all nodes belonging to the crack boundary
!     
            a=1.d30
            do k=istartcrackbou(i),iendcrackbou(i)
              n1=k
              if(k.eq.iendcrackbou(i)) then
                n2=istartcrackbou(i)
              else
                n2=k+1
              endif
!     
!     segments adjacent to the node should not be considered
!     
              if((n1.eq.jrel).or.(n2.eq.jrel)) cycle
              x1=ca*costruc(1,n1)+cb*costruc(2,n1)+cc*costruc(3,n1)+cd
              x2=ca*costruc(1,n2)+cb*costruc(2,n2)+cc*costruc(3,n2)+cd
              if((x1*x2.le.0.d0).and.(dabs(x2-x1).gt.0)) then
!     
!     line segment is cut by plane
!     
                do m=1,3
                  xsec(m)=(x2*costruc(m,n1)-x1*costruc(m,n2))/(x2-x1)
                enddo
!     
!     calculate the crack length
!     
                a=min(a,dsqrt((xsec(1)-costruc(1,jrel))**2+
     &               (xsec(2)-costruc(2,jrel))**2+
     &               (xsec(3)-costruc(3,jrel))**2))
              endif
            enddo
            acrack(j)=a
c     write(*,*) 'cracklength ',j,acrack(j)
          enddo
        enddo
      endif
!     
!     calculate the relative length along the front
!     may be needed for the shape factor calculation in shapefactor.f        
!     
      do i=1,nnfront
        dlength=0.d0
        posfront(istartfront(i))=0.d0
        do j=istartfront(i),iendfront(i)-1
          n1=ifront(j)
          n2=ifront(j+1)
          dist=dsqrt((co(1,n2)-co(1,n1))**2+
     &         (co(2,n2)-co(2,n1))**2+
     &         (co(3,n2)-co(3,n2))**2)
          dlength=dlength+dist
          posfront(j+1)=dlength
        enddo
        do j=istartfront(i),iendfront(i)            
          posfront(j)=posfront(j)/dlength
        enddo
      enddo
!     
!     final adjustments:
!     1. for a subsurface crack the crack length has to be halved
!     2. for a surface crack the locations outside the structure
!     take the neighboring value inside the structure
!     
      do i=1,nnfront
        if(isubsurffront(i).eq.1) then
          do j=istartfront(i),iendfront(i)
            acrack(j)=acrack(j)/2.d0
          enddo
        else
!
!     if the crack length for the node on or just inside the
!     structure was not found this may be due to the fact that
!     the node is exactly on the free surface. In that case the
!     value from a node just more inside is taken
!          
          if(acrack(istartfront(i)+1).eq.1.d30) then
            acrack(istartfront(i)+1)=acrack(istartfront(i)+2)
          endif
          if(acrack(iendfront(i)-1).eq.1.d30) then
            acrack(iendfront(i)-1)=acrack(iendfront(i)-2)
          endif
!     
!     crack length of nodes adjacent and external to the
!     structure is taken from the internal neighbors
!
          acrack(istartfront(i))=acrack(istartfront(i)+1)
          acrack(iendfront(i))=acrack(iendfront(i)-1)
!
!     check whether any crack length was not found
!
          do j=istartfront(i),iendfront(i)
            if(acrack(j).eq.1.d30) then
              write(*,*) '*ERROR in cracklength: crack'
              write(*,*) '       length could not be determined'
              write(*,*) '       for node ',ifront(j)
              ier=1
            endif
          enddo
c          do j=istartfront(i),iendfront(i)
c              write(*,*) 'cracklength ',j,acrack(j)
c          enddo
        endif
      enddo
!     
      return
      end
