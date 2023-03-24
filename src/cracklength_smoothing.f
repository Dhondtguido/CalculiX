!     
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
      subroutine cracklength_smoothing(nnfront,isubsurffront,
     &     istartfront,iendfront,ifrontrel,costruc,dist,a,
     &     istartcrackfro,iendcrackfro,acrack,amin,nfront,ncrack,
     &     idist,nstep)
!     
!     smoothes the crack length
!     calculates the minimum crack length for each crack      
!     
      implicit none
!     
      integer i,j,nnfront,isubsurffront(*),istart,iend,istartfront(*),
     &     iendfront(*),ifrontrel(*),jforward,jbackward,n,nmax,node,
     &     istartcrackfro(*),iendcrackfro(*),nfront,ncrack,k,idist(*),
     &     nstep,m
!     
      real*8 xp,yp,zp,costruc(3,*),distforward,distbackward,dist(*),
     &     radius,sum,temp,acrack(*),amin(*),alpha,a(*)
!     
!     copy the crack length in temporary field a
!     
      do i=1,nfront
          a(i)=acrack(i)
      enddo
!     
      do i=1,nnfront
!     
!     target number of nodes in the circle of influence
!     
        nmax=max(1,int((iendfront(i)-istartfront(i)+1)/1.1d0))
        if(nmax.le.2) cycle
!     
!       next line introduced on 27th of Nov 2022 to limit the radius 
!       of influence     
!     
        nmax=min(nmax,5)
!     
        if(isubsurffront(i).eq.1) then
          istart=istartfront(i)
          iend=iendfront(i)
        else
          istart=istartfront(i)+1
          iend=iendfront(i)-1
        endif
!     
        do j=istart,iend
          xp=costruc(1,ifrontrel(j))
          yp=costruc(2,ifrontrel(j))
          zp=costruc(3,ifrontrel(j))
!     
!     actual forward node and its distance from node
!     actual backward node and its distance from node
!     actual number of nodes n in the circle of influence          
!     
          jforward=j
          jbackward=j
          n=1
          dist(j)=0.d0
          idist(1)=j
!     
!     calculating one position forward
!     
          jforward=jforward+1
          if(isubsurffront(i).eq.1) then
            if(jforward.gt.iendfront(i)) jforward=istartfront(i)
          endif
!     
          if((isubsurffront(i).eq.0).and.
     &         (jforward.gt.iendfront(i))) then
            distforward=1.d30
          else
            distforward=
     &           dsqrt((costruc(1,ifrontrel(jforward))-xp)**2+
     &           (costruc(2,ifrontrel(jforward))-yp)**2+
     &           (costruc(3,ifrontrel(jforward))-zp)**2)
          endif
!     
!     calculating one position backward
!     
          jbackward=jbackward-1
          if(isubsurffront(i).eq.1) then
            if(jbackward.lt.istartfront(i)) jbackward=iendfront(i)
          endif
!     
          if((isubsurffront(i).eq.0).and.
     &         (jbackward.lt.istartfront(i))) then
            distbackward=1.d30
          else
            distbackward=
     &           dsqrt((costruc(1,ifrontrel(jbackward))-xp)**2+
     &           (costruc(2,ifrontrel(jbackward))-yp)**2+
     &           (costruc(3,ifrontrel(jbackward))-zp)**2)
          endif
!     
!     start infinite loop
!     
          do
            if(distforward.lt.distbackward) then
              n=n+1
              idist(n)=jforward
              dist(jforward)=distforward
              radius=distforward
              if(n.eq.nmax) exit
!     
              jforward=jforward+1
!     
!     calculate the next forward position
!     
              if(isubsurffront(i).eq.1) then
                if(jforward.gt.iendfront(i)) jforward=istartfront(i)
              endif
!     
              if((isubsurffront(i).eq.0).and.
     &             (jforward.gt.iendfront(i))) then
                distforward=1.d30
              else
                distforward=
     &               dsqrt((costruc(1,ifrontrel(jforward))-xp)**2+
     &               (costruc(2,ifrontrel(jforward))-yp)**2+
     &               (costruc(3,ifrontrel(jforward))-zp)**2)
              endif
            else
              n=n+1
              idist(n)=jbackward
              dist(jbackward)=distbackward
              radius=distbackward
              if(n.eq.nmax) exit
!     
              jbackward=jbackward-1
!     
!     calculate the next backward position
!     
              if(isubsurffront(i).eq.1) then
                if(jbackward.lt.istartfront(i)) jbackward=iendfront(i)
              endif
!     
              if((isubsurffront(i).eq.0).and.
     &             (jbackward.lt.istartfront(i))) then
                distbackward=1.d30
              else
                distbackward=
     &               dsqrt((costruc(1,ifrontrel(jbackward))-xp)**2+
     &               (costruc(2,ifrontrel(jbackward))-yp)**2+
     &               (costruc(3,ifrontrel(jbackward))-zp)**2)
              endif
            endif
          enddo
!     
!     calculate the weighted mean
!     
            sum=0.d0
            temp=0.d0
            do k=1,n
              node=idist(k)
              alpha=1.d0-dist(node)/radius
              sum=sum+alpha
              temp=temp+alpha*a(node)
            enddo
            acrack(j)=temp/sum
        enddo
!     
        if(isubsurffront(i).eq.0) then
            acrack(istartfront(i))=acrack(istartfront(i)+1)
            acrack(iendfront(i))=acrack(iendfront(i)-1)
        endif
      enddo
!     
!     determine the minimum crack length for each crack
!     
c      write(*,*)
c      write(*,*) 'cracklength after smoothing'
c      write(*,*)
c      do i=1,ncrack
c        amin(i)=1.d30
c        do j=istartcrackfro(i),iendcrackfro(i)
c            amin(i)=min(amin(i),acrack(j))
c            write(*,*) 'cracklength ',j,m,acrack(j)
c        enddo
c      enddo
!     
      return
      end
