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
      subroutine stressintensity_smoothing(nnfront,isubsurffront,
     &     istartfront,iendfront,ifrontrel,costruc,dist,
     &     istartcrackfro,iendcrackfro,xkeq,phi,nfront,ncrack,
     &     dk,p,idist,phimax,nstep)
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
     &     radius,sum,temp,xkeq(nstep,*),phi(nstep,*),alpha,dk(nstep,*),
     &     p(nstep,*),tempp,pi,phimax,phiactmax,phiscale,c1
!     
      pi=4.d0*datan(1.d0)
      phiactmax=0.d0
!     
!     copy the crack length in temporary field a
!     
      do i=1,nfront
        do m=1,nstep
          dk(m,i)=xkeq(m,i)
          p(m,i)=phi(m,i)
        enddo
      enddo
!     
      do i=1,nnfront
!     
!       target number of nodes in the circle of influence
!     
        nmax=max(1,int((iendfront(i)-istartfront(i)+1)/1.1d0))
        if(nmax.le.2) cycle
!     
!       next line introduced on 27th of Nov 2022 to limit the radius 
!       of influence     
!     
        nmax=min(nmax,5)
!     
        istart=istartfront(i)
        iend=iendfront(i)
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
          do m=1,nstep
            sum=0.d0
            temp=0.d0
            tempp=0.d0
            do k=1,n
              node=idist(k)
              alpha=1.d0-dist(node)/radius
              sum=sum+alpha
              temp=temp+alpha*dk(m,node)
              tempp=tempp+alpha*p(m,node)
            enddo
            xkeq(m,j)=temp/sum
            phi(m,j)=tempp/sum
            phiactmax=max(phiactmax,dabs(phi(m,j)))
          enddo
        enddo
      enddo
!     
!     if the actual maximum deflection angle exceeds the
!     maximum allowed angle defined by the user all angles
!     are reduced proportionally.
!     
      if(phiactmax.gt.phimax) then
        phiscale=phimax/phiactmax
        do j=1,nfront
          do m=1,nstep
            phi(m,j)=phi(m,j)*phiscale
          enddo
        enddo
      endif
!     
c      write(*,*)
c      write(*,*) 'keq/phi after smoothing'
c      write(*,*)
c      c1=1.d0/dsqrt(1000.d0)
c      do i=1,ncrack
c        do m=1,nstep
c          do j=istartcrackfro(i),iendcrackfro(i)
c            write(*,*) 'keq/phi ',m,j,c1*xkeq(m,j),phi(m,j)*180.d0/pi
c          enddo
c          write(*,*)
c        enddo
c      enddo
!     
      return
      end

