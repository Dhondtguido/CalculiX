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
      subroutine attach_3d(pneigh,pnode,nterms,ratio,dist,xil,etl,zel,
     &  loopa)
!
!     ataches node with coordinates in "pnode" to the face containing 
!     "nterms" nodes with coordinates in field "pneigh" (nterms < 20).
!     cave: the coordinates are stored in pneigh(1..3,*)
!
      implicit none
!
      integer nterms,i,j,k,m,imin,jmin,kmin,im,jm,km,n,loopa
!
      real*8 ratio(20),pneigh(3,20),pnode(3),a,xi(-1:1,-1:1,-1:1),
     &  et(-1:1,-1:1,-1:1),ze(-1:1,-1:1,-1:1),p(3),distmin,d1,dist,
     &  xil,etl,zel,dx(3),al,dummy
!
      n=3
!
      d1=1.d0
!
      xi(0,0,0)=0.d0
      et(0,0,0)=0.d0
      ze(0,0,0)=0.d0
      call distattach_3d(xi(0,0,0),et(0,0,0),ze(0,0,0),pneigh,pnode,a,p,
     &     ratio,nterms)
      distmin=a
      imin=0
      jmin=0
      kmin=0
!
      do m=1,loopa
!
!     initialisation
!
         d1=d1/10.d0
!     
         do i=-1,1
            do j=-1,1
               do k=-1,1
                  if((i.eq.0).and.(j.eq.0).and.(k.eq.0)) cycle
!
                  xi(i,j,k)=xi(0,0,0)+i*d1
                  et(i,j,k)=et(0,0,0)+j*d1
                  ze(i,j,k)=ze(0,0,0)+k*d1
!
!              check whether inside the (-1,1)x(-1,1)x(-1,1) domain
!
                  if((xi(i,j,k).le.1.d0).and.
     &               (xi(i,j,k).ge.-1.d0).and.
     &               (et(i,j,k).le.1.d0).and.
     &               (et(i,j,k).ge.-1.d0).and.
     &               (ze(i,j,k).le.1.d0).and.
     &               (ze(i,j,k).ge.-1.d0)) then
                     call distattach_3d(xi(i,j,k),et(i,j,k),ze(i,j,k),
     &                    pneigh,pnode,a,p,ratio,nterms)
!     
!     checking for smallest initial distance
!     
                     if(a.lt.distmin) then
                        distmin=a
                        imin=i
                        jmin=j
                        kmin=k
                     endif
                  endif
!
               enddo
            enddo
         enddo
!     
!     minimizing the distance from the face to the node
!     
         do
!     
!     exit if minimum found
!     
c            write(*,*) 'attach_3d',m,imin,jmin,kmin
            if((imin.eq.0).and.(jmin.eq.0).and.(kmin.eq.0)) exit
!
!           new center of 3x3x3 matrix
!
            xi(0,0,0)=xi(imin,jmin,kmin)
            et(0,0,0)=et(imin,jmin,kmin)
            ze(0,0,0)=ze(imin,jmin,kmin)
!
            im=imin
            jm=jmin
            km=kmin
!
            imin=0
            jmin=0
            kmin=0
!     
            do i=-1,1
               do j=-1,1
                  do k=-1,1
                     if((i+im.lt.-1).or.(i+im.gt.1).or.
     &                  (j+jm.lt.-1).or.(j+jm.gt.1).or.
     &                  (k+km.lt.-1).or.(k+km.gt.1)) then
!
                        xi(i,j,k)=xi(0,0,0)+i*d1
                        et(i,j,k)=et(0,0,0)+j*d1
                        ze(i,j,k)=ze(0,0,0)+k*d1
!
!              check whether inside the (-1,1)x(-1,1)x(-1,1) domain
!
                        if((xi(i,j,k).le.1.d0).and.
     &                     (xi(i,j,k).ge.-1.d0).and.
     &                     (et(i,j,k).le.1.d0).and.
     &                     (et(i,j,k).ge.-1.d0).and.
     &                     (ze(i,j,k).le.1.d0).and.
     &                     (ze(i,j,k).ge.-1.d0)) then
                           call distattach_3d(xi(i,j,k),et(i,j,k),
     &                          ze(i,j,k),pneigh,pnode,a,p,ratio,nterms)
!
!                       check for new minimum
!
                           if(a.lt.distmin) then
                              distmin=a
                              imin=i
                              jmin=j
                              kmin=k
                           endif
                        endif
!
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call distattach_3d(xi(0,0,0),et(0,0,0),ze(0,0,0),pneigh,pnode,a,p,
     &     ratio,nterms)
!
      do i=1,3
        pnode(i)=p(i)
      enddo
!
      dist=dsqrt(a)
!
      if((nterms.eq.20).or.(nterms.eq.8)) then
         xil=xi(0,0,0)
         etl=et(0,0,0)
         zel=ze(0,0,0)
      elseif((nterms.eq.4).or.(nterms.eq.10)) then
         xil=(xi(0,0,0)+1.d0)/2.d0
         etl=(et(0,0,0)+1.d0)/2.d0
         zel=(ze(0,0,0)+1.d0)/2.d0
         dx(1)=xil
         dx(2)=etl
         dx(3)=zel
         call insertsortd(dx,n)
c         call dsort(dx,iy,n,kflag)
         if(dx(3).gt.1.d-30) then
            al=dx(3)/(xil+etl+zel)
            xil=al*xil
            etl=al*etl
            zel=al*zel
         endif
      elseif((nterms.eq.6).or.(nterms.eq.15)) then
         xil=(xi(0,0,0)+1.d0)/2.d0
         etl=(et(0,0,0)+1.d0)/2.d0
         if(xil+etl.gt.1.d0) then
            dummy=xil
            xil=1.d0-etl
            etl=1.d0-dummy
         endif
         zel=ze(0,0,0)
      endif
!
      return
      end
      
