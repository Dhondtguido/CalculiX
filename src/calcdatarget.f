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
      subroutine calcdatarget(ifront,co,nnfront,istartfront,
     &     iendfront,isubsurffront,damax,datarget,acrack,nstep)
!     
!     calculate the crack propagation increment:
!     it is the minimum of:
!     - the user-defined increment
!     - one fifth of the minimum crack front curvature
!     - one fifth of the smallest crack length
!
      implicit none
!     
      integer i,j,ifront(*),nnfront,istartfront(*),nstep,m,
     &     isubsurffront(*),iendfront(*),istart,iend
!     
      real*8 datarget,damax,acrack(*),co(3,*),denominator,
     &     xp,yp,zp,xa,ya,za,xn,yn,zn,xlpa,xlan,xlnp,dd,rcur
!
      datarget=damax
!
!     loop over all fronts
!
      do i=1,nnfront
        if(isubsurffront(i).eq.1) then
          istart=istartfront(i)
          iend=iendfront(i)
        else
          istart=istartfront(i)+1
          iend=iendfront(i)-1
        endif
!
!     loop over nodes belonging to front
!
        do j=istart,iend
          if(j.eq.istart) then
!
!     previous node
!
            if(isubsurffront(i).eq.1) then
              xp=co(1,ifront(iend))
              yp=co(2,ifront(iend))
              zp=co(3,ifront(iend))
            else
              xp=co(1,ifront(j-1))
              yp=co(2,ifront(j-1))
              zp=co(3,ifront(j-1))
            endif
!
!     actual node
!
            xa=co(1,ifront(j))
            ya=co(2,ifront(j))
            za=co(3,ifront(j))
          else
!
!     new previous node is old actual node
!
            xp=xa
            yp=ya
            zp=za
!
!     new actual node is old next node
!
            xa=xn
            ya=yn
            za=zn
          endif
!
!     next node
!
          if((j.eq.iend).and.(isubsurffront(i).eq.1)) then
            xn=co(1,ifront(istart))
            yn=co(2,ifront(istart))
            zn=co(3,ifront(istart))
          else
            xn=co(1,ifront(j+1))
            yn=co(2,ifront(j+1))
            zn=co(3,ifront(j+1))
          endif
!
!     calculate the radius of a circle going through the previous,
!     actual and next node (formula of Heron)
!
          if(j.eq.istart) then
            xlpa=dsqrt((xa-xp)**2+(ya-yp)**2+(za-zp)**2)
          else
            xlpa=xlan
          endif
          xlan=dsqrt((xn-xa)**2+(yn-ya)**2+(zn-za)**2)
          xlnp=dsqrt((xp-xn)**2+(yp-yn)**2+(zp-zn)**2)
          dd=(xlpa+xlan+xlnp)/2.d0
!
!         radius of the circle
!
          denominator=4.d0*dsqrt(dd*dabs(dd-xlpa)*dabs(dd-xlan)*
     &         dabs(dd-xlnp))
          if(denominator.lt.1.d-10) then
            rcur=1.d10
          else
            rcur=xlpa*xlan*xlnp/denominator
          endif
c            
c          rcur=xlpa*xlan*xlnp/
c     &         (4.d0*dsqrt(dd*(dd-xlpa)*(dd-xlan)*(dd-xlnp)))
!
c     datarget=min(datarget,rcur/5.d0,acrack(j)/5.d0)
          datarget=min(datarget,rcur/5.d0)
          datarget=min(datarget,acrack(j)/5.d0)
        enddo
      enddo
!     
      return
      end
      
