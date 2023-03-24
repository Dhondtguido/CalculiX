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
      subroutine calcstress(objectset,iobject,stn,inode,stressval)
!     
!     calculates the stress measure defined in obejctset
!     
      implicit none
!     
      character*81 objectset(5,*)
!     
      integer iobject,inode,j
!     
      real*8 stn(6,*),stressval,q,maxs,p(3),p_(3),sd(6),mean,i1,i2,i3,
     &     val1,val2,phi,pi
!     
      pi=4.d0*datan(1.d0)
!     
!     Calculate Mises stress
!     
      if((objectset(1,iobject)(1:11).eq.'MISESSTRESS').or.
     &     (objectset(1,iobject)(1:11).eq.'MODALSTRESS')) then
        q=-(stn(1,inode)+stn(2,inode)+stn(3,inode))/3.d0
        stressval=dsqrt(1.5d0*
     &       ((stn(1,inode)+q)**2+(stn(2,inode)+q)**2+
     &       (stn(3,inode)+q)**2
     &       +2.d0*(stn(4,inode)**2+stn(5,inode)**2+stn(6,inode)**2)))
!     
!     Calculate Principal stresses
!     
      elseif((objectset(1,iobject)(1:9).eq.'PS1STRESS').or.
     &       (objectset(1,iobject)(1:9).eq.'PS3STRESS')) then   
!     
!     check if stress tensor contains only zeros
!     
        maxs=0.d0
        do j=1,6
          maxs=max(maxs,dabs(stn(j,inode)))
        enddo
        if(maxs.eq.0.d0) then
          p(1)=0.d0
          p(2)=0.d0
          p(3)=0.d0
        else
          do j=1,6
            sd(j)=stn(j,inode)/maxs
          enddo
          mean=0.d0
          do j=1,3
            mean=mean+sd(j)
          enddo
!     
!     calculation of invariants
!     
          i1=-(sd(1)+sd(2)+sd(3))
          i2=(sd(1)*sd(2))+(sd(2)*sd(3))+(sd(3)*sd(1))-(sd(4)*sd(4))
     &         -(sd(5)*sd(5))-(sd(6)*sd(6))
          i3=-(sd(1)*sd(2)*sd(3)+2*sd(4)*sd(5)*sd(6)-sd(1)*sd(5)
     &         *sd(5)-sd(2)*sd(6)*sd(6)-sd(3)*sd(4)*sd(4))
!     
          if((abs(i1).le.0.d0).and.
     &         (abs(i2).le.0.d0).and.
     &         (abs(i3).le.0.d0)) then   
!     
!     if all principals 0
!     
            p(1)=0.d0
            p(2)=0.d0
            p(3)=0.d0
          elseif((abs(i2).lt.0.d0).and.(abs(i3).le.0.d0)) then 
!     
!     if just one principal is not zero
!     
            p(1)=0.d0
            p(2)=0.d0
            p(3)=0.d0
            if(-i1.lt.0.d0) then
              p_(3)=-i1*maxs
            else
              p_(1)=-i1*maxs
            endif
          elseif((abs((sd(1)-sd(2))/mean).lt.1.0e-05).and.
     &           (abs((sd(1)-sd(3))/mean).lt.1.0e-05)) then     
!     
!     case of hydrstatic stress state
!     
            if((abs(sd(4)/mean).lt.1.0e-05).and.
     &           (abs(sd(5)/mean).lt.1.0e-05).and.
     &           (abs(sd(6)/mean).lt.1.0e-05)) then
!     
!     case that shear stresses are zero 
!     --> principal stresses already available
!     
              i1=0.d0
              do j=1,6
                i1=i1+sd(j)*maxs
              enddo
              p(1)=i1
              p(2)=i1
              p(3)=i1
            else
!     
!     case that shear stresses are to zero 
!     --> principal stresses have to be calculated
!     
              val1=(i2-i1*i1/3.)/3.
              val2=(i1*i1*i1/13.5-i1*i2/3.+i3)*.5
!     
              val1=dsqrt(abs(val1))*abs(val2)/val2
              phi=val2/(val1*val1*val1)  
!     
              if(phi.gt.1.d0) then
                phi=0.d0
              elseif(phi.lt.-1.d0) then
                phi=pi
              else
                phi=acos(phi)  
              endif
              p_(1)=-2*val1*cos(phi/3.)-i1/3.
              p_(2)=2*val1*cos(pi/3.-phi/3.)-i1/3.
              p_(3)=2*val1*cos(pi/3.+phi/3.)-i1/3.
            endif               
          else
!     
!     all other cases   
!     
            val1=(i2-i1*i1/3.)/3.
            val2=(i1*i1*i1/13.5-i1*i2/3.+i3)*.5
!     
            val1=dsqrt(abs(val1))*abs(val2)/val2
            phi=val2/(val1*val1*val1)  
!     
            if(phi.gt.1.d0) then
              phi=0.d0
            elseif(phi.lt.-1.d0) then
              phi=pi
            else
              phi=acos(phi)  
            endif
            p_(1)=-2*val1*cos(phi/3.)-i1/3.
            p_(2)=2*val1*cos(pi/3.-phi/3.)-i1/3.
            p_(3)=2*val1*cos(pi/3.+phi/3.)-i1/3.                       
          endif
        endif      
!     
!     sorting of principal stress
!     
        do j=1,3
          p_(j)=p_(j)*maxs
        enddo
!     
!     1st greatest   
!     
        if((p_(1).ge.p_(2)).and.(p_(1).ge.p_(3))) then
          p(1)=p_(1)
          if(p_(2).gt.p_(3)) then
            p(2)=p_(2)
            p(3)=p_(3)
          else
            p(2)=p_(3)
            p(3)=p_(2)
          endif
!     
!     2nd greatest
!     
        elseif((p_(2).ge.p_(1)).and.(p_(2).ge.p_(3))) then     
          p(1)=p_(2)
          if(p_(1).gt.p_(3)) then
            p(2)=p_(1)
            p(3)=p_(3)
          else
            p(2)=p_(3)
            p(3)=p_(1)
          endif
!     
!     3rd greatest 
!     
        elseif((p_(3).ge.p_(1)).and.(p_(3).ge.p_(2))) then     
          p(1)=p_(3)
          if(p_(1).gt.p_(2)) then
            p(2)=p_(1)
            p(3)=p_(2)
          else
            p(2)=p_(2)
            p(3)=p_(1)
          endif
        else
          write(*,*) 'ERROR: calcualtion of principal stresses as'
          write(*,*) '       input for the KS-function failed'
          call exit(201)
        endif
!     
!     return of PS1 or PS3 stress value
!     
        if(objectset(1,iobject)(1:9).eq.'PS1STRESS') then
          stressval=p(1)
        elseif(objectset(1,iobject)(1:9).eq.'PS3STRESS') then
          stressval=p(3)
        endif
      endif
!     
      return
      end
