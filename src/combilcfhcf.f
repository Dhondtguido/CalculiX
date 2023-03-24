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
      subroutine combilcfhcf(tempf,stressf,stress,hcfstress,temp,
     &     nbounnod,mei,nstep)
!
!     calculates LCF+HCF and LCF-HCF and stores the resulting fields
!     in tempf and stressf
!
      implicit none
!
      integer i,j,nbounnod,mei(*),nstep
!
      real*8 tempf(2,*),stressf(6,2,*),hcfstress(6,*),stress(6,nstep,*),
     &     temp(nstep,*)
!
      do i=1,nbounnod
        do j=1,6
          stressf(j,1,i)=stress(j,mei(2),i)+hcfstress(j,i)
          stressf(j,2,i)=stress(j,mei(2),i)-hcfstress(j,i)
        enddo
        tempf(1,i)=temp(mei(2),i)
        tempf(2,i)=temp(mei(2),i)
      enddo
!
      return
      end

