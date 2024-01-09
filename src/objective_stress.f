!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine objective_stress(nodeset,istartset,iendset,ialset,
     &  nk,idesvarc,iobject,mi,g0,nobject,stn,objectset,expks)
!
!     calculates the Kreisselmeier-Steinhauser function for the
!     chosen stress measure
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer nk,istartset(*),iendset(*),ialset(*),nodeset,
     & idesvarc,iobject,mi(*),j,k,nobject,idesvar
!
      real*8 g0(nobject),stn(6,*),p,rho,xstress,argument,expks,
     & stressval
!
      idesvar=idesvarc+1
!
!     reading rho and the mean stress for the Kreisselmeier-Steinhauser
!     function
!
      read(objectset(2,iobject)(41:60),'(f20.0)') rho
      read(objectset(2,iobject)(61:80),'(f20.0)') xstress
!
      g0(iobject)=0.d0
!
!     check for the existence of a set, else take the complete mesh
!
      if(nodeset.eq.0) then
         do j=1,nk
	    call calcstress(objectset,iobject,stn,j,stressval)
	    argument=rho*stressval/xstress
            if(argument.gt.600.d0) then
               write(*,*) '*ERROR in objective_stress: argument'
               write(*,*) '       ',argument,
     &            ' of exponential function is too big;'
               write(*,*) 
     &          '       choose other Kreisselmeier-Steinhauser'
               write(*,*) '       coefficients'
               call exit(201)
            endif
            g0(iobject)=g0(iobject)+dexp(argument)
         enddo
         expks=g0(iobject)
         g0(iobject)=dlog(g0(iobject))/rho
      else
         do j=istartset(nodeset),iendset(nodeset)
            if(ialset(j).gt.0) then
               k=ialset(j)
               call calcstress(objectset,iobject,stn,k,stressval)
	       argument=rho*stressval/xstress
               if(argument.gt.600.d0) then
                  write(*,*) '*ERROR in objective_stress: argument'
                  write(*,*) '       ',argument,
     &            ' of exponential function is too big;'
                  write(*,*) 
     &            '       choose other Kreisselmeier-Steinhauser'
                  write(*,*) '       coefficients'
                  call exit(201)
               endif
               g0(iobject)=g0(iobject)+dexp(argument)
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
		  call calcstress(objectset,iobject,stn,k,stressval)
	          argument=rho*stressval/xstress
                  if(argument.gt.600.d0) then
                     write(*,*) '*ERROR in objective_stress: argument'
                     write(*,*) '       ',argument,
     &               ' of exponential function is too big;'
                     write(*,*) 
     &               '       choose other Kreisselmeier-Steinhauser'
                     write(*,*) '       coefficients'
                     call exit(201)
                  endif
                  g0(iobject)=g0(iobject)+dexp(argument)
               enddo
            endif
         enddo
         expks=g0(iobject)
         g0(iobject)=dlog(g0(iobject))/rho
      endif
!     
      return
      end
