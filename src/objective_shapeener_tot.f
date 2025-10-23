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
      subroutine objective_shapeener_tot(ne,kon,ipkon,lakon,
     &   fint,vold,iperturb,mi,nactdof,dgdx,df,ndesi,iobject,
     &   jqs,irows,vec,nod1st,nactdofinv)
!
      implicit none
!
      character*8 lakon(*)
!
      integer ndesi,iobject,idesvar,i,j,jqs(*),irows(*),idof,
     &   ne,ipkon(*),ielem,iperturb(*),indexe,kon(*),mi(*),
     &   nactdof(0:mi(2),*),node,nod1st(*),nactdofinv(*),
     &   idir,inode,mt
!      
      real*8 dgdx(ndesi,*),df(*),vec(*),vold(0:mi(2),*),fint(*),val
!
!     ----------------------------------------------------------------
!     Calculation of the total differential:
!     non-linear:  dgdx = dgdx + fint^(T) * ( df )
!     linear:      dgdx = dgdx + vold^(T) * ( df )
!     ----------------------------------------------------------------
!
!      do idesvar=1,ndesi
!       dgdx(idesvar,iobject)=0.d0
!      enddo
      mt=mi(2)+1
      do idesvar=1,ndesi
!	write(5,*) '*********************************'
!	write(5,*) idesvar
!	write(5,*) '*********************************'
        do j=jqs(idesvar),jqs(idesvar+1)-1
           idof=irows(j)
	   if(iperturb(2).eq.1) then
	      val=fint(idof)
	   else
	      inode=nactdofinv(idof)/mt+1
	      idir=nactdofinv(idof)-mt*(nactdofinv(idof)/mt)
	      val=vold(idir,inode)
	   endif
           dgdx(idesvar,iobject)=dgdx(idesvar,iobject)+val*df(j)
!	   write(5,101) inode,idir,idof,val,df(j),dgdx(idesvar,iobject)
        enddo
      enddo     
! 
! 101  format(1(3x,i10,3x,i1,3x,i10,3x,e14.7,3x,e14.7,3x,e14.7))
!    
      return
      end
