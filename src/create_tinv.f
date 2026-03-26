!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
!
!     Generate local transformation matrix needed for quad-quad mortar method
!     see phd-thesis Sitzmann equation (4.3)
!     Author: Saskia Sitzmann
!
!    [out]    contr       field containing T_e contributions for current face
!    [out]    krow     (i)  row  of contribution(i)
!    [out]    kcol     (i)  column of contribution(i)
!    [out]    icounter    counter variable for contr
!    [in]     lface	 current slave face
!
      subroutine create_tinv(ipkon,kon,lakon,islavsurf,
     &     contr,krow,kcol,icounter,lface)
!     
!     Generate local transformation matrix T
!     
!     Author: Sitzmann,Saskia ;
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer ipkon(*),kon(*),konl(20),islavsurf(2,*),
     &     icounter,krow(*),kcol(*),j,nope,
     &     ifaces,nelems,jfaces,m,nopes,lface,
     &     ifac,getlocno,nodes(8),modf,idummy
!     
      real*8 contr(*),alpha
!     
      alpha=1.d0/5.d0
      icounter=0
      ifaces = islavsurf(1,lface)
      nelems = int(ifaces/10)
      jfaces = ifaces - nelems*10
!     
      call getnumberofnodes(nelems,jfaces,lakon,nope,
     &     nopes,idummy)
!     
      do j=1,nope
         konl(j)=kon(ipkon(nelems)+j)
      enddo
      do m=1,nopes
         ifac=getlocno(m,jfaces,nope)
         nodes(m)=konl(ifac)
      enddo
      if(nopes.eq.8) then
         do j=1,4
            icounter=icounter+1
            contr(icounter)=1.d0
            krow(icounter)=nodes(j)
            kcol(icounter)=nodes(j)
         enddo
         do j=5,8
            icounter=icounter+1
            contr(icounter)=1.d0/(1.d0-2.d0*alpha)
            krow(icounter)=nodes(j)
            kcol(icounter)=nodes(j)
        enddo
        do j=1,4
          icounter=icounter+1
          contr(icounter)=-alpha/(1-2.d0*alpha)
          krow(icounter)=nodes(j+4)
          kcol(icounter)=nodes(j)
          icounter=icounter+1
          contr(icounter)=-alpha/(1-2.d0*alpha)
          krow(icounter)=nodes(modf(4,j-1)+4)
          kcol(icounter)=nodes(j)
       enddo
       
      elseif(nopes.eq.4) then
         do j=1,4
            icounter=icounter+1
            contr(icounter)=1.d0
            krow(icounter)=nodes(j)
            kcol(icounter)=nodes(j)
         enddo                 
      elseif(nopes.eq.6) then
         do j=1,3
            icounter=icounter+1
            contr(icounter)=1.d0
            krow(icounter)=nodes(j)
            kcol(icounter)=nodes(j)
         enddo
         do j=4,6
            icounter=icounter+1
            contr(icounter)=1/(1.d0-2.d0*alpha)
            krow(icounter)=nodes(j)
            kcol(icounter)=nodes(j)
         enddo
         do j=1,3
            icounter=icounter+1
            contr(icounter)=-alpha/(1-2.d0*alpha)
            krow(icounter)=nodes(j+3)
            kcol(icounter)=nodes(j)
            icounter=icounter+1
            contr(icounter)=-alpha/(1-2.d0*alpha)
            krow(icounter)=nodes(modf(3,j-1)+3)
            kcol(icounter)=nodes(j)
         enddo
      else
         do j=1,3
            icounter=icounter+1
            contr(icounter)=1.d0
            krow(icounter)=nodes(j)
            kcol(icounter)=nodes(j)
         enddo
      endif 
!     
      return
      end
      
