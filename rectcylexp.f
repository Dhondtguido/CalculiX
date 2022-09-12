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
      subroutine rectcylexp(co,v,fn,stn,qfn,een,cs,nkt,icntrl,t,filab,
     &  imag,mi,iznode,nznode,nsectors,nk,emn)
!
!     special version of routine rectcyl for use in expand.c
!
!     icntrl=2:  rectangular to cylindrical coordinates for fields
!                v,fn,stn, een and emn
!     icntrl=-2: cylindrical to rectangular coordinates for fields
!                v,fn,stn, een and emn
!
!     for icntrl=2 the imaginary part is extra taken into account if
!     imag=1
!
!     nk: number of nodes in one segment
!     nkt: number of nodes in 360°
!
      implicit none
!
      character*87 filab(*)
      integer i,j,nkt,icntrl,imag,mi(*),iznode(*),nznode,nsectors,nk,
     &  ii,jj,node
      real*8 co(3,*),v(0:mi(2),*),fn(0:mi(2),*),stn(6,*),een(6,*),
     &  a(3,3),xr,xt,xz,b(3,3),cs(17,*),t(3),qfn(3,*),csab(7),emn(6,*)
!
      do i=1,7
         csab(i)=cs(5+i,1)
      enddo
!
      if(icntrl.eq.2) then
         do ii=1,nznode
            i=iznode(ii)
            j=i
            call transformatrix(csab,co(1,i),a)
!
            if((filab(1)(1:3).eq.'U  ').or.
     &         (filab(11)(1:4).eq.'PU'))  then 
               xr=v(1,j)*a(1,1)+v(2,j)*a(2,1)+v(3,j)*a(3,1)
               xt=v(1,j)*a(1,2)+v(2,j)*a(2,2)+v(3,j)*a(3,2)
               xz=v(1,j)*a(1,3)+v(2,j)*a(2,3)+v(3,j)*a(3,3)
               v(1,j)=xr
               v(2,j)=xt
               v(3,j)=xz
            endif
!
            if((filab(3)(1:4).eq.'S   ').or.
     &         (filab(18)(1:4).eq.'PHS ')) then
               b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(2,1)+stn(5,j)*a(3,1)
               b(1,2)=stn(1,j)*a(1,2)+stn(4,j)*a(2,2)+stn(5,j)*a(3,2)
               b(1,3)=stn(1,j)*a(1,3)+stn(4,j)*a(2,3)+stn(5,j)*a(3,3)
               b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(2,1)+stn(6,j)*a(3,1)
               b(2,2)=stn(4,j)*a(1,2)+stn(2,j)*a(2,2)+stn(6,j)*a(3,2)
               b(2,3)=stn(4,j)*a(1,3)+stn(2,j)*a(2,3)+stn(6,j)*a(3,3)
               b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(2,1)+stn(3,j)*a(3,1)
               b(3,2)=stn(5,j)*a(1,2)+stn(6,j)*a(2,2)+stn(3,j)*a(3,2)
               b(3,3)=stn(5,j)*a(1,3)+stn(6,j)*a(2,3)+stn(3,j)*a(3,3)
!
               stn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
               stn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
               stn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
               stn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
               stn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
               stn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
            endif
!
            if(filab(4)(1:4).eq.'E   ') then
               b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(2,1)+een(5,j)*a(3,1)
               b(1,2)=een(1,j)*a(1,2)+een(4,j)*a(2,2)+een(5,j)*a(3,2)
               b(1,3)=een(1,j)*a(1,3)+een(4,j)*a(2,3)+een(5,j)*a(3,3)
               b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(2,1)+een(6,j)*a(3,1)
               b(2,2)=een(4,j)*a(1,2)+een(2,j)*a(2,2)+een(6,j)*a(3,2)
               b(2,3)=een(4,j)*a(1,3)+een(2,j)*a(2,3)+een(6,j)*a(3,3)
               b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(2,1)+een(3,j)*a(3,1)
               b(3,2)=een(5,j)*a(1,2)+een(6,j)*a(2,2)+een(3,j)*a(3,2)
               b(3,3)=een(5,j)*a(1,3)+een(6,j)*a(2,3)+een(3,j)*a(3,3)
!
               een(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
               een(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
               een(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
               een(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
               een(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
               een(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
            endif
!
            if(filab(5)(1:4).eq.'RF  ') then
               xr=fn(1,j)*a(1,1)+fn(2,j)*a(2,1)+fn(3,j)*a(3,1)
               xt=fn(1,j)*a(1,2)+fn(2,j)*a(2,2)+fn(3,j)*a(3,2)
               xz=fn(1,j)*a(1,3)+fn(2,j)*a(2,3)+fn(3,j)*a(3,3)
               fn(1,j)=xr
               fn(2,j)=xt
               fn(3,j)=xz
            endif
!
            if(filab(9)(1:4).eq.'HFL ') then
               xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(2,1)+qfn(3,j)*a(3,1)
               xt=qfn(1,j)*a(1,2)+qfn(2,j)*a(2,2)+qfn(3,j)*a(3,2)
               xz=qfn(1,j)*a(1,3)+qfn(2,j)*a(2,3)+qfn(3,j)*a(3,3)
               qfn(1,j)=xr
               qfn(2,j)=xt
               qfn(3,j)=xz
            endif
!
            if(filab(32)(1:4).eq.'ME  ') then
               b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(2,1)+emn(5,j)*a(3,1)
               b(1,2)=emn(1,j)*a(1,2)+emn(4,j)*a(2,2)+emn(5,j)*a(3,2)
               b(1,3)=emn(1,j)*a(1,3)+emn(4,j)*a(2,3)+emn(5,j)*a(3,3)
               b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(2,1)+emn(6,j)*a(3,1)
               b(2,2)=emn(4,j)*a(1,2)+emn(2,j)*a(2,2)+emn(6,j)*a(3,2)
               b(2,3)=emn(4,j)*a(1,3)+emn(2,j)*a(2,3)+emn(6,j)*a(3,3)
               b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(2,1)+emn(3,j)*a(3,1)
               b(3,2)=emn(5,j)*a(1,2)+emn(6,j)*a(2,2)+emn(3,j)*a(3,2)
               b(3,3)=emn(5,j)*a(1,3)+emn(6,j)*a(2,3)+emn(3,j)*a(3,3)
!
               emn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
               emn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
               emn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
               emn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
               emn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
               emn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
            endif
!
!           imaginary part for cyclic symmetry frequency calculations
!
            if(imag.eq.1) then
!
               j=i+nkt
!
               if((filab(1)(1:3).eq.'U  ').or.
     &            (filab(11)(1:4).eq.'PU'))  then 
                  xr=v(1,j)*a(1,1)+v(2,j)*a(2,1)+v(3,j)*a(3,1)
                  xt=v(1,j)*a(1,2)+v(2,j)*a(2,2)+v(3,j)*a(3,2)
                  xz=v(1,j)*a(1,3)+v(2,j)*a(2,3)+v(3,j)*a(3,3)
                  v(1,j)=xr
                  v(2,j)=xt
                  v(3,j)=xz
               endif
!
               if((filab(3)(1:4).eq.'S   ').or.
     &            (filab(18)(1:4).eq.'PHS ')) then
                  b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(2,1)+stn(5,j)*a(3,1)
                  b(1,2)=stn(1,j)*a(1,2)+stn(4,j)*a(2,2)+stn(5,j)*a(3,2)
                  b(1,3)=stn(1,j)*a(1,3)+stn(4,j)*a(2,3)+stn(5,j)*a(3,3)
                  b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(2,1)+stn(6,j)*a(3,1)
                  b(2,2)=stn(4,j)*a(1,2)+stn(2,j)*a(2,2)+stn(6,j)*a(3,2)
                  b(2,3)=stn(4,j)*a(1,3)+stn(2,j)*a(2,3)+stn(6,j)*a(3,3)
                  b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(2,1)+stn(3,j)*a(3,1)
                  b(3,2)=stn(5,j)*a(1,2)+stn(6,j)*a(2,2)+stn(3,j)*a(3,2)
                  b(3,3)=stn(5,j)*a(1,3)+stn(6,j)*a(2,3)+stn(3,j)*a(3,3)
!
                  stn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
                  stn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
                  stn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
                  stn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
                  stn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
                  stn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
               endif
!
               if(filab(4)(1:4).eq.'E   ') then
                  b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(2,1)+een(5,j)*a(3,1)
                  b(1,2)=een(1,j)*a(1,2)+een(4,j)*a(2,2)+een(5,j)*a(3,2)
                  b(1,3)=een(1,j)*a(1,3)+een(4,j)*a(2,3)+een(5,j)*a(3,3)
                  b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(2,1)+een(6,j)*a(3,1)
                  b(2,2)=een(4,j)*a(1,2)+een(2,j)*a(2,2)+een(6,j)*a(3,2)
                  b(2,3)=een(4,j)*a(1,3)+een(2,j)*a(2,3)+een(6,j)*a(3,3)
                  b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(2,1)+een(3,j)*a(3,1)
                  b(3,2)=een(5,j)*a(1,2)+een(6,j)*a(2,2)+een(3,j)*a(3,2)
                  b(3,3)=een(5,j)*a(1,3)+een(6,j)*a(2,3)+een(3,j)*a(3,3)
!
                  een(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
                  een(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
                  een(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
                  een(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
                  een(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
                  een(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
               endif
!
               if(filab(5)(1:4).eq.'RF  ') then
                  xr=fn(1,j)*a(1,1)+fn(2,j)*a(2,1)+fn(3,j)*a(3,1)
                  xt=fn(1,j)*a(1,2)+fn(2,j)*a(2,2)+fn(3,j)*a(3,2)
                  xz=fn(1,j)*a(1,3)+fn(2,j)*a(2,3)+fn(3,j)*a(3,3)
                  fn(1,j)=xr
                  fn(2,j)=xt
                  fn(3,j)=xz
               endif
!
               if(filab(9)(1:4).eq.'HFL ') then
                  xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(2,1)+qfn(3,j)*a(3,1)
                  xt=qfn(1,j)*a(1,2)+qfn(2,j)*a(2,2)+qfn(3,j)*a(3,2)
                  xz=qfn(1,j)*a(1,3)+qfn(2,j)*a(2,3)+qfn(3,j)*a(3,3)
                  qfn(1,j)=xr
                  qfn(2,j)=xt
                  qfn(3,j)=xz
               endif
!
               if(filab(32)(1:4).eq.'ME  ') then
                  b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(2,1)+emn(5,j)*a(3,1)
                  b(1,2)=emn(1,j)*a(1,2)+emn(4,j)*a(2,2)+emn(5,j)*a(3,2)
                  b(1,3)=emn(1,j)*a(1,3)+emn(4,j)*a(2,3)+emn(5,j)*a(3,3)
                  b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(2,1)+emn(6,j)*a(3,1)
                  b(2,2)=emn(4,j)*a(1,2)+emn(2,j)*a(2,2)+emn(6,j)*a(3,2)
                  b(2,3)=emn(4,j)*a(1,3)+emn(2,j)*a(2,3)+emn(6,j)*a(3,3)
                  b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(2,1)+emn(3,j)*a(3,1)
                  b(3,2)=emn(5,j)*a(1,2)+emn(6,j)*a(2,2)+emn(3,j)*a(3,2)
                  b(3,3)=emn(5,j)*a(1,3)+emn(6,j)*a(2,3)+emn(3,j)*a(3,3)
!
                  emn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
                  emn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
                  emn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
                  emn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
                  emn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
                  emn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
               endif
            endif
         enddo
      elseif(icntrl.eq.-2) then
         do ii=1,nznode
            node=iznode(ii)
            do jj=1,nsectors
               i=node+(jj-1)*nk
               j=i
               call transformatrix(csab,co(1,i),a)
!     
               if((filab(1)(1:3).eq.'U  ').or.
     &              (filab(11)(1:4).eq.'PU'))  then 
                  xr=v(1,j)*a(1,1)+v(2,j)*a(1,2)+v(3,j)*a(1,3)
                  xt=v(1,j)*a(2,1)+v(2,j)*a(2,2)+v(3,j)*a(2,3)
                  xz=v(1,j)*a(3,1)+v(2,j)*a(3,2)+v(3,j)*a(3,3)
                  v(1,j)=xr
                  v(2,j)=xt
                  v(3,j)=xz
               endif
!     
               if((filab(3)(1:4).eq.'S   ').or.
     &              (filab(18)(1:4).eq.'PHS ')) then
                  b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(1,2)+stn(5,j)*a(1,3)
                  b(1,2)=stn(1,j)*a(2,1)+stn(4,j)*a(2,2)+stn(5,j)*a(2,3)
                  b(1,3)=stn(1,j)*a(3,1)+stn(4,j)*a(3,2)+stn(5,j)*a(3,3)
                  b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(1,2)+stn(6,j)*a(1,3)
                  b(2,2)=stn(4,j)*a(2,1)+stn(2,j)*a(2,2)+stn(6,j)*a(2,3)
                  b(2,3)=stn(4,j)*a(3,1)+stn(2,j)*a(3,2)+stn(6,j)*a(3,3)
                  b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(1,2)+stn(3,j)*a(1,3)
                  b(3,2)=stn(5,j)*a(2,1)+stn(6,j)*a(2,2)+stn(3,j)*a(2,3)
                  b(3,3)=stn(5,j)*a(3,1)+stn(6,j)*a(3,2)+stn(3,j)*a(3,3)
!     
                  stn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
                  stn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
                  stn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
                  stn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
                  stn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
                  stn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
               endif
!     
               if(filab(4)(1:4).eq.'E   ') then
                  b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(1,2)+een(5,j)*a(1,3)
                  b(1,2)=een(1,j)*a(2,1)+een(4,j)*a(2,2)+een(5,j)*a(2,3)
                  b(1,3)=een(1,j)*a(3,1)+een(4,j)*a(3,2)+een(5,j)*a(3,3)
                  b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(1,2)+een(6,j)*a(1,3)
                  b(2,2)=een(4,j)*a(2,1)+een(2,j)*a(2,2)+een(6,j)*a(2,3)
                  b(2,3)=een(4,j)*a(3,1)+een(2,j)*a(3,2)+een(6,j)*a(3,3)
                  b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(1,2)+een(3,j)*a(1,3)
                  b(3,2)=een(5,j)*a(2,1)+een(6,j)*a(2,2)+een(3,j)*a(2,3)
                  b(3,3)=een(5,j)*a(3,1)+een(6,j)*a(3,2)+een(3,j)*a(3,3)
!     
                  een(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
                  een(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
                  een(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
                  een(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
                  een(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
                  een(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
               endif
!     
               if(filab(5)(1:4).eq.'RF  ') then
                  xr=fn(1,j)*a(1,1)+fn(2,j)*a(1,2)+fn(3,j)*a(1,3)
                  xt=fn(1,j)*a(2,1)+fn(2,j)*a(2,2)+fn(3,j)*a(2,3)
                  xz=fn(1,j)*a(3,1)+fn(2,j)*a(3,2)+fn(3,j)*a(3,3)
                  fn(1,j)=xr
                  fn(2,j)=xt
                  fn(3,j)=xz
               endif
!     
               if(filab(9)(1:4).eq.'HFL ') then
                  xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(1,2)+qfn(3,j)*a(1,3)
                  xt=qfn(1,j)*a(2,1)+qfn(2,j)*a(2,2)+qfn(3,j)*a(2,3)
                  xz=qfn(1,j)*a(3,1)+qfn(2,j)*a(3,2)+qfn(3,j)*a(3,3)
                  qfn(1,j)=xr
                  qfn(2,j)=xt
                  qfn(3,j)=xz
               endif
!     
               if(filab(32)(1:4).eq.'ME  ') then
                  b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(1,2)+emn(5,j)*a(1,3)
                  b(1,2)=emn(1,j)*a(2,1)+emn(4,j)*a(2,2)+emn(5,j)*a(2,3)
                  b(1,3)=emn(1,j)*a(3,1)+emn(4,j)*a(3,2)+emn(5,j)*a(3,3)
                  b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(1,2)+emn(6,j)*a(1,3)
                  b(2,2)=emn(4,j)*a(2,1)+emn(2,j)*a(2,2)+emn(6,j)*a(2,3)
                  b(2,3)=emn(4,j)*a(3,1)+emn(2,j)*a(3,2)+emn(6,j)*a(3,3)
                  b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(1,2)+emn(3,j)*a(1,3)
                  b(3,2)=emn(5,j)*a(2,1)+emn(6,j)*a(2,2)+emn(3,j)*a(2,3)
                  b(3,3)=emn(5,j)*a(3,1)+emn(6,j)*a(3,2)+emn(3,j)*a(3,3)
!     
                  emn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
                  emn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
                  emn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
                  emn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
                  emn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
                  emn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
               endif
!     
            enddo
         enddo
      endif
!
      return
      end













