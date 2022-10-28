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
      subroutine printoutnode(prlab,v,t1,fn,ithermal,ii,node,
     &  rftot,trab,inotr,ntrans,co,mi,veold)
!
!     stores results in the .dat file
!
      implicit none
!
      character*1 local
      character*6 prlab(*)
!
      integer ithermal(*),node,ii,j,inotr(2,*),ntrans,mi(*)
!
      real*8 v(0:mi(2),*),t1(*),fn(0:mi(2),*),rftot(0:3),trab(7,*),
     &  co(3,*),a(3,3),veold(0:mi(2),*)
!
      local='L'
!
      if(prlab(ii)(1:4).eq.'U   ') then
        if((ntrans.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
            write(5,'(i10,1p,6(1x,e13.6))') node,
     &           (v(j,node),j=1,mi(2))
         elseif(inotr(1,node).eq.0) then
            write(5,'(i10,1p,6(1x,e13.6))') node,
     &           (v(j,node),j=1,mi(2))
         elseif(mi(2).eq.3) then
            call transformatrix(trab(1,inotr(1,node)),co(1,node),a)
            write(5,'(i10,1p,3(1x,e13.6),1x,a1)') node,
     &          v(1,node)*a(1,1)+v(2,node)*a(2,1)+v(3,node)*a(3,1),
     &          v(1,node)*a(1,2)+v(2,node)*a(2,2)+v(3,node)*a(3,2),
     &          v(1,node)*a(1,3)+v(2,node)*a(2,3)+v(3,node)*a(3,3),
     &          local
         else
            write(*,*) '*WARNING in printoutnode:'
            write(*,*) '         for output purposes only 4, 5 or 6'
            write(*,*) '         degrees of freedom are allowed'
            write(*,*) '         for generalized vectors;'
            write(*,*) '         actual degrees of freedom = ',mi(2)
            write(*,*) '         output request ist not performed;'
         endif
      elseif(prlab(ii)(1:4).eq.'V   ') then
         if((ntrans.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
            write(5,'(i10,1p,3(1x,e13.6))') node,
     &           (veold(j,node),j=1,3)
         elseif(inotr(1,node).eq.0) then
            write(5,'(i10,1p,3(1x,e13.6))') node,
     &           (veold(j,node),j=1,3)
         else
            call transformatrix(trab(1,inotr(1,node)),co(1,node),a)
            write(5,'(i10,1p,3(1x,e13.6),1x,a1)') node,
     &          veold(1,node)*a(1,1)+veold(2,node)*a(2,1)+
     &          veold(3,node)*a(3,1),
     &          veold(1,node)*a(1,2)+veold(2,node)*a(2,2)+
     &          veold(3,node)*a(3,2),
     &          veold(1,node)*a(1,3)+veold(2,node)*a(2,3)+
     &          veold(3,node)*a(3,3),
     &          local
         endif
      elseif((prlab(ii)(1:4).eq.'NT  ').or.
     &       (prlab(ii)(1:4).eq.'TS  ')) then
         if(ithermal(1).le.1) then
            write(5,'(i10,1x,1p,e13.6)') node,
     &           t1(node)
         else
            write(5,'(i10,1x,1p,e13.6)') node,
     &           v(0,node)
         endif
      elseif(prlab(ii)(1:4).eq.'PS  ') then
         write(5,'(i10,1x,1p,e13.6)') node,
     &        v(4,node)
      elseif(prlab(ii)(1:4).eq.'PN  ') then
            write(5,'(i10,1x,1p,e13.6)') node,
     &           v(2,node)
      elseif(prlab(ii)(1:4).eq.'MF  ') then
         write(5,'(i10,1x,1p,e13.6)') node,
     &        v(1,node)
      elseif(prlab(ii)(1:4).eq.'RF  ') then
         do j=1,3
            rftot(j)=rftot(j)+fn(j,node)
         enddo
         if(prlab(ii)(5:5).ne.'O') then
            if((ntrans.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
               write(5,'(i10,1p,3(1x,e13.6))') node,
     &              (fn(j,node),j=1,3)         
            elseif(inotr(1,node).eq.0) then
               write(5,'(i10,1p,3(1x,e13.6))') node,
     &              (fn(j,node),j=1,3) 
            else
               call transformatrix(trab(1,inotr(1,node)),co(1,node),a)
               write(5,'(i10,1p,3(1x,e13.6),1x,a1)') node,
     &            fn(1,node)*a(1,1)+fn(2,node)*a(2,1)+fn(3,node)*a(3,1),
     &            fn(1,node)*a(1,2)+fn(2,node)*a(2,2)+fn(3,node)*a(3,2),
     &            fn(1,node)*a(1,3)+fn(2,node)*a(2,3)+fn(3,node)*a(3,3),
     &            local
            endif
         endif
      elseif(prlab(ii)(1:4).eq.'RFL ') then
         rftot(0)=rftot(0)+fn(0,node)
         if(prlab(ii)(5:5).ne.'O') then
            write(5,'(i10,1p,3(1x,e13.6))') node,
     &           fn(0,node)
         endif
      endif
!
      flush(5)
!
      return
      end






