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
      subroutine effectivemodalmass(neq,nactdof,mi,adb,aub,jq,
     &  irow,nev,z,co,nk)
!
!     calculates the effective modal mass for frequency 
!     calculations
!
      implicit none
!
      integer i,j,k,nk,mi(*),nactdof(0:mi(2),*),neq(*),jq(*),irow(*),
     &  nev
!
      real*8 x(neq(2)),y(neq(2)),adb(*),aub(*),z(*),part(nev,6),
     &  toteffmass(6),effmodmass(nev,6),toteffmodmass(6),co(3,*)
!
!     translations in x, y and z
!
      do j=1,3
         do i=1,neq(2)
            x(i)=0.d0
         enddo
         do i=1,nk
            if(nactdof(j,i).gt.0) x(nactdof(j,i))=1.d0
         enddo
!
!        {y}=[M].{x}
!
         call op(neq(2),x,y,adb,aub,jq,irow)
!
!        participation factors {phi_k}^T.{y}
!
         do k=1,nev
            call multvec(neq(2),z((k-1)*neq(2)+1),y,part(k,j))
         enddo
!
!        total effective mass   {x}^T.{y}
!
         call multvec(neq(2),x,y,toteffmass(j))
!
!        effective modal mass and its total
!
         toteffmodmass(j)=0.d0
         do k=1,nev
            effmodmass(k,j)=(part(k,j)**2)
            toteffmodmass(j)=toteffmodmass(j)+effmodmass(k,j)
         enddo
      enddo
!
!     rotations about x, y and z
!
      do j=4,6
         do i=1,neq(2)
            x(i)=0.d0
         enddo
         do i=1,nk
            if(j.eq.4) then
               if(nactdof(2,i).gt.0) then
                  x(nactdof(2,i))=-co(3,i)
               endif
               if(nactdof(3,i).gt.0) then
                  x(nactdof(3,i))=co(2,i)
               endif
            elseif(j.eq.5) then
               if(nactdof(3,i).gt.0) then
                  x(nactdof(3,i))=-co(1,i)
               endif
               if(nactdof(1,i).gt.0) then
                  x(nactdof(1,i))=co(3,i)
               endif
            elseif(j.eq.6) then
               if(nactdof(1,i).gt.0) then
                  x(nactdof(1,i))=-co(2,i)
               endif
               if(nactdof(2,i).gt.0) then
                  x(nactdof(2,i))=co(1,i)
               endif
            endif
         enddo
!
!        {y}=[M].{x}
!
         call op(neq(2),x,y,adb,aub,jq,irow)
!
!        participation factors {phi_k}^T.{y}
!
         do k=1,nev
            call multvec(neq(2),z((k-1)*neq(2)+1),y,part(k,j))
         enddo
!
!        total effective mass   {x}^T.{y}
!
         call multvec(neq(2),x,y,toteffmass(j))
!
!        effective modal mass and its total
!
         toteffmodmass(j)=0.d0
         do k=1,nev
            effmodmass(k,j)=(part(k,j)**2)
            toteffmodmass(j)=toteffmodmass(j)+effmodmass(k,j)
         enddo
      enddo
!
!     writing the participation factors into the .dat-file
!
      write(5,*)
      write(5,*) '    P A R T I C I P A T I O N   F A C T O R S'
      write(5,*)
      write(5,100)
 100  format(   'MODE NO.   X-COMPONENT     Y-COMPONENT     Z-COMPONENT 
     &    X-ROTATION      Y-ROTATION      Z-ROTATION')
      write(5,*)
      do k=1,nev
         write(5,'(i7,6(2x,e14.7))') k,(part(k,j),j=1,6)
      enddo
!
!     writing the effective mass into the .dat-file
!
      write(5,*)
      write(5,*) '    E F F E C T I V E   M O D A L   M A S S'
      write(5,*)
      write(5,100)
      write(5,*)
      do k=1,nev
         write(5,'(i7,6(2x,e14.7))') k,(effmodmass(k,j),j=1,6)
      enddo
      write(5,'(a7,6(2x,e14.7))') 'TOTAL  ',(toteffmodmass(j),j=1,6)
!
!     writing the total effective mass into the .dat-file
!
      write(5,*)
      write(5,*) '    T O T A L   E F F E C T I V E   M A S S'
      write(5,*)
      write(5,100)
      write(5,*)
      write(5,'(a7,6(2x,e14.7))') '       ',(toteffmass(j),j=1,6)
      write(5,*)
!
      return
      end
