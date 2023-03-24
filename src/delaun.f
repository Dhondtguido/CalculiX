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
!     S.W. Sloan, Adv.Eng.Software,1987,9(1),34-55.
!     Permission for use with the GPL license granted by Prof. Scott
!     Sloan on 17. Nov. 2013
!
      subroutine delaun(numpts,n,x,y,list,stack,v,e,numtri)
!     
      implicit none
!     
      integer v(3,*),n,i,t,list(*),numtri,p,e(3,*),maxstk,topstk,
     &     v1,v2,v3,l,r,pop,a,b,c,erl,era,erb,edg,triloc,numpts,
     &     tstrt,tstop,stack(*)
!     
      real*8 x(*),y(*),xp,yp,c00000,c00100
!     
      logical swap
!     
      parameter(c00000=0.d0,
     &     c00100=100000.d0)
!     
      v1=numpts+1
      v2=numpts+2
      v3=numpts+3
      v(1,1)=v1
      v(2,1)=v2
      v(3,1)=v3
      e(1,1)=0
      e(2,1)=0
      e(3,1)=0
!     
      x(v1)=-c00100
      x(v2)= c00100
      x(v3)= c00000
      y(v1)=-c00100
      y(v2)=-c00100
      y(v3)= c00100
!     
      numtri=1
      topstk=0
      maxstk=numpts
      do 100 i=1,n
         p=list(i)
         xp=x(p)
         yp=y(p)
         t=triloc(xp,yp,x,y,v,e,numtri)
         a=e(1,t)
         b=e(2,t)
         c=e(3,t)
         v1=v(1,t)
         v2=v(2,t)
         v3=v(3,t)
         v(1,t)=p
         v(2,t)=v1
         v(3,t)=v2
         e(1,t)=numtri+2
         e(2,t)=a
         e(3,t)=numtri+1
!     
         numtri=numtri+1
         v(1,numtri)=p
         v(2,numtri)=v2
         v(3,numtri)=v3
         e(1,numtri)=t
         e(2,numtri)=b
         e(3,numtri)=numtri+1
         numtri=numtri+1
         v(1,numtri)=p
         v(2,numtri)=v3
         v(3,numtri)=v1
         e(1,numtri)=numtri-1
         e(2,numtri)=c
         e(3,numtri)=t
!     
         if(a.ne.0) then
            call push(t,maxstk,topstk,stack)
         end if
         if(b.ne.0) then
            e(edg(b,t,e),b)=numtri-1
            call push(numtri-1,maxstk,topstk,stack)
         end if
         if(c.ne.0) then
            e(edg(c,t,e),c)=numtri
            call push(numtri,maxstk,topstk,stack)
         end if
!     
 50      if(topstk.gt.0) then
            l=pop(topstk,stack)
            r=e(2,l)
!     
            erl=edg(r,l,e)
            era=mod(erl,3)+1
            erb=mod(era,3)+1
            v1=v(erl,r)
            v2=v(era,r)
            v3=v(erb,r)
            if(swap(x(v1),y(v1),x(v2),y(v2),x(v3),y(v3),xp,yp)) then
               a=e(era,r)
               b=e(erb,r)
               c=e(3,l)
               v(3,l)=v3
               e(2,l)=a
               e(3,l)=r
               v(1,r)=p
               v(2,r)=v3
               v(3,r)=v1
               e(1,r)=l
               e(2,r)=b
               e(3,r)=c
               if(a.ne.0) then
                  e(edg(a,r,e),a)=l
                  call push(l,maxstk,topstk,stack)
               end if
               if(b.ne.0) then
                  call push(r,maxstk,topstk,stack)
               end if
               if(c.ne.0) then
                  e(edg(c,l,e),c)=r
               end if
            end if
            goto 50
         end if
 100  continue
      if(numtri.ne.2*n+1) then
         write(6,'("o***error in subroutine delaun***")')
         write(6,'(" ***incorrect number of triangls formed***")')
         call exit(201)
      end if
      do 120 t=1,numtri
         if((v(1,t).gt.numpts).or.
     &        (v(2,t).gt.numpts).or.
     &        (v(3,t).gt.numpts))then
            do 110 i=1,3
               a=e(i,t)
               if(a.ne.0) then
                  e(edg(a,t,e),a)=0
               end if
 110        continue
            goto 125
         end if
 120  continue
 125  tstrt=t+1
      tstop=numtri
      numtri=t-1
      do 200 t=tstrt,tstop
         if((v(1,t).gt.numpts).or.
     &        (v(2,t).gt.numpts).or.
     &        (v(3,t).gt.numpts))then
            do 130 i=1,3
               a=e(i,t)
               if(a.ne.0) then
                  e(edg(a,t,e),a)=0
               end if
 130        continue
         else
            numtri=numtri+1
            do 140 i=1,3
               a=e(i,t)
               e(i,numtri)=a
               v(i,numtri)=v(i,t)
               if(a.ne.0) then
                  e(edg(a,t,e),a)=numtri
               end if
 140        continue
         endif
 200  continue
      end
      
