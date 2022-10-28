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
      subroutine beamintscheme(lakonl,mint3d,npropstart,prop,
     &  kk,xi,et,ze,weight)
!
!     provides the integration scheme for beams with a cross section
!     which is not rectangular nor elliptical
!
!     mint3d: number of integration points (returned if kk=0)
!     xi,et,ze: local coordinates of integration point kk
!     weight: weight for integration point kk
!
      implicit none
!
      character*8 lakonl
!
      integer mint3d,npropstart,jj,kk
!
      real*8 prop(*),xi,et,ze,weight,ratio,ratio2,dtheta,theta,r,
     &  t1,t2,t3,t4,a,b
!
!
!
      if(lakonl(8:8).eq.'P') then
!
!        pipe cross section
!
         if(kk.eq.0) then
            mint3d=16
            return
         endif
!
!        ratio of inner radius to outer radius
!
         ratio=(prop(npropstart+1)-prop(npropstart+2))/
     &              prop(npropstart+1)
         ratio2=ratio*ratio
!
         if(kk.gt.8) then
            jj=kk-8
            xi=1.d0/dsqrt(3.d0)
         else
            jj=kk
            xi=-1.d0/dsqrt(3.d0)
         endif
!
!        pi/4
!
         dtheta=datan(1.d0)
!
         theta=(jj-1)*dtheta
         r=dsqrt((ratio2+1.d0)/2.d0)
!
         et=r*dcos(theta)
         ze=r*dsin(theta)
         weight=dtheta*(1.d0-ratio2)/2.d0
c
c     Box cross section
      elseif(lakonl(8:8).eq.'B') then
         if(kk.eq.0) then
           mint3d=32
           return
         endif
!
!        2 pts in long direction xi
!
         if(kk.gt.16) then
            jj=kk-16
            xi=1.d0/dsqrt(3.d0)
         else
            jj=kk
            xi=-1.d0/dsqrt(3.d0)
         endif
!
!        pts in cross sections
!
         a=prop(npropstart+1)
         b=prop(npropstart+2)
         t1=prop(npropstart+3)
         t2=prop(npropstart+4)
         t3=prop(npropstart+5)
         t4=prop(npropstart+6)
!
         if(jj.eq.1)then 
            et = -(t4-b)/b
            ze = -(t1-a)/a
            weight = -((((-2*a)+2*t1+t3)*t4+t1*t2-2*b*t1)/(a*b))/6.0d+0
         elseif(jj.eq.2)then 
            et = -((3*t4-t2-2*b)/b)/4.0d+0
            ze = -(t1-a)/a
            weight = -((2*t1*t4+2*t1*t2-4*b*t1)/(a*b))/3.0d+0
         elseif(jj.eq.3)then 
            et = -((t4-t2)/b)/2.0d+0
            ze = -(t1-a)/a
            weight = -((t1*t4+t1*t2-2*b*t1)/(a*b))/3.0d+0
         elseif(jj.eq.4)then 
            et = -((t4-3*t2+2*b)/b)/4.0d+0
            ze = -(t1-a)/a
            weight = -((2*t1*t4+2*t1*t2-4*b*t1)/(a*b))/3.0d+0
         elseif(jj.eq.5)then 
            et = (t2-b)/b
            ze = -(t1-a)/a
            weight = -((t1*t4+t2*t3+(2*t1-2*a)*t2-2*b*t1)/(a*b))/6.0d+0
         elseif(jj.eq.6)then 
            et = (t2-b)/b
            ze = ((t3-3*t1+2*a)/a)/4.0d+0
            weight = -((2*t2*t3+(2*t1-4*a)*t2)/(a*b))/3.0d+0
         elseif(jj.eq.7)then 
            et = (t2-b)/b
            ze = ((t3-t1)/a)/2.0d+0
            weight = -((t2*t3+(t1-2*a)*t2)/(a*b))/3.0d+0
         elseif(jj.eq.8)then 
            et = (t2-b)/b
            ze = ((3*t3-t1-2*a)/a)/4.0d+0
            weight = -((2*t2*t3+(2*t1-4*a)*t2)/(a*b))/3.0d+0
         elseif(jj.eq.9)then 
            et = (t2-b)/b
            ze = (t3-a)/a
            weight = -((t3*t4+(2*t2-2*b)*t3+(t1-2*a)*t2)/(a*b))/6.0d+0
         elseif(jj.eq.10)then 
            et = -((t4-3*t2+2*b)/b)/4.0d+0
            ze = (t3-a)/a
            weight = -((2*t3*t4+(2*t2-4*b)*t3)/(a*b))/3.0d+0
         elseif(jj.eq.11)then 
            et = -((t4-t2)/b)/2.0d+0
            ze = (t3-a)/a
            weight = -((t3*t4+(t2-2*b)*t3)/(a*b))/3.0d+0
         elseif(jj.eq.12)then 
            et = -((3*t4-t2-2*b)/b)/4.0d+0
            ze = (t3-a)/a
            weight = -((2*t3*t4+(2*t2-4*b)*t3)/(a*b))/3.0d+0
         elseif(jj.eq.13)then 
            et = -(t4-b)/b
            ze = (t3-a)/a
            weight = -((((-2*a)+t1+2*t3)*t4+(t2-2*b)*t3)/(a*b))/6.0d+0
         elseif(jj.eq.14)then 
            et = -(t4-b)/b
            ze = ((3*t3-t1-2*a)/a)/4.0d+0
            weight = -(((2*t3+2*t1-4*a)*t4)/(a*b))/3.0d+0
         elseif(jj.eq.15)then 
            et = -(t4-b)/b
            ze = ((t3-t1)/a)/2.0d+0
            weight = -(((t3+t1-2*a)*t4)/(a*b))/3.0d+0
         elseif(jj.eq.16)then 
            et = -(t4-b)/b
            ze = ((t3-3*t1+2*a)/a)/4.0d+0
            weight = -(((2*t3+2*t1-4*a)*t4)/(a*b))/3.0d+0
         endif 
      endif
!     
      return
      end
