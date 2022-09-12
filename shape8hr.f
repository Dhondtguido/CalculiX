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
      subroutine shape8hr(xl,xsj,shp,gs,a)
!
!     shape functions and derivatives for a 8-node linear isoparametric
!     mean strain solid element with hourglass control
!
!     Reference: Flanagan, D.P., Belytschko, T.; "Uniform  strain hexahedron
!     and quadrilateral with orthogonal Hourglass control". Int. J. Num.
!     Meth. Engg., Vol. 17, 679-706, 1981. 
!
!     author: Otto-Ernst Bernhardi
!
      implicit none
!
      integer i,j,k,l
      real*8 shp(4,20),xl(3,20),xsj,vol
      real*8 x1,x2,x3,x4,x5,x6,x7,x8
      real*8 y1,y2,y3,y4,y5,y6,y7,y8
      real*8 z1,z2,z3,z4,z5,z6,z7,z8
      real*8 gb(8,4),gs(8,4),s0,a
!
!
!
      gb = reshape((
     &      / 1.0d0, 1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0, 1.0d0, 1.0d0,
     &        1.0d0,-1.0d0,-1.0d0, 1.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0,
     &        1.0d0,-1.0d0, 1.0d0,-1.0d0, 1.0d0,-1.0d0, 1.0d0,-1.0d0,
     &       -1.0d0, 1.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0, 1.0d0,-1.0d0 /),
     &        (/8,4/))
!
!     shape functions and their global derivatives
!
!     shape functions
!
      shp(4, 1)=1.0d0/8.0d0
      shp(4, 2)=1.0d0/8.0d0
      shp(4, 3)=1.0d0/8.0d0
      shp(4, 4)=1.0d0/8.0d0
      shp(4, 5)=1.0d0/8.0d0
      shp(4, 6)=1.0d0/8.0d0
      shp(4, 7)=1.0d0/8.0d0
      shp(4, 8)=1.0d0/8.0d0
!
!     extract node point coordinates of element
      x1=xl(1,1)
      x2=xl(1,2)
      x3=xl(1,3)
      x4=xl(1,4)
      x5=xl(1,5)
      x6=xl(1,6)
      x7=xl(1,7)
      x8=xl(1,8)
      y1=xl(2,1)
      y2=xl(2,2)
      y3=xl(2,3)
      y4=xl(2,4)
      y5=xl(2,5)
      y6=xl(2,6)
      y7=xl(2,7)
      y8=xl(2,8)
      z1=xl(3,1)
      z2=xl(3,2)
      z3=xl(3,3)
      z4=xl(3,4)
      z5=xl(3,5)
      z6=xl(3,6)
      z7=xl(3,7)
      z8=xl(3,8)
!
!     Average displacement derivative operator matrix, 
!     using eqn 16 in the reference above. 
!     Generated using maxima/wxmaxima and the following input lines. 
!     Note that shp array must be divided by the element volume.
!     h1:(1-r)*(1-s)*(1-t)/8;
!     h2:(1+r)*(1-s)*(1-t)/8;
!     h3:(1+r)*(1+s)*(1-t)/8;
!     h4:(1-r)*(1+s)*(1-t)/8;
!     h5:(1-r)*(1-s)*(1+t)/8;
!     h6:(1+r)*(1-s)*(1+t)/8;
!     h7:(1+r)*(1+s)*(1+t)/8;
!     h8:(1-r)*(1+s)*(1+t)/8;
!     H:matrix([h1,h2,h3,h4,h5,h6,h7,h8])$
!     Bx:diff(H,r);
!     By:diff(H,s);
!     Bz:diff(H,t);
!     B:matrix([Bx[1,1],Bx[1,2],Bx[1,3],Bx[1,4],Bx[1,5],Bx[1,6],Bx[1,7],Bx[1,8]],
!              [By[1,1],By[1,2],By[1,3],By[1,4],By[1,5],By[1,6],By[1,7],By[1,8]],
!              [Bz[1,1],Bz[1,2],Bz[1,3],Bz[1,4],Bz[1,5],Bz[1,6],Bz[1,7],Bz[1,8]]);
!     x:matrix([x1,x2,x3,x4,x5,x6,x7,x8],
!              [y1,y2,y3,y4,y5,y6,y7,y8],
!              [z1,z2,z3,z4,z5,z6,z7,z8]);
!     xr:B.transpose(x);
!     det:determinant(xr);
!     vol:factor(integrate(integrate(integrate(det,r,-1,1),s,-1,1),t,-1,1));
!     shp: matrix( [0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0])$
!     shp[1][1]:ratsimp(diff(vol,x1));
!     shp[1][2]:ratsimp(diff(vol,x2));
!     shp[1][3]:ratsimp(diff(vol,x3));
!     ... (more lines)
!     shp[3][8]:ratsimp(diff(vol,z8));
!     fortran(shp);
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1) = ((y5-y4)*z8+(y2-y5)*z6+(-y8+y6-y4+y2)*z5+(y8+y5-y3-y2)
     1   *z4+(y4-y2)*z3+(-y6-y5+y4+y3)*z2)/1.2d+1
      shp(1,2) = -((y6-y3)*z7+(-y7+y5-y3+y1)*z6+(y1-y6)*z5+(y3-y1)*z4+(
     1   y7+y6-y4-y1)*z3+(-y6-y5+y4+y3)*z1)/1.2d+1
      shp(1,3) = -((y7-y4)*z8+(-y8+y6-y4+y2)*z7+(y2-y7)*z6+(y8+y7-y2-y1
     1   )*z4+(-y7-y6+y4+y1)*z2+(y4-y2)*z1)/1.2d+1
      shp(1,4) = -((y7-y5+y3-y1)*z8+(y3-y8)*z7+(y8-y1)*z5+(-y8-y7+y2+y1
     1   )*z3+(y1-y3)*z2+(y8+y5-y3-y2)*z1)/1.2d+1
      shp(1,5) = ((y7+y6-y4-y1)*z8+(y6-y8)*z7+(-y8-y7+y2+y1)*z6+(y8-y1)
     1   *z4+(y1-y6)*z2+(y8-y6+y4-y2)*z1)/1.2d+1
      shp(1,6) = ((y7-y5)*z8+(-y8-y5+y3+y2)*z7+(y8+y7-y2-y1)*z5+(y2-y7)
     1   *z3+(-y7+y5-y3+y1)*z2+(y5-y2)*z1)/1.2d+1
      shp(1,7) = -((y6+y5-y4-y3)*z8+(-y8-y5+y3+y2)*z6+(y6-y8)*z5+(y8-y3
     1   )*z4+(y8-y6+y4-y2)*z3+(y3-y6)*z2)/1.2d+1
      shp(1,8) = ((y6+y5-y4-y3)*z7+(y5-y7)*z6+(-y7-y6+y4+y1)*z5+(y7-y5+
     1   y3-y1)*z4+(y7-y4)*z3+(y4-y5)*z1)/1.2d+1
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1) = -((x5-x4)*z8+(x2-x5)*z6+(-x8+x6-x4+x2)*z5+(x8+x5-x3-x2
     1   )*z4+(x4-x2)*z3+(-x6-x5+x4+x3)*z2)/1.2d+1
      shp(2,2) = ((x6-x3)*z7+(-x7+x5-x3+x1)*z6+(x1-x6)*z5+(x3-x1)*z4+(x
     1   7+x6-x4-x1)*z3+(-x6-x5+x4+x3)*z1)/1.2d+1
      shp(2,3) = ((x7-x4)*z8+(-x8+x6-x4+x2)*z7+(x2-x7)*z6+(x8+x7-x2-x1)
     1   *z4+(-x7-x6+x4+x1)*z2+(x4-x2)*z1)/1.2d+1
      shp(2,4) = ((x7-x5+x3-x1)*z8+(x3-x8)*z7+(x8-x1)*z5+(-x8-x7+x2+x1)
     1   *z3+(x1-x3)*z2+(x8+x5-x3-x2)*z1)/1.2d+1
      shp(2,5) = -((x7+x6-x4-x1)*z8+(x6-x8)*z7+(-x8-x7+x2+x1)*z6+(x8-x1
     1   )*z4+(x1-x6)*z2+(x8-x6+x4-x2)*z1)/1.2d+1
      shp(2,6) = -((x7-x5)*z8+(-x8-x5+x3+x2)*z7+(x8+x7-x2-x1)*z5+(x2-x7
     1   )*z3+(-x7+x5-x3+x1)*z2+(x5-x2)*z1)/1.2d+1
      shp(2,7) = ((x6+x5-x4-x3)*z8+(-x8-x5+x3+x2)*z6+(x6-x8)*z5+(x8-x3)
     1   *z4+(x8-x6+x4-x2)*z3+(x3-x6)*z2)/1.2d+1
      shp(2,8) = -((x6+x5-x4-x3)*z7+(x5-x7)*z6+(-x7-x6+x4+x1)*z5+(x7-x5
     1   +x3-x1)*z4+(x7-x4)*z3+(x4-x5)*z1)/1.2d+1
!
!     local derivatives of the shape functions: zeta-derivative
!
      shp(3,1) = ((x5-x4)*y8+(x2-x5)*y6+(-x8+x6-x4+x2)*y5+(x8+x5-x3-x2)
     1   *y4+(x4-x2)*y3+(-x6-x5+x4+x3)*y2)/1.2d+1
      shp(3,2) = -((x6-x3)*y7+(-x7+x5-x3+x1)*y6+(x1-x6)*y5+(x3-x1)*y4+(
     1   x7+x6-x4-x1)*y3+(-x6-x5+x4+x3)*y1)/1.2d+1
      shp(3,3) = -((x7-x4)*y8+(-x8+x6-x4+x2)*y7+(x2-x7)*y6+(x8+x7-x2-x1
     1   )*y4+(-x7-x6+x4+x1)*y2+(x4-x2)*y1)/1.2d+1
      shp(3,4) = -((x7-x5+x3-x1)*y8+(x3-x8)*y7+(x8-x1)*y5+(-x8-x7+x2+x1
     1   )*y3+(x1-x3)*y2+(x8+x5-x3-x2)*y1)/1.2d+1
      shp(3,5) = ((x7+x6-x4-x1)*y8+(x6-x8)*y7+(-x8-x7+x2+x1)*y6+(x8-x1)
     1   *y4+(x1-x6)*y2+(x8-x6+x4-x2)*y1)/1.2d+1
      shp(3,6) = ((x7-x5)*y8+(-x8-x5+x3+x2)*y7+(x8+x7-x2-x1)*y5+(x2-x7)
     1   *y3+(-x7+x5-x3+x1)*y2+(x5-x2)*y1)/1.2d+1
      shp(3,7) = -((x6+x5-x4-x3)*y8+(-x8-x5+x3+x2)*y6+(x6-x8)*y5+(x8-x3
     1   )*y4+(x8-x6+x4-x2)*y3+(x3-x6)*y2)/1.2d+1
      shp(3,8) = ((x6+x5-x4-x3)*y7+(x5-x7)*y6+(-x7-x6+x4+x1)*y5+(x7-x5+
     1   x3-x1)*y4+(x7-x4)*y3+(x4-x5)*y1)/1.2d+1
!
!     computation of element volume (eqn 15)
!
      vol=0.0d0
      do k=1,8
        vol=vol+xl(1,k)*shp(1,k)
      enddo
!
!     computation of the average jacobian determinant
!
      xsj=vol/8.0d0
!
!     hourglass control vectors(see appendix 2 from the reference above). 
!     divide shp array by element volume
      a=0.0d0
      do i=1,8
         do j=1,3
            a=a+shp(j,i)*shp(j,i)
            s0=shp(j,i)/vol
            shp(j,i)=s0
            do k=1,4
               gs(i,k)=gb(i,k)
               do l=1,8
                  gs(i,k)=gs(i,k)-s0*xl(j,l)*gb(l,k)
               enddo
            enddo
         enddo
      enddo
c
c     calculate hourglass control stiffness factor a 
c     (to be used in hgstiffness() and hgforce())
c
c     in ABAQUS, a 0.005 is used as default value.
      a=0.0005d0*a/vol
c
      return
      end
