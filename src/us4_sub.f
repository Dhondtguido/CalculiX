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
      SUBROUTINE us4_csys(xg,tm,tmg)
!
      IMPLICIT NONE                      
!
!     determines the material data for element iel
!     4-node shell element
!     author: Gil Rama
!
      REAL*8, INTENT(IN)  :: xg(4,3)         ! element coordinates in global csys
      REAL*8, INTENT(OUT) :: tm(3,3),tmg(24,24)         ! transformation matrix (e1,e2,e3)T
      REAL*8 :: e1(3),e2(3),e3(3),dl,dd,xno(3)
      REAL*8 :: xi,et,xl(3,8),xs(3,7),p(3),shp(7,8)
      REAL*8 :: a(3,3)    
      INTEGER :: iflag,j,i,l,k
      !
      tm(:,:)  = 0.d0
      tmg(:,:) = 0.d0
      !
      do i=1,3
        do j=1,4
      	 xl(i,j) = xg(j,i)
        enddo
      enddo
      xi = 0.d0
      et = 0.d0
      iflag = 2
      call shape4q(xi,et,xl,xno,xs,shp,iflag)
      dd = dsqrt(xno(1)*xno(1)+xno(2)*xno(2)+xno(3)*xno(3))
      do l = 1,3
        xno(l)=xno(l)/dd
      enddo
      ! coordinates of the point at (xi,et)=(0.,0.)      
      do l = 1,3
        p(l) = 0.d0
        do k = 1,4
          p(l) = p(l) + shp(4,k)*xl(l,k)
        enddo
      enddo      
      ! unit matrix
      do k = 1,3
        do l = 1,3
          a(k,l) = 0.d0
        enddo
        a(k,k) = 1.d0
      enddo
      dd = a(1,1)*xno(1)+a(2,1)*xno(2)+a(3,1)*xno(3)
      if(dabs(dd).gt.0.999999999536d0) then
        ! project the z-axis
        dd = a(1,3)*xno(1)+a(2,3)*xno(2)+a(3,3)*xno(3)
         do l = 1,3
           e1(l) = a(l,3)-dd*xno(l)
         enddo
      else
        ! project the x-axis
        do l =1 ,3
          e1(l) = a(l,1)-dd*xno(l)
        enddo
      endif
      !
      dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do l = 1,3
        e1(l) = e1(l)/dd
      enddo
      ! e2 = n x e1
      e2(1) = xno(2)*e1(3) - e1(2)*xno(3)
      e2(2) = xno(3)*e1(1) - e1(3)*xno(1)
      e2(3) = xno(1)*e1(2) - e1(1)*xno(2)
      do l = 1,3
        e3(l) = xno(l)
      enddo
      ! put in tm and out
      do j = 1,3
         tm(1,j) = e1(j)
         tm(2,j) = e2(j)
         tm(3,j) = e3(j)
      enddo 
      !
      do i = 1,3
        do j = 1,3
          tmg(i,j)       = tm(i,j) ! f1
          tmg(i+3,j+3)   = tm(i,j) ! r1
          tmg(i+6,j+6)   = tm(i,j) ! f2
          tmg(i+9,j+9)   = tm(i,j) ! r2
          tmg(i+12,j+12) = tm(i,j) ! f3
          tmg(i+15,j+15) = tm(i,j) ! r3
          tmg(i+18,j+18) = tm(i,j) ! f4
          tmg(i+21,j+21) = tm(i,j) ! r4          
        enddo
      enddo      
      RETURN
      END
      !
      SUBROUTINE us4_Km(X,Dm,Ki)      
      REAL*8, INTENT(IN)  :: X(4,3),Dm(3,3)           
      REAL*8, INTENT(OUT) :: Ki(24,24)
      !
      REAL*8 :: Nrs(4),dNr(4),dNs(4),Jm(2,2),invJm(2,2)
      REAL*8 :: ri,si,bb(3,24),g_p(4,3),detJm,detinvJm  
      REAL*8 :: dNx(4),dNy(4)
      !
      INTEGER :: i,k
      !
      ! Gauß points(4):
      g_p(1,1) = -0.577350269189626
      g_p(2,1) = +0.577350269189626       
      g_p(3,1) = +0.577350269189626 
      g_p(4,1) = -0.577350269189626  
      !   
      g_p(1,2) = -0.577350269189626 
      g_p(2,2) = -0.577350269189626       
      g_p(3,2) = +0.577350269189626 
      g_p(4,2) = +0.577350269189626
      !    
      g_p(1,3) = +1.000000000000000 
      g_p(2,3) = +1.000000000000000       
      g_p(3,3) = +1.000000000000000 
      g_p(4,3) = +1.000000000000000 
      !
      Ki(:,:) = 0.d0
      do k = 1,4    
        ri = g_p(k,1)
        si = g_p(k,2)               
      	call us4_Ni(ri,si,X,Nrs,dNr,dNs,Jm,invJm,detJm,detinvJm,dNx,dNy) 
      	call us4_Bmi(dNx,dNy,bb)
      	Ki=Ki+matmul(matmul(transpose(bb),Dm),bb)*detJm
      enddo 
      !
      RETURN
      END
      !
      SUBROUTINE us4_Kb(X,Db,Ki)      
      REAL*8, INTENT(IN)  :: X(4,3),Db(3,3)           
      REAL*8, INTENT(OUT) :: Ki(24,24)
      !
      REAL*8 :: Nrs(4),dNr(4),dNs(4),Jm(2,2),invJm(2,2)
      REAL*8 :: ri,si,bb(3,24),g_p(4,3),detJm,detinvJm  
      REAL*8 :: dNx(4),dNy(4)
      !
      INTEGER :: i,k
      !
      ! Gauß points(4):
      g_p(1,1) = -0.577350269189626
      g_p(2,1) = +0.577350269189626       
      g_p(3,1) = +0.577350269189626 
      g_p(4,1) = -0.577350269189626  
      !   
      g_p(1,2) = -0.577350269189626 
      g_p(2,2) = -0.577350269189626       
      g_p(3,2) = +0.577350269189626 
      g_p(4,2) = +0.577350269189626
      !
      Ki(:,:) = 0.d0
      do k = 1,4    
        ri = g_p(k,1)
        si = g_p(k,2)               
      	call us4_Ni(ri,si,X,Nrs,dNr,dNs,Jm,invJm,detJm,detinvJm,dNx,dNy) 
      	call us4_Bbi(dNx,dNy,bb)
      	Ki=Ki+matmul(matmul(transpose(bb),Db),bb)*detJm!*g_p(k,3)
      enddo 
      !
      RETURN
      END
      !      
      SUBROUTINE us4_Ks(X,Ds,Ki)      
      REAL*8, INTENT(IN)  :: X(4,3),Ds(2,2)           
      REAL*8, INTENT(OUT) :: Ki(24,24)
      !
      REAL*8 :: Nrs(4),dNr(4),dNs(4),Jm(2,2),invJm(2,2)
      REAL*8 :: ri,si,bs(2,24),g_p(4,3),detJm,detinvJm  
      REAL*8 :: dNx(4),dNy(4)
      !
      INTEGER :: i,k
      !
      ! Gauß points(4):
      g_p(1,1) = -0.577350269189626
      g_p(2,1) = +0.577350269189626       
      g_p(3,1) = +0.577350269189626 
      g_p(4,1) = -0.577350269189626  
      !   
      g_p(1,2) = -0.577350269189626 
      g_p(2,2) = -0.577350269189626       
      g_p(3,2) = +0.577350269189626 
      g_p(4,2) = +0.577350269189626
      !    
      g_p(1,3) = +1.000000000000000 
      g_p(2,3) = +1.000000000000000       
      g_p(3,3) = +1.000000000000000 
      g_p(4,3) = +1.000000000000000 
      !
      Ki(:,:) = 0.d0
      !
      do k = 1,4    
        ri = g_p(k,1)
        si = g_p(k,2)               
      	call us4_Ni(ri,si,X,Nrs,dNr,dNs,Jm,invJm,detJm,detinvJm,dNx,dNy) 
      	call us4_Bsi_ANS(ri,si,X,bs)  ! shear strain displacement matrix
      	Ki=Ki+matmul(matmul(transpose(bs),Ds),bs)*detJm*g_p(k,3)
      enddo 
      !
      !call us4_Ni(0.0,0.0,X,Nrs,dNr,dNs,Jm,invJm,detJm,detinvJm,dNx,dNy)
      !call us4_Bsi(dNx,dNy,bs,Nrs)
      !Ki = matmul(matmul(transpose(bs),Ds),bs)*detJm*2.d0
      !
      RETURN
      END      
      ! Ni  element interpolation
      SUBROUTINE us4_Ni(r,s,x,N,dNr,dNs,Jm,invJm,detJm,detinvJm,dx,dy)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: r         ! natural coordinate xi
      REAL*8, INTENT(IN)  :: s         ! natural coordinate eta
      REAL*8, INTENT(OUT) :: N(4)      ! array shape functions (4) 
      REAL*8, INTENT(OUT) :: dNr(4)    ! array dN/dxi  (4)
      REAL*8, INTENT(OUT) :: dNs(4)    ! array dN/deta (4) 
      REAL*8, INTENT(IN)  :: x(4,3)    ! nodal coordinates in element frame
      REAL*8, INTENT(OUT) :: Jm(2,2)   ! jacobian matrix
      REAL*8, INTENT(OUT) :: invJm(2,2)! inverse jacobian matrix   
      REAL*8, INTENT(OUT) :: detJm     !      
      REAL*8, INTENT(OUT) :: detinvJm  !
      REAL*8, INTENT(OUT) :: dx(4)     ! array dN/dx'  (4)
      REAL*8, INTENT(OUT) :: dy(4)     ! array dN/dy'  (4) 
      INTEGER :: i,k
      !
      ! shape functions (bi-linear)
      N(1) = 0.25d0*(1.d0-r)*(1.d0-s)
      N(2) = 0.25d0*(1.d0+r)*(1.d0-s)
      N(3) = 0.25d0*(1.d0+r)*(1.d0+s)            
      N(4) = 0.25d0*(1.d0-r)*(1.d0+s)
      ! dN/dxi:
      dNr(1) = -0.25d0*(1.d0-s)
      dNr(2) = +0.25d0*(1.d0-s)
      dNr(3) = +0.25d0*(1.d0+s)            
      dNr(4) = -0.25d0*(1.d0+s)
      ! dN/deta:      
      dNs(1) = -0.25d0*(1.d0-r)
      dNs(2) = -0.25d0*(1.d0+r)
      dNs(3) = +0.25d0*(1.d0+r)            
      dNs(4) = +0.25d0*(1.d0-r)    
      !  J-matrix
      do i=1,2
         do k=1,2
            Jm(i,k) = 0.d0
         enddo
      enddo
      !
      ! 2D J-matrix
      do i=1,4
         Jm(1,1) = Jm(1,1) + dNr(i)*x(i,1)
         Jm(1,2) = Jm(1,2) + dNr(i)*x(i,2)
         Jm(2,1) = Jm(2,1) + dNs(i)*x(i,1)
         Jm(2,2) = Jm(2,2) + dNs(i)*x(i,2)       
      enddo
      ! inverse J-matrix (direct):
      ! Calculate the inverse determinant of the matrix
      detJm    = (Jm(1,1)*Jm(2,2) - Jm(1,2)*Jm(2,1))
      detinvJm = 1/(detJm)
      !
      ! Calculate the inverse of the matrix
      invJm(1,1) = +detinvJm * Jm(2,2)
      invJm(2,1) = -detinvJm * Jm(2,1)
      invJm(1,2) = -detinvJm * Jm(1,2)
      invJm(2,2) = +detinvJm * Jm(1,1)
      !
      ! shape function derivatives dx', dy':    
      do i=1,4
         dx(i) = invJm(1,1)*dNr(i) + invJm(1,2)*dNs(i)
         dy(i) = invJm(2,1)*dNr(i) + invJm(2,2)*dNs(i)
      enddo
      !       
      RETURN
      END
      !
      ! Bmi computes membrane strain-displacement-matrix
      SUBROUTINE us4_Bmi(dNx,dNy,bm)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: dNx(4)    ! array dN/dx'  (4)
      REAL*8, INTENT(IN)  :: dNy(4)    ! array dN/dy'  (4) 
      REAL*8, INTENT(OUT) :: bm(3,24)  ! membrane strain displacement matrix  (3,12)
      INTEGER :: i,k,i1,i2
      !
      do i = 1,3
      	do k = 1,24
         bm(i,k) = 0.d0
        enddo  
      enddo    
      !  
      do i = 1,4
         i1 = i  + 5*(i-1)  
         i2 = i1 + 1
         bm(1,i1) = dNx(i)
         bm(2,i2) = dNy(i)
         bm(3,i1) = dNy(i)
         bm(3,i2) = dNx(i)  
      enddo
      !
      RETURN
      END
      !
      SUBROUTINE us4_Bbi(dNx,dNy,bb)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)  :: dNx(4)    ! array dN/dx'  (4)
      REAL*8, INTENT(IN)  :: dNy(4)    ! array dN/dy'  (4) 
      REAL*8, INTENT(OUT) :: bb(3,24)  
      INTEGER :: i,k,i4,i5
      !
      DO i = 1,3
      	DO k = 1,24
         bb(i,k) = 0.d0
        ENDDO  
      ENDDO    
      ! 
      do i = 1,4
         i4 = i + 3 + 5*(i-1)  
         i5 = i4 + 1
         !
         bb(1,i5) = +dNx(i)
         bb(2,i4) = -dNy(i)
         bb(3,i4) = -dNx(i)
         bb(3,i5) = +dNy(i)  
      enddo
      !
      RETURN
      END  
      ! Bsi computes shear strain-displacement-matrix
      SUBROUTINE us4_Bsi(dNx,dNy,bs,Nrs)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)   :: dNx(4)    ! array dN/dx'  (4)
      REAL*8, INTENT(IN)   :: dNy(4)    ! array dN/dy'  (4) 
      REAL*8, INTENT(OUT)  :: bs(2,24)  ! membrane strain displacement matrix  (3,12)
      REAL*8, INTENT(IN)   :: Nrs(4)    ! array shape functions (4)           
      INTEGER :: i,k,i3,i4,i5
      !
      do i = 1,2
      	do k = 1,24
         bs(i,k) = 0.d0
        enddo  
      enddo    
      !  
      do i = 1,4
         i3 = i + 2  + 5*(i-1)  
         i4 = i3 + 1
         i5 = i4 + 1
         bs(2,i3) = +dNy(i)
         bs(2,i4) = -Nrs(i)
         bs(1,i3) = +dNx(i)
         bs(1,i5) = +Nrs(i)
      enddo
      !
      RETURN
      END      
      !
      ! Bsi_ANS computes shear strain-displacement-matrix
      SUBROUTINE us4_Bsi_ANS(r,s,x,bs_ANS)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)   :: r         !
      REAL*8, INTENT(IN)   :: s         !
      REAL*8, INTENT(IN)   :: x(4,3)    ! 
      REAL*8, INTENT(OUT)  :: bs_ANS(2,24)    ! 
      !
      REAL*8 ::  N(4),dNr(4),dNs(4),Jm(2,2),invJm(2,2),detJm,detinvJm,
     &  dNx(4),dNy(4),bs_A(2,24),bs_B(2,24),bs_C(2,24),bs_D(2,24)
      !
      INTEGER :: i,k
      !
      do i = 1,2
         do k = 1,24
           bs_ANS(i,k) = 0.d0
           bs_A(i,k)   = 0.d0
           bs_B(i,k)   = 0.d0
           bs_C(i,k)   = 0.d0
           bs_D(i,k)   = 0.d0
         enddo
      enddo
      !
      ! assumed natural strains interpolated due to sampling points (A,B,C,D) values
      ! _______ ________
      ! *      A        *
      ! |      ^        |
      ! |               |
      ! B      x        D
      ! |               |
      ! |               |
      ! *_______C_______*   
      ! 
      ! strain displacement matrix @ A (0,+1)
      call us4_Ni(0.d0,+1.d0,x,N,dNr,dNs,Jm,
     &  invJm,detJm,detinvJm,dNx,dNy)
      call us4_Bsi_1(dNx,dNy,bs_A,N)        
      ! strain displacement matrix @ C (0,-1)
      call us4_Ni(0.d0,-1.d0,x,N,dNr,dNs,Jm,invJm,
     &  detJm,detinvJm,dNx,dNy)
      call us4_Bsi_1(dNx,dNy,bs_C,N) 
      ! strain displacement matrix @ B (-1,0)
      call us4_Ni(-1.d0,0.d0,x,N,dNr,dNs,Jm,invJm,
     &  detJm,detinvJm,dNx,dNy)
      call us4_Bsi_0(dNx,dNy,bs_B,N) 
      ! strain displacement matrix @ D (1,0)
      call us4_Ni(+1.d0,0.d0,x,N,dNr,dNs,Jm,invJm,
     &  detJm,detinvJm,dNx,dNy)
      call us4_Bsi_0(dNx,dNy,bs_D,N) 
      !
      ! assumed natural strains:
      bs_ANS = 0.5d0*((1.d0-s)*bs_A+(1.d0+s)*bs_C)
     &  + 0.5d0*((1.d0-r)*bs_B+(1.d0+r)*bs_D)                  
      RETURN
      END 
      ! Bsi computes shear strain-displacement-matrix
      SUBROUTINE us4_Bsi_1(dNx,dNy,bs,Nrs)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)   :: dNx(4)    ! array dN/dx'  (4)
      REAL*8, INTENT(IN)   :: dNy(4)    ! array dN/dy'  (4) 
      REAL*8, INTENT(OUT)  :: bs(2,24)  ! membrane strain displacement matrix  (3,12)
      REAL*8, INTENT(IN)   :: Nrs(4)    ! array shape functions (4)           
      INTEGER :: i,k,i1,i2,i3,i4,i5
      !
      do i = 1,2
      	do k = 1,24
         bs(i,k) = 0.d0
        enddo  
      enddo    
      !  
      do i = 1,4
         i3 = i + 2 + 5*(i-1)  
         i4 = i3 + 1
         i5 = i4 + 1
         bs(1,i3) = +dNx(i)
         bs(1,i5) = +Nrs(i) 
      enddo
      !
      RETURN
      END
      ! -------------------------------------------------------------
      ! Bsi computes shear strain-displacement-matrix
      SUBROUTINE us4_Bsi_0(dNx,dNy,bs,Nrs)
      IMPLICIT NONE                      
      REAL*8, INTENT(IN)   :: dNx(4)    ! array dN/dx'  (4)
      REAL*8, INTENT(IN)   :: dNy(4)    ! array dN/dy'  (4) 
      REAL*8, INTENT(OUT)  :: bs(2,24)  ! membrane strain displacement matrix  (3,12)
      REAL*8, INTENT(IN)   :: Nrs(4)    ! array shape functions (4)           
      INTEGER :: i,k,i1,i2,i3,i4,i5
      !
      do i = 1,2
      	do k = 1,24
         bs(i,k) = 0.d0
        enddo  
      enddo    
      !  
      do i = 1,4
         i1 = i  + 5*(i-1)  
         i3 = i + 2 + 5*(i-1)  
         i4 = i3 + 1
         i5 = i4 + 1
         bs(2,i3) = +dNy(i)
         bs(2,i4) = -Nrs(i) 
      enddo
      !
      RETURN
      END      
      !
      SUBROUTINE us4_M(X,h,rho,M)  
      REAL*8, INTENT(IN)  :: X(4,3),rho,h             
      REAL*8, INTENT(OUT) :: M(24,24) 
      REAL*8 :: ri,si,Nrs(4),dNr(4),dNs(4),Jm(2,2)
      REAL*8 :: invJm(2,2),detJm,detinvJm,dNx(4),dNy(4),q1
      REAL*8 :: m_3t(6,6), N_u(6,24),g_p(4,3)
      INTEGER :: k,j
      !
      ! Gauß points
      !
      M(:,:)    = 0.d0
      N_u(:,:)  = 0.d0 
      m_3t(:,:) = 0.d0
      !
      g_p(1,1) = -0.577350269189626
      g_p(2,1) = +0.577350269189626       
      g_p(3,1) = +0.577350269189626 
      g_p(4,1) = -0.577350269189626  
      !   
      g_p(1,2) = -0.577350269189626 
      g_p(2,2) = -0.577350269189626       
      g_p(3,2) = +0.577350269189626 
      g_p(4,2) = +0.577350269189626
      !    
      g_p(1,3) = +1.000000000000000 
      g_p(2,3) = +1.000000000000000       
      g_p(3,3) = +1.000000000000000 
      g_p(4,3) = +1.000000000000000 
      !
      q1 = rho*h
      m_3t(1,1) = q1
      m_3t(2,2) = q1
      m_3t(3,3) = q1
      q1 = (rho*h**3)/12.d0
      m_3t(4,4) = q1
      m_3t(5,5) = q1
      m_3t(6,6) = q1
      !
      do k = 1,4    
        ri = g_p(k,1)
        si = g_p(k,2)               
      	call us4_Ni(ri,si,X,Nrs,dNr,dNs,Jm,invJm,detJm,detinvJm,dNx,dNy) ! s4 interpolation
      	do j = 1,6
      	  N_u(j,j)    = Nrs(1)
      	  N_u(j,j+6)  = Nrs(2)
      	  N_u(j,j+12) = Nrs(3)
      	  N_u(j,j+18) = Nrs(4)
      	enddo
      	M=M+matmul(matmul(transpose(N_u),m_3t),N_u)*detJm*g_p(k,3)
      enddo
      !
      RETURN
      END      
      !
      SUBROUTINE us4_xu(x,ushell,ueg,xg,tm,vl)      
      REAL*8, INTENT(IN)  :: xg(4,3),tm(3,3),vl(6,4)            
      REAL*8, INTENT(OUT) :: x(4,3),ushell(24),ueg(24)
      !
      x(:,:) = 0.d0
      x(1,:) = matmul(tm,xg(1,:))
      x(2,:) = matmul(tm,xg(2,:))
      x(3,:) = matmul(tm,xg(3,:))
      x(4,:) = matmul(tm,xg(4,:))
      !
      ushell(1:3)   = matmul(tm,vl(1:3,1))
      ushell(4:6)   = matmul(tm,vl(4:6,1))
      ushell(7:9)   = matmul(tm,vl(1:3,2))
      ushell(10:12) = matmul(tm,vl(4:6,2))
      ushell(13:15) = matmul(tm,vl(1:3,3))
      ushell(16:18) = matmul(tm,vl(4:6,3))
      ushell(19:21) = matmul(tm,vl(1:3,4))
      ushell(22:24) = matmul(tm,vl(4:6,4))      
      !
      ueg(:) = 0.d0
      ueg(1:3)   = vl(1:3,1)
      ueg(4:6)   = vl(4:6,1)
      ueg(7:9)   = vl(1:3,2)
      ueg(10:12) = vl(4:6,2)
      ueg(13:15) = vl(1:3,3)
      ueg(16:18) = vl(4:6,3)
      ueg(19:21) = vl(1:3,3)
      ueg(22:24) = vl(4:6,3)      
      !
      RETURN
      END
