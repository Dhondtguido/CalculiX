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
      subroutine umat_abaqusnl_total(amat,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,elas,
     &        iorien,pgauss,orab,istep,kinc,pnewdt,nmethod,iperturb)
!
!     converting nonlinear fields used by CalculiX into nonlinear
!     fields used by Abaqus and vice versa.
!
!     before the umat call: convert the Lagrange strain into
!                           corotational logarithmic strain
!     after the umat call: convert corotational Cauchy stress
!                          into PK2 stress and convert the derivative
!                          of the corotational Cauchy stress w.r.t. the
!                          corotational logarithmic strain into the
!                          derivative of the PK2 stress w.r.t. the
!                          Lagrange strain
! 
      implicit none
!     
      character*80 amat
!     
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),i,
     &     iorien,nmethod,iperturb(*),istep,nprops,jstep(4),kinc,
     &     kel(4,21),j1,j2,j3,j4,j5,j6,j7,j8,jj,n,ier,j,matz,kal(2,6),
     &     keltot(4,36),iflag
!     
      real*8 elconloc(*),elas(21),emec(6),emec0(6),beta(6),stre(6),
     &     vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &     time,ttime,skl(3,3),xa(3,3),ya(3,3,3,3),z(3,3),
     &     xstateini(nstate_,mi(1),*),w(3),fv1(3),fv2(3),d(6),c(6),
     &     v1,v2,v3,eln(6),e(3,3),tkl(3,3),u(6),c2(6),dd,um1(3,3),
     &     expansion,ctot(3,3),ddsdde(6,6),spd,rpl,pnewdt,stran(6),temp,
     &     xstate(nstate_,mi(1),*),dp(21),ciicp(21),xm1(6),xm2(6),
     &     dum1dc(21),factor,al(3),xm1m1(21),xm2m2(21),xm3m3(21),aa(21),
     &     bb(21),d1,d2,d3,f1,f2,f3,g1,g2,g3,h1,h2,h3,
     &     dum1dcglob(3,3,3,3),dsigdc(3,3,3,3),cm1(3,3),elastot(36),
     &     xm3(6),delndc(21),delndcglob(3,3,3,3)
!     
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!     
      keltot=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3,2,2,1,1,3,3,1,1,3,3,2,2,
     &     1,2,1,1,1,2,2,2,1,2,3,3,1,3,1,1,1,3,2,2,1,3,3,3,
     &     1,3,1,2,2,3,1,1,2,3,2,2,2,3,3,3,2,3,1,2,2,3,1,3/),(/4,36/))
!     
      d=(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
      dp=(/1.d0,0.d0,1.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.5d0,0.d0,
     &     0.d0,0.d0,0.d0,0.5d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.5d0/)
!     
!     calculating the logarithmic mechanical strain at the
!     end of the increment
!     
      e(1,1)=emec(1)
      e(2,2)=emec(2)
      e(3,3)=emec(3)
      e(1,2)=emec(4)
      e(1,3)=emec(5)
      e(2,3)=emec(6)
      e(2,1)=emec(4)
      e(3,1)=emec(5)
      e(3,2)=emec(6)
!     
!     calculating the eigenvalues and eigenvectors
!     
      n=3
      matz=1
!     
      call rs(n,n,e,w,matz,z,fv1,fv2,ier)
!     
      if(ier.ne.0) then
        write(*,*) '
     &*ERROR calculating the eigenvalues/vectors in umat_abaqusnl'
        call exit(201)
      endif
!     
!     calculating the principal stretches at the end of the increment
!     
      do i=1,3
        if(w(i).lt.-0.5d0) then
          write(*,*) '*ERROR in umat_abaqusnl: negative eigenvalue'
          write(*,*) '         of the Cauchy-Green tensor;'
          call exit(201)
        else
          w(i)=dsqrt(2.d0*w(i)+1.d0)
        endif
      enddo
!     
!     calculating the invariants at the end of the increment
!     
      v1=w(1)+w(2)+w(3)
      v2=w(1)*w(2)+w(2)*w(3)+w(3)*w(1)
      v3=w(1)*w(2)*w(3)
!     
!     calculating the right Cauchy-Green tensor at the end of the
!     increment
!     
      do i=1,3
        c(i)=2.d0*emec(i)+1.d0
      enddo
      do i=4,6
        c(i)=2.d0*emec(i)
      enddo
!     
!     calculating the square of the right Cauchy-Green tensor at the
!     end of the increment
!     
      c2(1)=c(1)*c(1)+c(4)*c(4)+c(5)*c(5)
      c2(2)=c(4)*c(4)+c(2)*c(2)+c(6)*c(6)
      c2(3)=c(5)*c(5)+c(6)*c(6)+c(3)*c(3)
      c2(4)=c(1)*c(4)+c(4)*c(2)+c(5)*c(6)
      c2(5)=c(1)*c(5)+c(4)*c(6)+c(5)*c(3)
      c2(6)=c(4)*c(5)+c(2)*c(6)+c(6)*c(3)
!     
!     calculating the right stretch tensor at the end of the increment
!     (cf. Simo and Hughes, Computational Inelasticity)
!     
!     This is the mechanical right stretch tensor, i.e. without
!     thermal expansion (since it is based on emec)
!     
      dd=v1*v2-v3
      do i=1,6
        u(i)=(-c2(i)+(v1*v1-v2)*c(i)+v1*v3*d(i))/dd
      enddo
!     
!     calculating the inverse of the right stretch tensor at the end
!     of the increment
!     
      um1(1,1)=(c(1)-v1*u(1)+v2)/v3
      um1(2,2)=(c(2)-v1*u(2)+v2)/v3
      um1(3,3)=(c(3)-v1*u(3)+v2)/v3
      um1(1,2)=(c(4)-v1*u(4))/v3
      um1(1,3)=(c(5)-v1*u(5))/v3
      um1(2,3)=(c(6)-v1*u(6))/v3
      um1(2,1)=um1(1,2)
      um1(3,1)=um1(1,3)
      um1(3,2)=um1(2,3)
!     
!     calculating the logarithmic strain at the end of the increment
!     Elog=Z.ln(w).Z^T
!     
      do i=1,3
        w(i)=dlog(w(i))
      enddo
!     
!     logarithmic strain in global coordinates at the end of the
!     increment = ln(lambda) N diadic N;
!     
!     the true logarithmic strain is ln(lambda) n diadic n;
!     the rotation tensor R = n diadic N;
!     therefore, ln(lambda) N diadic N is R^T true logarithmic strain R, 
!     in other words the corotational true logarithmic strain (required
!     by the abaqus routine)
!     
      xm1(1)=z(1,1)*z(1,1)
      xm1(2)=z(2,1)*z(2,1)
      xm1(3)=z(3,1)*z(3,1)
      xm1(4)=z(1,1)*z(2,1)
      xm1(5)=z(1,1)*z(3,1)
      xm1(6)=z(2,1)*z(3,1)
!
      xm2(1)=z(1,2)*z(1,2)
      xm2(2)=z(2,2)*z(2,2)
      xm2(3)=z(3,2)*z(3,2)
      xm2(4)=z(1,2)*z(2,2)
      xm2(5)=z(1,2)*z(3,2)
      xm2(6)=z(2,2)*z(3,2)
!
      xm3(1)=z(1,3)*z(1,3)
      xm3(2)=z(2,3)*z(2,3)
      xm3(3)=z(3,3)*z(3,3)
      xm3(4)=z(1,3)*z(2,3)
      xm3(5)=z(1,3)*z(3,3)
      xm3(6)=z(2,3)*z(3,3)
!
      eln(1)=xm1(1)*w(1)+xm2(1)*w(2)+xm3(1)*w(3)          
      eln(2)=xm1(2)*w(1)+xm2(2)*w(2)+xm3(2)*w(3)          
      eln(3)=xm1(3)*w(1)+xm2(3)*w(2)+xm3(3)*w(3)          
      eln(4)=xm1(4)*w(1)+xm2(4)*w(2)+xm3(4)*w(3)          
      eln(5)=xm1(5)*w(1)+xm2(5)*w(2)+xm3(5)*w(3)          
      eln(6)=xm1(6)*w(1)+xm2(6)*w(2)+xm3(6)*w(3)  
c      eln(1)=z(1,1)*z(1,1)*w(1)+z(1,2)*z(1,2)*w(2)+
c     &     z(1,3)*z(1,3)*w(3)          
c      eln(2)=z(2,1)*z(2,1)*w(1)+z(2,2)*z(2,2)*w(2)+
c     &     z(2,3)*z(2,3)*w(3)          
c      eln(3)=z(3,1)*z(3,1)*w(1)+z(3,2)*z(3,2)*w(2)+
c     &     z(3,3)*z(3,3)*w(3)          
c      eln(4)=z(1,1)*z(2,1)*w(1)+z(1,2)*z(2,2)*w(2)+
c     &     z(1,3)*z(2,3)*w(3)          
c      eln(5)=z(1,1)*z(3,1)*w(1)+z(1,2)*z(3,2)*w(2)+
c     &     z(1,3)*z(3,3)*w(3)          
c      eln(6)=z(2,1)*z(3,1)*w(1)+z(2,2)*z(3,2)*w(2)+
c     &     z(2,3)*z(3,3)*w(3)  
!     
      do i=1,nstate_
        xstate(i,iint,iel)=xstateini(i,iint,iel)
      enddo
!     
      temp=t1l
!     
      nprops=-kode-100
!     
!     taking local material orientations into account
!     
      if(iorien.ne.0) then
        call transformatrix(orab(1,iorien),pgauss,skl)
!     
!     rotating the strain at the end of the increment
!     into the local system
!     
        xa(1,1)=eln(1)
        xa(1,2)=eln(4)
        xa(1,3)=eln(5)
        xa(2,1)=eln(4)
        xa(2,2)=eln(2)
        xa(2,3)=eln(6)
        xa(3,1)=eln(5)
        xa(3,2)=eln(6)
        xa(3,3)=eln(3)
!     
        do jj=1,6
          stran(jj)=0.d0
          j1=kal(1,jj)
          j2=kal(2,jj)
          do j3=1,3
            do j4=1,3
              stran(jj)=stran(jj)+
     &             xa(j3,j4)*skl(j3,j1)*skl(j4,j2)
            enddo
          enddo
        enddo
      else
        do jj=1,6
          stran(jj)=eln(jj)
        enddo
      endif
!     
      if(amat(1:11).eq.'JOHNSONCOOK') then
!     
!     changing physical strain into engineering strain
!     ABAQUS uses the engineering strain!
!     
        do i=4,6
          stran(i)=2.d0*stran(i)
        enddo
!     
        call umat_johnson_cook(stre,xstate(1,iint,iel),ddsdde,
     &       spd,rpl,stran,dtime,temp,elconloc,nprops,
     &       iel,iint,jstep,kinc,icmd)
!     
      elseif(kode.eq.-50) then
!     
!     deformation plasticity
!     
        call umat_def_plas(elconloc,elas,stran,icmd,stre,
     &       xstate(1,iint,iel),iel,iint)
      endif
!     
!     rotating the stress into the global system
!     S=J.U^(-1).s'.U^(-T)  (no orientation card) or
!     S=J.T.U^(-1).s'.U^(-T).T^T (orientation card)
!     
!     
!     Here, U is the total right stretch tensor, i.e. the
!     thermal effects have to be included.
!     
!     Isotropic thermal expansion is assumed. Indeed, only
!     in that case the rotation tensor with expansion equals
!     the rotation tensor without thermal expansion, and      
!     the right stretch tensor with expansion equals the right
!     stretch tensor without expansion times (1+alpha*Delta T).
!     
!     calculate the Cauchy-Green tensor at the beginning of 
!     the increment
!     
      do i=1,3
        do j=1,3
          ctot(i,j)=xkl(1,i)*xkl(1,i)+xkl(2,i)*xkl(2,i)+
     &         xkl(3,i)*xkl(3,i)
        enddo
      enddo
!     
!     Ctot=C*(1+alpha*Delta T)**2
!     expansion:=1+alpha*Delta T
!     
      expansion=dsqrt((ctot(1,1)/c(1)+ctot(2,2)/c(2)+
     &     ctot(3,3)/c(3))/3.d0)
!     
!     calculate the inverse total right stretch tensor at the end
!     of the increment (= inverse mechanical right stretch tensor
!     divided by the expansion at the end of the increment)
!     
      do i=1,3
        do j=1,3
          um1(i,j)=um1(i,j)/expansion
        enddo
      enddo
!     
      if(iorien.ne.0) then
        do i=1,3
          do j=1,3
            tkl(i,j)=skl(i,1)*um1(1,j)+skl(i,2)*um1(2,j)+
     &           skl(i,3)*um1(3,j)
          enddo
        enddo
      else
        do i=1,3
          do j=1,3
            tkl(i,j)=um1(i,j)
          enddo
        enddo
      endif
!     
      xa(1,1)=stre(1)
      xa(1,2)=stre(4)
      xa(1,3)=stre(5)
      xa(2,1)=stre(4)
      xa(2,2)=stre(2)
      xa(2,3)=stre(6)
      xa(3,1)=stre(5)
      xa(3,2)=stre(6)
      xa(3,3)=stre(3)
!     
      do jj=1,6
        stre(jj)=0.d0
        j1=kal(1,jj)
        j2=kal(2,jj)
        do j3=1,3
          do j4=1,3
            stre(jj)=stre(jj)+
     &           xa(j3,j4)*tkl(j1,j3)*tkl(j2,j4)
          enddo
        enddo
        stre(jj)=stre(jj)*vj
      enddo
!     
!     calculate the stiffness matrix (the matrix is symmetrized)
!     
      if(icmd.ne.3) then
        if(amat(1:11).eq.'JOHNSONCOOK') then
          elas(1)=ddsdde(1,1)
          elas(2)=ddsdde(2,1)
          elas(3)=ddsdde(2,2)
          elas(4)=ddsdde(3,1)
          elas(5)=ddsdde(3,2)
          elas(6)=ddsdde(3,3)
          elas(7)=ddsdde(4,1)
          elas(8)=ddsdde(4,2)
          elas(9)=ddsdde(4,3)
          elas(10)=ddsdde(4,4)
          elas(11)=ddsdde(5,1)
          elas(12)=ddsdde(5,2)
          elas(13)=ddsdde(5,3)
          elas(14)=ddsdde(5,4)
          elas(15)=ddsdde(5,5)
          elas(16)=ddsdde(6,1)
          elas(17)=ddsdde(6,2)
          elas(18)=ddsdde(6,3)
          elas(19)=ddsdde(6,4)
          elas(20)=ddsdde(6,5)
          elas(21)=ddsdde(6,6)
        endif
!
!        calculate the eigenvalues of the Green-Cauchy tensor
!        (= the square of the stretches)
!
        iflag=0
        if(iflag.eq.1) then
          do i=1,3
            al(i)=(dexp(w(i)))**2
          enddo
!     
          if((al(1)-al(2).lt.1.d-10).or.(al(1)-al(3).lt.1.d-10)) then
            write(*,*) 'three equal eigenvalues'
!     
!     three equal eigenvalues
!     
            factor=1.d0/(2.d0*al(1))
            do jj=1,21
              delndc(jj)=factor*dp(jj)
            enddo
            factor=-1.d0/(2.d0*al(1)*dsqrt(al(1)))
            do jj=1,21
              dum1dc(jj)=factor*dp(jj)
            enddo
          elseif((al(1)-al(2).lt.1.d-10).or.
     &           (al(2)-al(3).lt.1.d-10).or.
     &           (al(3)-al(1).lt.1.d-10)) then
!     
!     two equal eigenvalues: TO DO
!     
            write(*,*) 'two equal eigenvalues'
          else
!     
!     three different eigenvalues
!     
!     calculating (C diadic I + I diadic C)'
!     
            ciicp(1)=2.d0*c(1)
            ciicp(2)=0.d0
            ciicp(3)=2.d0*c(2)
            ciicp(4)=0.d0
            ciicp(5)=0.d0
            ciicp(6)=2.d0*c(3)
            ciicp(7)=c(4)
            ciicp(8)=c(4)
            ciicp(9)=0.d0
            ciicp(10)=(c(1)+c(2))/2.d0
            ciicp(11)=c(5)
            ciicp(12)=0.d0
            ciicp(13)=c(5)
            ciicp(14)=c(6)/2.d0
            ciicp(15)=(c(1)+c(3))/2.d0
            ciicp(16)=0.d0
            ciicp(17)=c(6)
            ciicp(18)=c(6)
            ciicp(19)=c(5)/2.d0
            ciicp(20)=c(4)/2.d0
            ciicp(21)=(c(2)+c(3))/2.d0
!     
!     calculating M_i diadic M_i (i=1,2,3)
!     
            jj=0
            do j2=1,6
              do j1=1,j2
                jj=jj+1
                xm1m1(jj)=xm1(j1)*xm1(j2)
                xm2m2(jj)=xm2(j1)*xm2(j2)
                xm3m3(jj)=xm3(j1)*xm3(j2)
              enddo
            enddo
!     
!     calculating tensor A and B (Dhondt and Hackenberg, 2021)
!     
            do jj=1,21
              aa(jj)=dp(jj)-xm1m1(jj)-xm2m2(jj)-xm3m3(jj)
              bb(jj)=ciicp(jj)-2.d0*(al(1)*xm1m1(jj)+
     &             al(2)*xm2m2(jj)+al(3)*xm3m3(jj))
            enddo
!     
!     calculating the derivative of the logarithmic strain               
!     tensor w.r.t. the Cauchy-Green tensor               
!     
            f1=1.d0/(2.d0*al(1))
            f2=1.d0/(2.d0*al(2))
            f3=1.d0/(2.d0*al(3))
            d1=(al(1)-al(2))*(al(1)-al(3))
            d2=(al(2)-al(3))*(al(2)-al(1))
            d3=(al(3)-al(1))*(al(3)-al(2))
            g1=w(1)/d1
            g2=w(2)/d2
            g3=w(3)/d3
            h1=(al(2)+al(3))
            h2=(al(3)+al(1))
            h3=(al(1)+al(2))
            do jj=1,21
              delndc(jj)=f1*xm1m1(jj)+f2*xm2m2(jj)+f3*xm3m3(jj)+
     &             (g1+g2+g3)*bb(jj)-(h1*g1+h2*g2+h3*g3)*aa(jj)
            enddo
     !             
            g1=1.d0/dsqrt(al(1))
            g2=1.d0/dsqrt(al(2))
            g3=1.d0/dsqrt(al(3))
            f1=-g1/(2.d0*al(1))
            f2=-g2/(2.d0*al(2))
            f3=-g3/(2.d0*al(3))
            g1=g1/d1
            g2=g2/d2
            g3=g3/d3
            do jj=1,21
              dum1dc(jj)=(f1*xm1m1(jj)+f2*xm2m2(jj)+f3*xm3m3(jj)+
     &             (g1+g2+g3)*bb(jj)-(h1*g1+h2*g2+h3*g3)*aa(jj))/
     &             expansion
            enddo
!     
          endif
!     
!     transform delndc, dum1dc and stiff into a (3,3,3,3)-maxtrix
!     these matrices are symmetric
!     
          call anisotropic(delndc,delndcglob)
          call anisotropic(dum1dc,dum1dcglob)
          call anisotropic(elas,ya)
!     
!     calculate dsigdc(3,3,3,3); this matrix may be asymmetric
!     
          do j1=1,3
            do j2=1,3
              do j3=1,3
                do j4=1,3
                  dsigdc(j1,j2,j3,j4)=0.d0
                  do j5=1,3
                    do j6=1,3
                      dsigdc(j1,j2,j3,j4)=dsigdc(j1,j2,j3,j4)+
     &                     ya(j1,j2,j5,j6)*delndcglob(j5,j6,j3,j4)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
!     
!     calculate the inverse of the Cauchy-Green tensor
!     C^{-1}=U^{-1}.U^{-1}
!     
          do j1=1,3
            do j2=j1,3
              cm1(j1,j2)=um1(j1,1)*um1(1,j2)+um1(j1,2)*um1(2,j2)+
     &             um1(j1,3)*um1(3,j2)
            enddo
          enddo
          cm1(2,1)=cm1(1,2)
          cm1(3,1)=cm1(1,3)
          cm1(3,2)=cm1(2,3)
!     
!     calculate the stiffness matrix (at first without transformation)
!     
          do jj=1,36
            j1=keltot(1,jj)
            j2=keltot(2,jj)
            j3=keltot(3,jj)
            j4=keltot(4,jj)
            elastot(jj)=0.d0
            do j5=1,3
              do j6=1,3
                elastot(jj)=elastot(jj)+
     &               um1(j1,j5)*xa(j5,j6)*um1(j6,j2)*cm1(j3,j4)/2.d0+
     &               (dum1dcglob(j1,j5,j3,j4)*xa(j5,j6)*um1(j6,j2)+
     &               um1(j1,j5)*dsigdc(j5,j6,j3,j4)*um1(j6,j2)+
     &               um1(j1,j5)*xa(j5,j6)*dum1dcglob(j6,j2,j3,j4))/
     &               expansion**2
              enddo
            enddo
            elastot(jj)=elastot(jj)*vj
          enddo
!     
!     symmetrize the stiffness matrix and multiply by 2 (derivative
!     w.r.t. E is needed, not w.r.t. C)
!     
          elas(1)=2.d0*elastot(1)
          elas(2)=(elastot(2)+elastot(22))
          elas(3)=2.d0*elastot(3)
          elas(4)=(elastot(4)+elastot(23))
          elas(5)=(elastot(5)+elastot(24))
          elas(6)=2.d0*elastot(6)
          elas(7)=(elastot(7)+elastot(25))
          elas(8)=(elastot(8)+elastot(26))
          elas(9)=(elastot(9)+elastot(27))
          elas(10)=2.d0*elastot(10)
          elas(11)=(elastot(11)+elastot(28))
          elas(12)=(elastot(12)+elastot(29))
          elas(13)=(elastot(13)+elastot(30))
          elas(14)=(elastot(14)+elastot(31))
          elas(15)=2.d0*elastot(15)
          elas(16)=(elastot(16)+elastot(32))
          elas(17)=(elastot(17)+elastot(33))
          elas(18)=(elastot(18)+elastot(34))
          elas(19)=(elastot(19)+elastot(35))
          elas(20)=(elastot(20)+elastot(36))
          elas(21)=2.d0*elastot(21)
        else
!     
!     rotating the stiffness coefficients into the global system
!     
          call anisotropic(elas,ya)
!     
          do jj=1,21
            j1=kel(1,jj)
            j2=kel(2,jj)
            j3=kel(3,jj)
            j4=kel(4,jj)
            elas(jj)=0.d0
            do j5=1,3
              do j6=1,3
                do j7=1,3
                  do j8=1,3
                    elas(jj)=elas(jj)+ya(j5,j6,j7,j8)*
     &                   tkl(j1,j5)*tkl(j2,j6)*tkl(j3,j7)*tkl(j4,j8)
                  enddo
                enddo
              enddo
            enddo
            elas(jj)=elas(jj)*vj
          enddo
        endif
      endif
!     
      return
      end
