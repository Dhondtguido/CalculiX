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
      subroutine rubber(elconloc,elas,emec,kode,didc,
     &  d2idc2,dibdc,d2ibdc2,dudc,d2udc2,dldc,d2ldc2,dlbdc,d2lbdc2,
     &  ithermal,icmd,beta,stre)
!
!     calculates stiffness and stresses for rubber and elastomeric
!     foam materials
!
!     icmd=3: stress at mechanical strain
!     else: stress and stiffness matrix at mechanical strain
!
      implicit none
!
      integer ogden,hyperfoam,taylor
!
      integer nelconst,kode,kk(84),i,j,k,l,m,nt,icmd,istart,iend,
     &  nc,n,ithermal(*),ii,jj,mm,neig
!
      real*8 elconloc(*),elas(*),emec(*),didc(3,3,3),
     &  d2idc2(3,3,3,3,3),dibdc(3,3,3),d2ibdc2(3,3,3,3,3),dudc(3,3),
     &  d2udc2(3,3,3,3),dldc(3,3,3),d2ldc2(3,3,3,3,3),dlbdc(3,3,3),
     &  d2lbdc2(3,3,3,3,3),v1,v2,v3,c(3,3),cinv(3,3),d(3,3),djth,
     &  coef,bb,cc,cm,cn,tt,pi,dd,al(3),v1b,v2b,v3b,
     &  alb(3),beta(*),v33,v36,all(3),term,stre(*),total,coefa,
     &  coefb,coefd,coefm,constant(21)
!
      kk=(/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &  1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,3,3,1,3,
     &  1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,1,2,2,3,1,3,2,3,
     &  2,3,2,3/)
!
!     copy the elastic constants into a new field, such that
!     they can be mixed without influencing the field in the
!     calling program
!
      do i=1,21
         constant(i)=elconloc(i)
      enddo
!
!     type of hyperelastic law; taylor stands for everything
!     which involves parts of a taylor expansion in terms of the
!     reduced Green deformation invariants
!
      ogden=0
      hyperfoam=0
      taylor=0
      if((kode.lt.-3).and.(kode.gt.-7)) then
         ogden=1
      elseif((kode.lt.-14).and.(kode.gt.-18)) then
         hyperfoam=1
      else
         taylor=1
      endif
!
c      if(icmd.eq.1) then
         istart=1
         iend=1
!
!     reclassifying some classes of hyperelastic materials as
!     subclasses of the polynomial model
!
      if(((kode.lt.-1).and.(kode.gt.-4)).or.
     &   ((kode.lt.-6).and.(kode.gt.-13)).or.
     &    (kode.eq.-14)) then
         if(kode.eq.-2) then
            kode=-7
            nelconst=1
         elseif((kode.eq.-3).or.(kode.eq.-10)) then
            constant(3)=constant(2)
            constant(2)=0.d0
            kode=-7
            nelconst=1
         elseif(kode.eq.-11) then
            constant(7)=constant(4)
            constant(6)=constant(3)
            constant(5)=0.d0
            constant(4)=0.d0
            constant(3)=constant(2)
            constant(2)=0.d0
            kode=-8
            nelconst=2
         elseif((kode.eq.-12).or.(kode.eq.-14)) then
            constant(12)=constant(6)
            constant(11)=constant(5)
            constant(10)=constant(4)
            constant(9)=0.d0
            constant(8)=0.d0
            constant(7)=0.d0
            constant(6)=constant(3)
            constant(5)=0.d0
            constant(4)=0.d0
            constant(3)=constant(2)
            constant(2)=0.d0
            kode=-9
            nelconst=3
         elseif(kode.eq.-7) then
            nelconst=1
         elseif(kode.eq.-8) then
            nelconst=2
         elseif(kode.eq.-9) then
            nelconst=3
         endif
      endif
!
!     major loop
!
      do ii=istart,iend
!
!        calculation of the Green deformation tensor for the total
!        strain and the thermal strain
!
         do i=1,3
            c(i,i)=emec(i)*2.d0+1.d0
         enddo
         c(1,2)=2.d0*emec(4)
         c(1,3)=2.d0*emec(5)
         c(2,3)=2.d0*emec(6)
!
!        calculation of the invariants of c
!
         v1=c(1,1)+c(2,2)+c(3,3)
         v2=c(2,2)*c(3,3)+c(1,1)*c(3,3)+c(1,1)*c(2,2)-
     &     (c(2,3)*c(2,3)+c(1,3)*c(1,3)+c(1,2)*c(1,2))
         v3=c(1,1)*(c(2,2)*c(3,3)-c(2,3)*c(2,3))
     &     -c(1,2)*(c(1,2)*c(3,3)-c(1,3)*c(2,3))
     &     +c(1,3)*(c(1,2)*c(2,3)-c(1,3)*c(2,2))
         v33=v3**(-1.d0/3.d0)
         v36=v3**(-1.d0/6.d0)
!
!        calculation of the thermal strain jacobian 
!        (not really needed)
!
         djth=1.d0
!
!        inversion of c
!
         cinv(1,1)=(c(2,2)*c(3,3)-c(2,3)*c(2,3))/v3
         cinv(2,2)=(c(1,1)*c(3,3)-c(1,3)*c(1,3))/v3
         cinv(3,3)=(c(1,1)*c(2,2)-c(1,2)*c(1,2))/v3
         cinv(1,2)=(c(1,3)*c(2,3)-c(1,2)*c(3,3))/v3
         cinv(1,3)=(c(1,2)*c(2,3)-c(2,2)*c(1,3))/v3
         cinv(2,3)=(c(1,2)*c(1,3)-c(1,1)*c(2,3))/v3
         cinv(2,1)=cinv(1,2)
         cinv(3,1)=cinv(1,3)
         cinv(3,2)=cinv(2,3)
!
!        creation of the delta Dirac matrix d
!
         do j=1,3
            do i=1,3
               d(i,j)=0.d0
            enddo
         enddo
         do i=1,3
            d(i,i)=1.d0
         enddo
!
!        derivative of the c-invariants with respect to c(k,l)
!
         do l=1,3
            do k=1,l
               didc(k,l,1)=d(k,l)
               didc(k,l,2)=v1*d(k,l)-c(k,l)
               didc(k,l,3)=v3*cinv(k,l)
            enddo
         enddo
!
!        second derivative of the c-invariants w.r.t. c(k,l) 
!        and c(m,n)
!
         if(icmd.ne.3) then
            nt=0
            do i=1,21
               k=kk(nt+1)
               l=kk(nt+2)
               m=kk(nt+3)
               n=kk(nt+4)
               nt=nt+4
               d2idc2(k,l,m,n,1)=0.d0
               d2idc2(k,l,m,n,2)=d(k,l)*d(m,n)-
     &                           (d(k,m)*d(l,n)+d(k,n)*d(l,m))/2.d0
               d2idc2(k,l,m,n,3)=v3*(cinv(m,n)*cinv(k,l)-
     &            (cinv(k,m)*cinv(n,l)+cinv(k,n)*cinv(m,l))/2.d0)
            enddo
         endif
!
!     derivatives for the reduced invariants used in rubber materials
!
         v1b=v1*v33
         v2b=v2*v33*v33
         v3b=dsqrt(v3)/djth
!
!        first derivative of the reduced c-invariants w.r.t. c(k,l)
!
         do l=1,3
            do k=1,l
               if(taylor.eq.1) then
                  dibdc(k,l,1)=-v33**4*v1*didc(k,l,3)/3.d0
     &                 +v33*didc(k,l,1)
                  dibdc(k,l,2)=-2.d0*v33**5*v2*didc(k,l,3)/3.d0
     &                 +v33**2*didc(k,l,2)
               endif
               dibdc(k,l,3)=didc(k,l,3)/(2.d0*dsqrt(v3)*djth)
            enddo
         enddo
!
!        second derivative of the reduced c-invariants w.r.t. c(k,l)
!        and c(m,n)
!
         if(icmd.ne.3) then
            nt=0
            do i=1,21
               k=kk(nt+1)
               l=kk(nt+2)
               m=kk(nt+3)
               n=kk(nt+4)
               nt=nt+4
               if(taylor.eq.1) then
                  d2ibdc2(k,l,m,n,1)=4.d0/9.d0*v33**7*v1*didc(k,l,3)
     &                 *didc(m,n,3)-v33**4/3.d0*(didc(m,n,1)*didc(k,l,3)
     &                 +didc(k,l,1)*didc(m,n,3))-v33**4/3.d0*v1*
     &                 d2idc2(k,l,m,n,3)+v33*d2idc2(k,l,m,n,1)
                  d2ibdc2(k,l,m,n,2)=10.d0*v33**8/9.d0*v2*didc(k,l,3)
     &                 *didc(m,n,3)-2.d0*v33**5/3.d0*(didc(m,n,2)
     &                 *didc(k,l,3)
     &                 +didc(k,l,2)*didc(m,n,3))-2.d0*v33**5/3.d0*v2*
     &                 d2idc2(k,l,m,n,3)+v33**2*d2idc2(k,l,m,n,2)
               endif
               d2ibdc2(k,l,m,n,3)=-didc(k,l,3)*didc(m,n,3)/
     &              (4.d0*djth*v3**1.5d0)+d2idc2(k,l,m,n,3)/
     &              (2.d0*dsqrt(v3)*djth)
            enddo
         endif
!
!     calculation of the principal stretches for the Ogden model and 
!     hyperfoam materials
!
         if((ogden.eq.1).or.(hyperfoam.eq.1)) then
!
!        taking the thermal jacobian into account        
!
            if((kode.lt.-14).and.(kode.gt.-18)) then
               dd=djth**(1.d0/3.d0)
            else
               dd=1.d0
            endif
!
            pi=4.d0*datan(1.d0)
!
!        determining the eigenvalues of c (Simo & Hughes) and taking
!        the square root to obtain the principal stretches
!
!           neig is the number of different eigenvalues
!
            neig=3
!
            bb=v2-v1*v1/3.d0
            cc=-2.d0*v1**3/27.d0+v1*v2/3.d0-v3
            if(dabs(bb).le.1.d-10) then
               if(dabs(cc).gt.1.d-10) then
                  al(1)=-cc**(1.d0/3.d0)
               else
                  al(1)=0.d0
               endif
               al(2)=al(1)
               al(3)=al(1)
               neig=1
            else
               cm=2.d0*dsqrt(-bb/3.d0)
               cn=3.d0*cc/(cm*bb)
               if(dabs(cn).gt.1.d0) then
                  if(cn.gt.1.d0) then
                     cn=1.d0
                  else
                     cn=-1.d0
                  endif
               endif
               tt=datan2(dsqrt(1.d0-cn*cn),cn)/3.d0
               al(1)=dcos(tt)
               al(2)=dcos(tt+2.d0*pi/3.d0)
               al(3)=dcos(tt+4.d0*pi/3.d0)
!
!              check for two equal eigenvalues
!
               if((dabs(al(1)-al(2)).lt.1.d-5).or.
     &            (dabs(al(1)-al(3)).lt.1.d-5).or.
     &            (dabs(al(2)-al(3)).lt.1.d-5)) neig=2
               al(1)=cm*al(1)
               al(2)=cm*al(2)
               al(3)=cm*al(3)
            endif
            do i=1,3
               al(i)=dsqrt(al(i)+v1/3.d0)
               all(i)=(6.d0*al(i)**5-4.d0*v1*al(i)**3+2.d0*al(i)*v2)*dd
            enddo
!
!        first derivative of the principal stretches w.r.t. c(k,l)
!
            if(neig.eq.3) then
!
!              three different principal stretches
!
               do i=1,3
                  do l=1,3
                     do k=1,l
                        dldc(k,l,i)=(al(i)**4*didc(k,l,1)
     &                       -al(i)**2*didc(k,l,2)+didc(k,l,3))/all(i)
                     enddo
                  enddo
               enddo
            elseif(neig.eq.1) then
!
!              three equal principal stretches
!
               do i=1,3
                  do l=1,3
                     do k=1,l
                        dldc(k,l,i)=didc(k,l,1)/(6.d0*al(i))
                     enddo
                  enddo
               enddo
            else
!
!              two equal principal stretches
!
               do i=1,3
                  do l=1,3
                     do k=1,l
                        dldc(k,l,i)=(dcos(tt+(i-1)*2.d0*pi/3.d0)*
     &                    (2.d0*v1*didc(k,l,1)-3.d0*didc(k,l,2))/
     &                    (3.d0*dsqrt(v1*v1-3.d0*v2))+didc(k,l,1)/3.d0)/
     &                    (2.d0*al(i))
                     enddo
                  enddo
               enddo
            endif
!
!        second derivative of the principal stretches w.r.t. c(k,l)
!        and c(m,n)
!
            if(icmd.ne.3) then
               if(neig.eq.3) then
!
!              three different principal stretches
!
                  do i=1,3
                     nt=0
                     do j=1,21
                        k=kk(nt+1)
                        l=kk(nt+2)
                        m=kk(nt+3)
                        n=kk(nt+4)
                        nt=nt+4
                        d2ldc2(k,l,m,n,i)=(-30.d0*al(i)**4
     &                    *dldc(k,l,i)*dldc(m,n,i)+al(i)**4
     &                    *d2idc2(k,l,m,n,1)
     &                    +4.d0*al(i)**3*(didc(k,l,1)*dldc(m,n,i)
     &                    +didc(m,n,1)
     &                    *dldc(k,l,i))+12.d0*v1*al(i)**2*dldc(k,l,i)*
     &                    dldc(m,n,i)-d2idc2(k,l,m,n,2)*al(i)**2-2.d0
     &                    *al(i)*
     &                    didc(k,l,2)*dldc(m,n,i)-2.d0*v2*dldc(k,l,i)*
     &                    dldc(m,n,i)-2.d0*al(i)*didc(m,n,2)*dldc(k,l,i)
     &                    +d2idc2(k,l,m,n,3))/all(i)
                     enddo
                  enddo
               elseif(neig.eq.1) then
!
!              three equal principal stretches
!
                  do i=1,3
                     nt=0
                     do j=1,21
                        k=kk(nt+1)
                        l=kk(nt+2)
                        m=kk(nt+3)
                        n=kk(nt+4)
                        nt=nt+4
                        d2ldc2(k,l,m,n,i)=(d2idc2(k,l,m,n,1)/6.d0
     &                    -dldc(k,l,i)*dldc(m,n,i))/al(i)
                     enddo
                  enddo
               else
!
!              two equal principal stretches
!
                  do i=1,3
                     nt=0
                     do j=1,21
                        k=kk(nt+1)
                        l=kk(nt+2)
                        m=kk(nt+3)
                        n=kk(nt+4)
                        nt=nt+4
                        d2ldc2(k,l,m,n,i)=(dcos(tt+(i-1)*2.d0*pi/3.d0)*
     &                   (-(2.d0*v1*didc(k,l,1)-3.d0*didc(k,l,2))*
     &                    (2.d0*v1*didc(m,n,1)-3.d0*didc(m,n,2))/
     &                    (6.d0*(v1*v1-3.d0*v2)**1.5d0)+
     &                    (2.d0*didc(k,l,1)*didc(m,n,1)+2.d0*v1*
     &                     d2idc2(k,l,m,n,1)-3.d0*d2idc2(k,l,m,n,2))/
     &                     (3.d0*dsqrt(v1*v1-3.d0*v2)))
     &                     +d2idc2(k,l,m,n,1)/3.d0)/(2.d0*al(i))-
     &                     dldc(k,l,i)*dldc(m,n,i)/al(i)
                     enddo
                  enddo
               endif
            endif
!
!        reduced principal stretches (Ogden model)
!
            if(ogden.eq.1) then
!
!           calculation of the reduced principal stretches
!
               do i=1,3
                  alb(i)=al(i)*v36
               enddo
!
!           first derivative of the reduced principal stretches
!           w.r.t. c(k,l)
!
               do i=1,3
                  do l=1,3
                     do k=1,l
                        dlbdc(k,l,i)=-v36**7*al(i)*didc(k,l,3)/6.d0
     &                       +v36*dldc(k,l,i)
                     enddo
                  enddo
               enddo
!
!           second derivative of the reduced principal stretches w.r.t.
!           c(k,l) and c(m,n)
!
               if(icmd.ne.3) then
                  do i=1,3
                     nt=0
                     do j=1,21
                        k=kk(nt+1)
                        l=kk(nt+2)
                        m=kk(nt+3)
                        n=kk(nt+4)
                        nt=nt+4
                        d2lbdc2(k,l,m,n,i)=7.d0*v36**13*al(i)
     &                       *didc(k,l,3)*didc(m,n,3)/36.d0-v36**7/6.d0
     &                       *(dldc(m,n,i)*didc(k,l,3)+al(i)
     &                       *d2idc2(k,l,m,n,3)+dldc(k,l,i)*didc(m,n,3))
     &                       +v36*d2ldc2(k,l,m,n,i)
                     enddo
                  enddo
               endif
!
            endif
         endif
!
!     calculation of the local stiffness matrix, and, if appropriate,
!     the stresses
!                     
!     Polynomial model
!
         if((kode.lt.-6).and.(kode.gt.-10)) then
!
!        first derivative of U w.r.t. c(k,l)
!
            do l=1,3
               do k=1,l
                  dudc(k,l)=0.d0
               enddo
            enddo
!
            nc=0
            do m=1,nelconst
               do j=0,m
                  i=m-j
                  nc=nc+1
                  coef=constant(nc)
                  if(dabs(coef).lt.1.d-20) cycle
                  do l=1,3
                     do k=1,l
                        total=0.d0
                        if(i.gt.0) then
                           term=dibdc(k,l,1)
                           if(i.gt.1) term=i*term*(v1b-3.d0)**(i-1)
                           if(j.gt.0) term=term*(v2b-3.d0)**j
                           total=total+term
                        endif
                        if(j.gt.0) then
                           term=dibdc(k,l,2)
                           if(i.gt.0) term=term*(v1b-3.d0)**i
                           if(j.gt.1) term=j*term*(v2b-3.d0)**(j-1)
                           total=total+term
                        endif
                        dudc(k,l)=dudc(k,l)+total*coef
                     enddo
                  enddo
               enddo
            enddo
            do m=1,nelconst
               nc=nc+1
               coef=constant(nc)
               do l=1,3
                  do k=1,l
                     dudc(k,l)=dudc(k,l)+2.d0*m*(v3b-1.d0)**
     &                    (2*m-1)*dibdc(k,l,3)/coef
                  enddo
               enddo
            enddo
!
!        tangent stiffness matrix
!        second derivative of U w.r.t. c(k,l) and c(m,n)
!
            if(icmd.ne.3) then
               nt=0
               do i=1,21
                  k=kk(nt+1)
                  l=kk(nt+2)
                  m=kk(nt+3)
                  n=kk(nt+4)
                  nt=nt+4
                  d2udc2(k,l,m,n)=0.d0
               enddo
               nc=0
               do mm=1,nelconst
                  do j=0,mm
                     i=mm-j
                     nc=nc+1
                     coef=constant(nc)
                     if(dabs(coef).lt.1.d-20) cycle
                     nt=0
                     do jj=1,21
                        k=kk(nt+1)
                        l=kk(nt+2)
                        m=kk(nt+3)
                        n=kk(nt+4)
                        nt=nt+4
                        total=0.d0
                        if(i.gt.1) then
                           term=dibdc(k,l,1)*dibdc(m,n,1)*i*(i-1)
                           if(i.gt.2) term=term*(v1b-3.d0)**(i-2)
                           if(j.gt.0) term=term*(v2b-3.d0)**j
                           total=total+term
                        endif
                        if((i.gt.0).and.(j.gt.0)) then
                           term=dibdc(k,l,1)*dibdc(m,n,2)+
     &                          dibdc(m,n,1)*dibdc(k,l,2)
                           if(i.gt.1) term=i*term*(v1b-3.d0)**(i-1)
                           if(j.gt.1) term=j*term*(v2b-3.d0)**(j-1)
                           total=total+term
                        endif
                        if(i.gt.0) then
                           term=d2ibdc2(k,l,m,n,1)
                           if(i.gt.1) term=i*term*(v1b-3.d0)**(i-1)
                           if(j.gt.0) term=term*(v2b-3.d0)**j
                           total=total+term
                        endif
                        if(j.gt.1) then
                           term=dibdc(k,l,2)*dibdc(m,n,2)*j*(j-1)
                           if(i.gt.0) term=term*(v1b-3.d0)**i
                           if(j.gt.2) term=term*(v2b-3.d0)**(j-2)
                           total=total+term
                        endif
                        if(j.gt.0) then
                           term=d2ibdc2(k,l,m,n,2)
                           if(i.gt.0) term=term*(v1b-3.d0)**i
                           if(j.gt.1) term=j*term*(v2b-3.d0)**(j-1)
                           total=total+term
                        endif
                        d2udc2(k,l,m,n)=d2udc2(k,l,m,n)+total*coef
                     enddo
                  enddo
               enddo
!
               do mm=1,nelconst
                  nc=nc+1
                  coef=constant(nc)
                  nt=0
                  do i=1,21
                     k=kk(nt+1)
                     l=kk(nt+2)
                     m=kk(nt+3)
                     n=kk(nt+4)
                     nt=nt+4
                     if(mm.eq.1) then
                        term=(2.d0*dibdc(k,l,3)*dibdc(m,n,3)+
     &                    2.d0*(v3b-1.d0)*d2ibdc2(k,l,m,n,3))/coef
                     else
                        term=
     &                    2.d0*mm*(v3b-1.d0)**(2*mm-2)/coef*
     &                    ((2*mm-1)*dibdc(k,l,3)*dibdc(m,n,3)
     &                    +(v3b-1.d0)*d2ibdc2(k,l,m,n,3))
                     endif
                     d2udc2(k,l,m,n)=d2udc2(k,l,m,n)+term
                  enddo
               enddo
            endif
         endif
!
!     Ogden form
!
         if((kode.lt.-3).and.(kode.gt.-7)) then
            if(kode.eq.-4) then
               nelconst=1
            elseif(kode.eq.-5) then
               nelconst=2
            elseif(kode.eq.-6) then
               nelconst=3
            endif
!
!        first derivative of U w.r.t. c(k,l)
!
            do l=1,3
               do k=1,l
                  dudc(k,l)=0.d0
               enddo
            enddo
!
            do m=1,nelconst
               coefa=constant(2*m)
               coefd=constant(2*nelconst+m)
               coefm=constant(2*m-1)
               do l=1,3
                  do k=1,l
                     term=0.d0
                     do i=1,3
                        term=term+alb(i)**(coefa-1.d0)*dlbdc(k,l,i)
                     enddo
                     dudc(k,l)=dudc(k,l)+2.d0*coefm/coefa
     &                    *term+2.d0*m/coefd*
     &                    (v3b-1.d0)**(2*m-1)*dibdc(k,l,3)
                  enddo
               enddo
            enddo
!
!        tangent stiffness matrix
!        second derivative of U w.r.t. c(k,l) and c(m,n)
!
            if(icmd.ne.3) then
               nt=0
               do i=1,21
                  k=kk(nt+1)
                  l=kk(nt+2)
                  m=kk(nt+3)
                  n=kk(nt+4)
                  nt=nt+4
                  d2udc2(k,l,m,n)=0.d0
               enddo
               do mm=1,nelconst
                  coefa=constant(2*mm)
                  coefd=constant(2*nelconst+mm)
                  coefm=constant(2*mm-1)
                  nt=0
                  do jj=1,21
                     k=kk(nt+1)
                     l=kk(nt+2)
                     m=kk(nt+3)
                     n=kk(nt+4)
                     nt=nt+4
                     term=0.d0
                     do i=1,3
                        term=term+alb(i)**(coefa-2.d0)*dlbdc(k,l,i)*
     &                       dlbdc(m,n,i)
                     enddo
                     term=term*(coefa-1.d0)
                     do i=1,3
                        term=term+alb(i)**(coefa-1.d0)
     &                        *d2lbdc2(k,l,m,n,i)
                     enddo
                     term=term*2.d0*coefm/coefa
                     d2udc2(k,l,m,n)=d2udc2(k,l,m,n)+term+(2*mm)*
     &                    (2*mm-1)/coefd*(v3b-1.d0)**(2*mm-2)*
     &                    dibdc(k,l,3)*dibdc(m,n,3)+2*mm/coefd
     &                    *(v3b-1.d0)**(2*mm-1)*d2ibdc2(k,l,m,n,3)
                  enddo
               enddo
            endif
         endif
!
!     Arruda-Boyce model
!
         if(kode.eq.-1) then
            coef=constant(2)
!
!        first derivative of U w.r.t. c(k,l)
!
            do l=1,3
               do k=1,l
                  dudc(k,l)=constant(1)*(0.5d0+v1b/(10.d0*
     &                 coef**2)+33.d0*v1b*v1b/(1050.d0*
     &                 coef**4)+76.d0*v1b**3/(7000.d0*
     &                 coef**6)+2595.d0*v1b**4/(673750.d0*
     &                 coef**8))*dibdc(k,l,1)+(v3b-1.d0/v3b)
     &                 *dibdc(k,l,3)/constant(3)
               enddo
            enddo
!
!        tangent stiffness matrix
!        second derivative of U w.r.t. c(k,l) and c(m,n)
!
            if(icmd.ne.3) then
               nt=0
               do jj=1,21
                  k=kk(nt+1)
                  l=kk(nt+2)
                  m=kk(nt+3)
                  n=kk(nt+4)
                  nt=nt+4
                  d2udc2(k,l,m,n)=constant(1)*(1.d0/(10.d0*
     &                 coef**2)+66.d0*v1b/(1050.d0*coef**4)+228.d0
     &                 *v1b**2/(7000.d0*coef**6)+10380.d0*v1b**3/
     &                 (673750.d0*coef**8))*dibdc(k,l,1)*dibdc(m,n,1)
     &                 +constant(1)*(0.5d0+v1b/(10.d0*coef**2)
     &                 +33.d0*v1b**2/
     &                 (1050.d0*coef**4)+76.d0*v1b**3/(7000.d0*coef**6)+
     &                 2595.d0*v1b**4/(673750.d0*coef**8))
     &                 *d2ibdc2(k,l,m,n,1)
     &                 +(1.d0+1.d0/v3b**2)*dibdc(k,l,3)*dibdc(m,n,3)/
     &                 constant(3)+(v3b-1.d0/v3b)*d2ibdc2(k,l,m,n,3)
     &                 /constant(3)
               enddo
            endif
         endif
!
!     elastomeric foam behavior
!
         if((kode.lt.-15).and.(kode.gt.-18)) then
            if(kode.eq.-15) then
               nelconst=1
            elseif(kode.eq.-16) then
               nelconst=2
            elseif(kode.eq.-17) then
               nelconst=3
            endif
!
!        first derivative of U w.r.t. c(k,l)
!
            do l=1,3
               do k=1,l
                  dudc(k,l)=0.d0
               enddo
            enddo
!
            do m=1,nelconst
               coefa=constant(2*m)
               coefb=constant(2*nelconst+m)/(1.d0-2.d0
     &              *constant(2*nelconst+m))
               coefm=constant(2*m-1)
               do l=1,3
                  do k=1,l
                     term=0.d0
                     do i=1,3
                        term=term+al(i)**(coefa-1.d0)*dldc(k,l,i)
                     enddo
                     dudc(k,l)=dudc(k,l)+2.d0*coefm/coefa
     &                    *(term-v3b**(-coefa*coefb-1.d0)*
     &                    dibdc(k,l,3))
                  enddo
               enddo
            enddo
!
!        tangent stiffness matrix
!        second derivative of U w.r.t. c(k,l) and c(m,n)
!
            if(icmd.ne.3) then
               nt=0
               do i=1,21
                  k=kk(nt+1)
                  l=kk(nt+2)
                  m=kk(nt+3)
                  n=kk(nt+4)
                  nt=nt+4
                  d2udc2(k,l,m,n)=0.d0
               enddo
               do mm=1,nelconst
                  coefa=constant(2*mm)
                  coefb=constant(2*nelconst+mm)/(1.d0-2.d0
     &                      *constant(2*nelconst+mm))
                  coefm=constant(2*mm-1)
                  nt=0
                  do jj=1,21
                     k=kk(nt+1)
                     l=kk(nt+2)
                     m=kk(nt+3)
                     n=kk(nt+4)
                     nt=nt+4
                     term=0.d0
                     do i=1,3
                        term=term+(coefa-1.d0)*al(i)**(coefa-2.d0)
     &                       *dldc(k,l,i)*dldc(m,n,i)
     &                       +al(i)**(coefa-1.d0)*d2ldc2(k,l,m,n,i)
                     enddo
                     d2udc2(k,l,m,n)=d2udc2(k,l,m,n)
     &                    +2.d0*coefm/
     &                    coefa*(term+(coefa*coefb+1.d0)*v3b
     &                    **(-coefa*coefb-2.d0)*dibdc(k,l,3)
     &                    *dibdc(m,n,3)-v3b**(-coefa*coefb-1.d0)
     &                    *d2ibdc2(k,l,m,n,3))
                  enddo
               enddo
            endif
         endif
!
!        storing the stiffness matrix and/or the stress
!
         if(icmd.ne.3) then
!
!           storing the stiffness matrix
!
            nt=0
            do i=1,21
               k=kk(nt+1)
               l=kk(nt+2)
               m=kk(nt+3)
               n=kk(nt+4)
               nt=nt+4
               elas(i)=4.d0*d2udc2(k,l,m,n)
            enddo
         endif
!
!           store the stress at mechanical strain
!
         stre(1)=2.d0*dudc(1,1)
         stre(2)=2.d0*dudc(2,2)
         stre(3)=2.d0*dudc(3,3)
         stre(4)=2.d0*dudc(1,2)
         stre(5)=2.d0*dudc(1,3)
         stre(6)=2.d0*dudc(2,3)
!               
      enddo
!
      return
      end
