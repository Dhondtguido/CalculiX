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
      subroutine linel(kode,mattyp,beta,emec,stre,elas,elconloc,
     &  iorien,orab,pgauss)
!
!     calculates stresses for linear elastic materials
!
      implicit none
!
      integer mattyp,j1,j2,j3,j4,j5,j6,j7,j8,j,jj,kel(4,21),
     &  iorien,i,kode
!
      real*8 beta(6),elas(21),stre(6),fxx,fyy,fzz,fxy,fxz,fyz,
     &  elconloc(*),emax,ya(3,3,3,3),orab(7,*),skl(3,3),e,un,
     &  um,um2,al,am1,pgauss(3),emec(6)
!
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &          1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &          3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &          1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!
!     engineering strain
!
      fxx=emec(1)
      fyy=emec(2)
      fzz=emec(3)
      fxy=2.d0*emec(4)
      fxz=2.d0*emec(5)
      fyz=2.d0*emec(6)
!
      if(kode.eq.2) then
!
!        isotropic
!
         do i=1,2
            elas(i)=elconloc(i)
         enddo
!
         e=elas(1)
         un=elas(2)
         um2=e/(1.d0+un)
         al=un*um2/(1.d0-2.d0*un)
         um=um2/2.d0
         am1=al+um2
!
         stre(1)=am1*fxx+al*(fyy+fzz)-beta(1)
         stre(2)=am1*fyy+al*(fxx+fzz)-beta(2)
         stre(3)=am1*fzz+al*(fxx+fyy)-beta(3)
         stre(4)=um*fxy-beta(4)
         stre(5)=um*fxz-beta(5)
         stre(6)=um*fyz-beta(6)
!
         mattyp=1
!
      elseif((kode.eq.9).or.(kode.eq.21)) then
!
         if((kode.eq.9).and.(iorien.eq.0)) then
!        
!           orthotropic
!
            do i=1,9
               elas(i)=elconloc(i)
            enddo
            do i=10,21
               elas(i)=0.d0
            enddo
!
            stre(1)=elas(1)*fxx+elas(2)*fyy+
     &           elas(4)*fzz-beta(1)
            stre(2)=elas(2)*fxx+elas(3)*fyy+
     &           elas(5)*fzz-beta(2)
            stre(3)=elas(4)*fxx+elas(5)*fyy+
     &           elas(6)*fzz-beta(3)
            stre(4)=elas(7)*fxy-beta(4)
            stre(5)=elas(8)*fxz-beta(5)
            stre(6)=elas(9)*fyz-beta(6)
!
            mattyp=2
!
         else
!
            do i=1,21
               elas(i)=elconloc(i)
            enddo
!
            mattyp=3
!
            if(iorien.ne.0) then
!
!              calculating the transformation matrix
!
               call transformatrix(orab(1,iorien),pgauss,skl)
!
!              transforming the elastic coefficients
!
               if(kode.eq.9) then
                  call orthotropic(elas,ya)
               else
                  call anisotropic(elas,ya)
               endif
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
     &                             skl(j1,j5)*skl(j2,j6)*skl(j3,j7)*
     &                             skl(j4,j8)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
!
!              determining the type: orthotropic or anisotropic
!
               emax=0.d0
               do j=1,21
                  emax=max(emax,dabs(elas(j)))
               enddo
               do j=7,9
                  if(dabs(elas(j)).gt.emax*1.d-10) then
                     emax=-1.d0
                     exit
                  endif
               enddo
               if(emax.ge.0.d0) then
                  do j=11,14
                     if(dabs(elas(j)).gt.emax*1.d-10) then
                        emax=-1.d0
                        exit
                     endif
                  enddo
               endif
               if(emax.ge.0.d0) then
                  do j=16,20
                     if(dabs(elas(j)).gt.emax*1.d-10) then
                        emax=-1.d0
                        exit
                     endif
                  enddo
               endif
               if(emax.ge.0.d0) then
                  elas(7)=elas(10)
                  elas(8)=elas(15)
                  elas(9)=elas(21)
!
                  do j=10,21
                     elas(j)=0.d0
                  enddo
c                  elas(10)=0.d0
c                  elas(15)=0.d0
c                  elas(21)=0.d0
!
                  mattyp=2
               endif
            endif
!
            if(mattyp.eq.2) then
!
!              orthotropic
!
               stre(1)=elas(1)*fxx+elas(2)*fyy+
     &              elas(4)*fzz-beta(1)
               stre(2)=elas(2)*fxx+elas(3)*fyy+
     &              elas(5)*fzz-beta(2)
               stre(3)=elas(4)*fxx+elas(5)*fyy+
     &              elas(6)*fzz-beta(3)
               stre(4)=elas(7)*fxy-beta(4)
               stre(5)=elas(8)*fxz-beta(5)
               stre(6)=elas(9)*fyz-beta(6)
            else
!
!              fully anisotropic
!
               stre(1)=elas(1)*fxx+
     &              elas(2)*fyy+
     &              elas(4)*fzz+
     &              elas(7)*fxy+
     &              elas(11)*fxz+
     &              elas(16)*fyz-beta(1)
               stre(2)=elas(2)*fxx+
     &              elas(3)*fyy+
     &              elas(5)*fzz+
     &              elas(8)*fxy+
     &              elas(12)*fxz+
     &              elas(17)*fyz-beta(2)
               stre(3)=elas(4)*fxx+
     &              elas(5)*fyy+
     &              elas(6)*fzz+
     &              elas(9)*fxy+
     &              elas(13)*fxz+
     &              elas(18)*fyz-beta(3)
               stre(4)=elas(7)*fxx+
     &              elas(8)*fyy+
     &              elas(9)*fzz+
     &              elas(10)*fxy+
     &              elas(14)*fxz+
     &              elas(19)*fyz-beta(4)
               stre(5)=elas(11)*fxx+
     &              elas(12)*fyy+
     &              elas(13)*fzz+
     &              elas(14)*fxy+
     &              elas(15)*fxz+
     &              elas(20)*fyz-beta(5)
               stre(6)=elas(16)*fxx+
     &              elas(17)*fyy+
     &              elas(18)*fzz+
     &              elas(19)*fxy+
     &              elas(20)*fxz+
     &              elas(21)*fyz-beta(6)
!     
            endif
         endif
      endif
!     
      return
      end
      
