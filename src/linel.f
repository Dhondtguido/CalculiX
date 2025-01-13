!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine linel(kode,mattyp,beta,emec,stre,stiff,elconloc,
     &     iorien,orab,pgauss,ncmat_,nalcon,imat)
!     
!     calculates stresses for linear elastic materials
!     
      implicit none
!     
      integer mattyp,j1,j2,j3,j4,j5,j6,j7,j8,j,jj,kel(4,21),
     &     iorien,i,kode,ncmat_,nalcon(2,*),imat
!     
      real*8 beta(6),stiff(21),stre(6),fxx,fyy,fzz,fxy,fxz,fyz,
     &     elconloc(*),emax,ya(3,3,3,3),orab(7,*),skl(3,3),e,un,
     &     um,um2,al,am1,pgauss(3),emec(6)
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
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
!     isotropic
!     
        e=elconloc(1)
        un=elconloc(2)
!     
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
!     anisotropic expansion makes the tangent matrix anisotropic
!     
        if(nalcon(1,imat).le.1) then
!     
!     isotropic expansion
!     
          stiff(1)=elconloc(1)
          stiff(2)=elconloc(2)
          mattyp=1
        elseif(nalcon(1,imat).eq.3) then
!     
!     orthotropic expansion
!     
          stiff(1)=am1
          stiff(2)=al
          stiff(3)=am1
          stiff(4)=al
          stiff(5)=al
          stiff(6)=am1
          stiff(7)=um2
          stiff(8)=um2
          stiff(9)=um2
          mattyp=2
        else
!     
!     anorthotropic expansion
!     
          do i=7,9
            stiff(i)=0.d0
          enddo
          do i=11,14
            stiff(i)=0.d0
          enddo
          do i=16,20
            stiff(i)=0.d0
          enddo
          stiff(1)=am1
          stiff(2)=al
          stiff(3)=am1
          stiff(4)=al
          stiff(5)=al
          stiff(6)=am1
          stiff(10)=um2
          stiff(15)=um2
          stiff(21)=um2
          mattyp=3
        endif
!     
      elseif((kode.eq.9).or.(kode.eq.21)) then
!     
        if((kode.eq.9).and.(iorien.eq.0)) then
!     
!         mechanically orthotropic
!     
          stre(1)=stiff(1)*fxx+stiff(2)*fyy+
     &         stiff(4)*fzz-beta(1)
          stre(2)=stiff(2)*fxx+stiff(3)*fyy+
     &         stiff(5)*fzz-beta(2)
          stre(3)=stiff(4)*fxx+stiff(5)*fyy+
     &         stiff(6)*fzz-beta(3)
          stre(4)=stiff(7)*fxy-beta(4)
          stre(5)=stiff(8)*fxz-beta(5)
          stre(6)=stiff(9)*fyz-beta(6)
!
          if(nalcon(1,imat).le.6) then
!
!           isotropic or orthotropic expansion
!
            do i=1,9
              stiff(i)=elconloc(i)
            enddo
            do i=10,21
              stiff(i)=0.d0
            enddo
            mattyp=2
          else
!
!           anorthotropic expansion
!
            do i=1,6
              stiff(i)=elconloc(i)
            enddo
            do i=7,9
              stiff(i)=0.d0
            enddo
            stiff(10)=elconloc(7)
            do i=11,14
              stiff(i)=0.d0
            enddo
            stiff(15)=elconloc(8)
            do i=16,20
              stiff(i)=0.d0
            enddo
            stiff(21)=elconloc(9)
            mattyp=3
          endif
!     
        else
!     
          do i=1,ncmat_
            stiff(i)=elconloc(i)
          enddo
!     
          mattyp=3
!     
          if(iorien.ne.0) then
!     
!     calculating the transformation matrix
!     
            call transformatrix(orab(1,iorien),pgauss,skl)
!     
!     transforming the elastic coefficients
!     
            if(kode.eq.9) then
              call orthotropic(stiff,ya)
            else
              call anisotropic(stiff,ya)
            endif
!     
            do jj=1,21
              j1=kel(1,jj)
              j2=kel(2,jj)
              j3=kel(3,jj)
              j4=kel(4,jj)
              stiff(jj)=0.d0
              do j5=1,3
                do j6=1,3
                  do j7=1,3
                    do j8=1,3
                      stiff(jj)=stiff(jj)+ya(j5,j6,j7,j8)*
     &                     skl(j1,j5)*skl(j2,j6)*skl(j3,j7)*
     &                     skl(j4,j8)
                    enddo
                  enddo
                enddo
              enddo
            enddo
!     
!     determining the type: orthotropic or anisotropic
!
            if(nalcon(1,imat).le.3) then
!     
!             at most orthotropic expansion
!
              emax=0.d0
              do j=1,21
                emax=max(emax,dabs(stiff(j)))
              enddo
              do j=7,9
                if(dabs(stiff(j)).gt.emax*1.d-10) then
                  emax=-1.d0
                  exit
                endif
              enddo
              if(emax.ge.0.d0) then
                do j=11,14
                  if(dabs(stiff(j)).gt.emax*1.d-10) then
                    emax=-1.d0
                    exit
                  endif
                enddo
              endif
              if(emax.ge.0.d0) then
                do j=16,20
                  if(dabs(stiff(j)).gt.emax*1.d-10) then
                    emax=-1.d0
                    exit
                  endif
                enddo
              endif
            else
!     
!             anorthotropic expansion
!
              emax=-1.d0
            endif
!
            if(emax.ge.0.d0) then
              stiff(7)=stiff(10)
              stiff(8)=stiff(15)
              stiff(9)=stiff(21)
!     
              do j=10,21
                stiff(j)=0.d0
              enddo
!     
              mattyp=2
            endif
          endif
!     
          if(mattyp.eq.2) then
!     
!     orthotropic
!     
            stre(1)=stiff(1)*fxx+stiff(2)*fyy+
     &           stiff(4)*fzz-beta(1)
            stre(2)=stiff(2)*fxx+stiff(3)*fyy+
     &           stiff(5)*fzz-beta(2)
            stre(3)=stiff(4)*fxx+stiff(5)*fyy+
     &           stiff(6)*fzz-beta(3)
            stre(4)=stiff(7)*fxy-beta(4)
            stre(5)=stiff(8)*fxz-beta(5)
            stre(6)=stiff(9)*fyz-beta(6)
          else
!     
!     fully anisotropic
!     
            stre(1)=stiff(1)*fxx+
     &           stiff(2)*fyy+
     &           stiff(4)*fzz+
     &           stiff(7)*fxy+
     &           stiff(11)*fxz+
     &           stiff(16)*fyz-beta(1)
            stre(2)=stiff(2)*fxx+
     &           stiff(3)*fyy+
     &           stiff(5)*fzz+
     &           stiff(8)*fxy+
     &           stiff(12)*fxz+
     &           stiff(17)*fyz-beta(2)
            stre(3)=stiff(4)*fxx+
     &           stiff(5)*fyy+
     &           stiff(6)*fzz+
     &           stiff(9)*fxy+
     &           stiff(13)*fxz+
     &           stiff(18)*fyz-beta(3)
            stre(4)=stiff(7)*fxx+
     &           stiff(8)*fyy+
     &           stiff(9)*fzz+
     &           stiff(10)*fxy+
     &           stiff(14)*fxz+
     &           stiff(19)*fyz-beta(4)
            stre(5)=stiff(11)*fxx+
     &           stiff(12)*fyy+
     &           stiff(13)*fzz+
     &           stiff(14)*fxy+
     &           stiff(15)*fxz+
     &           stiff(20)*fyz-beta(5)
            stre(6)=stiff(16)*fxx+
     &           stiff(17)*fyy+
     &           stiff(18)*fzz+
     &           stiff(19)*fxy+
     &           stiff(20)*fxz+
     &           stiff(21)*fyz-beta(6)
!     
          endif
        endif
      endif
!     
      return
      end
      
