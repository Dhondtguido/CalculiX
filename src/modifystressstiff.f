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
      subroutine modifystressstiff(stre,stiff,mattyp,eth,nalcon,imat,
     &     xthi,vthj)
!
      implicit none
!
      integer mattyp,i,nalcon(2,*),imat,kal(2,6),kel(4,21),jj,j1,j2,
     &     j3,j4,j5,j6,j7,j8
!
      real*8 stre(6),stiff(21),eth(6),fth1,fth2,fth3,skl(3,3),vthj,
     &     xthi(3,3),ya(3,3,3,3)
!     
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!
      if(mattyp.eq.1) then
!
!       elastically isotropic and isotropic expansion
!
        fth1=1.d0+eth(1)
        do i=1,6
          stre(i)=stre(i)*fth1
        enddo
        stiff(1)=stiff(1)/fth1
      elseif(mattyp.eq.2) then
!
!       elastically orthotropic or orthotropic expansion or both 
!
        if(nalcon(1,imat).eq.1) then
          fth1=1.d0+eth(1)
          do i=1,6
            stre(i)=stre(i)*fth1
          enddo
          do i=1,9
            stiff(i)=stiff(i)/fth1
          enddo
        else
          fth1=1.d0+eth(1)
          fth2=1.d0+eth(2)
          fth3=1.d0+eth(3)
          stre(1)=stre(1)*fth2*fth3/fth1
          stre(2)=stre(2)*fth1*fth3/fth2
          stre(3)=stre(3)*fth1*fth2/fth3
          stre(4)=stre(4)*fth3
          stre(5)=stre(5)*fth2
          stre(6)=stre(6)*fth1
          stiff(1)=stiff(1)*fth2*fth3/fth1**3
          stiff(2)=stiff(2)*fth3/(fth1*fth2)
          stiff(3)=stiff(3)*fth1*fth3/fth2**3
          stiff(4)=stiff(4)*fth2/(fth1*fth3)
          stiff(5)=stiff(5)*fth1/(fth2*fth3)
          stiff(6)=stiff(6)*fth1*fth2/fth3**3
          stiff(7)=stiff(7)*fth3/(fth1*fth2)
          stiff(8)=stiff(8)*fth2/(fth1*fth3)
          stiff(9)=stiff(9)*fth1/(fth2*fth3)
        endif
      else
!
!       elastically anisotropic or anisotropic expansion or both 
!
        if(nalcon(1,imat).eq.1) then
          fth1=1.d0+eth(1)
          do i=1,6
            stre(i)=stre(i)*fth1
          enddo
          do i=1,21
            stiff(i)=stiff(i)/fth1
          enddo
        elseif(nalcon(1,imat).eq.3) then
          fth1=1.d0+eth(1)
          fth2=1.d0+eth(2)
          fth3=1.d0+eth(3)
          stre(1)=stre(1)*fth2*fth3/fth1
          stre(2)=stre(2)*fth1*fth3/fth2
          stre(3)=stre(3)*fth1*fth2/fth3
          stre(4)=stre(4)*fth3
          stre(5)=stre(5)*fth2
          stre(6)=stre(6)*fth1
          stiff(1)=stiff(1)*fth2*fth3/fth1**3
          stiff(2)=stiff(2)*fth3/(fth1*fth2)
          stiff(3)=stiff(3)*fth1*fth3/fth2**3
          stiff(4)=stiff(4)*fth2/(fth1*fth3)
          stiff(5)=stiff(5)*fth1/(fth2*fth3)
          stiff(6)=stiff(6)*fth1*fth2/fth3**3
          stiff(7)=stiff(7)*fth3/(fth1*fth1)
          stiff(8)=stiff(8)*fth3/(fth2*fth2)
          stiff(9)=stiff(9)/fth3
          stiff(10)=stiff(10)*fth3/(fth1*fth2)
          stiff(11)=stiff(11)*fth2/(fth1*fth1)
          stiff(12)=stiff(12)/fth2
          stiff(13)=stiff(13)*fth2/(fth3*fth3)
          stiff(14)=stiff(14)/fth1
          stiff(15)=stiff(15)*fth2/(fth1*fth3)
          stiff(16)=stiff(16)/fth1
          stiff(17)=stiff(17)*fth1/(fth2*fth2)
          stiff(18)=stiff(18)*fth1/(fth3*fth3)
          stiff(19)=stiff(19)/fth2
          stiff(20)=stiff(20)/fth3
          stiff(21)=stiff(21)*fth1/(fth2*fth3)
        elseif(nalcon(1,imat).eq.6) then  
          skl(1,1)=stre(1)
          skl(2,2)=stre(2)
          skl(3,3)=stre(3)
          skl(1,2)=stre(4)
          skl(1,3)=stre(5)
          skl(2,3)=stre(6)
          skl(2,1)=stre(4)
          skl(3,1)=stre(5)
          skl(3,2)=stre(6)
!
!         modifying the stress
!          
          do jj=1,6
            stre(jj)=0.d0
            j1=kal(1,jj)
            j2=kal(2,jj)
            do j3=1,3
              do j4=1,3
                stre(jj)=stre(jj)+xthi(j1,j3)*skl(j3,j4)*xthi(j4,j2)
              enddo
            enddo
            stre(jj)=stre(jj)*vthj
          enddo
!
!         modifying the stiffness
!          
          call anisotropic(stiff,ya)
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
                    stiff(jj)=stiff(jj)+xthi(j1,j5)*xthi(j2,j6)
     &                   *xthi(j3,j7)*xthi(j4,j8)*ya(j5,j6,j7,j8)
                  enddo
                enddo
              enddo
            enddo
            stiff(jj)=stiff(jj)*vthj
          enddo
        endif
      endif
!       
      return
      end
