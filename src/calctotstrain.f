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
      subroutine calctotstrain(vkl,vokl,eloc,elineng,iperturb)
!
!     calculates the total strain from the displacement gradients
!
      implicit none
!
      integer iperturb(*)
!
      real*8 elineng(6),vkl(0:3,3),vokl(3,3),eloc(6)
!
!     calculating the strain
!
!     attention! elineng(4),elineng(5) and elineng(6) are engineering strains!
!     
      elineng(1)=vkl(1,1)
      elineng(2)=vkl(2,2)
      elineng(3)=vkl(3,3)
      elineng(4)=vkl(1,2)+vkl(2,1)
      elineng(5)=vkl(1,3)+vkl(3,1)
      elineng(6)=vkl(2,3)+vkl(3,2)
!     
      if(iperturb(2).eq.1) then
!     
!     Lagrangian strain
!     
         elineng(1)=elineng(1)+
     &        (vkl(1,1)**2+vkl(2,1)**2+vkl(3,1)**2)/2.d0
         elineng(2)=elineng(2)+
     &        (vkl(1,2)**2+vkl(2,2)**2+vkl(3,2)**2)/2.d0
         elineng(3)=elineng(3)+
     &        (vkl(1,3)**2+vkl(2,3)**2+vkl(3,3)**2)/2.d0
         elineng(4)=elineng(4)+vkl(1,1)*vkl(1,2)+vkl(2,1)*vkl(2,2)+
     &        vkl(3,1)*vkl(3,2)
         elineng(5)=elineng(5)+vkl(1,1)*vkl(1,3)+vkl(2,1)*vkl(2,3)+
     &        vkl(3,1)*vkl(3,3)
         elineng(6)=elineng(6)+vkl(1,2)*vkl(1,3)+vkl(2,2)*vkl(2,3)+
     &        vkl(3,2)*vkl(3,3)
!     
!     for frequency analysis or buckling with preload the
!     strains are calculated with respect to the deformed
!     configuration
!     
      elseif(iperturb(1).eq.1) then
         elineng(1)=elineng(1)+vokl(1,1)*vkl(1,1)+vokl(2,1)*vkl(2,1)+
     &        vokl(3,1)*vkl(3,1)
         elineng(2)=elineng(2)+vokl(1,2)*vkl(1,2)+vokl(2,2)*vkl(2,2)+
     &        vokl(3,2)*vkl(3,2)
         elineng(3)=elineng(3)+vokl(1,3)*vkl(1,3)+vokl(2,3)*vkl(2,3)+
     &        vokl(3,3)*vkl(3,3)
         elineng(4)=elineng(4)+vokl(1,1)*vkl(1,2)+vokl(1,2)*vkl(1,1)+
     &        vokl(2,1)*vkl(2,2)+vokl(2,2)*vkl(2,1)+
     &        vokl(3,1)*vkl(3,2)+vokl(3,2)*vkl(3,1)
         elineng(5)=elineng(5)+vokl(1,1)*vkl(1,3)+vokl(1,3)*vkl(1,1)+
     &        vokl(2,1)*vkl(2,3)+vokl(2,3)*vkl(2,1)+
     &        vokl(3,1)*vkl(3,3)+vokl(3,3)*vkl(3,1)
         elineng(6)=elineng(6)+vokl(1,2)*vkl(1,3)+vokl(1,3)*vkl(1,2)+
     &        vokl(2,2)*vkl(2,3)+vokl(2,3)*vkl(2,2)+
     &        vokl(3,2)*vkl(3,3)+vokl(3,3)*vkl(3,2)
      endif
!     
!     storing the local strains
!     
      if(iperturb(1).ne.-1) then
         eloc(1)=elineng(1)
         eloc(2)=elineng(2)
         eloc(3)=elineng(3)
         eloc(4)=elineng(4)/2.d0
         eloc(5)=elineng(5)/2.d0
         eloc(6)=elineng(6)/2.d0
      else
!     
!     linear iteration within a nonlinear increment:
!     
         eloc(1)=vokl(1,1)+
     &        (vokl(1,1)**2+vokl(2,1)**2+vokl(3,1)**2)/2.d0
         eloc(2)=vokl(2,2)+
     &        (vokl(1,2)**2+vokl(2,2)**2+vokl(3,2)**2)/2.d0
         eloc(3)=vokl(3,3)+
     &        (vokl(1,3)**2+vokl(2,3)**2+vokl(3,3)**2)/2.d0
         eloc(4)=(vokl(1,2)+vokl(2,1)+vokl(1,1)*vokl(1,2)+
     &        vokl(2,1)*vokl(2,2)+vokl(3,1)*vokl(3,2))/2.d0
         eloc(5)=(vokl(1,3)+vokl(3,1)+vokl(1,1)*vokl(1,3)+
     &        vokl(2,1)*vokl(2,3)+vokl(3,1)*vokl(3,3))/2.d0
         eloc(6)=(vokl(2,3)+vokl(3,2)+vokl(1,2)*vokl(1,3)+
     &        vokl(2,2)*vokl(2,3)+vokl(3,2)*vokl(3,3))/2.d0
      endif
!     
      return
      end
