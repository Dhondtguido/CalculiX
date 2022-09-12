!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine detectactivecont(gapnorm,gapdisp,auw,iroww,jqw,
     &     nslavs,springarea,iacti,nacti)
!     
!     computing g_Npre=g0+Wb^T*gapdisp
!     
      implicit none
!     
      integer i,inorm,j,icol,iroww(*),jqw(*),nslavs,iacti(*),nacti
!     
      real*8 gapnorm(*),gapdisp(*),auw(*),springarea(2,*),value
!
!     premultiply gapdisp with Wb^T (taking only normal directions
!     into account, i.e. nslav entries)
!
      do i=1,nslavs
        inorm=3*(i-1)+1 ! only normal node DOF
        do j=jqw(inorm),jqw(inorm+1)-1
          value=auw(j)
          icol=iroww(j)
          gapnorm(i)=gapnorm(i)+value*gapdisp(icol) ! TODOCMT friction check!
        enddo
      enddo
!     
      nacti=0
      do i=1,nslavs
!
      j = 3*(i-1)
!     
!     contact evaluation: active degrees of freedom are those for
!     which there is overlap (with added initial clearance at time 0)
!     and for the NON-zero columns (this are nodes which have no master face).
!     
      if((gapnorm(i)+springarea(2,i).le.0.d0).and.
     &     (jqw(j+1).ne.jqw(j+2)) )then
!
        if(jqw(j+1).eq.jqw(j+2))then
          write(*,*) 'Zero column detected!!! Singular contact matrix'
          stop
        endif
!
!     identifying the indices only of the active normals.
!
            iacti(j+1)=nacti+1
            iacti(j+2)=nacti+2
            iacti(j+3)=nacti+3
            nacti=nacti+3
        endif
      enddo
!     
      return
      end
