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
      subroutine feasibledirections(inpc,textpart,istat,n,key,iline,
     &     ipol,inl,ipoinp,inp,ipoinpc,nmethod,istep,ier,tmax)          
!     
!     reading the input deck: *FEASIBLE DIRECTION
!     
      implicit none
!     
      character*1 inpc(*)
      character*132 textpart(16)
!     
      integer nmethod,nobject,i,n,key,istat,istep,iline,ipol,inl,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier
!     
      real*8 tmax
!     
      tmax=0.d0
!
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *FEASIBLE DIRECTION:'
        write(*,*) '       *FEASIBLE DIRECTION can only be used'
        write(*,*) '       within a STEP'    
        ier=1
        return
      endif
!     
      nmethod=16 
!     
!     read optimization method
!     tmax=1.5: Gradient Descent Akin Method (GDAM, default)
!     tmax=2.5: Gradient Projection Method (GPM)
!     
      do i=2,n
         if(textpart(i)(1:7).eq.'METHOD=') then
            if(textpart(i)(8:22).eq.'GRADIENTDESCENT') then
               tmax=1.5
            elseif(textpart(i)(8:25).eq.'GRADIENTPROJECTION') then
               tmax=2.5
            else
               write(*,*) 
               write(*,*) '*WARNING reading *FEASIBLE DIRECTION; '     
               write(*,*) '         Method for computation of '     
               write(*,*) '         *FEASIBLE DIRECTION not valid;'
               write(*,*) '         Gradient Descent taken as default' 
               write(*,*) ' '
               call inputwarning(inpc,ipoinpc,iline,
     &           "*FEASIBLEDIRECTION%")
            endif
         endif
      enddo
!
      if(tmax.lt.1.d0) then
         tmax=1.5
         write(*,*)
         write(*,*) '*WARNING reading *FEASIBLE DIRECTION; '     
         write(*,*) '         Method for computation of '     
         write(*,*) '         *FEASIBLE DIRECTION not specified;'
         write(*,*) '         Gradient Descent taken as default'
         write(*,*) ' '
         call inputwarning(inpc,ipoinpc,iline,"*FEASIBLEDIRECTION%")
      endif          
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)     
!      
      return
      end      
      
