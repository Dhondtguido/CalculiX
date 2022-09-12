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
      subroutine objectives(inpc,textpart,istat,n,iline,
     &   ipol,inl,ipoinp,inp,ipoinpc,nobject,objectset,
     &   ier,nmethod,objective_flag)        
!
!     reading the input deck: *OBJECTIVE 
!
      implicit none
!
      logical copy,objective_flag
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 objectset(5,*)
!
      integer istat,n,key,i,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*),nobject,ier,nmethod,iobject,j,icopy   
!
      if(nmethod.ne.16) then
         write(*,*) '*ERROR reading *OBJECTIVE'
         write(*,*) '       *OBJECTIVE can only be defined'
         write(*,*) '       within a *FEASIBLE DIRECTION STEP'
         call inputerror(inpc,ipoinpc,iline,"*OBJECTIVE%",ier)
         return
       endif
!
       if(objective_flag) then
         write(*,*) '*ERROR reading *OBJECTIVE'
         write(*,*) '       *OBJECTIVE can only be defined'
         write(*,*) '       once within a *FEASIBLE DIRECTION STEP'
         call inputerror(inpc,ipoinpc,iline,"*OBJECTIVE%",ier)
         return
       endif
!
!     check if it is a minimization or maximization problem
!        
      if(textpart(2)(1:7).eq.'TARGET=') then
         if(textpart(2)(8:10).eq.'MIN') then
            objectset(2,1)(17:19)='MIN'
         elseif(textpart(2)(8:10).eq.'MAX') then 
            objectset(2,1)(17:19)='MAX'  
         else
            write(*,*) '*WARNING optimization TARGET not specified.'
            write(*,*) '         Minimization problem assumed as'
            write(*,*) '         default.'
            objectset(2,1)(17:19)='MIN' 
         endif              
      else
         write(*,*) '*WARNING optimization TARGET not specified.'
         write(*,*) '         Minimization problem assumed as'
         write(*,*) '         default.'
         objectset(2,1)(17:19)='MIN' 
      endif        
!
!     reading the design response which should be used for this objective
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      if((textpart(1)(1:1).eq.'*').or.(istat.lt.0).or.(key.eq.1)) then        
         write(*,*) '*ERROR reading *OBJECTIVE'
         write(*,*) '       no design response specified'
         call inputerror(inpc,ipoinpc,iline,"*OBJECTIVE%",ier)
         return
      endif     
!      
!     check if design response exists. if dr is already used, create a new
!     entry   
!
      iobject=0
      copy=.false.
      icopy=0
      do i=1,nobject
         if(objectset(5,i)(1:80).eq.textpart(1)(1:80)) then
            if(objectset(5,i)(81:81).ne.' ') then
               icopy=i
               copy=.true.
            else
               copy=.false.
               iobject=i
               exit
            endif
         endif
      enddo
!      
      if(copy) then
         nobject=nobject+1
         iobject=nobject
         do j=1,5
            objectset(j,iobject)(1:81)=objectset(j,icopy)(1:81)
         enddo
      endif
!     
      if(iobject.eq.0) then
         write(*,*) '*ERROR reading *OBJECTIVE'
         write(*,*) '       given name of design '
         write(*,*) '       response does not exist.'
         call inputerror(inpc,ipoinpc,iline,"*OBJECTIVE%",ier)
         return;
      endif
!      
      objectset(5,iobject)(81:81)='O'
      objectset(1,iobject)(19:20)='  '
!                 
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!      
      return
      end
