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
      subroutine filters(inpc,textpart,istep,istat,n,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc,objectset,ier,nobject,nmethod)        
!
!     reading the input deck: *FILTER
!
!     options: TYPE
!              BOUNDARY WEIGHTING
!              EDGE PRESERVATION
!              DIRECTION WEIGHTING
!            
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 objectset(5,*)
!
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),nmethod,
     &  inp(3,*),ipoinpc(0:*),ipos,boundact,filteract,ier,nobject  
!
      real*8 radius
!     
c      if(istep.lt.1) then
      if(nmethod.ne.12) then
        write(*,*) '*ERROR reading *FILTER: *FILTER can
     &only be used within a SENSITIVITY STEP'     
        ier=1
        return
      endif
!     
      if(nobject.eq.0) then
        write(*,*) '*ERROR reading *FILTER: at least one'
        write(*,*) '       *DESIGN RESPONSE must have been'
        write(*,*) '       defined before the definition of'
        write(*,*) '       a filter'
        ier=1
        return
      endif
!     
      boundact=0
      filteract=0
!
      do i=2,n
!  
!        reading filter options:
!     
!        type of filter  
!      
         if(textpart(i)(1:5).eq.'TYPE=') then
            objectset(2,1)(1:5)=textpart(i)(6:10)
            filteract=1      
!     
!        boundary weighting activated
!
         elseif(textpart(i)(1:18).eq.'BOUNDARYWEIGHTING=') then
            if(textpart(i)(19:21).eq.'YES') then
               boundact=1
               objectset(2,1)(6:8)='BOU'
            else         
               boundact=0
            endif
!     
!        edge weighting activated
!
         elseif(textpart(i)(1:17).eq.'EDGEPRESERVATION=') then
            if(textpart(i)(18:20).eq.'YES') then
               objectset(2,1)(10:12)='EDG'
            endif
!     
!        direction weighting activated
!
         elseif(textpart(i)(1:19).eq.'DIRECTIONWEIGHTING=') then
            if(textpart(i)(20:22).eq.'YES') then
               objectset(2,1)(14:16)='DIR'
            endif
!     
!        parameter not recognized
!
         else
            write(*,*) 
     &        '*WARNING reading *FILTER: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &           "*FILTER%")
         endif
      enddo 
!
!     reading the radii
!     
      if((filteract.eq.1).or.(boundact.eq.1)) then
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!     
!     reading in the filter radius
!
         if(filteract.eq.1) then
            read(textpart(1)(1:20),'(f20.0)',iostat=istat) radius
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &                       "*FILTER%",ier)
               return
            endif
            if(radius.lt.0.d0) then
               write(*,*) '*ERROR reading *FILTER'
               write(*,*) '       Radius of the sensitivity'
               write(*,*) '       filter cannot be less than 0' 
               write(*,*)  
               call inputerror(inpc,ipoinpc,iline,
     &                           "*FILTER%",ier)
               return
            endif
            objectset(2,1)(21:40)=textpart(1)(1:20)
         endif
!
!     reading in the radius for boundary weighting
!
         if((n.eq.2).and.(boundact.eq.1)) then
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) radius
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &                       "*FILTER%",ier)
               return
            endif
            if(radius.lt.0.d0) then
               write(*,*) '*ERROR reading *FILTER'
               write(*,*) '       Radius for the boundary'
               write(*,*) '       weighting cannot be less' 
               write(*,*) '       than 0'
               write(*,*)   
               call inputerror(inpc,ipoinpc,iline,
     &                            "*FILTER%",ier)
               return
            endif
            objectset(1,1)(21:40)=textpart(2)(1:20)
         elseif((n.eq.1).and.(boundact.eq.1)) then
            write(*,*) '*WARNING reading *FILTER:'
            write(*,*) '         boundary weighting activated'
            write(*,*) '         but no radius defined'
            write(*,*) '         The radius of the sensitivity' 
            write(*,*) '         filter will be taken'
            write(*,*)   
            call inputwarning(inpc,ipoinpc,iline,
     &           "*FILTER%")
            objectset(1,1)(21:40)=objectset(2,1)(21:40)
         endif
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
      else
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      endif
!     
      return
      end
      
      
