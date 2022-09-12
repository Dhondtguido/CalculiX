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
      subroutine designvariabless(inpc,textpart,tieset,tietol,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,ipoinpc,
     &       set,nset,ier)
!
!     reading the input deck: *DESIGN VARIABLES
!
      implicit none
!
      character*1 inpc(*)
      character*81 tieset(3,*),set(*)
      character*132 textpart(16)
!
      integer istep,istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ntie,ntie_,ipoinpc(0:*),nset,itype,ier,id
!
      real*8 tietol(4,*)
!
!     Check of correct position in Inputdeck
!
      if(istep.gt.0) then
        write(*,*)
     &       '*ERROR reading *DESIGN VARIABLES: *DESIGN VARIABLES'
         write(*,*) ' should be placed before all step definitions'
         ier=1
         return
      endif
!
!     Check of correct number of ties
!
      ntie=ntie+1
      if(ntie.gt.ntie_) then
         write(*,*) '*ERROR reading *DESIGN VARIABLES: increase ntie_'
         ier=1
         return
      endif
!
!     Read in *DESIGN VARIABLES
!
      itype=0
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            read(textpart(i)(6:85),'(a80)',iostat=istat) 
     &           tieset(1,ntie)(1:80)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DESIGN VARIABLES%",ier)
               return
            endif
            itype=1
        else
          write(*,*) 
     &   '*WARNING reading *DESIGN VARIABLES: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*DESIGN VARIABLES%")
         endif
       enddo  
!
      if(itype.eq.0) then
         write(*,*) 
     &'*ERROR reading *DESIGN VARIABLES: type is lacking'
         call inputerror(inpc,ipoinpc,iline,
     &        "*DESIGN VARIABLES%",ier)
         return
      endif
!
!     Add "D" at the end of the name of the designvariable keyword
!      
      tieset(1,ntie)(81:81)='D' 
!
      if(tieset(1,ntie)(1:10).eq.'COORDINATE') then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            write(*,*)
     &'*ERROR reading *DESIGN VARIABLES: definition'
            write(*,*) '      is not complete.'
            ier=1
            return
         endif
!
!        Read the name of the design variable node set
!
         tieset(2,ntie)(1:81)=textpart(1)(1:81)
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='N'
!
!        Check existence of the node set
!
         call cident81(set,tieset(2,ntie),nset,id)
         i=nset+1
         if(id.gt.0) then
           if(tieset(2,ntie).eq.set(id)) then
             i=id
           endif
         endif
         if(i.gt.nset) then
            write(*,*) '*ERROR reading *DESIGN VARIABLES'
            write(*,*) 'node set ',tieset(2,ntie)(1:ipos-1),
     &           'does not exist. Card image:'
            call inputerror(inpc,ipoinpc,iline,
     &           "*DESIGN VARIABLES%",ier)
            return
         endif
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end



