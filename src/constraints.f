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
      subroutine constraints(inpc,textpart,istat,n,
     &     iline,ipol,inl,ipoinp,inp,ipoinpc,nobject,objectset,
     &     ier,nmethod)
!     
!     reading the input deck: *CONSTRAINT
!     
      implicit none
!     
      logical copy
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 objectset(5,*),drname
!     
      integer istat,n,key,i,iline,ipol,inl,ipoinp(2,*),
     &     inp(3,*),ipoinpc(0:*),nobject,ier,nmethod,iobject,
     &     j,icopy
!      
      real*8 relval,absval
!     
!     constraints can only be defined within feasible direction steps.
!     We ensure that's the case below
!     
      if(nmethod.ne.16) then
        write(*,*) '*ERROR reading *CONSTRAINT'
        write(*,*) '       *CONSTRAINT can only be specified'
        write(*,*) '       within a *FEASIBILE DIRECTION step.'   
        call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
        return
      endif
!     
!     one can define multiple constraints under a single *CONSTRAINT flag which
!     means that we need to loop over all of them.
!     We do not know how many there are so we run until the line does not
!     contain a key or we reached the end of the file.
!     
      do
!     
!     reading the design response which should be used for this constraint
!     we first get a new line and then check if the given name is not longer
!     than 80 chars and the design response exists
!     
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
!        
        drname=textpart(1)(1:80)
        if((istat.lt.0).or.(key.eq.1)) then        
          return
        else
          if(textpart(1)(81:81).ne.' ') then
            write(*,*) '*ERROR reading *CONSTRAINT'
            write(*,*) '       name of design response must not '
            write(*,*) '       be longer than 80 chars.'
            call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
            return
          endif
        endif
!     
!     check if design response exists and store the id in iobject;
!     if the design response is already used for
!     either a (geometric) constraint or an objective    
!     create a new entry in the objectset with the same name.
!     
        iobject=0
        copy=.false.
        icopy=0
!
        do i=1,nobject
          if(objectset(5,i)(1:80).eq.drname) then
            if(objectset(5,i)(81:81).ne.' ') then
              icopy=i
              copy=.true.
            else
              iobject=i
              exit
            endif
          endif
        enddo
!     
!     no empty entry with the same name has been found
!     -> create a new entry and increment nobject
!     
        if(copy) then
          nobject=nobject+1
          iobject=nobject
          do j=1,5
            objectset(j,iobject)(1:81)=objectset(j,icopy)(1:81)
          enddo
        endif
!     
!     if iobject=0, there is no entry with the same name
!     
        if(iobject.eq.0) then
          write(*,*) '*ERROR reading *CONSTRAINT'
          write(*,*) '       name of given design'
          write(*,*) '       response does not exist:', drname
          call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
          return
        endif
!     
!     mark this entry as a constraint 
!     
        objectset(5,iobject)(81:81)='C'   
!     
!     check the constraint mode (LE, GE)
!     throw an error if either its not given or not recognised.
!     
        if(n.ge.2) then
          if(textpart(2)(1:2).eq.'LE') then
            objectset(1,iobject)(19:20)='LE'
          elseif (textpart(2)(1:2).eq.'GE') then
            objectset(1,iobject)(19:20)='GE'
          else
            write(*,*) '*ERROR reading *CONSTRAINT'
            write(*,*) '       mode of constraint must be'
            write(*,*) '       either LE or GE'
            call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
            return
          endif
        else
          write(*,*) '*ERROR reading *CONSTRAINT'
          write(*,*) '       mode of constraint must be'
          write(*,*) '       either LE or GE'
          call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
          return
        endif  
!     
!     The K-S functions only need relative constraints as the absolute one
!     is already incorporated in the function (e.g. as target stress)
!     All the other functions constraints can have both
!     
!        if((objectset(1,iobject)(1:8).eq.'ALL-DISP').or.
!     &       (objectset(1,iobject)(1:6).eq.'X-DISP').or.
!     &       (objectset(1,iobject)(1:6).eq.'Y-DISP').or.
!     &       (objectset(1,iobject)(1:6).eq.'Z-DISP').or.
!     &       (objectset(1,iobject)(1:4).eq.'MASS')) then
!     
!     relative constraint value
!     
          if(n.ge.3) then
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) relval
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
              return
            endif
            if(istat.le.0) then
              objectset(1,iobject)(41:60)=textpart(3)(1:20)
            endif
          endif 
!     
!     absolute constraint value
!     
          if(n.ge.4) then
            read(textpart(4)(1:20),'(f20.0)',iostat=istat) absval
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
              return
            endif
            objectset(1,iobject)(61:80)=textpart(4)(1:20)
          endif
!
!        elseif((objectset(1,iobject)(1:11).eq.'MODALSTRESS').or.
!     &       (objectset(1,iobject)(1:12).eq.'STRAINENERGY').or.
!     &       (objectset(1,iobject)(1:11).eq.'MISESSTRESS').or.
!     &       (objectset(1,iobject)(1:9).eq.'PS1STRESS').or.
!     &       (objectset(1,iobject)(1:9).eq.'PS3STRESS').or.
!     &       (objectset(1,iobject)(1:17).eq.'EQUIVALENTPLASTICSTRAIN')) 
!     &       then
!     
!     relative constraint value
!     
!          if(n.ge.3) then
!            read(textpart(3)(1:20),'(f20.0)',iostat=istat) absval
!            if(istat.gt.0) then
!              call inputerror(inpc,ipoinpc,iline,"*CONSTRAINT%",ier)
!              return
!            endif
!            objectset(1,iobject)(41:60)=textpart(3)(1:20)
!          endif
!        endif  
      enddo
!      
      return
      end
