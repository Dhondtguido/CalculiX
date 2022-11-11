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
      subroutine designresponses(inpc,textpart,istat,n,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc,nobject,objectset,
     &     set,nset,ntie,tieset,ier,nmethod)        
!     
!     reading the input deck: *DESIGN RESPONSE
!     
!     criteria: 
!     ALL-DISP
!     X-DISP
!     Y-DISP
!     Z-DISP
!     EIGENFREQUENCY
!     GREEN
!     MASS
!     MODALSTRESS
!     STRAIN ENERGY
!     MISESSTRESS
!     PS1STRESS
!     PS3STRESS
!     EQUIVALENT PLASTIC STRAIN
!     
      implicit none
!     
      character*1 inpc(*),settype
      character*132 textpart(16)
      character*81 objectset(5,*),set(*),tieset(3,*),desrespname
!     
      integer istat,n,key,i,iline,ipol,inl,ipoinp(2,*),nset,j,id,
     &     inp(3,*),ipoinpc(0:*),nobject,ipos,icoordinate,
     &     ntie,ier,nmethod,nsets,m
!     
      real*8 rho,stress
!     
!     check whether the design variables are the coordinates
!     
      icoordinate=0
      do i=1,ntie
        if(tieset(1,i)(81:81).eq.'D') then
          if(tieset(1,i)(1:10).eq.'COORDINATE') then
            icoordinate=1
            exit
          elseif(tieset(1,i)(1:11).eq.'ORIENTATION') then
            exit
          endif
        endif
      enddo
!     
!     check that this is a sensitivity step
!     
      if(nmethod.ne.12) then
        write(*,*) '*ERROR reading *DESIGN RESPONSE'
        write(*,*) '       *DESIGN RESPONSE can only be '
        write(*,*) '       used within a SENSITIVITY STEP.'     
        call inputerror(inpc,ipoinpc,iline,
     &       "*DESIGN RESPONSE%",ier)
        return
      endif
!     
!     keeping track of the name for the objective to reference
!     it within the optimization step. Names must not be longer than 80
!     characters.
!     
      do i=2,n
        if(textpart(i)(1:5).eq.'NAME=') then
          if(textpart(i)(86:86).ne.' ') then
            write(*,*) '*ERROR in *DESIGN RESPONSE'
            write(*,*) '       reading argument NAME of'
            write(*,*) '       *DESIGN RESPONSE, NAME must'
            write(*,*) '       not be longer than 80 digits.'
            call inputerror(inpc,ipoinpc,iline,
     &           "*DESIGN RESPONSE%",ier)
            return
          endif
!     
          do j=1,nobject
            if(objectset(5,j)(1:80).eq.textpart(i)(6:85)) then
              write(*,*) '*ERROR reading *DESIGN RESPONSE'
              write(*,*) '       NAME has already'
              write(*,*) '       been used for a different'
              write(*,*) '       design response:'
              write(*,*) textpart(i)(6:85)
              call inputerror(inpc,ipoinpc,iline,
     &             "*DESIGN RESPONSE%",ier)
              return
            endif
          enddo       
          desrespname=textpart(2)(6:85)
        else
          write(*,*) 
     &    '*WARNING reading *DESIGN RESPONSE: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*DESIGN RESPONSE%")
        endif
      enddo
!     
!     check that a name was defined for a design response with
!     coordinate design variables
!     
      if(icoordinate.eq.1) then
        if(desrespname(1:1).eq.' ') then
c        if(objectset(5,nobject)(1:1).eq.' ') then
          write(*,*) '*ERROR reading *DESIGN RESPONSE:'
          write(*,*) '       no name given.'
          call inputerror(inpc,ipoinpc,iline,
     &         "*DESIGN RESPONSE%",ier)
          return
        endif
      endif
!     
      do
!     
!     reading the design response type and the set to which it applies
!     
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
!     
        if((istat.lt.0).or.(key.eq.1)) then
          if(icoordinate.eq.1) then
            write(*,*) '*ERROR reading *DESIGN RESPONSE'
            write(*,*) '       no information about sets specified'
            call inputerror(inpc,ipoinpc,iline,
     &           "*DESIGN RESPONSE%",ier)
            return
          endif
          exit
        endif
!     
!     incrementing nobject
!     
        nobject=nobject+1
!     
        if(icoordinate.eq.1) then
          objectset(5,nobject)(1:81)=desrespname
        endif
!     
!     reading type of design response and defining the set
!     on which it acts
!     
        nsets=1
        settype=' '
        if(textpart(1)(1:8).eq.'ALL-DISP') then 
          settype='N'
          objectset(1,nobject)(1:8)='ALL-DISP'
        elseif(textpart(1)(1:6).eq.'X-DISP') then
          settype='N'
          objectset(1,nobject)(1:6)='X-DISP'
        elseif(textpart(1)(1:6).eq.'Y-DISP') then
          settype='N'
          objectset(1,nobject)(1:6)='Y-DISP'
        elseif(textpart(1)(1:6).eq.'Z-DISP') then
          settype='N'
          objectset(1,nobject)(1:6)='Z-DISP'
        elseif(textpart(1)(1:14).eq.'EIGENFREQUENCY') then
          settype=' '
          nsets=0
          objectset(1,nobject)(1:14)='EIGENFREQUENCY'
        elseif(textpart(1)(1:5).eq.'GREEN') then
          settype=' '
          nsets=0
          objectset(1,nobject)(1:5)='GREEN'
        elseif(textpart(1)(1:4).eq.'MASS') then
          settype='E'
          objectset(1,nobject)(1:4)='MASS'
        elseif(textpart(1)(1:12).eq.'STRAINENERGY') then
          settype='E'
          objectset(1,nobject)(1:12)='STRAINENERGY'
        elseif(textpart(1)(1:11).eq.'MISESSTRESS') then
          settype='N'
          objectset(1,nobject)(1:11)='MISESSTRESS'
        elseif(textpart(1)(1:9).eq.'PS1STRESS') then
          settype='N'
          objectset(1,nobject)(1:9)='PS1STRESS'
        elseif(textpart(1)(1:9).eq.'PS3STRESS') then
          settype='N'
          objectset(1,nobject)(1:9)='PS3STRESS'
        elseif(textpart(1)(1:23).eq.'EQUIVALENTPLASTICSTRAIN') then
          settype='N'
          objectset(1,nobject)(1:15)='EQPLASTICSTRAIN'
        elseif(textpart(1)(1:11).eq.'MODALSTRESS') then
          settype='N'
          objectset(1,nobject)(1:11)='MODALSTRESS'
        else
          write(*,*) '*ERROR reading *DESIGN RESPONSE'
          write(*,*) '       unknown type', textpart(1)
          call inputerror(inpc,ipoinpc,iline,
     &         "*DESIGN RESPONSE%",ier)
          return
        endif
!     
        do m=1,nsets
          read(textpart(1+m)(1:80),'(a80)',iostat=istat) 
     &         objectset(2+m,nobject)(1:80) 
!     
!     appending the set character (N or E for nodal or element sets)
!     
          ipos=index(objectset(2+m,nobject),' ')
          if(ipos.eq.1) exit
          objectset(2+m,nobject)(ipos:ipos)=settype
!     
!     check if the given set exists. for this, we iterate over all existing
!     sets and compare the names. if the set does not exist, we throw an error.
!     
          i=0
          call cident81(set,objectset(2+m,nobject),nset,id)
          if(id.gt.0) then
            if(set(id).eq.objectset(2+m,nobject)) then
              i=id
            endif
          endif
!     
          if(i.gt.nset) then
            write(*,*) '*ERROR reading *DESIGN RESPONSE'
            write(*,*) '       unknown set name.'
            write(*,*) objectset(2+m,nobject)
            call inputerror(inpc,ipoinpc,iline,
     &           "*DESIGN RESPONSE%",ier)
            return
          endif
        enddo
!     
!     for stresses: parsing for the Kreisselmeier Steinhauser 
!     parameters rho and the target stress
!     
        if((objectset(1,nobject)(1:11).eq.'MISESSTRESS').or.
     &       (objectset(1,nobject)(1:9).eq.'PS1STRESS').or.
     &       (objectset(1,nobject)(1:9).eq.'PS3STRESS').or.
     &       (objectset(1,nobject)(1:15).eq.'EQPLASTICSTRAIN').or.
     &       (objectset(1,nobject)(1:11).eq.'MODALSTRESS'))then
          rho=0.d0
          stress=0.d0
!     
!     rho
!     
          if(n.ge.3) then
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) rho
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DESIGN RESPONSE%",ier)
              return
            endif
            objectset(2,nobject)(41:60)=textpart(3)(1:20)
          endif
!     
          if(icoordinate.eq.1) then
            if(rho.lt.1.d0) then
              write(*,*) '*ERROR reading *DESIGN RESPONSE'
              write(*,*) '       first Kreisselmeier-Steinhauser'
              write(*,*) '       parameter rho cannot be less'
              write(*,*) '       than 1.'
              call inputerror(inpc,ipoinpc,iline,
     &             "*DESIGN RESPONSE%",ier)
              return
            endif
          endif
!     
!     the target stress
!     
          if(n.ge.4) then
            read(textpart(4)(1:20),'(f20.0)',iostat=istat) stress
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DESIGN RESPONSE%",ier)
              return
            endif
            objectset(2,nobject)(61:80)=textpart(4)(1:20)
          endif
!     
          if((stress.le.0.d0).and.(objectset(1,nobject)(1:9).ne.
     &          'PS3STRESS')) then
            if(icoordinate.eq.1) then
              write(*,*) '*ERROR reading *DESIGN RESPONSE'
              write(*,*) '       the target stress in the'
              write(*,*) '       Kreisselmeier-Steinhauser function'
              write(*,*) '       must be strictly positive.'
              call inputerror(inpc,ipoinpc,iline,
     &             "*DESIGN RESPONSE%",ier)
              return
            endif
          elseif((stress.gt.0.d0).and.(objectset(1,nobject)(1:9).eq.
     &          'PS3STRESS')) then
            if(icoordinate.eq.1) then
              write(*,*) '*ERROR reading *DESIGN RESPONSE'
              write(*,*) '       the target stress for PS3STRESS in'
              write(*,*) '       the Kreisselmeier-Steinhauser'
              write(*,*) '       function must be strictly negative.'
              call inputerror(inpc,ipoinpc,iline,
     &             "*DESIGN RESPONSE%",ier)
              return
            endif
          endif
        endif
!     
        if(icoordinate.eq.1) exit
      enddo
!     
      return
      end
