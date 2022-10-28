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
      subroutine contactpairs(inpc,textpart,tieset,istep,
     &                istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,
     &                iperturb,matname,nmat,ipoinpc,tietol,set,nset,
     &                mortar,ncmat_,ntmat_,elcon,ier)
!
!     reading the input deck: *CONTACT PAIR
!
      implicit none
!
      logical linear
!
      character*1 inpc(*)
      character*80 matname(*),material
      character*81 tieset(3,*),noset,set(*)
      character*132 textpart(16)
!
      integer istep,istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ntie,ntie_,iperturb(*),nmat,ipoinpc(0:*),nset,j,
     &  mortar,ncmat_,ntmat_,ier,id
!
      real*8 tietol(4,*),adjust,elcon(0:ncmat_,ntmat_,*)
!
!     tietol contains information on:
!            - small (tietol<0) or large (tietol>0) sliding
!            - the adjust value (only if dabs(tietol)>=2,
!                 adjust=dabs(tietol)-2
!
      if(istep.gt.0) then
         write(*,*) '*ERROR reading *CONTACT PAIR: *CONTACT PAIR should'
         write(*,*) '  be placed before all step definitions'
         ier=1
         return
      endif
!
      mortar=-2
      linear=.false.
!
      ntie=ntie+1
      if(ntie.gt.ntie_) then
         write(*,*) '*ERROR reading *CONTACT PAIR: increase ntie_'
         ier=1
         return
      endif
      tietol(1,ntie)=1.d0
!
!     default for "no clearance"
!
      tietol(3,ntie)=1.2357111317d0
      do j=1,80
         tieset(1,ntie)(j:j)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:12).eq.'INTERACTION=') then
            material=textpart(i)(13:92)
         elseif(textpart(i)(1:12).eq.'SMALLSLIDING') then
            tietol(1,ntie)=-tietol(1,ntie)
         elseif(textpart(i)(1:6).eq.'LINEAR') then
            linear=.true.
         elseif(textpart(i)(1:7).eq.'ADJUST=') then
            read(textpart(i)(8:25),'(f20.0)',iostat=istat) adjust
            if(istat.gt.0) then
               noset(1:80)=textpart(i)(8:87)
               noset(81:81)=' '
               ipos=index(noset,' ')
               noset(ipos:ipos)='N'
               call cident81(set,noset,nset,id)
               j=nset+1
               if(id.gt.0) then
                 if(noset.eq.set(id)) then
                   j=id
                 endif
               endif
               if(j.gt.nset) then
                  noset(ipos:ipos)=' '
                  write(*,*) 
     &               '*ERROR reading *CONTACT PAIR: adjust node set',
     &                   noset
                  write(*,*) '       has not been defined'
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*CONTACT PAIR%",ier)
                  return
               endif
               do j=1,ipos-1
                  tieset(1,ntie)(j:j)=noset(j:j)
               enddo
               do j=ipos,80
                  tieset(1,ntie)(j:j)=' '
               enddo
            else
               tietol(1,ntie)=dsign(1.d0,tietol(1,ntie))*(2.d0+adjust)
            endif
         elseif(textpart(i)(1:13).eq.'TYPE=MASSLESS') then
            mortar=-1
         elseif(textpart(i)(1:18).eq.'TYPE=NODETOSURFACE') then
            mortar=0
         elseif(textpart(i)(1:21).eq.'TYPE=SURFACETOSURFACE') then
            mortar=1
         elseif(textpart(i)(1:11).eq.'TYPE=MORTAR') then
            mortar=2
         elseif(textpart(i)(1:13).eq.'TYPE=PGMORTAR') then
            mortar=5
         elseif(textpart(i)(1:14).eq.'TYPE=LINMORTAR') then
            mortar=3
         elseif(textpart(i)(1:16).eq.'TYPE=PGLINMORTAR') then
            mortar=4
         else
            write(*,*) 
     &       '*WARNING reading *CONTACT PAIR: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CONTACT PAIR%")
         endif
      enddo
!
      if(mortar.eq.-2) then
         write(*,*) '*ERROR reading *CONTACT PAIR'
         write(*,*) '       no TYPE specified'
         call inputerror(inpc,ipoinpc,iline,
     &        "*CONTACT PAIR%",ier)
         return
      endif
!
!     SMALL SLIDING significates that the number of contact elements
!     within one increment is frozen. This is not allowed for
!     SURFACE TO SURFACE contact.
!
      if((tietol(1,ntie).lt.0.d0).and.(mortar.eq.1)) then
         write(*,*) '*WARNING reading *CONTACT PAIR'
         write(*,*) '         The option SMALL SLIDING cannot be'
         write(*,*) '         used with SURFACE TO SURFACE contact'
         write(*,*)
         tietol(1,ntie)=-tietol(1,ntie)
      endif
!
!     check for the existence of the surface interaction
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) '*ERROR reading *CONTACT PAIR: nonexistent surface'
         write(*,*) '       interaction; '
         call inputerror(inpc,ipoinpc,iline,
     &        "*CONTACT PAIR%",ier)
         return
      endif
      tietol(2,ntie)=i+0.5d0
!
!     check whether sigma_at_infinity is given for node-to-face penalty 
!     contact with a linear pressure-overclosure relationship
!
      if(ncmat_<2) then
         write(*,*) 
     &     '*ERROR reading *CONTACT PAIR: no PRESSURE-OVERCLOSURE'
         write(*,*) 
     &   '       has been defined for at least one *SURFACE INTERACTION'
         ier=1
         return
      elseif((mortar.eq.0).or.(mortar.eq.1)) then
         if(ncmat_<3) then
            write(*,*) 
     &           '*ERROR reading *CONTACT PAIR: no PRESSURE-OVERCLOSURE'
            write(*,*) 
     &   '       has been defined for at least one *SURFACE INTERACTION'
            ier=1
            return
         elseif(int(elcon(3,1,i)).le.0) then
            write(*,*) 
     &           '*ERROR reading *CONTACT PAIR: no PRESSURE-OVERCLOSURE'
            write(*,*) 
     &   '       has been defined for at least one *SURFACE INTERACTION'
            ier=1
            return
         endif
      endif
!
      if(ncmat_.ge.3) then
         if(int(elcon(3,1,i)).eq.2) then
            if(mortar.eq.0) then
               if(elcon(1,1,i).lt.1.d-30) then
                  write(*,*) '*ERROR reading *CONTACT PAIR:'
                  write(*,*) '       for node-to-face penalty contact'
                  write(*,*) '       with linear pressure-overclosure'
                  write(*,*) '       relationship, the'
                  write(*,*) '       tension at large clearances'
                  write(*,*) '       must exceed 1.e-30'
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*CONTACT PAIR%",ier)
                  return
               endif
            endif
         endif
      endif
!
      tieset(1,ntie)(81:81)='C'
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)'*ERROR reading *CONTACT PAIR: definition of the '
         write(*,*) '      contact pair is not complete.'
         ier=1
         return
      endif
!
!     storing the slave surface
!
      if(mortar.eq.1) then
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='T'
      elseif(mortar.eq.2) then
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='M'
      elseif(mortar.eq.3) then
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='O'
      elseif(mortar.eq.4) then
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='P'      
      elseif(mortar.eq.5) then
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='G'
      else
         tieset(2,ntie)(1:80)=textpart(1)(1:80)
         tieset(2,ntie)(81:81)=' '
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='S'
      endif
!
      tieset(3,ntie)(1:80)=textpart(2)(1:80)
      tieset(3,ntie)(81:81)=' '
      ipos=index(tieset(3,ntie),' ')
      tieset(3,ntie)(ipos:ipos)='T'
!
!     the definition of a contact pair triggers a call to
!     nonlingeo (for static calculations) but not automatically
!     to the nonlinear calculation of strains (i.e.
!     iperturb(2) should be zero unless NLGEOM is activated)
!
      if((iperturb(1).eq.0).and.(.not.linear)) then
         iperturb(1)=2
      endif
!
!     check for further contact pairs
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
!
         ntie=ntie+1
!
!        copying the information from the previous tie
!
         do i=1,3
            tietol(i,ntie)=tietol(i,ntie-1)
         enddo
         tieset(1,ntie)=tieset(1,ntie-1)
!
!        storing slave and master surface
!
         if(mortar.eq.1) then
            tieset(2,ntie)(1:80)=textpart(1)(1:80)
            tieset(2,ntie)(81:81)=' '
            ipos=index(tieset(2,ntie),' ')
            tieset(2,ntie)(ipos:ipos)='T'
         elseif(mortar.eq.2) then
            tieset(2,ntie)(1:80)=textpart(1)(1:80)
            tieset(2,ntie)(81:81)=' '
            ipos=index(tieset(2,ntie),' ')
            tieset(2,ntie)(ipos:ipos)='M'
         elseif(mortar.eq.3) then
            tieset(2,ntie)(1:80)=textpart(1)(1:80)
            tieset(2,ntie)(81:81)=' '
            ipos=index(tieset(2,ntie),' ')
            tieset(2,ntie)(ipos:ipos)='O'
         elseif(mortar.eq.4) then
            tieset(2,ntie)(1:80)=textpart(1)(1:80)
            tieset(2,ntie)(81:81)=' '
            ipos=index(tieset(2,ntie),' ')
            tieset(2,ntie)(ipos:ipos)='P'      
         elseif(mortar.eq.5) then
            tieset(2,ntie)(1:80)=textpart(1)(1:80)
            tieset(2,ntie)(81:81)=' '
            ipos=index(tieset(2,ntie),' ')
            tieset(2,ntie)(ipos:ipos)='G'
         else
            tieset(2,ntie)(1:80)=textpart(1)(1:80)
            tieset(2,ntie)(81:81)=' '
            ipos=index(tieset(2,ntie),' ')
            tieset(2,ntie)(ipos:ipos)='S'
         endif
!
         tieset(3,ntie)(1:80)=textpart(2)(1:80)
         tieset(3,ntie)(81:81)=' '
         ipos=index(tieset(3,ntie),' ')
         tieset(3,ntie)(ipos:ipos)='T'
      enddo
!
      return
      end



