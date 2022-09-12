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
      subroutine printout(set,nset,istartset,iendset,ialset,nprint,
     &     prlab,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener,
     &     mi,nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &     orab,ielorien,norien,nk,ne,inum,filab,vold,ikin,ielmat,
     &     thicke,eme,islavsurf,mortar,time,ielprop,prop,veold,orname,
     &     nelemload,nload,sideload,xload,rhcon,nrhcon,ntmat_,ipobody,
     &     ibody,xbody,nbody,nmethod)
!     
!     stores results in the .dat file
!     
      implicit none
!     
      logical force
!     
      character*1 cflag
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 noset,elset,orname(*)
      character*81 set(*),prset(*)
      character*87 filab(*)
!     
      integer nset,istartset(*),iendset(*),ialset(*),nprint,ipkon(*),
     &     mi(*),nstate_,ii,jj,iset,l,limit,node,ipos,ithermal(*),ielem,
     &     nelem,kon(*),inotr(2,*),ntrans,ielorien(mi(3),*),norien,nk,
     &     inum(*),nfield,ikin,nodes,ne0,nope,mt,ielmat(mi(3),*),iface,
     &     jfaces,mortar,islavsurf(2,*),ielprop(*),nload,i,ntmat_,id,
     &     nelemload(2,*),nrhcon(*),ipobody(2,*),ibody(3,*),nbody,
     &     nmethod,ne
!     
      real*8 v(0:mi(2),*),t1(*),fn(0:mi(2),*),stx(6,mi(1),*),bhetot,
     &     eei(6,mi(1),*),xstate(nstate_,mi(1),*),ener(2,mi(1),*),
     &     volumetot,co(3,*),qfx(3,mi(1),*),rftot(0:3),ttime,time,
     &     trab(7,*),orab(7,*),vold(0:mi(2),*),enerkintot,
     &     eme(6,mi(1),*),prop(*),veold(0:mi(2),*),xload(2,*),xmasstot,
     &     xinertot(6),cg(3),rhcon(0:1,ntmat_,*),xbody(7,*),energytot,
     &     thicke(mi(3),*)
!     
      mt=mi(2)+1
!     
!     interpolation in the original nodes of 1d and 2d elements
!     
      do ii=1,nprint
        if((prlab(ii)(1:4).eq.'U   ').or.
     &       ((prlab(ii)(1:4).eq.'NT  ').and.(ithermal(1).gt.1))) then
          if(filab(1)(5:5).ne.' ') then
            nfield=mt
            cflag=' '
            force=.false.
            call map3dto1d2d(v,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,force,mi,ielprop,prop)
          endif
          exit
        endif
      enddo
      do ii=1,nprint
        if((prlab(ii)(1:4).eq.'NT  ').and.(ithermal(1).le.1)) then
          if(filab(2)(5:5).ne.' ') then
            nfield=1
            cflag=' '
            force=.false.
            call map3dto1d2d(t1,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,force,mi,ielprop,prop)
          endif
          exit
        endif
      enddo
      do ii=1,nprint
        if(prlab(ii)(1:2).eq.'RF') then
          if(filab(1)(5:5).ne.' ') then
            nfield=mt
            cflag=' '
            force=.true.
            call map3dto1d2d(fn,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,force,mi,ielprop,prop)
          endif
          exit
        endif
      enddo
!     
      do ii=1,nprint
!     
!     nodal values
!     
        if((prlab(ii)(1:4).eq.'U   ').or.(prlab(ii)(1:4).eq.'NT  ').or.
     &      (prlab(ii)(1:4).eq.'RF  ').or.(prlab(ii)(1:4).eq.'RFL ').or. 
     &      (prlab(ii)(1:4).eq.'PS  ').or.(prlab(ii)(1:4).eq.'PN  ').or.
     &      (prlab(ii)(1:4).eq.'MF  ').or.(prlab(ii)(1:4).eq.'V   ').or.
     &       (prlab(ii)(1:4).eq.'TS  ')) 
     &       then
!     
          ipos=index(prset(ii),' ')
          noset='                    '
          noset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
!     printing the header
!     
          if(prlab(ii)(1:4).eq.'U   ') then
            write(5,*)
            if(mi(2).eq.3) then
              write(5,100) noset(1:ipos-2),ttime+time
            else
              write(5,133) noset(1:ipos-2),ttime+time
            endif
 100        format(' displacements (vx,vy,vz) for set ',A,
     &           ' and time ',e14.7)
 133        format(' displacements (v(i),i=1..ndof) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:4).eq.'NT  ').or.
     &           (prlab(ii)(1:4).eq.'TS  ')) then
            write(5,*)
            write(5,101) noset(1:ipos-2),ttime+time
 101        format(' temperatures for set ',A,' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'RF   ').or.
     &           (prlab(ii)(1:5).eq.'RF  T')) then
            write(5,*)
            write(5,102) noset(1:ipos-2),ttime+time
 102        format(' forces (fx,fy,fz) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'RFL ').or.
     &           (prlab(ii)(1:5).eq.'RFL T')) then
            write(5,*)
            write(5,103) noset(1:ipos-2),ttime+time
 103        format(' heat generation for set ',A,' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'PS  ') then
            write(5,*)
            write(5,115) noset(1:ipos-2),ttime+time
 115        format(' static pressures for set ',A,' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'PN  ') then
            write(5,*)
            write(5,117) noset(1:ipos-2),ttime+time
 117        format(   ' network pressures (total pressure for gases, sta
     &tic pressure for liquids and fluid depth for channels) for set ',
     &           A,' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'MF  ') then
            write(5,*)
            write(5,118) noset(1:ipos-2),ttime+time
 118        format(' mass flows for set ',A,' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'V   ') then
            write(5,*)
            write(5,119) noset(1:ipos-2),ttime+time
 119        format(' velocities (vx,vy,vz) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          endif
!     
!     printing the data
!     
          call cident81(set,prset(ii),nset,id)
          iset=nset+1
          if(id.gt.0) then
            if(prset(ii).eq.set(id)) then
              iset=id
            endif
          endif
          do jj=0,3
            rftot(jj)=0.d0
          enddo
          do jj=istartset(iset),iendset(iset)
            if(ialset(jj).lt.0) cycle
            if(jj.eq.iendset(iset)) then
              node=ialset(jj)
              call printoutnode(prlab,v,t1,fn,ithermal,ii,node,
     &             rftot,trab,inotr,ntrans,co,mi,veold)
            elseif(ialset(jj+1).gt.0) then
              node=ialset(jj)
              call printoutnode(prlab,v,t1,fn,ithermal,ii,node,
     &             rftot,trab,inotr,ntrans,co,mi,veold)
            else
              do node=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &             -ialset(jj+1)
                call printoutnode(prlab,v,t1,fn,ithermal,ii,node,
     &               rftot,trab,inotr,ntrans,co,mi,veold)
              enddo
            endif
          enddo
!     
!     writing total values to file
!     
          if((prlab(ii)(1:5).eq.'RF  O').or.
     &         (prlab(ii)(1:5).eq.'RF  T')) then
            write(5,*)
            write(5,104) noset(1:ipos-2),ttime+time
 104        format(' total force (fx,fy,fz) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,3(1x,e13.6))') rftot(1),rftot(2),rftot(3)
          elseif((prlab(ii)(1:5).eq.'RFL O').or.
     &           (prlab(ii)(1:5).eq.'RFL T')) then
            write(5,*)
            write(5,105)noset(1:ipos-2),ttime+time
 105        format(' total heat generation for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') rftot(0)
          endif
!     
!     integration point values
!     
        elseif((prlab(ii)(1:4).eq.'S   ').or.
     &         (prlab(ii)(1:4).eq.'E   ').or.
     &         (prlab(ii)(1:4).eq.'ME  ').or.
     &         (prlab(ii)(1:4).eq.'PEEQ').or.
     &         (prlab(ii)(1:4).eq.'ENER').or.
     &         (prlab(ii)(1:4).eq.'SDV ').or.
     &         (prlab(ii)(1:4).eq.'COOR').or.
     &         (prlab(ii)(1:4).eq.'HFL ')) then
!     
          ipos=index(prset(ii),' ')
          elset='                    '
          elset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
          limit=1
!     
          do l=1,limit
!     
!     printing the header
!     
            if(prlab(ii)(1:4).eq.'S   ') then
              write(5,*)
              write(5,106) elset(1:ipos-2),ttime+time
 106          format(    ' stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sx
     &z,syz) for set ',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'E   ') then
              write(5,*)
              write(5,107) elset(1:ipos-2),ttime+time
 107          format(    ' strains (elem, integ.pnt.,exx,eyy,ezz,exy,exz
     &,eyz) for set ',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'PEEQ') then
              write(5,*)
              write(5,108) elset(1:ipos-2),ttime+time
 108          format(    ' equivalent plastic strain (elem, integ.pnt.,p 
     &e)for set ',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'ENER') then
              write(5,*)
              write(5,109) elset(1:ipos-2),ttime+time
 109          format(    ' internal energy density (elem, integ.pnt.,ene
     &rgy) for set ',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'SDV ') then
              write(5,*)
              write(5,111) elset(1:ipos-2),ttime+time
 111          format
     &           (' internal state variables (elem, integ.pnt.,values) f
     &or set ',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'HFL ') then
              write(5,*)
              write(5,112) elset(1:ipos-2),ttime+time
 112          format(    ' heat flux (elem, integ.pnt.,qx,qy,qz) for set 
     &',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'ME  ') then
              write(5,*)
              write(5,130) elset(1:ipos-2),ttime+time
 130          format(    ' mechanical strains (elem, integ.pnt.,exx,eyy,
     &ezz,exy,exz,eyz) for set ',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'COOR') then
              write(5,*)
              write(5,139) elset(1:ipos-2),ttime+time
 139          format(    ' global coordinates (elem, integ.pnt.,x,y,z) f
     &or set ',A,' and time ',e14.7)
              write(5,*)
            endif
!     
!     printing the data
!     
            call cident81(set,prset(ii),nset,id)
            iset=nset+1
            if(id.gt.0) then
              if(prset(ii).eq.set(id)) then
                iset=id
              endif
            endif
            do jj=istartset(iset),iendset(iset)
              if(ialset(jj).lt.0) cycle
              if(jj.eq.iendset(iset)) then
                nelem=ialset(jj)
                call printoutint(prlab,ipkon,lakon,stx,eei,xstate,
     &               ener,mi(1),nstate_,ii,nelem,qfx,
     &               orab,ielorien,norien,co,kon,ielmat,thicke,eme,
     &               ielprop,prop,nelem,ithermal,orname)
              elseif(ialset(jj+1).gt.0) then
                nelem=ialset(jj)
                call printoutint(prlab,ipkon,lakon,stx,eei,xstate,
     &               ener,mi(1),nstate_,ii,nelem,qfx,orab,
     &               ielorien,norien,co,kon,ielmat,thicke,eme,
     &               ielprop,prop,nelem,ithermal,orname)
              else
                do nelem=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &               -ialset(jj+1)
                  call printoutint(prlab,ipkon,lakon,stx,eei,
     &                 xstate,ener,mi(1),nstate_,ii,nelem,
     &                 qfx,orab,ielorien,norien,co,kon,ielmat,
     &                 thicke,eme,ielprop,prop,nelem,ithermal,
     &                 orname)
                enddo
              endif
            enddo
!     
          enddo
!     
!     whole element values
!     
        elseif((prlab(ii)(1:4).eq.'ELSE').or.
     &         (prlab(ii)(1:4).eq.'ELKE').or.
     &         (prlab(ii)(1:4).eq.'EVOL').or.
     &         (prlab(ii)(1:4).eq.'EMAS').or.
     &         (prlab(ii)(1:4).eq.'EBHE').or.
     &         (prlab(ii)(1:4).eq.'CENT').or.
     &         (prlab(ii)(1:4).eq.'CSTR').or.
     &         (prlab(ii)(1:4).eq.'CDIS').or.
     &         (prlab(ii)(1:4).eq.'CNUM').or.
     &         (prlab(ii)(1:4).eq.'CELS')) then
!     
          ipos=index(prset(ii),' ')
          elset='                    '
          elset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
!     printing the header
!     
          if((prlab(ii)(1:5).eq.'ELSE ').or.
     &         (prlab(ii)(1:5).eq.'ELSET')) then
            write(5,*)
            write(5,113) elset(1:ipos-2),ttime+time
 113        format(' internal energy (element, energy) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'ELKE ').or.
     &           (prlab(ii)(1:5).eq.'ELKET')) then
            write(5,*)
            write(5,110) elset(1:ipos-2),ttime+time
 110        format(' kinetic energy (elem, energy) for set '
     &           ,A,' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'EVOL ').or.
     &           (prlab(ii)(1:5).eq.'EVOLT')) then
            write(5,*)
            write(5,114) elset(1:ipos-2),ttime+time
 114        format(' volume (element, volume) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'EMAS ').or.
     &           (prlab(ii)(1:5).eq.'EMAST')) then
            write(5,*)
            write(5,136) elset(1:ipos-2),ttime+time
 136        format(    ' mass (element, mass) and mass moment of itertia 
     &(xx,yy,zz,xy,xz,yz) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'EBHE ').or.
     &           (prlab(ii)(1:5).eq.'EBHET')) then
            if(nmethod.eq.2) cycle
            write(5,*)
            write(5,131) elset(1:ipos-2),ttime+time
 131        format(' body heating (element, volume) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'CSTR ').or.
     &           (prlab(ii)(1:5).eq.'CSTRT')) then
            write(5,*)
            if(mortar.eq.0) then
              write(5,122) ttime+time
 122          format(' contact stress (slave node,press,'
     &             'tang1,tang2) for all contact elements and time',
     &             e14.7)
            elseif(mortar.eq.1) then
              write(5,126) ttime+time
 126          format(' contact stress (slave element+face,press,'
     &             'tang1,tang2) for all contact elements and time',
     &             e14.7)
            endif
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'CDIS ').or.
     &           (prlab(ii)(1:5).eq.'CDIST')) then
            write(5,*)
            if(mortar.eq.0) then
              write(5,123) ttime+time
 123          format(' relative contact displacement (slave node,'
     &             'normal,tang1,tang2) for all contact elements and '
     &             'time',e14.7)
            elseif(mortar.eq.1) then
              write(5,127) ttime+time
 127          format(
     &             ' relative contact displacement (slave element+face,'
     &             'normal,tang1,tang2) for all contact elements and '
     &             'time',e14.7)
            endif
            write(5,*)
          elseif((prlab(ii)(1:5).eq.'CELS ').or.
     &           (prlab(ii)(1:5).eq.'CELST')) then
            write(5,*)
            if(mortar.eq.0) then
              write(5,124) ttime+time
 124          format(' contact spring energy (slave node,energy) for '  
     &             'all contact elements and time',e14.7)
            elseif(mortar.eq.1) then
              write(5,128) ttime+time
 128          format(
     &             ' contact spring energy (slave element+face,energy) '  
     &             'for all contact elements and time',e14.7)
            endif
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'CENT') then
            write(5,*)
            write(5,140)  elset(1:ipos-2),ttime+time
 140        format(' centrifugal force(element,omega square) for '  
     &           'set ',A,' and time ',e14.7)
            write(5,*)
          endif
!     
!     printing the data
!     
          volumetot=0.d0
          bhetot=0.d0
          energytot=0.d0
          enerkintot=0.d0
          xmasstot=0.d0
          do jj=1,6
            xinertot(jj)=0.d0
          enddo
          do jj=1,3
            cg(jj)=0.d0
          enddo
!     
          if ((prlab(ii)(1:4).eq.'CSTR').or.
     &         (prlab(ii)(1:4).eq.'CDIS').or.
     &         (prlab(ii)(1:4).eq.'CNUM').or.
     &         (prlab(ii)(1:4).eq.'CELS')) then
!     
!     ne0 is the number of the first contact element
!     
            do jj=ne,1,-1
              if((lakon(jj)(2:2).ne.'S').or.
     &             (lakon(jj)(7:7).ne.'C')) then
                ne0=jj+1
                exit
              endif
            enddo
!     
            if(prlab(ii)(1:4).ne.'CNUM') then
              if(mortar.eq.0) then
                do nelem=ne0,ne
                  read(lakon(nelem)(8:8),'(i1)') nope
                  nope=nope+1
                  nodes=kon(ipkon(nelem)+nope)
                  call printoutelem(prlab,ipkon,lakon,kon,co,
     &                 ener,mi(1),ii,nelem,energytot,volumetot,
     &                 enerkintot,ne,stx,nodes,thicke,ielmat,
     &                 ielem,iface,mortar,ielprop,prop,
     &                 sideload,nload,nelemload,xload,bhetot,
     &                 xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &                 ntmat_,t1,vold,ipobody,ibody,xbody,nbody)
                enddo
              elseif(mortar.eq.1) then
                do nelem=ne0,ne
                  jfaces=
     &              islavsurf(1,kon(ipkon(nelem)+kon(ipkon(nelem))+2))
                  ielem=int(jfaces/10.d0)
                  iface=jfaces-10*ielem
                  call printoutelem(prlab,ipkon,lakon,kon,co,
     &                 ener,mi(1),ii,nelem,energytot,volumetot,
     &                 enerkintot,ne,stx,nodes,thicke,ielmat,
     &                 ielem,iface,mortar,ielprop,prop,
     &                 sideload,nload,nelemload,xload,bhetot,
     &                 xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &                 ntmat_,t1,vold,ipobody,ibody,xbody,nbody)
                enddo
              endif
            endif
          else
            call cident81(set,prset(ii),nset,id)
            iset=nset+1
            if(id.gt.0) then
              if(prset(ii).eq.set(id)) then
                iset=id
              endif
            endif
            do jj=istartset(iset),iendset(iset)
              if(ialset(jj).lt.0) cycle
              if(jj.eq.iendset(iset)) then
                nelem=ialset(jj)
                call printoutelem(prlab,ipkon,lakon,kon,co,
     &               ener,mi(1),ii,nelem,energytot,volumetot,
     &               enerkintot,ne,stx,nodes,thicke,ielmat,
     &               ielem,iface,mortar,ielprop,prop,
     &               sideload,nload,nelemload,xload,bhetot,
     &               xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &               ntmat_,t1,vold,ipobody,ibody,xbody,nbody)
              elseif(ialset(jj+1).gt.0) then
                nelem=ialset(jj)
                call printoutelem(prlab,ipkon,lakon,kon,co,
     &               ener,mi(1),ii,nelem,energytot,volumetot,
     &               enerkintot,ne,stx,nodes,thicke,ielmat,
     &               ielem,iface,mortar,ielprop,prop,
     &               sideload,nload,nelemload,xload,bhetot,
     &               xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &               ntmat_,t1,vold,ipobody,ibody,xbody,nbody)
              else
                do nelem=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &               -ialset(jj+1)
                  call printoutelem(prlab,ipkon,lakon,kon,co,
     &                 ener,mi(1),ii,nelem,energytot,volumetot,
     &                 enerkintot,ne,stx,nodes,thicke,ielmat,
     &                 ielem,iface,mortar,ielprop,prop,
     &                 sideload,nload,nelemload,xload,bhetot,
     &                 xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &                 ntmat_,t1,vold,ipobody,ibody,xbody,nbody)
                enddo
              endif
            enddo
          endif
!     
!     writing total values to file
!     
          if((prlab(ii)(1:5).eq.'ELSEO').or.
     &         (prlab(ii)(1:5).eq.'ELSET')) then
            write(5,*)
            write(5,116) elset(1:ipos-2),ttime+time
 116        format(' total internal energy for set ',A,' and time ',
     &           e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') energytot
          elseif((prlab(ii)(1:5).eq.'ELKEO').or.
     &           (prlab(ii)(1:5).eq.'ELKET')) then
            write(5,*)
            write(5,120) elset(1:ipos-2),ttime+time
 120        format(' total kinetic energy for set ',A,' and time ',
     &           e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') enerkintot
          elseif((prlab(ii)(1:5).eq.'EVOLO').or.
     &           (prlab(ii)(1:5).eq.'EVOLT')) then
            write(5,*)
            write(5,121) elset(1:ipos-2),ttime+time
 121        format(' total volume for set ',A,' and time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') volumetot
          elseif((prlab(ii)(1:5).eq.'EMASO').or.
     &           (prlab(ii)(1:5).eq.'EMAST')) then
            write(5,*)
            write(5,137) elset(1:ipos-2),ttime+time
 137        format(' total mass for set ',A,' and time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') xmasstot
            write(5,*)
            write(5,134) elset(1:ipos-2),ttime+time
 134        format(
     &     ' total mass moment of inertia (xx,yy,zz,xy,xz,yz) for set ',
     &A,' and time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,6(1x,e13.6))') (xinertot(i),i=1,6)
            write(5,*)
            write(5,135) elset(1:ipos-2),ttime+time
 135        format(' center of gravity for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,3(1x,e13.6))') (cg(i)/xmasstot,i=1,3)
            write(5,*)
            write(5,138) elset(1:ipos-2),ttime+time
 138        format(   ' total mass moment of inertia about the center of
     &gravity (xx,yy,zz,xy,xz,yz) for set ',
     &           A,' and time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,6(1x,e13.6))')
     &           xinertot(1)-cg(1)*cg(1)/xmasstot,
     &           xinertot(2)-cg(2)*cg(2)/xmasstot,
     &           xinertot(3)-cg(3)*cg(3)/xmasstot,
     &           xinertot(4)-cg(1)*cg(2)/xmasstot,
     &           xinertot(5)-cg(1)*cg(3)/xmasstot,
     &           xinertot(6)-cg(2)*cg(3)/xmasstot
          elseif((prlab(ii)(1:5).eq.'EBHEO').or.
     &           (prlab(ii)(1:5).eq.'EBHET')) then
            if(nmethod.eq.2) cycle
            write(5,*)
            write(5,132) elset(1:ipos-2),ttime+time
 132        format(' total body heating for set ',A,' and time ',
     &           e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') bhetot
          elseif((prlab(ii)(1:5).eq.'CELSO').or.
     &           (prlab(ii)(1:5).eq.'CELST')) then
            write(5,*)
            write(5,125) ttime+time
 125        format(' total contact spring energy for time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') energytot
          elseif(prlab(ii)(1:4).eq.'CNUM') then
            write(5,*)
            write(5,129) ttime+time
 129        format
     &           (' total number of contact elements for time ',e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,i10)') ne-ne0+1
!     
          endif
        endif
      enddo
!     
      return
      end
