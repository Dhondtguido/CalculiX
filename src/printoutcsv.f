!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine printoutcsv(set,nset,istartset,iendset,ialset,nprint,
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
      character*1 cflag
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 noset,elset,orname(*)
      character*81 set(*),prset(*)
      character*87 filab(*),filename
!     
      integer nset,istartset(*),iendset(*),ialset(*),nprint,ipkon(*),
     &     mi(*),nstate_,ii,jj,iset,l,limit,node,ipos,ithermal(*),ielem,
     &     nelem,kon(*),inotr(2,*),ntrans,ielorien(mi(3),*),norien,nk,
     &     inum(*),nfield,ikin,nodes,ne0,nope,mt,ielmat(mi(3),*),iface,
     &     jfaces,mortar,islavsurf(2,*),ielprop(*),nload,i,ntmat_,id,
     &     nelemload(2,*),nrhcon(*),ipobody(2,*),ibody(3,*),nbody,
     &     nmethod,ne,iforce,u_unit,f_unit,s_unit,e_unit,pe_unit,
     &     ener_unit,svd_unit,me_unit,coor_unit,else_unit,elke_unit,
     &     evol_unit,emas_unit,cent_unit
!     
      real*8 v(0:mi(2),*),t1(*),fn(0:mi(2),*),stx(6,mi(1),*),bhetot,
     &     eei(6,mi(1),*),xstate(nstate_,mi(1),*),ener(2,mi(1),*),
     &     volumetot,co(3,*),qfx(3,mi(1),*),rftot(0:3),ttime,time,
     &     trab(7,*),orab(7,*),vold(0:mi(2),*),enerkintot,
     &     eme(6,mi(1),*),prop(*),veold(0:mi(2),*),xload(2,*),xmasstot,
     &     xinertot(6),cg(3),rhcon(0:1,ntmat_,*),xbody(7,*),energytot,
     &     thicke(mi(3),*), random_unit
!     
      mt=mi(2)+1
      call random_number(random_unit)
      u_unit=nint(random_unit*1000)
      f_unit=u_unit+1
      s_unit=f_unit+1
      e_unit=s_unit+1
      pe_unit=e_unit+1
      ener_unit=pe_unit+1
      svd_unit=ener_unit+1
      me_unit=svd_unit+1
      coor_unit=me_unit+1
      else_unit=coor_unit+1
      elke_unit=else_unit+1
      evol_unit=elke_unit+1
      emas_unit=evol_unit+1
      cent_unit=emas_unit+1

!     
!     interpolation in the original nodes of 1d and 2d elements
!     
      do ii=1,nprint
        if((prlab(ii)(1:4).eq.'U   ').or.
     &       ((prlab(ii)(1:4).eq.'NT  ').and.(ithermal(1).gt.1))) then
          if(filab(1)(5:5).ne.' ') then
            nfield=mt
            cflag=' '
            iforce=0
            call map3dto1d2d(v,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,iforce,mi,ielprop,prop)
          endif
          exit
        endif
      enddo
      do ii=1,nprint
        if((prlab(ii)(1:4).eq.'NT  ').and.(ithermal(1).le.1)) then
          if(filab(2)(5:5).ne.' ') then
            nfield=1
            cflag=' '
            iforce=0
            call map3dto1d2d(t1,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,iforce,mi,ielprop,prop)
          endif
          exit
        endif
      enddo
      do ii=1,nprint
        if(prlab(ii)(1:2).eq.'RF') then
          if(filab(1)(5:5).ne.' ') then
            nfield=mt
            cflag=' '
            iforce=1
            call map3dto1d2d(fn,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,iforce,mi,ielprop,prop)
          endif
          exit
        endif
      enddo
!     
      do ii=1,nprint
!     
!     nodal values
!     
        if((prlab(ii)(1:4).eq.'U   ').or.(prlab(ii)(1:4).eq.'RF  '))
     &       then
!     
          ipos=index(prset(ii),' ')
          noset='                    '
          noset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
!     printing the header
!     
1000      format(A,'_'A,'_',e14.7,'.csv')
100       format('# displacements (vx,vy,vz) for set ',A,
     &           ' and time ',e14.7)
101       format('node,vx,vy,vz')
102       format('# forces (fx,fy,fz) for set ',A,
     &           ' and time ',e14.7)
103       format('node,fx,fy,fz')

          if(prlab(ii)(1:4).eq.'U   ') then
            write(filename, 1000) 'U',noset(1:ipos-2),ttime+time
            open(unit=u_unit,file=filename, form='formatted')
            write(u_unit,100) noset(1:ipos-2),ttime+time
            write(u_unit,101)
          elseif((prlab(ii)(1:5).eq.'RF   ').or.
     &           (prlab(ii)(1:5).eq.'RF  T')) then
            write(filename, 1000) 'RF', noset(1:ipos-2),ttime+time
            open(unit=f_unit,file=filename, form='formatted')
            write(f_unit,102) noset(1:ipos-2),ttime+time
            write(f_unit,103)
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
              call printoutnodecsv(prlab,v,t1,fn,ithermal,ii,node,
     &             rftot,trab,inotr,ntrans,co,mi,veold,u_unit,f_unit)
            elseif(ialset(jj+1).gt.0) then
              node=ialset(jj)
              call printoutnodecsv(prlab,v,t1,fn,ithermal,ii,node,
     &             rftot,trab,inotr,ntrans,co,mi,veold,u_unit,f_unit)
            else
              do node=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &             -ialset(jj+1)
                call printoutnodecsv(prlab,v,t1,fn,ithermal,ii,node,
     &               rftot,trab,inotr,ntrans,co,mi,veold,u_unit,f_unit)
              enddo
            endif
          enddo

!     
!     integration point values
!     
        elseif((prlab(ii)(1:4).eq.'S   ').or.
     &         (prlab(ii)(1:4).eq.'E   ').or.
     &         (prlab(ii)(1:4).eq.'ME  ').or.
     &         (prlab(ii)(1:4).eq.'PEEQ').or.
     &         (prlab(ii)(1:4).eq.'ENER').or.
     &         (prlab(ii)(1:4).eq.'SDV ').or.
     &         (prlab(ii)(1:4).eq.'COOR')) then
!     
          ipos=index(prset(ii),' ')
          elset='                    '
          elset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
          limit=1
          write(filename,1000) trim(prlab(ii)(1:4)),
     &                         elset(1:ipos-2),ttime+time
          
!     
          do l=1,limit

!     
!     printing the header
!     
104       format('# stresses (element, integ.pnt.,sxx,syy,szz,sxy,sxz'
     &',syz) for set ',A,' and time ',e14.7)
1041       format('element,integ.pnt,sxx,syy,szz,sxy,sxz,syz')

105       format('# strains (element, integ.pnt.,exx,eyy,ezz,exy,exz'
     &',eyz) for set ',A,' and time ',e14.7)
1051       format('element,integ.pnt,exx,eyy,ezz,exy,exz,eyz')

106       format('# equivalent plastic strain (element, integ.pnt.,pe)' 
     &' for set ',A,' and time ',e14.7)
1061       format('element, integ.pnt,pe')

107       format('# internal energy density (element, integ.pnt.'
     &',energy), for set ',A,' and time ',e14.7)
1071       format('element,integ.pnt,energy')

108       format('# internal state variables (element, integ.pnt.'
     &',values), for set ',A,' and time ',e14.7)
1081       format('element,integ.pnt,values')

109       format('# mechanical strains (element, integ.pnt.,exx,eyy,'
     &' ezz,exy,exz,eyz) for set ',A,' and time ',e14.7)
1091       format('element,integ.pnt,exx,eyy,ezz,exy,exz,eyz')

110       format('# global coordinates (element, integ.pnt.,x,y,z) '
     &' for set ',A,' and time ',e14.7)
1101       format('element,integ.pnt,x,y,z')

            if(prlab(ii)(1:4).eq.'S   ') then
              open(unit=s_unit,file=filename, form='formatted')
              write(s_unit,104) elset(1:ipos-2),ttime+time
              write(s_unit,1041)
            elseif(prlab(ii)(1:4).eq.'E   ') then
              open(unit=e_unit,file=filename, form='formatted')
              write(e_unit,105) elset(1:ipos-2),ttime+time
              write(e_unit,1051)
            elseif(prlab(ii)(1:4).eq.'PEEQ') then
              open(unit=pe_unit,file=filename, form='formatted')
              write(pe_unit,106) elset(1:ipos-2),ttime+time
              write(pe_unit,1061)
            elseif(prlab(ii)(1:4).eq.'ENER') then
              open(unit=ener_unit,file=filename, form='formatted')
              write(ener_unit,107) elset(1:ipos-2),ttime+time
              write(ener_unit,1071)
            elseif(prlab(ii)(1:4).eq.'SDV ') then
              open(unit=svd_unit,file=filename, form='formatted')
              write(svd_unit,108) elset(1:ipos-2),ttime+time
              write(svd_unit,1081)
            elseif(prlab(ii)(1:4).eq.'ME  ') then
              open(unit=me_unit,file=filename, form='formatted')
              write(me_unit,109) elset(1:ipos-2),ttime+time
              write(me_unit,1091)
            elseif(prlab(ii)(1:4).eq.'COOR') then
              open(unit=coor_unit,file=filename, form='formatted')
              write(coor_unit,110) elset(1:ipos-2),ttime+time
              write(coor_unit,1101)
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
                call printoutintcsv(prlab,ipkon,lakon,stx,eei,xstate,
     &               ener,mi(1),nstate_,ii,nelem,qfx,
     &               orab,ielorien,norien,co,kon,ielmat,thicke,eme,
     &               ielprop,prop,nelem,ithermal,orname,s_unit,e_unit,
     &               pe_unit,ener_unit,svd_unit,me_unit,coor_unit)
              elseif(ialset(jj+1).gt.0) then
                nelem=ialset(jj)
                call printoutintcsv(prlab,ipkon,lakon,stx,eei,xstate,
     &               ener,mi(1),nstate_,ii,nelem,qfx,orab,
     &               ielorien,norien,co,kon,ielmat,thicke,eme,
     &               ielprop,prop,nelem,ithermal,orname,s_unit,e_unit,
     &               pe_unit,ener_unit,svd_unit,me_unit,coor_unit)
              else
                do nelem=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &               -ialset(jj+1)
                  call printoutintcsv(prlab,ipkon,lakon,stx,eei,
     &                 xstate,ener,mi(1),nstate_,ii,nelem,
     &                 qfx,orab,ielorien,norien,co,kon,ielmat,
     &                 thicke,eme,ielprop,prop,nelem,ithermal,
     &                 orname,s_unit,e_unit,pe_unit,ener_unit,svd_unit,
     &                 me_unit,coor_unit)
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
     &         (prlab(ii)(1:4).eq.'EMAS')) then
!     
          ipos=index(prset(ii),' ')
          elset='                    '
          elset(1:ipos-1)=prset(ii)(1:ipos-1)

          write(filename,1000) trim(prlab(ii)(1:5)),
     &                         elset(1:ipos-2),ttime+time
     

!     
!     printing the header
!

111       format('# internal energy (element, energy) for set ',A,
     &           ' and time ',e14.7)
1111      format('element,energy')

112       format('# kinetic energy (element, energy) for set '
     &           ,A,' and time ',e14.7)
1121      format('element,energy')

113       format('# volume (element, volume) for set ',A,
     &           ' and time ',e14.7)
1131      format('element,volume')

114       format('# mass (element, mass) and mass moment of itertia '
     &,          '(xx,yy,zz,xy,xz,yz) for set ',A,
     &           ' and time ',e14.7)
1141      format('element,mass,xx,yy,zz,xy,xz,yz')

115       format('# centrifugal force(element,omega square) for '  
     &           'set ',A,' and time ',e14.7)
1151      format('element,omega^2')


          if((prlab(ii)(1:5).eq.'ELSE ').or.
     &         (prlab(ii)(1:5).eq.'ELSET')) then
            open(unit=else_unit,file=filename, form='formatted')
            write(else_unit,111) elset(1:ipos-2),ttime+time
            write(else_unit,1111)
          elseif((prlab(ii)(1:5).eq.'ELKE ').or.
     &           (prlab(ii)(1:5).eq.'ELKET')) then
            open(unit=elke_unit,file=filename, form='formatted')
            write(elke_unit,112) elset(1:ipos-2),ttime+time
            write(elke_unit,1121)
          elseif((prlab(ii)(1:5).eq.'EVOL ').or.
     &           (prlab(ii)(1:5).eq.'EVOLT')) then
            open(unit=evol_unit,file=filename, form='formatted')
            write(evol_unit,113) elset(1:ipos-2),ttime+time
            write(evol_unit,1131)
          elseif((prlab(ii)(1:5).eq.'EMAS ').or.
     &           (prlab(ii)(1:5).eq.'EMAST')) then
            open(unit=emas_unit,file=filename, form='formatted')
            write(emas_unit,114) elset(1:ipos-2),ttime+time
            write(emas_unit,1141)
          elseif(prlab(ii)(1:4).eq.'CENT') then
            open(unit=cent_unit,file=filename, form='formatted')
            write(emas_unit,115) elset(1:ipos-2),ttime+time
            write(emas_unit,1151)
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
              call printoutelemcsv(prlab,ipkon,lakon,kon,co,
     &               ener,mi(1),ii,nelem,energytot,volumetot,
     &               enerkintot,ne,stx,nodes,thicke,ielmat,
     &               ielem,iface,mortar,ielprop,prop,
     &               sideload,nload,nelemload,xload,bhetot,
     &               xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &               ntmat_,t1,vold,ipobody,ibody,xbody,nbody,
     &               else_unit,elke_unit,evol_unit,emas_unit,
     &               cent_unit)
            elseif(ialset(jj+1).gt.0) then
              nelem=ialset(jj)
              call printoutelemcsv(prlab,ipkon,lakon,kon,co,
     &               ener,mi(1),ii,nelem,energytot,volumetot,
     &               enerkintot,ne,stx,nodes,thicke,ielmat,
     &               ielem,iface,mortar,ielprop,prop,
     &               sideload,nload,nelemload,xload,bhetot,
     &               xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &               ntmat_,t1,vold,ipobody,ibody,xbody,nbody,
     &               else_unit,elke_unit,evol_unit,emas_unit,
     &               cent_unit)
            else
              do nelem=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &               -ialset(jj+1)
                call printoutelemcsv(prlab,ipkon,lakon,kon,co,
     &                 ener,mi(1),ii,nelem,energytot,volumetot,
     &                 enerkintot,ne,stx,nodes,thicke,ielmat,
     &                 ielem,iface,mortar,ielprop,prop,
     &                 sideload,nload,nelemload,xload,bhetot,
     &                 xmasstot,xinertot,cg,ithermal,rhcon,nrhcon,
     &                 ntmat_,t1,vold,ipobody,ibody,xbody,nbody,
     &                 else_unit,elke_unit,evol_unit,emas_unit,
     &                 cent_unit)
              enddo
            endif
          enddo
          
        endif
      enddo
      
      close(u_unit)
      close(f_unit)
      close(s_unit)
      close(e_unit)
      close(pe_unit)
      close(ener_unit)
      close(svd_unit)
      close(me_unit)
      close(coor_unit)
      close(else_unit)
      close(elke_unit)
      close(evol_unit)
      close(emas_unit)
      close(cent_unit)
      return
      end
