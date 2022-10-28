!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine printoutfluidfem(set,nset,istartset,iendset,ialset,
     &     nprint,prlab,prset,v,ipkon,lakon,sti,
     &     mi,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &     orab,ielorien,norien,vold,ielmat,
     &     thicke,physcon,ielprop,prop,orname,vcon,nk,nknew,
     &     nelnew)
!     
!     stores results in the .dat file
!     
      implicit none
!     
      character*6 prlab(*)
      character*8 lakon(*)
      character*80 noset,elset,orname(*)
      character*81 set(*),prset(*)
!     
      integer nset,istartset(*),iendset(*),ialset(*),nprint,ipkon(*),
     &     mi(*),ii,jj,iset,l,limit,node,ipos,ithermal(*),id,
     &     nelem,kon(*),inotr(2,*),ntrans,ielorien(mi(3),*),norien,
     &     mt,ielmat(mi(3),*),nk,ielprop(*),nknew(*),nelnew(*)
!     
      real*8 sti(6,mi(1),*),co(3,*),qfx(3,mi(1),*),ttime,
     &     trab(7,*),orab(7,*),vold(0:mi(2),*),v(nk,0:mi(2)),
     &     thicke(mi(3),*),physcon(*),prop(*),
     &     vcon(nk,0:mi(2))
!     
      mt=mi(2)+1
!     
      do ii=1,nprint
!     
!     nodal values
!     
        if((prlab(ii)(1:4).eq.'VF  ').or.(prlab(ii)(1:4).eq.'PSF ').or.
     &       (prlab(ii)(1:4).eq.'TSF ').or.
     &       (prlab(ii)(1:4).eq.'PTF ').or. 
     &       (prlab(ii)(1:4).eq.'TTF ').or.
     &       (prlab(ii)(1:4).eq.'CP  ').or.
     &       (prlab(ii)(1:4).eq.'TURB').or.
     &       (prlab(ii)(1:4).eq.'DEPF').or.
     &       (prlab(ii)(1:4).eq.'MACH')) then
!     
          ipos=index(prset(ii),' ')
          noset='                    '
          noset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
!     printing the header
!     
          if(prlab(ii)(1:4).eq.'VF  ') then
            write(5,*)
            write(5,100) noset(1:ipos-2),ttime
 100        format(' velocities (vx,vy,vz) for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'PSF ') then
            write(5,*)
            write(5,101) noset(1:ipos-2),ttime
 101        format(' static pressures for set ',A,' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:5).eq.'TSF ') then
            write(5,*)
            write(5,102) noset(1:ipos-2),ttime
 102        format(' static temperatures for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:5).eq.'PTF ') then
            write(5,*)
            write(5,103) noset(1:ipos-2),ttime
 103        format(' total pressures for set ',A,' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'TTF ') then
            write(5,*)
            write(5,115) noset(1:ipos-2),ttime
 115        format(' total temperatures for set ',A,' and time ',
     &           e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'CP  ') then
            write(5,*)
            write(5,117) noset(1:ipos-2),ttime
 117        format(' pressure coefficients for set ',
     &           A,' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'TURB') then
            write(5,*)
            write(5,118) noset(1:ipos-2),ttime
 118        format(' turbulence conservative variables for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'MACH') then
            write(5,*)
            write(5,119) noset(1:ipos-2),ttime
 119        format(' Mach number for set ',A,
     &           ' and time ',e14.7)
            write(5,*)
          elseif(prlab(ii)(1:4).eq.'DEPF') then
            write(5,*)
            write(5,120) noset(1:ipos-2),ttime
 120        format(' Change in depth for set ',A,
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
          do jj=istartset(iset),iendset(iset)
            if(ialset(jj).lt.0) cycle
            if(jj.eq.iendset(iset)) then
              node=ialset(jj)
              call printoutnodefluidfem(prlab,v,vold,physcon,
     &             ii,node,trab,inotr,ntrans,co,mi,vcon,nk,nknew)
            elseif(ialset(jj+1).gt.0) then
              node=ialset(jj)
              call printoutnodefluidfem(prlab,v,vold,physcon,
     &             ii,node,trab,inotr,ntrans,co,mi,vcon,nk,nknew)
            else
              do node=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &             -ialset(jj+1)
                call printoutnodefluidfem(prlab,v,vold,physcon,
     &               ii,node,trab,inotr,ntrans,co,mi,vcon,nk,nknew)
              enddo
            endif
          enddo
!     
!     integration point values
!     
        elseif((prlab(ii)(1:4).eq.'SVF ').or.
     &         (prlab(ii)(1:4).eq.'HFLF')) then
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
            if(prlab(ii)(1:4).eq.'SVF ') then
              write(5,*)
              write(5,106) elset(1:ipos-2),ttime
 106          format(' viscous stresses (elem, integ.pnt.,sxx,syy,sz
     &z,sxy,sxz,syz) for set ',A,' and time ',e14.7)
              write(5,*)
            elseif(prlab(ii)(1:4).eq.'HFLF') then
              write(5,*)
              write(5,112) elset(1:ipos-2),ttime
 112          format(' heat flux (elem, integ.pnt.,qx,qy,qz) for set 
     &',A,' and time ',e14.7)
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
                call printoutintfluidfem(prlab,ipkon,lakon,sti,
     &               mi(1),ii,nelem,qfx,
     &               orab,ielorien,norien,co,kon,ielmat,thicke,
     &               ielprop,prop,nelem,ithermal,orname,nelnew)
              elseif(ialset(jj+1).gt.0) then
                nelem=ialset(jj)
                call printoutintfluidfem(prlab,ipkon,lakon,sti,
     &               mi(1),ii,nelem,qfx,orab,
     &               ielorien,norien,co,kon,ielmat,thicke,
     &               ielprop,prop,nelem,ithermal,orname,nelnew)
              else
                do nelem=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &               -ialset(jj+1)
                  call printoutintfluidfem(prlab,ipkon,lakon,sti,
     &                 mi(1),ii,nelem,
     &                 qfx,orab,ielorien,norien,co,kon,ielmat,
     &                 thicke,ielprop,prop,nelem,ithermal,
     &                 orname,nelnew)
                enddo
              endif
            enddo
!     
          enddo
        endif
      enddo
!     
      return
      end
