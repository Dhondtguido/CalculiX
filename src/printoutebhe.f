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
      subroutine printoutebhe(set,nset,istartset,iendset,ialset,nprint,
     &     prlab,prset,t1,ipkon,lakon,stx,ener,mi,ithermal,co,kon,ttime,
     &     ne,vold,ielmat,thicke,mortar,time,ielprop,prop,
     &     nelemload,nload,sideload,xload,rhcon,nrhcon,ntmat_,ipobody,
     &     ibody,xbody,nbody,nmethod)
!     
!     stores results in the .dat file
!     
      implicit none
!     
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 elset
      character*81 set(*),prset(*)
!     
      integer nset,istartset(*),iendset(*),ialset(*),nprint,ipkon(*),
     &     mi(*),ii,jj,iset,ipos,ithermal(*),ielem,
     &     nelem,kon(*),nodes,mt,ielmat(mi(3),*),iface,
     &     mortar,ielprop(*),nload,ntmat_,id,
     &     nelemload(2,*),nrhcon(*),ipobody(2,*),ibody(3,*),nbody,
     &     nmethod,ne
!     
      real*8 t1(*),stx(6,mi(1),*),bhetot,ener(2,mi(1),*),
     &     volumetot,co(3,*),ttime,time,vold(0:mi(2),*),enerkintot,
     &     prop(*),xload(2,*),xmasstot,xinertot(6),cg(3),
     &     rhcon(0:1,ntmat_,*),xbody(7,*),energytot,thicke(mi(3),*)
!     
      mt=mi(2)+1
!     
      do ii=1,nprint
!     
!     whole element values
!     
        if(prlab(ii)(1:4).eq.'EBHE') then
!     
          ipos=index(prset(ii),' ')
          elset='                    '
          elset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
!     printing the header
!     
          if((prlab(ii)(1:5).eq.'EBHE ').or.
     &           (prlab(ii)(1:5).eq.'EBHET')) then
            write(5,*)
            write(5,131) elset(1:ipos-2),ttime+time
 131        format(' body heating (element, volume) for set ',A,
     &           ' and time ',e14.7)
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
!     
!     writing total values to file
!     
          if((prlab(ii)(1:5).eq.'EBHEO').or.
     &           (prlab(ii)(1:5).eq.'EBHET')) then
            write(5,*)
            write(5,132) elset(1:ipos-2),ttime+time
 132        format(' total body heating for set ',A,' and time ',
     &           e14.7)
            write(5,*)
            write(5,'(6x,1p,1x,e13.6)') bhetot
          endif
        endif
      enddo
!     
      return
      end
