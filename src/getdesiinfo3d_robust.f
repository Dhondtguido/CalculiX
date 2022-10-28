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
      subroutine getdesiinfo3d_robust(set,istartset,iendset,ialset,
     &  nset,mi,nactdof,ndesi,nodedesi,ntie,tieset,itmp,nmpc,nodempc,
     &  ipompc,nodedesiinv,iponoel,inoel,lakon,ipkon,kon,iregion,
     &  ipoface,nodface,nk,irandomtype,jobnamef)    
!
!     storing the design variables in nodedesi
!     marking which nodes are design variables in nodedesiinv
!
!     a node is a design variable for the random field if:
!     1) it belongs to the design variable set AND
!     2) not all dofs in the node are defined by SPC's AND
!     3) no MPC is applied to any of its dofs AND
!     4) it belongs to at least one face whose number of
!        design variables exceeds half its nodes
!
      implicit none
!
      character*8 lakon(*)
!
      character*81 setname
      character*81 set(*)
      character*81 tieset(3,*)
      character*132 jobnamef
      character*256 fn
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ndesi,
     &  node,nodedesi(*),nset,ntie,i,j,k,l,m,nmpc,nodempc(3,*),
     &  nactdof(0:mi(2),*),itmp(*),ntmp,index,id,ipompc(*),
     &  nodedesiinv(*),iponoel(*),inoel(2,*),nelem,nope,nopedesi,
     &  ipkon(*),nnodes,kon(*),iregion,konl(26),iaux,kflag,
     &  ipoface(*),nodface(5,*),jfacem,nopesurf(9),ifaceq(8,6),
     &  ifacet(6,4),ifacew1(4,5),ifacew2(8,5),nopem,nk,ishift,
     &  indexe,jface,konl2d(26),expandhex(20),expandwed(15),
     &  idelta,idesi,irandomtype(*),iwrite,ilen
!
      setname(1:1)=' '
      ndesi=0
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
!     opening a file to store the nodes which are rejected as
!     design variables
!
      iwrite=0
      ilen=index(jobnamef,' ')-1
      fn=jobnamef(1:ilen)//'_WarnNodeDesignReject.nam'
      open(40,file=fn,status='unknown')
      write(40,*) '*NSET,NSET=WarnNodeDesignReject'
!
!------------------------------------------------------------------------
!     assign the nodes of the set to the appropriate variables
!------------------------------------------------------------------------
!
      do i=1,nk
         if(irandomtype(i).gt.0) then
            ndesi=ndesi+1
            nodedesi(ndesi)=i
            nodedesiinv(i)=-2      
         endif
      enddo 
!     
!------------------------------------------------------------------------
!     Remove interior nodes from the set of design variables    
!------------------------------------------------------------------------
!
      do j=1,nk
!        
         if(ipoface(j).eq.0) cycle
         indexe=ipoface(j)
!        
         do
!
            nelem=nodface(3,indexe)
            jface=nodface(4,indexe)
!
            if((lakon(nelem)(7:7).eq.'A').or.
     &         (lakon(nelem)(7:7).eq.'S').or.      
     &         (lakon(nelem)(7:7).eq.'E')) then
!
!              for plane stress/strain/axi only faces 
!              3 and higher are taken into account for the normal
!
               if(jface.le.2) then
                  indexe=nodface(5,indexe)
                  if(indexe.eq.0) then
                     exit
                  else
                     cycle
                  endif
               endif
            elseif(lakon(nelem)(7:7).eq.'L') then
!
!              for shells only faces 2 and lower
!              taken into account for the normal
!
               if(jface.gt.2) then
                  indexe=nodface(5,indexe)
                  if(indexe.eq.0) then
                     exit
                  else
                     cycle
                  endif
               endif
            endif
!
!           nopem: # of nodes in the surface
!           nope: # of nodes in the element
!     
            if(lakon(nelem)(4:4).eq.'8') then
               nopem=4
               nope=8
               nopedesi=3
               ishift=8
            elseif(lakon(nelem)(4:5).eq.'20') then
               nopem=8
               nope=20
               nopedesi=5
               ishift=20
            elseif(lakon(nelem)(4:5).eq.'10') then
               nopem=6
               nope=10
               nopedesi=4
            elseif(lakon(nelem)(4:4).eq.'4') then
               nopem=3
               nope=4
               nopedesi=3
!     
!     treatment of wedge faces
!     
            elseif(lakon(nelem)(4:4).eq.'6') then
               nope=6
               if(jface.le.2) then
                  nopem=3
               else
                  nopem=4
               endif
               nopedesi=3
               ishift=6
            elseif(lakon(nelem)(4:5).eq.'15') then
               nope=15
               if(jface.le.2) then
                  nopem=6
                  nopedesi=4
               else
                  nopem=8
                  nopedesi=5
               endif
               ishift=15
            endif
            if(iregion.eq.0) then
         nopedesi=0
      endif
!     
!     actual position of the nodes belonging to the surface
!     
            if((lakon(nelem)(7:7).eq.'A').or.
     &         (lakon(nelem)(7:7).eq.'S').or.      
     &         (lakon(nelem)(7:7).eq.'E')) then
               if((lakon(nelem)(4:5).eq.'20').or.
     &            (lakon(nelem)(4:5).eq.'8 ')) then
                  do k=1,nope
                     konl(k)=kon(ipkon(nelem)+k)
                     konl2d(k)=kon(ipkon(nelem)+ishift+expandhex(k))
                  enddo
               elseif((lakon(nelem)(4:5).eq.'15').or.
     &                 (lakon(nelem)(4:5).eq.'6 ')) then
                  do k=1,nope
                     konl(k)=kon(ipkon(nelem)+k)
                     konl2d(k)=kon(ipkon(nelem)+ishift+expandwed(k))
                  enddo
               endif
            else
               do k=1,nope
                  konl(k)=kon(ipkon(nelem)+k)
               enddo
            endif
!     
            if((nope.eq.20).or.(nope.eq.8)) then
               do m=1,nopem
                  nopesurf(m)=konl(ifaceq(m,jface))
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then 
               do m=1,nopem
                  nopesurf(m)=konl(ifacet(m,jface))
               enddo
            elseif(nope.eq.15) then
               do m=1,nopem
                  nopesurf(m)=konl(ifacew2(m,jface))
               enddo
            else
               do m=1,nopem
                  nopesurf(m)=konl(ifacew1(m,jface))
               enddo
            endif
!    
!     sum up how many designvariables are on that surface
!
            nnodes=0
            do m=1,nopem
               if(nodedesiinv(nopesurf(m)).ne.0) then
                  nnodes=nnodes+1
               endif
            enddo
!     
!     set design variables located on external surfaces 
!     equal to 1 in nodedesiinv
!     if sufficient design variables are defined on this surface
!     
            nnodes=0
            do m=1,nopem
               if((nodedesiinv(nopesurf(m)).lt.0).or.
     &            (nodedesiinv(nopesurf(m)).eq.1)) then
                  nnodes=nnodes+1
               endif
            enddo               
!
            if(nnodes.ge.nopedesi) then
               do m=1,nopem
                  if(nodedesiinv(nopesurf(m)).lt.0) then
                     nodedesiinv(nopesurf(m))=1
                  endif
               enddo
      else
               do m=1,nopem
                  if(nodedesiinv(nopesurf(m)).eq.-2) then
                     nodedesiinv(nopesurf(m))=-1
                  endif
               enddo
            endif

            indexe=nodface(5,indexe)
            if(indexe.eq.0) exit
!            
         enddo
      enddo
!
      kflag=1
      call isortii(nodedesi,iaux,ndesi,kflag)
!
!------------------------------------------------------------------------
!     if node i in nodedesiinv(i) is -2 or -1 --> delete node i from 
!     set of designvariables
!------------------------------------------------------------------------
!
      do i=1,nk
         if(nodedesiinv(i).eq.-2) then
!
            write(*,*) '*WARNING in getdesiinfo:'
            write(*,*) '          node ',i,' is removed'
            write(*,*) '          from the set of design'
            write(*,*) '          variables as it is an' 
            write(*,*) '          internal node'
            write(40,*) i
            iwrite=1
!
         elseif(nodedesiinv(i).eq.-1) then
!
            write(*,*) '*WARNING in getdesiinfo:'
            write(*,*) '          node ',i,' is removed'
            write(*,*) '          from the set of design'
            write(*,*) '          variables as not sufficient '
            write(*,*) '          other variables are on the  '
            write(*,*) '          surrounding element faces  '
            write(40,*) i
            iwrite=1
!
            nodedesiinv(i)=0
            call nident(nodedesi,i,ndesi,id)
            do k=id+1,ndesi
               nodedesi(k-1)=nodedesi(k)
            enddo
            ndesi=ndesi-1    
         endif
      enddo 
!
!------------------------------------------------------------------------
!     write remaining design variables in the dat file
!------------------------------------------------------------------------
!
      write(5,*)'NUMBER OF DESIGN VARIABLES TAKEN INTO ACCOUNT '
      write(5,*)'FOR RANDOM FIELD ASSESSMENT:'
      write(5,*) ndesi
      write(5,*)''
      write(5,*)'NODE NUMBERS OF DESIGN VARIABLES: '
      do i=1,ndesi
         write(5,*) nodedesi(i),','
      enddo
!
      if(iwrite.eq.1) then
        write(*,*) '*INFO in getdesiinfo3d_robust:'
        write(*,*) '      rejected design nodes are stored in'
        write(*,*) '      file ',fn(1:ilen+25)
        write(*,*) '      This file can be loaded into'
        write(*,*) '      an active cgx-session by typing'
        write(*,*) 
     &       '      read ',fn(1:ilen+25),' inp'
        write(*,*)
        close(40)
      else
        close(40,status='delete')
      endif
!
      return
      end

