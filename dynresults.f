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
      subroutine dynresults(nk,v,ithermal,nactdof,vold,nodeboun,
     &  ndirboun,xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
     &  b,bp,veold,dtime,mi,imdnode,nmdnode,imdboun,nmdboun,imdmpc,
     &  nmdmpc,nmethod,time)
!
!     calculates the displacements or temperatures in a modal dynamics 
!     calculation
!
      implicit none
!
      character*20 labmpc(*)
!
      integer nodeboun(*),ndirboun(*),ipompc(*),imdnode(*),nmdnode,
     &  nodempc(3,*),nk,ithermal(*),i,j,index,mi(*),nactdof(0:mi(2),*),
     &  nboun,nmpc,ist,ndir,node,incrementalmpc,jmin,jmax,
     &  imdboun(*),nmdboun,imdmpc(*),nmdmpc,nmethod
!
      real*8 v(0:mi(2),*),vold(0:mi(2),*),xboun(*),coefmpc(*),
     &  fixed_disp,b(*),veold(0:mi(2),*),bp(*),dtime,time,omega
!
!     omega is needed to calculated the velocity in steady state
!     dynamics calculations
!
      omega=8.d0*datan(1.d0)*time
!
      if(ithermal(1).le.1) then
         jmin=1
         jmax=3
      elseif(ithermal(1).eq.2) then
         jmin=0
         jmax=min(2,mi(2))
      else
         jmin=0
         jmax=3
      endif
!     
!     output for all nodes 
!     
      if(nmdnode.eq.0) then
!
!     extracting the displacement information from the solution
!
         do i=1,nk
            do j=jmin,jmax
               if(nactdof(j,i).gt.0) then
                  v(j,i)=b(nactdof(j,i))
c                  vold(j,i)=b(nactdof(j,i))
               else
                  v(j,i)=0.d0
c                  vold(j,i)=0.d0
               endif
            enddo
         enddo
!
!     inserting the boundary conditions
!     
         do i=1,nboun
            if(ndirboun(i).gt.3) cycle
            fixed_disp=xboun(i)
            v(ndirboun(i),nodeboun(i))=fixed_disp
c            vold(ndirboun(i),nodeboun(i))=fixed_disp
         enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
         do i=1,nmpc
            if((labmpc(i)(1:20).eq.'                    ').or.
     &           (labmpc(i)(1:6).eq.'CYCLIC').or.
     &           (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
               incrementalmpc=0
            else
               incrementalmpc=1
            endif
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.eq.0) then
               if(ithermal(1).lt.2) cycle
            else
               if(ithermal(1).eq.2) cycle
            endif
            index=nodempc(3,ist)
            fixed_disp=0.d0
            if(index.ne.0) then
               do
                  if(incrementalmpc.eq.0) then
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    v(nodempc(2,index),nodempc(1,index))
                  else
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    (v(nodempc(2,index),nodempc(1,index))-
     &                    vold(nodempc(2,index),nodempc(1,index)))
                  endif
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
            fixed_disp=fixed_disp/coefmpc(ist)
            if(incrementalmpc.eq.1) then
               fixed_disp=fixed_disp+vold(ndir,node)
            endif
            v(ndir,node)=fixed_disp
c            vold(ndir,node)=fixed_disp
         enddo
!     
!     extracting the velocity information from the solution
!     
         if(nmethod.eq.4) then
            do i=1,nk
               do j=jmin,jmax
                  if(nactdof(j,i).gt.0) then
                     veold(j,i)=bp(nactdof(j,i))
                  else
                     veold(j,i)=0.d0
                  endif
               enddo
            enddo
!     
!     inserting the boundary conditions
!     
            do i=1,nboun
               if(ndirboun(i).gt.3) cycle
               fixed_disp=xboun(i)
               veold(ndirboun(i),nodeboun(i))=
     &              (fixed_disp-vold(ndirboun(i),nodeboun(i)))/dtime
            enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
            do i=1,nmpc
               if((labmpc(i)(1:20).eq.'                    ').or.
     &              (labmpc(i)(1:6).eq.'CYCLIC').or.
     &              (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
                  incrementalmpc=0
               else
                  incrementalmpc=1
               endif
               ist=ipompc(i)
               node=nodempc(1,ist)
               ndir=nodempc(2,ist)
               if(ndir.eq.0) then
                  if(ithermal(1).lt.2) cycle
               else
                  if(ithermal(1).eq.2) cycle
               endif
               index=nodempc(3,ist)
               fixed_disp=0.d0
               if(index.ne.0) then
                  do
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    veold(nodempc(2,index),nodempc(1,index))
                     index=nodempc(3,index)
                     if(index.eq.0) exit
                  enddo
               endif
               veold(ndir,node)=fixed_disp/coefmpc(ist)
            enddo
!     
!     extracting the velocity information from the solution
!     
         elseif(nmethod.eq.5) then
            do i=1,nk
               do j=jmin,jmax
                  if(nactdof(j,i).gt.0) then
                     veold(j,i)=bp(nactdof(j,i))*omega
                  else
                     veold(j,i)=0.d0
                  endif
               enddo
            enddo
!     
!     inserting the boundary conditions
!     
            do i=1,nboun
               if(ndirboun(i).gt.3) cycle
               veold(ndirboun(i),nodeboun(i))=xboun(i)*omega
            enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
            do i=1,nmpc
               if((labmpc(i)(1:20).eq.'                    ').or.
     &              (labmpc(i)(1:6).eq.'CYCLIC').or.
     &              (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
                  incrementalmpc=0
               else
                  incrementalmpc=1
               endif
               ist=ipompc(i)
               node=nodempc(1,ist)
               ndir=nodempc(2,ist)
               if(ndir.eq.0) then
                  if(ithermal(1).lt.2) cycle
               else
                  if(ithermal(1).eq.2) cycle
               endif
               index=nodempc(3,ist)
               fixed_disp=0.d0
               if(index.ne.0) then
                  do
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    veold(nodempc(2,index),nodempc(1,index))
                     index=nodempc(3,index)
                     if(index.eq.0) exit
                  enddo
               endif
               veold(ndir,node)=fixed_disp/coefmpc(ist)
            enddo
         endif
!
c         do i=1,nk
c            do j=jmin,jmax
c               vold(j,i)=v(j,i)
c            enddo
c         enddo
!     
!     output for a selected number of nodes (fields imdnode,
!     imdboun and imdmpc)  
!     
      else
!     
!     extracting the displacement information from the solution
!     
         do i=1,nmdnode
            do j=jmin,jmax
               if(nactdof(j,imdnode(i)).gt.0) then
                  v(j,imdnode(i))=b(nactdof(j,imdnode(i)))
c                  vold(j,imdnode(i))=b(nactdof(j,imdnode(i)))
               else
                  v(j,imdnode(i))=0.d0
c                  vold(j,imdnode(i))=0.d0
               endif
            enddo
         enddo
!     
!     inserting the boundary conditions
!     
         do j=1,nmdboun
            i=imdboun(j)
            if(ndirboun(i).gt.3) cycle
            fixed_disp=xboun(i)
            v(ndirboun(i),nodeboun(i))=fixed_disp
c            vold(ndirboun(i),nodeboun(i))=fixed_disp
         enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
         do j=1,nmdmpc
            i=imdmpc(j)
            if((labmpc(i)(1:20).eq.'                    ').or.
     &           (labmpc(i)(1:6).eq.'CYCLIC').or.
     &           (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
               incrementalmpc=0
            else
               incrementalmpc=1
            endif
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.eq.0) then
               if(ithermal(1).lt.2) cycle
            else
               if(ithermal(1).eq.2) cycle
            endif
            index=nodempc(3,ist)
            fixed_disp=0.d0
            if(index.ne.0) then
               do
                  if(incrementalmpc.eq.0) then
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    v(nodempc(2,index),nodempc(1,index))
                  else
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    (v(nodempc(2,index),nodempc(1,index))-
     &                    vold(nodempc(2,index),nodempc(1,index)))
                  endif
                  index=nodempc(3,index)
                  if(index.eq.0) exit
               enddo
            endif
            fixed_disp=fixed_disp/coefmpc(ist)
            if(incrementalmpc.eq.1) then
               fixed_disp=fixed_disp+vold(ndir,node)
            endif
            v(ndir,node)=fixed_disp
c            vold(ndir,node)=fixed_disp
         enddo
!     
!     extracting the velocity information from the solution
!  
         if(nmethod.eq.4) then
            do i=1,nmdnode
               do j=jmin,jmax
                  if(nactdof(j,imdnode(i)).gt.0) then
                     veold(j,imdnode(i))=bp(nactdof(j,imdnode(i)))
                  else
                     veold(j,imdnode(i))=0.d0
                  endif
               enddo
            enddo
!     
!     inserting the boundary conditions
!     
            do j=1,nmdboun
               i=imdboun(j)
               if(ndirboun(i).gt.3) cycle
               fixed_disp=xboun(i)
               veold(ndirboun(i),nodeboun(i))=
     &              (fixed_disp-vold(ndirboun(i),nodeboun(i)))/dtime
            enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
            do j=1,nmdmpc
               i=imdmpc(j)
               if((labmpc(i)(1:20).eq.'                    ').or.
     &              (labmpc(i)(1:6).eq.'CYCLIC').or.
     &              (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
                  incrementalmpc=0
               else
                  incrementalmpc=1
               endif
               ist=ipompc(i)
               node=nodempc(1,ist)
               ndir=nodempc(2,ist)
               if(ndir.eq.0) then
                  if(ithermal(1).lt.2) cycle
               else
                  if(ithermal(1).eq.2) cycle
               endif
               index=nodempc(3,ist)
               fixed_disp=0.d0
               if(index.ne.0) then
                  do
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    veold(nodempc(2,index),nodempc(1,index))
                     index=nodempc(3,index)
                     if(index.eq.0) exit
                  enddo
               endif
               veold(ndir,node)=fixed_disp/coefmpc(ist)
            enddo
         elseif(nmethod.eq.5) then
            do i=1,nmdnode
               do j=jmin,jmax
                  if(nactdof(j,imdnode(i)).gt.0) then
                     veold(j,imdnode(i))=bp(nactdof(j,imdnode(i)))*omega
                  else
                     veold(j,imdnode(i))=0.d0
                  endif
               enddo
            enddo
!     
!     inserting the boundary conditions
!     
            do j=1,nmdboun
               i=imdboun(j)
               if(ndirboun(i).gt.3) cycle
               veold(ndirboun(i),nodeboun(i))=xboun(i)*omega
            enddo
!     
!     inserting the mpc information
!     the parameter incrementalmpc indicates whether the
!     incremental displacements enter the mpc or the total 
!     displacements (incrementalmpc=0)
!     
            do j=1,nmdmpc
               i=imdmpc(j)
               if((labmpc(i)(1:20).eq.'                    ').or.
     &              (labmpc(i)(1:6).eq.'CYCLIC').or.
     &              (labmpc(i)(1:9).eq.'SUBCYCLIC')) then
                  incrementalmpc=0
               else
                  incrementalmpc=1
               endif
               ist=ipompc(i)
               node=nodempc(1,ist)
               ndir=nodempc(2,ist)
               if(ndir.eq.0) then
                  if(ithermal(1).lt.2) cycle
               else
                  if(ithermal(1).eq.2) cycle
               endif
               index=nodempc(3,ist)
               fixed_disp=0.d0
               if(index.ne.0) then
                  do
                     fixed_disp=fixed_disp-coefmpc(index)*
     &                    veold(nodempc(2,index),nodempc(1,index))
                     index=nodempc(3,index)
                     if(index.eq.0) exit
                  enddo
               endif
               veold(ndir,node)=fixed_disp/coefmpc(ist)
            enddo
         endif
!
c         do i=1,nmdnode
c            do j=jmin,jmax
c               vold(j,imdnode(i))=v(j,imdnode(i))
c            enddo
c         enddo
      endif
!     
      return
      end
      
