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
      subroutine readforce(zc,neq,nev,nactdof,ikmpc,nmpc,ipompc,
     &  nodempc,mi,coefmpc,jobnamef,a,igeneralizedforce)
!
!     reads a complex force (e.g. response of a fluid to a harmonic
!     structural excitation)
!     the force has to be stored in file 'dummy' in the form:
!
!     node,Fx-real,Fx-imag,Fy-real,Fy-imag,Fz-real,Fz-imag
!
!     only nodes in which a force is applied have to be stored. Modes
!     have to be separated by line starting with **
!
!     for cyclic symmetric structures the eigenmodes come in pairs
!     forces must be given for the first mode of each pair only
!
      implicit none
!
      logical exi
!
      character*132 jobnamef
      character*144 name
!
      integer mi(*),neq,nev,i,j,k,nactdof(0:mi(2),*),ikmpc(*),nmpc,
     &  jdof,id,ist,ipompc(*),index,nodempc(3,*),node,istat,
     &  igeneralizedforce
!
      real*8 coefmpc(*),comp(6)
!
      complex*16 zc(neq,*),force(3),a(nev,*)
!
      igeneralizedforce=0
!
!     creating name for force file
!
      do i=1,132
         if(jobnamef(i:i).eq.' ') exit
         name(i:i)=jobnamef(i:i)
      enddo
      i=i-1
      name(i+1:i+6)='_force'
      do j=i+7,144
         name(j:j)=' '
      enddo
!
!     if a force file exists, it is read. If it does not exist,
!     a generalized force file is looked for (= product of the
!     force due to eigenmode i with eigenmode j, leading to the
!     nev x nev a-matrix).
!
      inquire(file=name,exist=exi)
!
      if(exi) then
         open(27,file=name,status='unknown')
!     
         do i=1,nev
            do
               read(27,*,iostat=istat) node,(comp(k),k=1,6)
               if(istat.ne.0) then
                  exit
               endif
!     
               do k=1,3
                  force(k)=comp(2*k-1)*(1.d0,0.d0)+comp(2*k)*(0.d0,1.d0)
               enddo
!     
               do k=1,3
                  jdof=nactdof(k,node)
                  if(jdof.gt.0) then
                     zc(jdof,i)=zc(jdof,i)-force(k)
                  else
!     
!     node is a dependent node of a MPC: distribute
!     the forces among the independent nodes
!     (proportional to their coefficients)
!     
                     jdof=8*(node-1)+k
                     call nident(ikmpc,jdof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.jdof) then
                           ist=ipompc(id)
                           index=nodempc(3,ist)
                           if(index.eq.0) cycle
                           do
                              jdof=nactdof(nodempc(2,index),
     &                             nodempc(1,index))
                              if(jdof.gt.0) then
                                 zc(jdof,i)=zc(jdof,i)-
     &                              coefmpc(index)*force(k)/coefmpc(ist)
                              endif
                              index=nodempc(3,index)
                              if(index.eq.0) exit
                           enddo
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo
!     
         close(27)
      else
!
!     creating name for generalized force file
!
         do i=1,132
            if(jobnamef(i:i).eq.' ') exit
            name(i:i)=jobnamef(i:i)
         enddo
         i=i-1
         name(i+1:i+9)='_genforce'
         do j=i+10,144
            name(j:j)=' '
         enddo
!
         inquire(file=name,exist=exi)
!
         if(exi) then
!
            igeneralizedforce=1
!
            open(27,file=name,status='unknown')
            do
               read(27,*,iostat=istat)i,j,a(i,j)
               if(istat.ne.0) exit
            enddo
            close(27)
         else
            write(*,*) '*ERROR in readforce: neither a force file'
            write(*,*) '       nor a generalized force file exists'
            call exit(201)
         endif
      endif
!     
      return
      end
      
