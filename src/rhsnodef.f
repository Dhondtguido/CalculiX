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
      subroutine rhsnodef(co,kon,ne,
     &  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &  nforc,fext,nactdof,nmethod,ikmpc,ntmat_,iperturb,
     &  mi,ikactmech,nactmech,ntrans,inotr,trab,fnext)
!
!     filling the right hand side load vector b
!
!     b contains the contributions due to mechanical forces only
!
      implicit none
!
      integer kon(*),ipompc(*),nodempc(3,*),nodeforc(2,*),ndirforc(*),
     &     ikmpc(*),mi(*),nactdof(0:mi(2),*),konl(20),rhsi,
     &     ntrans,inotr(2,*),ne,nmpc,nforc,nmethod,m,idm,idir,itr,
     &     iperturb(*),i,j,k,idist,jj,node,id,ist,index,jdof1,jdof,
     &     node1,ntmat_,icalccg,ikactmech(*),nactmech
!
      real*8 co(3,*),coefmpc(*),xforc(*),fnext(0:mi(2),*),
     &  fext(*),a(3,3),trab(7,*)
!
      icalccg=0
!
!     point forces
!
!     modal dynamics and steady state dynamics: 
!     location of nonzeros is stored
!
      if((nmethod.ge.4).and.(iperturb(1).lt.2)) then
         do i=1,nforc
            if(ndirforc(i).gt.3) cycle
            if(dabs(xforc(i)).lt.1.d-30) cycle
! 
!           checking for local coordinate system
!    
            node=nodeforc(1,i)
            if(ntrans.eq.0) then
               itr=0
            else
               itr=inotr(1,node)
            endif
!     
            if(itr.eq.0) then
!
!              no local coordinate system
!
               jdof=nactdof(ndirforc(i),node)
               if(jdof.gt.0) then
                  fext(jdof)=fext(jdof)+xforc(i)
                  call nident(ikactmech,jdof-1,nactmech,
     &                 idm)
                  do
                     if(idm.gt.0) then
                        if(ikactmech(idm).eq.jdof-1) exit
                     endif
                     nactmech=nactmech+1
                     do m=nactmech,idm+2,-1
                        ikactmech(m)=ikactmech(m-1)
                     enddo
                     ikactmech(idm+1)=jdof-1
                     exit
                  enddo
               else
!     
!                 node is a dependent node of a MPC: distribute
!                 the forces among the independent nodes
!                 (proportional to their coefficients)
!     
                  jdof=8*(node-1)+ndirforc(i)
                  call nident(ikmpc,jdof,nmpc,id)
                  if(id.gt.0) then
                     if(ikmpc(id).eq.jdof) then
                        ist=ipompc(id)
                        index=nodempc(3,ist)
                        if(index.eq.0) cycle
                        do
                           jdof=nactdof(nodempc(2,index),
     &                                  nodempc(1,index))
                           if(jdof.gt.0) then
                              fext(jdof)=fext(jdof)-
     &                             coefmpc(index)*xforc(i)/coefmpc(ist)
                              call nident(ikactmech,jdof-1,nactmech,
     &                             idm)
                              do
                                 if(idm.gt.0) then
                                    if(ikactmech(idm).eq.jdof-1) exit
                                 endif
                                 nactmech=nactmech+1
                                 do m=nactmech,idm+2,-1
                                    ikactmech(m)=ikactmech(m-1)
                                 enddo
                                 ikactmech(idm+1)=jdof-1
                                 exit
                              enddo
                           endif
                           index=nodempc(3,index)
                           if(index.eq.0) exit
                        enddo
                     endif
                  endif
               endif
            else
!
!              local coordinate system
!
               call transformatrix(trab(1,itr),co(1,node),a)
               idir=ndirforc(i)
!
!              loop over all dofs
!
               do j=1,3
                  jdof=nactdof(j,node)
                  if(jdof.gt.0) then
                     fext(jdof)=fext(jdof)+a(j,idir)*xforc(i)
                     call nident(ikactmech,jdof-1,nactmech,
     &                    idm)
                     do
                        if(idm.gt.0) then
                           if(ikactmech(idm).eq.jdof-1) exit
                        endif
                        nactmech=nactmech+1
                        do m=nactmech,idm+2,-1
                           ikactmech(m)=ikactmech(m-1)
                        enddo
                        ikactmech(idm+1)=jdof-1
                        exit
                     enddo
                  else
!     
!                    node is a dependent node of a MPC: distribute
!                    the forces among the independent nodes
!                    (proportional to their coefficients)
!     
                     jdof=8*(node-1)+j
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
                                 fext(jdof)=fext(jdof)-
     &                                coefmpc(index)*a(j,idir)*xforc(i)/
     &                                coefmpc(ist)
                                 call nident(ikactmech,jdof-1,nactmech,
     &                                idm)
                                 do
                                    if(idm.gt.0) then
                                       if(ikactmech(idm).eq.jdof-1) exit
                                    endif
                                    nactmech=nactmech+1
                                    do m=nactmech,idm+2,-1
                                       ikactmech(m)=ikactmech(m-1)
                                    enddo
                                    ikactmech(idm+1)=jdof-1
                                    exit
                                 enddo
                              endif
                              index=nodempc(3,index)
                              if(index.eq.0) exit
                           enddo
                        endif
                     endif
                  endif
               enddo


            endif
         enddo
      else
!
!       other procedures
!
         rhsi=1
         call mafillsmforc(nforc,ndirforc,nodeforc,xforc,nactdof,
     &        fext,ipompc,nodempc,coefmpc,mi,rhsi,fnext,nmethod,
     &        ntrans,inotr,trab,co)
!         
      endif
!
      return
      end
