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
      subroutine checkconstraint(nobject,objectset,g0,nactive,
     &     nnlconst,ipoacti,ndesi,dgdxglob,nk,nodedesi,iconstacti,
     &     objnorm,inameacti)               
!     
!     check which constraints are active on the basis of the 
!     function values of the constraints     
!     
      implicit none
!     
      character*81 objectset(5,*)
      character*20 empty
!     
      integer iobject,nobject,istat,nactive,nnlconst,inameacti(*),
     &   ipoacti(*),ifree,i,ndesi,nk,node,nodedesi(*),
     &   iconstacti(*),nconst
!     
      real*8 g0(*),bounds(nobject),scale,bound,objnorm(*),
     &   dgdxglob(2,nk,*)
!
      empty='                    '
!   
!     
!     header written in dat file
!     
      write(5,*)
      write(5,*)
      write(5,'(a113)') '  ################################################
     &#################################################################'
      write(5,*) '  A S S E M B L Y   O F   A C T I V E   S E T'
      write(5,*)
      write(5,101)
     &'NUMBER OF    ','CONSTRAINT      ','LE/     ','FUNCTION         ',
     &'FUNCTION         ','FUNCTION      ','  ACTIVE/ ','   NAME OF' 
      write(5,101)
     &'CONSTRAINT   ','FUNCTION        ','GE      ','VALUE            ',
     &'BOUND            ','VIOLATION     ','  INACTIVE','   CONSTRAINT' 
      write(5,'(a113)') '  ################################################
     &#################################################################'
      write(5,*)
!     
!     determine bounds of constraints
!     
      do iobject=2,nobject
         if(objectset(5,iobject)(81:81).eq.'G') then
            read(objectset(1,iobject)(61:80),'(f20.0)',
     &      iostat=istat) bound
            bounds(iobject)=bound
         elseif(objectset(5,iobject)(81:81).eq.'C') then
            if(objectset(1,iobject)(61:80).ne.empty) then
               read(objectset(1,iobject)(61:80),'(f20.0)',
     &         iostat=istat) bound
            else
               write(*,*) '*WARNING in checkconstraint'
               write(*,*) '         no absolute constraint boundary'
               write(*,*) '         defined, system value taken' 
               write(*,*)
               bound=g0(iobject)
            endif
            if(objectset(1,iobject)(41:60).ne.empty) then
               read(objectset(1,iobject)(41:60),'(f20.0)',
     &         iostat=istat) scale
            else
               write(*,*) '*WARNING in checkconstraint'
               write(*,*) '         no relative constraint boundary'
               write(*,*) '         defined, 1.0 taken' 
               write(*,*)
               scale=1.0d0
            endif
            bounds(iobject)=bound*scale
         endif
      enddo
!     
!     determine all active constraints
!     
      nconst=0
      nactive=0
      nnlconst=0
      ifree=1
!
      do iobject=2,nobject
!     
!        determine all active nonlinear constraints
!     
         if(objectset(5,iobject)(81:81).eq.'C') then
            nconst=nconst+1
            if(objectset(1,iobject)(19:20).eq.'LE') then
               objnorm(ifree)=g0(iobject)/bounds(iobject)-1
               if(objnorm(ifree).ge.0.d0) then
                  nactive=nactive+1
                  nnlconst=nnlconst+1
                  ipoacti(ifree)=iobject
                  inameacti(ifree)=iobject
                  iconstacti(ifree)=-1
                  write(5,102) nconst,objectset(1,iobject),'LE  ',
     &               g0(iobject),bounds(iobject),objnorm(ifree),
     &               'ACTIVE  ',objectset(5,iobject)
                  ifree=ifree+1
               else
                  write(5,102) nconst,objectset(1,iobject),'LE  ',
     &               g0(iobject),bounds(iobject),objnorm(ifree),
     &               'INACTIVE',objectset(5,iobject)            
               endif
            elseif(objectset(1,iobject)(19:20).eq.'GE') then
               objnorm(ifree)=-1*(g0(iobject)/bounds(iobject))+1
               if(objnorm(ifree).ge.0.d0) then
                  nactive=nactive+1
                  nnlconst=nnlconst+1
                  ipoacti(ifree)=iobject
                  inameacti(ifree)=iobject
                  iconstacti(ifree)=1
                  write(5,102) nconst,objectset(1,iobject),'GE  ',
     &               g0(iobject),bounds(iobject),objnorm(ifree),
     &               'ACTIVE  ',objectset(5,iobject)
                  ifree=ifree+1
               else
                  write(5,102) nconst,objectset(1,iobject),'GE  ',
     &               g0(iobject),bounds(iobject),objnorm(ifree),
     &               'INACTIVE',objectset(5,iobject)            
               endif
            endif
!     
!        determine all active linear constraints
!     
         elseif(objectset(5,iobject)(81:81).eq.'G') then
            nconst=nconst+1
            if(objectset(1,iobject)(19:20).eq.'LE') then
               if(g0(iobject).gt.0d0) then
                  do i=1,ndesi
                     node=nodedesi(i)
                     if(i.eq.1) then
                        objnorm(iobject)=dgdxglob(2,node,iobject)
                     else
                        objnorm(iobject)=
     &                  max(objnorm(iobject),dgdxglob(2,node,iobject))
                     endif
                     if(dgdxglob(2,node,iobject).ge.0.d0) then
                        ipoacti(ifree)=i
                        inameacti(ifree)=iobject
                        iconstacti(ifree)=-1
                        ifree=ifree+1
                        nactive=nactive+1
                     endif
                  enddo
                  write(5,102) nconst,objectset(1,iobject),'LE  ',
     &               g0(iobject),bounds(iobject),objnorm(iobject),
     &               'ACTIVE  ',objectset(5,iobject)    
               else
                  do i=1,ndesi
                     node=nodedesi(i)
                     if(i.eq.1) then
                        objnorm(iobject)=dgdxglob(2,node,iobject)
                     else
                        objnorm(iobject)=
     &                  max(objnorm(iobject),dgdxglob(2,node,iobject))
                     endif
                  enddo
                  write(5,102) nconst,objectset(1,iobject),'LE  ',
     &               g0(iobject),bounds(iobject),objnorm(iobject),
     &               'INACTIVE',objectset(5,iobject)    
               endif
            elseif(objectset(1,iobject)(19:20).eq.'GE') then
               if(g0(iobject).gt.0.d0) then
                  do i=1,ndesi
                     node=nodedesi(i)
                     if(i.eq.1) then
                        objnorm(iobject)=dgdxglob(2,node,iobject)
                     else
                        objnorm(iobject)=
     &                  max(objnorm(iobject),dgdxglob(2,node,iobject))
                     endif
                     if(dgdxglob(2,node,iobject).ge.0.d0) then
                        ipoacti(ifree)=i
                        inameacti(ifree)=iobject
                        iconstacti(ifree)=1
                        ifree=ifree+1
                        nactive=nactive+1
                     endif
                  enddo
                  write(5,102) nconst,objectset(1,iobject),'GE  ',
     &               g0(iobject),bounds(iobject),objnorm(iobject),
     &               'ACTIVE  ',objectset(5,iobject)    
               else
                  do i=1,ndesi
                     node=nodedesi(i)
                     if(i.eq.1) then
                        objnorm(iobject)=dgdxglob(2,node,iobject)
                     else
                        objnorm(iobject)=
     &                  max(objnorm(iobject),dgdxglob(2,node,iobject))
                     endif
                  enddo
                  write(5,102) nconst,objectset(1,iobject),'GE  ',
     &               g0(iobject),bounds(iobject),objnorm(iobject),
     &               'INACTIVE',objectset(5,iobject)    
               endif
            endif
         endif
      enddo
!     
      return        
!
 101  format(3x,13a,3x,a16,a11,3x,a11,3x,a11,3x,a8,3x,a10,3x,a10)
 102  format(3x,i2,8x,3x,a16,a4,3x,e14.7,3x,e14.7,3x,e14.7,3x,a8,3x,a80)
!
      end
