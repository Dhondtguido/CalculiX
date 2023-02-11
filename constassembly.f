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
      subroutine constassembly(nobject,objectset,g0,ndesi,dgdxglob,
     &   nk,nodedesi,gradproj,set,nset,nodedesiboun,istartset,
     &   iendset,ialset,nodedesiinv)               
!     
!     assesmbly of all gradients to the assembled gradient vector
!     
      implicit none
!     
      character*81 objectset(5,*),set(*)
      character*20 empty
!     
      integer nobject,istat,i,j,k,ndesi,nk,nodedesi(*),inode,iset,
     &   nset,ndesiboun,nodedesiboun(*),id,istartset(*),iendset(*),
     &   ialset(*),nodedesiinv(*)
!     
      real*8 g0(nobject),scale,bound,dgdxglob(2,nk,*),gradproj(3,*),
     &   dd,funcvalabs,funcvalnorm,objnorm,obj
!     
      empty='                    '
!
      write(5,*)
      write(5,*)
      write(5,'(a113)') '  ################################################
     &#################################################################'
      write(5,*) '  A S S E M B L Y   O F   B A R R I E R   F U N C T I 
     &O N'
      write(5,*)
      write(5,101)
     &'NUMBER OF    ','CONSTRAINT      ','LE/     ','FUNCTION         ',
     &'FUNCTION         ','FUNCTION      ','  ACTIVE/ ','   NAME OF' 
      write(5,101)
     &'CONSTRAINT   ','FUNCTION        ','GE      ','VALUE            ',
     &'BOUND            ','VALUE NORM.   ','  INACTIVE','   CONSTRAINT' 
      write(5,'(a113)') '  ################################################
     &#################################################################'
      write(5,*)
!   
!     Initialize field of assembled constraint gradient
!
      do i=1,ndesi
         inode=nodedesi(i)
         gradproj(1,inode)=0.d0
      enddo
!
!     starting loop from the first constraint 
!     which is the 2nd entry in objectset
!
      do i=2,nobject
!
!        in case of a geometric constraint, determine the design
!        variables which have this constraint
!
         if(objectset(5,i)(81:81).eq.'G') then
            call cident81(set,objectset(3,i),nset,id)
            iset=nset+1
            if(id.gt.0) then
               if(objectset(3,i).eq.set(id)) then
                  iset=id
               endif
            endif
!
            if(iset.le.nset) then
               ndesiboun=0
               do j=istartset(iset),iendset(iset)
                  if(ialset(j).gt.0) then
                     k=ialset(j)
                     if(nodedesiinv(k).ne.1) cycle
                     ndesiboun=ndesiboun+1
                     nodedesiboun(ndesiboun)=k
                  else
                     k=ialset(j-2)
                     do
                        k=k-ialset(j)
                        if(k.ge.ialset(j-1)) exit
                        if(nodedesiinv(k).ne.1) cycle
                        ndesiboun=ndesiboun+1
                        nodedesiboun(ndesiboun)=k
                     enddo
                  endif
               enddo
            endif 
         endif       
!
!        compute entries of the assembled gradient
!
         if(objectset(1,i)(1:13).eq.'MAXMEMBERSIZE') then      
            read(objectset(1,i)(61:80),'(f20.0)',iostat=istat) bound
            do j=1,ndesiboun
               inode=nodedesiboun(j)
               obj=dgdxglob(1,inode,i)
               objnorm=dgdxglob(2,inode,i)         
               if(obj.ge.0.d0) then
                  gradproj(1,inode)=gradproj(1,inode)+1/(-1*objnorm)
                  if(j.eq.1) then
                    funcvalnorm=objnorm
                    funcvalabs=obj
                  else
                     funcvalnorm=max(funcvalnorm,objnorm )
                     funcvalabs=max(funcvalabs,obj)
                  endif
               endif          
            enddo
            write(5,102) i-1,objectset(1,i),objectset(1,i)(19:22),
     &         funcvalabs,bound,funcvalnorm,'ACTIVE  '      
         elseif(objectset(1,i)(1:13).eq.'MINMEMBERSIZE') then     
            read(objectset(1,i)(61:80),'(f20.0)',iostat=istat) bound
            do j=1,ndesiboun
               inode=nodedesiboun(j)
               obj=dgdxglob(1,inode,i)
               objnorm=dgdxglob(2,inode,i)
               if(obj.ge.0.d0) then
                  gradproj(1,inode)=gradproj(1,inode)-1/(-1*objnorm)
                  if(j.eq.1) then
                     funcvalnorm=objnorm
                     funcvalabs=obj
                  else
                     funcvalnorm=max(funcvalnorm,objnorm)
                     funcvalabs=min(funcvalabs,obj)
                  endif
               endif          
            enddo
            write(5,102) i-1,objectset(1,i),objectset(1,i)(19:22),
     &         funcvalabs,bound,funcvalnorm,'ACTIVE  '      
         elseif(objectset(1,i)(1:12).eq.'MAXSHRINKAGE') then     
            read(objectset(1,i)(61:80),'(f20.0)',iostat=istat) bound
            do j=1,ndesiboun
               inode=nodedesiboun(j)
               obj=dgdxglob(1,inode,i)
               objnorm=dgdxglob(2,inode,i)
               gradproj(1,inode)=gradproj(1,inode)-1/(-1*objnorm)
               if(j.eq.1) then
                  funcvalnorm=objnorm
                  funcvalabs=obj
               else
                  funcvalnorm=max(funcvalnorm,objnorm)
                  funcvalabs=min(funcvalabs,obj)
               endif          
            enddo
            write(5,102) i-1,objectset(1,i),objectset(1,i)(19:22),
     &         funcvalabs,bound,funcvalnorm,'ACTIVE  '      
         elseif(objectset(1,i)(1:9).eq.'MAXGROWTH') then     
            read(objectset(1,i)(61:80),'(f20.0)',iostat=istat) bound
            do j=1,ndesiboun
               inode=nodedesiboun(j)
               obj=dgdxglob(1,inode,i)
               objnorm=dgdxglob(2,inode,i)
               gradproj(1,inode)=gradproj(1,inode)+1/(-1*objnorm)
               if(j.eq.1) then
                  funcvalnorm=objnorm
                  funcvalabs=obj
               else
                  funcvalnorm=max(funcvalnorm,objnorm)
                  funcvalabs=max(funcvalabs,obj)
               endif          
            enddo
            write(5,102) i-1,objectset(1,i),objectset(1,i)(19:22),
     &         funcvalabs,bound,funcvalnorm,'ACTIVE  '      
         elseif(objectset(1,i)(1:9).eq.'PACKAGING') then     
            bound=0.d0
            do j=1,ndesiboun
               inode=nodedesiboun(j)
               obj=dgdxglob(1,inode,i)
               objnorm=dgdxglob(2,inode,i)
               gradproj(1,inode)=gradproj(1,inode)+1/(-1*objnorm)
               if(j.eq.1) then
                  funcvalnorm=objnorm
                  funcvalabs=obj
               else
                  funcvalnorm=max(funcvalnorm,objnorm)
                  funcvalabs=min(funcvalabs,obj)
               endif        
            enddo
            write(5,102) i-1,objectset(1,i),objectset(1,i)(19:22),
     &         funcvalabs,bound,funcvalnorm,'ACTIVE  '      
!
!        all nonlinear constraints
!
         else
            if(objectset(1,i)(61:80).ne.empty) then
               read(objectset(1,i)(61:80),'(f20.0)',
     &            iostat=istat) bound
            else
               write(*,*) '*WARNING in checkconstraint'
               write(*,*) '         no absolute constraint boundary'
               write(*,*) '         defined, system value taken' 
               write(*,*)
               bound=g0(i)
            endif
            if(objectset(1,i)(41:60).ne.empty) then
               read(objectset(1,i)(41:60),'(f20.0)',
     &           iostat=istat) scale
            else
               write(*,*) '*WARNING in checkconstraint'
               write(*,*) '         no relative constraint boundary'
               write(*,*) '         defined, 1.0 taken' 
               write(*,*)
               scale=1.0d0
            endif
            bound=bound*scale
!
            if(objectset(1,i)(19:20).eq.'LE') then
               objnorm=g0(i)/bound-1
               do j=1,ndesi
                  inode=nodedesi(j)
                  gradproj(1,inode)=gradproj(1,inode)+
     &               dgdxglob(2,inode,i)/(-1*objnorm)
               enddo
            elseif(objectset(1,i)(19:20).eq.'GE') then
               objnorm=g0(i)/bound+1
               do j=1,ndesi
                  inode=nodedesi(j)
                  gradproj(1,inode)=gradproj(1,inode)+
     &               dgdxglob(2,inode,i)/objnorm       
               enddo
            endif
            write(5,102) i-1,objectset(1,i),objectset(1,i)(19:22),
     &         g0(i),bound,objnorm,'ACTIVE  ',objectset(5,i)    
         endif
      enddo
!   
!     normalization of assembled constraint gradient vector   
!
      if(nobject.gt.1) then   
         dd=0.d0
         do i=1,ndesi
            inode=nodedesi(i)
            dd=dd+gradproj(1,inode)**2
         enddo
         if(dd.le.0.d0) then
            dd=1.0
         else
            dd=dsqrt(dd) 
         endif
         do i=1,ndesi
            inode=nodedesi(i)
            gradproj(1,inode)=gradproj(1,inode)/dd
         enddo
      endif
!      
 101  format(3x,13a,3x,a16,a11,3x,a11,3x,a11,3x,a8,3x,a10,3x,a10)
 102  format(3x,i2,8x,3x,a16,a4,3x,e14.7,3x,e14.7,3x,e14.7,3x,a8,3x,a80)
!
      return        
      end
