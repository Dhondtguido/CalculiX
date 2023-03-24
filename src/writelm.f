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
      subroutine writelm(iter,xlambd,nactive,nnlconst,objectset,
     &   nobject,ipoacti,iconstacti,inameacti,nodedesi,dgdxglob,nk)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer nactive,nobject,nnlconst,iter,ipos,i,ipoacti(*),
     &   iconstacti(*),inameacti(*),nodedesi(*),node,nk
!
      real*8 xlambd(*),dgdxglob(2,nk,*),val
!
      write(5,*)
      write(5,*)
      write(5,*) '  #######################################
     &#####################################'
      if(iter.eq.1) then      
         write(5,*) '  L A G R A N G E   M U L T I P L I E R S
     &   1ST   I T E R A T I O N' 
      elseif(iter.eq.2) then
         write(5,*) '  L A G R A N G E   M U L T I P L I E R S
     &   2ND   I T E R A T I O N' 
      elseif(iter.eq.3) then
         write(5,*) '  L A G R A N G E   M U L T I P L I E R S
     &   3RD   I T E R A T I O N' 
      elseif((iter.gt.3).and.(iter.lt.10)) then
         write(5,'(a42,i1,a22)') '  L A G R A N G E
     &   M U L T I P L I E R S   ',iter,'TH   I T E R A T I O N' 
      else
         write(5,'(a42,i3,a22)') '  L A G R A N G E
     &   M U L T I P L I E R S   ',iter,'TH   I T E R A T I O N' 
      endif
      write(5,*)
      write(5,103) 'NUMBER OF
     &    ','CONSTRAINT      ','LE/     ','LAGRANGE      ','  ACTIVE/
     & ','   NAME OF' 
      write(5,103) 'CONSTRAINT
     &   ','FUNCTION        ','GE      ','MULTIPLIER    ','  INACTIVE'
     &,'   CONSTRAINT' 
      write(5,*) '  #######################################
     &#####################################'
      write(5,*)
!
      do i=1,nactive
         ipos=ipoacti(i)
!
!        writing of all nonlinear constraints
!
         if(i.le.nnlconst) then
            if(iconstacti(i).eq.-1) then
               if(xlambd(i).lt.0.d0) then            
                  write(5,101)
     &            ipos-1,objectset(1,ipos),'LE  ',xlambd(i),'ACTIVE  ',
     &            objectset(5,ipos)
               else
                  write(5,101)
     &            ipos-1,objectset(1,ipos),'LE  ',xlambd(i),'INACTIVE',
     &            objectset(5,ipos)
               endif
            else
               if(xlambd(i).lt.0.d0) then
                  write(5,101)
     &            ipos-1,objectset(1,ipos),'GE  ',xlambd(i),'INACTIVE', 
     &            objectset(5,ipos)
               else
                  write(5,101)
     &            ipos-1,objectset(1,ipos),'GE  ',xlambd(i),'ACTIVE  ',
     &            objectset(5,ipos)
               endif
            endif   
!
!        writing of all linear (geometric) constraints
!
         else
!
!           MAXMEMBERSIZE
!
            if(objectset(1,inameacti(i))(1:13).eq.'MAXMEMBERSIZE') then
               node=nodedesi(ipoacti(i))
               val=dgdxglob(2,node,inameacti(i))     
               if((xlambd(i).lt.0.d0).and.(val.gt.0.d0)) then
                  write(5,102)
     &            inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &            xlambd(i),'ACTIVE  ',nodedesi(ipos)              
               else
                  write(5,102)
     &            inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &            xlambd(i),'INACTIVE',nodedesi(ipos)          
               endif
!
!           MINMEMBERSIZE
!
            else if(objectset(1,inameacti(i))(1:13).eq.
     &'MINMEMBERSIZE') then
               node=nodedesi(ipoacti(i))
               val=dgdxglob(2,node,inameacti(i))     
               if((xlambd(i).gt.0.d0).and.(val.gt.0.d0)) then
                  write(5,102)
     &            inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &            xlambd(i),'ACTIVE  ',nodedesi(ipos)              
               else
                  write(5,102)
     &            inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &            xlambd(i),'INACTIVE',nodedesi(ipos)          
               endif
!
!           MAXSHRINKAGE
!
            elseif(objectset(1,inameacti(i))(4:12).eq.'SHRINKAGE') then
              node=nodedesi(ipoacti(i))
              val=dgdxglob(2,node,inameacti(i))    
              if((xlambd(i).gt.0.d0).and.(val.ge.0.d0)) then  
                 write(5,102)
     &           inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &           xlambd(i),'ACTIVE  ',nodedesi(ipos)
              else
                 write(5,102)
     &           inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &           xlambd(i),'INACTIVE',nodedesi(ipos) 
              endif                  
!
!           MAXGROWTH
!
            elseif(objectset(1,inameacti(i))(4:9).eq.'GROWTH') then
              node=nodedesi(ipoacti(i))
              val=dgdxglob(2,node,inameacti(i))    
              if((xlambd(i).lt.0.d0).and.(val.ge.0.d0)) then  
                 write(5,102)
     &           inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &           xlambd(i),'ACTIVE  ',nodedesi(ipos)
              else
                 write(5,102)
     &           inameacti(i)-1,objectset(1,inameacti(i)),'LE  ',
     &           xlambd(i),'INACTIVE',nodedesi(ipos) 
              endif                  
!
!           PACKAGING
!
            elseif(objectset(1,inameacti(i))(1:9).eq.'PACKAGING') then
              node=nodedesi(ipoacti(i))
              val=dgdxglob(2,node,inameacti(i))    
              if((xlambd(i).lt.0.d0).and.(val.ge.0.d0)) then  
                 write(5,102)
     &           inameacti(i)-1,objectset(1,inameacti(i)),'GE  ',
     &           xlambd(i),'ACTIVE  ',nodedesi(ipos)
              else
                 write(5,102)
     &           inameacti(i)-1,objectset(1,inameacti(i)),'GE  ',
     &           xlambd(i),'INACTIVE',nodedesi(ipos) 
              endif                  
            endif
         endif
      enddo  
      write(5,*)
!
      return        
!
 101  format(1(3x,i2,8x,3x,a16,a4,3x,e14.7,3x,a8,3x,a80))
 102  format(1(3x,i2,8x,3x,a16,a4,3x,e14.7,3x,a8,3x,i6))
 103  format(1(3x,13a,3x,a16,a8,3x,a14,5x,a10,3x,a10))
!
      end
