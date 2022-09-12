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
      subroutine flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,rho,physcon,g,co,dvi,numf,
     &     vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,ider,
     &     ttime,time,iaxial,iplausi)
!
!     gas element routines 
!     
!     mass flow input for all gas element routines is the gas
!     flow for a 2 degrees segment with the correct sign
!     (positive if from node 1 to node2 of the element,
!      negative if from node 2 to node1 of the element)
!
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!      
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(8),idirf(8),kflag,ipkon(*),kon(*),
     &     nshcon(*), nrhcon(*),ntmat_,mi(*),ider,iaxial,iplausi
!      
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(8),R,cp,physcon(*),rho,
     &     g(3),co(3,*),dvi,vold(0:mi(2),*),shcon(0:3,ntmat_,*),
     &     rhcon(0:1,ntmat_,*),ttime,time
!
      if((lakon(nelem)(2:4).eq.'ATR')
     &        .or.(lakon(nelem)(2:4).eq.'RTA')) then
!
!        absolute to relative system or vice versa
!     
         call absolute_relative(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,numf,set,mi,ttime,time,iaxial,
     &        iplausi)
!     
      elseif(lakon(nelem)(2:8).eq.'ACCTUBO') then
!
!        code not available in the internet distribution of CalculiX
!         
         call acctube_one(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,mi,ider,
     &        ttime,time,iaxial,iplausi)
!
      elseif(lakon(nelem)(2:8).eq.'ACCTUBE') then
!
!        code not available in the internet distribution of CalculiX
!
         call acctube(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,mi,ider,
     &        ttime,time,iaxial,iplausi)
!
!	 code not available in the internet distribution of CalculiX
!
      elseif(lakon(nelem)(2:5).eq.'AVLV') then 
!
         call air_valve(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,
     &     ttime,time,iaxial,iplausi)
!
      elseif(lakon(nelem)(2:6).eq.'CARBS') then  
!
!        carbon seal
!
         call carbon_seal(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,R,physcon,dvi,numf,set,mi,ttime,time,
     &     iaxial,iplausi)
!     
      elseif(lakon(nelem)(2:5).eq.'CHAR') then 
!
!        characteristic
!     
         call characteristic(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,
     &        mi,ttime,time,iaxial,iplausi)
!
      elseif(lakon(nelem)(2:5).eq.'CROS') then 
!         
!        cross split
!
         call cross_split(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,ttime,time,
     &     iaxial,iplausi)
!
!     code not available in the internet distribution of CalculiX
!
      elseif(lakon(nelem)(2:5).eq.'FDPF') then 
         call free_disc_pumping(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!  
!     code not available in the internet distribution of CalculiX
! 
      elseif(lakon(nelem)(2:5).eq.'FCVF') then 
         call free_convection(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!
!     gas pipe fanno
!
      elseif(lakon(nelem)(2:5).eq.'GAPF') then 
!
         call gaspipe_fanno(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!
!     rotating gas pipe
!
      elseif(lakon(nelem)(2:5).eq.'GAPR') then 
!
         call gaspipe_rot(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!
!     straight and stepped labyrinth
!
      elseif(lakon(nelem)(2:4).eq.'LAB') then 
!         
         call labyrinth(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,co,dvi,numf,vold,set,
     &     kon,ipkon,mi,ttime,time,iaxial,iplausi)
!
!     liquid pipes including loss elements (hydraulic elements)
!         
      elseif(lakon(nelem)(2:5).eq.'LIPI') then
!         
         call liquidpipe(node1,node2,nodem,nelem,lakon,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,dvi,numf,vold,mi,ipkon,kon,set,ttime,time,
     &           iaxial,iplausi)
!
!     liquid pump
!
      elseif(lakon(nelem)(2:5).eq.'LIPU') then
!        
         call liquidpump(node1,node2,nodem,nelem,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,numf,mi,ttime,time,iaxial,iplausi)
!
!     liquid pipes including loss elements (types derived from their
!     compressible equivalent)
!
      elseif(lakon(nelem)(2:3).eq.'LP') then
!         
         call liquidpipe(node1,node2,nodem,nelem,lakon,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,dvi,numf,vold,mi,ipkon,kon,set,ttime,time,
     &           iaxial,iplausi)
!
!     element that fixes the mass flow as a specific percentage of the
!     sum of the massflow of up to 10 other elements
!
      elseif(lakon(nelem)(2:5).eq.'MFPC') then 
         call massflow_percent(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!     
!     Moehring
! 
      elseif(lakon(nelem)(2:4).eq.'MRG') then 
!     
         call moehring(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,dvi,numf,set,mi,ttime,time,
     &        iaxial,iplausi)
!
!     Bleed tapping, orifice and pre-swirl nozzle
!
      elseif(lakon(nelem)(2:3).eq.'OR') then 
!         
         call orifice(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,co,vold,mi,
     &        ttime,time,iaxial,iplausi)
! 
!     code not available in the internet distribution of CalculiX
!
      elseif(lakon(nelem)(2:4).eq.'RCV') then   
!
         call rcavi(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,dvi,numf,set,mi,ttime,time,
     &     iaxial,iplausi)
!
!     code not available in the internet distribution of CalculiX
!
      elseif(lakon(nelem)(2:3).eq.'RO') then   
!
         call rcavi2(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,dvi,numf,set,mi,ttime,time,
     &     iaxial,iplausi)
!
!     restrictors
!
      elseif(((lakon(nelem)(2:3).eq.'RE').or.
     &        (lakon(nelem)(2:3).eq.'RB')).and.
     &       (lakon(nelem)(2:8).ne.'REBRSI1').and.
     &       (lakon(nelem)(2:8).ne.'REBRSI2')) then 
!          
         call restrictor(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,mi,ttime,time,iaxial,
     &        co,vold,iplausi)
!
!     code not available in the internet distribution of CalculiX
!
      elseif((lakon(nelem)(2:5).eq.'RIMS').or.
     &       (lakon(nelem)(2:8).eq.'RIMFLEX')) then   
!
         call rimseal(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,mi,
     &        ttime,time,iaxial,co,vold,iplausi)
!     
!     code not available in the internet distribution of CalculiX
!
      elseif(lakon(nelem)(2:6).eq.'SPUMP') then 
!
        call scavenge_pump(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,ntmat_,mi,
     &        ttime,time,iaxial,iplausi)   
!   
!     branch split Idelchik2
!  
      elseif(lakon(nelem)(2:8).eq.'REBRSI2') then 
!         
         call tee(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,ttime,time,
     &     iaxial,iplausi)
!
!     user element
!
      elseif(lakon(nelem)(2:2).eq.'U') then
!         
         call user_network_element(node1,node2,nodem,nelem,lakon,kon,
     &        ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,co,vold,mi,
     &        ttime,time,iaxial,iplausi)
!
!     vortex
!
      elseif(lakon(nelem)(2:3).eq.'VO') then 
!     
         call vortex(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,numf,set,mi,ttime,time,iaxial,
     &        iplausi) 
!
!     branch split Idelchik1
!  
      elseif(lakon(nelem)(2:8).eq.'REBRSI1') then 
!         
         call wye(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,ttime,time,
     &     iaxial,iplausi,dvi)
!
      else
         identity=.true.
!          
      endif
!    
      return
      end
      
