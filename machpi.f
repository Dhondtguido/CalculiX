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
!
!     author: Yannick Muller
!
      subroutine machpi (MACH, PI,kappa, rgas)
!
!-----------------------------------------------------------------------
!                                                                      |
!     Dieses Unterprogramm berechnet die Mach-Zahl fuer das            |
!     eingegebene Druckverhaeltnis PI.                                 |
!                                                                      |
!     Eingabe-Groessen:                                                |
!       PI     = Druckverhaeltnis PS/PT                                |
!                                                                      |
!     Ausgabe-Groessen:                                                |
!       MACH   = Mach-Zahl                                             |
!                                                                      |
!-----------------------------------------------------------------------
!
      IMPLICIT CHARACTER*1 (A-Z)
      real*8    PI, MACH, MA2, kappa, rgas, kappam,kappax,pikrit
!
!-----------------------------------------------------------------------
!
      kappax = (kappa-1)/kappa
      KAPPAM = 2.d0 / (KAPPA - 1.d0)
      PIKRIT = (2.d0/(KAPPA+1.d0)) ** (KAPPA/(KAPPA-1.d0))
!
      IF (PI.GE.1.d0) THEN
!       Druckverhaeltnis groesser gleich 1
        MACH = 0.d0
      ELSEIF (PI.GT.PIKRIT) THEN
!       Druckverhaeltnis unterkritisch
        MA2  = KAPPAM * (PI**(-KAPPAX) - 1.d0)
        IF (MA2.GT.0.d0) THEN
          MACH = SQRT (MA2)
        ELSE
          MACH = 0.d0
        ENDIF
      ELSEIF (PI.GT.0.d0) THEN
!       Druckverhaeltnis ueberkritisch
        MACH = 1.d0
      ELSE
!       Druckverhaeltnis ungueltig
        MACH = 1.d20
      ENDIF
!
      RETURN
      END
!
