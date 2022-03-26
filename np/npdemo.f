C*********************************************************************
C*                                                                   *
C*  This is demo program for calculation of differential and total   *
C*  cross sections  for  neutrinoproduction of pions off nuclei in   *
C*  the approach of B.Z. Kopeliovich et al.                          *
C*                                                                   *
C*********************************************************************
C
      PROGRAM    NPDEMO
C
C     Kinematic variables
C
      REAL*8     E                 !! Neutrino energy
      REAL*8     x,xmin,xmax       !! x, x min., x max.
      REAL*8     y,ymin,ymax       !! y, y min., y max.
      REAL*8     t,tmin,tmax       !! t, t min., t max.
C
C     Functions for total and differential cross sections
C
      REAL*8     SigNA,SigNAx,SigNAxy,SigNAxyt
      EXTERNAL   SigNA,SigNAx,SigNAxy,SigNAxyt
C
C     First step is to set a reaction. One should specify:
C
C       Lepton = -1 (  e-), 1 (  e+)
C                -2 ( mu-), 2 ( mu+)
C                -3 (tau-), 3 (tau+)
C                 0 (nu)
C
C       Pion   = -1 (pi-)
C                 1 (pi+)
C                 0 (pi0)
C
C       Nucleus Z and A
C
C     For example, for charge current neutrino scattering off neon 
C
C       nu Ne --> mu- pi+ Ne
C
C     one should call (for neon Z=10 and A=20)
C
C                  mu- pi+  Z   A
C                   |   |   |   |
C                   v   v   v   v
      CALL SETREAC(-2,  1, 10, 20)
C
C     All variables, if not specified exactly, are in GeV units. For
C     cross sections one can specify unit value. Default value is 1.
C     It corresponds to cross secton in GeV^{-2} units. Here we want
C     to have results for cross sections in 10^{-40} cm^2. As far as
C     1 fm = 5.067728853 GeV^{-1} and 1 cm = 10^{13} fm one can set
C
      CALL SETSIGU(1D-40*(1D13*5.067728853D0)**2)
C
C     Next step is to find kinematical regions for all variables. In
C     this program we use variables "x", "y" and "t". Instead of "x"
C     and "y" one can use "nu" and "Q^2".  For example, set neutrino
C     energy to some value E = 2 GeV
C
      E = 2
C
C     To find allowed interval for variable "x" (i.e. xmin and xmax)
C     one can call
C
      CALL KINEX(E,xmin,xmax)
C
C     Now we select some "x" in this interval
C
      x = xmin+(xmax-xmin)*0.2D0
C
C     And find interval for variable "y"
C
      CALL KINEXY(E,x,ymin,ymax)
C
C     Select some "y"
C
      y = ymin+(ymax-ymin)*0.4D0
C
C     And find interval for "t"
C
      CALL KINEXYT(E,x,y,tmin,tmax)
C
C     Finally set some "t"
C
      t = tmin+(tmax-tmin)*0.8D0
C
C     Print all kinematic variables
C
      PRINT *
      PRINT *,' --- Kinematic variables ---'
      PRINT *
      PRINT 1,'E',E,'GeV  '
      PRINT 1,'x',x
      PRINT 1,'y',y
      PRINT 1,'t',t,'GeV^2'
C
C     Now select cross section type:
C
C       Type = 0 - coherent (default)
C              1 - incoherent
C
      CALL SETSIGT(0)
C
C     And finally differential and total cross sections
C
      PRINT *
      PRINT *,' --- Coherent cross sections ---'
      PRINT *
      PRINT 2,'d sigma / dx dy dt',SigNAxyt(E,x,y,t),'/GeV^2'
      PRINT 2,'d sigma / dx dy   ',SigNAxy (E,x,y)
      PRINT 2,'d sigma / dx      ',SigNAx  (E,x)
      PRINT 2,'  sigma total     ',SigNA   (E)

    1 FORMAT(1X,A3 ,' = ',F8.3,' '             ,A5)
    2 FORMAT(1X,A20,' = ',F8.3,' 10^{-40} cm^2',A6)

      END
