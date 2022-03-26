C-----------------------------------------------------------
C     Neutrinoproduction of Pions: Nucleus Functions
C-----------------------------------------------------------
C
C --- Set nucleus parameters and functions
C
C     Z - nucleus charge
C     A - nucleus weight
C
      SUBROUTINE SETNUCL(Z,A)
      INTEGER            Z,A
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.inc'
      INCLUDE 'npmath.fun'
      INTEGER  j0,j1,j2,i,j,k
      REAL*8   f0,f1,f2,T0,T1,T2
      REAL*8   RR,RA,RB
      SAVE     RR,RA,RB
      DATA     RR /0D0/
      DATA     RA /0D0/
      DATA     RB /0D0/
      REAL*8   Rho,RhoRR,RhoBZB
      EXTERNAL Rho,RhoRR,RhoBZB
      NuclZ = Z
      NuclA = A
      RhoB  = 0D0
C
C ... C12
C
      IF     ((Z .EQ.  6) .AND. (A .EQ.  12)) THEN
        Nucl = 'C12'
        RhoR = 2.300D0*fermi
        RhoA = 0.455D0*fermi
C
C ... N14
C
      ELSEIF ((Z .EQ.  7) .AND. (A .EQ.  14)) THEN
        Nucl = 'N14'
        RhoR = 2.5700D0*fermi
        RhoA = 0.5052D0*fermi
        RhoB =-0.1800D0
C
C ... O16
C
      ELSEIF ((Z .EQ.  8) .AND. (A .EQ.  16)) THEN
        Nucl = 'O16'
        RhoR = 2.608D0*fermi
        RhoA = 0.513D0*fermi
        RhoB =-0.051D0
C
C ... F19
C
      ELSEIF ((Z .EQ.  9) .AND. (A .EQ.  19)) THEN
        Nucl = 'F19'
        RhoR = 2.580D0*fermi
        RhoA = 0.567D0*fermi
C
C ... Ne20
C
      ELSEIF ((Z .EQ. 10) .AND. (A .EQ.  20)) THEN
        Nucl = 'Ne20'
        RhoR = 2.805D0*fermi
        RhoA = 0.571D0*fermi
C
C ... Ne22
C
      ELSEIF ((Z .EQ. 10) .AND. (A .EQ.  22)) THEN
        Nucl = 'Ne22'
        RhoR = 2.782D0*fermi
        RhoA = 0.549D0*fermi
C
C ... Mg24
C
      ELSEIF ((Z .EQ. 12) .AND. (A .EQ.  24)) THEN
        Nucl = 'Mg24'
        RhoR = 3.108D0*fermi
        RhoA = 0.607D0*fermi
        RhoB =-0.163D0
C
C ... Al27
C
      ELSEIF ((Z .EQ. 13) .AND. (A .EQ.  27)) THEN
        Nucl = 'Al27'
        RhoR = 2.840D0*fermi
        RhoA = 0.569D0*fermi
C
C ... Si28
C
      ELSEIF ((Z .EQ. 14) .AND. (A .EQ.  28)) THEN
        Nucl = 'Si28'
        RhoR = 3.340D0*fermi
        RhoA = 0.580D0*fermi
        RhoB =-0.233D0
C
C ... P31
C
      ELSEIF ((Z .EQ. 15) .AND. (A .EQ.  31)) THEN
        Nucl = 'P31'
        RhoR = 3.369D0*fermi
        RhoA = 0.582D0*fermi
        RhoB =-0.173D0
C
C ... Cl35
C
      ELSEIF ((Z .EQ. 17) .AND. (A .EQ.  35)) THEN
        Nucl = 'Cl35'
        RhoR = 3.476D0*fermi
        RhoA = 0.599D0*fermi
        RhoB =-0.100D0
C
C ... Cl37
C
      ELSEIF ((Z .EQ. 17) .AND. (A .EQ.  37)) THEN
        Nucl = 'Cl37'
        RhoR = 3.554D0*fermi
        RhoA = 0.588D0*fermi
        RhoB =-0.130D0
C
C ... Ar40
C
      ELSEIF ((Z .EQ. 18) .AND. (A .EQ.  40)) THEN
        Nucl = 'Ar40'
        RhoR = 3.73D0*fermi
        RhoA = 0.62D0*fermi
        RhoB =-0.19D0
C
C ... Ca40
C
      ELSEIF ((Z .EQ. 20) .AND. (A .EQ.  40)) THEN
        Nucl = 'Ca40'
        RhoR = 3.766D0*fermi
        RhoA = 0.586D0*fermi
        RhoB =-0.161D0
C
C ... Ca48
C
      ELSEIF ((Z .EQ. 20) .AND. (A .EQ.  48)) THEN
        Nucl = 'Ca48'
        RhoR = 3.7369D0*fermi
        RhoA = 0.5245D0*fermi
        RhoB =-0.0300D0
C
C ... Ti48
C
      ELSEIF ((Z .EQ. 22) .AND. (A .EQ.  48)) THEN
        Nucl = 'Ti48'
        RhoR = 3.843D0*fermi
        RhoA = 0.588D0*fermi
C
C ... Cr50
C
      ELSEIF ((Z .EQ. 24) .AND. (A .EQ.  50)) THEN
        Nucl = 'Cr50'
        RhoR = 3.979D0*fermi
        RhoA = 0.520D0*fermi
C
C ... Cr52
C
      ELSEIF ((Z .EQ. 24) .AND. (A .EQ.  52)) THEN
        Nucl = 'Cr52'
        RhoR = 4.010D0*fermi
        RhoA = 0.497D0*fermi
C
C ... Cr53
C
      ELSEIF ((Z .EQ. 24) .AND. (A .EQ.  53)) THEN
        Nucl = 'Cr53'
        RhoR = 4.000D0*fermi
        RhoA = 0.557D0*fermi
C
C ... Cr54
C
      ELSEIF ((Z .EQ. 24) .AND. (A .EQ.  54)) THEN
        Nucl = 'Cr54'
        RhoR = 4.021D0*fermi
        RhoA = 0.524D0*fermi
C
C ... Mn55
C
      ELSEIF ((Z .EQ. 25) .AND. (A .EQ.  55)) THEN
        Nucl = 'Mn55'
        RhoR = 3.890D0*fermi
        RhoA = 0.567D0*fermi
C
C ... Ni58
C
      ELSEIF ((Z .EQ. 28) .AND. (A .EQ.  58)) THEN
        Nucl = 'Ni58'
        RhoR = 4.3092D0*fermi
        RhoA = 0.5169D0*fermi
        RhoB =-0.1308D0
C
C ... Ni60
C
      ELSEIF ((Z .EQ. 28) .AND. (A .EQ.  60)) THEN
        Nucl = 'Ni60'
        RhoR = 4.4891D0*fermi
        RhoA = 0.5369D0*fermi
        RhoB =-0.2668D0
C
C ... Cu63
C
      ELSEIF ((Z .EQ. 29) .AND. (A .EQ.  63)) THEN
        Nucl = 'Cu63'
        RhoR = 4.218D0*fermi
        RhoA = 0.595D0*fermi
C
C ... Cu65
C
      ELSEIF ((Z .EQ. 29) .AND. (A .EQ.  65)) THEN
        Nucl = 'Cu65'
        RhoR = 4.218D0*fermi
        RhoA = 0.595D0*fermi
C
C ... Ge70
C
      ELSEIF ((Z .EQ. 32) .AND. (A .EQ.  70)) THEN
        Nucl = 'Ge70'
        RhoR = 4.440D0*fermi
        RhoA = 0.585D0*fermi
C
C ... Ge72
C
      ELSEIF ((Z .EQ. 32) .AND. (A .EQ.  72)) THEN
        Nucl = 'Ge72'
        RhoR = 4.450D0*fermi
        RhoA = 0.573D0*fermi
C
C ... Nb93
C
      ELSEIF ((Z .EQ. 41) .AND. (A .EQ.  93)) THEN
        Nucl = 'Nb93'
        RhoR = 4.870D0*fermi
        RhoA = 0.573D0*fermi
C
C ... Pd110
C
      ELSEIF ((Z .EQ. 46) .AND. (A .EQ. 110)) THEN
        Nucl = 'Pd110'
        RhoR = 5.301D0*fermi
        RhoA = 0.581D0*fermi
C
C ... Sn120
C
      ELSEIF ((Z .EQ. 50) .AND. (A .EQ. 120)) THEN
        Nucl = 'Sn120'
        RhoR = 5.315D0*fermi
        RhoA = 0.576D0*fermi
C
C ... La139
C
      ELSEIF ((Z .EQ. 57) .AND. (A .EQ. 139)) THEN
        Nucl = 'La139'
        RhoR = 5.710D0*fermi
        RhoA = 0.535D0*fermi
C
C ... Gd156
C
      ELSEIF ((Z .EQ. 64) .AND. (A .EQ. 156)) THEN
        Nucl = 'Gd156'
        RhoR = 5.930D0*fermi
        RhoA = 0.576D0*fermi
C
C ... Ta181
C
      ELSEIF ((Z .EQ. 73) .AND. (A .EQ. 181)) THEN
        Nucl = 'Ta181'
        RhoR = 6.38D0*fermi
        RhoA = 0.64D0*fermi
C
C ... Pb207
C
      ELSEIF ((Z .EQ. 82) .AND. (A .EQ. 207)) THEN
        Nucl = 'Pb207'
        RhoR = 6.620D0*fermi
        RhoA = 0.546D0*fermi
C
C ... U238
C
      ELSEIF ((Z .EQ. 92) .AND. (A .EQ. 238)) THEN
        Nucl = 'U238'
        RhoR = 6.8054D0*fermi
        RhoA = 0.6050D0*fermi
      ELSE
C
C ... Generic
C
        WRITE (Nucl,'(2HA(,I3,1H,,I3,1H))') A,Z
        RhoR = 1.120D0*fermi*A**(1D0/3D0)
        RhoA = 0.530D0*fermi
      ENDIF
C
C ... m_N
C
      KinMN = (Z*Massp+(A-Z)*Massn)/A
C
C ... Compare with previous parameters
C
      IF ((RhoR .EQ. RR) .AND.
     .    (RhoA .EQ. RA) .AND.
     .    (RhoB .EQ. RB)) RETURN
      RR = RhoR
      RA = RhoA
      RB = RhoB
C
C ... rho(r)
C
      RhoRmax = RhoR+15*RhoA
      IF (RhoB .LT. 0D0) THEN
        RhoRmax = MIN(RhoRmax,SQRT(-1D0/RhoB)*RhoR)
      ENDIF
      RhoC = 1D0
      RhoC = A/(4D0*Pi*DINTEG(RhoRR,0D0,RhoRmax,-1D-10,50))
C
C ... BpN
C
      DO k = 1,RHOBNUM
        RhoBdat(k) = 16D0/(RHOBNUM-1)*(k-1)
      ENDDO
C
C ... rho(b,z,BpN)
C
      IF ((RHORNUM .LT. 3) .OR. (RHORNUM/2*2 .EQ. RHORNUM)) THEN
        CALL NPwarn('SETNUCL: RHORNUM should be odd and >= 3')
        STOP
      ENDIF
      IF ((RHOBNUM .LT. 3) .OR. (RHOBNUM/2*2 .EQ. RHOBNUM)) THEN
        CALL NPwarn('SETNUCL: RHOBNUM should be odd and >= 3')
        STOP
      ENDIF
      RhoH = RhoRmax/(RHORNUM-1)
      DO i = 1,RHORNUM
        RhoRdat(i) = RhoH*(i-1)
      ENDDO
      DO i = 1,RHORNUM
      DO j = 1,RHORNUM
        RhoMdat(i,j,1) = Rho(SQRT(RhoRdat(i)**2+RhoRdat(j)**2))
      ENDDO
      ENDDO
      DO k = 2,RHOBNUM
      DO i = 1,RHORNUM
      DO j = 1,RHORNUM
        RhoMdat(i,j,k) = RhoBZB(RhoRdat(i),RhoRdat(j),RhoBdat(k))
      ENDDO
      ENDDO
      ENDDO
C
C ... T_z(b,z,BpN)
C
      DO k = 1,RHOBNUM
      DO i = 1,RHORNUM
        j  =   RHORNUM
        T2             = 0D0
        RhoTdat(i,j,k) = 0D0
        DO WHILE (j .GE. 3)
          j0 = j-2
          j1 = j-1
          j2 = j
          f0 = RhoMdat(i,j0,k)
          f1 = RhoMdat(i,j1,k)*0.666666666667D0
          f2 = RhoMdat(i,j2,k)
          T1 = T2+RhoH*(-0.083333333333D0*f0+0.416666666667D0*f2+f1)
          T0 = T1+RhoH*(-0.083333333333D0*f2+0.416666666667D0*f0+f1)
          RhoTdat(i,j1,k) = T1
          RhoTdat(i,j0,k) = T0
          T2 = T0
          j  = j0
        ENDDO
      ENDDO
      ENDDO
      END
C
C --- rho(r)
C
      FUNCTION Rho(r)
      REAL*8   Rho,r
      INCLUDE 'npmain.inc'
      Rho = RhoC/(1D0+EXP((r-RhoR)/RhoA))*(1D0+RhoB*(r/RhoR)**2)
      END
C
C --- rho(r) r**2
C
      FUNCTION RhoRR(r)
      REAL*8   RhoRR,r
      REAL*8   Rho
      EXTERNAL Rho
      RhoRR = r*r*Rho(r)
      END
C
C --- rho(b,z,BpN), integral
C
      FUNCTION RhoBZB(b,z,BpN)
      REAL*8   RhoBZB,b,z,BpN
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      REAL*8          bb,b2,z2,BN,BN2
      COMMON  /RhoBZ/ bb,b2,z2,BN,BN2
      REAL*8   RhoBZBb
      EXTERNAL RhoBZBb
      bb     = b
      b2     = b*b
      z2     = z*z
      BN     = 1.0D0/BpN
      BN2    = 0.5D0*BN
      RhoBZB = BN*DINTEG(RhoBZBb,0D0,RhoRmax,RhoC*1D-7,10)
      END

      FUNCTION RhoBZBb(bt)
      REAL*8   RhoBZBb,bt
      INCLUDE 'npmath.fun'
      REAL*8          bb,b2,z2,BN,BN2
      COMMON  /RhoBZ/ bb,b2,z2,BN,BN2
      REAL*8   Rho,bt2,bte
      EXTERNAL Rho
      bt2     = bt*bt
      bte     = DLOGI0(bb*bt*BN)-(b2+bt2)*BN2
      RhoBZBb = bt*Rho(SQRT(z2+bt2))*EXP(bte)
      END
C
C --- rho(b,z,BpN), approximation
C
      FUNCTION RhoM(b,z,BpN)
      REAL*8   RhoM,b,z,BpN
      INCLUDE 'npmain.inc'
      REAL*8   DINTER3
      EXTERNAL DINTER3
      RhoM = DINTER3(
     .  RHORNUM,RHORNUM,RHOBNUM,
     .  RhoRdat,RhoRdat,RhoBdat,RhoMdat,b,ABS(z),BpN)
      END
C
C --- T_z(b,z,BpN), approximation
C
      FUNCTION RhoT(b,z,BpN)
      REAL*8   RhoT,b,z,BpN
      INCLUDE 'npmain.inc'
      REAL*8   DINTER3
      EXTERNAL DINTER3
      IF (z .GE. 0D0) THEN
        RhoT = DINTER3(
     .    RHORNUM,RHORNUM,RHOBNUM,
     .    RhoRdat,RhoRdat,RhoBdat,RhoTdat,b,  z,BpN)
      ELSE
        RhoT = DINTER3(
     .    RHORNUM,RHORNUM,RHOBNUM,
     .    RhoRdat,RhoRdat,RhoBdat,RhoTdat,b,0D0,BpN)*2D0
     .       - DINTER3(
     .    RHORNUM,RHORNUM,RHOBNUM,
     .    RhoRdat,RhoRdat,RhoBdat,RhoTdat,b, -z,BpN)
      ENDIF
      END
C
C --- Set BpN for RhoMB(b,z) and RhoTB(b,z)
C
      SUBROUTINE SETBPN(BpN)
      REAL*8            BpN
      INCLUDE 'npmain.inc'
      INTEGER  k1,k2,k3,i,j
      REAL*8   v1,v2,v3,v12,v13,v23,vv1,vv2,vv3
      INTEGER  DINTERind
      EXTERNAL DINTERind
      k1  = DINTERind(RHOBNUM,RhoBdat,BpN)
      k2  = k1+1
      k3  = k2+1
      v1  = BpN-RhoBdat(k1)
      v2  = BpN-RhoBdat(k2)
      v3  = BpN-RhoBdat(k3)
      v12 =  v2-v1
      v13 =  v3-v1
      v23 =  v3-v2
      vv1 =  v2*v3/(v12*v13)
      vv2 = -v1*v3/(v12*v23)
      vv3 =  v1*v2/(v13*v23)
      DO i = 1,RHORNUM
      DO j = 1,RHORNUM
        RhoMdatB(i,j) = vv1*RhoMdat(i,j,k1)
     .                + vv2*RhoMdat(i,j,k2)
     .                + vv3*RhoMdat(i,j,k3)
        RhoTdatB(i,j) = vv1*RhoTdat(i,j,k1)
     .                + vv2*RhoTdat(i,j,k2)
     .                + vv3*RhoTdat(i,j,k3)
      ENDDO
      ENDDO
      END
C
C --- rho(b,z), approximation
C
      FUNCTION RhoMB(b,z)
      REAL*8   RhoMB,b,z
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      RhoMB = DINTER2(
     .  RHORNUM,RHORNUM,
     .  RhoRdat,RhoRdat,RhoMdatB,b,ABS(z))
      END
C
C --- T_z(b,z), approximation
C
      FUNCTION RhoTB(b,z)
      REAL*8   RhoTB,b,z
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      IF (z .GE. 0D0) THEN
        RhoTB = DINTER2(
     .    RHORNUM,RHORNUM,
     .    RhoRdat,RhoRdat,RhoTdatB,b,  z)
      ELSE
        RhoTB = DINTER2(
     .    RHORNUM,RHORNUM,
     .    RhoRdat,RhoRdat,RhoTdatB,b,0D0)*2D0
     .        - DINTER2(
     .    RHORNUM,RHORNUM,
     .    RhoRdat,RhoRdat,RhoTdatB,b, -z)
      ENDIF
      END
