C-----------------------------------------------------------
C     Neutrinoproduction of Pions: Cross Sections
C-----------------------------------------------------------
C
C --- Default values
C
      BLOCK DATA NPMAIN
      INCLUDE   'npmain.inc'
      DATA OutFile   /6/
      DATA OutLevel  /1/
      DATA SigType   /0/
      DATA SigUnit   /1D0/
      DATA SigEps    /1D-2/
      END
C
C --- Output File and printout Level
C
      SUBROUTINE SETOUT(F,L)
      INTEGER :: F
      INTEGER :: L
      INCLUDE   'npmain.inc'
      OutFile = F
      IF (F .GT. 0) THEN
        OutLevel = L
      ELSE
        OutLevel = 0
      ENDIF
      END
C
C --- Warning messages
C
      SUBROUTINE NPwarn(Text)
      CHARACTER*(*)     Text
      INCLUDE   'npmain.inc'
      IF (OutLevel .GE. 1) WRITE (OutFile,*) '*** ',Text
      END

      SUBROUTINE NPwarnI (Text,V)
      CHARACTER*(*)       Text
      INTEGER                  V
      CHARACTER*60 Msg
      WRITE       (Msg,*) Text,V
      CALL NPwarn (Msg)
      END

      SUBROUTINE NPwarnR (Text,V)
      CHARACTER*(*)       Text
      REAL                     V
      CHARACTER*60 Msg
      WRITE       (Msg,*) Text,V
      CALL NPwarn (Msg)
      END

      SUBROUTINE NPwarnD (Text,V)
      CHARACTER*(*)       Text
      REAL*8                   V
      CHARACTER*60 Msg
      WRITE       (Msg,*) Text,V
      CALL NPwarn (Msg)
      END
C
C --- Output parameters
C
      SUBROUTINE OUTPAR(F,P)
      INTEGER           F
      CHARACTER*(*)       P
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.inc'
      IF (OutLevel .GE. 1) THEN
        WRITE (F,1) P
        WRITE (F,*) P
        WRITE (F,*) P,' ',Reac
        WRITE (F,*) P
        IF (OutLevel .GE. 2) THEN
          WRITE (F,*) P,' KinML    = ',KinML
          WRITE (F,*) P,' KinMP    = ',KinMP
          WRITE (F,*) P,' KinMN    = ',KinMN
          WRITE (F,*) P,' RhoA     = ',RhoA   /fermi,' [fm]'
          WRITE (F,*) P,' RhoR     = ',RhoR   /fermi,' [fm]'
          WRITE (F,*) P,' RhoRmax  = ',RhoRmax/fermi,' [fm]'
          WRITE (F,*) P,' RhoH     = ',RhoH   /fermi,' [fm]'
          WRITE (F,*) P,' SigUnit  = ',SigUnit
          WRITE (F,*) P,' SigFact  = ',SigFact
          WRITE (F,*) P,' SigEps   = ',SigEps
          WRITE (F,*) P
        ENDIF
        WRITE (F,1) P
      ENDIF
    1 FORMAT(1X,1A,50(1H-))
      END
C
C --- Set reaction
C
C     L - Lepton  charge (0 - nu , +1/-1 - e, +2/-2 - mu, +3/-3 -tau)
C     P - Pion    charge
C     Z - nucleus charge
C     A - nucleus weight
C
      SUBROUTINE SETREAC(L,P,Z,A) BIND(C)
      USE ISO_C_BINDING
      INTEGER (C_INT), VALUE :: L,P,Z,A
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.inc'
      INTEGER      C
      CHARACTER*40 proc
      CHARACTER*5  neut
      CHARACTER*4  lept(-3:3)
      CHARACTER*3  pion(-1:1)
      DATA         lept /'tau-','mu-','e-','nu','e+','mu+','tau+'/
      DATA         pion /'pi-','pi0','pi+'/
      INTEGER      i,j,m
C
C ... Check process
C
      IF (L .EQ. 0) THEN
        C = 0
      ELSE
        C = L/ABS(L)
      ENDIF
      IF (C+P .NE. 0) THEN
        CALL NPwarn('SETREAC: invalid process')
        STOP
      ENDIF
C
C ... Lepton
C
      IF (L .EQ. 0) THEN
        KinML = 0D0
      ELSEIF (ABS(L) .EQ. 1) THEN
        KinML = Masse
      ELSEIF (ABS(L) .EQ. 2) THEN
        KinML = Massmu
      ELSEIF (ABS(L) .EQ. 3) THEN
        KinML = Masstau
      ELSE
        CALL NPwarn('SETREAC: invalid lepton')
        STOP
      ENDIF
C
C ... Pion
C
      IF (ABS(P) .GT. 1) THEN
        CALL NPwarn('SETREAC: invalid pion')
        STOP
      ENDIF
C
C ... Charge factor
C
      IF (L .EQ. 0) THEN
        SigFact = (ConstG*Constfpi/Pi)**2/2D0
      ELSE
        SigFact = (ConstG*Constfpi/Pi)**2
      ENDIF
      CALL SETNUCL(Z,A)
      CALL SETPION(P)
C
C ... Neutrino threshold
C
      KinEM = (KinML+KinMP)*(1D0+(KinML+KinMP)/(2D0*KinMN))
C
C ... Reaction name
C
      IF (L .GT. 0) THEN
        neut = 'nubar'
      ELSE
        neut = 'nu'
      ENDIF
      proc = neut   //' '              //Nucl//' --> '
     .     //lept(L)//' '//pion(P)//' '//Nucl
      j = 1
      m = 1
      DO i = 1,40
        Reac(i:i) = ' '
        IF (proc(i:i) .NE. ' ') THEN
            Reac(j:j) = proc(i:i)
            j = j+1
            m = 1
        ELSE
          IF (m .NE. 0) THEN
            j = j+1
            m = 0
          ENDIF
        ENDIF
      ENDDO
C
C ... Info
C
      CALL OUTPAR(OutFile,' ')
      END
C
C --- Set cross section type (0 - coherent, 1 - incoherent)
C
      SUBROUTINE SETSIGT(T) BIND(C)
      USE ISO_C_BINDING
      INTEGER(C_INT), VALUE :: T
      INCLUDE 'npmain.inc'
      IF (T .EQ. 0) THEN
        SigType = 0
      ELSE
        SigType = 1
      ENDIF
      END
C
C --- Set cross section unit (deafult is 1, i.e. in GeV^{-2})
C
      SUBROUTINE SETSIGU(U) BIND(C)
      USE ISO_C_BINDING
      REAL (C_DOUBLE), VALUE :: U
      INCLUDE 'npmain.inc'
      SigUnit = U
      END
C
C --- Set cross section accuracy (deafult is 0.01)
C
      SUBROUTINE SETSIGA(eps) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE :: eps
      INCLUDE 'npmain.inc'
      SigEps = eps
      END
C
C --- Set kinematic variables
C
C     SETKINE  : KinE              --> init. misc. variables
C
C     SETKINS  : KinE              --> KinSNM, KinSNX
C     SETKINSM : KinE,KinSN        --> KinMuM, KinMuX
C     SETNQXY  : KinE,KinSN,KinMu  --> KinNu , KinQ2, KinX, KinY
C
C     SETKINN  : KinE              --> KinNuM, KinNuX
C     SETKINNQ : KinE,KinNu        --> KinQ2M, KinQ2X
C
C     SETKINQ  : KinE              --> KinQ2M, KinQ2X
C     SETKINQN : KinE,KinQ2        --> KinNuM, KinNuX
C
C     SETKINX  : KinE              --> KinXM , KinXX
C     SETKINXY : KinE,KinX         --> KinYM , KinYX
C
C     SETKINY  : KinE              --> KinYM , KinYX
C     SETKINYX : KinE,KinY         --> KinXM , KinXX
C
C     SETTMX   : KinX,KinY         --> KinTM , KinTX, KinTC, KinPC, KinPNC
C     SETKTL   : ..., KinT         --> KinKT , KinKL, KinKLC
C
C ... init. misc. variables
C
      SUBROUTINE SETKINE
      INCLUDE 'npmain.inc'
      REAL*8   vc,gc,Ec,Ecn,Eclmax,Q20,Q2A,nu0,nuA,nuB,x0
      SAVE     vc,gc,Ec,Ecn,Eclmax,Q20,Q2A,nu0,nuA,nuB,x0
      REAL*8   Fpl,FEm,FEp,FQm,FQp,El,Ecl,ksi,mu,x
      REAL*8   a,b,c,aa,bb,cc,dd
      Fpl(x) = SQRT(MAX(x**2-KinML**2,0D0))
      FEm(x) = gc*(x-vc*Fpl(x))
      FEp(x) = gc*(x+vc*Fpl(x))
      FQm(x) = 2*Ecn*(x-Fpl(x))-KinML**2
      FQp(x) = 2*Ecn*(x+Fpl(x))-KinML**2
      IF (KinE .LE. KinEM) THEN
        CALL NPwarnD('SETKINE: Too low neutrino energy = ',KinE)
        STOP
      ENDIF
      vc     = KinE/(KinE+KinMN)
      gc     = 1/SQRT(1-vc*vc)
      Ec     = (KinE+KinMN)/gc
      Ecn    =  KinE*KinMN /Ec
      Eclmax = (Ec**2+KinML**2-(KinMP+KinMN)**2)/(2*Ec)
      Q20    = (2*Ecn*Eclmax-KinML**2)/(2*KinMN)
      Q2A    = FQm(KinML)
      nu0    = KinE-gc*Eclmax
      nuA    = KinE-FEp(KinML)
      nuB    = KinE-FEm(Eclmax)
      x0     = KinML*(2*Ecn-KinML)/(2*KinMN*(KinE-gc*KinML))
      RETURN
C
C ... s_N range
C
      ENTRY      SETKINS
      KinSNM = (KinMN+KinMP)**2
      KinSNX = (Ec   -KinML)**2
      RETURN
C
C ... s_N - fixed, mu range
C
      ENTRY      SETKINSM
      Ecl    = (Ec**2+KinML**2-KinSN)/(2*Ec)
      KinMuX =  KinE/Ec*Fpl(Ecl)
      KinMuM = -KinMuX
      RETURN
C
C ... s_N, mu --> nu, Q^2, x, y
C
      ENTRY      SETNQXY
      Ecl   = (Ec**2+KinML**2-KinSN)/(2*Ec)
      KinNu = KinE-gc*Ecl-KinMu
      KinQ2 = 2*( Ecn*Ecl-KinMu*KinMN)-KinML**2
      KinX  = KinQ2/(2*KinMN*KinNu)
      KinY  = KinNu/KinE
      RETURN
C
C ... Q^2 range
C
      ENTRY      SETKINQ
      KinQ2M = FQm(Eclmax)
      KinQ2X = FQp(Eclmax)
      RETURN
C
C ... Q^2 - fixed, nu range
C
      ENTRY      SETKINQN
      ksi    = (KinQ2 +KinML**2)/(2*Ecn)
      Ecl    = (ksi**2+KinML**2)/(2*ksi)
      KinNuM = KinE-gc*(Eclmax+vc*(Eclmax-ksi))
      IF (KinQ2 .LE. Q2A) THEN
        KinNuX = KinE-FEp(Ecl)
      ELSE
        KinNuX = KinE-FEm(Ecl)
      ENDIF
      RETURN
C
C ... nu range
C
      ENTRY      SETKINN
      KinNuM = KinE-FEp(    Eclmax          )
      KinNuX = KinE-FEm(MIN(Eclmax,gc*KinML))
      RETURN
C
C ... nu - fixed, Q^2 range
C
      ENTRY      SETKINNQ
      El = KinE-KinNu
      IF (KinNu .LE. nuA) THEN
        KinQ2M = FQm(FEm(El))
      ELSE
        KinQ2M = FQp(FEm(El))
      ENDIF
      IF ((KinNu .LE. nuB) .OR. (nuA .LE. nuB)) THEN
        KinQ2X = 2*Ecn*(Eclmax-(El/gc-Eclmax)/vc)-KinML**2
      ELSE
        KinQ2X = FQp(FEp(El))
      ENDIF
      RETURN
C
C ... x range
C
      ENTRY      SETKINX
      mu    = KinE/Ec*SQRT(Eclmax**2-KinML**2)
      KinXM = (Q20-mu)/(nu0-mu)
      KinXX = (Q20+mu)/(nu0+mu)
      RETURN
C
C ... x - fixed, y range
C
      ENTRY      SETKINXY
      mu    = (Q20-nu0*KinX)/(1-KinX)
      KinYM = (KinE-gc*Eclmax-mu)/KinE
      a     = 2*KinMN*KinX*gc+2*Ecn
      b     = 2*KinMN*KinX*KinE+KinML**2
      c     = 2*KinMN*(1-KinX)*KinE/Ec
      aa    = (a-c)*(a+c)
      bb    = 2*a*b
      cc    = b*b+(c*KinML)**2
      dd    = SQRT(bb*bb-4*aa*cc)
      Ecl   = (bb-dd)/(2*aa)
      mu    = KinE/Ec*SQRT(Ecl**2-KinML**2)
      IF (KinX .LE. x0) THEN
        KinYX = (KinE-gc*Ecl-mu)/KinE
      ELSE
        KinYX = (KinE-gc*Ecl+mu)/KinE
      ENDIF
      RETURN
C
C ... y range
C
      ENTRY      SETKINY
      KinYM = 1-FEp(    Eclmax          )/KinE
      KinYX = 1-FEm(MIN(Eclmax,gc*KinML))/KinE
      RETURN
C
C ... y - fixed, x range
C
      ENTRY      SETKINYX
      KinNu = KinE*KinY
      El    = KinE-KinNu
      IF (KinNu .LE. nuA) THEN
        KinQ2M = FQm(FEm(El))
      ELSE
        KinQ2M = FQp(FEm(El))
      ENDIF
      IF ((KinNu .LE. nuB) .OR. (nuA .LE. nuB)) THEN
        KinQ2X = 2*Ecn*(Eclmax-(El/gc-Eclmax)/vc)-KinML**2
      ELSE
        KinQ2X = FQp(FEp(El))
      ENDIF
      KinXM = KinQ2M/(2*KinMN*KinNu)
      KinXX = KinQ2X/(2*KinMN*KinNu)
      RETURN
      END
C
C ... t range, t_c, p_c, p_N^c
C
      SUBROUTINE SETTMX
      INCLUDE 'npmain.inc'
      REAL*8   vc,gc,ENcp,nu,Q2,q,E,Ec,ENc,dt
      SAVE     vc,gc,ENcp
      nu     = KinE*KinY
      Q2     = 2*KinMN*nu*KinX
      E      = nu+KinMN
      q      = SQRT(nu**2+Q2)
      vc     = q/E
      gc     = 1/SQRT(1-vc**2)
      Ec     = E/gc
      ENc    = gc*KinMN
      KinPNC = vc*ENc
      ENcp   = (Ec**2+KinMN**2-KinMP**2)/(2*Ec)
      KinPC  = SQRT(ABS(ENcp**2-KinMN**2))
      dt     = 2*KinPNC*KinPC
      KinTC  = 2*(KinMN**2-ENc*ENcp)
      KinTM  = KinTC-dt
      KinTX  = KinTC+dt
      RETURN
C
C ... k_T, k_L, k_L^c
C
      ENTRY      SETKTL
      KinKLC = (KinTC-KinT)/(2*KinPNC)
      KinKL  = gc*(KinKLC+vc*ENcp)
      KinKT  = SQRT(ABS(KinPC**2-KinKLC**2))
      END
C
C --- User's functions for kinematic
C
C     KINEMIN  :            --> Emin
C
C     KINEN    : E          --> numin, numax
C     KINENQ   : E, nu      --> Q2min, Q2max
C     KINEQ    : E          --> Q2min, Q2max
C     KINEQN   : E, Q2      --> numin, numax
C
C     KINEX    : E          --> xmin , xmax
C     KINEXY   : E, x       --> ymin , ymax
C     KINEY    : E          --> ymin , ymax
C     KINEYX   : E, y       --> xmin , xmax
C
C     KINENQT  : E, nu, Q2  --> tmin , tmax
C     KINEXYT  : E, x , y   --> tmin , tmax
C
C     KINEXYNQ : E, x , y   --> nu, Q2
C     KINENQXY : E, nu, Q2  --> x , y
C
      SUBROUTINE KINEMIN(Emin) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE) :: Emin
      INCLUDE 'npmain.inc'
      Emin = KinEM
      END

      SUBROUTINE KINEN(E,numin,numax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E
      REAL(C_DOUBLE), INTENT(OUT) ::  numin,numax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINN
      numin = KinNuM
      numax = KinNuX
      END

      SUBROUTINE KINENQ(E,nu,Q2min,Q2max) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E,nu
      REAL(C_DOUBLE), INTENT(OUT) :: Q2min,Q2max
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINN
      IF ((nu .LT. KinNuM) .OR. (KinNuX .LT. nu)) THEN
        CALL NPwarnD('KINENQ: Out of region: nu = ',nu)
        STOP
      ENDIF
      KinNu = nu
      CALL SETKINNQ
      Q2min = KinQ2M
      Q2max = KinQ2X
      END

      SUBROUTINE KINEQ(E,Q2min,Q2max) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) ::  E
      REAL(C_DOUBLE), INTENT(OUT) :: Q2min,Q2max
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINQ
      Q2min = KinQ2M
      Q2max = KinQ2X
      END

      SUBROUTINE KINEQN(E,Q2,numin,numax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E,Q2
      REAL(C_DOUBLE), INTENT(OUT) :: numin,numax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINQ
      IF ((Q2 .LT. KinQ2M) .OR. (KinQ2X .LT. Q2)) THEN
        CALL NPwarnD('KINEQN: Out of region: Q2 = ',Q2)
        STOP
      ENDIF
      KinQ2 = Q2
      CALL SETKINQN
      numin = KinNuM
      numax = KinNuX
      END

      SUBROUTINE KINENQT(E,nu,Q2,tmin,tmax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) ::  E,nu,Q2
      REAL(C_DOUBLE), INTENT(OUT) :: tmin,tmax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINN
      IF ((nu .LT. KinNuM) .OR. (KinNuX .LT. nu)) THEN
        CALL NPwarnD('KINENQT: Out of region: nu = ',nu)
        STOP
      ENDIF
      KinNu = nu
      CALL SETKINNQ
      IF ((Q2 .LT. KinQ2M) .OR. (KinQ2X .LT. Q2)) THEN
        CALL NPwarnD('KINENQT: Out of region: Q2 = ',Q2)
        STOP
      ENDIF
      KinQ2 = Q2
      KinX  = KinQ2/(2*KinMN*KinNu)
      KinY  = KinNu/KinE
      CALL SETTMX
      tmin  = KinTM
      tmax  = KinTX
      END

      SUBROUTINE KINEX(E,xmin,xmax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E
      REAL(C_DOUBLE), INTENT(OUT) :: xmin,xmax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINX
      xmin = KinXM
      xmax = KinXX
      END

      SUBROUTINE KINEXY(E,x,ymin,ymax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E,x
      REAL(C_DOUBLE), INTENT(OUT) :: ymin,ymax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINX
      IF ((x .LT. KinXM) .OR. (KinXX .LT. x)) THEN
        CALL NPwarnD('KINEXY: Out of region: x = ',x)
        STOP
      ENDIF
      KinX = x
      CALL SETKINXY
      ymin = KinYM
      ymax = KinYX
      END

      SUBROUTINE KINEY(E,ymin,ymax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E
      REAL(C_DOUBLE), INTENT(OUT) :: ymin,ymax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINY
      ymin = KinYM
      ymax = KinYX
      END

      SUBROUTINE KINEYX(E,y,xmin,xmax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E,y
      REAL(C_DOUBLE), INTENT(OUT) :: xmin,xmax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINY
      IF ((y .LT. KinYM) .OR. (KinYX .LT. y)) THEN
        CALL NPwarnD('KINEYX: Out of region: y = ',y)
        STOP
      ENDIF
      KinY = y
      CALL SETKINYX
      xmin = KinXM
      xmax = KinXX
      END

      SUBROUTINE KINEXYT(E,x,y,tmin,tmax) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E,x,y
      REAL(C_DOUBLE), INTENT(OUT) :: tmin,tmax
      INCLUDE 'npmain.inc'
      KinE = E
      CALL SETKINE
      CALL SETKINX
      IF ((x .LT. KinXM) .OR. (KinXX .LT. x)) THEN
        CALL NPwarnD('KINEXYT: Out of region: x = ',x)
        STOP
      ENDIF
      KinX = x
      CALL SETKINXY
      IF ((y .LT. KinYM) .OR. (KinYX .LT. y)) THEN
        CALL NPwarnD('KINEXYT: Out of region: y = ',y)
        STOP
      ENDIF
      KinY = Y
      CALL SETTMX
      tmin = KinTM
      tmax = KinTX
      END

      SUBROUTINE KINEXYNQ(E,x,y,nu,Q2) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E,x,y
      REAL(C_DOUBLE), INTENT(OUT) :: nu,Q2
      INCLUDE 'npmain.inc'
      nu = y*E
      Q2 = 2*KinMN*nu*x
      END

      SUBROUTINE KINENQXY(E,nu,Q2,x,y) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E,nu,Q2
      REAL(C_DOUBLE), INTENT(OUT) :: x,y
      INCLUDE 'npmain.inc'
      x = Q2/(2*KinMN*nu)
      y = nu/E
      END
C
C --- SigNA - total cross section sigma(nu A --> l pi A)
C
      FUNCTION SigNA(E) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E
      REAL*8   SigNA
      INCLUDE 'npmain.inc'
      INCLUDE 'nppion.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAs
      EXTERNAL SigNAs
      KinE = E
      CALL SETKINE
      CALL SETKINS
      IF ((KinSNM .LT. sarn) .AND.
     .    (KinSNX .GT. sarn)) THEN
        SigNA = DINTAD(SigNAs,KinSNM,sarn  ,-SigEps,10)
     .        + DINTAD(SigNAs,sarn  ,KinSNX,-SigEps,10)
      ELSE
        SigNA = DINTAD(SigNAs,KinSNM,KinSNX,-SigEps,10)
      ENDIF
      END

      FUNCTION SigNAs(sN) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) ::   sN
      REAL*8   SigNAs
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAm
      EXTERNAL SigNAm
      KinSN  = sN
      CALL SETKINSM
      SigNAs = DINTAD2(SigNAm,KinMuM,KinMuX,-SigEps,10)
      END

      
      FUNCTION SigNAm(mu) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: mu
      REAL*8   SigNAm
      INCLUDE 'npmain.inc'
      REAL*8   SigNAtt
      EXTERNAL SigNAtt
      KinMu  = mu
      CALL SETNQXY
      SigNAm = SigNAtt()/(2*KinMN*KinE*KinNu)
      END

      FUNCTION SigNAtt() BIND(C)
      USE ISO_C_BINDING
      REAL*8   SigNAtt
      INCLUDE 'npmain.inc'
      INCLUDE 'nppion.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAt,tm
      EXTERNAL SigNAt
      CALL SETTMX
      IF ((SigType .EQ. 1) .AND. (KinSN .GE. sarn)) THEN
        tm = KinTC
      ELSE
        tm = KinTM
      ENDIF
      SigNAtt = DINTAD3(SigNAt,tm,KinTX,-SigEps,10)
      END
      
      FUNCTION SigNAt(t) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: t
      REAL*8   SigNAt
      INCLUDE 'npmain.inc'
      REAL*8   SigNAsmt
      EXTERNAL SigNAsmt
      SigNAt = SigNAsmt(KinE,KinX,KinY,t)
      END
C
C --- SigNAnu - differential cross section d sigma / d nu
C
      FUNCTION SigNAnu(E,nu) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E, nu
      REAL*8   SigNAnu
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAnQ
      EXTERNAL SigNAnQ
      KinE  = E
      KinNu = nu
      CALL SETKINE
      CALL SETKINN
      IF ((KinNu .LT. KinNuM)  .OR.
     .    (KinNu .GT. KinNuX)) THEN
        SigNAnu = 0D0
        RETURN
      ENDIF
      CALL SETKINNQ
      SigNAnu = DINTAD2(SigNAnQ,KinQ2M,KinQ2X,-SigEps,10)
      END

      FUNCTION SigNAnQ(Q2) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: Q2
      REAL*8   SigNAnQ
      INCLUDE 'npmain.inc'
      INCLUDE 'nppion.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAtt
      EXTERNAL SigNAtt
      KinQ2   = Q2
      KinSN   = (KinMN+2*KinNu)*KinMN-KinQ2
      KinX    = KinQ2/(2*KinMN*KinNu)
      KinY    = KinNu/KinE
      SigNAnQ = SigNAtt()/(2*KinMN*KinE*KinNu)
      END
C
C --- SigNAQ2 - differential cross section d sigma / d Q^2
C
      FUNCTION SigNAQ2(E,Q2) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E, Q2
      REAL*8   SigNAQ2
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAQn
      EXTERNAL SigNAQn
      KinE  = E
      KinQ2 = Q2
      CALL SETKINE
      CALL SETKINQ
      IF ((KinQ2 .LT. KinQ2M)  .OR.
     .    (KinQ2 .GT. KinQ2X)) THEN
        SigNAQ2 = 0D0
        RETURN
      ENDIF
      CALL SETKINQN
      SigNAQ2 = DINTAD2(SigNAQn,KinNuM,KinNuX,-SigEps,10)
      END

      FUNCTION SigNAQn(nu) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: nu
      REAL*8   SigNAQn
      INCLUDE 'npmain.inc'
      INCLUDE 'nppion.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAtt
      EXTERNAL SigNAtt
      KinNu   = nu
      KinSN   = (KinMN+2*KinNu)*KinMN-KinQ2
      KinX    = KinQ2/(2*KinMN*KinNu)
      KinY    = KinNu/KinE
      SigNAQn = SigNAtt()/(2*KinMN*KinE*KinNu)
      END
C
C --- SigNAx - differential cross section d sigma / d x
C
      FUNCTION SigNAx(E,x) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E, x
      REAL*8   SigNAx
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAxx
      EXTERNAL SigNAxx
      KinE = E
      KinX = x
      CALL SETKINE
      CALL SETKINX
      IF ((KinX .LT. KinXM)  .OR.
     .    (KinX .GT. KinXX)) THEN
        SigNAx = 0D0
        RETURN
      ENDIF
      CALL SETKINXY
      SigNAx = DINTAD2(SigNAxx,KinYM,KinYX,-SigEps,10)
      END

      FUNCTION SigNAxx(y) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: y
      REAL*8   SigNAxx
      INCLUDE 'npmain.inc'
      INCLUDE 'nppion.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAtt
      EXTERNAL SigNAtt
      KinY    = y
      KinNu   = KinY*KinE
      KinQ2   = KinX*(2*KinNu*KinMN)
      KinSN   = (KinMN+2*KinNu)*KinMN-KinQ2
      SigNAxx = SigNAtt()
      END
C
C --- SigNAy - differential cross section d sigma / d y
C
      FUNCTION SigNAy(E,y) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E, y
      REAL*8   SigNAy
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAyy
      EXTERNAL SigNAyy
      KinE = E
      KinY = y
      CALL SETKINE
      CALL SETKINY
      IF ((KinY .LT. KinYM)  .OR.
     .    (KinY .GT. KinYX)) THEN
        SigNAy = 0D0
        RETURN
      ENDIF
      CALL SETKINYX
      SigNAy = DINTAD2(SigNAyy,KinXM,KinXX,-SigEps,10)
      END

      FUNCTION SigNAyy(x) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: x
      REAL*8   SigNAyy
      INCLUDE 'npmain.inc'
      INCLUDE 'nppion.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAtt
      EXTERNAL SigNAtt
      KinX    = x
      KinNu   = KinY*KinE
      KinQ2   = KinX*(2*KinNu*KinMN)
      KinSN   = (KinMN+2*KinNu)*KinMN-KinQ2
      SigNAyy = SigNAtt()
      END
C
C --- SigNAxy - differential cross section d sigma / dx dy
C
      FUNCTION SigNAxy(E,x,y) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E, x, y
      REAL*8   SigNAxy
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.fun'
      REAL*8   SigNAtt
      EXTERNAL SigNAtt
      KinE = E
      KinX = x
      KinY = y
      CALL SETKINE
      CALL SETKINX
      IF ((KinX .LT. KinXM)  .OR.
     .    (KinX .GT. KinXX)) THEN
        SigNAxy = 0D0
        RETURN
      ENDIF
      CALL SETKINXY
      IF ((KinY .LT. KinYM)  .OR.
     .    (KinY .GT. KinYX)) THEN
        SigNAxy = 0D0
        RETURN
      ENDIF
      KinNu   = KinY*KinE
      KinQ2   = KinX*(2*KinNu*KinMN)
      KinSN   = (KinMN+2*KinNu)*KinMN-KinQ2
      SigNAxy = SigNAtt()
      END
C
C --- SigNAxyt - differential cross section:
C
C       d sigma (nu A --> l pi A)
C       -------------------------
C               dx dy dt
C
C     E      - neutrino energy
C     x,y,t  - kinematic variables
C
      FUNCTION SigNAxyt(E,x,y,t)

      REAL*8   E, x, y, t
      REAL*8   SigNAxyt
      REAL*8   SigNAsmt
      INCLUDE 'npmain.inc'
      INCLUDE 'npmath.inc'
      REAL*8   sN,Sig,Phi
      SAVE     sN,Sig
      DATA     sN /0D0/
      REAL*8   SigPNtot,SigPNel,DSigPNel,AlpPN,BslPN
      EXTERNAL SigPNtot,SigPNel,DSigPNel,AlpPN,BslPN

      KinE = E
      KinX = x
      KinY = y

      CALL SETKINE
      CALL SETKINY
      IF     (KinY .LT. KinYM) THEN
              KinY  =   KinYM
      ELSEIF (KinY .GT. KinYX) THEN
              KinY  =   KinYX
      ENDIF
      CALL SETKINYX
      IF     (KinX .LT. KinXM) THEN
              KinX  =   KinXM
      ELSEIF (KinX .GT. KinXX) THEN
              KinX  =   KinXX
      ENDIF

      KinNu = KinY*KinE
      KinQ2 = KinX*(2*KinMN*KinNu)
      KinSN = (KinMN+2*KinNu)*KinMN-KinQ2

      CALL SETTMX

      ENTRY SigNAsmt(E,x,y,t)
      KinT = t

      IF     (KinT .LT. KinTM) THEN
              KinT  =   KinTM
      ELSEIF (KinT .GT. KinTX) THEN
              KinT  =   KinTX
      ENDIF

      CALL SETKTL

      IF (KinPC .LT. 1D-7) THEN
        SigNAxyt = 0D0
        RETURN
      ENDIF

      IF (ABS(sN-KinSN) .GE. 1D-7) THEN
        sN      = KinSN
        SigBpN  = BslPN(sN)
        CALL SETBPN(SigBpN)
        SigTot  = SigPNtot(sN)
        SigEl   = SigPNel (sN)
        SigTotM = SigTot*(1D0-AlpPN(sN)*(0D0,1D0))
        SigIn   = SigTotM-sigEl
        Sig     = (DREAL(SigTotM)**2+DIMAG(SigTotM)**2)/(16*Pi)
      ENDIF

      IF ((NuclZ .EQ. 1) .AND. (NuclA .EQ. 1)) THEN
        Phi = DSigPNel(KinSN,-KinKLC/KinPC)/ABS(2*KinPC*KinPNC)
      ELSE
        IF (SigType .EQ. 0) THEN
          CALL PhiCoh(Phi)
          Phi = Phi*Sig*ABS(KinKLC/KinPNC)
        ELSE
          CALL PhiInc(Phi)
          Phi = Phi*DSigPNel(KinSN,-KinKLC/KinPC)/ABS(2*KinPC*KinPNC)
        ENDIF
      ENDIF
      SigNAxyt = SigFact*KinMN*(KinE-KinNu)/(1+KinQ2/Massa1**2)**2*Phi
      SigNAxyt = SigNAxyt/SigUnit
      END

      FUNCTION SigNAnQt(E,nu,Q2,t) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: E, nu, Q2, t
      REAL*8   SigNAnQt
      INCLUDE 'npmain.inc'
      REAL*8   SigNAxyt,x,y
      EXTERNAL SigNAxyt
      x        = Q2/(2*KinMN*nu)
      y        = nu/E
      SigNAnQt = SigNAxyt(E,x,y,t)/(2*KinMN*E*nu)
      END
C
C --- Coherent
C
      SUBROUTINE PhiCoh(PhiC) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), INTENT(OUT) :: PhiC
      INCLUDE   'npmain.inc'
      INCLUDE   'npmath.inc'
      INCLUDE   'npmath.fun'
      COMPLEX*16     Phi,Phi0,alp
      COMPLEX*16     PhiP(RHORNUM)
      COMPLEX*16     PhiN(RHORNUM)
      COMPLEX*16     PhiB(RHORNUM)
      COMMON /PhiCB/ PhiB
      COMPLEX*16     PhiCohB
      EXTERNAL       PhiCohB
      INTEGER        i,j,f
      REAL*8         R
      alp = -0.5D0*SigTotM
      DO i = 1,RHORNUM
        Phi0 = EXP(-SigTotM*RhoTdatB(i,1))
        DO j = 1,RHORNUM
          Phi     = EXP(alp*RhoTdatB(i,j))
          R       =         RhoMdatB(i,j)
          PhiP(j) = R*Phi
          PhiN(j) = R/Phi*Phi0
        ENDDO
        PhiB(i) = CSIFOU(RHORNUM,PhiP,RhoH, KinKL)
     .          + CSIFOU(RHORNUM,PhiN,RhoH,-KinKL)
      ENDDO
      f = RhoH*KinKT*2
      IF (f .GE. 1) THEN
        Phi = CINTEG(PhiCohB,0D0,RhoRmax,1D-4*RhoC,f*(RHORNUM-1)/2)
      ELSE
        DO i = 1,RHORNUM
           PhiB(i) = RhoRdat(i)*DFUNJ0(RhoRdat(i)*KinKT)*PhiB(i)
        ENDDO
        Phi = CSISUM(RHORNUM,PhiB,RhoH)
      ENDIF
      PhiC = (2*Pi*ABS(Phi))**2
      END

      FUNCTION   PhiCohB(b) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: b
      COMPLEX*16 PhiCohB
      INCLUDE   'npmain.inc'
      INCLUDE   'npmath.fun'
      COMPLEX*16     PhiB(RHORNUM)
      COMMON /PhiCB/ PhiB
      PhiCohB = b*DFUNJ0(b*KinKT)*CINTER(RHORNUM,RhoRdat,PhiB,b)
      END

C
C --- Incoherent
C
      SUBROUTINE PhiInc(PhiI) BIND(C)
      USE ISO_C_BINDING
      REAL(C_DOUBLE), INTENT(OUT) :: PhiI

      INCLUDE   'npmain.inc'
      INCLUDE   'npmath.inc'
      INCLUDE   'npmath.fun'
      COMPLEX*16 Phi,c1
      REAL*8     R,T,c2,cT
      COMPLEX*16 alp1,p1       ,alp2,p2       ,alp3,p3
      COMPLEX*16 Phi1P(RHORNUM),Phi2P(RHORNUM),Phi3P(RHORNUM)
      COMPLEX*16 Phi1N(RHORNUM),Phi2N(RHORNUM),Phi3N(RHORNUM)
      COMPLEX*16 Phi1(-RHORNUM+1:RHORNUM-1)
      COMPLEX*16 Phi2(-RHORNUM+1:RHORNUM-1)
      COMPLEX*16 PhiS(-RHORNUM+1:RHORNUM-1)
      COMPLEX*16 PhiZ (RHORNUM)
      DATA       PhiZ (RHORNUM)/(0D0,0D0)/
      REAL*8     PhiB (RHORNUM)
      INTEGER    i,j,k,j1,j1n,j1x,j2,j2n,NS
      alp1 = - SigTotM     *0.5D0
      alp2 = -(SigIn-SigEl)*0.5D0
      alp3 = - SigIn
      cT   = ABS(SigTotM)**2
      c1   = (SigIn-SigEl)*SigTotM/SigEl
      c2   = -cT/(4*SigEl)
      DO i = 1,RHORNUM
        T  = RhoTdatB(i,1)
        p1 = EXP(2*alp1*T)
        p2 = EXP(2*alp2*T)
        p3 = EXP(2*alp3*T)
        DO j = 1,RHORNUM
          R        = RhoMdatB(i,j)
          T        = RhoTdatB(i,j)
          Phi      = EXP(alp1*T)
          Phi1P(j) = R*Phi
          Phi1N(j) = R/Phi*p1
          Phi      = EXP(alp2*T)
          Phi2P(j) = R*Phi
          Phi2N(j) = R/Phi*p2
          Phi      = EXP(alp3*T)
          Phi3P(j) = R*Phi
          Phi3N(j) = R/Phi*p3
          k        = j-1
          Phi1( k) = Phi1P(j)
          Phi1(-k) = Phi1N(j)
          Phi2( k) = Phi2P(j)
          Phi2(-k) = Phi2N(j)
        ENDDO
        j1n  = -(RHORNUM-1)
        j1x  =   RHORNUM-1
        j2n  = j1n
        NS   = 2*RHORNUM-1
        DO j = 1,RHORNUM-1
          j2 = j2n
          DO j1 = j1n,j1x
             PhiS(j1) = Phi1(j1)*Phi2(j2)
             j2       = j2+1
          ENDDO
          PhiZ(j) = CSISUM(NS,PhiS,RhoH)
          j1x     = j1x-2
          j2n     = j2n+2
          NS      = NS -2
        ENDDO
        p1 = CSIFOU(RHORNUM,PhiZ ,RhoH, KinKL*2)
        p2 = CSIFOU(RHORNUM,Phi1P,RhoH, KinKL)
     .     + CSIFOU(RHORNUM,Phi1N,RhoH,-KinKL)
        p3 = CSISUM(RHORNUM,Phi3P,RhoH)
     .     + CSISUM(RHORNUM,Phi3N,RhoH)
        PhiB(i) = RhoRdat(i)*(ABS(p1*c1)+c2*ABS(p2)**2+ABS(p3))
      ENDDO
      PhiI = 2*Pi*DSISUM(RHORNUM,PhiB,RhoH)
      END
