C-----------------------------------------------------------
C     Neutrinoproduction of Pions: Math Library
C-----------------------------------------------------------
C
C --- Default values
C
      BLOCK DATA NPMATH
      INCLUDE   'npmath.inc'
      DATA DINTEGnum /20/, DINTEG2num /20/, DINTEG3num /20/
      DATA DINTEGmax /16/, DINTEG2max /16/, DINTEG3max /16/
      DATA CINTEGnum /20/, CINTEG2num /20/, CINTEG3num /20/
      DATA CINTEGmax /16/, CINTEG2max /16/, CINTEG3max /16/
      DATA DINTADnum / 5/, DINTAD2num / 5/, DINTAD3num / 5/
      END
C
C --- Spline initialization
C
C  Input:  N  - Number of points
C          X  - array x(i), i = 1..N
C          Y  - array y(i), i = 1..N
C          Y1 - first derivative y'(x(1))
C          Yn - first derivative y'(x(N))
C
C  If first derivative (Y1 and/or YN) >= 10^100 then
C  natural spline is used for corresponding boundary
C
C  Output: W  - array w(i), i = 1..N
C          U  - array u(i), i = 1..N - temporary
C
      SUBROUTINE DSPLINEI(N,X,Y,Y1,YN,W,U)
      INTEGER    N
      REAL*8     X(N)
      REAL*8     Y(N),Y1,YN
      REAL*8     W(N)
      REAL*8     U(N)
      REAL*8     s,p,x0m,xpm,xp0
      INTEGER    i
C
C ... Check input
C
      IF (N .LT. 3) THEN
        CALL NPwarnI('DSPLINEI: not enough points N = ',N)
        STOP
      ENDIF
C
C ... Calculate W(i)
C
      IF (Y1 .GE. 1D100) THEN
        W(1) = 0D0
        U(1) = 0D0
      ELSE
        U(1) = (3D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-Y1)
        W(1) = -0.5D0
      ENDIF
      DO i = 2,N-1
        x0m  = X(i  )-X(i-1)
        xpm  = X(i+1)-X(i-1)
        xp0  = X(i+1)-X(i  )
        s    = x0m/xpm
        p    = s*W(i-1)+2D0
        W(i) = (s-1D0)/p
        U(i) = (Y(i+1)-Y(i))/xp0-(Y(i)-Y(i-1))/x0m
        U(i) = (6.0*U(i)/xpm-s*U(i-1))/p
      ENDDO
      IF (YN .GE. 1D100) THEN
        s = 0D0
        p = 0D0
      ELSE
        s = 0.5D0
        p = (3D0/(X(N)-X(N-1)))*(YN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      W(N) = (p-s*U(N-1))/(s*W(N-1)+1D0)
      DO i = N-1,1,-1
        W(i) = W(i)*W(i+1)+U(i)
      ENDDO
      END
C
C --- Spline
C
C  Input:  N  - Number of points
C          X  - array x(i), i = 1..N
C          Y  - array y(i), i = 1..N
C          W  - array w(i), i = 1..N
C          XX - point x
C
      FUNCTION DSPLINE(N,X,Y,W,XX)
      REAL*8   DSPLINE
      INTEGER  N
      REAL*8   X(N)
      REAL*8   Y(N)
      REAL*8   W(N)
      REAL*8   XX
      INTEGER  l,m,u
      REAL*8   a,b,h
C
C ... Check input
C
      IF (N .LT. 3) THEN
        CALL NPwarnI('DSPLINE: not enough points N = ',N)
        STOP
      ENDIF
C
C ... Search index on x
C
      l = 0
      u = N+1
      DO WHILE (u-l .GT. 1)
        m = (l+u)/2
        IF (XX .LT. X(m)) THEN
          u = m
        ELSE
          l = m
        ENDIF
      ENDDO
      IF (l .LE. 0) THEN
        IF (X(1)-XX .GT. 1D-7) THEN
          CALL NPwarnD('DSPLINE: value below range = ',X)
        ENDIF
        l = 1
      ENDIF
      IF (l .GE. N) THEN
        IF (XX-X(N) .GT. 1D-7) THEN
          CALL NPwarnD('DSPLINE: value above range = ',X)
        ENDIF
        l = N-1
      ENDIF
      u = l+1
C
C ... Result
C
      h       =  X(u)-X(l)
      a       = (X(u)-XX)/h
      b       = 1D0-a
      DSPLINE = a*Y(l)+b*Y(u)+((a*a*a-a)*W(l)+(b*b*b-b)*W(u))*(h*h)/6D0
      END
C
C --- Interpolation, search index
C
C     N  - number of x points
C     XX - array  of x_i (i = 1..N, monotonically increasing)
C     X  - x point
C
      FUNCTION DINTERind(N,XX,X)
      INTEGER  DINTERind
      INTEGER     N
      REAL*8   XX(N)
      REAL*8      X
      INTEGER  l,m,u
      l = 0
      u = N+1
      DO WHILE (u-l .GT. 1)
        m = (l+u)/2
        IF (X .LT. XX(m)) THEN
          u = m
        ELSE
          l = m
        ENDIF
      ENDDO
      IF (l .LE. 0) THEN
        IF (XX(1)-X .GT. 1D-7) THEN
          CALL NPwarnD('DINTERind: value below range = ',X)
        ENDIF
        l = 1
      ENDIF
      IF (l .GE. N) THEN
        IF (X-XX(N) .GT. 1D-7) THEN
          CALL NPwarnD('DINTERind: value above range = ',X)
        ENDIF
        l = N-1
      ENDIF
      IF (MOD(l,2) .EQ. 0  ) l = l-1
      IF (l        .GE. N-1) l = l-1
      DINTERind = l
      END
C
C --- Interpolation function (one variable)
C
C      N  - number of   x points
C      XX - array  of   x_i (i = 1..N, monotonically increasing)
C      YY - array  of y[x_i]
C      X  - x point
C
      FUNCTION DINTER(N,XX,YY,X)
      REAL*8   DINTER
      INTEGER     N
      REAL*8   XX(N)
      REAL*8   YY(N)
      REAL*8      X
      INCLUDE 'npmath.inc'
      INTEGER  i1,i2,i3
      REAL*8   v1,v2,v3,v12,v13,v23,vv1,vv2,vv3
      INTEGER  DINTERind
      EXTERNAL DINTERind
C
C ... Check input
C
      IF (N .LT. 3) THEN
        CALL NPwarnI('DINTER: not enough points N = ',N)
        STOP
      ENDIF
C
C ... Search index on x
C
      i1 = DINTERind(N,XX,X)
      i2 = i1+1
      i3 = i2+1
C
C ... Interpolation
C
      v1     =   X-XX(i1)
      v2     =   X-XX(i2)
      v3     =   X-XX(i3)
      v12    =  v2-v1
      v13    =  v3-v1
      v23    =  v3-v2
      vv1    =  v2*v3/(v12*v13)
      vv2    = -v1*v3/(v12*v23)
      vv3    =  v1*v2/(v13*v23)
      DINTER = vv1*YY(i1)+vv2*YY(i2)+vv3*YY(i3)
      END
C
C --- Interpolation function (two variables)
C
C      N1  - number of x1 points
C      N2  - number of x2 points
C      XX1 - array  of x1_i (i = 1..N1, monotonically increasing)
C      XX2 - array  of x2_i (i = 1..N2, monotonically increasing)
C      YY  - array  of y[x1_i,x2_j]
C      X1  - x1 point
C      X2  - x2 point
C
      FUNCTION DINTER2(N1,N2,XX1,XX2,YY,X1,X2)
      REAL*8   DINTER2
      INTEGER      N1 ,    N2
      REAL*8   XX1(N1),XX2(N2)
      REAL*8   YY (N1 ,    N2)
      REAL*8       X1 ,    X2
      INCLUDE 'npmath.inc'
      INTEGER  i1,i2,i3
      INTEGER  j1,j2,j3
      REAL*8   v1,v2,v3,v12,v13,v23,vv1,vv2,vv3
      REAL*8   y1,y2,y3
      INTEGER  DINTERind
      EXTERNAL DINTERind
C
C ... Check input
C
      IF (N1 .LT. 3) THEN
        CALL NPwarnI('DINTER2: not enough points N1 = ',N1)
        STOP
      ENDIF
      IF (N2 .LT. 3) THEN
        CALL NPwarnI('DINTER2: not enough points N2 = ',N2)
        STOP
      ENDIF
C
C ... Search index on x1
C
      i1 = DINTERind(N1,XX1,X1)
      i2 = i1+1
      i3 = i2+1
C
C ... Search index on x2
C
      j1 = DINTERind(N2,XX2,X2)
      j2 = j1+1
      j3 = j2+1
C
C ... Interpolation
C
      v1      =  X2-XX2(j1)
      v2      =  X2-XX2(j2)
      v3      =  X2-XX2(j3)
      v12     =  v2-v1
      v13     =  v3-v1
      v23     =  v3-v2
      vv1     =  v2*v3/(v12*v13)
      vv2     = -v1*v3/(v12*v23)
      vv3     =  v1*v2/(v13*v23)
      y1      = vv1*YY(i1,j1)+vv2*YY(i1,j2)+vv3*YY(i1,j3)
      y2      = vv1*YY(i2,j1)+vv2*YY(i2,j2)+vv3*YY(i2,j3)
      y3      = vv1*YY(i3,j1)+vv2*YY(i3,j2)+vv3*YY(i3,j3)
      v1      =  X1-XX1(i1)
      v2      =  X1-XX1(i2)
      v3      =  X1-XX1(i3)
      v12     =  v2-v1
      v13     =  v3-v1
      v23     =  v3-v2
      vv1     =  v2*v3/(v12*v13)
      vv2     = -v1*v3/(v12*v23)
      vv3     =  v1*v2/(v13*v23)
      DINTER2 = vv1*y1+vv2*y2+vv3*y3
      END
C
C --- Interpolation function (three variables)
C
C      N1  - number of x1 points
C      N2  - number of x2 points
C      N3  - number of x3 points
C      XX1 - array  of x1_i (i = 1..N1, monotonically increasing)
C      XX2 - array  of x2_i (i = 1..N2, monotonically increasing)
C      XX3 - array  of x3_i (i = 1..N3, monotonically increasing)
C      YY  - array  of y[x1_i,x2_j,x3_k]
C      X1  - x1 point
C      X2  - x2 point
C      X3  - x3 point
C
      FUNCTION DINTER3(N1,N2,N3,XX1,XX2,XX3,YY,X1,X2,X3)
      REAL*8   DINTER3
      INTEGER      N1 ,    N2 ,    N3
      REAL*8   XX1(N1),XX2(N2),XX3(N3)
      REAL*8   YY (N1 ,    N2 ,    N3)
      REAL*8       X1 ,    X2 ,    X3
      INCLUDE 'npmath.inc'
      INTEGER  i1,i2,i3
      INTEGER  j1,j2,j3
      INTEGER  k1,k2,k3
      REAL*8   v1,v2,v3,v12,v13,v23,vv1,vv2,vv3
      REAL*8   y1,y2,y3,y11,y21,y31,y12,y22,y32,y13,y23,y33
      INTEGER  DINTERind
      EXTERNAL DINTERind
C
C ... Check input
C
      IF (N1 .LT. 3) THEN
        CALL NPwarnI('DINTER3: not enough points N1 = ',N1)
        STOP
      ENDIF
      IF (N2 .LT. 3) THEN
        CALL NPwarnI('DINTER3: not enough points N2 = ',N2)
        STOP
      ENDIF
      IF (N3 .LT. 3) THEN
        CALL NPwarnI('DINTER3: not enough points N3 = ',N3)
        STOP
      ENDIF
C
C ... Search index on x1
C
      i1 = DINTERind(N1,XX1,X1)
      i2 = i1+1
      i3 = i2+1
C
C ... Search index on x2
C
      j1 = DINTERind(N2,XX2,X2)
      j2 = j1+1
      j3 = j2+1
C
C ... Search index on x3
C
      k1 = DINTERind(N3,XX3,X3)
      k2 = k1+1
      k3 = k2+1
C
C ... Interpolation
C
      v1      =  X3-XX3(k1)
      v2      =  X3-XX3(k2)
      v3      =  X3-XX3(k3)
      v12     =  v2-v1
      v13     =  v3-v1
      v23     =  v3-v2
      vv1     =  v2*v3/(v12*v13)
      vv2     = -v1*v3/(v12*v23)
      vv3     =  v1*v2/(v13*v23)
      y11     = vv1*YY(i1,j1,k1)+vv2*YY(i1,j1,k2)+vv3*YY(i1,j1,k3)
      y21     = vv1*YY(i2,j1,k1)+vv2*YY(i2,j1,k2)+vv3*YY(i2,j1,k3)
      y31     = vv1*YY(i3,j1,k1)+vv2*YY(i3,j1,k2)+vv3*YY(i3,j1,k3)
      y12     = vv1*YY(i1,j2,k1)+vv2*YY(i1,j2,k2)+vv3*YY(i1,j2,k3)
      y22     = vv1*YY(i2,j2,k1)+vv2*YY(i2,j2,k2)+vv3*YY(i2,j2,k3)
      y32     = vv1*YY(i3,j2,k1)+vv2*YY(i3,j2,k2)+vv3*YY(i3,j2,k3)
      y13     = vv1*YY(i1,j3,k1)+vv2*YY(i1,j3,k2)+vv3*YY(i1,j3,k3)
      y23     = vv1*YY(i2,j3,k1)+vv2*YY(i2,j3,k2)+vv3*YY(i2,j3,k3)
      y33     = vv1*YY(i3,j3,k1)+vv2*YY(i3,j3,k2)+vv3*YY(i3,j3,k3)
      v1      =  X2-XX2(j1)
      v2      =  X2-XX2(j2)
      v3      =  X2-XX2(j3)
      v12     =  v2-v1
      v13     =  v3-v1
      v23     =  v3-v2
      vv1     =  v2*v3/(v12*v13)
      vv2     = -v1*v3/(v12*v23)
      vv3     =  v1*v2/(v13*v23)
      y1      = vv1*y11+vv2*y12+vv3*y13
      y2      = vv1*y21+vv2*y22+vv3*y23
      y3      = vv1*y31+vv2*y32+vv3*y33
      v1      =  X1-XX1(i1)
      v2      =  X1-XX1(i2)
      v3      =  X1-XX1(i3)
      v12     =  v2-v1
      v13     =  v3-v1
      v23     =  v3-v2
      vv1     =  v2*v3/(v12*v13)
      vv2     = -v1*v3/(v12*v23)
      vv3     =  v1*v2/(v13*v23)
      DINTER3 = vv1*y1+vv2*y2+vv3*y3
      END
C
C --- Interpolation function (one variable), complex
C
C      N  - number of   x points
C      XX - array  of   x_i (i = 1..N, monotonically increasing)
C      YY - array  of y[x_i] (complex)
C      X  - x point
C
      FUNCTION   CINTER(N,XX,YY,X)
      COMPLEX*16 CINTER
      INTEGER       N
      REAL*8     XX(N)
      COMPLEX*16 YY(N)
      REAL*8        X
      INCLUDE   'npmath.inc'
      INTEGER    i1,i2,i3
      REAL*8     v1,v2,v3,v12,v13,v23,vv1,vv2,vv3
      INTEGER    DINTERind
      EXTERNAL   DINTERind
C
C ... Check input
C
      IF (N .LT. 3) THEN
        CALL NPwarnI('CINTER: not enough points N = ',N)
        STOP
      ENDIF
C
C ... Search index on x
C
      i1 = DINTERind(N,XX,X)
      i2 = i1+1
      i3 = i2+1
C
C ... Interpolation
C
      v1     =   X-XX(i1)
      v2     =   X-XX(i2)
      v3     =   X-XX(i3)
      v12    =  v2-v1
      v13    =  v3-v1
      v23    =  v3-v2
      vv1    =  v2*v3/(v12*v13)
      vv2    = -v1*v3/(v12*v23)
      vv3    =  v1*v2/(v13*v23)
      CINTER = vv1*YY(i1)+vv2*YY(i2)+vv3*YY(i3)
      END
C
C --- Simpson Summation
C
C     N should be odd and >= 3
C
C ... Double
C
      FUNCTION DSISUM(N,F,H)
      REAL*8   DSISUM
      INTEGER  N,i
      REAL*8   F(N),S1,S2
      REAL*8   H
      S1 = 0D0
      S2 = 0D0
      DO i = 2,N-1,2
        S1 = S1+F(i)
      ENDDO
      IF (N .GT. 3) THEN
        DO i = 3,N-2,2
          S2 = S2+F(i)
        ENDDO
      ENDIF
      DSISUM = (F(1)+F(N)+4D0*S1+2D0*S2)*H/3D0
      END
C
C ... Complex
C
      FUNCTION   CSISUM(N,F,H)
      COMPLEX*16 CSISUM
      INTEGER    N,i
      COMPLEX*16 F(N),S1,S2
      REAL*8     H
      S1 = (0D0,0D0)
      S2 = (0D0,0D0)
      DO i = 2,N-1,2
        S1 = S1+F(i)
      ENDDO
      IF (N .GT. 3) THEN
        DO i = 3,N-2,2
          S2 = S2+F(i)
        ENDDO
      ENDIF
      CSISUM = (F(1)+F(N)+4D0*S1+2D0*S2)*H/3D0
      END
C
C  --- Simpson summation for Fourier integral:
C
C      xmax   iwx
C       | dx e   f(x), where xmax = h*(N-1)
C       0
C
      FUNCTION   CSIFOU(N,F,h,w)
      COMPLEX*16 CSIFOU
      INTEGER    N
      COMPLEX*16 F(N),G(N)
      REAL*8     h,w
      COMPLEX*16 FC,FS,Fm,Fp,F0,Fa,Fb
      REAL*8     hw,p,p1,p2,p3
      REAL*8     Ch,Ch0,Ch2
      REAL*8     Sh,Sh0,Sh2
      COMPLEX*16 CSISUM
      EXTERNAL   CSISUM
      INTEGER    i
      hw = h*w
      Ch = COS(hw)
      Sh = SIN(hw)
      IF (ABS(hw) .LT. 1D-2) THEN
        Ch0 = 1
        Sh0 = 0
        DO i = 1,N
          G(i) = F(i)*(Ch0+(0,1)*Sh0)
          p    = Ch0*Ch-Sh0*Sh
          Sh0  = Sh0*Ch+Ch0*Sh
          Ch0  = p
        ENDDO
        CSIFOU = CSISUM(N,G,h)
        RETURN
      ENDIF
      Ch0 = Ch
      Sh0 = Sh
      Ch2 = Ch*Ch-Sh*Sh
      Sh2 = Sh*Ch*2
      p1  = 2*hw*Ch
      p2  =   hw*hw-2
      p3  =  (hw*Ch-Sh)*hw
      FC  = (0,0)
      FS  = (0,0)
      DO i = 1,N-2,2
        Fm  = F(i)
        F0  = F(i+1)
        Fp  = F(i+2)
        Fa  = Fp+Fm
        Fa  = (Fa-2*F0)*p1+(Fa*p2+4*F0)*Sh
        Fb  = (Fp-Fm)*p3
        FC  = FC+Fa*Ch0+Fb*Sh0
        FS  = FS+Fa*Sh0-Fb*Ch0
        p   = Ch0*Ch2-Sh0*Sh2
        Sh0 = Sh0*Ch2+Ch0*Sh2
        Ch0 = p
      ENDDO
      CSIFOU = (FC+(0,1)*FS)/(hw*hw*w)
      END
C
C --- Integration (Simpson)
C
C     FUN(x) - function
C     A      - low   limit
C     B      - upper limit
C     EPS    - accuracy
C     NUM    - number of intervals
C
C ... Double, first
C
      FUNCTION DINTEG(FUN,A,B,EPS,NUM)
      REAL*8   DINTEG,FUN,A,B,EPS
      EXTERNAL FUN
      INTEGER  NUM
      INCLUDE 'npmath.inc'
      INTEGER  n,i,j
      REAL*8   h,a1,a2,f0,f1,f2,u,In,Io
      IF (NUM .LT. 1) THEN
        n = DINTEGnum
      ELSE
        n = NUM
      ENDIF
      h  = (B-A)/n
      a1 = A+h*0.5D0
      a2 = A
      f0 = FUN(A)+FUN(B)
      f1 = FUN(a1)
      f2 = 0D0
      DO i = 2,N
        a1 = a1+h
        a2 = a2+h
        f1 = f1+FUN(a1)
        f2 = f2+FUN(a2)
      ENDDO
      j  = 0
      Io = 0
 1    In = h/6D0*(f0+4D0*f1+2D0*f2)
      IF (EPS .EQ. 0D0) GOTO 2
      IF (j   .GT. 0  ) THEN
        u = ABS(In-Io)
        IF (u   .EQ. 0D0) GOTO 2
        IF (EPS .GT. 0D0) THEN
          u = u/(ABS(In)+1D0)
        ELSE
          u = u/ ABS(In)
        ENDIF
        IF (u .LE. ABS(EPS) ) GOTO 2
        IF (j .GE. DINTEGmax) THEN
          CALL NPwarnI('DINTEG: iteration limit = ',DINTEGmax)
          GOTO 2
        ENDIF
      ENDIF
      j  = j+1
      Io = In
      n  = n*2
      h  = h*0.5D0
      a1 = h*0.5D0+A
      f2 = f2+f1
      f1 = 0D0
      DO i = 1,N
        f1 = f1+FUN(a1)
        a1 = a1+h
      ENDDO
      GOTO 1
 2    DINTEGitr = j
      DINTEG    = In
      END
C
C ... Double, second
C
      FUNCTION DINTEG2(FUN,A,B,EPS,NUM)
      REAL*8   DINTEG2,FUN,A,B,EPS
      EXTERNAL FUN
      INTEGER  NUM
      INCLUDE 'npmath.inc'
      INTEGER  n,i,j
      REAL*8   h,a1,a2,f0,f1,f2,u,In,Io
      IF (NUM .LT. 1) THEN
        n = DINTEG2num
      ELSE
        n = NUM
      ENDIF
      h  = (B-A)/n
      a1 = A+h*0.5D0
      a2 = A
      f0 = FUN(A)+FUN(B)
      f1 = FUN(a1)
      f2 = 0D0
      DO i = 2,N
        a1 = a1+h
        a2 = a2+h
        f1 = f1+FUN(a1)
        f2 = f2+FUN(a2)
      ENDDO
      j  = 0
      Io = 0
 1    In = h/6D0*(f0+4D0*f1+2D0*f2)
      IF (EPS .EQ. 0D0) GOTO 2
      IF (j   .GT. 0  ) THEN
        u = ABS(In-Io)
        IF (u   .EQ. 0D0) GOTO 2
        IF (EPS .GT. 0D0) THEN
          u = u/(ABS(In)+1D0)
        ELSE
          u = u/ ABS(In)
        ENDIF
        IF (u .LE. ABS(EPS)  ) GOTO 2
        IF (j .GE. DINTEG2max) THEN
          CALL NPwarnI('DINTEG2: iteration limit = ',DINTEG2max)
          GOTO 2
        ENDIF
      ENDIF
      j  = j+1
      Io = In
      n  = n*2
      h  = h*0.5D0
      a1 = h*0.5D0+A
      f2 = f2+f1
      f1 = 0D0
      DO i = 1,N
        f1 = f1+FUN(a1)
        a1 = a1+h
      ENDDO
      GOTO 1
 2    DINTEG2itr = j
      DINTEG2    = In
      END
C
C ... Double, third
C
      FUNCTION DINTEG3(FUN,A,B,EPS,NUM)
      REAL*8   DINTEG3,FUN,A,B,EPS
      EXTERNAL FUN
      INTEGER  NUM
      INCLUDE 'npmath.inc'
      INTEGER  n,i,j
      REAL*8   h,a1,a2,f0,f1,f2,u,In,Io
      IF (NUM .LT. 1) THEN
        n = DINTEG3num
      ELSE
        n = NUM
      ENDIF
      h  = (B-A)/n
      a1 = A+h*0.5D0
      a2 = A
      f0 = FUN(A)+FUN(B)
      f1 = FUN(a1)
      f2 = 0D0
      DO i = 2,N
        a1 = a1+h
        a2 = a2+h
        f1 = f1+FUN(a1)
        f2 = f2+FUN(a2)
      ENDDO
      j  = 0
      Io = 0
 1    In = h/6D0*(f0+4D0*f1+2D0*f2)
      IF (EPS .EQ. 0D0) GOTO 2
      IF (j   .GT. 0  ) THEN
        u = ABS(In-Io)
        IF (u   .EQ. 0D0) GOTO 2
        IF (EPS .GT. 0D0) THEN
          u = u/(ABS(In)+1D0)
        ELSE
          u = u/ ABS(In)
        ENDIF
        IF (u .LE. ABS(EPS)  ) GOTO 2
        IF (j .GE. DINTEG3max) THEN
          CALL NPwarnI('DINTEG3: iteration limit = ',DINTEG3max)
          GOTO 2
        ENDIF
      ENDIF
      j  = j+1
      Io = In
      n  = n*2
      h  = h*0.5D0
      a1 = h*0.5D0+A
      f2 = f2+f1
      f1 = 0D0
      DO i = 1,N
        f1 = f1+FUN(a1)
        a1 = a1+h
      ENDDO
      GOTO 1
 2    DINTEG3itr = j
      DINTEG3    = In
      END
C
C ... Complex, first
C
      FUNCTION   CINTEG(FUN,A,B,EPS,NUM)
      COMPLEX*16 CINTEG
      COMPLEX*16 FUN
      EXTERNAL   FUN
      REAL*8     A,B,EPS
      INTEGER    NUM
      INCLUDE   'npmath.inc'
      INTEGER    n ,i ,j
      REAL*8     h ,a1,a2,u
      COMPLEX*16 f0,f1,f2,In,Io
      IF (NUM .LT. 1) THEN
        n = CINTEGnum
      ELSE
        n = NUM
      ENDIF
      h  = (B-A)/n
      a1 = A+h*0.5D0
      a2 = A
      f0 = FUN(A)+FUN(B)
      f1 = FUN(a1)
      f2 = (0D0,0D0)
      DO i = 2,N
        a1 = a1+h
        a2 = a2+h
        f1 = f1+FUN(a1)
        f2 = f2+FUN(a2)
      ENDDO
      j  = 0
 1    In = h/6D0*(f0+4D0*f1+2D0*f2)
      IF (EPS .EQ. 0D0) GOTO 2
      IF (j   .GT. 0  ) THEN
        u = ABS(In-Io)
        IF (u   .EQ. 0D0) GOTO 2
        IF (EPS .GT. 0D0) THEN
          u = u/(ABS(In)+1D0)
        ELSE
          u = u/ ABS(In)
        ENDIF
        IF (u .LE. ABS(EPS) ) GOTO 2
        IF (j .GE. CINTEGmax) THEN
          CALL NPwarnI('CINTEG: iteration limit = ',CINTEGmax)
          GOTO 2
        ENDIF
      ENDIF
      j  = j+1
      Io = In
      n  = n*2
      h  = h*0.5D0
      a1 = h*0.5D0+A
      f2 = f2+f1
      f1 = (0D0,0D0)
      DO i = 1,N
        f1 = f1+FUN(a1)
        a1 = a1+h
      ENDDO
      GOTO 1
 2    CINTEGitr = j
      CINTEG    = In
      END
C
C ... Complex, second
C
      FUNCTION   CINTEG2(FUN,A,B,EPS,NUM)
      COMPLEX*16 CINTEG2
      COMPLEX*16 FUN
      EXTERNAL   FUN
      REAL*8     A,B,EPS
      INTEGER    NUM
      INCLUDE   'npmath.inc'
      INTEGER    n ,i ,j
      REAL*8     h ,a1,a2,u
      COMPLEX*16 f0,f1,f2,In,Io
      IF (NUM .LT. 1) THEN
        n = CINTEG2num
      ELSE
        n = NUM
      ENDIF
      h  = (B-A)/n
      a1 = A+h*0.5D0
      a2 = A
      f0 = FUN(A)+FUN(B)
      f1 = FUN(a1)
      f2 = (0D0,0D0)
      DO i = 2,N
        a1 = a1+h
        a2 = a2+h
        f1 = f1+FUN(a1)
        f2 = f2+FUN(a2)
      ENDDO
      j  = 0
 1    In = h/6D0*(f0+4D0*f1+2D0*f2)
      IF (EPS .EQ. 0D0) GOTO 2
      IF (j   .GT. 0  ) THEN
        u = ABS(In-Io)
        IF (u   .EQ. 0D0) GOTO 2
        IF (EPS .GT. 0D0) THEN
          u = u/(ABS(In)+1D0)
        ELSE
          u = u/ ABS(In)
        ENDIF
        IF (u .LE. ABS(EPS)  ) GOTO 2
        IF (j .GE. CINTEG2max) THEN
          CALL NPwarnI('CINTEG2: iteration limit = ',CINTEG2max)
          GOTO 2
        ENDIF
      ENDIF
      j  = j+1
      Io = In
      n  = n*2
      h  = h*0.5D0
      a1 = h*0.5D0+A
      f2 = f2+f1
      f1 = (0D0,0D0)
      DO i = 1,N
        f1 = f1+FUN(a1)
        a1 = a1+h
      ENDDO
      GOTO 1
 2    CINTEG2itr = j
      CINTEG2    = In
      END
C
C ... Complex, third
C
      FUNCTION   CINTEG3(FUN,A,B,EPS,NUM)
      COMPLEX*16 CINTEG3
      COMPLEX*16 FUN
      EXTERNAL   FUN
      REAL*8     A,B,EPS
      INTEGER    NUM
      INCLUDE   'npmath.inc'
      INTEGER    n ,i ,j
      REAL*8     h ,a1,a2,u
      COMPLEX*16 f0,f1,f2,In,Io
      IF (NUM .LT. 1) THEN
        n = CINTEG3num
      ELSE
        n = NUM
      ENDIF
      h  = (B-A)/n
      a1 = A+h*0.5D0
      a2 = A
      f0 = FUN(A)+FUN(B)
      f1 = FUN(a1)
      f2 = (0D0,0D0)
      DO i = 2,N
        a1 = a1+h
        a2 = a2+h
        f1 = f1+FUN(a1)
        f2 = f2+FUN(a2)
      ENDDO
      j  = 0
 1    In = h/6D0*(f0+4D0*f1+2D0*f2)
      IF (EPS .EQ. 0D0) GOTO 2
      IF (j   .GT. 0  ) THEN
        u = ABS(In-Io)
        IF (u   .EQ. 0D0) GOTO 2
        IF (EPS .GT. 0D0) THEN
          u = u/(ABS(In)+1D0)
        ELSE
          u = u/ ABS(In)
        ENDIF
        IF (u .LE. ABS(EPS)  ) GOTO 2
        IF (j .GE. CINTEG3max) THEN
          CALL NPwarnI('CINTEG3: iteration limit = ',CINTEG3max)
          GOTO 2
        ENDIF
      ENDIF
      j  = j+1
      Io = In
      n  = n*2
      h  = h*0.5D0
      a1 = h*0.5D0+A
      f2 = f2+f1
      f1 = (0D0,0D0)
      DO i = 1,N
        f1 = f1+FUN(a1)
        a1 = a1+h
      ENDDO
      GOTO 1
 2    CINTEG3itr = j
      CINTEG3    = In
      END
C
C --- Integration (Adaptive Simpson)
C
C     FUN(x) - function
C     A      - low   limit
C     B      - upper limit
C     EPS    - accuracy
C     NUM    - number of intervals
C
C ... Double, first
C
      FUNCTION DINTAD(FUN,A,B,EPS,NUM)
      REAL*8   DINTAD,FUN,A,B,EPS
      EXTERNAL FUN
      INTEGER  NUM
      INCLUDE 'npmath.inc'
      INTEGER  p(0:DINTADMAX)
      REAL*8   x(0:DINTADMAX)
      REAL*8   f(0:DINTADMAX)
      REAL*8   V(  DINTADMAX)
      REAL*8   D(  DINTADMAX)
      LOGICAL  M(  DINTADMAX)
      INTEGER  n,i,k
      REAL*8   h,u,vi,L,In,Io
      INTEGER  i0,i1
      REAL*8   x0,x1
      Io = 0D0
      In = 0D0
      IF (NUM .LT. 1) THEN
        n = 2*DINTADnum
      ELSE
        n = 2*NUM
      ENDIF
      IF (n .GE. DINTADMAX) GOTO 9
      L = B-A
      h = L/n
      DO i = 0,n
        p(i) = i
        x(i) = i*h+A
        f(i) = FUN(x(i))
      ENDDO
      h = h/12
      DO i = 1,n,2
        V(i)   = h*(5*f(i-1)+8*f(i)-f(i+1))
        V(i+1) = h*(5*f(i+1)+8*f(i)-f(i-1))
        In     = In+V(i)+V(i+1)
        M(i)   = .TRUE.
        M(i+1) = .TRUE.
      ENDDO
      IF ((EPS .EQ. 0D0) .OR. (In .EQ. 0D0)) GOTO 2
 1    Io = In
      i  = 1
      DO WHILE (i .LE. n)
        i0 = p(i-1)
        i1 = p(i)
        x0 = x(i0)
        x1 = x(i1)
        vi = V(i1)
        h  = x1-x0
        u  = ABS(In*EPS*h/L)
        IF ((ABS(vi).GE.u) .AND. (M(i1) .OR. (ABS(D(i1)).GE.u))) THEN
          IF (n .EQ. DINTADMAX) GOTO 9
          n     = n+1
          h     = h/24
          x(n)  = (x1+x0)/2
          f(n)  = FUN(x(n))
          V(n)  = h*(5*f(i0)+8*f(n)-f(i1))
          V(i1) = h*(5*f(i1)+8*f(n)-f(i0))
          u     = V(i1)+V(n)-vi
          In    = In+u
          D(i1) = u/2
          D(n)  = u/2
          M(i1) = .FALSE.
          M(n)  = .FALSE.
          DO k = n,i+1,-1
            p(k) = p(k-1)
          ENDDO
          p(i) = n
          i    = i+1 
        ENDIF
        i = i+1
      ENDDO
      u = ABS(In-Io)
      IF (u   .EQ. 0D0) GOTO 2
      IF (EPS .GT. 0D0) THEN
        u = u/(ABS(In)+1D0)
      ELSE
        u = u/ ABS(In)
      ENDIF
      IF (u .LE. ABS(EPS)) GOTO 2
      GOTO 1
 9    CALL NPwarn('DINTAD: not enough iterations...')
 2    DINTAD = In
      END
C
C --- Integration (Adaptive Simpson)
C
C     FUN(x) - function
C     A      - low   limit
C     B      - upper limit
C     EPS    - accuracy
C     NUM    - number of intervals
C
C ... Double, second
C
      FUNCTION DINTAD2(FUN,A,B,EPS,NUM)
      REAL*8   DINTAD2,FUN,A,B,EPS
      EXTERNAL FUN
      INTEGER  NUM
      INCLUDE 'npmath.inc'
      INTEGER  p(0:DINTADMAX)
      REAL*8   x(0:DINTADMAX)
      REAL*8   f(0:DINTADMAX)
      REAL*8   V(  DINTADMAX)
      REAL*8   D(  DINTADMAX)
      LOGICAL  M(  DINTADMAX)
      INTEGER  n,i,k
      REAL*8   h,u,vi,L,In,Io
      INTEGER  i0,i1
      REAL*8   x0,x1
      Io = 0D0
      In = 0D0
      IF (NUM .LT. 1) THEN
        n = 2*DINTADnum
      ELSE
        n = 2*NUM
      ENDIF
      IF (n .GE. DINTADMAX) GOTO 9
      L = B-A
      h = L/n
      DO i = 0,n
        p(i) = i
        x(i) = i*h+A
        f(i) = FUN(x(i))
      ENDDO
      h = h/12
      DO i = 1,n,2
        V(i)   = h*(5*f(i-1)+8*f(i)-f(i+1))
        V(i+1) = h*(5*f(i+1)+8*f(i)-f(i-1))
        In     = In+V(i)+V(i+1)
        M(i)   = .TRUE.
        M(i+1) = .TRUE.
      ENDDO
      IF ((EPS .EQ. 0D0) .OR. (In .EQ. 0D0)) GOTO 2
 1    Io = In
      i  = 1
      DO WHILE (i .LE. n)
        i0 = p(i-1)
        i1 = p(i)
        x0 = x(i0)
        x1 = x(i1)
        vi = V(i1)
        h  = x1-x0
        u  = ABS(In*EPS*h/L)
        IF ((ABS(vi).GE.u) .AND. (M(i1) .OR. (ABS(D(i1)).GE.u))) THEN
          IF (n .EQ. DINTADMAX) GOTO 9
          n     = n+1
          h     = h/24
          x(n)  = (x1+x0)/2
          f(n)  = FUN(x(n))
          V(n)  = h*(5*f(i0)+8*f(n)-f(i1))
          V(i1) = h*(5*f(i1)+8*f(n)-f(i0))
          u     = V(i1)+V(n)-vi
          In    = In+u
          D(i1) = u/2
          D(n)  = u/2
          M(i1) = .FALSE.
          M(n)  = .FALSE.
          DO k = n,i+1,-1
            p(k) = p(k-1)
          ENDDO
          p(i) = n
          i    = i+1
        ENDIF
        i = i+1
      ENDDO
      u = ABS(In-Io)
      IF (u   .EQ. 0D0) GOTO 2
      IF (EPS .GT. 0D0) THEN
        u = u/(ABS(In)+1D0)
      ELSE
        u = u/ ABS(In)
      ENDIF
      IF (u .LE. ABS(EPS)) GOTO 2
      GOTO 1
 9    CALL NPwarn('DINTAD2: not enough iterations...')
 2    DINTAD2 = In
      END
C
C --- Integration (Adaptive Simpson)
C
C     FUN(x) - function
C     A      - low   limit
C     B      - upper limit
C     EPS    - accuracy
C     NUM    - number of intervals
C
C ... Double, third
C
      FUNCTION DINTAD3(FUN,A,B,EPS,NUM)
      REAL*8   DINTAD3,FUN,A,B,EPS
      EXTERNAL FUN
      INTEGER  NUM
      INCLUDE 'npmath.inc'
      INTEGER  p(0:DINTADMAX)
      REAL*8   x(0:DINTADMAX)
      REAL*8   f(0:DINTADMAX)
      REAL*8   V(  DINTADMAX)
      REAL*8   D(  DINTADMAX)
      LOGICAL  M(  DINTADMAX)
      INTEGER  n,i,k
      REAL*8   h,u,vi,L,In,Io
      INTEGER  i0,i1
      REAL*8   x0,x1
      Io = 0D0
      In = 0D0
      IF (NUM .LT. 1) THEN
        n = 2*DINTADnum
      ELSE
        n = 2*NUM
      ENDIF
      IF (n .GE. DINTADMAX) GOTO 9
      L = B-A
      h = L/n
      DO i = 0,n
        p(i) = i
        x(i) = i*h+A
        f(i) = FUN(x(i))
      ENDDO
      h = h/12
      DO i = 1,n,2
        V(i)   = h*(5*f(i-1)+8*f(i)-f(i+1))
        V(i+1) = h*(5*f(i+1)+8*f(i)-f(i-1))
        In     = In+V(i)+V(i+1)
        M(i)   = .TRUE.
        M(i+1) = .TRUE.
      ENDDO
      IF ((EPS .EQ. 0D0) .OR. (In .EQ. 0D0)) GOTO 2
 1    Io = In
      i  = 1
      DO WHILE (i .LE. n)
        i0 = p(i-1)
        i1 = p(i)
        x0 = x(i0)
        x1 = x(i1)
        vi = V(i1)
        h  = x1-x0
        u  = ABS(In*EPS*h/L)
        IF ((ABS(vi).GE.u) .AND. (M(i1) .OR. (ABS(D(i1)).GE.u))) THEN
          IF (n .EQ. DINTADMAX) GOTO 9
          n     = n+1
          h     = h/24
          x(n)  = (x1+x0)/2
          f(n)  = FUN(x(n))
          V(n)  = h*(5*f(i0)+8*f(n)-f(i1))
          V(i1) = h*(5*f(i1)+8*f(n)-f(i0))
          u     = V(i1)+V(n)-vi
          In    = In+u
          D(i1) = u/2
          D(n)  = u/2
          M(i1) = .FALSE.
          M(n)  = .FALSE.
          DO k = n,i+1,-1
            p(k) = p(k-1)
          ENDDO
          p(i) = n
          i    = i+1
        ENDIF
        i = i+1
      ENDDO
      u = ABS(In-Io)
      IF (u   .EQ. 0D0) GOTO 2
      IF (EPS .GT. 0D0) THEN
        u = u/(ABS(In)+1D0)
      ELSE
        u = u/ ABS(In)
      ENDIF
      IF (u .LE. ABS(EPS)) GOTO 2
      GOTO 1
 9    CALL NPwarn('DINTAD3: not enough iterations...')
 2    DINTAD3 = In
      END
C
C --- Modified Bessel function of the first kind: ln(I_0(x))
C
      FUNCTION DLOGI0(X)
      REAL*8   DLOGI0,X
      REAL*8   r,t,u
      IF (X .GE. 0D0) THEN
        r =  X
      ELSE
        r = -X
      ENDIF
      IF (r .LE. 3.75D0) THEN
        t = (r / 3.75D0)**2
        u = 1D0+
     .  t*(3.5156229+t*(3.0899424+t*(1.2067492+
     .  t*(0.2659732+t*(0.0360768+t* 0.0045813)))))
        DLOGI0 = LOG(u)
      ELSE
        t = 3.75D0/r
        u = 0.39894228+t*( 0.01328592+t*( 0.00225319+
     .  t*(-0.00157565+t*( 0.00916281+t*(-0.02057706+
     .  t*( 0.02635537+t*(-0.01647633+t*  0.00392377)))))))
        DLOGI0 = r+LOG(u/SQRT(r))
      ENDIF                                                             
      END                                                               
C
C --- Bessel functions J0, J1, Y0, Y1 (from CERN library)
C
      FUNCTION   DFUNJ0(X)
      REAL*8     DFUNJ0,X
      REAL*8     DFUNY0
      REAL*8     DFUNJ1
      REAL*8     DFUNY1
      LOGICAL    LJ0,LY0,LJ1,LY1
      REAL*8     A(0:16),B(0:16),D(0:16),E(0:16)
      COMPLEX*16 C(0:19),F(0:19),I
      REAL*8     ALFA,CE,H,HF,P,PI,PI1,PI3,PI4,R,R8,R32,V,Y,Z1
      REAL*8      B0, B1, B2
      COMPLEX*16 CB0,CB1,CB2
      INTEGER    IT

      PARAMETER (I   = (0,1))
      PARAMETER (Z1  = 1, HF = Z1/2, R8 = Z1/8, R32 = Z1/32)
      PARAMETER (CE  = 0.57721 56649 01532 861D0)
      PARAMETER (PI  = 3.14159 26535 89793 238D0)
      PARAMETER (PI1 = 2/PI, PI3 = 3*PI/4, PI4 = PI/4)

      DATA P /0D0/

      DATA A( 0) /+0.15772 79714 74890 120D0/
      DATA A( 1) /-0.00872 34423 52852 221D0/
      DATA A( 2) /+0.26517 86132 03336 810D0/
      DATA A( 3) /-0.37009 49938 72649 779D0/
      DATA A( 4) /+0.15806 71023 32097 261D0/
      DATA A( 5) /-0.03489 37694 11408 885D0/
      DATA A( 6) /+0.00481 91800 69467 605D0/
      DATA A( 7) /-0.00046 06261 66206 275D0/
      DATA A( 8) /+0.00003 24603 28821 005D0/
      DATA A( 9) /-0.00000 17619 46907 762D0/
      DATA A(10) /+0.00000 00760 81635 924D0/
      DATA A(11) /-0.00000 00026 79253 531D0/
      DATA A(12) /+0.00000 00000 78486 963D0/
      DATA A(13) /-0.00000 00000 01943 835D0/
      DATA A(14) /+0.00000 00000 00041 253D0/
      DATA A(15) /-0.00000 00000 00000 759D0/
      DATA A(16) /+0.00000 00000 00000 012D0/

      DATA B( 0) /-0.02150 51114 49657 551D0/
      DATA B( 1) /-0.27511 81330 43518 791D0/
      DATA B( 2) /+0.19860 56347 02554 156D0/
      DATA B( 3) /+0.23425 27461 09021 802D0/
      DATA B( 4) /-0.16563 59817 13650 413D0/
      DATA B( 5) /+0.04462 13795 40669 282D0/
      DATA B( 6) /-0.00693 22862 91523 188D0/
      DATA B( 7) /+0.00071 91174 03752 303D0/
      DATA B( 8) /-0.00005 39250 79722 939D0/
      DATA B( 9) /+0.00000 30764 93288 108D0/
      DATA B(10) /-0.00000 01384 57181 230D0/
      DATA B(11) /+0.00000 00050 51054 369D0/
      DATA B(12) /-0.00000 00001 52582 850D0/
      DATA B(13) /+0.00000 00000 03882 867D0/
      DATA B(14) /-0.00000 00000 00084 429D0/
      DATA B(15) /+0.00000 00000 00001 587D0/
      DATA B(16) /-0.00000 00000 00000 026D0/

      DATA
     + C( 0)/ (+0.99898 80898 58965 153D0, -0.01233 15205 78544 144D0)/,
     1 C( 1)/ (-0.00133 84285 49971 856D0, -0.01224 94962 81259 475D0)/,
     2 C( 2)/ (-0.00031 87898 78061 893D0, +0.00009 64941 84993 423D0)/,
     3 C( 3)/ (+0.00000 85112 32210 657D0, +0.00001 36555 70490 357D0)/,
     4 C( 4)/ (+0.00000 06915 42349 139D0, -0.00000 08518 06644 426D0)/,
     5 C( 5)/ (-0.00000 00907 70101 537D0, -0.00000 00272 44053 414D0)/,
     6 C( 6)/ (+0.00000 00014 54928 079D0, +0.00000 00096 46421 338D0)/,
     7 C( 7)/ (+0.00000 00009 26762 487D0, -0.00000 00006 83347 518D0)/,
     8 C( 8)/ (-0.00000 00001 39166 198D0, -0.00000 00000 60627 380D0)/,
     9 C( 9)/ (+0.00000 00000 03237 975D0, +0.00000 00000 21695 716D0)/,
     * C(10)/ (+0.00000 00000 02535 357D0, -0.00000 00000 02304 899D0)/
      DATA
     A C(11)/ (-0.00000 00000 00559 090D0, -0.00000 00000 00122 554D0)/,
     B C(12)/ (+0.00000 00000 00041 919D0, +0.00000 00000 00092 314D0)/,
     C C(13)/ (+0.00000 00000 00008 733D0, -0.00000 00000 00016 778D0)/,
     D C(14)/ (-0.00000 00000 00003 619D0, +0.00000 00000 00000 754D0)/,
     E C(15)/ (+0.00000 00000 00000 594D0, +0.00000 00000 00000 462D0)/,
     F C(16)/ (-0.00000 00000 00000 010D0, -0.00000 00000 00000 159D0)/,
     G C(17)/ (-0.00000 00000 00000 024D0, +0.00000 00000 00000 025D0)/,
     H C(18)/ (+0.00000 00000 00000 008D0, +0.00000 00000 00000 000D0)/,
     I C(19)/ (-0.00000 00000 00000 001D0, -0.00000 00000 00000 001D0)/

      DATA D( 0) /+0.64835 87706 05264 921D0/
      DATA D( 1) /-1.19180 11605 41216 873D0/
      DATA D( 2) /+1.28799 40988 57677 620D0/
      DATA D( 3) /-0.66144 39341 34543 253D0/
      DATA D( 4) /+0.17770 91172 39728 283D0/
      DATA D( 5) /-0.02917 55248 06154 208D0/
      DATA D( 6) /+0.00324 02701 82683 857D0/
      DATA D( 7) /-0.00026 04443 89348 581D0/
      DATA D( 8) /+0.00001 58870 19239 932D0/
      DATA D( 9) /-0.00000 07617 58780 540D0/
      DATA D(10) /+0.00000 00294 97070 073D0/
      DATA D(11) /-0.00000 00009 42421 298D0/
      DATA D(12) /+0.00000 00000 25281 237D0/
      DATA D(13) /-0.00000 00000 00577 740D0/
      DATA D(14) /+0.00000 00000 00011 386D0/
      DATA D(15) /-0.00000 00000 00000 196D0/
      DATA D(16) /+0.00000 00000 00000 003D0/

      DATA E( 0) /-0.04017 29465 44414 076D0/
      DATA E( 1) /-0.44444 71476 30558 063D0/
      DATA E( 2) /-0.02271 92444 28417 736D0/
      DATA E( 3) /+0.20664 45410 17490 520D0/
      DATA E( 4) /-0.08667 16970 56948 524D0/
      DATA E( 5) /+0.01763 67030 03163 134D0/
      DATA E( 6) /-0.00223 56192 94485 095D0/
      DATA E( 7) /+0.00019 70623 02701 541D0/
      DATA E( 8) /-0.00001 28858 53299 241D0/
      DATA E( 9) /+0.00000 06528 47952 359D0/
      DATA E(10) /-0.00000 00264 50737 175D0/
      DATA E(11) /+0.00000 00008 78030 117D0/
      DATA E(12) /-0.00000 00000 24343 279D0/
      DATA E(13) /+0.00000 00000 00572 612D0/
      DATA E(14) /-0.00000 00000 00011 578D0/
      DATA E(15) /+0.00000 00000 00000 203D0/
      DATA E(16) /-0.00000 00000 00000 003D0/

      DATA
     + F( 0)/ (+1.00170 22348 53820 996D0, +0.03726 17150 00537 654D0)/,
     1 F( 1)/ (+0.00225 55728 46561 180D0, +0.03714 53224 79807 690D0)/,
     2 F( 2)/ (+0.00054 32164 87508 013D0, -0.00013 72632 38201 907D0)/,
     3 F( 3)/ (-0.00001 11794 61895 408D0, -0.00001 98512 94687 597D0)/,
     4 F( 4)/ (-0.00000 09469 01382 392D0, +0.00000 10700 14057 386D0)/,
     5 F( 5)/ (+0.00000 01110 32677 121D0, +0.00000 00383 05261 714D0)/,
     6 F( 6)/ (-0.00000 00012 94398 927D0, -0.00000 00116 28723 277D0)/,
     7 F( 7)/ (-0.00000 00011 14905 944D0, +0.00000 00007 59733 092D0)/,
     8 F( 8)/ (+0.00000 00001 57637 232D0, +0.00000 00000 75476 075D0)/,
     9 F( 9)/ (-0.00000 00000 02830 457D0, -0.00000 00000 24752 781D0)/
      DATA
     * F(10)/ (-0.00000 00000 02932 169D0, +0.00000 00000 02493 893D0)/,
     A F(11)/ (+0.00000 00000 00617 809D0, +0.00000 00000 00156 198D0)/,
     B F(12)/ (-0.00000 00000 00043 162D0, -0.00000 00000 00103 385D0)/,
     C F(13)/ (-0.00000 00000 00010 133D0, +0.00000 00000 00018 129D0)/,
     D F(14)/ (+0.00000 00000 00003 973D0, -0.00000 00000 00000 709D0)/,
     E F(15)/ (-0.00000 00000 00000 632D0, -0.00000 00000 00000 520D0)/,
     F F(16)/ (+0.00000 00000 00000 006D0, +0.00000 00000 00000 172D0)/,
     G F(17)/ (+0.00000 00000 00000 027D0, -0.00000 00000 00000 026D0)/,
     H F(18)/ (-0.00000 00000 00000 008D0, -0.00000 00000 00000 000D0)/,
     I F(19)/ (+0.00000 00000 00000 001D0, +0.00000 00000 00000 001D0)/

      LJ0 = .TRUE.
      LY0 = .FALSE.
      GOTO 11

      ENTRY DFUNY0(X)

      LJ0 = .FALSE.
      LY0 = .TRUE.
      IF (X .LE. 0) THEN
        P = 0
        CALL NPwarn('DFUNY0: argument < 0')
        GOTO 9
      ENDIF

   11 V = ABS(X)
      IF (V .LT. 8) THEN
        H    = R32*V**2-1
        ALFA = H+H
        B1   = 0
        B2   = 0
        DO 1 IT = 16,0,-1
        B0 = A(IT)+ALFA*B1-B2
        B2 = B1
    1   B1 = B0
        P  = B0-H*B2
        IF (LY0) THEN
          B1 = 0
          B2 = 0
          DO 2 IT = 16,0,-1
          B0 = B(IT)+ALFA*B1-B2
          B2 = B1
    2     B1 = B0
          P  = PI1*(CE+LOG(HF*X))*P+B0-H*B2
        ENDIF
      ELSE
        R    = 1/V
        H    = 10*R-1
        ALFA = H+H
        CB1  = 0
        CB2  = 0
        DO 3 IT = 19,0,-1
        CB0 = C(IT)+ALFA*CB1-CB2
        CB2 = CB1
    3   CB1 = CB0
        CB0 = SQRT(PI1*R)*EXP(I*(V-PI4))*(CB0-H*CB2)
        IF (LJ0) P =   CB0
        IF (LY0) P =-I*CB0
      ENDIF
      GOTO 9

      ENTRY DFUNJ1(X)

      LJ1 = .TRUE.
      LY1 = .FALSE.
      GOTO 12

      ENTRY DFUNY1(X)

      LJ1 = .FALSE.
      LY1 = .TRUE.

      IF (X .LE. 0) THEN
        P = 0
        CALL NPwarn('DFUNY1: argument < 0')
        GOTO 9
      ENDIF

   12 V = ABS(X)
      IF (V .LT. 8) THEN
        Y    = R8*V
        H    = 2*Y**2-1
        ALFA = H+H
        B1   = 0
        B2   = 0
        DO 4 IT = 16,0,-1
        B0 = D(IT)+ALFA*B1-B2
        B2 = B1
    4   B1 = B0
        P  = Y*(B0-H*B2)
        IF (LY1) THEN
          B1 = 0
          B2 = 0
          DO 5 IT = 16,0,-1
          B0 = E(IT)+ALFA*B1-B2
          B2 = B1
    5     B1 = B0
          P  = PI1*((CE+LOG(HF*X))*P-1/X)+Y*(B0-B2)
        ENDIF
      ELSE
        R    = 1/V
        H    = 10*R-1
        ALFA = H+H
        CB1  = 0
        CB2  = 0
        DO 6 IT = 19,0,-1
        CB0 = F(IT)+ALFA*CB1-CB2
        CB2 = CB1
    6   CB1 = CB0
        CB0 = SQRT(PI1*R)*EXP(I*(V-PI3))*(CB0-H*CB2)
        IF (LJ1) P =   CB0
        IF (LY1) P =-I*CB0
      ENDIF
      IF (X .LT. 0) P = -P
    9 CONTINUE
      DFUNJ0 = P
      END
