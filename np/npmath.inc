C
C --- Constants
C
      REAL*8     Pi
      PARAMETER (Pi = 3.14159265358979D0)
C
C --- Integration
C
      INTEGER    DINTADMAX
      PARAMETER (DINTADMAX = 20000)

      INTEGER          DINTEGnum,DINTEG2num,DINTEG3num
      INTEGER          DINTEGitr,DINTEG2itr,DINTEG3itr
      INTEGER          DINTEGmax,DINTEG2max,DINTEG3max

      INTEGER          CINTEGnum,CINTEG2num,CINTEG3num
      INTEGER          CINTEGitr,CINTEG2itr,CINTEG3itr
      INTEGER          CINTEGmax,CINTEG2max,CINTEG3max

      INTEGER          DINTADnum,DINTAD2num,DINTAD3num

      COMMON /NPINTEG/ DINTEGnum,DINTEG2num,DINTEG3num
     .                ,DINTEGitr,DINTEG2itr,DINTEG3itr
     .                ,DINTEGmax,DINTEG2max,DINTEG3max
     .                ,CINTEGnum,CINTEG2num,CINTEG3num
     .                ,CINTEGitr,CINTEG2itr,CINTEG3itr
     .                ,CINTEGmax,CINTEG2max,CINTEG3max
     .                ,DINTADnum,DINTAD2num,DINTAD3num
