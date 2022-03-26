C
C --- Output file
C
      INTEGER          OutFile,OutLevel
      COMMON /NPOUT/   OutFile,OutLevel
C
C --- Constants
C
      REAL*8     alphaem                    !! fine-structure constant
      PARAMETER (alphaem  = 0.00729735308D0)
      REAL*8     fermi                      !! fermi GeV^{-1}
      PARAMETER (fermi    = 5.067728853D0)
      REAL*8     mbarn                      !! millibarn GeV^{-2}
      PARAMETER (mbarn    = 2.56818756D0)
      REAL*8     ConstG                     !! weak coupling G GeV^{-2}
      PARAMETER (ConstG   = 1.16639D-5)
      REAL*8     Constfpi                   !! pion coupling f_pi GeV
      PARAMETER (Constfpi = 0.1307D0)
C
C --- Masses
C
      REAL*8     Masse                      !! electron
      PARAMETER (Masse    = 0.511D-3)
      REAL*8     Massmu                     !! muon
      PARAMETER (Massmu   = 0.10566D0)
      REAL*8     Masstau                    !! tau
      PARAMETER (Masstau  = 1.77699D0)
      REAL*8     Masspi                     !! pion (charged)
      PARAMETER (Masspi   = 0.13957D0)
      REAL*8     Masspi0                    !! pion (neutral)
      PARAMETER (Masspi0  = 0.13498D0)
      REAL*8     Massp                      !! proton
      PARAMETER (Massp    = 0.93827D0)
      REAL*8     Massn                      !! neutron
      PARAMETER (Massn    = 0.93957D0)
      REAL*8     Massa1                     !! a_1
      PARAMETER (Massa1   = 1.230D0)
C
C --- Reaction
C
      CHARACTER*40     Reac                 !! Reaction
      CHARACTER*10     Nucl                 !! Nucleus
      INTEGER          NuclZ,NuclA          !! Nucleus charge and weight
      COMMON /NPREAC/  NuclZ,NuclA,Nucl,Reac
C
C --- Nucleus density and thickness
C
      REAL*8           RhoA,RhoB,RhoC,RhoR,RhoRmax,RhoH
      COMMON /NPRHO/   RhoA,RhoB,RhoC,RhoR,RhoRmax,RhoH

      INTEGER          RHORNUM,RHOBNUM
      PARAMETER       (RHORNUM = 41)        !! 2*n+1, n >= 1
      PARAMETER       (RHOBNUM = 7)         !! 2*n+1, n >= 1

      REAL*8           RhoRdat(RHORNUM)
      REAL*8           RhoBdat(RHOBNUM)
      REAL*8           RhoMdat(RHORNUM,RHORNUM,RHOBNUM)
      REAL*8           RhoTdat(RHORNUM,RHORNUM,RHOBNUM)
      COMMON /NPRHORB/ RhoRdat,RhoBdat
      COMMON /NPRHOM/  RhoMdat
      COMMON /NPRHOT/  RhoTdat

      REAL*8           RhoMdatB(RHORNUM,RHORNUM)
      REAL*8           RhoTdatB(RHORNUM,RHORNUM)
      COMMON /NPRHOMB/ RhoMdatB
      COMMON /NPRHOTB/ RhoTdatB
C
C --- Kinematics
C
      REAL*8           KinE ,KinEM          !! E  , E_min (neutrino)
      REAL*8           KinML,KinMP ,KinMN   !! m_l, m_pi, m_N
      REAL*8           KinSN,KinSNM,KinSNX  !! s_N, s_N_min, s_N_max
      REAL*8           KinMu,KinMuM,KinMuX  !!  mu,  mu_min,  mu_max
      REAL*8           KinNu,KinNuM,KinNuX  !!  nu,  nu_min,  nu_max
      REAL*8           KinQ2,KinQ2M,KinQ2X  !!  Q2,  Q2_min,  Q2_max
      REAL*8           KinX ,KinXM ,KinXX   !!   x,   x_min,   x_max
      REAL*8           KinY ,KinYM ,KinYX   !!   y,   y_min,   y_max
      REAL*8           KinT ,KinTM ,KinTX   !!   t,   t_min,   t_max
      REAL*8           KinTC,KinPC ,KinPNC  !! t_c, p_c, |p_N^c|
      REAL*8           KinKT,KinKL ,KinKLC  !! k_T, k_L,  k_L^c
      COMMON /NPKIN/   KinE ,KinEM
     .                ,KinML,KinMP ,KinMN
     .                ,KinSN,KinSNM,KinSNX
     .                ,KinMu,KinMuM,KinMuX
     .                ,KinNu,KinNuM,KinNuX
     .                ,KinQ2,KinQ2M,KinQ2X
     .                ,KinX ,KinXM ,KinXX
     .                ,KinY ,KinYM ,KinYX
     .                ,KinT ,KinTM ,KinTX
     .                ,KinTC,KinPC ,KinPNC
     .                ,KinKT,KinKL ,KinKLC
C
C --- Cross sections
C
      REAL*8           SigUnit              !! unit value (in GeV)
      REAL*8           SigFact              !! factor (CC,NC in units)
      INTEGER          SigType              !! type (0 - coh., 1 - inc.)
      REAL*8           SigBpN               !! slope BpN
      REAL*8           SigTot               !! sig_tot
      COMPLEX*16       SigTotM              !! sig_tot modified
      REAL*8           SigEl                !! sig_el
      COMPLEX*16       SigIn                !! sig_in
      REAL*8           SigEps               !! accuracy
      COMMON /NPSIG/   SigTotM,SigIn
     .                ,SigTot ,SigEl
     .                ,SigBpN
     .                ,SigEps
     .                ,SigUnit
     .                ,SigFact
     .                ,SigType
