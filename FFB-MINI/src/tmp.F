      SUBROUTINE VEL3D1
     * (MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *  IFORM,BTDCOE,DT,
     *  ITIME,DEVLP1,ACCELX,ACCELY,ACCELZ,
     *  NMAX,EPS,RESU,RESV,RESW,NITRU,NITRV,NITRW,
     *  ME,N,N1,N2,NE,NP,NEX,NODE,
     *  U,V,W,VISC,UE,VE,WE,
     *  MELM,E,EX,EY,EZ,EXX,EYY,EZZ,EXY,EXZ,EYZ,
     *  SN,
     *  NUMIP,
     *  A,NPP,NCRS,IPCRS,APCRS,LTAB,
     *  NPINLT,LPINLT,UINLT,VINLT,WINLT,
     *  NPWALL,LPWALL,UWALL,VWALL,WWALL,
     *  NPSYMT,NPSYM2,LPSYMT,XPSYMT,YPSYMT,ZPSYMT,
     *  IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *  LPFIX,LFIX3D,
     *  UG,VG,WG,WRK01,WRK02,WRK03,WRK04,A0,AR,RHSU,RHSV,RHSW,
     *  RX,RY,
     *  JUNROL,NPPMAX,NCRS2,TS,TACRS,ITPCRS,
     *  IUT0,IERR)
C
#include "timing.h"
      IMPLICIT NONE
C
      INTEGER*4 MCOLOR,MCPART
      INTEGER*4 NCOLOR(4),NCPART(MCOLOR,4),LLOOP(MCPART,MCOLOR,4)
C
      INTEGER*4 IFORM,ITIME,IPART,
     *          NMAX,NUMIP,
     *          ME,N,NEX,N1,N2,NE,NP,MELM,
     *          NODE,
     *          LPFIX,NCRS,
     *          NPWALL,LPWALL,
     *          NPINLT,LPINLT,
     *          NPSYMT,NPSYM2,LPSYMT,
     *          NDOM,LDOM,MBPDOM,NBPDOM,IPSLF,IPSND,
     *          LFIX3D,
     *          NPP,IPCRS,LTAB,IUT0,IERR
      REAL*4    BTDCOE(4),DT,EPS,RESU,RESV,RESW,
     *          DEVLP1,ACCELX,ACCELY,ACCELZ,
     *          U,V,W,VISC,
     *          SN,E,EX,EY,EZ,EXX,EYY,EZZ,EXY,EXZ,EYZ,
     *          UINLT,VINLT,WINLT,
     *          UWALL,VWALL,WWALL,XPSYMT,YPSYMT,ZPSYMT,
     *          RX,RY,
     *          A,UG,VG,WG,UE,VE,WE,
     *          WRK01,WRK02,WRK03,WRK04,A0,AR,
     *          RHSU,RHSV,RHSW,APCRS
C
      DIMENSION NEX(12),
     *          U(NP),V(NP),W(NP),UE(NE),VE(NE),WE(NE),
     *          NODE(N1,NE),
     *          VISC(NE),
     *          SN(N1,NE),E(MELM),
     *          EX (MELM),EY (MELM),EZ (MELM),
     *          EXX(MELM),EYY(MELM),EZZ(MELM),
     *          EXY(MELM),EXZ(MELM),EYZ(MELM) 
C         
      DIMENSION NUMIP(NP),LPFIX(NP),LPINLT(NPINLT),
     *          LPWALL(NPWALL),LPSYMT(NPSYM2),
     *          UINLT(NPINLT),VINLT(NPINLT),WINLT (NPINLT),
     *          UWALL(NPWALL),VWALL(NPWALL),WWALL(NPWALL),
     *          XPSYMT(NPSYM2),YPSYMT(NPSYM2),ZPSYMT(NPSYM2),
     *          LDOM(NDOM),NBPDOM(NDOM),
     *          IPSLF(MBPDOM,NDOM),IPSND(MBPDOM,NDOM) 
C
      DIMENSION RX(0:N,ME),RY(0:N,ME),LFIX3D(NP),
     *          A(N1,N2,NE),UG(*),VG(*),WG(*),
     *          WRK01(*),WRK02(*),WRK03(*),WRK04(*),A0(NP),AR(NP),
     *          RHSU(NP),RHSV(NP),RHSW(NP)
C
      DIMENSION NPP(NP),IPCRS(NCRS),APCRS(NCRS),LTAB(N1,N2,NE)
C
      REAL*4  DIJ(8,8)
      DATA DIJ / 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &           0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &           0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &           0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
     &           0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
     &           0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     &           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
     &           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 /
C
      INTEGER*4 NEHEX,NHEX,NSKIP4,
     *          IP,IE,I,J,IITRE,
     *          IES,IEE,
     *          IP1,IP2,IP3,IP4,IP5,IP6,IP7,IP8,IBP,
     *          IELM0,
     *          IELM1A,IELM2A,IELM3A,IELM4A,
     *          IELM5A,IELM6A,IELM7A,IELM8A,
     *          IELM1B,IELM2B,IELM3B,IELM4B,
     *          IELM5B,IELM6B,IELM7B,IELM8B,
     *          MAXBUF,NPFIX,IDUM,NITRE,
     *          NITRU,NITRV,NITRW,IRESU,IRESV,IRESW,
     *          IERR1,IERR2,IERR3,ICOLOR,ICPART
      REAL*4    DTH,CEB1,CEB2,UU,VV,WW,AT,
     *          AT1,AT2,AT3,AT4,AT5,AT6,AT7,AT8,
     *          AC1,AC2,AC3,AC4,AC5,AC6,AC7,AC8,
     *          AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8,
     *          CRHS1,CRHS2,CRHS3,CRHS4,CRHS5,CRHS6,CRHS7,CRHS8,
     *          ABTD1,ABTD2,ABTD3,ABTD4,ABTD5,ABTD6,ABTD7,ABTD8,
     *          FX0,FY0,FZ0,
     *          COF
C
C     [FULL UNROOL]
      INTEGER*4 JUNROL
      INTEGER*4 NPPMAX,NCRS2,ITPCRS(NCRS2)
      REAL*4    TS(0:NP),TACRS(NCRS2)
C      
C     [IN:MID NODE COLORING]
C
      CHARACTER*60 ERMSGC
     & /' ## SUBROUTINE VEL3D1: FATAL      ERROR REPORT   ; RETURNED' /
C
      INTEGER*4 IBCGS
      DATA IBCGS  / 0 /
      REAL*4 WRK05(NP)
C
C
C      CALCULATE VELOCITY PREDICTOR
C         ( 3-D CALCULATION : SINGLE WORD & MULTI ELEMENT VERSION )
C                                              2009.01.30 Y.YAMADE
C
C ************ COMPUTATIONAL COST EXCEPT FOR MATRIX SOLVER *******
C =============================HEX======================================
C          OPERATION COUNTS:   2656 FLOP /ELEMENT
C          DATA LOADINGS   :   1417 WORDS/ELEMENT
C                           (    41 WORDS CONTIGUOUSLY,
C                               978 WORDS BY STRIDE, AND
C                               398 WORDS BY LIST )
C
C     ARGUMENT LISTINGS
C       (1) INPUT
C          NLOOP       ;NUMBER OF LOOPS
C          LLOOP       ;POINTER FOR SPLITTED ELEMENT LIST
C          IFORM       ; SPECIFIES MOMENTUM EQUATIONS METHOD
C
C          TIMER       ; PRESENT TIME OF OVERSET CONDITIONS DATA
C           NOTES ; 'TIMER' WILL BE REFERED TO FOR INTER-FLAME OVERSET.
C          DT          ; TIME INCTREMENT
C
C          NMAX        ; NUMBER OF MATRIX SOLVER ITERATIONS
C          EPS         ; MAXIMUM ALLOWABLE ERROR
C
C          ITIME       ; CUREENT TIME STEP
C
C          DEVLP1      ; DEVELOPMENT FUNCTION FOR INLET VELOCITIES
C          DEVLP2      ; DEVELOPMENT FUNCTION FOR ALL THE OTHER VALUES
C          ACCELX      ; X-DIR. ACCELERATION TERMS ADDED TO ALL FRAMES
C          ACCELY      ; Y-DIR. ACCELERATION TERMS ADDED TO ALL FRAMES
C          ACCELZ      ; Z-DIR. ACCELERATION TERMS ADDED TO ALL FRAMES
C
C          ME          ; MAX. NUMBER OF TOTAL ELEMENTS
C          N           ; NUMBER OF NODES ASSIGNED TO ONE ELEMENT
C          NE          ; NUMBER OF TOTAL ELEMENTS
C          NP          ; NUMBER OF TOTAL    NODES
C          NEX(I)      ; INCLUDES NUMBER OF ELEMENTS AND NUMBER OF LOCAL NODES
C                        AS FOLOOWS
C          NEX(1)      ; NUMBER OF TET.    ELEMENTS
C          NEX(2)      ; NUMBER OF PYRAMID ELEMENTS
C          NEX(3)      ; NUMBER OF WEGDE   ELEMENTS
C          NEX(4)      ; NUMBER OF HEX.    ELEMENTS
C          NEX(5)      ; NUMBER OF LOCAL NODES IN A TET.    ELEMENT (=4)
C          NEX(6)      ; NUMBER OF LOCAL NODES IN A PYRAMID ELEMENT (=5)
C          NEX(7)      ; NUMBER OF LOCAL NODES IN A WEGDE   ELEMENT (=6)
C          NEX(8)      ; NUMBER OF LOCAL NODES IN A HEX.    ELEMENT (=8)
C
C          NODE  (I,IE); NODE NO. TABLE BASED ON ELEMENT
C          VISC    (IE); ELEMENT VISCOSITY
C
C          SN    (I,IE); INTEGRATED ELEMENT VECTOR OF N
C          E   (I,J,IE); INTEGRATED ELEMENT MATRIX OF N*NT
C          EX  (I,J,IE); INTEGRATED ELEMENT MATRIX OF N*NXT
C          EY  (I,J,IE); INTEGRATED ELEMENT MATRIX OF N*NYT
C          EZ  (I,J,IE); INTEGRATED ELEMENT MATRIX OF N*NZT
C          EXX (I,J,IE); INTEGRATED ELEMENT MATRIX OF NX*NXT
C          EYY (I,J,IE); INTEGRATED ELEMENT MATRIX OF NY*NYT
C          EZZ (I,J,IE); INTEGRATED ELEMENT MATRIX OF NZ*NZT
C          EXY (I,J,IE); INTEGRATED ELEMENT MATRIX OF NX*NYT
C          EXZ (I,J,IE); INTEGRATED ELEMENT MATRIX OF NX*NZT
C          EYZ (I,J,IE); INTEGRATED ELEMENT MATRIX OF NY*NZT
C
C          NUMIP   (IP); NUMBER OF NEIGHBORING DOMAINS THAT NODE
C                        'IP' BELONG TO
C
C        A. INLET BOUNDARY
C          NPINLT      ; NUMBER OF INLET BOUNDARY NODES
C          LPINLT (IBP); INLET BOUNDARY NODES
C          UINLT  (IBP); INLET BOUNDARY U-VELOCITIES
C          VINLT  (IBP); INLET BOUNDARY V-VELOCITIES
C          WINLT  (IBP); INLET BOUNDARY W-VELOCITIES
C
C        B. WALL BOUNDARY
C          NPWALL      ; NUMBER OF WALL BOUNDARY NODES
C          LPWALL (IBP); WALL BOUNDARY NODES
C          UWALL  (IBP); WALL BOUNDARY U-VELOCITIES
C          VWALL  (IBP); WALL BOUNDARY V-VELOCITIES
C          WWALL  (IBP); WALL BOUNDARY W-VELOCITIES
C
C        C. SYMMETRIC BOUNDARY
C          NPSYMT      ; NUMBER OF SYMMETRIC BOUNDARY NODES
C          LPSYMT (IBP); SYMMETRIC BOUNDARY NODES
C          XPSYMT (IBP); X-DIR COMPONENT OF SYMMETRIC NODE NORMAL VECTOR
C          YPSYMT (IBP); Y-DIR COMPONENT OF SYMMETRIC NODE NORMAL VECTOR
C          ZPSYMT (IBP); Z-DIR COMPONENT OF SYMMETRIC NODE NORMAL VECTOR
C
C        D. INTER-CONNECT BOUNDARY
C          IPART       ; SUB-DOMAIN NUMBER THAT THIS TASK SHOULD TAKE/IS
C                       TAKING CARE OF. IPART BEING SET ZERO MEANS THAT
C                       THE PROGRAM SHOULD RUN/IS RUNNING IN SERIAL 
C                       MODE.
C
C          LDOM  (IDOM); NEIBERING SUB-DOMAIN NUMBER
C          NBPDOM(IDOM); NUMBER OF INTER-CONNECT BOUNDARY NODES
C                       SHARING WITH THE IDOM'TH NEIBERING SUB-DOMAIN,
C                       LDOM(IDOM)
C          NDOM        ; NUMBER OF THE NERIBERING SUB-DOMAINS
C          IPSLF (IBP,IDOM); INTER-CONNECT BOUNDARY NODE NUMBER IN THE
C                           CALLING TASK'S SUB-DOMAIN, FOR THE IDOM'TH
C                           NEIBERING SUB-DOMAIN, LDOM(IDOM)
C          IPSND (IBP,IDOM); INTER-CONNECT BOUNDARY NODE NUMBER IN THE
C                           SUB-DOMAIN THAT IS RECEIVING THE CALLING
C                           TASK'S RESIDUALS.
C          MBPDOM      ; THE MAXIMUM NUMBER OF THE INTER-CONNECT 
C                       BOUNDARY NODES FOR ONE NEIBERING SUB-DOMAIN
C
C
C          NPP      (IP); NUMBER OF ADJACENT NODES    TO NODE    IP
C          NCRS    ; NUMBER OF NONZERO ELEMENTS IN MATRIX OF CRS FORMAT
C          IPCRS  (ICRS); NODE NO. TABLE BASED ON CRS FORMAT
C          LTAB(J1,J2,IE); CRS INDEX TABLE FOR NODE-BASE MATRIX
C                          COEFFICIENT
C
C          IUT0        ; FILE NUMBER TO REPORT ERROR OCCURENCE
C
C       (2) OUTPUT
C          RESU        ;RESIDUAL OF U-EQUATION
C          RESV        ;RESIDUAL OF V-EQUATION 
C          RESW        ;RESIDUAL OF W-EQUATION
C
C          NRNU        ;NUMBER OF U-EQUATION ITERATIONS 
C          NRNV        ;NUMBER OF V-EQUATION ITERATIONS  
C          NRNW        ;NUMBER OF W-EQUATION ITERATIONS 
C
C          IERR        ; RETURN CODE TO REPORT ERROR OCCURENCE
C                   0 --- NORMAL TERMINATION
C                   1 --- A FATAL ERROR HAS OCCURED
C
C       (3) INPUT-OUTPUT
C          U       (IP); X-DIR. VELOCITY COMPONENT
C          V       (IP); Y-DIR. VELOCITY COMPONENT
C          W       (IP); Z-DIR. VELOCITY COMPONENT
C
C       (4) WORK
C          LPFIX   (IB); VELOCITY FIX BOUNDARY CONDITION NODES WHICH ARE
C                        COMPOSED OF INLET, WALL, SYMMETRIC AND OVERSET
C                        BOUNDARY CONDITIONS NODES
C          LWROK   (IP); WORK REGION TO MAKE 'LPFIX'
C
C          A   (I,J,IE); ELEMENT-WISE MATRIX COEFFICIENT WHICH INCLUDES 
C                        ALL THE ERMS AND WILL BE PASSED TO THE MATROX SOLVER
C          APCRS  (ICRS); NODE-BASE MATRIX COEFFICIENT WHICH INCLUDES
C                  ALL THE ERMS AND WILL BE PASSED TO THE MATRIX SOLVER
C
C          UG      (IE); WORK REGION FOR U-VELOCITY AT ELEMENTS
C          VG      (IE); WORK REGION FOR V-VELOCITY AT ELEMENTS
C          WG      (IE); WORK REGION FOR W-VELOCITY AT ELEMENTS
C
C          RHSU    (IP); WORK REGION FOR RHS IN U-EQUATION
C          RHSV    (IP); WORK REGION FOR RHS IN V-EQUATION
C          RHSW    (IP); WORK REGION FOR RHS IN W-EQUATION
      IERR=0
C
C
CVEL3D1 [ 1.] SET CONTROL PARAMETER
C
      MAXBUF = NE*(N+1)
      DTH=DT*0.5E0
C
CC    CRANCK-NICOLSON FOR CONVECTION TERM
      CEB1 = 0.0E0
      IF(IFORM.EQ.5)  CEB1 = 1.0E0
C
      IF(IFORM.EQ.1.OR.IFORM.EQ.2) THEN 
CC    CRANCK-NICOLSON FOR DIFF. TERM
         CEB2 = 0.0E0
      ELSE
CC    EULAR BACKWARD DIFF. TERM 
         CEB2 = 1.0E0
      END IF  
C
      IF(IFORM.EQ.2 .OR. IFORM.EQ.4 .OR. IFORM.EQ.5) THEN 
C
      ELSE
         BTDCOE(1) = 0.0E0
         BTDCOE(2) = 0.0E0
         BTDCOE(3) = 0.0E0
         BTDCOE(4) = 0.0E0
      END IF   
C
      NEHEX=NEX(4)
      NHEX=NEX(8)
      NSKIP4=NEX(12)
C
CVEL3D1 [ 2.] MAKE FIX BOUNDARY NODES LIST (LPFIX)
C
      DO 1000 IP = 1 , NP
          LFIX3D(IP) = 0
          LPFIX(IP) = 0
 1000 CONTINUE
C
*POPTION INDEP(LFIX3D)
C*$*ASSERT PERMUTATION ( LPINLT )
!ocl norecurrence(LFIX3D)      
      DO 1010 IBP = 1 , NPINLT
          LFIX3D(LPINLT(IBP))=1
          U(LPINLT(IBP)) = DEVLP1*UINLT(IBP)
          V(LPINLT(IBP)) = DEVLP1*VINLT(IBP)
          W(LPINLT(IBP)) = DEVLP1*WINLT(IBP)
 1010 CONTINUE
C
*POPTION INDEP(LFIX3D)
C*$*ASSERT PERMUTATION ( LPWALL )
!ocl norecurrence(LFIX3D)      
      DO 1020 IBP = 1 , NPWALL
          LFIX3D(LPWALL(IBP))=1
 1020 CONTINUE
C
CCCCCC*POPTION INDEP(LFIX3D)
CCCCCCC*$*ASSERT PERMUTATION ( LPSYMT )
CCCCCC      DO 130 IBP = 1 , NPSYMT
CCCCCC          LFIX3D(LPSYMT(IBP))=1
CCCCCC  130 CONTINUE
C
*POPTION INDEP(LPFIX)
C*$*ASSERT PERMUTATION ( LFIX3D )
      NPFIX=0
!ocl norecurrence(LFIX3D)      
      DO 1070 IP = 1 , NP
          IF(LFIX3D(IP).EQ.0) GO TO 1070
          NPFIX=NPFIX+1
          LPFIX(NPFIX) = IP
 1070 CONTINUE
C
C INITIALIZE
      DO 1080 IP = 1 , NP
         AR   (IP)=0.0E0
         WRK04(IP)=0.0E0
 1080 CONTINUE
      IF(ITIME.EQ.0)THEN
          GO TO 3000
      ENDIF
C
CVEL3D1 [ 3.] CAL. TIME, ADV. AND VIS. TEAM IN LHS AND RHS
C
C*    call maprof_time_start(TM_CALUEL)
      CALL CALUEL(N1,NE,NP,NEX,NODE,U,V,W,UE,VE,WE,IERR)
C*    call maprof_time_stop(TM_CALUEL)
C
      DO 1100 IP=1,NP
          RHSU(IP)=0.0E0
          RHSV(IP)=0.0E0
          RHSW(IP)=0.0E0
          A0  (IP)=0.0E0
 1100 CONTINUE
C
      DO 1220 IE=1,NE
      DO 1210 I=1,N
      DO 1200 J=1,N
          A(J,I,IE)=0.0E0
 1200 CONTINUE
 1210 CONTINUE
 1220 CONTINUE
C
C          OPERATION COUNTS:  FLOP /ELEMENT
C          DATA LOADINGS   : 1792 WORDS/ELEMENT
C                           (  32 WORDS CONTIGUOUSLY,
C                             512 WORDS BY STRIDE, AND
C                             216 WORDS BY LIST )
C
      call maprof_time_start(TM_VEL3D1_OP1)
      DO 1343 ICOLOR=1,NCOLOR(4)
!ocl norecurrence(A,A0,AR,WRK04,RHSU,RHSV,RHSW)
      DO 1342 ICPART=1,NCPART(ICOLOR,4)
          IES=LLOOP(ICPART  ,ICOLOR,4)
          IEE=LLOOP(ICPART+1,ICOLOR,4)-1
!ocl nosimd
!ocl noswp
          DO 1341 IE=IES,IEE
          DO 1340 I=1,NHEX
              IP =NODE(I,IE)
              IP1=NODE(1,IE)
              IP2=NODE(2,IE)
              IP3=NODE(3,IE)
              IP4=NODE(4,IE)
              IP5=NODE(5,IE)
              IP6=NODE(6,IE)
              IP7=NODE(7,IE)
              IP8=NODE(8,IE)
C
              UU=UE(IE)
              VV=VE(IE)
              WW=WE(IE)
C
              IELM0 = 1+NSKIP4*(IE-1)
C
              IELM1A=IELM0+(I-1)+NHEX*0
              IELM2A=IELM0+(I-1)+NHEX*1
              IELM3A=IELM0+(I-1)+NHEX*2
              IELM4A=IELM0+(I-1)+NHEX*3
              IELM5A=IELM0+(I-1)+NHEX*4
              IELM6A=IELM0+(I-1)+NHEX*5
              IELM7A=IELM0+(I-1)+NHEX*6
              IELM8A=IELM0+(I-1)+NHEX*7
              IELM1B=IELM0+0+NHEX*(I-1)
              IELM2B=IELM0+1+NHEX*(I-1)
              IELM3B=IELM0+2+NHEX*(I-1)
              IELM4B=IELM0+3+NHEX*(I-1)
              IELM5B=IELM0+4+NHEX*(I-1)
              IELM6B=IELM0+5+NHEX*(I-1)
              IELM7B=IELM0+6+NHEX*(I-1)
              IELM8B=IELM0+7+NHEX*(I-1)
C
              ABTD1 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM1A)+VV*EXY(IELM1A)+WW*EXZ(IELM1A))
     &           +VV*(UU*EXY(IELM1B)+VV*EYY(IELM1A)+WW*EYZ(IELM1A))
     &           +WW*(UU*EXZ(IELM1B)+VV*EYZ(IELM1B)+WW*EZZ(IELM1A)))*DTH
              ABTD2 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM2A)+VV*EXY(IELM2A)+WW*EXZ(IELM2A))
     &           +VV*(UU*EXY(IELM2B)+VV*EYY(IELM2A)+WW*EYZ(IELM2A))
     &           +WW*(UU*EXZ(IELM2B)+VV*EYZ(IELM2B)+WW*EZZ(IELM2A)))*DTH
              ABTD3 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM3A)+VV*EXY(IELM3A)+WW*EXZ(IELM3A))
     &           +VV*(UU*EXY(IELM3B)+VV*EYY(IELM3A)+WW*EYZ(IELM3A))
     &           +WW*(UU*EXZ(IELM3B)+VV*EYZ(IELM3B)+WW*EZZ(IELM3A)))*DTH
              ABTD4 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM4A)+VV*EXY(IELM4A)+WW*EXZ(IELM4A))
     &           +VV*(UU*EXY(IELM4B)+VV*EYY(IELM4A)+WW*EYZ(IELM4A))
     &           +WW*(UU*EXZ(IELM4B)+VV*EYZ(IELM4B)+WW*EZZ(IELM4A)))*DTH
              ABTD5 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM5A)+VV*EXY(IELM5A)+WW*EXZ(IELM5A))
     &           +VV*(UU*EXY(IELM5B)+VV*EYY(IELM5A)+WW*EYZ(IELM5A))
     &           +WW*(UU*EXZ(IELM5B)+VV*EYZ(IELM5B)+WW*EZZ(IELM5A)))*DTH
              ABTD6 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM6A)+VV*EXY(IELM6A)+WW*EXZ(IELM6A))
     &           +VV*(UU*EXY(IELM6B)+VV*EYY(IELM6A)+WW*EYZ(IELM6A))
     &           +WW*(UU*EXZ(IELM6B)+VV*EYZ(IELM6B)+WW*EZZ(IELM6A)))*DTH
              ABTD7 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM7A)+VV*EXY(IELM7A)+WW*EXZ(IELM7A))
     &           +VV*(UU*EXY(IELM7B)+VV*EYY(IELM7A)+WW*EYZ(IELM7A))
     &           +WW*(UU*EXZ(IELM7B)+VV*EYZ(IELM7B)+WW*EZZ(IELM7A)))*DTH
              ABTD8 = 0.5E0*DT*BTDCOE(4)
     &          *(UU*(UU*EXX(IELM8A)+VV*EXY(IELM8A)+WW*EXZ(IELM8A))
     &           +VV*(UU*EXY(IELM8B)+VV*EYY(IELM8A)+WW*EYZ(IELM8A))
     &           +WW*(UU*EXZ(IELM8B)+VV*EYZ(IELM8B)+WW*EZZ(IELM8A)))*DTH
C
              AT =E(IELM1A)+E(IELM2A)+E(IELM3A)+E(IELM4A) 
     *           +E(IELM5A)+E(IELM6A)+E(IELM7A)+E(IELM8A)
              AT1=DIJ(1,I)*AT
              AT2=DIJ(2,I)*AT
              AT3=DIJ(3,I)*AT
              AT4=DIJ(4,I)*AT
              AT5=DIJ(5,I)*AT
              AT6=DIJ(6,I)*AT
              AT7=DIJ(7,I)*AT
              AT8=DIJ(8,I)*AT
C
              AC1=(UU*EX(IELM1A)+VV*EY(IELM1A)+WW*EZ(IELM1A))*DTH
              AC2=(UU*EX(IELM2A)+VV*EY(IELM2A)+WW*EZ(IELM2A))*DTH
              AC3=(UU*EX(IELM3A)+VV*EY(IELM3A)+WW*EZ(IELM3A))*DTH
              AC4=(UU*EX(IELM4A)+VV*EY(IELM4A)+WW*EZ(IELM4A))*DTH
              AC5=(UU*EX(IELM5A)+VV*EY(IELM5A)+WW*EZ(IELM5A))*DTH
              AC6=(UU*EX(IELM6A)+VV*EY(IELM6A)+WW*EZ(IELM6A))*DTH
              AC7=(UU*EX(IELM7A)+VV*EY(IELM7A)+WW*EZ(IELM7A))*DTH
              AC8=(UU*EX(IELM8A)+VV*EY(IELM8A)+WW*EZ(IELM8A))*DTH
C
              AD1=(EXX(IELM1A)+EYY(IELM1A)+EZZ(IELM1A))*VISC(IE)*DTH
              AD2=(EXX(IELM2A)+EYY(IELM2A)+EZZ(IELM2A))*VISC(IE)*DTH
              AD3=(EXX(IELM3A)+EYY(IELM3A)+EZZ(IELM3A))*VISC(IE)*DTH
              AD4=(EXX(IELM4A)+EYY(IELM4A)+EZZ(IELM4A))*VISC(IE)*DTH
              AD5=(EXX(IELM5A)+EYY(IELM5A)+EZZ(IELM5A))*VISC(IE)*DTH
              AD6=(EXX(IELM6A)+EYY(IELM6A)+EZZ(IELM6A))*VISC(IE)*DTH
              AD7=(EXX(IELM7A)+EYY(IELM7A)+EZZ(IELM7A))*VISC(IE)*DTH
              AD8=(EXX(IELM8A)+EYY(IELM8A)+EZZ(IELM8A))*VISC(IE)*DTH
C
              A(1,I,IE)=AT1+AC1*(1.0E0+CEB1)+AD1*(1.0E0+CEB2)+ABTD1
              A(2,I,IE)=AT2+AC2*(1.0E0+CEB1)+AD2*(1.0E0+CEB2)+ABTD2
              A(3,I,IE)=AT3+AC3*(1.0E0+CEB1)+AD3*(1.0E0+CEB2)+ABTD3
              A(4,I,IE)=AT4+AC4*(1.0E0+CEB1)+AD4*(1.0E0+CEB2)+ABTD4
              A(5,I,IE)=AT5+AC5*(1.0E0+CEB1)+AD5*(1.0E0+CEB2)+ABTD5
              A(6,I,IE)=AT6+AC6*(1.0E0+CEB1)+AD6*(1.0E0+CEB2)+ABTD6
              A(7,I,IE)=AT7+AC7*(1.0E0+CEB1)+AD7*(1.0E0+CEB2)+ABTD7
              A(8,I,IE)=AT8+AC8*(1.0E0+CEB1)+AD8*(1.0E0+CEB2)+ABTD8
              A0(IP)=A0(IP)+A(I,I,IE)
C
              WRK04(IP)=WRK04(IP)
     *                 +A(1,I,IE)+A(2,I,IE)+A(3,I,IE)+A(4,I,IE)
     *                 +A(5,I,IE)+A(6,I,IE)+A(7,I,IE)+A(8,I,IE)
C
              CRHS1=AT1-AC1*(1.0E0-CEB1)-AD1*(1.0E0-CEB2)-ABTD1
              CRHS2=AT2-AC2*(1.0E0-CEB1)-AD2*(1.0E0-CEB2)-ABTD2 
              CRHS3=AT3-AC3*(1.0E0-CEB1)-AD3*(1.0E0-CEB2)-ABTD3
              CRHS4=AT4-AC4*(1.0E0-CEB1)-AD4*(1.0E0-CEB2)-ABTD4
              CRHS5=AT5-AC5*(1.0E0-CEB1)-AD5*(1.0E0-CEB2)-ABTD5
              CRHS6=AT6-AC6*(1.0E0-CEB1)-AD6*(1.0E0-CEB2)-ABTD6
              CRHS7=AT7-AC7*(1.0E0-CEB1)-AD7*(1.0E0-CEB2)-ABTD7
              CRHS8=AT8-AC8*(1.0E0-CEB1)-AD8*(1.0E0-CEB2)-ABTD8
C
              RHSU(IP)=RHSU(IP)
     *              +CRHS1*U(IP1)+CRHS2*U(IP2)+CRHS3*U(IP3)+CRHS4*U(IP4)
     *              +CRHS5*U(IP5)+CRHS6*U(IP6)+CRHS7*U(IP7)+CRHS8*U(IP8)
              RHSV(IP)=RHSV(IP)
     *              +CRHS1*V(IP1)+CRHS2*V(IP2)+CRHS3*V(IP3)+CRHS4*V(IP4)
     *              +CRHS5*V(IP5)+CRHS6*V(IP6)+CRHS7*V(IP7)+CRHS8*V(IP8)
              RHSW(IP)=RHSW(IP)
     *              +CRHS1*W(IP1)+CRHS2*W(IP2)+CRHS3*W(IP3)+CRHS4*W(IP4)
     *              +CRHS5*W(IP5)+CRHS6*W(IP6)+CRHS7*W(IP7)+CRHS8*W(IP8)
 1340     CONTINUE
 1341     CONTINUE
 1342 CONTINUE
 1343 CONTINUE
      call maprof_time_stop(TM_VEL3D1_OP1)
C
      DO IP=1,NP
          AR(IP)=A0(IP)
          IF (AR(IP).EQ.0.0E0) AR(IP)=1.0E0
      ENDDO
C 
C
C
CVEL3D1 [ 5.] CAL. ACCELERATION TERMS IN RHS
C
C
      FX0=ACCELX
      FY0=ACCELY
      FZ0=ACCELZ
!ocl norecurrence(WRK01,WRK02,WRK03)      
      DO 1500 IE = 1,NE
          WRK01(IE)=FX0*DT
          WRK02(IE)=FY0*DT
          WRK03(IE)=FZ0*DT
 1500 CONTINUE
C
C          OPERATION COUNTS:  48 FLOP /ELEMENT
C          DATA LOADINGS   :  43 WORDS/ELEMENT
C                           (  3 WORDS CONTIGUOUSLY,
C                             16 WORDS BY STRIDE, AND
C                             24 WORDS BY LIST )
C
C*    call maprof_time_start(TM_VEL3D1_OP2)
#ifdef USE_DETAIL
      call start_collection('vel3d1_1642')
#endif      
      DO 1642 ICOLOR=1,NCOLOR(4)
!ocl norecurrence(RHSU,RHSV,RHSW)
      DO 1641 ICPART=1,NCPART(ICOLOR,4)
          IES=LLOOP(ICPART  ,ICOLOR,4)
          IEE=LLOOP(ICPART+1,ICOLOR,4)-1
      DO 1640 IE=IES,IEE
          IP1=NODE(1,IE)
          IP2=NODE(2,IE)
          IP3=NODE(3,IE)
          IP4=NODE(4,IE)
          IP5=NODE(5,IE)
          IP6=NODE(6,IE)
          IP7=NODE(7,IE)
          IP8=NODE(8,IE)
C
          RHSU(IP1)=RHSU(IP1)+WRK01(IE)*SN(1,IE)
          RHSU(IP2)=RHSU(IP2)+WRK01(IE)*SN(2,IE)
          RHSU(IP3)=RHSU(IP3)+WRK01(IE)*SN(3,IE)
          RHSU(IP4)=RHSU(IP4)+WRK01(IE)*SN(4,IE)
          RHSU(IP5)=RHSU(IP5)+WRK01(IE)*SN(5,IE)
          RHSU(IP6)=RHSU(IP6)+WRK01(IE)*SN(6,IE)
          RHSU(IP7)=RHSU(IP7)+WRK01(IE)*SN(7,IE)
          RHSU(IP8)=RHSU(IP8)+WRK01(IE)*SN(8,IE)
C
          RHSV(IP1)=RHSV(IP1)+WRK02(IE)*SN(1,IE)
          RHSV(IP2)=RHSV(IP2)+WRK02(IE)*SN(2,IE)
          RHSV(IP3)=RHSV(IP3)+WRK02(IE)*SN(3,IE)
          RHSV(IP4)=RHSV(IP4)+WRK02(IE)*SN(4,IE)
          RHSV(IP5)=RHSV(IP5)+WRK02(IE)*SN(5,IE)
          RHSV(IP6)=RHSV(IP6)+WRK02(IE)*SN(6,IE)
          RHSV(IP7)=RHSV(IP7)+WRK02(IE)*SN(7,IE)
          RHSV(IP8)=RHSV(IP8)+WRK02(IE)*SN(8,IE)
C
          RHSW(IP1)=RHSW(IP1)+WRK03(IE)*SN(1,IE)
          RHSW(IP2)=RHSW(IP2)+WRK03(IE)*SN(2,IE)
          RHSW(IP3)=RHSW(IP3)+WRK03(IE)*SN(3,IE)
          RHSW(IP4)=RHSW(IP4)+WRK03(IE)*SN(4,IE)
          RHSW(IP5)=RHSW(IP5)+WRK03(IE)*SN(5,IE)
          RHSW(IP6)=RHSW(IP6)+WRK03(IE)*SN(6,IE)
          RHSW(IP7)=RHSW(IP7)+WRK03(IE)*SN(7,IE)
          RHSW(IP8)=RHSW(IP8)+WRK03(IE)*SN(8,IE)
 1640 CONTINUE
 1641 CONTINUE
 1642 CONTINUE
#ifdef USE_DETAIL
      call stop_collection('vel3d1_1642')
#endif
C*    call maprof_time_stop(TM_VEL3D1_OP2)
C
C
CVEL3D1 [ 6.] MAKE CRS MATIRX
      call maprof_time_start(TM_E2PMATR)
      CALL E2PMATR(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *             N2,N1,NE,NCRS,A,APCRS,LTAB,IERR)
      call maprof_time_stop(TM_E2PMATR)
      IF (IERR.NE.0) THEN
          WRITE(IUT0,*)
          WRITE(IUT0,*) ERMSGC
          RETURN
      ENDIF
C
C 
C SUPERIMPOSE DIAGONAL TERM (AR)
C
      IDUM=1
C*    call maprof_time_start(TM_VEL3D1_COM)
      CALL DDCOMX(IPART,IDUM,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *            AR,AR,AR,NP,IUT0,IERR,RX,RY,MAXBUF)
C*    call maprof_time_stop(TM_VEL3D1_COM)
      IF(IERR.NE.0) THEN
          WRITE(IUT0,*)
          WRITE(IUT0,*) ERMSGC
          RETURN
      ENDIF
C
CVEL3D1 [10.] DIAGONAL SCALING OF CRS MATRIX AND RHS
C
      call maprof_time_start(TM_DGNSCL)
      CALL DGNSCL(APCRS,AR,NP,NCRS,NPP)
      call maprof_time_stop(TM_DGNSCL)
C
      DO 4100 IP=1, NP
          RHSU(IP) = RHSU(IP)/AR(IP)
          RHSV(IP) = RHSV(IP)/AR(IP)
          RHSW(IP) = RHSW(IP)/AR(IP)
 4100 CONTINUE
C
      IDUM = 3
C*    call maprof_time_start(TM_VEL3D1_COM)
      CALL DDCOMX(IPART,IDUM,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *            RHSU,RHSV,RHSW,NP,IUT0,IERR,RX,RY,MAXBUF)
C*    call maprof_time_stop(TM_VEL3D1_COM)
      IF(IERR.NE.0) THEN
          WRITE(IUT0,*)
          WRITE(IUT0,*) ERMSGC
          RETURN
      ENDIF
C
CVEL3D1 [11.] CLEAR CRS MATRIX FOR DIRICHLET B.C.
C
      call maprof_time_start(TM_CLRCRS)
      CALL CLRCRS(APCRS,NP,NCRS,IPCRS,NPP,
     *            LFIX3D,NUMIP,UG)
      call maprof_time_stop(TM_CLRCRS)
C
CVEL3D1 [12.] ITERATIONS FOR SOLVING MOMENTUM EQUATIONS
C
 3000 CONTINUE
C
      IF(ITIME.EQ.0) THEN
         NITRE=0
      ELSE
         NITRE=1
      END IF
C
      NITRU=0
      NITRV=0
      NITRW=0
C
      IRESU=0
      IRESV=0
      IRESW=0
C
      DO 5000 IITRE=0,NITRE
          IF(IITRE.EQ.0) GOTO 5100
C
C          OPERATION COUNTS:   4  FLOP /NODE/NPPAVE/ITERATION
C          DATA LOADINGS   :   6  WORDS/NODE/NPPAVE/ITERATION
C                           (  4  WORDS CONTIGUOUSLY,
C                              0  WORDS BY 4-WORD STRIDE, AND
C                              2  WORDS BY LIST )
C           NPPAVE: AVE. NUMBER OF ADJACENT NODES AT A NODE
C           --> 15 (TET), 20 (PRD,WED), 27(HEX.) 
C
C-----------------------------------------------------------------------
C
C
C
CVEL3D1 [12.1] CALL MATRIX SOLVER
C
      IERR1=1
      IERR2=1
      IERR3=1
C
      call maprof_time_start(TM_CRSCVA)
      IF(JUNROL.EQ.1) THEN
         CALL CRSCVA(NP,NPPMAX,NCRS,NCRS2,NPP,APCRS,TACRS)
      ENDIF
      call maprof_time_stop(TM_CRSCVA)
C
      call maprof_time_start(TM_BCGS3X)
      CALL BCGS3X(IBCGS,IPART,NMAX,EPS,ME,N,NE,NP,
     *            NPP,NCRS,IPCRS,APCRS,RHSU,U,NITRU,RESU,
     *            NDOM,MBPDOM,LDOM,NBPDOM,IPSLF,IPSND,NUMIP,
     *            UG,VG,WG,WRK01,WRK02,WRK03,WRK04,WRK05,RX,RY,
     *            IUT0,IERR1,
     *            JUNROL,NPPMAX,NCRS2,TS,TACRS,ITPCRS)
      IF(NITRU.LT.NMAX) IRESU=1
      CALL BCGS3X(IBCGS,IPART,NMAX,EPS,ME,N,NE,NP,
     *            NPP,NCRS,IPCRS,APCRS,RHSV,V,NITRV,RESV,
     *            NDOM,MBPDOM,LDOM,NBPDOM,IPSLF,IPSND,NUMIP,
     *            UG,VG,WG,WRK01,WRK02,WRK03,WRK04,WRK05,RX,RY,
     *            IUT0,IERR2,
     *            JUNROL,NPPMAX,NCRS2,TS,TACRS,ITPCRS)
      IF(NITRV.LT.NMAX) IRESV=1
      CALL BCGS3X(IBCGS,IPART,NMAX,EPS,ME,N,NE,NP,
     *            NPP,NCRS,IPCRS,APCRS,RHSW,W,NITRW,RESW,
     *            NDOM,MBPDOM,LDOM,NBPDOM,IPSLF,IPSND,NUMIP,
     *            UG,VG,WG,WRK01,WRK02,WRK03,WRK04,WRK05,RX,RY,
     *            IUT0,IERR3,
     *            JUNROL,NPPMAX,NCRS2,TS,TACRS,ITPCRS)
      IF(NITRW.LT.NMAX) IRESW=1
      call maprof_time_stop(TM_BCGS3X)
C
      IF(IERR1.NE.0 .OR. IERR2.NE.0 .OR. IERR3.NE.0) THEN
         WRITE(IUT0,*)
         WRITE(IUT0,*) ERMSGC
         RETURN
      ENDIF
C
       IF(IERR.NE.0) THEN
          WRITE(IUT0,*)
          WRITE(IUT0,*) ERMSGC
          RETURN
       ENDIF
C
CVEL3D1 [12.2] SET BOUNDARY CONDITIONS
C
 5100    CONTINUE
C
C
CCYY   A. INLET BOUNDARY CONDITIONS
C
*POPTION INDEP(U,V,W)
C*$*ASSERT PERMUTATION ( LPINLT )
!ocl norecurrence(U,V,W)
         DO 6010 IBP = 1 , NPINLT
            U(LPINLT(IBP)) = DEVLP1*UINLT(IBP)
            V(LPINLT(IBP)) = DEVLP1*VINLT(IBP)
            W(LPINLT(IBP)) = DEVLP1*WINLT(IBP)
 6010    CONTINUE
C
C
CCYY   B. WALL BOUNDARY CONDITIONS
C
C
*POPTION INDEP(U,V,W)
C*$*ASSERT PERMUTATION ( LPWALL )
!ocl norecurrence(U,V,W)
         DO 6020 IBP = 1 , NPWALL
            U(LPWALL(IBP)) = DEVLP1*UWALL(IBP)
            V(LPWALL(IBP)) = DEVLP1*VWALL(IBP)
            W(LPWALL(IBP)) = DEVLP1*WWALL(IBP)
 6020    CONTINUE
C
C
CCYY   C. SYMMETRIC BOUNDARY CONDITIONS
C
C
*POPTION INDEP(U,V,W)
C*$*ASSERT PERMUTATION ( LPSYMT )
!ocl norecurrence(U,V,W)
         DO 6040 IBP = 1 , NPSYMT
            COF = XPSYMT(IBP)*U(LPSYMT(IBP))
     *           +YPSYMT(IBP)*V(LPSYMT(IBP))
     *           +ZPSYMT(IBP)*W(LPSYMT(IBP))
            U(LPSYMT(IBP)) = U(LPSYMT(IBP))-COF*XPSYMT(IBP)
            V(LPSYMT(IBP)) = V(LPSYMT(IBP))-COF*YPSYMT(IBP)
            W(LPSYMT(IBP)) = W(LPSYMT(IBP))-COF*ZPSYMT(IBP)
 6040    CONTINUE
C
C
C
C        - RETURN IF ALL THREE EQUATIONS ARE CONVERGED -
C
         IF(IRESU*IRESV*IRESW.EQ.1) THEN
            RETURN
         ENDIF
C
C
CCC    SET R.H.S THAT CORRESPONDS TO THE DIRECLET BOUNDARY NODES
C
C
*POPTION INDEP(RHSU,RHSV,RHSW)
C*$*ASSERT PERMUTATION ( LPFIX )
!ocl norecurrence(RHSU,RHSV,RHSW)
         DO 6200 IBP = 1 , NPFIX
            IP=LPFIX(IBP)
            RHSU(IP) = U(IP)
            RHSV(IP) = V(IP)
            RHSW(IP) = W(IP)
 6200    CONTINUE
C
 5000 CONTINUE   
C
      RETURN
      END
