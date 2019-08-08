      SUBROUTINE SMATRIX1(P,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.6.1, 2017-12-12
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO for Madevent Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: c c~ > z > mu+ mu- mu+ mu- WEIGHTED<=8 @1
C     
      USE DISCRETESAMPLER
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'genps.inc'
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'maxamps.inc'
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=74)
      INTEGER    NDIAGS
      PARAMETER (NDIAGS=74)
      INTEGER    THEL
      PARAMETER (THEL=2*NCOMB)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
C     
C     global (due to reading writting) 
C     
      LOGICAL GOODHEL(NCOMB,2)
      INTEGER NTRY(2)
      COMMON/BLOCK_GOODHEL/NTRY,GOODHEL
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB)
      INTEGER ISHEL(2)
      REAL*8 T,MATRIX1
      REAL*8 R,SUMHEL,TS(NCOMB)
      INTEGER I,IDEN
      INTEGER JC(NEXTERNAL),II
      REAL*8 HWGT, XTOT, XTRY, XREJ, XR, YFRAC(0:NCOMB)
      INTEGER NGOOD(2), IGOOD(NCOMB,2)
      INTEGER JHEL(2), J, JJ
      INTEGER THIS_NTRY(2)
      SAVE THIS_NTRY
      DATA THIS_NTRY /0,0/
C     This is just to temporarily store the reference grid for
C      helicity of the DiscreteSampler so as to obtain its number of
C      entries with ref_helicity_grid%n_tot_entries
      TYPE(SAMPLEDDIMENSION) REF_HELICITY_GRID
C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2,       JAMP2

      CHARACTER*101         HEL_BUFF
      COMMON/TO_HELICITY/  HEL_BUFF

      INTEGER IMIRROR
      COMMON/TO_MIRROR/ IMIRROR

      REAL*8 POL(2)
      COMMON/TO_POLARIZATION/ POL

      INTEGER          ISUM_HEL
      LOGICAL                    MULTI_CHANNEL
      COMMON/TO_MATRIX/ISUM_HEL, MULTI_CHANNEL
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, ICONFIG
      INTEGER SUBDIAG(MAXSPROC),IB(2)
      COMMON/TO_SUB_DIAG/SUBDIAG,IB
      DATA XTRY, XREJ /0,0/
      DATA NGOOD /0,0/
      DATA ISHEL/0,0/
      SAVE YFRAC, IGOOD, JHEL
      DATA (NHEL(I,   1),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA IDEN/144/

C     To be able to control when the matrix<i> subroutine can add
C      entries to the grid for the MC over helicity configuration
      LOGICAL ALLOW_HELICITY_GRID_ENTRIES
      COMMON/TO_ALLOW_HELICITY_GRID_ENTRIES/ALLOW_HELICITY_GRID_ENTRIES

C     ----------
C     BEGIN CODE
C     ----------
      NTRY(IMIRROR)=NTRY(IMIRROR)+1
      THIS_NTRY(IMIRROR) = THIS_NTRY(IMIRROR)+1
      DO I=1,NEXTERNAL
        JC(I) = +1
      ENDDO

      IF (MULTI_CHANNEL) THEN
        DO I=1,NDIAGS
          AMP2(I)=0D0
        ENDDO
        JAMP2(0)=1
        DO I=1,INT(JAMP2(0))
          JAMP2(I)=0D0
        ENDDO
      ENDIF
      ANS = 0D0
      WRITE(HEL_BUFF,'(20I5)') (0,I=1,NEXTERNAL)
      DO I=1,NCOMB
        TS(I)=0D0
      ENDDO

        !   If the helicity grid status is 0, this means that it is not yet initialized.
        !   If HEL_PICKED==-1, this means that calls to other matrix<i> where in initialization mode as well for the helicity.
      IF ((ISHEL(IMIRROR).EQ.0.AND.ISUM_HEL.EQ.0).OR.(DS_GET_DIM_STATUS
     $('Helicity').EQ.0).OR.(HEL_PICKED.EQ.-1)) THEN
        DO I=1,NCOMB
          IF (GOODHEL(I,IMIRROR) .OR. NTRY(IMIRROR).LE.MAXTRIES.OR.(ISU
     $M_HEL.NE.0).OR.THIS_NTRY(IMIRROR).LE.2) THEN
            T=MATRIX1(P ,NHEL(1,I),JC(1))
            DO JJ=1,NINCOMING
              IF(POL(JJ).NE.1D0.AND.NHEL(JJ,I).EQ.INT(SIGN(1D0,POL(JJ))
     $         )) THEN
                T=T*ABS(POL(JJ))
              ELSE IF(POL(JJ).NE.1D0)THEN
                T=T*(2D0-ABS(POL(JJ)))
              ENDIF
            ENDDO
            IF (ISUM_HEL.NE.0.AND.DS_GET_DIM_STATUS('Helicity')
     $       .EQ.0.AND.ALLOW_HELICITY_GRID_ENTRIES) THEN
              CALL DS_ADD_ENTRY('Helicity',I,T)
            ENDIF
            ANS=ANS+DABS(T)
            TS(I)=T
          ENDIF
        ENDDO
        IF(NTRY(IMIRROR).EQ.(MAXTRIES+1)) THEN
          CALL RESET_CUMULATIVE_VARIABLE()  ! avoid biais of the initialization
        ENDIF
        IF (ISUM_HEL.NE.0) THEN
            !         We set HEL_PICKED to -1 here so that later on, the call to DS_add_point in dsample.f does not add anything to the grid since it was already done here.
          HEL_PICKED = -1
            !         For safety, hardset the helicity sampling jacobian to 0.0d0 to make sure it is not .
          HEL_JACOBIAN   = 1.0D0
            !         We don't want to re-update the helicity grid if it was already updated by another matrix<i>, so we make sure that the reference grid is empty.
          REF_HELICITY_GRID = DS_GET_DIMENSION(REF_GRID,'Helicity')
          IF((DS_GET_DIM_STATUS('Helicity').EQ.1).AND.(REF_HELICITY_GRI
     $D%N_TOT_ENTRIES.EQ.0)) THEN
              !           If we finished the initialization we can update the grid so as to start sampling over it.
              !           However the grid will now be filled by dsample with different kind of weights (including pdf, flux, etc...) so by setting the grid_mode of the reference grid to 'initialization' we make sure it will be overwritten (as opposed to 'combined') by the running grid at the next update.
            CALL DS_UPDATE_GRID('Helicity')
            CALL DS_SET_GRID_MODE('Helicity','init')
          ENDIF
        ELSE
          JHEL(IMIRROR) = 1
          IF(NTRY(IMIRROR).LE.MAXTRIES.OR.THIS_NTRY(IMIRROR).LE.2)THEN
            DO I=1,NCOMB
              IF (.NOT.GOODHEL(I,IMIRROR) .AND. (DABS(TS(I)).GT.ANS
     $         *LIMHEL/NCOMB)) THEN
                GOODHEL(I,IMIRROR)=.TRUE.
                NGOOD(IMIRROR) = NGOOD(IMIRROR) +1
                IGOOD(NGOOD(IMIRROR),IMIRROR) = I
                PRINT *,'Added good helicity ',I,TS(I)*NCOMB/ANS,' in'
     $           //' event ',NTRY(IMIRROR), 'local:',THIS_NTRY(IMIRROR)
              ENDIF
            ENDDO
          ENDIF
          IF(NTRY(IMIRROR).EQ.MAXTRIES)THEN
            ISHEL(IMIRROR)=MIN(ISUM_HEL,NGOOD(IMIRROR))
          ENDIF
        ENDIF
      ELSE  ! random helicity 
C       The helicity configuration was chosen already by genps and put
C        in a common block defined in genps.inc.
        I = HEL_PICKED

        T=MATRIX1(P ,NHEL(1,I),JC(1))

        DO JJ=1,NINCOMING
          IF(POL(JJ).NE.1D0.AND.NHEL(JJ,I).EQ.INT(SIGN(1D0,POL(JJ))))
     $      THEN
            T=T*ABS(POL(JJ))
          ELSE IF(POL(JJ).NE.1D0)THEN
            T=T*(2D0-ABS(POL(JJ)))
          ENDIF
        ENDDO
C       Always one helicity at a time
        ANS = T
C       Include the Jacobian from helicity sampling
        ANS = ANS * HEL_JACOBIAN

        WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
      ENDIF
      IF (ANS.NE.0D0.AND.(ISUM_HEL .NE. 1.OR.HEL_PICKED.EQ.-1)) THEN
        CALL RANMAR(R)
        SUMHEL=0D0
        DO I=1,NCOMB
          SUMHEL=SUMHEL+DABS(TS(I))/ANS
          IF(R.LT.SUMHEL)THEN
            WRITE(HEL_BUFF,'(20i5)')(NHEL(II,I),II=1,NEXTERNAL)
C           Set right sign for ANS, based on sign of chosen helicity
            ANS=DSIGN(ANS,TS(I))
            GOTO 10
          ENDIF
        ENDDO
 10     CONTINUE
      ENDIF
      IF (MULTI_CHANNEL) THEN
        XTOT=0D0
        DO I=1,NDIAGS
          XTOT=XTOT+AMP2(I)
        ENDDO
        IF (XTOT.NE.0D0) THEN
          ANS=ANS*AMP2(SUBDIAG(1))/XTOT
        ELSE IF(ANS.NE.0D0) THEN
          WRITE(*,*) 'Problem in the multi-channeling. All amp2 are'
     $     //' zero but not the total matrix-element'
          STOP 1
        ENDIF
      ENDIF
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX1(P,NHEL,IC)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.6.1, 2017-12-12
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: c c~ > z > mu+ mu- mu+ mu- WEIGHTED<=8 @1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=74)
      INCLUDE 'genps.inc'
      INCLUDE 'nexternal.inc'
      INCLUDE 'maxamps.inc'
      INTEGER    NWAVEFUNCS,     NCOLOR
      PARAMETER (NWAVEFUNCS=26, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=1, NSQAMPSO=1)
      LOGICAL CHOSEN_SO_CONFIGS(NSQAMPSO)
      DATA CHOSEN_SO_CONFIGS/.TRUE./
      SAVE CHOSEN_SO_CONFIGS
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,M,N
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR,NAMPSO)
      COMPLEX*16 W(6,NWAVEFUNCS)
C     Needed for v4 models
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     FUNCTION
C     
      INTEGER SQSOINDEX1
C     
C     GLOBAL VARIABLES
C     
      DOUBLE PRECISION AMP2(MAXAMPS), JAMP2(0:MAXFLOW)
      COMMON/TO_AMPS/  AMP2,       JAMP2
      INCLUDE 'coupl.inc'
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    3/
C     1 T(2,1)
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,1),MDL_MC,NHEL(1),+1*IC(1),W(1,1))
      CALL OXXXXX(P(0,2),MDL_MC,NHEL(2),-1*IC(2),W(1,2))
      CALL IXXXXX(P(0,3),MDL_MMU,NHEL(3),-1*IC(3),W(1,3))
      CALL OXXXXX(P(0,4),MDL_MMU,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),MDL_MMU,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX(P(0,6),MDL_MMU,NHEL(6),+1*IC(6),W(1,6))
      CALL FFS5_3(W(1,1),W(1,2),GC_105,MDL_MH,MDL_WH,W(1,7))
      CALL FFV2_3_3(W(1,3),W(1,4),GC_82,GC_78,MDL_MZ,MDL_WZ,W(1,8))
      CALL FFS5_2(W(1,5),W(1,7),GC_116,MDL_MMU,ZERO,W(1,9))
C     Amplitude(s) for diagram number 1
      CALL FFV2_3_0(W(1,9),W(1,6),W(1,8),GC_82,GC_78,AMP(1))
      CALL FFS5_1(W(1,6),W(1,7),GC_116,MDL_MMU,ZERO,W(1,10))
C     Amplitude(s) for diagram number 2
      CALL FFV2_3_0(W(1,5),W(1,10),W(1,8),GC_82,GC_78,AMP(2))
      CALL FFV1P0_3(W(1,1),W(1,2),GC_2,ZERO,ZERO,W(1,11))
      CALL FFV1_2(W(1,5),W(1,11),GC_3,MDL_MMU,ZERO,W(1,12))
C     Amplitude(s) for diagram number 3
      CALL FFV2_3_0(W(1,12),W(1,6),W(1,8),GC_82,GC_78,AMP(3))
      CALL FFV1_1(W(1,6),W(1,11),GC_3,MDL_MMU,ZERO,W(1,13))
C     Amplitude(s) for diagram number 4
      CALL FFV2_3_0(W(1,5),W(1,13),W(1,8),GC_82,GC_78,AMP(4))
      CALL FFV2_3_3(W(1,1),W(1,2),GC_81,GC_77,MDL_MZ,MDL_WZ,W(1,14))
      CALL FFS1_3(W(1,3),W(1,4),GC_14,MDL_MUB,MDL_WUB,W(1,15))
      CALL FFV2_3_2(W(1,5),W(1,14),GC_82,GC_78,MDL_MMU,ZERO,W(1,16))
C     Amplitude(s) for diagram number 5
      CALL FFS1_0(W(1,16),W(1,6),W(1,15),GC_14,AMP(5))
      CALL FFV2_3_1(W(1,6),W(1,14),GC_82,GC_78,MDL_MMU,ZERO,W(1,17))
C     Amplitude(s) for diagram number 6
      CALL FFS1_0(W(1,5),W(1,17),W(1,15),GC_14,AMP(6))
      CALL FFS5_3(W(1,3),W(1,4),GC_116,MDL_MH,MDL_WH,W(1,15))
C     Amplitude(s) for diagram number 7
      CALL FFS5_0(W(1,16),W(1,6),W(1,15),GC_116,AMP(7))
C     Amplitude(s) for diagram number 8
      CALL FFS5_0(W(1,5),W(1,17),W(1,15),GC_116,AMP(8))
      CALL FFV1P0_3(W(1,3),W(1,4),GC_3,ZERO,ZERO,W(1,18))
C     Amplitude(s) for diagram number 9
      CALL FFV1_0(W(1,16),W(1,6),W(1,18),GC_3,AMP(9))
C     Amplitude(s) for diagram number 10
      CALL FFV1_0(W(1,5),W(1,17),W(1,18),GC_3,AMP(10))
C     Amplitude(s) for diagram number 11
      CALL FFV2_3_0(W(1,16),W(1,6),W(1,8),GC_82,GC_78,AMP(11))
C     Amplitude(s) for diagram number 12
      CALL FFV2_3_0(W(1,5),W(1,17),W(1,8),GC_82,GC_78,AMP(12))
      CALL FFV2_3_3(W(1,5),W(1,6),GC_82,GC_78,MDL_MZ,MDL_WZ,W(1,19))
C     Amplitude(s) for diagram number 13
      CALL VVS1_0(W(1,8),W(1,19),W(1,7),GC_102,AMP(13))
C     Amplitude(s) for diagram number 14
      CALL VVS1_0(W(1,14),W(1,19),W(1,15),GC_102,AMP(14))
      CALL FFS5_3(W(1,5),W(1,6),GC_116,MDL_MH,MDL_WH,W(1,20))
C     Amplitude(s) for diagram number 15
      CALL VVS1_0(W(1,14),W(1,8),W(1,20),GC_102,AMP(15))
      CALL FFV2_3_3(W(1,3),W(1,6),GC_82,GC_78,MDL_MZ,MDL_WZ,W(1,21))
      CALL FFS5_1(W(1,4),W(1,7),GC_116,MDL_MMU,ZERO,W(1,22))
C     Amplitude(s) for diagram number 16
      CALL FFV2_3_0(W(1,5),W(1,22),W(1,21),GC_82,GC_78,AMP(16))
C     Amplitude(s) for diagram number 17
      CALL FFV2_3_0(W(1,9),W(1,4),W(1,21),GC_82,GC_78,AMP(17))
      CALL FFV1_1(W(1,4),W(1,11),GC_3,MDL_MMU,ZERO,W(1,9))
C     Amplitude(s) for diagram number 18
      CALL FFV2_3_0(W(1,5),W(1,9),W(1,21),GC_82,GC_78,AMP(18))
C     Amplitude(s) for diagram number 19
      CALL FFV2_3_0(W(1,12),W(1,4),W(1,21),GC_82,GC_78,AMP(19))
      CALL FFS1_3(W(1,3),W(1,6),GC_14,MDL_MUB,MDL_WUB,W(1,12))
      CALL FFV2_3_1(W(1,4),W(1,14),GC_82,GC_78,MDL_MMU,ZERO,W(1,23))
C     Amplitude(s) for diagram number 20
      CALL FFS1_0(W(1,5),W(1,23),W(1,12),GC_14,AMP(20))
C     Amplitude(s) for diagram number 21
      CALL FFS1_0(W(1,16),W(1,4),W(1,12),GC_14,AMP(21))
      CALL FFS5_3(W(1,3),W(1,6),GC_116,MDL_MH,MDL_WH,W(1,12))
C     Amplitude(s) for diagram number 22
      CALL FFS5_0(W(1,5),W(1,23),W(1,12),GC_116,AMP(22))
C     Amplitude(s) for diagram number 23
      CALL FFS5_0(W(1,16),W(1,4),W(1,12),GC_116,AMP(23))
      CALL FFV1P0_3(W(1,3),W(1,6),GC_3,ZERO,ZERO,W(1,24))
C     Amplitude(s) for diagram number 24
      CALL FFV1_0(W(1,5),W(1,23),W(1,24),GC_3,AMP(24))
C     Amplitude(s) for diagram number 25
      CALL FFV1_0(W(1,16),W(1,4),W(1,24),GC_3,AMP(25))
C     Amplitude(s) for diagram number 26
      CALL FFV2_3_0(W(1,5),W(1,23),W(1,21),GC_82,GC_78,AMP(26))
C     Amplitude(s) for diagram number 27
      CALL FFV2_3_0(W(1,16),W(1,4),W(1,21),GC_82,GC_78,AMP(27))
      CALL FFV2_3_3(W(1,5),W(1,4),GC_82,GC_78,MDL_MZ,MDL_WZ,W(1,16))
C     Amplitude(s) for diagram number 28
      CALL VVS1_0(W(1,21),W(1,16),W(1,7),GC_102,AMP(28))
C     Amplitude(s) for diagram number 29
      CALL VVS1_0(W(1,14),W(1,16),W(1,12),GC_102,AMP(29))
      CALL FFS5_3(W(1,5),W(1,4),GC_116,MDL_MH,MDL_WH,W(1,25))
C     Amplitude(s) for diagram number 30
      CALL VVS1_0(W(1,14),W(1,21),W(1,25),GC_102,AMP(30))
      CALL FFS5_2(W(1,3),W(1,7),GC_116,MDL_MMU,ZERO,W(1,26))
C     Amplitude(s) for diagram number 31
      CALL FFV2_3_0(W(1,26),W(1,6),W(1,16),GC_82,GC_78,AMP(31))
C     Amplitude(s) for diagram number 32
      CALL FFV2_3_0(W(1,3),W(1,10),W(1,16),GC_82,GC_78,AMP(32))
      CALL FFV1_2(W(1,3),W(1,11),GC_3,MDL_MMU,ZERO,W(1,10))
C     Amplitude(s) for diagram number 33
      CALL FFV2_3_0(W(1,10),W(1,6),W(1,16),GC_82,GC_78,AMP(33))
C     Amplitude(s) for diagram number 34
      CALL FFV2_3_0(W(1,3),W(1,13),W(1,16),GC_82,GC_78,AMP(34))
      CALL FFS1_3(W(1,5),W(1,4),GC_14,MDL_MUB,MDL_WUB,W(1,13))
      CALL FFV2_3_2(W(1,3),W(1,14),GC_82,GC_78,MDL_MMU,ZERO,W(1,11))
C     Amplitude(s) for diagram number 35
      CALL FFS1_0(W(1,11),W(1,6),W(1,13),GC_14,AMP(35))
C     Amplitude(s) for diagram number 36
      CALL FFS1_0(W(1,3),W(1,17),W(1,13),GC_14,AMP(36))
C     Amplitude(s) for diagram number 37
      CALL FFS5_0(W(1,11),W(1,6),W(1,25),GC_116,AMP(37))
C     Amplitude(s) for diagram number 38
      CALL FFS5_0(W(1,3),W(1,17),W(1,25),GC_116,AMP(38))
      CALL FFV1P0_3(W(1,5),W(1,4),GC_3,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 39
      CALL FFV1_0(W(1,11),W(1,6),W(1,13),GC_3,AMP(39))
C     Amplitude(s) for diagram number 40
      CALL FFV1_0(W(1,3),W(1,17),W(1,13),GC_3,AMP(40))
C     Amplitude(s) for diagram number 41
      CALL FFV2_3_0(W(1,11),W(1,6),W(1,16),GC_82,GC_78,AMP(41))
C     Amplitude(s) for diagram number 42
      CALL FFV2_3_0(W(1,3),W(1,17),W(1,16),GC_82,GC_78,AMP(42))
C     Amplitude(s) for diagram number 43
      CALL FFV2_3_0(W(1,26),W(1,4),W(1,19),GC_82,GC_78,AMP(43))
C     Amplitude(s) for diagram number 44
      CALL FFV2_3_0(W(1,3),W(1,22),W(1,19),GC_82,GC_78,AMP(44))
C     Amplitude(s) for diagram number 45
      CALL FFV2_3_0(W(1,10),W(1,4),W(1,19),GC_82,GC_78,AMP(45))
C     Amplitude(s) for diagram number 46
      CALL FFV2_3_0(W(1,3),W(1,9),W(1,19),GC_82,GC_78,AMP(46))
      CALL FFS1_3(W(1,5),W(1,6),GC_14,MDL_MUB,MDL_WUB,W(1,9))
C     Amplitude(s) for diagram number 47
      CALL FFS1_0(W(1,11),W(1,4),W(1,9),GC_14,AMP(47))
C     Amplitude(s) for diagram number 48
      CALL FFS1_0(W(1,3),W(1,23),W(1,9),GC_14,AMP(48))
C     Amplitude(s) for diagram number 49
      CALL FFS5_0(W(1,11),W(1,4),W(1,20),GC_116,AMP(49))
C     Amplitude(s) for diagram number 50
      CALL FFS5_0(W(1,3),W(1,23),W(1,20),GC_116,AMP(50))
      CALL FFV1P0_3(W(1,5),W(1,6),GC_3,ZERO,ZERO,W(1,9))
C     Amplitude(s) for diagram number 51
      CALL FFV1_0(W(1,11),W(1,4),W(1,9),GC_3,AMP(51))
C     Amplitude(s) for diagram number 52
      CALL FFV1_0(W(1,3),W(1,23),W(1,9),GC_3,AMP(52))
C     Amplitude(s) for diagram number 53
      CALL FFV2_3_0(W(1,11),W(1,4),W(1,19),GC_82,GC_78,AMP(53))
C     Amplitude(s) for diagram number 54
      CALL FFV2_3_0(W(1,3),W(1,23),W(1,19),GC_82,GC_78,AMP(54))
      CALL FFS5_2(W(1,1),W(1,15),GC_105,MDL_MC,ZERO,W(1,23))
C     Amplitude(s) for diagram number 55
      CALL FFV2_3_0(W(1,23),W(1,2),W(1,19),GC_81,GC_77,AMP(55))
      CALL FFV2_3_2(W(1,1),W(1,19),GC_81,GC_77,MDL_MC,ZERO,W(1,23))
C     Amplitude(s) for diagram number 56
      CALL FFS5_0(W(1,23),W(1,2),W(1,15),GC_105,AMP(56))
      CALL FFV1_2(W(1,1),W(1,18),GC_2,MDL_MC,ZERO,W(1,15))
C     Amplitude(s) for diagram number 57
      CALL FFV2_3_0(W(1,15),W(1,2),W(1,19),GC_81,GC_77,AMP(57))
C     Amplitude(s) for diagram number 58
      CALL FFV1_0(W(1,23),W(1,2),W(1,18),GC_2,AMP(58))
      CALL FFV2_3_2(W(1,1),W(1,8),GC_81,GC_77,MDL_MC,ZERO,W(1,18))
C     Amplitude(s) for diagram number 59
      CALL FFS5_0(W(1,18),W(1,2),W(1,20),GC_105,AMP(59))
      CALL FFS5_2(W(1,1),W(1,20),GC_105,MDL_MC,ZERO,W(1,15))
C     Amplitude(s) for diagram number 60
      CALL FFV2_3_0(W(1,15),W(1,2),W(1,8),GC_81,GC_77,AMP(60))
C     Amplitude(s) for diagram number 61
      CALL FFV1_0(W(1,18),W(1,2),W(1,9),GC_2,AMP(61))
      CALL FFV1_2(W(1,1),W(1,9),GC_2,MDL_MC,ZERO,W(1,15))
C     Amplitude(s) for diagram number 62
      CALL FFV2_3_0(W(1,15),W(1,2),W(1,8),GC_81,GC_77,AMP(62))
C     Amplitude(s) for diagram number 63
      CALL FFV2_3_0(W(1,18),W(1,2),W(1,19),GC_81,GC_77,AMP(63))
C     Amplitude(s) for diagram number 64
      CALL FFV2_3_0(W(1,23),W(1,2),W(1,8),GC_81,GC_77,AMP(64))
      CALL FFS5_2(W(1,1),W(1,12),GC_105,MDL_MC,ZERO,W(1,23))
C     Amplitude(s) for diagram number 65
      CALL FFV2_3_0(W(1,23),W(1,2),W(1,16),GC_81,GC_77,AMP(65))
      CALL FFV2_3_2(W(1,1),W(1,16),GC_81,GC_77,MDL_MC,ZERO,W(1,23))
C     Amplitude(s) for diagram number 66
      CALL FFS5_0(W(1,23),W(1,2),W(1,12),GC_105,AMP(66))
      CALL FFV1_2(W(1,1),W(1,24),GC_2,MDL_MC,ZERO,W(1,12))
C     Amplitude(s) for diagram number 67
      CALL FFV2_3_0(W(1,12),W(1,2),W(1,16),GC_81,GC_77,AMP(67))
C     Amplitude(s) for diagram number 68
      CALL FFV1_0(W(1,23),W(1,2),W(1,24),GC_2,AMP(68))
      CALL FFV2_3_2(W(1,1),W(1,21),GC_81,GC_77,MDL_MC,ZERO,W(1,24))
C     Amplitude(s) for diagram number 69
      CALL FFS5_0(W(1,24),W(1,2),W(1,25),GC_105,AMP(69))
      CALL FFS5_2(W(1,1),W(1,25),GC_105,MDL_MC,ZERO,W(1,12))
C     Amplitude(s) for diagram number 70
      CALL FFV2_3_0(W(1,12),W(1,2),W(1,21),GC_81,GC_77,AMP(70))
C     Amplitude(s) for diagram number 71
      CALL FFV1_0(W(1,24),W(1,2),W(1,13),GC_2,AMP(71))
      CALL FFV1_2(W(1,1),W(1,13),GC_2,MDL_MC,ZERO,W(1,12))
C     Amplitude(s) for diagram number 72
      CALL FFV2_3_0(W(1,12),W(1,2),W(1,21),GC_81,GC_77,AMP(72))
C     Amplitude(s) for diagram number 73
      CALL FFV2_3_0(W(1,24),W(1,2),W(1,16),GC_81,GC_77,AMP(73))
C     Amplitude(s) for diagram number 74
      CALL FFV2_3_0(W(1,23),W(1,2),W(1,21),GC_81,GC_77,AMP(74))
C     JAMPs contributing to orders ALL_ORDERS=1
      JAMP(1,1)=+AMP(1)+AMP(2)+AMP(3)+AMP(4)+AMP(5)+AMP(6)+AMP(7)
     $ +AMP(8)+AMP(9)+AMP(10)+AMP(11)+AMP(12)+AMP(13)+AMP(14)+AMP(15)
     $ -AMP(16)-AMP(17)-AMP(18)-AMP(19)-AMP(20)-AMP(21)-AMP(22)-AMP(23)
     $ -AMP(24)-AMP(25)-AMP(26)-AMP(27)-AMP(28)-AMP(29)-AMP(30)-AMP(31)
     $ -AMP(32)-AMP(33)-AMP(34)-AMP(35)-AMP(36)-AMP(37)-AMP(38)-AMP(39)
     $ -AMP(40)-AMP(41)-AMP(42)+AMP(43)+AMP(44)+AMP(45)+AMP(46)+AMP(47)
     $ +AMP(48)+AMP(49)+AMP(50)+AMP(51)+AMP(52)+AMP(53)+AMP(54)+AMP(55)
     $ +AMP(56)+AMP(57)+AMP(58)+AMP(59)+AMP(60)+AMP(61)+AMP(62)+AMP(63)
     $ +AMP(64)-AMP(65)-AMP(66)-AMP(67)-AMP(68)-AMP(69)-AMP(70)-AMP(71)
     $ -AMP(72)-AMP(73)-AMP(74)

      MATRIX1 = 0.D0
      DO M = 1, NAMPSO
        DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
            ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
          ENDDO
          DO N = 1, NAMPSO
            IF (CHOSEN_SO_CONFIGS(SQSOINDEX1(M,N))) THEN
              MATRIX1 = MATRIX1 + ZTEMP*DCONJG(JAMP(I,N))/DENOM(I)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      AMP2(1)=AMP2(1)+AMP(1)*DCONJG(AMP(1))
      AMP2(2)=AMP2(2)+AMP(2)*DCONJG(AMP(2))
      AMP2(3)=AMP2(3)+AMP(3)*DCONJG(AMP(3))
      AMP2(4)=AMP2(4)+AMP(4)*DCONJG(AMP(4))
      AMP2(5)=AMP2(5)+AMP(5)*DCONJG(AMP(5))
      AMP2(6)=AMP2(6)+AMP(6)*DCONJG(AMP(6))
      AMP2(7)=AMP2(7)+AMP(7)*DCONJG(AMP(7))
      AMP2(8)=AMP2(8)+AMP(8)*DCONJG(AMP(8))
      AMP2(9)=AMP2(9)+AMP(9)*DCONJG(AMP(9))
      AMP2(10)=AMP2(10)+AMP(10)*DCONJG(AMP(10))
      AMP2(11)=AMP2(11)+AMP(11)*DCONJG(AMP(11))
      AMP2(12)=AMP2(12)+AMP(12)*DCONJG(AMP(12))
      AMP2(13)=AMP2(13)+AMP(13)*DCONJG(AMP(13))
      AMP2(14)=AMP2(14)+AMP(14)*DCONJG(AMP(14))
      AMP2(15)=AMP2(15)+AMP(15)*DCONJG(AMP(15))
      AMP2(16)=AMP2(16)+AMP(16)*DCONJG(AMP(16))
      AMP2(17)=AMP2(17)+AMP(17)*DCONJG(AMP(17))
      AMP2(18)=AMP2(18)+AMP(18)*DCONJG(AMP(18))
      AMP2(19)=AMP2(19)+AMP(19)*DCONJG(AMP(19))
      AMP2(20)=AMP2(20)+AMP(20)*DCONJG(AMP(20))
      AMP2(21)=AMP2(21)+AMP(21)*DCONJG(AMP(21))
      AMP2(22)=AMP2(22)+AMP(22)*DCONJG(AMP(22))
      AMP2(23)=AMP2(23)+AMP(23)*DCONJG(AMP(23))
      AMP2(24)=AMP2(24)+AMP(24)*DCONJG(AMP(24))
      AMP2(25)=AMP2(25)+AMP(25)*DCONJG(AMP(25))
      AMP2(26)=AMP2(26)+AMP(26)*DCONJG(AMP(26))
      AMP2(27)=AMP2(27)+AMP(27)*DCONJG(AMP(27))
      AMP2(28)=AMP2(28)+AMP(28)*DCONJG(AMP(28))
      AMP2(29)=AMP2(29)+AMP(29)*DCONJG(AMP(29))
      AMP2(30)=AMP2(30)+AMP(30)*DCONJG(AMP(30))
      AMP2(31)=AMP2(31)+AMP(31)*DCONJG(AMP(31))
      AMP2(32)=AMP2(32)+AMP(32)*DCONJG(AMP(32))
      AMP2(33)=AMP2(33)+AMP(33)*DCONJG(AMP(33))
      AMP2(34)=AMP2(34)+AMP(34)*DCONJG(AMP(34))
      AMP2(35)=AMP2(35)+AMP(35)*DCONJG(AMP(35))
      AMP2(36)=AMP2(36)+AMP(36)*DCONJG(AMP(36))
      AMP2(37)=AMP2(37)+AMP(37)*DCONJG(AMP(37))
      AMP2(38)=AMP2(38)+AMP(38)*DCONJG(AMP(38))
      AMP2(39)=AMP2(39)+AMP(39)*DCONJG(AMP(39))
      AMP2(40)=AMP2(40)+AMP(40)*DCONJG(AMP(40))
      AMP2(41)=AMP2(41)+AMP(41)*DCONJG(AMP(41))
      AMP2(42)=AMP2(42)+AMP(42)*DCONJG(AMP(42))
      AMP2(43)=AMP2(43)+AMP(43)*DCONJG(AMP(43))
      AMP2(44)=AMP2(44)+AMP(44)*DCONJG(AMP(44))
      AMP2(45)=AMP2(45)+AMP(45)*DCONJG(AMP(45))
      AMP2(46)=AMP2(46)+AMP(46)*DCONJG(AMP(46))
      AMP2(47)=AMP2(47)+AMP(47)*DCONJG(AMP(47))
      AMP2(48)=AMP2(48)+AMP(48)*DCONJG(AMP(48))
      AMP2(49)=AMP2(49)+AMP(49)*DCONJG(AMP(49))
      AMP2(50)=AMP2(50)+AMP(50)*DCONJG(AMP(50))
      AMP2(51)=AMP2(51)+AMP(51)*DCONJG(AMP(51))
      AMP2(52)=AMP2(52)+AMP(52)*DCONJG(AMP(52))
      AMP2(53)=AMP2(53)+AMP(53)*DCONJG(AMP(53))
      AMP2(54)=AMP2(54)+AMP(54)*DCONJG(AMP(54))
      AMP2(55)=AMP2(55)+AMP(55)*DCONJG(AMP(55))
      AMP2(56)=AMP2(56)+AMP(56)*DCONJG(AMP(56))
      AMP2(57)=AMP2(57)+AMP(57)*DCONJG(AMP(57))
      AMP2(58)=AMP2(58)+AMP(58)*DCONJG(AMP(58))
      AMP2(59)=AMP2(59)+AMP(59)*DCONJG(AMP(59))
      AMP2(60)=AMP2(60)+AMP(60)*DCONJG(AMP(60))
      AMP2(61)=AMP2(61)+AMP(61)*DCONJG(AMP(61))
      AMP2(62)=AMP2(62)+AMP(62)*DCONJG(AMP(62))
      AMP2(63)=AMP2(63)+AMP(63)*DCONJG(AMP(63))
      AMP2(64)=AMP2(64)+AMP(64)*DCONJG(AMP(64))
      AMP2(65)=AMP2(65)+AMP(65)*DCONJG(AMP(65))
      AMP2(66)=AMP2(66)+AMP(66)*DCONJG(AMP(66))
      AMP2(67)=AMP2(67)+AMP(67)*DCONJG(AMP(67))
      AMP2(68)=AMP2(68)+AMP(68)*DCONJG(AMP(68))
      AMP2(69)=AMP2(69)+AMP(69)*DCONJG(AMP(69))
      AMP2(70)=AMP2(70)+AMP(70)*DCONJG(AMP(70))
      AMP2(71)=AMP2(71)+AMP(71)*DCONJG(AMP(71))
      AMP2(72)=AMP2(72)+AMP(72)*DCONJG(AMP(72))
      AMP2(73)=AMP2(73)+AMP(73)*DCONJG(AMP(73))
      AMP2(74)=AMP2(74)+AMP(74)*DCONJG(AMP(74))
      DO I = 1, NCOLOR
        DO M = 1, NAMPSO
          DO N = 1, NAMPSO
            IF (CHOSEN_SO_CONFIGS(SQSOINDEX1(M,N))) THEN
              JAMP2(I)=JAMP2(I)+DABS(DBLE(JAMP(I,M)*DCONJG(JAMP(I,N))))
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END

C     Set of functions to handle the array indices of the split orders


      INTEGER FUNCTION SQSOINDEX1(ORDERINDEXA, ORDERINDEXB)
C     
C     This functions plays the role of the interference matrix. It can
C      be hardcoded or 
C     made more elegant using hashtables if its execution speed ever
C      becomes a relevant
C     factor. From two split order indices, it return the
C      corresponding index in the squared 
C     order canonical ordering.
C     
C     CONSTANTS
C     

      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERINDEXA, ORDERINDEXB
C     
C     LOCAL VARIABLES
C     
      INTEGER I, SQORDERS(NSO)
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      DATA (AMPSPLITORDERS(  1,I),I=  1,  1) /    1/
      COMMON/AMPSPLITORDERS1/AMPSPLITORDERS
C     
C     FUNCTION
C     
      INTEGER SOINDEX_FOR_SQUARED_ORDERS1
C     
C     BEGIN CODE
C     
      DO I=1,NSO
        SQORDERS(I)=AMPSPLITORDERS(ORDERINDEXA,I)+AMPSPLITORDERS(ORDERI
     $NDEXB,I)
      ENDDO
      SQSOINDEX1=SOINDEX_FOR_SQUARED_ORDERS1(SQORDERS)
      END

      INTEGER FUNCTION SOINDEX_FOR_SQUARED_ORDERS1(ORDERS)
C     
C     This functions returns the integer index identifying the squared
C      split orders list passed in argument which corresponds to the
C      values of the following list of couplings (and in this order).
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO, NAMPSO
      PARAMETER (NSO=1, NSQSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER SQSPLITORDERS(NSQSO,NSO)
      DATA (SQSPLITORDERS(  1,I),I=  1,  1) /    2/
      COMMON/SQPLITORDERS1/SQPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NSQSO
        DO J=1,NSO
          IF (ORDERS(J).NE.SQSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_SQUARED_ORDERS1 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping in function'
      WRITE(*,*) 'SOINDEX_FOR_SQUARED_ORDERS1'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

      SUBROUTINE GET_NSQSO_BORN1(NSQSO)
C     
C     Simple subroutine returning the number of squared split order
C     contributions returned when calling smatrix_split_orders 
C     

      INTEGER    NSQUAREDSO
      PARAMETER  (NSQUAREDSO=1)

      INTEGER NSQSO

      NSQSO=NSQUAREDSO

      END

C     This is the inverse subroutine of SOINDEX_FOR_SQUARED_ORDERS.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_SQUARED_ORDERS_FOR_SOINDEX1(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the squared
C      split order index in argument. Order values correspond to
C      following list of couplings (and in this order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO
      PARAMETER (NSO=1, NSQSO=1)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER SQPLITORDERS(NSQSO,NSO)
      COMMON/SQPLITORDERS1/SQPLITORDERS
C     
C     BEGIN CODE
C     
      IF (SOINDEX.GT.0.AND.SOINDEX.LE.NSQSO) THEN
        DO I=1,NSO
          ORDERS(I) =  SQPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function GET_SQUARED_ORDERS_FOR_SOIN'
     $ //'DEX1'
      WRITE(*,*) 'Could not find squared orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This is the inverse subroutine of getting amplitude SO orders.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_ORDERS_FOR_AMPSOINDEX1(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the split order
C      index in argument. Order values correspond to following list of
C      couplings (and in this order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS1/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      IF (SOINDEX.GT.0.AND.SOINDEX.LE.NAMPSO) THEN
        DO I=1,NSO
          ORDERS(I) =  AMPSPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function GET_ORDERS_FOR_AMPSOINDEX1'
      WRITE(*,*) 'Could not find amplitude split orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This function is not directly useful, but included for
C      completeness
      INTEGER FUNCTION SOINDEX_FOR_AMPORDERS1(ORDERS)
C     
C     This functions returns the integer index identifying the
C      amplitude split orders passed in argument which correspond to
C      the values of the following list of couplings (and in this
C      order):
C     []
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=1, NAMPSO=1)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS1/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NAMPSO
        DO J=1,NSO
          IF (ORDERS(J).NE.AMPSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_AMPORDERS1 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping function SOINDEX_FOR_AMPORDERS1'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

