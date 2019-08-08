ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1 = -(MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_2 = (2.000000D+00*MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_3 = -(MDL_EE*MDL_COMPLEXI)
      GC_76 = (MDL_EE*MDL_COMPLEXI*MDL_SW)/(3.000000D+00*MDL_CW)
      GC_77 = (-2.000000D+00*MDL_EE*MDL_COMPLEXI*MDL_SW)/(3.000000D+00
     $ *MDL_CW)
      GC_78 = (MDL_EE*MDL_COMPLEXI*MDL_SW)/MDL_CW
      GC_80 = -(MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
     $ -(MDL_EE*MDL_COMPLEXI*MDL_SW)/(6.000000D+00*MDL_CW)
      GC_81 = (MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
     $ -(MDL_EE*MDL_COMPLEXI*MDL_SW)/(6.000000D+00*MDL_CW)
      GC_82 = -(MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
     $ +(MDL_EE*MDL_COMPLEXI*MDL_SW)/(2.000000D+00*MDL_CW)
      GC_102 = MDL_EE__EXP__2*MDL_COMPLEXI*MDL_VEV+(MDL_CW__EXP__2
     $ *MDL_EE__EXP__2*MDL_COMPLEXI*MDL_VEV)/(2.000000D+00
     $ *MDL_SW__EXP__2)+(MDL_EE__EXP__2*MDL_COMPLEXI*MDL_SW__EXP__2
     $ *MDL_VEV)/(2.000000D+00*MDL_CW__EXP__2)
      GC_105 = -((MDL_COMPLEXI*MDL_YC)/MDL_SQRT__2)
      GC_108 = -((MDL_COMPLEXI*MDL_YDO)/MDL_SQRT__2)
      GC_116 = -((MDL_COMPLEXI*MDL_YM)/MDL_SQRT__2)
      GC_118 = -((MDL_COMPLEXI*MDL_YS)/MDL_SQRT__2)
      GC_125 = -((MDL_COMPLEXI*MDL_YUP)/MDL_SQRT__2)
      GC_14 = MDL_COMPLEXI*MDL_GUBMU
      END
