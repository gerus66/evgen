C########################################################################
C#                                                                     ## 
C#                                                                     ##
C#   ##   ##  ##      #####  #####    ####    ##        THIS  FILE     ##
C#   ##  ##   ##      ##     ##  ##  ##       ##      CONTAINES  THE   ##
C#   ## ##    ##      ####   #####    #####   ##         FOLLOWING     ##
C#   #####    ##      ##     ##  ##       ##  ##      SUBROUTINES AND  ##
C#   ## ###   ##  ##  ##     ##  ##  ##   ##  ##         FUNCTIONS :   ##
C#   ##  ###  ######  #####  #####    #####   ##                       ##
C#                                                                     ##
C#                                                                     ##
C#  1)  Function  H(AJ-REAL*8)          H=Sqrt(2AJ+1)                  ##
C#  2)  Function  Z(J-INTEGER)          Z=(-)^J                        ##
C#  3)  Function  Y(AJ-REAL*8)          Y=(-)^Int(AJ)                  ##
C#  4)  Function  CLB(All-REAL*8)                                      ##
C#  5)  Function  S3J(All-REAL*8)                                      ##
C#  6)  Function  S6J(All-REAL*8)                                      ##
C#  7)  Function  S9J(All-REAL*8)                                      ##
C#  8)  Function  SFEFUN(L,M,TET,FI)    Re( Y  (TET,FI) )              ##
C#  9)  Function  PLM(L,IABS(M),COS(TETA))                             ##
C#  10) Function  YNORM(L,IABS(M))                                     ##
C#      Y(L,M) = (-)**M * YNORM(L,IABS(M)) * PLM(L,IABS(M),COS(TETA))  ##
C#                                                                     ##
C#       date : 30.10.94    7.05.01                                    ##
C#                                                                     ##
C########################################################################
C
      FUNCTION H(AJ)
      USE MYDATA,  ONLY: KIND, CR1, CR2
      REAL(KIND) H, AJ
      H = SQRT( CR2 * AJ + CR1 )
      END
C
C########################################################################
C
      FUNCTION Z(J)
      USE MYDATA,  ONLY: KIND, CR1
      REAL(KIND) Z
	INTEGER J
      IF(J.EQ.0) THEN
        Z = CR1
      ELSE
        Z = ( -CR1 )**IABS(J)
      END IF
      END
C
C########################################################################
C
      FUNCTION Y(AJ)
      USE MYDATA,  ONLY: KIND, CR1
      REAL(KIND) Y, AJ
      J = INT(ABS(AJ)+0.1)
      IF(J.EQ.0) THEN
        Y = CR1
      ELSE
        Y = ( -CR1 )**IABS(J)
      END IF
      END
C
C########################################################################
C
	FUNCTION CLB(J1,M1,J2,M2,J,M)
      USE MYDATA,  ONLY: KIND, CR0, CR1, CR2
C
C-----------------------------------------------------------------------
C...program name
	CHARACTER(5) NAME
	PARAMETER (NAME='CLB')
C...input scalar arguments
	REAL(KIND) J1,M1,J2,M2,J,M
C...output function
	REAL(KIND) CLB
C...external functions
      REAL(KIND) S3J
      EXTERNAL S3J
C
C...executable statements
C
      CLB=SQRT(CR2*J+CR1)*S3J(J1,M1,J2,M2,J,-M)
	IF(CLB.NE.CR0) THEN
	  IF(CR1*INT(J1-J2+M).EQ.J1-J2+M) THEN
	    IF(MOD(INT(J1-J2+M),2).NE.0) CLB=-CLB
          ELSE
	    PRINT *,'  error(',NAME,'):  J1-J2+M  is not integer'
	    PRINT *,'  J1 =  ',J1
	    PRINT *,'  J2 =  ',J2
	    PRINT *,'  M =  ',M
	    STOP
	  END IF
	END IF
C
      RETURN
	END
C
C########################################################################
C
      FUNCTION S3J(J1,M1,J2,M2,J3,M3)
      USE MYDATA,  ONLY: KIND, CR0, CR0D5, CR1, CR2
C
C-----------------------------------------------------------------------
C...program name
        CHARACTER*5 NAME
        PARAMETER (NAME='S3J')
C...parameters
	INTEGER*4 L_M,NLF_M
	PARAMETER (L_M=100,NLF_M=3*L_M+1)
C...input scalar arguments
	REAL(KIND) J1,M1,J2,M2,J3,M3
C...output function
      REAL(KIND) S3J
C...save variables
	LOGICAL*4 INIT
	REAL(KIND) LF(0:NLF_M)
	SAVE INIT,LF
C...local scalars
	INTEGER*4 DJ1,DM1,DJ2,DM2,DJ3,DM3,ILF,NLF,Z,ZMIN,ZMAX
	REAL(KIND) LCF,SGN
C...local arrays
	INTEGER*4 WK(13)
C...constants
	DATA INIT/.TRUE./
C
C...executable statements
C
	IF(INIT) THEN
	  LF(0)=CR0
	  LF(1)=CR0
	  DO ILF=2,NLF_M
	    LF(ILF)=LF(ILF-1)+LOG(CR1*ILF)
	  END DO
	  INIT=.FALSE.
	END IF
C
	IF(J1.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J1 < 0'
	  PRINT *,'  J1 =  ',J1
	  STOP
	END IF
	IF(J2.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J2 < 0'
	  PRINT *,'  J2 =  ',J2
	  STOP
	END IF
	IF(J3.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J3 < 0'
	  PRINT *,'  J3 =  ',J3
	  STOP
	END IF
C
	IF(DBLE(INT(J1)).EQ.J1.AND.CR1*INT(M1).EQ.M1) THEN
	  DJ1=2*INT(J1)
	  DM1=2*INT(M1)
	ELSE
	  IF(CR1*INT(J1)+CR0D5.EQ.J1.AND.
     &       CR1*INT(M1)+SIGN(CR0D5,M1).EQ.M1) THEN
	    DJ1=2*INT(J1)+1
	    DM1=2*INT(M1)+INT(SIGN(CR1,M1))
	  ELSE
	    PRINT *,'  error(',NAME,'):  J1 M1  are wrong'
	    PRINT *,'  J1 =  ',J1
	    PRINT *,'  M1 =  ',M1
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J2).EQ.J2.AND.CR1*INT(M2).EQ.M2) THEN
	  DJ2=2*INT(J2)
	  DM2=2*INT(M2)
	ELSE
	  IF(CR1*INT(J2)+CR0D5.EQ.J2.AND.
     &       CR1*INT(M2)+SIGN(CR0D5,M2).EQ.M2) THEN
	    DJ2=2*INT(J2)+1
	    DM2=2*INT(M2)+INT(SIGN(CR1,M2))
	  ELSE
	    PRINT *,'  error(',NAME,'):  J2 M2  are wrong'
	    PRINT *,'  J2 =  ',J2
	    PRINT *,'  M2 =  ',M2
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J3).EQ.J3.AND.CR1*INT(M3).EQ.M3) THEN
	  DJ3=2*INT(J3)
	  DM3=2*INT(M3)
	ELSE
	  IF(CR1*INT(J3)+CR0D5.EQ.J3.AND.
     &       CR1*INT(M3)+SIGN(CR0D5,M3).EQ.M3) THEN
	    DJ3=2*INT(J3)+1
	    DM3=2*INT(M3)+INT(SIGN(CR1,M3))
	  ELSE
	    PRINT *,'  error(',NAME,'):  J3 M3  are wrong'
	    PRINT *,'  J3 =  ',J3
	    PRINT *,'  M3 =  ',M3
	    STOP
	  END IF
	END IF
C
	IF(IABS(DM1).GT.DJ1) THEN
	  PRINT *,'  warning(',NAME,'):  ABS(M1) > J1'
	  PRINT *,'  J1 =  ',J1
	  PRINT *,'  M1 =  ',M1
          S3J=CR0
	  RETURN
	END IF
	IF(IABS(DM2).GT.DJ2) THEN
	  PRINT *,'  warning(',NAME,'):  ABS(M2) > J2'
	  PRINT *,'  J2 =  ',J2
	  PRINT *,'  M2 =  ',M2
          S3J=CR0
	  RETURN
	END IF
	IF(IABS(DM3).GT.DJ3) THEN
	  PRINT *,'  warning(',NAME,'):  ABS(M3) > J3'
	  PRINT *,'  J3 =  ',J3
	  PRINT *,'  M3 =  ',M3
          S3J=CR0
	  RETURN
	END IF
C
	IF(MOD(DJ1+DJ2+DJ3,2).NE.0) THEN
	  PRINT *,'  warning(',NAME,'):  J1+J2+J3 is not integer'
	  PRINT *,'  J1 =  ',J1
	  PRINT *,'  J2 =  ',J2
	  PRINT *,'  J3 =  ',J3
          S3J=CR0
	  RETURN
	END IF
C
	IF(DM1+DM2+DM3.NE.0) THEN
          S3J=CR0
	  RETURN
	END IF
C
	IF(( DJ1+DJ2-DJ3).LT.0.OR.
     &     ( DJ1-DJ2+DJ3).LT.0.OR.
     &     (-DJ1+DJ2+DJ3).LT.0) THEN
          S3J=CR0
	  RETURN
	END IF
C
	WK(1)=(DJ1+DJ2-DJ3)/2
	WK(2)=(DJ1-DJ2+DJ3)/2
	WK(3)=(-DJ1+DJ2+DJ3)/2
	WK(4)=(DJ1+DJ2+DJ3)/2+1
	WK(5)=(DJ1+DM1)/2
	WK(6)=(DJ1-DM1)/2
	WK(7)=(DJ2+DM2)/2
	WK(8)=(DJ2-DM2)/2
	WK(9)=(DJ3+DM3)/2
	WK(10)=(DJ3-DM3)/2
	WK(11)=(-DJ1-DM2+DJ3)/2
	WK(12)=(DM1-DJ2+DJ3)/2
	WK(13)=(DJ1-DJ2-DM3)/2
C
	NLF=MAX0(WK(4),DJ1,DJ2,DJ3)
	IF(NLF.GT.NLF_M) THEN
	  PRINT *,'  error(',NAME,'):  NLF > NLF_M'
	  PRINT *,'  NLF =  ',NLF
	  PRINT *,'  NLF_M =  ',NLF_M
	  STOP
	END IF
C
	LCF=0.5D0*(LF(WK(1))+LF(WK(2))+LF(WK(3))-LF(WK(4))+
     &	   LF(WK(5))+LF(WK(6))+LF(WK(7))+LF(WK(8))+LF(WK(9))+LF(WK(10)))
C
	ZMIN=MAX0(0,-WK(11),-WK(12))
	ZMAX=MIN0(WK(1),WK(6),WK(7))
C
	IF(MOD(WK(13)+ZMIN,2).EQ.0) THEN
	  SGN=-CR1
	ELSE
	  SGN=CR1
	END IF
        S3J=CR0
	DO Z=ZMIN,ZMAX
	  SGN=-SGN
          S3J=S3J+SGN*EXP(LCF-LF(Z)-LF(WK(11)+Z)-LF(WK(12)+Z)-
     &                LF(WK(1)-Z)-LF(WK(6)-Z)-LF(WK(7)-Z))
        END DO
C
	END
C
C########################################################################
C
      FUNCTION S6J(J1,J2,J3,J4,J5,J6)
      USE MYDATA,  ONLY: KIND, CR0, CR0D5, CR1, CR2
C
C-----------------------------------------------------------------------
C...program name
        CHARACTER*5 NAME
        PARAMETER (NAME='S6J')
C...parameters
	INTEGER*4 L_M,NLF_M
	PARAMETER (L_M=100,NLF_M=4*L_M+1)
C...input scalar arguments
	REAL(KIND) J1,J2,J3,J4,J5,J6
C...output function
      REAL(KIND) S6J
C...save variables
	LOGICAL*4 INIT
	REAL(KIND) LF(0:NLF_M)
	SAVE INIT,LF
C...local scalars
	INTEGER*4 DJ1,DJ2,DJ3,DJ4,DJ5,DJ6,ILF,NLF,Z,ZMIN,ZMAX
	REAL(KIND) LCF,SGN
C...local arrays
	INTEGER*4 WK(23)
C...constants
	DATA INIT/.TRUE./
C
C...executable statements
C
	IF(INIT) THEN
	  LF(0)=CR0
	  LF(1)=CR0
	  DO ILF=2,NLF_M
	    LF(ILF)=LF(ILF-1)+LOG(CR1*ILF)
	  END DO
	  INIT=.FALSE.
	END IF
C
	IF(J1.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J1 < 0'
	  PRINT *,'  J1 =  ',J1
	  STOP
	END IF
	IF(J2.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J2 < 0'
	  PRINT *,'  J2 =  ',J2
	  STOP
	END IF
	IF(J3.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J3 < 0'
	  PRINT *,'  J3 =  ',J3
	  STOP
	END IF
	IF(J4.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J4 < 0'
	  PRINT *,'  J4 =  ',J4
	  STOP
	END IF
	IF(J5.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J5 < 0'
	  PRINT *,'  J5 =  ',J5
	  STOP
	END IF
	IF(J6.LT.CR0) THEN
	  PRINT *,'  error(',NAME,'):  J6 < 0'
	  PRINT *,'  J6 =  ',J6
	  STOP
	END IF
C
	IF(CR1*INT(J1).EQ.J1) THEN
	  DJ1=2*INT(J1)
	ELSE
	  IF(CR1*INT(J1)+0.5D0.EQ.J1) THEN
	    DJ1=2*INT(J1)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J1 is wrong'
	    PRINT *,'  J1 =  ',J1
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J2).EQ.J2) THEN
	  DJ2=2*INT(J2)
	ELSE
	  IF(CR1*INT(J2)+CR0D5.EQ.J2) THEN
	    DJ2=2*INT(J2)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J2 is wrong'
	    PRINT *,'  J2 =  ',J2
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J3).EQ.J3) THEN
	  DJ3=2*INT(J3)
	ELSE
	  IF(CR1*INT(J3)+CR0D5.EQ.J3) THEN
	    DJ3=2*INT(J3)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J3 is wrong'
	    PRINT *,'  J3 =  ',J3
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J4).EQ.J4) THEN
	  DJ4=2*INT(J4)
	ELSE
	  IF(CR1*INT(J4)+CR0D5.EQ.J4) THEN
	    DJ4=2*INT(J4)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J4 is wrong'
	    PRINT *,'  J4 =  ',J4
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J5).EQ.J5) THEN
	  DJ5=2*INT(J5)
	ELSE
	  IF(CR1*INT(J5)+CR0D5.EQ.J5) THEN
	    DJ5=2*INT(J5)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J5 is wrong'
	    PRINT *,'  J5 =  ',J5
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J6).EQ.J6) THEN
	  DJ6=2*INT(J6)
	ELSE
	  IF(CR1*INT(J6)+CR0D5.EQ.J6) THEN
	    DJ6=2*INT(J6)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J6 is wrong'
	    PRINT *,'  J6 =  ',J6
	    STOP
	  END IF
	END IF
C
	IF(MOD(DJ1+DJ2+DJ3,2).NE.0) THEN
	  PRINT *,'  warning(',NAME,'):  J1+J2+J3 is not integer'
	  PRINT *,'  J1 =  ',J1
	  PRINT *,'  J2 =  ',J2
	  PRINT *,'  J3 =  ',J3
          S6J=CR0
	  RETURN
	END IF
	IF(MOD(DJ1+DJ5+DJ6,2).NE.0) THEN
	  PRINT *,'  warning(',NAME,'):  J1+J5+J6 is not integer'
	  PRINT *,'  J1 =  ',J1
	  PRINT *,'  J5 =  ',J5
	  PRINT *,'  J6 =  ',J6
          S6J=CR0
	  RETURN
	END IF
	IF(MOD(DJ4+DJ2+DJ6,2).NE.0) THEN
	  PRINT *,'  warning(',NAME,'):  J4+J2+J6 is not integer'
	  PRINT *,'  J4 =  ',J4
	  PRINT *,'  J2 =  ',J2
	  PRINT *,'  J6 =  ',J6
          S6J=CR0
	  RETURN
	END IF
	IF(MOD(DJ4+DJ5+DJ3,2).NE.0) THEN
	  PRINT *,'  warning(',NAME,'):  J4+J5+J3 is not integer'
	  PRINT *,'  J4 =  ',J4
	  PRINT *,'  J5 =  ',J5
	  PRINT *,'  J3 =  ',J3
          S6J=CR0
	  RETURN
	END IF
C
	IF(( DJ1+DJ2-DJ3).LT.0.OR.
     &     ( DJ1-DJ2+DJ3).LT.0.OR.
     &     (-DJ1+DJ2+DJ3).LT.0.OR.
     &     ( DJ1+DJ5-DJ6).LT.0.OR.
     &     ( DJ1-DJ5+DJ6).LT.0.OR.
     &     (-DJ1+DJ5+DJ6).LT.0.OR.
     &     ( DJ4+DJ2-DJ6).LT.0.OR.
     &     ( DJ4-DJ2+DJ6).LT.0.OR.
     &     (-DJ4+DJ2+DJ6).LT.0.OR.
     &     ( DJ4+DJ5-DJ3).LT.0.OR.
     &     ( DJ4-DJ5+DJ3).LT.0.OR.
     &     (-DJ4+DJ5+DJ3).LT.0) THEN
          S6J=CR0
	  RETURN
	END IF
C
	WK(1)=(DJ1+DJ2-DJ3)/2
	WK(2)=(DJ1-DJ2+DJ3)/2
	WK(3)=(-DJ1+DJ2+DJ3)/2
	WK(4)=(DJ1+DJ2+DJ3)/2+1
	WK(5)=(DJ1+DJ5-DJ6)/2
	WK(6)=(DJ1-DJ5+DJ6)/2
	WK(7)=(-DJ1+DJ5+DJ6)/2
	WK(8)=(DJ1+DJ5+DJ6)/2+1
	WK(9)=(DJ4+DJ2-DJ6)/2
	WK(10)=(DJ4-DJ2+DJ6)/2
	WK(11)=(-DJ4+DJ2+DJ6)/2
	WK(12)=(DJ4+DJ2+DJ6)/2+1
	WK(13)=(DJ4+DJ5-DJ3)/2
	WK(14)=(DJ4-DJ5+DJ3)/2
	WK(15)=(-DJ4+DJ5+DJ3)/2
	WK(16)=(DJ4+DJ5+DJ3)/2+1
	WK(17)=(DJ1+DJ2+DJ4+DJ5)/2
	WK(18)=(DJ2+DJ3+DJ5+DJ6)/2
	WK(19)=(DJ3+DJ1+DJ6+DJ4)/2
	WK(20)=(DJ1+DJ2+DJ3)/2
	WK(21)=(DJ1+DJ5+DJ6)/2
	WK(22)=(DJ4+DJ2+DJ6)/2
	WK(23)=(DJ4+DJ5+DJ3)/2
C
	ZMIN=MAX0(WK(20),WK(21),WK(22),WK(23))
	ZMAX=MIN0(WK(17),WK(18),WK(19))
C
	NLF=ZMAX+1
	IF(NLF.GT.NLF_M) THEN
	  PRINT *,'  error(',NAME,'):  NLF > NLF_M'
	  PRINT *,'  NLF =  ',NLF
	  PRINT *,'  NLF_M =  ',NLF_M
	  STOP
	END IF
C
	LCF=CR0D5*(LF(WK(1))+LF(WK(2))+LF(WK(3))-LF(WK(4))+
     &             LF(WK(5))+LF(WK(6))+LF(WK(7))-LF(WK(8))+
     &             LF(WK(9))+LF(WK(10))+LF(WK(11))-LF(WK(12))+
     &             LF(WK(13))+LF(WK(14))+LF(WK(15))-LF(WK(16)))
C
	IF(MOD(ZMIN,2).EQ.0) THEN
	  SGN=-CR1
	ELSE
	  SGN=CR1
	END IF
        S6J=CR0
	DO Z=ZMIN,ZMAX
	  SGN=-SGN
          S6J=S6J+SGN*EXP(LCF+LF(Z+1)-
     &          LF(WK(17)-Z)-LF(WK(18)-Z)-LF(WK(19)-Z)-
     &          LF(Z-WK(20))-LF(Z-WK(21))-LF(Z-WK(22))-LF(Z-WK(23)))
      END DO
C
	END
C
C#######################################################################
C
      FUNCTION S9J(J11,J12,J13,J21,J22,J23,J31,J32,J33)
      USE MYDATA,  ONLY: KIND, CR0, CR0D5, CR1, CR2
C
C-----------------------------------------------------------------------
C...program name
      CHARACTER*5 NAME
      PARAMETER (NAME='S9J')
C...input scalar arguments
	REAL(KIND) J11,J12,J13,J21,J22,J23,J31,J32,J33
C...output function
        REAL(KIND) S9J
C...local scalars
	INTEGER*4 DJ11,DJ12,DJ21,DJ23,DJ32,DJ33,DJ,DJMIN,DJMAX
	REAL(KIND) J
C...external functions
      REAL(KIND) S6J
      EXTERNAL S6J
C
C...executable statements
C
	IF(CR1*INT(J11).EQ.J11) THEN
	  DJ11=2*INT(J11)
	ELSE
	  IF(CR1*INT(J11)+CR0D5.EQ.J11) THEN
	    DJ11=2*INT(J11)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J11 is wrong'
	    PRINT *,'  J11 =  ',J11
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J12).EQ.J12) THEN
	  DJ12=2*INT(J12)
	ELSE
	  IF(CR1*INT(J12)+CR0D5.EQ.J12) THEN
	    DJ12=2*INT(J12)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J12 is wrong'
	    PRINT *,'  J12 =  ',J12
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J21).EQ.J21) THEN
	  DJ21=2*INT(J21)
	ELSE
	  IF(CR1*INT(J21)+CR0D5.EQ.J21) THEN
	    DJ21=2*INT(J21)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J21 is wrong'
	    PRINT *,'  J21 =  ',J21
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J23).EQ.J23) THEN
	  DJ23=2*INT(J23)
	ELSE
	  IF(CR1*INT(J23)+CR0D5.EQ.J23) THEN
	    DJ23=2*INT(J23)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J23 is wrong'
	    PRINT *,'  J23 =  ',J23
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J32).EQ.J32) THEN
	  DJ32=2*INT(J32)
	ELSE
	  IF(CR1*INT(J32)+CR0D5.EQ.J32) THEN
	    DJ32=2*INT(J32)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J32 is wrong'
	    PRINT *,'  J32 =  ',J32
	    STOP
	  END IF
	END IF
	IF(CR1*INT(J33).EQ.J33) THEN
	  DJ33=2*INT(J33)
	ELSE
	  IF(CR1*INT(J33)+CR0D5.EQ.J33) THEN
	    DJ33=2*INT(J33)+1
	  ELSE
	    PRINT *,'  error(',NAME,'):  J33 is wrong'
	    PRINT *,'  J33 =  ',J33
	    STOP
	  END IF
	END IF
C
	DJMIN=MAX0(IABS(DJ11-DJ33),IABS(DJ12-DJ23),IABS(DJ21-DJ32))
	DJMAX=MIN0(DJ11+DJ33,DJ12+DJ23,DJ21+DJ32)
C
        S9J=CR0
	DO DJ=DJMIN,DJMAX,2
	  J=CR0D5*DJ
          S9J=S9J+(CR2*J+CR1)*S6J(J11,J21,J31,J32,J33,J)*
     &                              S6J(J12,J22,J32,J21,J,J23)*
     &                              S6J(J13,J23,J33,J,J11,J12)
	END DO
        IF(MOD(DJMIN,2).NE.0) S9J=-S9J
C
	END
C
C####################################################################
C
      FUNCTION SFEFUN(LURA,MURA,TETURA,FIURA)
C
C **************************************************************
C ***                                                        ***
C ***    IT  GIVES  Re( Y         (TETA,FIURA)  )            ***
C ***                    LURA,MURA                           ***
C ***                                                        ***
C **************************************************************
C
      USE MYDATA,  ONLY: KIND, CHEPI, CR0, CR1, CR2
      IMPLICIT REAL(KIND) (A-F,H,O-Z)
      COMMON/NJ/IT,G(302)
      DIMENSION I(4)
      EQUIVALENCE (I(1),I1),(I(2),I2),(I(3),I3),(I(4),I4)
C
      COE = CR1
      COE = H(CR1*LURA)*SQRT(COE/CHEPI)
C
      AJ=LURA
      AM1=IABS(MURA)
      AM2=CR0
      ANGL=TETURA
C **************************************************************
C ***  CALCULATION OF D-FUNCTIONS (GRATE)                    ***
C ***                                                        ***
C ***       NOTE: GRATE - BUT WITH FI = 0 OR 180 DEGREES     ***
C **************************************************************
      ROLL=CR0
      IF(IT-20431) 10,30,10
 10   IT=20431
      G(1)=CR0
      G(2)=CR0
      DO 20 J=3,302
      X=J-1
      GXXXX=X
 20   G(J)=G(J-1)+LOG(GXXXX)
 30   J1=CR2*AJ+0.1D0
      M1=CR2*AM1+0.1D0
      M2=CR2*AM2+0.1D0
      I(1)=J1+M1
      I(2)=J1-M1
      I(3)=J1+M2
      I(4)=J1-M2
      DO 50 J=1,4
      K=I(J)/2
      IF(I(J)-2*K) 300,40,300
 40   IF(K) 300,50,50
 50   I(J)=K+1
      Y=ANGL/2.0D0
      X=COS(Y)
      Y=SIN(Y)
      IL=1
      IF(I3-I1) 150,150,140
 140  IL=I3-I1+1
 150  IN=I2
      IF(I3-I2) 160,170,170
 160  IN=I3
 170  C=(G(I1)+G(I2)+G(I3)+G(I4))/CR2
      S=CR0
      ST=CR1
      DO 180 J=IL,IN
      J1=I2-J+1
      J2=I3-J+1
      J3=J
      J4=I1-I3+J
      J5=I2+I3-2*J
      J6=I1-I3+2*(J-1)
      A1=CR1
      IF(J5) 300,173,171
 171  DO 172 L=1,J5
 172  A1=A1*X
 173  A2=CR1
      IF(J6) 300,176,174
 174  DO 175 L=1,J6
 175  A2=A2*Y
 176  CONTINUE
      S=S+ST*A1*A2*EXP(C-G(J1)-G(J2)-G(J3)-G(J4))
 180  ST=-ST
      K=(IL-1)/2+0.1D0
      IF(IL-1-2*K) 190,200,190
 190  S=-S
 200  ROLL=S
 300  CONTINUE
      LEMURA=IABS(MURA)
      AKOT=CR1
      IF(MURA.GT.0) AKOT=(-CR1)**LEMURA
      SFEFUN=ROLL*COS(FIURA*MURA)*AKOT*COE
      RETURN
      END
C
C##############################################################################
C#     MODYFIED P_L^M FUNCTION: P_L^M(X)/(2M-1)!!   NORMALIZED FOR           ##
C#     INTGRATION OVER d(COS(TETA)) SPHERICAL FUNCTION IS OBTAINED           ##
C#     Y(L,M) = (-1)**M * YNORM(L,IABS(M)) * PLM(L,IABS(M),COS(TETA))        ##
C##############################################################################

      FUNCTION PLM(L,M,X)
      USE MYDATA,  ONLY: KIND, CR1, CR2

      IMPLICIT REAL(KIND) (A-H,O-Z)     
      INTEGER, INTENT(IN) :: L,M
      REAL(KIND), INTENT(IN) :: X
      Y=-SQRT(CR1-X*X) 
      A=Y**M
      IF(L.EQ.M)THEN
        PLM=A
        RETURN
      ENDIF
      B=X*(CR2*M+CR1)*A
      IF(L.EQ.M+1)THEN
        PLM=B
        RETURN
      END IF
      DO I=M+2,L
        PLM=(X*(CR2*I-CR1)*B-(I+M-CR1)*A)/(I-M)
        A=B
        B=PLM
      END DO
      RETURN
      END FUNCTION PLM

