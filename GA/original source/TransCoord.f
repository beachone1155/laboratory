      PARAMETER (M=1000)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HIGH(M),HOK(M),TOK(M)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND
      CHARACTER KIND*1

      PI=3.1415926D0
      RAD=180./PI

      KIND="W"  ! World Coordinate
      IND=2     ! Degree Format (not Do-Fun-Byou Format)

c     Origin of the coordinate (Do not change!)
      OLAT1=36.000000   ! Do (Degree)
      OLAT2=0.0000000   ! Fun (Minute)
      OLAT3=0.0000000   ! Byou (Second)
      OLON =139.83333   ! Degree Format Only

c     Data Info.
      N=1               ! Number of Data Point
      HIGH(1)=0.        ! Height (m)
      HOK(1) =37.000000 ! Latitude (Degree)
      TOK(1) =139.83333 ! Longitude (Degree)

      FA1=OLAT1
      FA2=OLAT2
      FA3=OLAT3
      CALL CONST(KIND)
      CALL SHIGO(S0)

      DO 100 I=1,N
      IF(IND.NE.1) GOTO 10
      FA1=REAL(INT(HOK(I)/10000.))
      FA2=REAL(INT((HOK(I)-FA1*10000.)/100.))
      FA3=( HOK(I)-FA1*10000.-FA2*100.)
      R1 =REAL(INT(TOK(I)/10000.))
      R2 =REAL(INT((TOK(I)- R1*10000.)/100.))
      R3 =( TOK(I)- R1*10000.- R2*100.)
      GOTO 20
   10 FA1=HOK(I)
      R1 =TOK(I)

c     Transformation: (Latitude, Longitude) -> (X, Y)
   20 CALL SHIGO(S)
      CALL COMPA(OLON)
      CALL ZAHYO(S,S0)
      WRITE(*,'(4A,8A,8A)') 'HIGH','   X(+N)','   Y(+E)'
      WRITE(*,'(I3,2(A,F8.0))') INT(HIGH(I)),',',X,',',Y

c     Transformation: (X, Y) -> (Latitude, Longitude)     
      IF(IND.EQ.1) FA=(OLAT1*3600.+OLAT2*60.+OLAT3)/(3600.*RAD)
      IF(IND.EQ.2) FA=OLAT1/RAD
      CALL DEGFA(S,S0,FA)
      CALL COOD(FA,HOK_I,TOK_I)
      HOK_I=HOK_I*RAD
      TOK_I=TOK_I*RAD
      WRITE(*,'(5A,8A,8A)') 'HIGH ','Lat.(N) ','Lon.(E) '
      WRITE(*,'(I3,A,F8.5,A,F9.5)') INT(HIGH(I)),',',HOK_I,',',TOK_I

  100 CONTINUE
 
      STOP
      END


      SUBROUTINE CONST(KIND)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND
      CHARACTER KIND*1
c     ---------------------------------- World/Japan Coordinate
      IF(KIND.EQ.'J'.OR.KIND.EQ.'j') THEN
      A =6.377397155D6
      F =2.99152813D2
      ELSE
      A =6.378137D6
      F =2.98257222101D2
      ENDIF

      E =SQRT(2.*F-1.)/F
      ED=SQRT(2.*F-1.)/(F-1.)

      RETURN
      END


      SUBROUTINE SHIGO (S)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND

      AA=1.+3.*E**2/4.+45.*E**4/64.+175.*E**6/256.+11025.*E**8/16384.
     *  +43659.*E**10/65536.+693693.*E**12/1048576.
     *  +19324305.*E**14/29360128.+4927697775.*E**16/7516192768.
      BB=3.*E**2/4.+15.*E**4/16.+525.*E**6/512.+2205.*E**8/2048.
     *  +72765.*E**10/65536.+297297.*E**12/262144.
     *  +135270135.*E**14/117440512.+547521975.*E**16/469762048.
      CC=15.*E**4/64.+105.*E**6/256.+2205.*E**8/4096.
     *  +10395.*E**10/16384.+1486485.*E**12/2097152.
     *  +45090045.*E**14/58720256.+766530765.*E**16/939524096.
      DD=35.*E**6/512.+315.*E**8/2048.+31185.*E**10/131072.
     *  +165165.*E**12/524288.+45090045.*E**14/117440512.
     *  +209053845.*E**16/469762048.
      EE=315.*E**8/16384.+3465.*E**10/65536.+99099.*E**12/1048576.
     *  +4099095.*E**14/29360128.+348423075.*E**16/1879048192.
      FF=693.*E**10/131072.+9009.*E**12/524288.
     *  +4099095.*E**14/117440512.+26801775.*E**16/469762048.
      GG=3003.*E**12/2097152.+315315.*E**14/58720256.
     *  +11486475.*E**16/939524096.
      HH=45045.*E**14/117440512.+765765.*E**16/469762048.
      RI=765765.*E**16/7516192768.

      B1=A*(1.-E**2)*AA
      B2=A*(1.-E**2)*(-BB/2.)
      B3=A*(1.-E**2)*(CC/4.)
      B4=A*(1.-E**2)*(-DD/6.)
      B5=A*(1.-E**2)*(EE/8.)
      B6=A*(1.-E**2)*(-FF/10.)
      B7=A*(1.-E**2)*(GG/12.)
      B8=A*(1.-E**2)*(-HH/14.)
      B9=A*(1.-E**2)*(RI/16.)

      IF(IND.EQ.1) FAI=(FA1*3600.+FA2*60.+FA3)/(3600.*RAD)
      IF(IND.EQ.2) FAI=FA1/RAD

      S02=SIN(FAI* 2.)
      S04=SIN(FAI* 4.)
      S06=SIN(FAI* 6.)
      S08=SIN(FAI* 8.)
      S10=SIN(FAI*10.)
      S12=SIN(FAI*12.)
      S14=SIN(FAI*14.)
      S16=SIN(FAI*16.)
      S=B1*FAI+B2*S02+B3*S04+B4*S06+B5*S08+B6*S10+B7*S12+B8*S14+B9*S16

      RETURN
      END


      SUBROUTINE SHIGO2 (S,FAI)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND

      AA=1.+3.*E**2/4.+45.*E**4/64.+175.*E**6/256.+11025.*E**8/16384.
     *  +43659.*E**10/65536.+693693.*E**12/1048576.
     *  +19324305.*E**14/29360128.+4927697775.*E**16/7516192768.
      BB=3.*E**2/4.+15.*E**4/16.+525.*E**6/512.+2205.*E**8/2048.
     *  +72765.*E**10/65536.+297297.*E**12/262144.
     *  +135270135.*E**14/117440512.+547521975.*E**16/469762048.
      CC=15.*E**4/64.+105.*E**6/256.+2205.*E**8/4096.
     *  +10395.*E**10/16384.+1486485.*E**12/2097152.
     *  +45090045.*E**14/58720256.+766530765.*E**16/939524096.
      DD=35.*E**6/512.+315.*E**8/2048.+31185.*E**10/131072.
     *  +165165.*E**12/524288.+45090045.*E**14/117440512.
     *  +209053845.*E**16/469762048.
      EE=315.*E**8/16384.+3465.*E**10/65536.+99099.*E**12/1048576.
     *  +4099095.*E**14/29360128.+348423075.*E**16/1879048192.
      FF=693.*E**10/131072.+9009.*E**12/524288.
     *  +4099095.*E**14/117440512.+26801775.*E**16/469762048.
      GG=3003.*E**12/2097152.+315315.*E**14/58720256.
     *  +11486475.*E**16/939524096.
      HH=45045.*E**14/117440512.+765765.*E**16/469762048.
      RI=765765.*E**16/7516192768.

      B1=A*(1.-E**2)*AA
      B2=A*(1.-E**2)*(-BB/2.)
      B3=A*(1.-E**2)*(CC/4.)
      B4=A*(1.-E**2)*(-DD/6.)
      B5=A*(1.-E**2)*(EE/8.)
      B6=A*(1.-E**2)*(-FF/10.)
      B7=A*(1.-E**2)*(GG/12.)
      B8=A*(1.-E**2)*(-HH/14.)
      B9=A*(1.-E**2)*(RI/16.)

      S02=SIN(FAI* 2.)
      S04=SIN(FAI* 4.)
      S06=SIN(FAI* 6.)
      S08=SIN(FAI* 8.)
      S10=SIN(FAI*10.)
      S12=SIN(FAI*12.)
      S14=SIN(FAI*14.)
      S16=SIN(FAI*16.)
      S=B1*FAI+B2*S02+B3*S04+B4*S06+B5*S08+B6*S10+B7*S12+B8*S14+B9*S16

      RETURN
      END


      SUBROUTINE COMPA(OLON)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND

      OLON=139.83333
      RAM=OLON/RAD
      IF(IND.EQ.1) DR=((R1*3600.+ R2*60.+ R3)/3600.-OLON)/RAD
      IF(IND.EQ.2) DR=(R1-OLON)/RAD

      IF(IND.EQ.1) FAI=(FA1*3600.+FA2*60.+FA3)/(3600.*RAD)
      IF(IND.EQ.2) FAI=FA1/RAD
      CS=COS(FAI)
      TA=TAN(FAI)

      V=SQRT(1.+(ED*CS)**2)
      C=A*F/(F-1.)
      RN=C/V
      ET=ED*CS
 
      RETURN
      END


      SUBROUTINE COMPA2(FAI)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND

      SS=SIN(FAI)
      CS=COS(FAI)
      TA=TAN(FAI)

      V=SQRT(1.+(ED*CS)**2)
      C=A*F/(F-1.)
      RN=C/V
      ET=ED*CS

      RETURN
      END


      SUBROUTINE ZAHYO(S,S0)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND
      DATA RM0/0.9999/

      E2=ET**2
      T2=TA**2
      T4=TA**4
      T6=TA**6

      X=((S-S0)+RN*CS**2*TA*DR**2/2.
     * +RN*CS**4*TA*(5.-T2+9.*E2+4.*E2**2)*DR**4/24.
     * -RN*CS**6*TA*(-61.+58.*T2-T4-270.*E2+330.*T2*E2)*DR**6/720.
     * -RN*CS**8*TA*(-1385.+3111.*T2-543.*T4+T6)*DR**8/40320.)*RM0

      Y=(RN*CS*DR-RN*CS**3*(-1.+T2-E2)*DR**3/6.
     * -RN*CS**5*(-5.+18.*T2-T4-14.*E2+58.*T2*E2)*DR**5/120.
     * -RN*CS**7*(-61.+479.*T2-179.*T4+T6)*DR**7/5040.)*RM0

      RETURN
      END


      SUBROUTINE DEGFA(S,S0,FAI)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND
      DATA RM0/0.9999/

      RM=S0+X/RM0
      PFA=FAI
c     ---------------------------------- Repeate eps(10^-5)
   10 CONTINUE

      CALL COMPA2(PFA)
      CALL SHIGO2(S,PFA)
c     ---------------------------------- Calc. FAI 
      FAI=PFA+2.*(S-RM)*(1.-(E*SS)**2)**1.5
     &   /(3.*E**2*(S-RM)*SS*CS*(1.-(E*SS)**2)**0.5-2.*A*(1.-E**2))
      SFA=FAI*3600.*RAD
      SPF=PFA*3600.*RAD
      IF(ABS(SFA-SPF).LT.2.D-5) RETURN
      PFA=FAI
      GOTO 10

      END


      SUBROUTINE COOD(FAI,HOK,TOK)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CST/FA1,FA2,FA3,R1,R2,R3,DR,CS,TA,RN,
     &        X,Y,RAM,SS,ET,E,ED,A,F,PI,RAD,IND
      DATA RM0/0.9999/

      E2=ET**2
      T2=TA**2
      T4=TA**4
      T6=TA**6
      YM=Y/RM0

      HOK=FAI-0.5*TA*(1.+E2)*YM**2/RN**2
     *   +TA*(5.+3*T2+6.*E2-6.*T2*E2-3.*E2**2-9.*T2*E2**2)*YM**4
     *   /(24.*RN**4)
     *   -TA*(61.+90.*T2+45.*T4+107.*E2-162.*T2*E2-45.*T4*E2)*YM**6
     *   /(720.*RN**6)
     *   +TA*(1385.+3633.*T2+4095.*T4+1575*T6)*YM**8/(40320.*RN**8)

      DR=YM/(RN*CS)-(1.+2.*T2+E2)*YM**3/(6.*RN**3*CS)
     *  +(5.+28.*T2+24.*T4+6.*E2+8.*T2*E2)*YM**5/(120.*RN**5*CS)
     *  -(61.+662.*T2+1320.*T4+720.*T6)*YM**7/(5040.*RN**7*CS)
      TOK=RAM+DR

      RETURN
      END

