c   --------------------------------------------------------------
c                          SUBROUTINE  EOSI
c   --------------------------------------------------------------
c
c              ECUACION DE ESTADO IONICA
C***************************************************************

      subroutine eosi(ro,t,mui,ye,p,dpdt,dpdr,e,dedt,dedr,xn,a,z,n)
C
C  EOS OF THE IONS (COMPLETELY IONIZED, EXCLUDING NUCLEAR INTERACTIONS)
C    Bravo&Garcia 1992
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MUE,MUI,LAMBDA
      
c      PARAMETER (n=6)
      
      DIMENSION z(n),a(n),AZ(30),BZ(30),CZ(30),x(n),xn(n)
c      common/xnuc2/xn(n)

      DATA R,UMA/8.31451D7,1.660540D-24/
      DATA A1,B1,C1,D1,E1/-0.898004D0,0.96786D0,0.220703D0,-0.86097D0,
     & 2.5269D0/
      DATA BET,GAM/0.29561D0,1.9885D0/
      DATA AZ/  0.01900d0,  0.17563d0,  0.55068d0,  1.19420d0,2.15100d0,
     &3.46289d0,5.16330d0,7.29083d0,  9.86272d0, 12.91897d0, 16.46465d0,
     &20.53046d0,25.14461d0,30.30823d0,36.05356d0,42.37249d0,49.32911d0,
     &56.91062d0,65.12891d0,73.98992d0,83.55745d0,93.80458d0,
     &104.73189d0,116.37763d0,128.75405d0,141.88308d0,155.78658d0,
     &170.41906d0,185.83520d0,202.04560d0/
      DATA BZ/  0.30480d0,0.59284d0,  0.91871d0,  1.27392d0,1.64098d0,
     &2.01314d0,2.38638d0,2.76091d0,  3.13770d0,  3.51750d0,  3.89974d0,
     &4.28643d0,4.67456d0,5.06763d0,  5.46386d0,  5.86207d0,  6.26504d0,
     &6.66890d0,7.07476d0,7.48181d0,  7.89211d0,  8.30193d0,  8.71471d0,
     &9.12949d0,9.54329d0,9.96024d0, 10.37804d0, 10.79736d0, 11.21502d0,
     &11.64054d0/
      DATA CZ/ -0.10403d0,-0.12354d0,-0.15982d0, -0.20283d0, -0.24347d0,
     &-0.28024d0,-0.31301d0,-0.34274d0,-0.37057d0,-0.39686d0,-0.42249d0,
     &-0.44745d0,-0.47209d0,-0.49634d0,-0.52022d0,-0.54365d0,-0.56705d0,
     &-0.58976d0,-0.61214d0,-0.63407d0,-0.65578d0,-0.67722d0,-0.69802d0,
     &-0.71876d0,-0.73889d0,-0.75865d0,-0.77830d0,-0.79757d0,-0.81668d0,
     &-0.83569d0/

c      if(n.ne.nP)then
c         write(*,*) 'WARNING!! Repasar n (niso) en eosi'
c         write(*,*) '     eosi','     Main'
c         write(*,*)
c         write(*,*) 'n (niso)',n,nP
c         stop
c      endif


C  CHEMICAL COMPOSITION MEANS
C
      Z1=0.D0
      Z2=0.D0
      Z53=0.D0
      Z512=0.D0
      ZM512=0.D0
      ZL=0.D0
      ZA=0.D0
      ZB=0.D0
      ZC=0.D0
c      MUI=0.D0

c      do i=1,n
c         xn(i)=0.d0
c      enddo
c      xn(2)=.5d0
c      xn(3)=.5d0

c      do 3 j=1,n
c         MUI=MUI+X(J)/A(J)
c 3    continue
c      MUI=1.D0/MUI


C.. Pasa a fracciones en numero (en lugar de en masa) pq son las que se
C.. necesitan para las correcciones coulombianas.
      do 4 j=1,n
         x(j)=xn(j)/a(j)*mui
 4    continue
      DO 1 J=1,N
         IF(Z(J).EQ.0.) GOTO 1
         Z1=Z1+X(J)*Z(J)
         Z2=Z2+X(J)*Z(J)*Z(J)
         Z53=Z53+X(J)*Z(J)**(5.d0/3.d0)
         Z512=Z512+X(J)*Z(J)**(5.d0/12.d0)
         ZM512=ZM512+X(J)/Z(J)**(5.d0/12.d0)
         ZL=ZL+X(J)*DLOG(Z(J))
         ipos=idint(z(j))
         ZA=ZA+X(J)*AZ(ipos)
         ZB=ZB+X(J)*BZ(ipos)
         ZC=ZC+X(J)*CZ(ipos)
 1    CONTINUE
      RMUI=R/MUI
C     
C  CALCULATION OF THE PLASMA PARAMETERS
C
      MUE=1.D0/YE
      XX=1.0088D-2*(RO*YE)**(1.d0/3.d0)
      GAME=2.255347D7*XX/T
      ETA=3.423034D-2/XX*DSQRT(1.D0+XX*XX)
      GAMMA=GAME*Z53
      LAMBDA=7.833114D3*DSQRT(RO*YE*Z1/MUI)/T
C +++ 09/01/04
c      IF(LAMBDA.GT.500.) RETURN
      IF(LAMBDA.GT.500.d0) then
         write(*,*) 'Lambda>500'
         write(*,*) 'Too high density and/or too low Temperature'
         write(*,*) 'Density,Temp,Lambda'
c         write(*,*) ro,t,lambda
c         stop
      endif
C +++
C
C  SUBROUTINE FOR THE SOLID, AND CALCULATION OF THE QUANTUM CORRECTION
C
      CALL EOSIS(RO,T,MUE,Z1,MUI,LAMBDA,GAMMA,P,E,S,DPDT,DPDR,DEDT,DEDR,
     &     PQE,PQM,EQM,SQM,DPDTQM,DPDRQM,DEDTQM,DEDRQM,PQEQM,PT,ET,ST,
     &     DPDTT,DPDRT,DEDTT,DEDRT,PQET)
      IF(GAMMA.GE.180.D0) RETURN
C     
C  CALCULATION OF THE LIQUID EOS
C
C  IDEAL GAS TERMS
C
      PID=RMUI*RO*T
      EID=1.5D0*RMUI*T
      SI=0.D0
      DO 2 J=1,N
         IF(X(J).EQ.0.D0) GOTO 2
         ARI=3.120707D-4*A(J)**2.5d0/X(J)*T**1.5d0/RO
         SI=SI+X(J)/A(J)*DLOG(ARI)
    2 CONTINUE
      SID=2.5D0*RMUI*RO+R*RO*SI
      PQEID=0.D0
      DPDTID=RMUI*RO
      DPDRID=RMUI*T
      DEDTID=1.5D0*RMUI
      DEDRID=0.D0
C
C  COULOMB CORRECTION: 1) UNIFORM ELECTRON BACKGROUND
C
      IF(GAMMA.LT.1.D0) THEN
         GAM32=GAMMA**1.5d0
         GAMG=GAMMA**GAM
         PUB=-RMUI*RO*T*(0.288675D0*GAM32-BET/3.D0*GAMG)
         EUB=3.D0*PUB/RO
         SUB=-RMUI*RO*(0.288675D0*GAM32-BET*(GAM-1.D0)/GAM*GAMG)
         PQEUB=UMA*PUB/YE/RO
         DPDTUB=RMUI*RO*(0.144337D0*GAM32-BET*(GAM-1.D0)/3.D0*GAMG)
         DPDRUB=-RMUI*T*(0.4330127D0*GAM32-BET*(GAM+3.D0)/9.D0*GAMG)
         DEDTUB=3.D0*DPDTUB/RO
         DEDRUB=-RMUI*T/RO*(0.4330127D0*GAM32-BET*GAM/3.D0*GAMG)
      ELSE
         GAME0=GAME*Z53
         GAME14=GAME**.25d0*Z512
         GAMM14=ZM512*Z512/GAME14
         PUB=RMUI*RO*T/3.D0*(A1*GAME0+B1*GAME14+C1*GAMM14+D1)
         EUB=3.D0*PUB/RO
c         SUB=-RMUI*RO*(3.D0*B1*GAME14-5.D0*C1*GAMM14+D1*DLOG(GAME)+
c     &        1.6666667D0*D1*ZL-(D1+E1))
         fivethird=5.d0/3.d0
         SUB=-RMUI*RO*(3.D0*B1*GAME14-5.D0*C1*GAMM14+D1*DLOG(GAME)+
     &        fivethird*D1*ZL-(D1+E1))
         PQEUB=UMA*PUB/YE/RO
c         DPDTUB=RMUI*RO*(.25D0*B1*GAME14+0.41666667D0*C1*GAMM14+D1/3.D0)
         fivetwelve=5.d0/12.d0
         DPDTUB=RMUI*RO*(.25D0*B1*GAME14+fivetwelve*C1*GAMM14+D1/3.D0)
c         DPDRUB=RMUI*T*(0.4444444D0*A1*GAME0+0.36111111D0*B1*GAME14+
c     &        0.305555556D0*C1*GAMM14+D1/3.D0)
         fourninth=4.d0/9.d0
         thirteenthirtysix=13.d0/36.d0
         eleventhirtysix=11.d0/36.d0
         DPDRUB=RMUI*T*(fourninth*A1*GAME0+thirteenthirtysix*B1*GAME14+
     &        eleventhirtysix*C1*GAMM14+D1/3.D0)
         DEDTUB=3.D0*DPDTUB/RO
         DEDRUB=RMUI*T/RO*(A1/3.D0*GAME0+(B1*GAME14-C1*GAMM14)/12.D0)
      ENDIF
C
C  COULOMB CORRECTION: 2) POLARIZATION
C
      IF(GAMMA.LT.1.D0) THEN
         X2=XX*XX
         X12=DSQRT(1.D0+X2)
         FPOL=-3.953662D6*DSQRT(Z2/MUI*T*RO**3)*YE*X12/X2
         PPO=(X2-1.D0)/6.D0/(1.D0+X2)*FPOL
         EPO=.5D0*FPOL/RO
         SPO=-.5D0*FPOL/T
         PQEPO=(1.D0+2.D0*X2)/3.D0/(1.D0+X2)*UMA*FPOL/YE/RO
         DPDTPO=.5D0*PPO/T
         DPDRPO=(7.D0*X2*X2+6.D0*X2-5.D0)/36.D0/(1.D0+X2)/(1.D0+X2)
     *        *FPOL/RO
         DEDTPO=.25D0*FPOL/RO/T
         DEDRPO=.5D0*PPO/RO/RO
      ELSE
         X2=XX*XX
         X4=X2*X2
         X21=1.D0+X2
         GAME0=GAME*ZA
         GAME14=GAME**.25d0*ZB
         PPO=-RMUI*RO*T*ETA/3.D0/X21*(X2*GAME0+.25D0*(X2-3.D0)
     *        *GAME14-ZC)
         EPO=-RMUI*T*ETA*(GAME0+.25D0*GAME14)
         SPO=RMUI*RO*ETA*(.75D0*GAME14+ZC)
         PQEPO=UMA*PPO/RO/YE
         DPDTPO=RMUI*RO*ETA/X21*(.0625D0*(3.D0-X2)*GAME14+ZC/3.D0)
         DPDRPO=-RMUI*T*ETA/9.D0/X21/X21*((4.D0*X4+5.D0*X2)*GAME0
     *        +.0625D0*(13.D0*X4+2.D0*X2-27.D0)*GAME14-(X2+2.D0)*ZC)
         DEDTPO=-.1875D0*RMUI*ETA*GAME14
         DEDRPO=-RMUI/RO*T*ETA/X21*(X2/3.D0*GAME0+(X2-3.D0)/48.D0
     *        *GAME14)
      ENDIF
C
C  CALCULATION OF TOTAL MAGNITUDES
C
      P=PID+PUB+PPO+PQM
      E=EID+EUB+EPO+EQM
      S=SID+SUB+SPO+SQM
      PQE=PQEID+PQEUB+PQEPO+PQEQM
      DPDT=DPDTID+DPDTUB+DPDTPO+DPDTQM
      DPDR=DPDRID+DPDRUB+DPDRPO+DPDRQM
      DEDT=DEDTID+DEDTUB+DEDTPO+DEDTQM
      DEDR=DEDRID+DEDRUB+DEDRPO+DEDRQM
c       p=pid
c       dpdt=dpdtid
c       dpdr=dpdrid
c      e=eid
c      dedt=dedtid
c      dedr=dedrid
      RETURN
      END


C      ***************************************************************
C                   *************************
      SUBROUTINE EOSIS(RHO,T,XMUE,ZM,AM,LAMBDA,GAMMA,P,U,ENTRO,DPT,
     &     DPD,DUT,DUD,PQE,PQM,UQM,SQM,DPTQM,DPDQM,DUTQM,DUDQM,PQEQM,
     &     PT,UT,ST,DPTT,DPDT,DUTT,DUDT,PQET)
C  PROGRAMA QUE GENERA LAS FUNCIONES TERMODINAMICAS DE UN SOLIDO
C  BCC A TRAVES DE UNA APROXIMACION ANALITICA.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMBDA
      COMMON/COEF/X(3),XMAX(3)
      COMMON/FUN/D(6,3)
      DATA E,A/-1612.5d0,-0.895929d0/
      DATA R/8.31451D7/
      RMUI=R/AM
      X(1)=LAMBDA
      IF(X(1) .LE. 4.d0) THEN
         CALL FTER1 (T,RHO,ZM,AM,XMUE,ENTRO,U,P,DUT,DUD,DPT,DPD)
      ELSE
         XMAX(1)=X(1)*0.84957279d0
         XMAX(2)=0.30353d0*X(1)
         XMAX(3)=0.6756d0*X(1)
         CALL FTER2 (T,RHO,ZM,AM,XMUE,ENTRO,U,P,DUT,DUD,DPT,DPD)
      ENDIF
      IF(GAMMA.LT.180.D0) THEN
         UQM=U-3.D0*RMUI*T
         PQM=P-1.5D0*RMUI*RHO*T
         SQM=ENTRO+RMUI*(3.D0*DLOG(LAMBDA)-5.498183D0)
         DUTQM=DUT-3.D0*RMUI
         DUDQM=DUD
         DPTQM=DPT-1.5D0*RMUI*RHO
         DPDQM=DPD-1.5D0*RMUI*T
         PQEQM=0.D0
      ELSE
C     CORRECCION ENARMONICA EN EL SOLIDO
         D(1,1)=E/GAMMA/GAMMA
         D(2,1)=2.*D(1,1)/T
c         D(3,1)=-0.6666*D(1,1)/RHO
         twothird=2.d0/3.d0
         D(3,1)=-twothird*D(1,1)/RHO
         D(4,1)=D(2,1)/T
         D(5,1)=2.d0/T*D(3,1)
c         D(6,1)=-1.66667/RHO*D(3,1)
         fivethird=5.d0/3.d0
         D(6,1)=-fivethird/RHO*D(3,1)
         CALL FUNTER (T,RHO,1,AM,ENTROE,UE,PE,DUTE,DUDE,DPTE,DPDE)
C     PARTE COULOMBIANA
         D(1,1)=A*GAMMA
         D(2,1)=-D(1,1)/T
         D(3,1)=D(1,1)/3.d0/RHO
         D(4,1)=-2./T*D(2,1)
         D(5,1)=-D(3,1)/T
c         D(6,1)=-0.66667/RHO*D(3,1)
         D(6,1)=-twothird/RHO*D(3,1)
         CALL FUNTER (T,RHO,1,AM,ENTROC,UC,PC,DUTC,DUDC,DPTC,DPDC)
C     MAGNITUDES TOTALES
         U=U+UE+UC
         P=P+PE+PC
         ENTRO=ENTRO+ENTROE+ENTROC
         DUT=DUT+DUTE+DUTC
         DUD=DUD+DUDE+DUDC
         DPD=DPD+DPDE+DPDC
         DPT=DPT+DPTE+DPTC
         PQE=0.D0
      ENDIF
      RETURN
      END
C     ***************************************************************
C                   *************************
      SUBROUTINE FTER1(T,RHO,ZM,AM,XMUE,ENTRO,U,P,DUT,DUD,DPT,DPD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X1(12)
      COMMON/SCR/D11,D21,D31,D41,D51,D61,D111,D211,D311,D411,D511,D611
      COMMON/COEF/X(3),XMAX(3)
      COMMON/FUN/D(6,3)
      DATA AK1,AK2,AK4,AK6,AK8,AK10,AK12/-2.498183d0,4.166667D-02,
     $     -2.1128D-04,2.46426D-06,-3.6368D-08,5.98375D-10,1.050742D-11/
      X1(1)=X(1)
      DO 1 I=2,12
         X1(I)=X(1)*X1(I-1)
 1    ENDDO
      NMODES=1
      D(1,1)=3.*DLOG(X(1))+AK1+AK2*X1(2)+AK4*X1(4)+AK6*X1(6)+AK8*
     $     X1(8)+AK10*X1(10)+AK12*X1(12)
C     
      D1=3./X(1)+2.d0*AK2*X(1)+4.d0*AK4*X1(3)+6.d0*AK6*X1(5)+
     $     8.d0*AK8*X1(7)+10.d0*AK10*X1(9)+12.d0*AK12*X1(11)
C     
      D2=-3.d0/X1(2)+2.d0*AK2+12.d0*AK4*X1(2)+30.d0*AK6*X1(4)+
     $     56.d0*AK8*X1(6)+90.d0*AK10*X1(8)+132.d0*AK12*X1(10)
C     
C     CALCULA LAS DERIVADAS RESPECTO RHO,T
      D(2,1)=-D1*X(1)/T
      D(3,1)=D1*X(1)/RHO/2.d0
      D(4,1)=D2*X1(2)/T/T-D(2,1)*2.d0/T
      D(5,1)=-D2*X1(2)/RHO/T/2.d0+D(2,1)/RHO/2.d0
      D(6,1)=D2*X1(2)/4.d0/RHO/RHO-D(3,1)/2.d0/RHO
      IF(RHO.LE.1.D+05) THEN
C     CORRIGE DE SCREENING EL MODO LONGITUDINAL
         CALL SCREEN (T,RHO,XMUE,ZM,AM,X(1))
         D(1,1)=D(1,1)+D11
         D(2,1)=D(2,1)+D21
         D(3,1)=D(3,1)+D31
         D(4,1)=D(4,1)+D41
         D(5,1)=D(5,1)+D51
         D(6,1)=D(6,1)+D61
      ENDIF
      CALL FUNTER (T,RHO,NMODES,AM,ENTRO,U,P,DUT,DUD,DPT,DPD)
      RETURN
      END
C     ***************************************************************
C               *******************************
      SUBROUTINE FTER2(T,RHO,ZM,AM,XMUE,ENTRO,U,P,DUT,DUD,DPT,DPD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D1(3),D2(3)
      COMMON/SCR/D11,D21,D31,D41,D51,D61,D111,D211,D311,D411,D511,D611
      COMMON/COEF/X(3),XMAX(3)
      COMMON/FUN/D(6,3)
      NMODES=3
C     CALCULO DEL MODO LONGITUDINAL
      A=DEXP(-XMAX(1))
      D0=DLOG(1.d0-A)
      D1(1)=A/(1.d0-A)
      D2(1)=-A/(1.d0-A)/(1.d0-A)
      DT=-D1(1)*XMAX(1)/T
      DR=D1(1)*XMAX(1)/2.d0/RHO
      DDT=D2(1)*XMAX(1)*XMAX(1)/T/T-2.d0*DT/T
      DDTD=-D2(1)*XMAX(1)*XMAX(1)/2.d0/T/RHO+DT/RHO/2.d0
      DDD=D2(1)*XMAX(1)*XMAX(1)/4.d0/RHO/RHO-DR/2.d0/RHO
C  CORRECCION SCREENING MODO LONGITUDINAL
C  EL PRIMER TERMINO ES EL LONGITUDINAL CORREGIDO DE SCREENING
C  SI ES NECESARIO. EL SEGUNDO TERMINO ES EL GROUND STATE DE LOS
C  TRES MODOS Y EL TERCERO LA CORRECCION DE SCREENING AL GROUND ST.
      IF (RHO .LE. 1.D+05) THEN
         CALL SCREEN (T,RHO,XMUE,ZM,AM,X(1))
         D(1,1)=D111*D0+0.5d0*1.5344d0*X(1)+D11
         D(2,1)=D211*D0+D111*DT-0.5d0*1.5344d0*X(1)/T+D21
         D(3,1)=D311*D0+D111*DR+0.5d0*1.5344d0*X(1)/2.d0/RHO+D31
         D(4,1)=D411*D0+2.d0*D211*DT+D111*DDT+1.5344d0*X(1)/T/T+D41
         D(5,1)=D511*D0+D211*DR+D311*DT+D111*DDTD
     $        -0.25d0*1.5344d0*X(1)/RHO/T+D51
         D(6,1)=D611*D0+2.d0*D311*DR+D111*DDD
     $        -0.125d0*1.5344d0*X(1)/RHO/RHO+D61
      ELSE
         D(1,1)=D0+0.5d0*1.53444d0*X(1)
         D(2,1)=DT-0.5d0*1.53444d0*X(1)/T
         D(3,1)=DR+0.5d0*1.5344d0*X(1)/2.d0/RHO
         D(4,1)=DDT+1.53444d0*X(1)/T/T
         D(5,1)=DDTD-0.25d0*1.53444d0*X(1)/RHO/T
         D(6,1)=DDD-0.125d0*1.53444d0*X(1)/RHO/RHO
      ENDIF
C  CALCULO MODOS TRANSVERSALES
      DO 1 I=2,NMODES
         A=DEXP(-XMAX(I))
         XM2=XMAX(I)*XMAX(I)
         XM3=XMAX(I)*XM2
c         D(1,I)=1./XM3*(6.4939477d0-A*(XM3+3.*XM2+6.*XMAX(I)+6.)
c     $        -0.5d0*A*A*(XM3+1.5*XM2+1.5*XMAX(I)+0.75)-A**3/3*(XM3+
c     $        XM2+0.66666*XMAX(I)+2./9.))
         twothird=2.d0/3.d0
         D(1,I)=1.d0/XM3*(6.4939477d0-A*(XM3+3.d0*XM2+6.d0*XMAX(I)+6.d0)
     $        -0.5d0*A*A*(XM3+1.5*XM2+1.5d0*XMAX(I)+0.75)-A**3/3.d0*
     $        (XM3+XM2+twothird*XMAX(I)+2.d0/9.d0))
         D1(I)=-3.d0/XMAX(I)*D(1,I)+A/(1.d0-A)
         D2(I)=3.d0/XM2*D(1,I)-3.d0/XMAX(I)*D1(I)-A/(1.d0-A)**2
         D(1,I)=DLOG(1.d0-A)-D(1,I)
         D1(I)=A/(1.d0-A)-D1(I)
         D2(I)=-A/(1.d0-A)/(1.d0-A)-D2(I)
 1    ENDDO
      DO 2 I=2,NMODES
         D(2,I)=-D1(I)*XMAX(I)/T
         D(3,I)=D1(I)*XMAX(I)/RHO/2.d0
         D(4,I)=D2(I)*XMAX(I)*XMAX(I)/T/T-D(2,I)*2.d0/T
         D(5,I)=-D2(I)*XMAX(I)*XMAX(I)/RHO/T/2.d0+D(2,I)/RHO/2.d0
         D(6,I)=D2(I)*XMAX(I)*XMAX(I)/4.d0/RHO/RHO-D(3,I)/2.d0/RHO
 2    ENDDO
      CALL FUNTER (T,RHO,NMODES,AM,ENTRO,U,P,DUT,DUD,DPT,DPD)
      RETURN
      END
C     ****************************************************************
C         CALCULA LAS FUNCIONES TERMODINAMICAS
C     ***************************************************************
      SUBROUTINE FUNTER (T,RHO,NMODES,AM,ENTRO,U,P,DUT,DUD,DPT,DPD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(6)
      COMMON/COEF/X(3),XMAX(3)
      COMMON/FUN/D(6,3)
      C=8.31451D+07*T/AM
      DO 1 I=1,6
         S(I)=0.d0
 1    ENDDO
      DO 3 I=1,6
       DO 2 J=1,NMODES
          S(I)=S(I)+D(I,J)
 2     ENDDO
 3    ENDDO
      F=C*S(1)
      ENTRO=-F/T-C*S(2)
      U=F+T*ENTRO
      P=C*RHO*RHO*S(3)
      DUT=2.d0*U/T-C*T*S(4)
      DUD=-C*T*S(5)
      DPD=2.d0*P/RHO+C*RHO*RHO*S(6)
      DPT=1.d0/T*(P-DUD*RHO*RHO)
      RETURN
      END
C      **************************************************************
C           CALCULA LA CORRECCION DE APANTALLAMIENTO
C      *************************************************************
      SUBROUTINE SCREEN (T,RHO,XMUE,ZM,AM,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ninth
      DIMENSION E(3),C(6),H(6),G(3)
      COMMON/SCR/D11,D21,D31,D41,D51,D61,D111,D211,D311,D411,D511,D611
      DATA E(1),E(2),E(3)/-0.032633d0,1.96414D-03,-1.1159D-04/
      DATA (C(I),I=1,6)/-4.7618D-02,4.7768D-03,-6.5133D-04,
     $  7.7568D-05,-6.6271D-06,2.7384D-07/
c      DATA G(1),G(2),G(3),AK2/-0.02772,1.6687D-03,-9.48D-05,4.1667D-02/
c      H(1)=5.496*ZM*XMUE**0.6666/AM**0.3333/RHO**0.33333
      DATA G(1),G(2),G(3),AK2/-0.02772d0,1.6687D-03,-9.48D-05,
     &     4.1667D-02/

c      H(1)=5.496*ZM*XMUE**0.6666/AM**0.3333/RHO**0.33333
      third=1.d0/3.d0
      twothird=2.d0*third
      ninth=third*third
      H(1)=5.496d0*ZM*XMUE**twothird/AM**third/RHO**third
      DO 1 I=1,5
         H(I+1)=H(I)*H(1)
 1    ENDDO
      SLW=0.d0
      DSLW=0.d0
      DDSLW=0.d0
      SW2=0.d0
      DSW2=0.d0
      DDSW2=0.d0
      IF (X .LE. 4.d0) THEN
         X12=X*X
         DO 2 K=1,3
            SLW=SLW+E(K)*H(K)
            DSLW=DSLW+K*E(K)*H(K)
            DDSLW=DDSLW+K*K*E(K)*H(K)
 2       ENDDO
         DO 3 K=1,6
            SW2=SW2+C(K)*H(K)
            DSW2=DSW2+K*C(K)*H(K)
            DDSW2=DDSW2+K*K*C(K)*H(K)
 3       ENDDO
         DSLW=-third/RHO*DSLW
         DDSLW=-DSLW/RHO+ninth/RHO/RHO*DDSLW
         DSW2=-third/RHO*DSW2
         DDSW2=-DSW2/RHO+ninth/RHO/RHO*DDSW2
         D11=SLW+AK2*SW2*X12
         D21=-AK2*SW2*X12/T*2.d0
         D31=DSLW+AK2*DSW2*X12+AK2*SW2*X12/RHO
         D41=6.d0*AK2*SW2*X12/T/T
         D51=-AK2*DSW2*X12/T*2.d0-AK2*SW2*X12/T/RHO*2.d0
         D61=DDSLW+AK2*DDSW2*X12+2.d0*AK2*DSW2*X12/RHO
      ELSE
         DO 4 K=1,3
            SW2=SW2+G(K)*H(K)
            DSW2=DSW2+K*G(K)*H(K)
            DDSW2=DDSW2+K*K*G(K)*H(K)
 4       ENDDO
         DSW2=-third/RHO*DSW2
         DDSW2=-DSW2/RHO+ninth/RHO/RHO*DDSW2
         D11=0.5d0*SW2*X
         D21=-D11/T
         D31=0.5d0*DSW2*X+D11/2.d0/RHO
         D41=-D21/T+D11/T/T
         D51=-D31/T
         D61=D31/2.d0/RHO-D11/2.d0/RHO/RHO+0.5d0*DDSW2*X+0.5d0*DSW2*X/
     1        2.d0/RHO
         D111=DEXP(-SW2*X)
         D211=D111*SW2*X/T
         D311=-D111*(DSW2+SW2/2.d0/RHO)*X
         D411=D211*SW2*X/T-2.d0*D211/T
         D511=X/T*(D411*SW2+D111*DSW2)+D211/RHO/2.d0
         D611=D311*D311/D111-D111*(DDSW2+DSW2/RHO-SW2/4.d0/RHO/RHO)*X
      ENDIF
      RETURN
      END
