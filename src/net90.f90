! ***********************************************************************
!                       GARCIA-SENZ, et al. (2024)                        
!                       CABEZON el al. (ApJS, 2004)
! ***********************************************************************
!
!     INPUT:
!
!     *** tempin    (D): temperature
!     *** xin       (D): mass fractions
!     *** rho       (D): density (g/cm^3) of the current time-step
!     *** delta     (D): time-step in seconds
!     *** screen    (L): Turns the screening corrections to the rates on/off
!     *** cv        (D): heat capacity.
!     *** val1      (D): (dp/dTemp)*(deltarho/rho**2)
!     *** val2      (D): (dp/dTemp)*(deltarho/rho**2)*Temp
!     *** dUedYein  (D): dUe/dYe (from EOS)
!     *** theta     (D): =1 for fully implicit method.
!
!     OUTPUT:
!     *** tempout   (D): temperature
!     *** xout      (D): mass fractions
!     *** sumtot    (D): total mass fraction (=1 if conserved)
!     *** nucenergy (D): total nuclear energy generation rate (erg/s)
!     *** k_iter    (I): Total number of Newton-Raphson iterations
!
! ***************
! ***********************************************************************

  SUBROUTINE net90(tempin,xin,rho,delta,screen,ecapture,tabulated,cv,val1,val2,dUedYein,theta,&
                 & tempout,xout,sumtot,nucenergy,k_iter)

      USE nuclear90_module

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(in)  :: tempin,delta,rho,cv,val1,val2,theta,dUedYein
      DOUBLE PRECISION,INTENT(out) :: tempout,sumtot,nucenergy
      DOUBLE PRECISION,INTENT(in),DIMENSION(niso)  :: xin
      DOUBLE PRECISION,INTENT(out),DIMENSION(niso) :: xout
      
      LOGICAL,INTENT(in)   :: screen,ecapture,tabulated
      INTEGER, INTENT(out) :: k_iter
      
      DOUBLE PRECISION EF,maxrow,temp,correction
      DOUBLE PRECISION,DIMENSION(iso) :: vindep,norm,vindepold,vindepold2,y_in,y_out

! **********************
! *** Initialization ***
! **********************
      dUedYe=dUedYein

      if(theta.eq.0.d0) then
         print *,'Wrong theta value!'
         print *,'0 < Theta <= 1'
         stop
      endif

      do i=1,7
        q(i)=0.d0
        do j=1,8
           fit(i,j)=0.d0
        enddo
      enddo

      do i=1,6
         do j=1,24
            choose(i,j)=0.d0
         enddo
      enddo

      y_in(iso)=ye
      y_in(1)=0.d0
      do i=2,iso-1
        y_in(i)=xin(i-1)/a(i)
        vindepold(i)=1.d50
        vindep(i)=1.d50
      enddo
      y_out(:)=y_in(:)
      tempout=tempin

      k_iter=0

! ***********************************************************************
! ***********************************************************************

! OOOOOOOOOOOOOOOOOOO START RUN OOOOOOOOOOOOOOOOOOO

      ! Newton-Raphson iterations
      ! Minimum 2 NR iterations as mass conservation is always slightly better
      do while((k_iter.lt.2.or.(correction.gt.NR_limit)) .and. (k_iter.lt.NR_iter))
         k_iter=k_iter+1
         vindepold2(:)=vindepold(:)
         vindepold(:)=vindep(:)

         temp=tempin*(1.d0-theta)+tempout*theta   !Temperature at which we calculate
         y(:)=y_in(:)*(1.d0-theta)+y_out(:)*theta !Molar fractions at which we calculate
         ye=y(iso)
         
         if(screen)call chempot(temp,rho)         !Chemical potentials

         ! Calculation of the rates.
         ! Until 16O(a,g)20Ne, including heavy ion reactions: calculated with WWFHZ78
         !    using subroutines 'effect' and 'deffcalc'.
         ! From 20Ne(a,g)24Mg until the end of the network: calculated with RT00
         !    using subroutines 'direct' and 'invers'.
         ! Subroutine direct also calculates nu and nubar energy losses
         call effect(temp)
         call deffcalc(temp)
         call direct(temp,rho,ecapture,tabulated) !We need rho here for electronic rates, but can be removed when using table.
         call inverse(temp)
         do i=1,rates
            eff(i)=eff(i)*rho
            deff(i)=deff(i)*rho
         enddo
         eff(5)=eff(5)*rho     !triple alpha
         deff(5)=deff(5)*rho
         do i=138,rates        
            if(i.eq. (rates-2) .or. i.eq.(rates-1)) cycle !rates-2, rates-1 are photodisintegrations
            l(i)=l(i)*rho      !x dens because these are not photodisintegrations 
            dl(i)=dl(i)*rho
         enddo

         ! Rates corrections due to the chemical potentials.
         ! It only affects the direct reaction rates and their derivatives.
         if(screen)then
            do i=1,rates
               EF=exp(deltamukbt(i))
               eff(i)=eff(i)*EF
               deff(i)=deff(i)*EF-2.d0*eff(i)*deltamukbt(i)/temp
            enddo
            do i=138,rates
               if(i.eq. (rates-2) .or. i.eq.(rates-1)) cycle !rates-2, rates-1 are photodisintegrations
               EF=exp(deltamukbt(i))
               l(i)=l(i)*EF
               dl(i)=dl(i)*EF-2.d0*l(i)*deltamukbt(i)/temp
            enddo
         endif
         
         ! Function and matrix calculation
         ! myfunc: Function evaluation that will become the indepent term.
         ! myfuncprima: Same as myfunc but evaluated with the rates derivatives instead
         !              of the rates themselves.
         ! phicalc: matrix of coefficients
         call myfunc()

         !Discretizacion correction
         do i=2,iso
            f(i)=f(i)-(y_out(i)-y_in(i))/delta
         enddo
         call myfuncprima()

         call phicalc(delta,temp,rho,theta)

         ! Calculates the temperature term of the energy equation.
         ! Includes neutrino energy rates for emmited enu and enubar.
         phi(1,1)=cv-theta*(val1-delta*(dEneutr*y(2)+dEaneutr*y(3)))

         do j=1,iso
            do i=2,iso
               phi(i,j)=phi(i,j)*theta
            enddo
         enddo

         ! Diagonal corrections
         do i=2,iso
            phi(i,i)=phi(i,i)+1.d0-theta
         enddo 
         vindep(1)=val2-cv*(tempout-tempin)-delta*(Eneutr*y(2)+Eaneutr*y(3))!+dif*deltadif
         do i=3,iso-1
            vindep(1)=vindep(1)-phi(1,i)*(y_out(i)-y_in(i))
         enddo
         vindep(1)=vindep(1)-(phi(1,iso)-theta*delta*(dEneutrdYe*y(2)+dEaneutrdYe*y(3)))*(y_out(iso)-y_in(iso))
         vindep(1)=vindep(1)-(phi(1,2)-theta*delta*(Eneutr+Eaneutr))*(y_out(2)-y_in(2))
         f(1)=vindep(1)

! OOOOOOOOOOOOOOOOOOO SYSTEM RESOLUTION OOOOOOOOOOOOOOOOOOO

         ! Matrix normalization.
         ! norm   = array with the normalization values for every equation of the system.
         ! vindep = array containing the independent terms of the equations.
         do i=1,iso
            maxrow=abs(phi(i,1))
            do j=1,iso
               if(abs(phi(i,j)).gt.maxrow) maxrow=abs(phi(i,j))
            enddo
            norm(i)=maxrow
            if(maxrow.eq.0.d0)then
               write(*,*) 'Zero row in the matrix...!!'
               stop
            endif
            do j=1,iso
               phi(i,j)=phi(i,j)/maxrow
            enddo
         enddo
         vindep(1)=vindep(1)/norm(1)

         do i=2,iso
            vindep(i)=delta*f(i)/norm(i)
         enddo
         ! An additional column in the matrix is needed for subroutine 'eigen' to run correctly.
         ! In this column is stored the array of the independent terms of the system.
         ! (PRANTZOS ET AL. returns result in vindep)
         do i=1,iso
            phi(i,iso+1)=vindep(i)
         enddo
         call eigen(iso,vindep,vindep) !Solves system of equations

         correction=abs(vindep(1)/temp)

         !If vindep is increasing, NR is not converging. We reduce DeltaY (i.e. vindep)
         !and force to pass the limit in these cases.
         !Equivalent to reducing the timestep to allow convergence
         if(k_iter.gt.2.and.k_iter.lt.9.and.(abs(vindep(1)).gt.abs(vindepold(1))))then
            vindep(:)=vindep(:)*0.21d0
            correction=NR_limit/1.d6 !1.d-16
         endif

         tempout=tempout+vindep(1)          !Update Tempout

         do i=2,iso                         !Update Y
            y_out(i)=y_out(i)+vindep(i)
            if(y_out(i).le.ymin)y_out(i)=0.d0
          enddo

         if(isnan(vindep(1)))print *,val2,norm(1)
      enddo

! OOOOOOOOOOOOOOOOOOO SYSTEM RESOLUTION END OOOOOOOOOOOOOOOOOOO
      sumtot=0.d0

      !Only until iso-1 because electrons have X=0 due to A=0
      do i=2,iso-1
         xout(i-1)=y_out(i)*a(i)
         sumtot=sumtot+xout(i-1)
      enddo

      ye=y_out(iso)   !Update ye

      nucenergy=0.d0  !Nuclear and neutrino energy calculation
      do i=2,iso
         nucenergy=nucenergy+(y_out(i)-y_in(i))*BE(i)*MeVperparticle2erg
      enddo
      Eneutr=Eneutr*delta*y(2)
      Eaneutr=Eaneutr*delta*y(3)
      nucenergy=nucenergy!+Eneutr
      return
    end subroutine net90

!*******************************************************************
! *** SUBROUTINE direct.
! *** Calculates the direct nuclear reaction rates using RT00
! *** From 20Ne(a,g)24Mg to 56Ni(a,g)60Zn
      SUBROUTINE direct(temp,rho,ecapture,tabulated)

      USE nuclear90_module

      IMPLICIT NONE

      DOUBLE PRECISION,INTENT(in)::temp,rho
      LOGICAL,INTENT(in)::ecapture,tabulated

      DOUBLE PRECISION t9,t913,t923,t953,t9i,t9i2,t9i13,t9i23,t9i43,lt9
      DOUBLE PRECISION ue,ue2,ue3,ue4,ue5,ue6,due,ueaux,t92,t93,t94,t95,t96,t9i3
      DOUBLE PRECISION dUdYe
      !DOUBLE PRECISION aue,b,bx,x,x2,aux,invaux,dxdYe,dUedx
      DOUBLE PRECISION,DIMENSION(nc):: interpolated

      t9=temp/1.d9
      t913=t9**(third)
      t923=t913*t913
      t953=t9**(5.d0*third)
      t9i=1.d0/t9
      t9i2=t9i*t9i
      t9i13=1.d0/t913
      t9i23=t9i13*t9i13
      t9i43=t9i23*t9i23
      lt9=log(t9)

      do i=8,rates
        coef(i)=fit(i,2)*t9i+fit(i,3)*t9i13+fit(i,4)*t913+fit(i,5)*t9&
     &         +fit(i,6)*t953+fit(i,7)*lt9
        dcoef(i)=(-fit(i,2)*t9i2+(-fit(i,3)*t9i43+fit(i,4)*t9i23+5.d0&
     &          *fit(i,6)*t923)*(third)+fit(i,5)+fit(i,7)*t9i)*1.d-9
        eff(i)=exp(fit(i,1)+coef(i))
        deff(i)=eff(i)*dcoef(i)
      enddo
          
      ! From most of the (a,p) and (a,n) we only have the inverse reaction.
      ! For that reason we use the inverse of the inverse as direct reaction.
      ! So the direct reactions from the tables are the inverse of our network.
      do i=138,154
        l(i)=eff(i)
        dl(i)=deff(i)
      enddo
      l(157)=eff(157)
      dl(157)=deff(157)
      l(158)=eff(158)
      dl(158)=deff(158)

      ! Electron capture Rate (effe) and Derivatives (deffe=deffe/dT and deffedYe=deffe/dYe)
      ! Energy of neutrinos and dUe/dYe:
      if(ecapture) then
        if(tabulated) then
          call interp(temp,rho*y(iso),interpolated)
          effe=interpolated(1)
          deffe=interpolated(2)
          deffedYe=interpolated(3)*rho
          Eneutr=interpolated(4)

          dEneutr=interpolated(5)
          dEneutrdYe=interpolated(6)*rho

          effp=interpolated(7)
          deffp=interpolated(8)
          deffpdYe=interpolated(9)*rho
          Eaneutr=interpolated(10)
          dEaneutr=interpolated(11)
          dEaneutrdYe=interpolated(12)*rho

        else

          ueaux=1.d-2*(rho*y(iso))**third
          ue=min(ueaux,6.d-6*rho*y(iso)*t9i2) !Chemical Potential
          t92=t9*t9
          t93=t92*t9
          t94=t93*t9
          t95=t94*t9
          t96=t95*t9
          t9i3=t9i2*t9i
          ue2=ue*ue
          ue3=ue2*ue
          ue4=ue3*ue
          ue5=ue4*ue
          ue6=ue5*ue

          effe=9.4d-5*(ue5+4.d-2*(ue3*t92+4.d-1*t95)) !1/s
          Eneutr=3.8d-4*(ue6*sixth+6.d-2*ue4*t92+3.2d-3*t96) !erg/s
          if(ue.eq.ueaux)then
            deffe=7.52d-6*(ue3*t9+t94)
            dEneutr=4.56d-7*ue4*t9+7.296d-8*t95
            dUdYe=third*1.d-2*rho**third*y(iso)**(-2*third)
            deffedYe=4.7d-4*dUdYe*(ue4+2.4d-2*ue2*t92)
            dEneutrdYe=3.8d-4*dUdYe*(ue5+2.4d-3*ue3*t92)
          else
            dUe=-1.2d-5*rho*y(iso)*t9i3
            dUdYe=6.d-6*rho*t9i2
            deffe=4.7d-4*(dUe*ue4+1.6d-2*(1.5d0*dUe*ue2*t92+ue3*t9+t94))
            deffedYe=4.7d-4*dUdYe*(ue4+2.4d-2*ue2*t92)
            dEneutr=3.8d-4*(dUe*(ue5+2.4d-3*ue3*t92)+19.2d-3*t95)
            dEneutrdYe=3.8d-4*dUdYe*(ue5+2.4d-3*ue3*t92)
          endif
          deffe=deffe*1.d-9
          Eneutr=Eneutr*4.93d17
          dEneutr=dEneutr*4.93d17*1.d-9
          dEneutrdYe=dEneutrdYe*4.93d17

          
          ! Calculation of dUe/dYe from Chandrasekhar expression
          ! Not necessay if dUedYe comes from EOS
          
          !aue=6.01d22
          !b=1.1218d-18
          !bx=1.006d-02
          !x=bx*(rho*ye)**third
          !x2=x*x
          !aux=sqrt(x2+1.d0)
          !invaux=1.d0/aux
          !dUedx=aue*(24.d0*x2*x2-24.d0*(aux-1.d0)*x2+6.d0)*invaux
          !dUedx=dUedx+aue*b*temp*temp*(-2.d0+(1.d0+3.d0*x)*invaux-(aux-1.d0)/x2)
          !dxdYe=bx*third*rho**third*ye**(-2.d0*third)
          !dUedYe=dUedx*dxdYe/rho
        endif
      else !Possibility to deactivate electron captures if ecapture=.false.
        effe=0.d0
        deffe=0.d0
        deffedYe=0.d0
        Eneutr=0.d0
        dUedYe=0.d0
        dEneutr=0.d0
        dEneutrdYe=0.d0
        effp=0.d0
        deffp=0.d0
        deffpdYe=0.d0
        Eaneutr=0.d0
        dEaneutr=0.d0
        dEaneutrdYe=0.d0
      endif
     return

    end subroutine direct


!*******************************************************************
! *** SUBROUTINE inverse.
! *** Calculates the inverse nuclear reaction rates using RT00
      SUBROUTINE inverse(temp)

      USE nuclear90_module

      IMPLICIT NONE
      INTEGER fin,ini
      DOUBLE PRECISION,INTENT(in)::temp
      DOUBLE PRECISION t9,t9i,t932,val1,val2,val3,val4,part
      DOUBLE PRECISION,DIMENSION(19)::aux1,aux2

      t9=temp/1.d9
      t9i=1.d0/t9
      t932=t9**1.5d0

      if(t9.ge.0.01d0.and.t9.lt.0.15d0) then
        k=1
      else if(t9.ge.0.15d0.and.t9.lt.0.2d0) then
        k=2
      else if(t9.ge.0.2d0.and.t9.lt.0.3d0) then
        k=3
      else if(t9.ge.0.3d0.and.t9.lt.0.4d0) then
        k=4
      else if(t9.ge.0.4d0.and.t9.lt.0.5d0) then
        k=5
      else if(t9.ge.0.5d0.and.t9.lt.0.6d0) then
        k=6
      else if(t9.ge.0.6d0.and.t9.lt.0.7d0) then
        k=7
      else if(t9.ge.0.7d0.and.t9.lt.0.8d0) then
        k=8
      else if(t9.ge.0.8d0.and.t9.lt.0.9d0) then
        k=9
      else if(t9.ge..9d0.and.t9.lt.1.d0) then
        k=10
      else if(t9.ge.1.d0.and.t9.lt.1.5d0) then
        k=11
      else if(t9.ge.1.5d0.and.t9.lt.2.d0) then
        k=12
      else if(t9.ge.2.d0.and.t9.lt.2.5d0) then
        k=13
      else if(t9.ge.2.5d0.and.t9.lt.3.d0) then
        k=14
      else if(t9.ge.3.d0.and.t9.lt.3.5d0) then
        k=15
      else if(t9.ge.3.5d0.and.t9.lt.4.d0) then
        k=16
      else if(t9.ge.4.d0.and.t9.lt.4.5d0) then
        k=17
      else if(t9.ge.4.5d0.and.t9.lt.5.d0) then
        k=18
      else if(t9.ge.5.d0.and.t9.lt.6.d0) then
        k=19
      else if(t9.ge.6.d0.and.t9.lt.7.d0) then
        k=20
      else if(t9.ge.7.d0.and.t9.lt.8.d0) then
        k=21
      else if(t9.ge.8.d0.and.t9.lt.9.d0) then
        k=22
      else if(t9.ge.9.d0.and.t9.lt.10.d0) then
        k=23
      else if(t9.ge.10.d0) then
        k=24
      endif

      do i=1,17
        aux1(i)=l(i+137)
        aux2(i)=dl(i+137)
      enddo
      aux1(18)=l(157)
      aux2(18)=dl(157)
      aux1(19)=l(158)
      aux2(19)=dl(158)

      val1=11.6045d0*t9i
      val2=1.5d0*log(t9)
      val3=val1*t9i*1.d-9
      val4=1.5d-9*t9i

!  *** target=target of the direct reaction not of the inverse ***
      do i=8,rates  
        if(i.le.137 .or.i.eq.(rates-2) .or. i.eq.rates-1) then
           fin=final(i)
           ini=target(i)
           if(i.eq.rates-1) then 
              part=1.d0
           else
              part=choose(ini,k)/choose(fin,k)
           endif
           l(i)=part*exp(fit(i,8)+coef(i)-val1*q(i)+val2)
           dl(i)=l(i)*(dcoef(i)+val3*q(i)+val4)
        endif
      enddo
      
      do i=138,rates   !These are not photodesintegrations so they don't have val2. Except rates-2 and rates-1
        if(i.eq. (rates-2) .or. i.eq.(rates-1))cycle
        fin=final(i)
        ini=target(i)
        if(i.eq.rates) then 
           part=1.d0
        else
           part=choose(ini,k)/choose(fin,k)
        endif
        l(i)=part*exp(fit(i,8)+coef(i)-val1*q(i))
        dl(i)=l(i)*(dcoef(i)+val3*q(i))
      enddo 
! --------------------------------------------------------------

      ! For most of the (a,p) and (a,n) we only have the inverse reaction.
      ! For that reason we use the inverse of the inverse as direct reaction.
      ! So the direct reactions from the tables are the inverse of our network.
      do i=138,154
        eff(i)=l(i)
        deff(i)=dl(i)
      enddo
      eff(157)=l(157)
      deff(157)=dl(157)
      eff(158)=l(158)
      deff(158)=dl(158)

      do i=1,17
        l(i+137)=aux1(i)
        dl(i+137)=aux2(i)
      enddo
      l(157)=aux1(18)
      dl(157)=aux2(18)
      l(158)=aux1(19)
      dl(158)=aux2(19)
      return
    end subroutine inverse


! *******************************************************************

      SUBROUTINE myfunc()

      USE nuclear90_module

      IMPLICIT NONE
      INTEGER a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a88,&
            & a89,a90,a91,a128,a129,a138,a139,a148,a149,ant

      do i=1,iso
        f(i)=0.d0
      enddo

      a8=8
      a9=9
      a10=10
      a11=11
      a12=12
      a13=13
      a14=14
      a15=15
      a16=16
      a17=17
      a18=18
      a19=19
      a20=20
      a88=88
      a89=89
      a90=90
      a91=91
      a128=128
      a129=129
      a138=138
      a139=139
      a148=148
      a149=149
      ant=7

      do while(ant.le.81)
        f(a8)=l(a13)*y(a13)+l(a14)*y(a14)+y(ant)*y(2)*eff(a8)&
          &   -y(a8)*(l(a8)+eff(a13)*y(2)+eff(a14)*y(3))
        f(a9)=y(ant)*y(4)*eff(a128)&
          &   +y(a10)*y(3)*eff(a9)+y(a12)*y(2)*eff(a91)&
          &   -y(a9)*(l(a9)+l(a91)+l(a128))
        f(a10)=l(a9)*y(a9)+eff(a10)*y(a13)*y(3)+eff(a89)*y(a14)*y(2)&
          &   +eff(a138)*y(ant)*y(4)&
          &   -y(a10)*(l(a10)+l(a89)+y(3)*(l(a138)+eff(a9)))
        f(a11)=l(a88)*y(a14)+l(a15)*y(a15)+eff(a11)*y(ant)*y(3)&
          &   -y(a11)*(l(a11)+eff(a88)*y(2)+eff(a15)*y(3))
        f(a12)=l(a91)*y(a9)+eff(a12)*y(a14)*y(3)+eff(a90)*y(a15)*y(2)&
          &   +eff(a148)*y(4)*y(ant)&
          &   -y(a12)*(l(a12)+l(a90)+y(2)*(l(a148)+eff(a91)))
        f(a13)=l(a10)*y(a10)+eff(a13)*y(a8)*y(2)&
          &   -y(a13)*(l(a13)+eff(a10)*y(3))
        f(a14)=l(a89)*y(a10)+l(a12)*y(a12)&
          &   +eff(a14)*y(a8)*y(3)+eff(a88)*y(a11)*y(2)&
          &   -y(a14)*(l(a88)+l(a14)+eff(a89)*y(2)+eff(a12)*y(3))
        f(a15)=l(a90)*y(a12)+eff(a15)*y(a11)*y(3)&
          &   -y(a15)*(l(a15)+eff(a90)*y(2))
        f(4)=f(4)+l(a128)*y(a9)+l(a138)*y(a10)*y(3)+l(a148)*y(a12)*y(2)&
          &   -y(4)*y(ant)*(eff(a128)+eff(a138)+eff(a148))
        f(2)=f(2)+l(a8)*y(a8)+l(a13)*y(a13)+l(a88)*y(a14)+l(a89)*y(a10)&
          &   +l(a90)*y(a12)+l(a91)*y(a9)+eff(a148)*y(ant)*y(4)&
          &   -y(2)*(eff(a8)*y(ant)+eff(a13)*y(a8)+eff(a88)*y(a11)&
          &   +eff(a89)*y(a14)+eff(a90)*y(a15)+y(a12)*(eff(a91)+l(a148)))
        f(3)=f(3)+l(a11)*y(a11)+l(a14)*y(a14)+l(a15)*y(a15)+l(a10)*y(a10)&
          &   +l(a12)*y(a12)+l(a9)*y(a9)+eff(a138)*y(ant)*y(4)&
          &   -y(3)*(eff(a11)*y(ant)+eff(a14)*y(a8)+eff(a15)*y(a11)&
          &   +eff(a12)*y(a14)+eff(a10)*y(a13)+y(a10)*(eff(a9)+l(a138)))
        ant=a9
        if(ant.eq.81)exit
        f(a9)=f(a9)+l(a129)*y(a17)+l(a139)*y(a18)*y(3)&
          &   +l(a149)*y(a20)*y(2)+l(a16)*y(a16)+l(a19)*y(a19)&
          &   -y(a9)*(eff(a16)*y(2)+eff(a19)*y(3)&
          &   +y(4)*(eff(a129)+eff(a139)+eff(a149)))
        a8=a16
        a9=a17
        a10=a18
        a11=a19
        a12=a20
        a13=a13+8
        a14=a14+8
        a15=a15+8
        a16=a16+8
        a17=a17+8
        a18=a18+8
        a19=a19+8
        a20=a20+8
        a88=a88+4
        a89=a89+4
        a90=a90+4
        a91=a91+4
        a128=a129
        a129=a129+1
        a138=a139
        a139=a139+1
        a148=a149
        a149=a149+1
      end do
      f(2)=f(2)+cp*.5d0*eff(1)*y(5)**2.d0+op*.5d0*eff(3)*y(6)**2.d0-effe*y(2)+effp*y(3)
      f(2)=f(2)+y(79)*y(4)*eff(158)-y(2)*y(88)*l(158)-&
         & y(5)*y(2)*eff(160)+y(90)*l(160)+y(4)*y(90)*eff(161)-y(2)*y(6)*l(161)
      f(3)=f(3)+effe*y(2)-effp*y(3)-y(3)*y(88)*eff(159)+y(89)*l(159)
      f(4)=f(4)+3.d0*l(5)*y(5)+l(6)*y(6)+l(7)*y(7)&
        &   -y(4)*(.5d0*eff(5)*y(4)**2.d0+eff(6)*y(5)+eff(7)*y(6))&
        &   +.5d0*(ca*eff(1)*y(5)**2.d0+oa*eff(3)*y(6)**2.d0)&
        &   +eff(2)*y(5)*y(6)-y(4)*y(79)*eff(158)+y(2)*y(88)*l(158)-&
        &   y(4)*y(90)*eff(161)+y(2)*y(6)*l(161)
      f(5)=f(5)+l(6)*y(6)+sixth*eff(5)*y(4)**3.d0-y(5)*(l(5)+eff(6)*y(4)&
        &   +eff(1)*y(5)+eff(2)*y(6))-y(2)*y(5)*eff(160)+y(90)*l(160)
      f(6)=f(6)+l(7)*y(7)+eff(6)*y(5)*y(4)-y(6)*(l(6)+eff(7)*y(4)&
        &   +eff(2)*y(5)+eff(3)*y(6))
      f(7)=f(7)+l(128)*y(9)+l(138)*y(10)*y(3)+l(148)*y(12)*y(2)&
        &   +l(8)*y(8)+l(11)*y(11)+y(6)*y(4)*eff(7)&
        &   -y(7)*(l(7)+eff(8)*y(2)+eff(11)*y(3)&
        &   +y(4)*(eff(128)+eff(138)+eff(148)))&
        &   +ca*.5d0*eff(1)*y(5)**2.d0
      f(9)=f(9)+eff(2)*y(5)*y(6)
      f(12)=f(12)+cp*.5d0*eff(1)*y(5)**2.d0
      f(6)=f(6)+y(4)*y(90)*eff(161)-y(2)*y(6)*l(161)
      f(17)=f(17)+oa*.5d0*eff(3)*y(6)**2.d0
      f(28)=f(28)+op*.5d0*eff(3)*y(6)**2.d0
      f(79)=f(79)-y(79)*y(4)*eff(158)+y(2)*y(88)*l(158)
      f(88)=y(4)*y(79)*eff(158)-y(2)*y(88)*l(158)-y(3)*y(88)*eff(159)+&
              & y(89)*l(159)
      f(89)=y(3)*y(88)*eff(159)-y(89)*l(159)
      f(90)=y(2)*y(5)*eff(160)-y(90)*l(160)-y(4)*y(90)*eff(161)+y(2)*y(6)*l(161)
      f(iso)=-effe*y(2)
      RETURN
      END

! *******************************************************************

      SUBROUTINE myfuncprima()

      USE nuclear90_module

      IMPLICIT NONE
      INTEGER a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a88,&
          &   a89,a90,a91,a128,a129,a138,a139,a148,a149,ant

      do i=1,iso
        fp(i)=0.d0
      enddo

      a8=8
      a9=9
      a10=10
      a11=11
      a12=12
      a13=13
      a14=14
      a15=15
      a16=16
      a17=17
      a18=18
      a19=19
      a20=20
      a88=88
      a89=89
      a90=90
      a91=91
      a128=128
      a129=129
      a138=138
      a139=139
      a148=148
      a149=149
      ant=7

      do while(ant.le.81)
          fp(a8)=dl(a13)*y(a13)+dl(a14)*y(a14)+y(ant)*y(2)*deff(a8)&
          &     -y(a8)*(dl(a8)+deff(a13)*y(2)+deff(a14)*y(3))
          fp(a9)=y(ant)*y(4)*deff(a128)&
          &     +y(a10)*y(3)*deff(a9)+y(a12)*y(2)*deff(a91)&
          &     -y(a9)*(dl(a9)+dl(a91)+dl(a128))
          fp(a10)=dl(a9)*y(a9)+deff(a10)*y(a13)*y(3)+deff(a89)*y(a14)*y(2)&
          &      +deff(a138)*y(ant)*y(4)&
          &      -y(a10)*(dl(a10)+dl(a89)+y(3)*(dl(a138)+deff(a9)))
          fp(a11)=dl(a88)*y(a14)+dl(a15)*y(a15)+deff(a11)*y(ant)*y(3)&
          &      -y(a11)*(dl(a11)+deff(a88)*y(2)+deff(a15)*y(3))
          fp(a12)=dl(a91)*y(a9)+deff(a12)*y(a14)*y(3)+deff(a90)*y(a15)*y(2)&
          &      +deff(a148)*y(4)*y(ant)&
          &      -y(a12)*(dl(a12)+dl(a90)+y(2)*(dl(a148)+deff(a91)))
          fp(a13)=dl(a10)*y(a10)+deff(a13)*y(a8)*y(2)&
          &      -y(a13)*(dl(a13)+deff(a10)*y(3))
          fp(a14)=dl(a89)*y(a10)+dl(a12)*y(a12)&
          &      +deff(a14)*y(a8)*y(3)+deff(a88)*y(a11)*y(2)&
          &      -y(a14)*(dl(a88)+dl(a14)+deff(a89)*y(2)+deff(a12)*y(3))
          fp(a15)=dl(a90)*y(a12)+deff(a15)*y(a11)*y(3)&
          &      -y(a15)*(dl(a15)+deff(a90)*y(2))
          fp(4)=fp(4)+dl(a128)*y(a9)+dl(a138)*y(a10)*y(3)&
          &    +dl(a148)*y(a12)*y(2)&
          &    -y(4)*y(ant)*(deff(a128)+deff(a138)+deff(a148))
          fp(2)=fp(2)+dl(a8)*y(a8)+dl(a13)*y(a13)+dl(a88)*y(a14)&
          &    +dl(a89)*y(a10)&
          &    +dl(a90)*y(a12)+dl(a91)*y(a9)+deff(a148)*y(ant)*y(4)&
          &    -y(2)*(deff(a8)*y(ant)+deff(a13)*y(a8)+deff(a88)*y(a11)&
          &    +deff(a89)*y(a14)+deff(a90)*y(a15)+y(a12)*(deff(a91)&
          &    +dl(a148)))
          fp(3)=fp(3)+dl(a11)*y(a11)+dl(a14)*y(a14)+dl(a15)*y(a15)&
          &    +dl(a10)*y(a10)&
          &    +dl(a12)*y(a12)+dl(a9)*y(a9)+deff(a138)*y(ant)*y(4)&
          &    -y(3)*(deff(a11)*y(ant)+deff(a14)*y(a8)+deff(a15)*y(a11)&
          &    +deff(a12)*y(a14)+deff(a10)*y(a13)+y(a10)*(deff(a9)&
          &    +dl(a138)))

      ant=a9
      if(ant.eq.81)exit
      fp(a9)=fp(a9)+dl(a129)*y(a17)+dl(a139)*y(a18)*y(3)&
      &     +dl(a149)*y(a20)*y(2)+dl(a16)*y(a16)+dl(a19)*y(a19)&
      &     -y(a9)*(deff(a16)*y(2)+deff(a19)*y(3)&
      &     +y(4)*(deff(a129)+deff(a139)+deff(a149)))
      a8=a16
      a9=a17
      a10=a18
      a11=a19
      a12=a20
      a13=a13+8
      a14=a14+8
      a15=a15+8
      a16=a16+8
      a17=a17+8
      a18=a18+8
      a19=a19+8
      a20=a20+8
      a88=a88+4
      a89=a89+4
      a90=a90+4
      a91=a91+4
      a128=a129
      a129=a129+1
      a138=a139
      a139=a139+1
      a148=a149
      a149=a149+1
    end do

    fp(2)=fp(2)+cp*.5d0*deff(1)*y(5)**2.d0+op*.5d0*deff(3)*y(6)**2.d0-deffe*y(2)+deffp*y(3)
    fp(2)=fp(2)+y(79)*y(4)*deff(158)-y(2)*y(88)*dl(158)-y(5)*y(2)*deff(160)+&
        &   y(90)*dl(160)+y(4)*y(90)*deff(161)-y(2)*y(6)*dl(161)
    fp(3)=fp(3)+deffe*y(2)-deffp*y(3)-y(3)*y(88)*deff(159)+y(89)*dl(159)
    fp(4)=fp(4)+3.d0*dl(5)*y(5)+dl(6)*y(6)+dl(7)*y(7)&
    &    -y(4)*(.5d0*deff(5)*y(4)**2.d0+deff(6)*y(5)+deff(7)*y(6))&
    &    +.5d0*(ca*deff(1)*y(5)**2.d0+oa*deff(3)*y(6)**2.d0)&
    &    +deff(2)*y(5)*y(6)-y(4)*y(79)*deff(158)+y(2)*y(88)*dl(158)-&
    &    y(4)*y(90)*deff(161)+y(2)*y(6)*dl(161)
    fp(5)=dl(6)*y(6)+sixth*deff(5)*y(4)**3.d0-y(5)*(dl(5)+deff(6)*y(4)&
    &    +deff(1)*y(5)+deff(2)*y(6))-y(2)*y(5)*deff(160)+y(90)*dl(160)
    fp(6)=dl(7)*y(7)+deff(6)*y(5)*y(4)-y(6)*(dl(6)+deff(7)*y(4)&
    &    +deff(2)*y(5)+deff(3)*y(6))
    fp(6)=fp(6)+y(4)*y(90)*deff(161)-y(2)*y(6)*dl(161)
    fp(7)=dl(128)*y(9)+dl(138)*y(10)*y(3)+dl(148)*y(12)*y(2)&
    &    +dl(8)*y(8)+dl(11)*y(11)+y(6)*y(4)*deff(7)&
    &    -y(7)*(dl(7)+deff(8)*y(2)+deff(11)*y(3)&
    &    +y(4)*(deff(128)+deff(138)+deff(148)))&
    &    +ca*.5d0*deff(1)*y(5)**2.d0
    fp(9)=fp(9)+deff(2)*y(5)*y(6)
    fp(12)=fp(12)+cp*.5d0*deff(1)*y(5)**2.d0
    fp(17)=fp(17)+oa*.5d0*deff(3)*y(6)**2.d0
    fp(28)=fp(28)+op*.5d0*deff(3)*y(6)**2.d0
    fp(79)=fp(79)-y(79)*y(4)*deff(158)+y(2)*y(88)*dl(158)
    fp(88)=y(4)*y(79)*deff(158)-y(2)*y(88)*dl(158)-y(3)*y(88)*deff(159)+&
            & y(89)*dl(159)
    fp(89)=y(3)*y(88)*deff(159)-y(89)*dl(159)
    fp(90)=y(2)*y(5)*deff(160)-y(90)*dl(160)-y(4)*y(90)*deff(161)+y(2)*y(6)*dl(161)

    fp(iso)=-deffe*y(2)
    RETURN
    END

! *******************************************************************

      SUBROUTINE phicalc(delta,temp,rho,theta)

      USE nuclear90_module

      IMPLICIT NONE
      INTEGER a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a88,&
            & a89,a90,a91,a128,a129,a138,a139,a148,a149,ant

      DOUBLE PRECISION,INTENT(in)::delta,temp,rho,theta
      DOUBLE PRECISION kbt,nakbt,ne,ae,gam,gamma
      DOUBLE PRECISION funcion,ggt1,glt1

      phi(:,:)=0.d0

      a8=8
      a9=9
      a10=10
      a11=11
      a12=12
      a13=13
      a14=14
      a15=15
      a16=16
      a17=17
      a18=18
      a19=19
      a20=20
      a88=88
      a89=89
      a90=90
      a91=91
      a128=128
      a129=129
      a138=138
      a139=139
      a148=148
      a149=149
      ant=7

      do while(ant.le.81)
        phi(a8,2)=-delta*(eff(a8)*y(ant)-eff(a13)*y(a8))
        phi(a8,3)=-delta*(-eff(a14)*y(a8))
        phi(a8,ant)=-delta*eff(a8)*y(2)
        phi(a8,a8)=1.d0-delta*(-(l(a8)+eff(a13)*y(2)+eff(a14)*y(3)))
        phi(a8,a13)=-delta*l(a13)
        phi(a8,a14)=-delta*l(a14)

        phi(a9,2)=-delta*eff(a91)*y(a12)
        phi(a9,3)=-delta*eff(a9)*y(a10)
        phi(a9,4)=-delta*eff(a128)*y(ant)
        phi(a9,ant)=-delta*eff(a128)*y(4)
        phi(a9,a9)=1.d0-delta*(-(l(a128)+l(a9)+l(a91)))
        phi(a9,a10)=-delta*eff(a9)*y(3)
        phi(a9,a12)=-delta*eff(a91)*y(2)

        phi(a10,2)=-delta*eff(a89)*y(a14)
        phi(a10,3)=-delta*(-y(a10)*(eff(a9)+l(a138))+eff(a10)*y(a13))
        phi(a10,4)=-delta*eff(a138)*y(ant)
        phi(a10,ant)=-delta*eff(a138)*y(4)
        phi(a10,a9)=-delta*l(a9)
        phi(a10,a10)=1.d0-delta*(-(y(3)*(l(a138)+eff(a9))+l(a10)+l(a89)))
        phi(a10,a13)=-delta*eff(a10)*y(3)
        phi(a10,a14)=-delta*eff(a89)*y(2)

        phi(a11,2)=-delta*(-eff(a88)*y(a11))
        phi(a11,3)=-delta*(eff(a11)*y(ant)-eff(a15)*y(a11))
        phi(a11,ant)=-delta*eff(a11)*y(3)
        phi(a11,a11)=1.d0-delta*(-(l(a11)+eff(a88)*y(2)+eff(a15)*y(3)))
        phi(a11,a14)=-delta*l(a88)
        phi(a11,a15)=-delta*l(a15)

        phi(a12,2)=-delta*(eff(a90)*y(a15)-y(a12)*(eff(a91)+l(a148)))
        phi(a12,3)=-delta*eff(a12)*y(a14)
        phi(a12,4)=-delta*eff(a148)*y(ant)
        phi(a12,ant)=-delta*eff(a148)*y(4)
        phi(a12,a9)=-delta*l(a91)
        phi(a12,a12)=1.d0-delta*(-(l(a12)+l(a90)+y(2)*(eff(a91)+l(a148))))
        phi(a12,a14)=-delta*eff(a12)*y(3)
        phi(a12,a15)=-delta*eff(a90)*y(2)

        phi(a13,2)=-delta*eff(a13)*y(a8)
        phi(a13,3)=-delta*(-eff(a10)*y(a13))
        phi(a13,a8)=-delta*eff(a13)*y(2)
        phi(a13,a10)=-delta*l(a10)
        phi(a13,a13)=1.d0-delta*(-(l(a13)+eff(a10)*y(3)))

        phi(a14,2)=-delta*(eff(a88)*y(a11)-eff(a89)*y(a14))
        phi(a14,3)=-delta*(eff(a14)*y(a8)-eff(a12)*y(a14))
        phi(a14,a8)=-delta*eff(a14)*y(3)
        phi(a14,a10)=-delta*l(a89)
        phi(a14,a11)=-delta*eff(a88)*y(2)
        phi(a14,a12)=-delta*l(a12)
        phi(a14,a14)=1.d0-delta*(-(l(a14)+l(a88)+eff(a89)*y(2)+eff(a12)*y(3)))


        phi(a15,2)=-delta*(-eff(a90)*y(a15))
        phi(a15,3)=-delta*eff(a15)*y(a11)
        phi(a15,a11)=-delta*eff(a15)*y(3)
        phi(a15,a12)=-delta*l(a90)
        phi(a15,a15)=1.d0-delta*(-(l(a15)+eff(a90)*y(2)))

        phi(4,2)=phi(4,2)-delta*l(a148)*y(a12)
        phi(4,3)=phi(4,3)-delta*l(a138)*y(a10)
        phi(4,4)=phi(4,4)-delta*(-y(ant)*(eff(a128)+eff(a138)+eff(a148)))
        phi(4,ant)=phi(4,ant)-delta*(-y(4)*(eff(a128)+eff(a138)+eff(a148)))
        phi(4,a9)=-delta*l(a128)
        phi(4,a10)=-delta*l(a138)*y(3)
        phi(4,a12)=-delta*l(a148)*y(2)

        phi(2,2)=phi(2,2)+delta*((eff(a8)*y(ant)+eff(a13)*y(a8)&
        &        +eff(a88)*y(a11)+eff(a90)*y(a15)+eff(a89)*y(a14)&
        &        +y(a12)*(l(a148)+eff(a91))))
        phi(2,4)=phi(2,4)-delta*eff(a148)*y(ant)
        phi(2,ant)=phi(2,ant)+delta*(eff(a8)*y(2)-eff(a148)*y(4))
        phi(2,a8)=-delta*(l(a8)-eff(a13)*y(2))
        phi(2,a9)=-delta*l(a91)
        phi(2,a10)=-delta*l(a89)
        phi(2,a11)=-delta*(-eff(a88)*y(2))
        phi(2,a12)=-delta*(l(a90)-y(2)*(l(a148)+eff(a91)))
        phi(2,a13)=-delta*l(a13)
        phi(2,a14)=-delta*(l(a88)-eff(a89)*y(2))
        phi(2,a15)=-delta*(-eff(a90)*y(2))

        phi(3,3)=phi(3,3)+delta*((eff(a11)*y(ant)+eff(a14)*y(a8)&
        &        +eff(a15)*y(a11)+eff(a10)*y(a13)+eff(a12)*y(a14)&
        &        +y(a10)*(eff(a9)+l(a138))))
        phi(3,4)=phi(3,4)-delta*eff(a138)*y(ant)
        phi(3,ant)=phi(3,ant)+delta*(eff(a11)*y(3)-eff(a138)*y(4))
        phi(3,a8)=-delta*(-eff(a14)*y(3))
        phi(3,a9)=-delta*l(a9)
        phi(3,a10)=-delta*(l(a10)-y(3)*(eff(a9)+l(a138)))
        phi(3,a11)=-delta*(l(a11)-eff(a15)*y(3))
        phi(3,a12)=-delta*l(a12)
        phi(3,a13)=-delta*(-eff(a10)*y(3))
        phi(3,a14)=-delta*(l(a14)-eff(a12)*y(3))
        phi(3,a15)=-delta*l(a15)

        ant=a9
        if(ant.eq.81)exit

        phi(a9,2)=phi(a9,2)+delta*(eff(a16)*y(a9)-l(a149)*y(a20))
        phi(a9,3)=phi(a9,3)+delta*(eff(a19)*y(a9)-l(a139)*y(a18))
        phi(a9,4)=phi(a9,4)+delta*(y(a9)*(eff(a129)+eff(a139)+eff(a149)))
        phi(a9,a9)=phi(a9,a9)+delta*((eff(a16)*y(2)+eff(a19)*y(3)&
        &          +y(4)*(eff(a129)+eff(a139)+eff(a149))))
        phi(a9,a16)=-delta*l(a16)
        phi(a9,a17)=-delta*l(a129)
        phi(a9,a18)=-delta*l(a139)*y(3)
        phi(a9,a19)=-delta*l(a19)
        phi(a9,a20)=-delta*l(a149)*y(2)

        a8=a16
        a9=a17
        a10=a18
        a11=a19
        a12=a20
        a13=a13+8
        a14=a14+8
        a15=a15+8
        a16=a16+8
        a17=a17+8
        a18=a18+8
        a19=a19+8
        a20=a20+8
        a88=a88+4
        a89=a89+4
        a90=a90+4
        a91=a91+4
        a128=a129
        a129=a129+1
        a138=a139
        a139=a139+1
        a148=a149
        a149=a149+1
      end do
      !Electron and positron captures ************************************
      phi(2,2)=phi(2,2)+1.d0+delta*effe
      phi(2,3)=-delta*effp
      phi(3,2)=-delta*effe
      phi(3,3)=phi(3,3)+1.d0+delta*effp
      phi(2,iso)=-delta*(deffedYe*y(3)-deffpdYe*y(2))
      phi(3,iso)=-delta*(deffedYe*y(2)-deffpdYe*y(3))
      phi(iso,2)=delta*effe
      phi(iso,iso)=1.d0+delta*deffedYe*y(2)
      !*******************************************************************

      phi(2,2)=phi(2,2)+delta*y(88)*l(158)+delta*y(5)*eff(160)+&
              & delta*y(6)*l(161)
      phi(2,4)=phi(2,4)-delta*y(79)*eff(158)-delta*y(90)*eff(161)
      phi(2,5)=-delta*cp*eff(1)*y(5)+delta*y(2)*eff(160)
      phi(2,6)=-delta*op*eff(3)*y(6)+delta*y(2)*l(161)
      phi(2,79)=phi(2,79)-delta*y(4)*eff(158)
      phi(2,88)=delta*y(2)*l(158)
      phi(2,90)=-delta*l(160)-delta*y(4)*eff(161)
      phi(3,3)=phi(3,3)+delta*y(88)*eff(159)
      phi(3,88)=phi(3,88)+delta*y(3)*eff(159)
      phi(3,89)=-delta*l(159)
      phi(4,2)=phi(4,2)-delta*y(88)*l(158)-delta*y(6)*l(161)
      phi(4,4)=phi(4,4)+ &
        & delta*((1.5d0*eff(5)*y(4)**2.d0+eff(6)*y(5)+eff(7)*y(6)))
      phi(4,4)=phi(4,4)+delta*y(79)*eff(158)+delta*y(90)*eff(161)
      phi(4,5)=-delta*(3.d0*l(5)-eff(6)*y(4)+ca*eff(1)*y(5)+eff(2)*y(6))
      phi(4,6)=-delta*(l(6)-eff(7)*y(4)+eff(2)*y(5)+oa*eff(3)*y(6))-delta*y(2)*l(161)
      phi(4,7)=phi(4,7)-delta*l(7)
      phi(4,79)=phi(4,79)+delta*y(4)*eff(158)
      phi(4,88)=-delta*y(2)*l(158)
      phi(4,90)=delta*y(4)*eff(161)

      phi(5,2)=phi(5,2)+delta*y(5)*eff(160)
      phi(5,4)=-delta*(.5d0*eff(5)*y(4)**2.d0-eff(6)*y(5))
      phi(5,5)=1.d0-delta*(-(l(5)+eff(6)*y(4)+2.d0*eff(1)*y(5)+eff(2)*y(6)))+&
                    delta*y(2)*eff(160)
      phi(5,6)=-delta*(l(6)-eff(2)*y(5)) 
      phi(5,90)=-delta*l(160)

      phi(6,2)=phi(6,2)+delta*y(6)*l(161)
      phi(6,4)=-delta*(eff(6)*y(5)-eff(7)*y(6))-delta*y(90)*eff(161)
      phi(6,5)=-delta*(eff(6)*y(4)-eff(2)*y(6))
      phi(6,6)=1.d0-delta*(-(l(6)+eff(7)*y(4)+eff(2)*y(5)+2.d0*eff(3)*y(6)))+&
              &  delta*y(2)*l(161)
      phi(6,7)=-delta*l(7)
      phi(6,90)=-delta*y(4)*eff(161)

      phi(7,2)=-delta*(l(148)*y(12)-eff(8)*y(7))
      phi(7,3)=-delta*(l(138)*y(10)-eff(11)*y(7))
      phi(7,4)=-delta*(eff(7)*y(6)-y(7)*(eff(128)+eff(138)+eff(148)))
      phi(7,5)=-delta*ca*eff(1)*y(5)
      phi(7,6)=-delta*eff(7)*y(4)
      phi(7,7)=1.d0-delta*(-(l(7)+y(2)*eff(8)+y(3)*eff(11)&
      &        +y(4)*(eff(128)+eff(138)+eff(148))))
      phi(7,8)=-delta*l(8)
      phi(7,9)=-delta*l(128)
      phi(7,10)=-delta*l(138)*y(3)
      phi(7,11)=-delta*l(11)
      phi(7,12)=-delta*l(148)*y(2)

      phi(9,5)=-delta*eff(2)*y(6)
      phi(9,6)=-delta*eff(2)*y(5)

      phi(12,5)=-delta*cp*eff(1)*y(5)

      phi(17,6)=-delta*oa*eff(3)*y(6)

      phi(28,6)=-delta*op*eff(3)*y(6)
      
      phi(79,2)=phi(79,2)-delta*y(88)*l(158)
      phi(79,4)=phi(79,4)+delta*y(79)*eff(158)
      phi(79,79)=phi(79,79)+delta*y(4)*eff(158)
      phi(79,88)=-delta*y(2)*l(158)
      phi(88,2)=delta*y(88)*l(158)
      phi(88,4)=-delta*y(79)*eff(158)
      phi(88,3)= delta*y(88)*eff(159)
      phi(88,79)=-delta*y(4)*eff(158)
      phi(88,88)=1.d0+delta*y(2)*l(158)+delta*y(3)*eff(159) 
      phi(88,89)=-delta*l(159)
      phi(89,3)=-delta*y(88)*eff(159)
      phi(89,88)=-delta*y(3)*eff(159)
      phi(89,89)=1.d0+delta*l(159)
      phi(90,2)=-delta*y(5)*eff(160)-delta*y(6)*l(161) 
      phi(90,4)= delta*y(90)*eff(161)
      phi(90,6)=-delta*y(2)*l(161)
      phi(90,5)=-delta*y(2)*eff(160)
      phi(90,90)=1.d0+delta*l(160)+delta*y(4)*eff(161)
      phi(4,4)=phi(4,4)+1.d0

      ! Energy term includes ideal gas correction.
      kbt=kb*temp
      nakbt=kbt*na
      do i=2,iso
         phi(i,1)=-delta*fp(i)
         phi(1,i)=-be(i)*MeVperparticle2erg+nakbt*1.5d0
      enddo
      phi(1,iso)=phi(1,iso)-nakbt*1.5d0+dUedYe+theta*delta*(dEneutrdYe*y(2)+dEaneutrdYe*y(3))
      phi(1,2)=phi(1,2)+theta*Eneutr*delta
      phi(1,3)=phi(1,3)+theta*Eaneutr*delta

      ! Corrections to the energy term due to coulombian energy.
      ! Always active because it's more coherent due to the fact that
      ! these corrections are included in the calculation of Cv.
      ne=rho*na*ye
      ae=(3.d0/(4.d0*pi*ne))**third
      gam=e2/(kbt*ae)
      do i=2,iso
         gamma=gam*z(i)**(5.d0/3.d0)
         if(gamma.gt.1.d0) then
            funcion=ggt1(gamma)
         else if(gamma.le.1.d0) then
            funcion=glt1(gamma)
         endif
         phi(1,i)=phi(1,i)+nakbt*funcion    
      enddo
      RETURN
      END

! *******************************************************************

      FUNCTION ggt1(gamma)
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(in)::gamma
      DOUBLE PRECISION a1,b1,c1,d1,sqroot2gamma,ggt1
      a1=-.898004d0
      b1=.96786d0
      c1=.220703d0
      d1=-.86097d0
      sqroot2gamma=sqrt(sqrt(gamma))
      ggt1=a1*gamma+b1*gamma**(.25d0)+c1/sqroot2gamma+d1
      return
      end function ggt1

! *******************************************************************

      FUNCTION glt1(gamma)
      implicit none
      DOUBLE PRECISION,INTENT(in)::gamma
      DOUBLE PRECISION a1,b1,c1,glt1
      a1=-.5d0*sqrt(3.d0)
      b1=.29561d0
      c1=1.9885d0
      glt1=a1*gamma**(1.5d0)+b1*gamma**c1
      return
      end function glt1


! *******************************************************************

      SUBROUTINE effect(temp)
      USE nuclear90_module
      IMPLICIT NONE

      DOUBLE PRECISION,INTENT(in)::temp

      DOUBLE PRECISION t9,t9r,t92,t93,t95,t912,t913,t923,t932,t943,t953
      DOUBLE PRECISION t9m1,t9m13,t9m23,t9m32,t9rm1,t9r32
      DOUBLE PRECISION t9a,t9a2,t9a13,t9a23,t9a56
      DOUBLE PRECISION r2abe,rbeac,r3a,rev,rg3a,r24,rc28,r1216,r32
      DOUBLE PRECISION rcag,roga,roag,rnega

      t9r=temp*1.0d-09
      t9=min(10.d00,t9r)
      t92=t9*t9
      t93=t92*t9
      t95=t92*t93
      t912=sqrt(t9)
      t913=t9**(third)
      t923=t913*t913
      t932=t9*t912
      t943=t9*t913
      t953=t9*t923
      t9m1=1.d0/t9
      t9m13=1.d0/t913
      t9m23=1.d0/t923
      t9m32=1.d0/t932
      t9rm1=1.d0/t9r
      t9r32=t9*sqrt(t9r)
      do i=1,rates
        eff(i)=0.d0
      enddo
      do i=1,rates
        l(i)=0.d0
      enddo

! **
! ** triple alpha reaction rate
! **
      r2abe=(7.40d+05*t9m32)*exp(-1.0663d0*t9m1)&
     &      +4.164d+09*t9m23*exp(-13.49d0*t9m13-t92/0.009604d0)&
     &      *(1.d0+0.031d0*t913+8.009d0*t923+1.732d0*t9+49.883d0&
     &      *t943+27.426d0*t953)
      rbeac=(130.d0*t9m32)*exp(-3.3364d0*t9m1)&
     &      +2.510d+07*t9m23*exp(-23.57d0*t9m13-t92/0.055225d0)&
     &      *(1.d0+0.018d0*t913+5.249d0*t923+0.650d0*t9+19.176d0&
     &      *t943+6.034d0*t953)
      if(t9.gt.0.08d0) then
            r3a=2.90d-16*(r2abe*rbeac)&
     &          +0.1d0*1.35d-07*t9m32*exp(-24.811d0*t9m1)
      else
            r3a=2.90d-16*(r2abe*rbeac)&
     &          *(0.01d0+0.2d0*(1.d0&
     &          +4.d0*exp(-(0.025d0*t9m1)**3.263d0))&
     &          /(1.d0+4.d0*exp(-(t9/0.025d0)**9.227d0)))&
     &          +0.1d0*1.35d-07*t9m32*exp(-24.811d0*t9m1)
      endif
      eff(5)=r3a
!      rev=2.d20*exp(-84.424035d0*t9m1)
!This one has a different Q than in net14 ?!?!
      rev=2.d20*exp(-84.426375d0*t9m1)


! ** 12c photodesintegration
      rg3a=rev*(t93)*r3a
      l(5)=rg3a

! **
! ** 12c+12c reaction rate
! **
      t9a=t9/(1.d0+0.0396d0*t9)
      t9a13=t9a**third
      t9a56=t9a**fsix
      r24=4.27d+26*t9a56*t9m32*exp(-84.165d0/t9a13-2.12d-03*t93)
      eff(1)=r24

! **
! ** c12+o16 reaction rate
! **
      if(t9.ge..5d0) then
            t9a=t9/(1.d0+0.055d0*t9)
            t9a2=t9a*t9a
            t9a13=t9a**third
            t9a23=t9a13*t9a13
            t9a56=t9a**fsix
            rc28=1.72d+31*t9a56*t9m32*exp(-106.594d0/t9a13)&
     &          /(exp(-0.18d0*t9a2)+1.06d-03*exp(2.562d0*t9a23))
            r1216=rc28
      else
            r1216=0.d0
      endif
      eff(2)=r1216

! **
! ** 16o+16o reaction rate
! **
      r32=7.10d+36*t9m23*exp(-135.93d0*t9m13-0.629d0*t923-0.445d0*t943&
     &    +0.0103d0*t92)
      eff(3)=r32

! **
! ** 12c(a,g)16o
! **
      rcag=(1.04d+08/(t92*(1.d0+0.0489d0*t9m23)**2.d0)&
     &     *exp(-32.120d0*t9m13-t92/12.222016d0)&
     &     +1.76d+08/(t92*(1.d0+0.2654d0*t9m23)**2.d0)&
     &     *exp(-32.120d0*t9m13)&
     &     +1.25d+03*t9m32*exp(-27.499d0*t9m1)+1.43d-02*t95&
     &     *exp(-15.541d0*t9m1))
      eff(6)=rcag
      roga=rcag*5.13d+10*t9r32*exp(-83.11501d0*t9rm1)
      l(6)=roga

! **
! ** 16o(a,g)20ne
! **
      roag=(9.37d+09*t9m23*exp(-39.757d0*t9m13-t92/2.515396d0)&
     &     +62.1d0*t9m32*exp(-10.297d0*t9m1)+538.d0*t9m32&
     &     *exp(-12.226d0*t9m1)+13.d0*t92*exp(-20.093d0*t9m1))
      eff(7)=roag

! ** 20Ne photodesintegration
      rnega=roag*5.65d+10*t9r32*exp(-54.93807d0*t9rm1)
      l(7)=rnega

      RETURN
      END

! *******************************************************************
! *** Analytical calculation of  deff/dtemp

      SUBROUTINE deffcalc(temp)

      USE nuclear90_module

      IMPLICIT NONE

      DOUBLE PRECISION,INTENT(in)::temp
      DOUBLE PRECISION t9,t92,t93,t94,t95,t912,t913,t923,t932,t943,t952,t953
      DOUBLE PRECISION t9i,t9i2,t9i3,t9i12,t9i13,t9i23,t9i32,t9i43,t9i52,t9i53
      DOUBLE PRECISION vA,vB,vC,vD,vE,vF,vG,val,r2abe,rbeac,vA56
      DOUBLE PRECISION dvA,dvB,dvC,dvD,dvE,dvF,dvG,dval,dr2abe,drbeac

      t9=temp*1.d-9
      t92=t9*t9
      t93=t92*t9
      t94=t93*t9
      t95=t94*t9
      t912=sqrt(t9)
      t913=t9**third
      t923=t913*t913
      t932=t9*t912
      t943=t9*t913
      t952=t9*t932
      t953=t9*t923
      t9i=1.d0/t9
      t9i2=t9i*t9i
      t9i3=t9i2*t9i
      t9i12=1.d0/t912
      t9i13=1.d0/t913
      t9i23=1.d0/t923
      t9i32=1.d0/t932
      t9i43=1.d0/t943
      t9i52=1.d0/t952
      t9i53=1.d0/t953
      do i=1,rates
        deff(i)=0.d0
      enddo
      do i=1,rates
        dl(i)=0.d0
      enddo

! *** dl(4) and deff(5). Alpha
      vA=-24.811d0*t9i
      vB=-1.0663d0*t9i
      vC=-13.49d0*t9i13-t92*104.123282d0
      vD=1.+.031d0*t913+8.009d0*t923+1.732d0*t9+49.883d0*t943+27.426d0*t953
      vE=-3.3364d0*t9i
      vF=-23.57d0*t9i13-18.10774106d0*t92
      vG=1.+.018d0*t913+5.249d0*t923+.650d0*t9+19.176d0*t943+6.034d0*t953
      dvA=24.811d0*t9i2
      dvB=1.0663d0*t9i2
      dvC=t9i43*13.49d0/3.-208.246564d0*t9
      dvD=(t9i23*.031d0+t9i13*16.018d0+5.196d0+t913*199.532d0+t923*137.13d0)/3.d0
      dvE=3.3364d0*t9i2
      dvF=t9i43*23.57d0/3.-36.21548212d0*t9
      dvG=(t9i23*.018d0+t9i13*10.498d0+1.950d0+t913*76.704d0+t923*30.17d0)/3.d0

      r2abe=7.4d5*t9i32*exp(vB)+4.164d9*t9i23*vD*exp(vC)
      rbeac=130.*t9i32*exp(vE)+2.510d7*t9i23*vG*exp(vF)
      dr2abe=exp(vB)*(-1.11d6*t9i52+7.4d5*t9i32*dvB)+&
     &       4.164d9*exp(vC)*(-t9i53*vD*2.d0/3.d0+t9i23*dvC*vD+t9i23*dvD)
      drbeac=exp(vE)*(-195.d0*t9i52+130.d0*t9i32*dvE)+&
     &       2.510d7*exp(vF)*(-2.d0*t9i53*vG/3.d0+t9i23*dvF*vG+t9i23*dvG)
      deff(5)=2.90d-16*(dr2abe*rbeac+r2abe*drbeac)+1.35d-8*exp(vA)*&
     &        (-1.5d0*t9i52+t9i32*dvA)


! *** dl(5) and deff(6). 12C
      vA=-84.424035d0*t9i
      vB=(1.d0+.0489d0*t9i23)**2.d0
      vC=-32.120d0*t9i13-(t9/3.496d0)**2.d0
      vD=(1.d0+.2654d0*t9i23)**2.d0
      vE=-32.120d0*t9i13
      vF=-27.499d0*t9i
      vG=-15.541d0*t9i
      dvA=84.424035d0*t9i2
      dvB=-(2.d0+.0978d0*t9i23)*(.0326d0*t9i53)
      dvC=32.120d0*t9i43/3.d0-2.d0*t9/(3.496d0)**2.d0
      dvD=-(2.d0+.5308d0*t9i23)*(.5308d0*t9i53/3.d0)
      dvE=32.120d0*t9i43/3.d0
      dvF=27.499d0*t9i2
      dvG=15.541d0*t9i2

      dl(5)=2.00d20*exp(vA)*t93*(dvA*eff(5)+3.d0*t9i*eff(5)+deff(5))
      deff(6)=1.04d8*exp(vC)*t9i2*(-2.d0*t9i+dvC-dvB/vB)/vB+&
     &        1.76d8*exp(vE)*t9i2*(-2.d0*t9i+dvE-dvD/vD)/vD+&
     &        1.25d3*exp(vF)*(-1.5d0*t9i52+dvF*t9i32)+&
     &        1.43d-2*exp(vG)*(5.d0*t94+dvG*t95)


! *** dl(6) and deff(7). 16O
      vA=-83.11501d0*t9i
      vB=-39.757d0*t9i13-(t9/1.586d0)**2.d0
      vC=-10.297d0*t9i
      vD=-12.226d0*t9i
      vE=-20.093d0*t9i
      dvA=83.11501d0*t9i2
      dvB=39.757d0*t9i43/3.d0-2.d0*t9/(1.586d0)**2.d0
      dvC=10.297d0*t9i2
      dvD=12.226d0*t9i2
      dvE=20.093d0*t9i2

      dl(6)=5.13d10*exp(vA)*(deff(6)*t932+eff(6)*1.5d0*t912+eff(6)*t932*dvA)
      deff(7)=9.37d9*exp(vB)*(-t9i53*2.d0/3.d0+t9i23*dvB)+&
     &        62.1d0*exp(vC)*(-1.5d0*t9i52+t9i32*dvC)+&
     &        538.d0*exp(vD)*(-1.5d0*t9i52+t9i32*dvD)+&
     &        13.d0*exp(vE)*(2.d0*t9+t92*dvE)


! *** dl(7). 20Ne
      vA=-54.93807d0*t9i
      dvA=54.93807d0*t9i2
      dl(7)=5.65d10*exp(vA)*(deff(7)*t932+1.5d0*eff(7)*t912+eff(7)*t932*dvA)


! *** deff(1) 12C+12C
      vA=t9/(1.d0+.0396d0*t9)
      vA56=vA**fsix
      vB=-84.165d0*vA**(-third)-2.12d-3*t93

      dvA=vA*vA*t9i2
      dvB=28.055d0*dvA*vA**(-4.d0/3.d0)-6.36d-3*t92

      deff(1)=4.27d26*t9i32*exp(vB)*(vA**(-1.d0/6.d0)*dvA*fsix&
     &        -1.5d0*vA56*t9i+vA56*dvB)

! *** deff(2) 12C+16O
      if(t9.ge..5d0) then
        vA=t9/(1.d0+.055d0*t9)
        vA56=vA**(5.d0/6.d0)
        vB=-106.594d0*vA**(-1.d0/3.d0)
        vC=-.18d0*vA*vA
        vD=2.562d0*vA**(2.d0/3.d0)
        val=exp(vC)+1.06d-3*exp(vD)

        dvA=vA*vA*t9i2
        dvB=106.594d0*vA**(-4.d0/3.d0)*dvA/3.d0
        dvC=-.36d0*vA*dvA
        dvD=1.708d0*dvA*vA**(-third)
        dval=dvC*exp(vC)+1.06d-3*dvD*exp(vD)


        deff(2)=1.72d31*t9i32*(vA**(-1.d0/6.d0)*dvA*fsix&
     &         -1.5d0*vA56*t9i&
     &         +vA56*(dvB-dval/val))*exp(vB)/val
      else
        deff(2)=0.d0
      endif


! *** deff(3) 16O+16O
      vA=-135.93d0*t9i13-.629d0*t923-.445d0*t943+.0103d0*t92

      dvA=45.31d0*t9i43-.629d0*t9i13*2.d0/3.d0-.445d0*t913*4.d0/3.d0+.0206d0*t9

      deff(3)=7.10d36*exp(vA)*t9i23*(-t9i*2.d0/3.d0+dvA)

      do i=1,7
        deff(i)=deff(i)/1.d9
      enddo
      do i=1,7
        dl(i)=dl(i)/1.d9
      enddo

      return
    end subroutine deffcalc


! ******************************************************************************
      SUBROUTINE chempot(temp,rho)

!     CALCULATION OF THE CHEMICAL POTENTIALS
!     OUTPUT: mukbt & deltamukbt through nuclear90_module

      USE nuclear90_module

      implicit none

      DOUBLE PRECISION,INTENT(in)::temp,rho

      DOUBLE PRECISION ne,kbt,ae,gam
      DOUBLE PRECISION a1,b1,c1,d1,e1
      DOUBLE PRECISION gamp,sqrootgamp,sqroot2gamp

      integer ini1,ini2,fin1,fin2

! In these ones final2 are neutrons because they have mukbt=0, just the same as photons.

      ne=rho*na*ye
      kbt=kb*temp
      ae=(3.d0/(4.d0*pi*ne))**third
      gam=e2/(kbt*ae)

      a1=-.898004d0
      b1=.96786d0
      c1=.220703d0
      d1=-.86097d0
      e1=2.520058332d0

      do i=1,91
         mukbt(i)=0.d0
         deltamukbt(i)=0.d0
      enddo
      do i=92,rates
         deltamukbt(i)=0.d0
      enddo

      do i=2,91
         if(i.eq.3)cycle
         gamp=gam*z(i)**(5.d0/3.d0)
         sqrootgamp=sqrt(gamp)
         sqroot2gamp=sqrt(sqrootgamp)
         if(gamp.le.1.d0) then
            mukbt(i)=-(1.d0/sqrt(3.d0))*gamp*sqrootgamp+&
            &  (gamp**1.9885d0)*.29561d0/1.9885d0
         else
            mukbt(i)=a1*gamp+4.d0*(b1*sqroot2gamp-c1/sqroot2gamp)&
            &  +d1*log(gamp)-e1
         endif
      enddo

! mu for neutrons must be zero
      mukbt(3)=0.d0
      do i=2,rates
         if(i.eq.3.or.i.eq.4)cycle
         ini1=target1(i)
         ini2=target2(i)
         fin1=final1(i)
         fin2=final2(i)
         deltamukbt(i)=mukbt(ini1)+mukbt(ini2)-(mukbt(fin1)+mukbt(fin2))
      enddo

! Triple alpha correction. We must sum a mukbt(alpha) because there are three alfas.
! The intermediate step of Berillium cancels, that's why is not taken into account.
      deltamukbt(5)=deltamukbt(5)+mukbt(4)

! Additional deltamukbt for the two extra channels of heavy ion reactions.
! We use 12C+12C=24Mg+G y 16O+16O=32S+G
! We use the gamma channels to calculate the deltamukbt of both channels
! because in fact the outcome of the heavy ion reactions is an excited nucleus
! that decays spontaneously to one of the channels used in the reactions.
! Nevertheless, this decay is not a nuclear reaction so it is not completely
! correct to use this coulombian corrections.
! In net14 it is handled as if it was truly the reaction 12C+12C->20Ne+alpha,
! for example, but it is more correct in this way.
! In any case, these corrections are not absolutely perfect anyway...
      deltamukbt(1)=2.d0*mukbt(5)-mukbt(9)
      deltamukbt(3)=2.d0*mukbt(6)-mukbt(25)
      return
      end


      ! *******************************************************************

            subroutine eigen(nisp,c,x)

              use nuclear90_module

              implicit none
      !-----------------------------------------------------------------------
      !c    Prantzos, Arnould & Arcoragi  (ApJ 315,209 - 1987)
      !c    inversion of the matrix (special sparse form)
      !-----------------------------------------------------------------------
      !c    solves the equation: phi2*x=c
      !c    warning: for computation convenience the array phi2 that is fed through
      !c             common/mat/ allocates an additional column, which is equal
      !c             to vector c
      !-----------------------------------------------------------------------
      !c nisp == dimension of the true matrix a (nisp x nisp), and number of
      !c         columns of array a in common/mat/. INPUT
      !c c == vector containing the independent terms in the equation system. INPUT
      !c x == vector containing the solution. OUTPUT
      !c
      !c phi2 == array of dimensions (nisp x nisp+1) containing matrix a and vector c.
      !c         INPUT through common/mat/
      !c
      !c nds == width of the upper diagonal. INPUT through common/ceign/
      !c ndi == width of the lower diagonal. INPUT through common/ceign/
      !c ijd == number of columns in the left box + 1. INPUT through common/ceign/
      !-----------------------------------------------------------------------

              DOUBLE PRECISION,dimension(iso)::c,x
              integer,intent(in)::nisp

              integer nn,idel,jdel,jbal,jmax,imax,jmin
              DOUBLE PRECISION divj,som
              DOUBLE PRECISION,dimension(iso)::div(iso)

              nn=nisp
              idel=ijd-1
              jdel=ijd-1
      !
      !---beginning of the gaussian elimination procedure
      !-(1)-elimination of the lower diagonals
      !
              do 1000 jbal=ijd,nn-1
                 div(jbal)=-1./phi(jbal,jbal)
                 jmax=jbal+ndi
                 if(jmax.gt.nn)jmax=nn
                 imax=jbal+nds
                 if(imax.gt.nn)imax=nn
      !*vdir ignore recrdeps
                 do 1001 j=jbal+1,jmax
                    divj=div(jbal)*phi(j,jbal)
                    if(divj.eq.0.)goto 1001
                    do 10 i=1,ijd-1
                       !if(divj.eq.0.)goto 10
                       phi(j,i)=divj*phi(jbal,i)+phi(j,i)
      10            enddo
                    do 20 i=jbal+1,imax
                       !if(divj.eq.0.)goto 20
                       phi(j,i)=divj*phi(jbal,i)+phi(j,i)
      20            enddo
      1001       enddo
      !*vdir ignore recrdeps
                 do 1002 j=jbal+1,jmax
                    if(phi(j,jbal).eq.0.)goto 1002
                    c(j)=div(jbal)*phi(j,jbal)*c(jbal)+c(j)
      1002       enddo
      1000    enddo
              div(nn)=-1./phi(nn,nn)
      !
      !-(2)-elimination of the upper diagonals and of the horizontal band
      !
              do 2000 jbal=nn,ijd+1,-1
                 jmin=jbal-ndi
                 if(jmin.lt.ijd)jmin=ijd
                 do 200 j=jmin,jbal-1
                    if(phi(j,jbal).eq.0.)goto 200
                    divj=div(jbal)*phi(j,jbal)
                    do 30 i=1,idel
                       phi(j,i)=divj*phi(jbal,i)+phi(j,i)
      30            enddo
                    c(j)=divj*c(jbal)+c(j)
      200        enddo
                 do 300 j=1,jdel
                    if(phi(j,jbal).eq.0.)goto 300
                    divj=div(jbal)*phi(j,jbal)
                    do 40 i=1,idel
                       phi(j,i)=divj*phi(jbal,i)+phi(j,i)
      40            enddo
                    c(j)=divj*c(jbal)+c(j)
      300        enddo
      2000    enddo
              do 400 j=1,jdel
                 if(phi(j,ijd).eq.0.)goto 400
                 divj=div(ijd)*phi(j,ijd)
                 do 50 i=1,idel
                    phi(j,i)=divj*phi(ijd,i)+phi(j,i)
      50         enddo
                 c(j)=divj*c(ijd)+c(j)
      400     enddo
      !
      !-(3)-gaussian elimination of the upper left square matrix
      !
              do 3000 jbal=1,ijd-2
                 div(jbal)=-1./phi(jbal,jbal)
                 do 3001 j=jbal+1,ijd-1
                    if(phi(j,jbal).eq.0.)goto 3000
                    divj=div(jbal)*phi(j,jbal)
                    do 60 i=jbal+1,ijd-1
                       phi(j,i)=divj*phi(jbal,i)+phi(j,i)
      60            enddo
                    c(j)=divj*c(jbal)+c(j)
      3001       enddo
      3000    enddo
              div(ijd-1)=-1./phi(ijd-1,ijd-1)
              x(ijd-1)=-div(ijd-1)*c(ijd-1)
              do 4000 jbal=ijd-2,1,-1
                 som=0.
                 do 70 i=ijd-1,jbal+1,-1
                    som=som+phi(jbal,i)*x(i)
      70         enddo
                 x(jbal)=div(jbal)*(som-c(jbal))
      4000    enddo
              som=0.
              do 80 i=1,idel
                 som=som+phi(nn,i)*x(i)
      80      enddo
              x(nn)=div(nn)*(som-c(nn))
              do 5000 jbal=nn-1,ijd,-1
                 som=0.
                 do 90 i=1,idel
                    som=som+phi(jbal,i)*x(i)
      90         enddo
                 x(jbal)=div(jbal)*(som-c(jbal))
      5000    enddo
              return
            end subroutine eigen

      ! **************************************************************

