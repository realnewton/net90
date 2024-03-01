!**********************************************************************
! Program that calls the nuclear network net90 for testing.
! See R.M. Cabezon et al. ApJ 151 (2004) 
! for further information.
!**********************************************************************
!
!************* Input parameters. Defined in parameters.in *************
!
! All relevant parameters are defined in parameters.in. This contains:
!
! - Mass fractions of each species (DP)
! - Initial_temp  (D): initial temperature in K (DP)
! - Initial_dens  (D): initial density in g/cm^3 (DP)
! - Screening     (L): Activates the screening corrections (L)
! - Adiabatic_exp (L): (T) an adiabatic expansion follows after NSE
!                      (F) the process is isochoric
! - Isotherm_test (L): (T) temperature is kept constant via an 
!                          artificially high cv
!                      (F) temperature is updated by net90
! - electron_capt (L): Activates e-/e+ reactions
! - tabul_rates   (L): (T) e-/e+ rates interpolating on a table
!                      (F) e-/e+ rates analytically
! - Theta         (D): Chooses between explicit (0), implicit(1), or mixed.
!                      Default value is 0.7d0
! - C(c_p)Na_prob (D): Probability of 12C+12C->23Na+p branch
!                      Default value is 0.5d0
! - C(c_a)Ne_prob (D): Probability of 12C+12C->20Ne+4He branch
!                      Default value is 0.5d0
! - O(O_p)P_prob  (D): Probability of 16O+16O->31P+p branch
!                      Default value is 0.6d0
! - O(O_a)S_prob  (D): Probability of 16O+16O->28S+a branch
!                      Default value is 0.4d0
! - NR_max_iter   (I): Maximum number of Newton-Raphson iterations 
!                      in net90 per timestep. Default value is 10
! - NR_limit      (D): Convergence limit for NR iterations
!                      Default value is 1.d-7  
! - type_eos      (I): Selects the EOS
!                      (1) relativistic electrons (Nadyozhin 1974),
!                          ions (Bravo & Garcia-Senz 1992), and radiation.
!                      (2) helmholtz EOS (Cox & Giuli chapter 24;
!                          Timmes & Swesty ApJ 1999)
!
! Nuclear species included in the network and its corresponding index
! in the mass fraction array:
!  1-p,     2-n,     3-He4,   4-C12,   5-O16,   6-Ne20,
!  7-Na21,  8-Mg24,  9-Mg23, 10-Ne21, 11-Na23, 12-Mg22, 13-Na22, 14-Ne22,
! 15-Al25, 16-Si28, 17-Si27, 18-Mg25, 19-Al27, 20-Si26, 21-Al26, 22-Mg26,
! 23-P29 , 24-S32 , 25-S31,  26-Si29, 27-P31,  28-S30,  29-P30,  30-Si30,
! 31-Cl33, 32-Ar36, 33-Ar35, 34-S33,  35-Cl35, 36-Ar34, 37-Cl34, 38-S34,
! 39-K37 , 40-Ca40, 41-Ca39, 42-Ar37, 43-K39,  44-Ca38, 45-K38,  46-Ar38,
! 47-Sc41, 48-Ti44, 49-Ti43, 50-Ca41, 51-Sc43, 52-Ti42, 53-Sc42, 54-Ca42,
! 55-V45 , 56-Cr48, 57-Cr47, 58-Ti45, 59-V47,  60-Cr46, 61-V46,  62-Ti46,
! 63-Mn49, 64-Fe52, 65-Fe51, 66-Cr49, 67-Mn51, 68-Fe50, 69-Mn50, 70-Cr50,
! 71-Co53, 72-Ni56, 73-Ni55, 74-Fe53, 75-Co55, 76-Ni54, 77-Co54, 78-Fe54,
! 79-Cu57, 80-Zn60, 81-Zn59, 82-Ni57, 83-Cu59, 84-Zn58, 85-Cu58, 86-Ni58,
! 87-Co57, 88-Co58, 89-N13
!
! Ruben M. Cabezon & D. García-Senz (2023)
!**********************************************************************

  PROGRAM testnuclear

    USE nuclear90_module, only:iso,niso,cp,ca,op,oa,ye,a,z,&
              & NR_iter,NR_limit,electroncapture,tabulated

    IMPLICIT NONE

    DOUBLE PRECISION,PARAMETER:: initial_dt = 1.d-16   !Initial timestep
    DOUBLE PRECISION,PARAMETER:: gammam1    = 5.d0/3.d0-1.d0 
    DOUBLE PRECISION,PARAMETER:: maxtime    = 2.d0
    DOUBLE PRECISION,PARAMETER:: maxdelta   = 2.d-3 

    !Program ends only when reching T=1.d7 after expansion
    !Or after certain time (maxtime), if there is no expansion

    INTEGER, PARAMETER:: nsteps = 1000000000
    LOGICAL,PARAMETER :: noprint=.false.   !To avoid screen verbosity set to true

    INTEGER max_abound,type_eos
    INTEGER i,j

    LOGICAL adiab,isotherm,screen,isExpand,expt

    DOUBLE PRECISION temp,tempout,dens,delta,cv,v1,v2,sumtot
    DOUBLE PRECISION densold,deltadens,pres,Utotal,dUdYe
    DOUBLE PRECISION time,deltatemp,densini,t_hd,deltatime,t0
    DOUBLE PRECISION theta,deltad,deltat,dold,nucenergy,Exptime
    DOUBLE PRECISION abar,zbar,yi,mui
    CHARACTER*13 admy

    DOUBLE PRECISION,DIMENSION(niso) :: x,xout
    INTEGER,DIMENSION(:),ALLOCATABLE::k_iter


!************ Read parameters ************
    open(1,file='parameters.in')
 1  format(115(1x,es17.10))
    do i=1,niso
       read(1,*) admy,x(i)
    enddo
    read(1,*) admy,temp
    read(1,*) admy,densini
    read(1,*) admy,screen
    read(1,*) admy,adiab
    read(1,*) admy,isotherm
    read(1,*) admy,electroncapture
    read(1,*) admy,tabulated
    read(1,*) admy,theta
    read(1,*) admy,cp
    read(1,*) admy,ca
    read(1,*) admy,op
    read(1,*) admy,oa
    read(1,*) admy,NR_iter
    read(1,*) admy,NR_limit
    read(1,*) admy,type_eos
    close(1)
    allocate(k_iter(10000))
    k_iter=0

    
!****************** Set EOS ******************
    if(type_eos.eq.1) THEN
       print *,'EOS: Nadyozhin + IONS + Radiation'
    else if (type_eos.eq.2) THEN
       print *,'EOS: Helmholtz'
       call read_helm_table
    else
       stop 'Wrong type_eos'
    endif

    !read neutrino reaction rates
    if(tabulated) call read_rates_table

!************ Checking parameters ************
    if(sum(x).ne.1.d0) stop 'Wrong total composition. Sum(X)!=1'

!************ Estimation of expansion time ************
    if(adiab) then
       max_abound=maxloc(x,DIM=1)
       if(max_abound.eq.3) then  !He ignition test
          t_hd    = 446.d0/sqrt(densini)
          Exptime = 2.d0*t_hd 
       else if (max_abound.eq.4 .or. max_abound.eq.5) then !CO ignition test
           t_hd = 446.d0/sqrt(densini)
           Exptime=50.d0*t_hd
       else
          t_hd    = 446.d0/sqrt(densini)
          Exptime = t_hd 
       endif
    endif
    expt=.true.

!************ Initializations ************
    isExpand  = .false.         !Fluid is not expanding yet
    delta     = initial_dt      !Timestep
    time      = delta           !Physical time
    dens      = densini         !Initial density
    densold   = dens            !Previous density
    deltatemp = 0.d0            !Variation of temperature
    deltadens = 0.d0            !Variation of density
    ye        = 0.d0            !Electron number
    do i=1,niso
      ye = ye+x(i)/a(i+1)*z(i+1)
    enddo

    open(1,file='result.d')

    if(.not.noprint)print '(10(1x,a13))','Temp(K)','Dens(g/cm3)','n','X(4He)','X(12C)','X(16O)',&
                        & 'X(28Si)','X(56Ni)','X(58Ni)','X(59Ni)'
    do i=1,nsteps
         if(mod(i,10).eq.0.and..not.noprint) print '(10(1x,es13.6))', &
                      & temp,dens,x(2),x(3),x(4),x(5),x(16),x(72),x(86)

         !Call equation of state
         call eos(type_eos,dens,temp,x,cv,v1,v2,pres,Utotal,dUdYe)
         v1=v1*deltadens
         v2=v2*deltadens

         !Artificially increase cv for isothermal test so that temperature doesn't change
         if(isotherm.and..not.isExpand) cv=cv*1.d20

         !Call nuclear network
         call net90(temp,x,dens,delta,screen,cv,v1,v2,dUdYe,theta,&
                  & tempout,xout,sumtot,nucenergy,k_iter(i))                  

         !Timestep control
         deltat    = 1.d60
         deltad    = 1.d60
         deltatemp = tempout-temp
         dold      = delta
         if(deltatemp.ne.0.d0)  deltat = 0.01d0*temp*delta/abs(deltatemp)
         if(deltadens.ne.0.d0)  deltad = 0.02d0*dens*delta/abs(deltadens)
         delta     = min(deltat,deltad,maxdelta)
         if(delta.gt.1.25d0*dold) delta = 1.25d0*dold

         !Update quantities
         x       = xout
         temp    = tempout
         time    = time+delta
         densold = dens
         if(temp.ge.1.6d9 .and. expt) then 
            Exptime=time+Exptime
            expt=.false.
          endif

         !Change density if adiabatic expansion
         if(adiab) then
            if(time.gt.Exptime) then
               if(temp.lt.1.d7) stop 'Min. density reached'
               if(.not.isExpand) then
                  t0       = time
                  isExpand = .true.
                  delta    = 1.d-4
                  if(.not.noprint)print *,'Expansion: ',i
               endif
               deltatime   = time-t0
               dens        = densini*exp(-deltatime/(2.d0*t_hd))
            endif
         else
            if(time.gt.maxtime) stop 'Max. physical time reached'
         endif

         deltadens = dens-densold

         abar=0.d0
         zbar=0.d0
         yi=0.d0
         do j=1,niso
            abar=abar+xout(j)*a(j+1)
            zbar=zbar+xout(j)*z(j+1)
            yi=yi+xout(j)/a(j+1)
         enddo
         mui=1.d0/yi

         !Write output
         write(1,1) time,tempout,sumtot,dens,pres,&
                  & max(0.,nucenergy)/delta,min(0.,nucenergy)/delta,(xout(j),j=1,niso),ye,mui
      enddo
      close(1)

  end PROGRAM testnuclear

!**********************************************

  SUBROUTINE eos(type_eos,dens,temp,x,cv,v1,v2,pres,Utotal,dUdYe)

    USE nuclear90_module, only: niso,a,z,ye

    IMPLICIT NONE

    INTEGER i
    INTEGER, intent(in):: type_eos
    LOGICAL error

    DOUBLE PRECISION, intent(in):: dens,temp
    DOUBLE PRECISION, DIMENSION(niso),intent(in):: x

    DOUBLE PRECISION, intent(out):: cv,v1,v2,pres,Utotal,dUdYe
    DOUBLE PRECISION, DIMENSION(niso)::anum,znum

    DOUBLE PRECISION mue,pe,dpdte,dpdre,ue,dudte,dudre,dudt,abar,zbar
    DOUBLE PRECISION mui,pii,dpdti,dpdri,ui,dudti,dudri
    DOUBLE PRECISION pr,dpdtr,dpdrr,ur,dudtr,dudrr
    DOUBLE PRECISION dens2,dpdt,Aavg,Zavg

    if(type_eos.eq.1) then
      mui=0.d0
      do i=1,niso
        mui=mui+x(i)/a(i+1)
      enddo
      mue=1.d0/ye
      mui=1.d0/mui
      call eose(dens,temp,mue,pe,dpdte,dpdre,ue,dudte,dudre,error)
      do i=1,niso
        anum(i)=a(i+1)
        znum(i)=z(i+1)
      enddo
      call eosi(dens,temp,mui,ye,pii,dpdti,dpdri,ui,dudti,dudri,x,anum,znum,niso)
      call eosr(dens,temp,pr,dpdtr,dpdrr,ur,dudtr,dudrr)

      dens2=dens*dens
      dpdt=dpdte+dpdti+dpdtr
      v1=dpdt/dens2
      v2=temp*v1
      cv=dudte+dudti+dudtr
      pres=pe+pii+pr
      Utotal=ui+ur+ue
      dUdYe=0.d0

    else if(type_eos.eq.2) THEN

      abar=0.d0
      zbar=0.d0
      Aavg=0.d0
      Zavg=0.d0
      do i=1,niso
         abar=abar+x(i)/a(i+1)
         zbar=zbar+x(i)*z(i+1)/a(i+1)
         Aavg=Aavg+x(i)*a(i+1)
         Zavg=Zavg+x(i)*z(i+1)
      end do
      abar=1.d0/abar
      zbar=abar*zbar

      call helmeos(temp,dens,abar,zbar,pres,Utotal,dudt,dpdt,dUdYe)
      dens2=dens*dens
      v1=dpdt/dens2
      v2=temp*v1
      cv=dudt
    else
      stop 'Wrong type_eos'
    endif

  END SUBROUTINE eos

!**********************************************
