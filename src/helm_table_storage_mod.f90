MODULE helm_table_storage_mod

      implicit none
! sizes of the tables
! normal table, big table, bigger table, denser bigger table

      integer          imax,jmax

! original
!      parameter        (imax = 211, jmax = 71)

! standard
!      parameter        (imax = 271, jmax = 101)

! twice as dense
      parameter        (imax = 541, jmax = 201)

! half as dense
!      parameter        (imax = 136, jmax = 51)



! for the electrons
! density and temperature
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi
      double precision d(imax),t(jmax)
      common /dttabc1/ d,t, &
                       tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi

! for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax), &
                       ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax), &
                       fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax), &
                       fddtt(imax,jmax)

      common /frtabc1/ f,fd, &
                       ft,fdd,ftt, &
                       fdt,fddt,fdtt, &
                       fddtt

! for the pressure derivative with density tables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax), &
                       dpdft(imax,jmax),dpdfdt(imax,jmax)

      common /dpdtab1/ dpdf,dpdfd, &
                       dpdft,dpdfdt


! for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax), &
                       eft(imax,jmax),efdt(imax,jmax)

      common /eftabc1/ ef,efd, &
                       eft,efdt


! for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax), &
                       xft(imax,jmax),xfdt(imax,jmax)

      common /xftabc1/ xf,xfd, &
                       xft,xfdt


! for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax), &
                       dti_sav(jmax),dt2i_sav(jmax),dt3i_sav(jmax), &
                       dd_sav(imax),dd2_sav(imax), &
                       ddi_sav(imax),dd2i_sav(imax),dd3i_sav(imax)

      common /diftabc1/dt_sav,dt2_sav, &
                       dti_sav,dt2i_sav,dt3i_sav, &
                       dd_sav,dd2_sav, &
                       ddi_sav,dd2i_sav,dd3i_sav




! for the ions
! density and temperature
      double precision tion_lo,tion_hi,tion_stp,tion_stpi, &
                       dion_lo,dion_hi,dion_stp,dion_stpi
      double precision dion(imax),tion(jmax)
      common /dttabc2/ dion,tion, &
                       tion_lo,tion_hi,tion_stp,tion_stpi, &
                       dion_lo,dion_hi,dion_stp,dion_stpi

! for the helmholtz free energy tables
      double precision fion(imax,jmax),fiond(imax,jmax), &
                       fiont(imax,jmax),fiondd(imax,jmax), &
                       fiontt(imax,jmax),fiondt(imax,jmax), &
                       fionddt(imax,jmax),fiondtt(imax,jmax), &
                       fionddtt(imax,jmax)

      common /frtabc2/ fion,fiond, &
                       fiont,fiondd,fiontt, &
                       fiondt,fionddt,fiondtt, &
                       fionddtt

! for the pressure derivative with density tables
      double precision dpiondf(imax,jmax),dpiondfd(imax,jmax), &
                       dpiondft(imax,jmax),dpiondfdt(imax,jmax)

      common /dpdtab2/ dpiondf,dpiondfd, &
                       dpiondft,dpiondfdt


! for chemical potential tables
      double precision efion(imax,jmax),efiond(imax,jmax), &
                       efiont(imax,jmax),efiondt(imax,jmax)

      common /eftabc2/ efion,efiond, &
                       efiont,efiondt


! for the number density tables
      double precision xfion(imax,jmax),xfiond(imax,jmax), &
                       xfiont(imax,jmax),xfiondt(imax,jmax)

      common /xftabc2/ xfion,xfiond, &
                       xfiont,xfiondt


! for storing the differences
      double precision dt_sav_ion(jmax),dt2_sav_ion(jmax), &
                       dti_sav_ion(jmax),dt2i_sav_ion(jmax), &
                       dt3i_sav_ion(jmax),dd_sav_ion(imax), &
                       dd2_sav_ion(imax),ddi_sav_ion(imax), &
                       dd2i_sav_ion(imax),dd3i_sav_ion(imax)

      common /diftabc2/dt_sav_ion,dt2_sav_ion, &
                       dti_sav_ion,dt2i_sav_ion, &
                       dt3i_sav_ion,dd_sav_ion, &
                       dd2_sav_ion,ddi_sav_ion, &
                       dd2i_sav_ion,dd3i_sav_ion

end module
