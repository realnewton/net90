c     ************************************************************
c    --------------------------------------------------------------
c                          SUBROUTINE  eosr
c    --------------------------------------------------------------
c     ************************************************************

      subroutine eosr(rho,t,pr,dpdtr,dpdrr,ur,dudtr,dudrr) 
      implicit none
      DOUBLE PRECISION,intent(in)::rho,t
      DOUBLE PRECISION,intent(out)::pr,dpdrr,ur
      DOUBLE PRECISION t3,t4,dudtr,dudrr,dpdtr,entror      
      t3=t*t*t
      t4=t3*t
      pr=2.521d-15*t4
      ur=7.564d-15*t4/rho
      dpdtr=1.009d-14*t3
      dpdrr=0.d0
      dudtr=3.026d-14*t3/rho
      dudrr=-7.56d-15*t4/rho/rho
      entror=4.d0/3.d0*7.564d-15*t*t*t/rho
      return
      end
