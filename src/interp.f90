subroutine interp(t,rhoye,interpolated)

  USE nuclear90_module,only:nt,nr,nc,logtemp_rates,&
                            &logrhoYe_rates,data_rate

  implicit none
  INTEGER i,it,ir,it1,it2,ir1,ir2,v1,v2,v3,v4
  DOUBLE PRECISION logtemp,logrhoye,x2x,y2y,xx1,yy1,x2x1,y2y1
  DOUBLE PRECISION, intent(in):: t,rhoye
  DOUBLE PRECISION, dimension(nc),intent(out)::interpolated

  logtemp=log10(T)
  logrhoye=log10(rhoye)

  i=1
  do while(logtemp.gt.logtemp_rates(1))
    if(logtemp.ge.logtemp_rates(i).and.logtemp.lt.logtemp_rates(i+1))exit
    i=i+1
    if(i.eq.nt) exit
  enddo
  if(i.eq.1) then !logtemp.le.logtemp_rates(1)
    it1=1
    it2=2
  else if(i.eq.nt) then !logtemp.ge.logtemp_rates(nt)
    it1=nt-1
    it2=nt
  else
    it1=i
    it2=i+1
  endif

  i=1
  do while(logrhoye.gt.logrhoye_rates(1))
    if(logrhoye.ge.logrhoYe_rates(i).and.logrhoye.lt.logrhoYe_rates(i+1))exit
    i=i+1
    if(i.eq.nr) exit
  enddo
  if(i.eq.1) then !logtemp.le.logtemp_rates(1)
    ir1=1
    ir2=2
  else if(i.eq.nr) then !logtemp.ge.logtemp_rates(nt)
    ir1=nr-1
    ir2=nr
  else
    ir1=i
    ir2=i+1
  endif

  v1=(it1-1)*nr+ir1
  v2=(it1-1)*nr+ir2
  v3=(it2-1)*nr+ir1
  v4=(it2-1)*nr+ir2
  x2x=logtemp_rates(it2)-logtemp
  xx1=logtemp-logtemp_rates(it1)
  y2y=logrhoYe_rates(ir2)-logrhoye
  yy1=logrhoye-logrhoYe_rates(ir1)
  x2x1=logtemp_rates(it2)-logtemp_rates(it1)
  y2y1=logrhoYe_rates(ir2)-logrhoYe_rates(ir1)
  !interpolate
  interpolated(:)=(data_rate(v1,:)*x2x*y2y+&
                  &data_rate(v3,:)*xx1*y2y+&
                  &data_rate(v2,:)*x2x*yy1+&
                  &data_rate(v4,:)*xx1*yy1)/&
                  &(x2x1*y2y1)

end subroutine
