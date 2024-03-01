  SUBROUTINE read_rates_table

    USE nuclear90_module, only:nt,nr,nc,data_rate,logtemp_rates,&
                              &logrhoYe_rates,rates_file

    implicit NONE
    INTEGER i,j,k

    open(18,file=rates_file,form='unformatted')
    do i=1,nt*nr
      read(18) logtemp_rates((i-1)/nr+1),logrhoYe_rates(mod(i-1,nr)+1),&
            &   (data_rate(i,j),j=1,nc)
    enddo
    close(18)
  END SUBROUTINE
