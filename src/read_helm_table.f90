subroutine read_helm_table

use helm_table_storage_mod
include 'implno.dek'

! this routine reads the helmholtz eos file, and
! must be called once before the helmeos routine is invoked.

! declare local variables
integer          i,j
double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                 dd,dd2,ddi,dd2i,dd3i


! open the file (use softlinks to input the desired table)
 open(unit=19,file='src/helm_table.dat',status='old')

! for standard table limits
 tlo   = 3.0d0
 thi   = 13.0d0
 tstp  = (thi - tlo)/float(jmax-1)
 tstpi = 1.0d0/tstp
 dlo   = -12.0d0
 dhi   = 15.0d0
 dstp  = (dhi - dlo)/float(imax-1)
 dstpi = 1.0d0/dstp

! read the helmholtz free energy and its derivatives
 do j=1,jmax
  tsav = tlo + (j-1)*tstp
  t(j) = 10.0d0**(tsav)
  do i=1,imax
   dsav = dlo + (i-1)*dstp
   d(i) = 10.0d0**(dsav)
   read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
            fddt(i,j),fdtt(i,j),fddtt(i,j)
  enddo
 enddo
!       write(6,*) 'read main table'


! read the pressure derivative with density table
 do j=1,jmax
  do i=1,imax
   read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
  enddo
 enddo

! read the electron chemical potential table
 do j=1,jmax
  do i=1,imax
   read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
  enddo
 enddo

! read the number density table
 do j=1,jmax
  do i=1,imax
   read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
  enddo
 enddo

! close the file
close(unit=19)

! construct the temperature and density deltas and their inverses
 do j=1,jmax-1
  dth          = t(j+1) - t(j)
  dt2         = dth * dth
  dti         = 1.0d0/dth
  dt2i        = 1.0d0/dt2
  dt3i        = dt2i*dti
  dt_sav(j)   = dth
  dt2_sav(j)  = dt2
  dti_sav(j)  = dti
  dt2i_sav(j) = dt2i
  dt3i_sav(j) = dt3i
 end do
 do i=1,imax-1
  dd          = d(i+1) - d(i)
  dd2         = dd * dd
  ddi         = 1.0d0/dd
  dd2i        = 1.0d0/dd2
  dd3i        = dd2i*ddi
  dd_sav(i)   = dd
  dd2_sav(i)  = dd2
  ddi_sav(i)  = ddi
  dd2i_sav(i) = dd2i
  dd3i_sav(i) = dd3i
 enddo

return
end
