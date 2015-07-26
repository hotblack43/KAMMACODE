data=get_data('kammas_folded_data.dat')
c=get_data('kamma_wrap/julian-dates/carbs.txt')
ba=get_data('kamma_wrap/julian-dates/basal.txt')
bo=get_data('kamma_wrap/julian-dates/bolus.txt')
fivemin=get_data('kamma_wrap/julian-dates/kamma_5_minute.dat')
openw,22,'jwptest_ca.dat'
openw,23,'jwptest_in.dat'
print,'Going from ',min(long(data(0,*)))+1,' to ',max(long(data(0,*)))-1
for ijd=min(long(data(0,*)))+1,max(long(data(0,*)))-1,1 do begin
  idxca=where(long(c(0,*)) eq ijd)
  idxba=where(long(ba(0,*)) eq ijd)
  idxbo=where(long(bo(0,*)) eq ijd)
  idxfm=where(long(fivemin(0,*)) eq ijd)
  jdx=where(long(data(0,*)) eq ijd)
  jd=data(0,jdx)
  foldedcarbs=data(2,jdx)
  foldedboba=data(3,jdx)
  int_foldedcarbs=int_tabulated(jd,foldedcarbs,/double)
  int_foldedboba=int_tabulated(jd,foldedboba,/double)
  printf,22,format='(1(1x,i7),3(1x,f12.7))',ijd,288.*int_foldedcarbs,total(c(1,idxca)),total(fivemin(5,idxfm))
  printf,23,format='(1(1x,i7),3(1x,f12.7))',ijd,288.*int_foldedboba,total(bo(1,idxbo))+total(ba(1,idxba)),total(fivemin(3,idxfm))+total(fivemin(4,idxfm))
endfor
close,22
close,23

end
