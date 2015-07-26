basal=get_data('kamma_wrap/julian-dates/basal.txt')
bolus=get_data('kamma_wrap/julian-dates/bolus.txt')
openw,22,'kammas_basal.dat'
for ijd=min(long(basal(0,*))),max(long(basal(0,*)))-2,1 do begin
idx=where(long(basal(0,*)) eq ijd)
print,ijd,total(basal(1,idx))
printf,22,ijd,total(basal(1,idx))
endfor
close,22
openw,22,'kammas_bolus.dat'
for ijd=min(long(bolus(0,*))),max(long(bolus(0,*)))-2,1 do begin
idx=where(long(bolus(0,*)) eq ijd)
print,ijd,total(bolus(1,idx))
printf,22,ijd,total(bolus(1,idx))
endfor
close,22
basal=get_data('kammas_basal.dat')
bolus=get_data('kammas_bolus.dat')
caldat,basal(0,*),mm,dd,yy,hh,mi
fracyear=yy+(mm-1.)/12.+dd/365.25+hh/24./365.25
!P.MULTI=[0,1,2]
plot,xstyle=3,fracyear,bolus(1,*),xtitle='Year',ytitle='basal and bolus',title='Basal is red'
oplot,fracyear,basal(1,*),color=fsc_color('red')
plot,xstyle=3,fracyear,bolus(1,*)+basal(1,*),xtitle='Year',ytitle='basal + bolus'
x=fracyear
y=bolus(1,*)+basal(1,*)
res=robust_linefit(x,y)
yhat=res(0)+res(1)*x
oplot,x,yhat,linestyle=2
end
