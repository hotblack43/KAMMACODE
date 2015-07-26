c=get_data('kamma_wrap/julian-dates/carbs.txt')
openw,22,'kammas_mad.dat'
for ijd=min(long(c(0,*))),max(long(c(0,*)))-2,1 do begin
idx=where(long(c(0,*)) eq ijd)
print,ijd,total(c(1,idx))
printf,22,ijd,total(c(1,idx))
endfor
close,22
m=get_data('kammas_mad.dat')
caldat,m(0,*),mm,dd,yy,hh,mi
fracyear=yy+(mm-1.)/12.+dd/365.25+hh/24./365.25
;plot,m(0,*),m(1,*),xtitle='JD',ytitle='Daily total carbs [g]'
plot,xstyle=3,fracyear,m(1,*),xtitle='Year',ytitle='Daily total carbs [g]'
end
