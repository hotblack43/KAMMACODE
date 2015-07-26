openw,33,'kammas_good_night_dist.dat'
openw,34,'kammas_bad_night_dist.dat'
data=get_data('kammas_minmax_night_latest.dat')
l=size(data,/dimensions)
for irow=0,l(1)-1,1 do begin
	if (data(6,irow) le median(data(6,*))) then printf,33,format='(f15.7,5(1x,f9.5))',data(0:5,irow)
	if (data(6,irow) gt median(data(6,*))) then printf,34,format='(f15.7,5(1x,f9.5))',data(0:5,irow)
endfor
close,33
close,34
data=get_data('kammas_good_night_dist.dat')
JD=reform(data(0,*))
col1=reform(data(1,*))
col2=reform(data(2,*))
col3=reform(data(3,*))
col4=reform(data(4,*))
col5=reform(data(5,*))
;
l=size(data,/dimensions)
nrows=l(1)
nvarsplus1=l(0)
nvars=nvarsplus1-1
openw,33,'data.dat'
!P.MULTI=[0,1,2]
for k=0,nrows-1,1 do begin
	printf,33,jd(k),total(data(1:nvars,k))
	print,jd(k),total(data(1:nvars,k))
	if (k eq 0) then plot,data(1:nvars,k),psym=7,xtitle='Hour since start',ytitle='ba+bo',ystyle=3,xstyle=3,yrange=[0,max(data(1:nvars,*))]
	if (k gt 0) then oplot,data(1:nvars,k),psym=7
endfor
close,33
; overplot the bad ones
data=get_data('kammas_bad_night_dist.dat')
JD=reform(data(0,*))
col1=reform(data(1,*))
col2=reform(data(2,*))
col3=reform(data(3,*))
col4=reform(data(4,*))
col5=reform(data(5,*))
;
l=size(data,/dimensions)
nrows=l(1)
nvarsplus1=l(0)
nvars=nvarsplus1-1
for k=0,nrows-1,1 do begin
	print,jd(k),total(data(1:nvars,k)),' bad'
	oplot,data(1:nvars,k),psym=6
endfor
;-------------
data=get_data('data.dat')
histo,data(1,*),0.4,1.2,0.07,xtitle='Sum of ba+bo for entire good nights',/abs
end
