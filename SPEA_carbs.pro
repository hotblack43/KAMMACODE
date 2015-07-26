!P.MULTI=[0,2,2]
!P.charsize=1.8
path='./kamma_wrap/julian-dates/'
data=get_data(path+'ig.txt')
ig_t=reform(data(0,*))
ig_y=reform(data(1,*))
print,'There are ',n_elements(ig_y),' ig values'
data=get_data(path+'carbs.txt')
c_t=reform(data(0,*))
c_y=reform(data(1,*))
caldat,c_t,mm,dd,yy,hh,mi,se
idx=where(c_y ge 16 and (hh gt 6 and hh lt 11))	; large morning meals only!
print,'There are ',n_elements(idx),' large morning meals.'
data=data(*,idx)
c_t=reform(data(0,*))
c_y=reform(data(1,*))
nc=n_elements(c_y)
print,'nc=',nc
w=0.25
openw,45,'large_meals.dat'
for ic=0,nc-1,1 do begin
ig_0=interpol(ig_y,ig_t,c_t(ic))
	kdx=where(24.*(ig_t-c_t(ic)) gt 0 and 24.*(ig_t-c_t(ic)) lt 4)
if (ic eq 0) then begin
	if (kdx(0) ne -1) then begin
		print,c_t(ic),'at carbs event # ',ic,' there are ',n_elements(kdx),' in the zone'
        	for klm=0,n_elements(kdx)-1,1 do printf,45,24.*(ig_t(kdx(klm))-c_t(ic)),ig_y(kdx(klm)),ig_y(kdx(klm))-ig_0
	endif
	if (kdx(0) eq -1) then print,c_t(ic),'at carbs event # ',ic,' there are ',0,' in the zone'
	plot,ytitle='ig',psym=3,24.*(ig_t-c_t(ic)),xtitle='Time [hours]',ig_y,xrange=[0,4],yrange=[0,20],xstyle=3,title='Large meals only (>16 g)'
endif 
if (ic ne 0) then begin
	if (kdx(0) ne -1) then begin
		print,c_t(ic),'at carbs event # ',ic,' there are ',n_elements(kdx),' in the zone'
        	for klm=0,n_elements(kdx)-1,1 do printf,45,24.*(ig_t(kdx(klm))-c_t(ic)),ig_y(kdx(klm)),ig_y(kdx(klm))-ig_0
	endif
	if (kdx(0) eq -1) then print,c_t(ic),'at carbs event # ',ic,' there are ',0,' in the zone'
	oplot,psym=3,24.*(ig_t-c_t(ic)),ig_y
endif
;a=get_kbrd()
endfor
for ic=0,nc-1,1 do begin
ig_0=interpol(ig_y,ig_t,c_t(ic))
if (ic eq 0) then begin
	plot,ytitle='!7D!3 ig',psym=3,24.*(ig_t-c_t(ic)),xtitle='Time [hours]',ig_y-ig_0,xrange=[0,4],yrange=[-20,20],xstyle=3,title='Large meals only (>16 g)'
endif 
if (ic ne 0) then begin
	oplot,psym=3,24.*(ig_t-c_t(ic)),ig_y-ig_0
endif
endfor
;-------------------------------------
data=get_data(path+'ig.txt')
ig_t=reform(data(0,*))
ig_y=reform(data(1,*))
data=get_data(path+'carbs.txt')
c_t=reform(data(0,*))
c_y=reform(data(1,*))
idx=where(c_y le 10)	; small meals only!
data=data(*,idx)
c_t=reform(data(0,*))
c_y=reform(data(1,*))
nc=n_elements(c_y)
for ic=0,nc-1,1 do begin
ig_0=interpol(ig_y,ig_t,c_t(ic))
if (ic eq 0) then begin
	plot,ytitle='ig',psym=3,24.*(ig_t-c_t(ic)),xtitle='Time [hours]',ig_y,xrange=[0,4],yrange=[0,20],xstyle=3,title='Small meals only (<10 g)'
endif 
if (ic ne 0) then begin
	oplot,psym=3,24.*(ig_t-c_t(ic)),ig_y
endif
endfor
for ic=0,nc-1,1 do begin
ig_0=interpol(ig_y,ig_t,c_t(ic))
if (ic eq 0) then begin
	plot,ytitle='!7D!3 ig',psym=3,24.*(ig_t-c_t(ic)),xtitle='Time [hours]',ig_y-ig_0,xrange=[0,4],yrange=[-20,20],xstyle=3,title='Small meals only (<10 g)'
endif 
if (ic ne 0) then begin
	oplot,psym=3,24.*(ig_t-c_t(ic)),ig_y-ig_0
endif
endfor
close,45
end

