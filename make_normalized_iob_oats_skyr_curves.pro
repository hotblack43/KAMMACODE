FUNCTION hat,minutes
if (n_elements(minutes) gt 1) then begin
     value = 3.8047663250396013E-21*minutes^8-3.405102229245495E-11*minutes^4+2.4986575697934149E-6*minutes^2-3.0645854679789398E-4*minutes+2.6513956407895902E-3
idx=where(minutes le 9 or minutes gt 234) 
if (idx(0) ne -1) then value(idx)=0.0
endif else begin
     value = 3.8047663250396013E-21*minutes^8-3.405102229245495E-11*minutes^4+2.4986575697934149E-6*minutes^2-3.0645854679789398E-4*minutes+2.6513956407895902E-3
if (minutes le 9) then value=0.0
if (minutes gt 234) then value=0.0
endelse
return,-value
end

FUNCTION iob,t
; t is time in minutesUTES
value = 0.992473793241444 + 0.00265139564078959*t + 8.32885856597805e-7*t^3 + 4.22751813893289e-22*t^9 - 6.81020445849099e-12*t^5 - 0.000153229273398947*t^2
idx=where(t gt 234)
if (idx(0) ne -1) then begin
	value(idx)=0.0
endif
return,value
end

PRO savenormalized,x,y,title
n=n_elements(x)
integral=int_tabulated(findgen(n),y,/double)
;integral=int_tabulated(x,y,/double)
y=y/integral
integral2=int_tabulated(findgen(n),y,/double)
;integral2=int_tabulated(x,y,/double)
print,'Integral is: ',integral2
openw,34,title
for i=0,n_elements(x)-1,1 do begin
printf,34,x(i),y(i)
endfor
close,34
return
end

; OATS
!P.charsize=2
!P.thick=4
data=get_data('oat_meals_breakfast.dat')
jd=reform(data(0,*))
dt=reform(data(1,*))/24.0d0	; now also in days
ig=reform(data(2,*))
dig=reform(data(3,*))
;
binwidth=5./60./24.0d0;0.1
openw,1,'oats_binned.dat'
for i=0.0,3.5/24.0d0-binwidth,binwidth do begin
;for i=0.0,4.-binwidth,binwidth do begin
idx=where(dt gt i and dt le i+binwidth)
printf,1,i+binwidth/2.,max([0.0,median(dig(idx))])
endfor
close,1
;
data_binned=get_data('oats_binned.dat')
savenormalized,data_binned(0,*),data_binned(1,*),'oats_normalized.dat'
plot,dt,dig,psym=7,xtitle='Days since meal start',ytitle='ig'
data_binned=get_data('oats_normalized.dat')
oplot,thick=4,data_binned(0,*),data_binned(1,*),color=fsc_color('red')
; SKYR
data=get_data('skyr_meals_breakfast.dat')
jd=reform(data(0,*))
dt=reform(data(1,*))/24.0d0	; now also in days
ig=reform(data(2,*))
dig=reform(data(3,*))
;
openw,1,'skyr_binned.dat'
for i=0.0,3.5/24.0-binwidth,binwidth do begin
;for i=0.0,4.-binwidth,binwidth do begin
idx=where(dt gt i and dt le i+binwidth)
printf,1,i+binwidth/2.,max([0.0,median(dig(idx))])
endfor
close,1
;
data_binned=get_data('skyr_binned.dat')
oplot,dt,dig,psym=1
savenormalized,data_binned(0,*),data_binned(1,*),'skyr_normalized.dat'
data_binned=get_data('skyr_normalized.dat')
oplot,data_binned(0,*),data_binned(1,*),color=fsc_color('green')

; and now iob
t=findgen(240)
y=hat(t)
openw,44,'iob_kernel.dat'
for i=0.,4.*60.,5. do begin
printf,44,i/60./24.,interpol(y,t,i*1.0)
endfor
close,44
data=get_data('iob_kernel.dat')
t=reform(data(0,*))
y=reform(data(1,*))
n=n_elements(y)
y=y/int_tabulated(findgen(n),y)
;y=y/int_tabulated(t,y)
openw,44,'iob_kernel.dat'
for i=0,n_elements(t)-1,1 do begin
printf,44,t(i),y(i)
endfor
close,44
data=get_data('iob_kernel.dat')
t=reform(data(0,*))
y=reform(data(1,*))
print,int_tabulated(findgen(n),y)
end
