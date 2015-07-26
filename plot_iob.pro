FUNCTION hat,minutes
if (n_elements(minutes) gt 1) then begin
     value = 3.8047663250396013E-21*minutes^8-3.405102229245495E-11*minutes^4+2.4986575697934149E-6*minutes^2-3.0645854679789398E-4*minutes+2.6513956407895902E-3
idx=where(minutes lt 9 or minutes gt 234) 
if (idx(0) ne -1) then value(idx)=0.0
endif else begin
     value = 3.8047663250396013E-21*minutes^8-3.405102229245495E-11*minutes^4+2.4986575697934149E-6*minutes^2-3.0645854679789398E-4*minutes+2.6513956407895902E-3
if (minutes lt 9) then value=0.0
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





t=findgen(300)
plot,t,-deriv(iob(t))
a=get_kbrd()
oplot,t,hat(t),color=fsc_color('red')
end
