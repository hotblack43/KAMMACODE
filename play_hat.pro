FUNCTION hat,minutes
 if (n_elements(minutes) gt 1) then begin
; if the argument is an array
     value = 3.8047663250396013d-21*minutes^8-3.405102229245495d-11*minutes^4+2.4986575697934149d-6*minutes^2-3.0645854679789398d-4*minutes+2.6513956407895902d-3
     idx=where(minutes le 9 or minutes gt 234)
     if (idx(0) ne -1) then value(idx)=0.0
     endif else begin
; if the argument is a scalas
     value = 3.8047663250396013d-21*minutes^8-3.405102229245495d-11*minutes^4+2.4986575697934149d-6*minutes^2-3.0645854679789398d-4*minutes+2.6513956407895902d-3
     if (minutes lt 9) then value=0.0
     if (minutes gt 234) then value=0.0
     endelse
     value=abs(value)
; NOTE: this 'hat' function happens to have area 1
 return,value
 end



minutes=findgen(400)
plot,minutes,hat(minutes)
end
