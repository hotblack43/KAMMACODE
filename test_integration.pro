 PRO make_holes,infile
; will make holæes in raw_table.dat just like the ones found in "infile"
; NOTE: this does require that the two files have overlapping time axes ....
data=get_data(infile)
twant=reform(data(0,*))
n=n_elements(twant)
data=get_data('raw_table.dat')
t_have=reform(data(0,*))
y1=reform(data(1,*))
y2=reform(data(2,*))
; interpolate
y1_out=interpol(y1,t_have,twant)
y2_out=interpol(y2,t_have,twant)
; clip to avoid end-effects
idx=where(twant ge min(t_have) and twant le max(t_have))
n=n_elements(y1_out)
openw,29,'temp.dat'
for i=0L,n-1,1 do begin
printf,29,format='(f15.7,2(1x,f12.8))', twant(i),y1_out(i),y2_out(i)
endfor
close,29
; now overwrite files ...
spawn,'mv temp.dat raw_table.dat'
 return
 end

 PRO gettheKERNEL,filename,kernel
 print,' Kernel file: ',filename
 data=get_data(filename)
 length=n_elements(data(1,*))
 time=reform(data(0,0:length-7-1))
 print,'stepsize in kernel: ',(time(5)-time(4))*60.*24.,' expressed in minutes'
 print,'stepsize in kernel: ',(time(5)-time(4)),' in file units [days]'
;
 kernel=reform(data(1,0:length-7-1))
 length=n_elements(kernel)
 kernel=[indgen(length)*0,kernel]
;
 n=n_elements(kernel)
 return
 end
 
 PRO foldwithKERNEL,data,kernel,time,folded_obs
 ; define a observation function
 time=reform(data(0,*))
 print,'stepsize in data: ',(time(5)-time(4))*60.*24.,'expressed in  minutes'
 print,'stepsize in data: ',(time(5)-time(4)),'expressed in file units [days]'
 observation=reform(data(1,*))
 ; perform convolution
 folded_obs=convol(observation,reverse(kernel),/edge_truncate,/center)
 return
 end
PRO setup_raw_table,ndays
 openw,23,'raw_table.dat'
 for jd=2456892.0451389d0,2456892.0451389d0+ndays,5./60./24. do begin
     intday=long(jd)
     frd=jd-intday
     meal=randomu(seed)*25.0	; max 25 g per meal
     bolus=randomu(seed)*2.2	; max 2.2 units insulin per bolus shot
     flip=randomu(seed)
     if (frd gt 0.75 and frd lt 0.985 and flip lt 0.1)then begin
         printf,23,format='(f15.7,2(1x,f9.5))',jd,meal,bolus
         endif else begin
         printf,23,format='(f15.7,2(1x,f9.5))',jd,0.0,0.0
         endelse
     endfor
 close,23
 return
 end
 
 ;---------------------------------------------
 ; Code to test how integration, folding and 
 ; kernel normalization works together
 ;---------------------------------------------
 ; first set up a list of howlarge the meals were on the days - each day is just one number
 ndays=365
 want_holes=0
 setup_raw_table,ndays
 data=get_data('raw_table.dat')
 if (want_holes eq 1) then make_holes,'./kamma_wrap/julian-dates/kamma_5_minute.dat'
 jd=reform(data(0,*))
 meals=reform(data(1,*))
 bolus=reform(data(2,*))
print,'A total of ',total(meals,/double),' g carbs was eaten in the period.'
; get the kernels
 gettheKERNEL,'oats_normalized.dat',oats_kernel
 gettheKERNEL,'iob_kernel.dat',iob_kernel
print,'===================================================================='
; now fold with the oats kernel
 foldwithKERNEL,data([0,1],*),oats_kernel,time1,folded_carbs
print,'Total of folded carbs: ',total(folded_carbs,/double),' g for the period.'
print,'288*Integral over folded carbs: ',288.0d0*int_tabulated(time1,folded_carbs,/double),' g for the period.'
print,'===================================================================='
print,'A total of ',total(bolus,/double),' units insulin was given in the period.'
; now fold with the iob kernel
 foldwithKERNEL,data([0,2],*),iob_kernel,time2,folded_bolus
print,'Total of folded bolus: ',total(folded_bolus,/double),' units for the period.'
print,'288*Integral over folded bolus: ',288.0d0*int_tabulated(time2,folded_bolus,/double),' units for the period.'
print,'===================================================================='
 end
 
