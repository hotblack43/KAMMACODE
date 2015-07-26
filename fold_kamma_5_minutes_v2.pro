 PRO gettheKERNEL,filename,kernel
 print,'File: ',filename
 data=get_data(filename)
 length=n_elements(data(1,*))
 time=reform(data(0,0:length-7-1))
 print,'stepsize in kernel: ',(time(5)-time(4))*60.*24.,' minutes'
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
 print,'stepsize in data: ',(time(5)-time(4))*60.*24.,' minutes'
 observation=reform(data(1,*))
 ; perform convolution
 folded_obs=convol(observation,reverse(kernel),/edge_truncate,/center)
 return
 end
 
;===============================================================
; Version 2 of code that folds carbs and insulin with kernels
; Fixes ig=-1 bug form Version 1
 idx=[30000]
 !P.MULTI=[0,1,3]
 !P.CHARSIZE=1.7
 alldata=get_data('./kamma_wrap/julian-dates/kamma_5_minute.dat')
 gettheKERNEL,'oats_normalized.dat',oats_kernel
 gettheKERNEL,'iob_kernel.dat',iob_kernel
 ; IG
 ig_data=alldata([0,1],*)
 ig=reform(ig_data(1,*))
; get rid of lines with ig=-1 by filling with interpolated values
 idx=where(ig ne -1)
 jdx=where(ig eq -1)
 alldata(1,jdx)=interpol(ig(idx),alldata(0,idx),alldata(0,jdx))
; and start over
 ig_data=alldata([0,1],*)
 ig=reform(ig_data(1,*))
;plot,xstyle=3,psym=-7,ig_data(0,idx),ig(idx),ytitle='IG'
 ; CARBS
 data=alldata([0,5],*)
 foldwithKERNEL,data,oats_kernel,time,folded_carbs
;plot,xstyle=3,psym=-7,time(idx),data(1,idx),ytitle='Carbs'
;oplot,time(idx),folded_carbs(idx),color=fsc_color('red')              
 ; BO
 data=alldata([0,4],*)
 foldwithKERNEL,data,iob_kernel,time,folded_bo
;plot,time(idx),data(1,idx),ytitle='BO (green), BA(blue)'              
;oplot,time(idx),folded_bo(idx),color=fsc_color('green')
 ; BA
 data=alldata([0,3],*)
 foldwithKERNEL,data,iob_kernel,time,folded_ba
;oplot,time(idx),data(1,idx)
;oplot,time(idx),folded_ba(idx),color=fsc_color('blue')
 ; FLAGS
 F1=reform(alldata(6,*))
 F2=reform(alldata(7,*))
 F3=reform(alldata(8,*))
 cag=reform(alldata(9,*))
 frd=reform(alldata(10,*))
 toc=reform(alldata(11,*))
 cat=reform(alldata(12,*))
 tsb=reform(alldata(13,*))
 igb=reform(alldata(14,*))
 bob=reform(alldata(15,*))
 isw=reform(alldata(16,*))
 fmt='(f15.7,3(1x,f9.5),3(1x,i2),3(1x,f15.7),1x,i4,4(1x,f9.5))'
; str='kammas_folded_data_'+string(hour,format='(f6.4)')+'.dat'
 str='kammas_folded_data.dat'
 openw,55,str
 for k=0L,n_elements(time)-1,1 do begin
         printf,55,format=fmt,time(k),ig(k),folded_carbs(k),folded_bo(k)+folded_ba(k),f1(k),f2(k),f3(k),cag(k),frd(k),toc(k),cat(k),tsb(k),igb(k),bob(k),isw(k)
 endfor
 close,55
 !P.MULTI=[0,1,1]
 data=get_data(str)
 t=reform(data(0,*))
 y=reform(data(1,*))
 x1=reform(data(2,*))
 x2=reform(data(3,*))
 jdx=where(x1 gt 0 and (f3 ne 1))
 res=regress([transpose(x1(jdx)),transpose(x2(jdx))],y(jdx),yfit=yhat,/double)
plot,y(jdx),yhat,psym=3,xtitle='ig',ytitle='Simple ig model',charsize=1.7
print,'R: ',correlate(y(jdx),yhat)
; test
 print,'test of file ',str
print,'Total and 288*Integral over ig: ',total(data(1,*),/double),288.*int_tabulated(data(0,*),data(1,*),/double)
print,'Total and 288*Integral over ca: ',total(data(2,*),/double),288.*int_tabulated(data(0,*),data(2,*),/double)
print,'Total and 288*Integral over bo+ba: ',total(data(3,*),/double),288.*int_tabulated(data(0,*),data(3,*),/double)
;endfor
 end
