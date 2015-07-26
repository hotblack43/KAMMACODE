FUNCTION hatcarbs,minutes
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
     value=value^3*10000.0d0
; NOTE: this 'hat' function happens to have area 1
 return,value
 end
 
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
 
 FUNCTION go_get_data,filename
 print,filename
 str="sed '1,$s/\// /g' "+filename+" | sed '1,$s/:/ /g' | sed '1,$s/\,/ /g' > out"
 spawn,str
 data=get_data('out')
 dd=data(0,*)
 mm=data(1,*)
 yy=data(2,*)
 hh=data(3,*)
 mi=data(4,*)
 val=data(5,*)
 jd=julday(mm,dd,yy,hh,mi)
 return,[jd,val]
 end
 
 FUNCTION logistic,t
 return,1.0d0/(1.0d0+exp(-t))
 end
 
 FUNCTION fn,hours
 days=hours/24.
 value=(-15.26+25.2*logistic(0.0983/(0.0033+days)))/10.0d0	; function is 1 at start and then decays
 return,value
 end
 
 
 openw,44,'kamma.data'
 path='./kamma_wrap/julian-dates/'
 carbs=get_data(path+'carbs.txt')
 ;ig=go_get_data('ig_alt.txt')
 ig=get_data(path+'ig.txt')
 l=size(ig,/dimensions)
 n=l(1)
 bolus=get_data(path+'bolus.txt')
 basal=get_data(path+'basal.txt')
 delta_1=4./24.
 for i=0,n-1,1 do begin
     tnow=ig(0,i)
     print,'tnow: ',tnow
     caldat,tnow,mm,dd,yyy,hh,mi,se
     if (hh ge 6 and hh lt 11) then begin
         idx=where(carbs(0,*)-tnow le 0.0 and abs(carbs(0,*)-tnow) lt delta_1)
         jdx=where(bolus(0,*)-tnow le 0.0 and abs(bolus(0,*)-tnow) lt delta_1)
         kdx=where(basal(0,*)-tnow le 0.0 and abs(basal(0,*)-tnow) lt delta_1)
        if (idx(0) ne -1 or jdx(0) ne -1) then begin
             ; straight average over delta_1 interval:
             ;            print,i,ig(1,i),mean(carbs(1,idx)),mean(bolus(1,jdx)),' mean'
             ;            printf,44,ig(1,i),mean(carbs(1,idx)),mean(bolus(1,jdx))
             ; bolus is weighted with the 'hat' function
	BB=0.0
        BA=0.0
	CC=0.0
	if(idx(0) ne -1) then begin
             dt_carbs=reform(-(carbs(0,idx)-tnow))
	print,'dt_carbs: ',dt_carbs
	CC=total(carbs(1,idx)*hatcarbs(60.*24.*dt_carbs))
endif
	if(jdx(0) ne -1) then begin
             dt_bolus=reform(-(bolus(0,jdx)-tnow))
	print,'dt_bolus: ',dt_bolus
	BB=total(bolus(1,jdx)*hat(60.*24.*dt_bolus))
endif
	if(kdx(0) ne -1) then begin
             dt_basal=reform(-(basal(0,kdx)-tnow))
	print,'dt_basal: ',dt_basal
	BA=total(basal(1,kdx)*hat(60.*24.*dt_basal))
endif
;             printf,44,format='(4(f9.5,1x),f15.7,1x,f8.5)',tnow-long(tnow),CC,BB,ig(1,i),tnow,BA
             printf,44,format='(4(f9.5,1x),f15.7,1x,f8.5)',tnow-long(tnow),total(carbs(1,idx)*hat(60.*24.*dt_carbs)),total(bolus(1,jdx)*hat(60.*24.*dt_bolus)),ig(1,i),tnow,BA
            endif
         endif
     endfor
 close,44
 data=get_data('kamma.data')
 y=reform(data(3,*))
 x=reform(data(1:2,*))
 jd=reform(data(4,*))
 basal=reform(data(5,*))
 x=[shift(x(0,*),4),shift(x(1,*),3),transpose(basal)]
 res=regress(x,y,/double,yfit=yhat,sigma=sigs,mcorr=mcorr,const=const)
 print,'Intercept:   ',const
 print,'carbs coeff: ',res(0),' +/- ',sigs(0)
 print,'bolus coeff: ',res(1),' +/- ',sigs(1)
 print,'basal coeff: ',res(2),' +/- ',sigs(2)
 plot,xstyle=3,ystyle=3,y,yhat,psym=7,xtitle='Observed BG',ytitle='Model'
 print,'R(x,Y) = ',mcorr
; write out kamma.data again
openw,78,'new_kamma.data'
 for klm=0,n_elements(y)-1,1 do begin
 printf,78,format='(4(f8.5,1x),f15.7,1x,f8.5)',data(0,klm),x(0,klm),x(1,klm),y(klm),jd(klm),basal(klm)
 endfor
 close,78
 spawn,'cp new_kamma.data kamma.data'

 
 ; now make a 'delta' version of kamma.data
 ; instead of ig use ig(t)-ig(t-delta) where delta=4 hrs
 data=get_data('kamma.data')
 t=reform(data(0,*))
 c=reform(data(1,*))
 b=reform(data(2,*))
 ig=reform(data(3,*))
 jd=reform(data(4,*))
 basal=reform(data(5,*))
 n=n_elements(t)
 openw,55,'delta_kamma.data'
 for i=0,n-1,1 do begin
;     ig_previous=interpol(ig,t,t(i)-4./24.)
     ig_previous=interpol(ig,t,0.5)
     if (c(i) ne 0 and b(i) ne 0) then printf,55,format='(4(f9.5,1x),f15.7,1x,f8.5)',t(i),c(i),b(i),ig(i)-ig_previous,jd(i),basal(i)
     endfor
 close,55
; plot
data=get_data('delta_kamma.data')
!P.MULTI=[0,1,2]
plot,data(1,*),data(3,*),psym=7,xtitle='Carbs',ytitle='Delta ig',xstyle=3,ystyle=3
res=robust_linefit(data(1,*),data(3,*),yfit)                                      
oplot,data(1,*),yfit,color=fsc_color('red') 
;
plot,data(2,*),data(3,*),psym=7,xtitle='Bolus',ytitle='Delta ig',xstyle=3,ystyle=3
res=robust_linefit(data(2,*),data(3,*),yfit)                                      
oplot,data(2,*),yfit,color=fsc_color('red') 
 end
