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
;mbg=go_get_data('ig_alt.txt')
 mbg=get_data(path+'ig.txt')
 l=size(mbg,/dimensions)
 n=l(1)
 bolus=get_data(path+'bolus.txt')
 delta_1=4./24.
 for i=0,n-1,1 do begin
     tnow=mbg(0,i)
     caldat,tnow,mm,dd,yyy,hh,mi,se
    if (hh ge 6 and hh lt 11) then begin
         idx=where(carbs(0,*)-tnow le 0.0 and abs(carbs(0,*)-tnow) lt delta_1)
         jdx=where(bolus(0,*)-tnow le 0.0 and abs(bolus(0,*)-tnow) lt delta_1)
         if (idx(0) ne -1 and jdx(0) ne -1) then begin
; straight average over delta_1 interval:
;            print,i,mbg(1,i),mean(carbs(1,idx)),mean(bolus(1,jdx)),' mean'
;            printf,44,mbg(1,i),mean(carbs(1,idx)),mean(bolus(1,jdx))
; bolus is weighted with the 'fn' function
		dt_bolus=-(bolus(0,jdx)-tnow)
             print,i,tnow-long(tnow),mean(carbs(1,idx)),total(bolus(1,jdx)*fn(24.*dt_bolus)),mbg(1,i)
;            printf,44,tnow-long(tnow),mean(carbs(1,idx)),total(bolus(1,jdx)*fn(24.*dt_bolus)),mbg(1,i)
             printf,44,tnow-long(tnow),total(carbs(1,idx)*fn(24.*dt_bolus)),total(bolus(1,jdx)*fn(24.*dt_bolus)),mbg(1,i)
             endif
        endif
     endfor
 close,44
 print,'Now split kamma.data into the training and testing files ...'
 data=get_data('kamma.data')
 y=reform(data(3,*))
 x=reform(data(0:2,*))
 res=regress(x,y,/double,yfit=yhat,sigma=sigs,mcorr=mcorr,const=const)
 print,const
 print,res(0),' +/- ',sigs(0)
 print,res(1),' +/- ',sigs(1)
 print,res(2),' +/- ',sigs(2)
 plot,xstyle=3,ystyle=3,y,yhat,psym=7,xtitle='Observed BG',ytitle='Model'
 print,'R(x,Y) = ',mcorr
; now make a 'delta' version of kamma.data
; instead of ig use ig(t)-ig(t-delta) where delta=4 hrs
data=get_data('kamma.data')
t=reform(data(0,*))
c=reform(data(1,*))
b=reform(data(2,*))
ig=reform(data(3,*))
n=n_elements(t)
openw,55,'delta_kamma.data'
for i=0,n-1,1 do begin
ig_previous=interpol(ig,t,t(i)-4./24.)
printf,55,t(i),c(i),b(i),ig(i)-ig_previous
endfor
close,55
 end
