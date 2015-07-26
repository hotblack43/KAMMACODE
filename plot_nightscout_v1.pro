PRO godofancyregress,bg,unfiltered
; applies formula for BG
xx=[transpose(shift(bg,-1)),transpose(shift(bg,-2)),transpose(shift(bg,-3)),transpose(shift(unfiltered,-1)),transpose(shift(unfiltered,-2)),transpose(shift(unfiltered,-3))]
res=regress(xx,bg,yfit=yhat,/double,sigma=sigs)
for i=0,n_elements(res)-1,1 do begin
print,i,res(i),' +/- ',sigs(i),' or, Z= ',abs(res(i)/sigs(i))
endfor
ac1=a_correlate(bg,1)
tau=(1.+ac1)/(1.-ac1)
print,'Decorr time in BG: ',tau,' steps.'
stop

return
end

PRO goregress,y,x
n=n_elements(x)
res=robust_linefit(x,y)
yhat=res(0)+res(1)*x
residuals=y-yhat
RMSE=sqrt(total(residuals^2)/float(n))
R=correlate(x,y)
print,' y = a + b*x'
print,' a =',res(0)
print,' b =',res(1)
print,'RMSE: ',RMSE
print,'R   : ',R
return
end

PRO goplot,x,y,noise,xstr,ystr,ifdiag
!P.MULTI=[0,2,2]
xra=[0,max([x])]
yra=[0,max([y])]
plot,title='Noise=1,2,3',psym=1,x,y,xtitle=xstr,ytitle=ystr,xstyle=3,ystyle=3,xrange=xra,yrange=yra
if (ifdiag eq 1) then oplot,[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],linestyle=2,color=fsc_color('red')
idx=where(noise eq 1)
plot,psym=1,x(idx),y(idx),title='Noise = 1',xtitle=xstr,ytitle=ystr,xstyle=3,ystyle=3,xrange=xra,yrange=yra
if (ifdiag eq 1) then oplot,[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],linestyle=2,color=fsc_color('red')
idx=where(noise eq 2)
plot,psym=1,x(idx),y(idx),title='Noise = 2',xtitle=xstr,ytitle=ystr,xstyle=3,ystyle=3,xrange=xra,yrange=yra
if (ifdiag eq 1) then oplot,[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],linestyle=2,color=fsc_color('red')
idx=where(noise eq 3)
plot,psym=1,x(idx),y(idx),title='Noise = 3',xtitle=xstr,ytitle=ystr,xstyle=3,ystyle=3,xrange=xra,yrange=yra
if (ifdiag eq 1) then oplot,[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],[max([!x.crange(0),!Y.crange(0)]),min([!x.crange(1),!Y.crange(1)])],linestyle=2,color=fsc_color('red')
return
end

PRO getscout,jd,noise,filtered,unfiltered,bg
data=get_data('nightscout_sgv.txt')
jd=reform(data(0,*))
noise=reform(data(2,*))
filtered=reform(data(4,*))
unfiltered=reform(data(5,*))
bg=reform(data(6,*))
return
end

;==============================================================
; Version 1
; PLotter nightscout data og finder relationer mellem variabler
getscout,jd,noise,filtered,unfiltered,bg
goplot,filtered,unfiltered,noise,'Filtered','Unfiltered',1
goplot,bg,unfiltered,noise,'BG','Unfiltered',0
goplot,bg,filtered,noise,'BG','Filtered',0
goplot,bg,unfiltered,noise,'BG','Unfiltered',0
;
print,'Regression BG vs UNfiltered:'
goregress,bg,unfiltered
print,'Regression BG vs Filtered:'
goregress,bg,filtered
; go do fancy regress
godofancyregress,bg,unfiltered
end
