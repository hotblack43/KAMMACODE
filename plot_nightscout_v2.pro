FUNCTION decorrlength,x,seriesname
 ac1=a_correlate(x,1)
 tau=(1.+ac1)/(1.-ac1)
 print,format='(a,f6.1,a)','Decorr time in series '+seriesname+' : ',tau,' steps.'
 return,tau
 end

PRO godofancyregress,bg,unfiltered,siglev
 ; applies formula for BG
 xx=[transpose(shift(bg,1)),transpose(shift(bg,2)),transpose(shift(bg,3)),transpose(shift(unfiltered,1)),transpose(shift(unfiltered,2)),transpose(shift(unfiltered,3))]
; givehere the names of the variables indicating shifts with X(-1), X(2) and so on
 varname=['bg(1)','bg(2)','bg(3)','UNF(1)','UNF(2)','UNF(3)']
 res=backw_elim(xx,bg,siglev,varlist=varlist,const=const,yfit=yfit,sigma=sigs)
 fmtstr='(a,1x,a7,1x,a,1x,g15.7,a,g12.7,a,f6.2)'
 print,'Backwards elimination applied. At S.L.= ',siglev
 print,'Regression intercept=',const
 for ivar=0,n_elements(varlist)-1,1 do begin
     print,format=fmtstr,' Coefficient for ',varname(varlist(ivar)),' is significant: ',res(ivar),' +/- ',sigs(ivar),', or Z= ',abs(res(ivar)/sigs(ivar))
     endfor
 tau=decorrlength(bg,'BG')
 R=correlate(bg,yfit)
 print,format='(a,f5.2)','R(BG,model)      : ',R
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
 
 ;==================================================================
 ; Version 2
 ; Bruger nightscout data og finder relationer mellem variabler
 ; kan finde relevante regressorer via 'backwards elimination'
 getscout,jd,noise,filtered,unfiltered,bg
 siglev=0.001	; Test significance level - the smaller the better
print,'-------------------------------------------'
 godofancyregress,bg,unfiltered,siglev
print,'-------------------------------------------'
 godofancyregress,bg,filtered,siglev
print,'-------------------------------------------'
 end
