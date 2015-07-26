PRO gofindpointers,nrows,idx,jdx
 ; will produce two lists of mutually exclusive integers
 idx=[]
 jdx=[]
 for k=0L,nrows-1,1 do begin
     rnd=randomn(seed)
     if (rnd gt 0.0) then idx=[idx,k]
     if (rnd le 0.0) then jdx=[jdx,k]
     endfor
 end
 
 PRO validation,x_in,y_in
 x=x_in
 y=y_in
 ; TEST 1
 print,'----------------------------------------------------'
 print,'TEST 1: R from completely random numbers vs y'
 ; will validate the regression of y on x by performing same regression on entirely random numbers
 l=size(x,/dimensions)
 nvars=l(0)
 nrows=l(1)
 nMC=1000
 R=[]
 for iMC=0,nMC-1,1 do begin
     xx=randomn(seed,l)
     res=regress(xx,y,/double,yfit=yhat)
     R=[R,correlate(y,yhat)]
     endfor
 print,format='(a,f6.2,a,f6.2)','R on random numbers is: ',mean(R),' +/- ',stddev(r)
 ; TEST 2
 print,'----------------------------------------------------'
 print,'TEST 2: R validated on withheld subsets of input data'
 nMC=1000
 R1=[]
 R2=[]
 for iMC=0,nMC-1,1 do begin
     gofindpointers,nrows,idx,jdx
     xx=x_in(*,idx)
     yy=y_in(idx)
     res=regress(xx,yy,/double,yfit=yhat,const=const)
     R1=[R1,correlate(yy,yhat)]
     xx=x_in(*,jdx)
     yy=y_in(jdx)
     yhat2=const
     for ivar=0,nvars-1,1 do yhat2=yhat2+res(ivar)*xx(ivar,*)
     R2=[R2,correlate(yy,yhat2)]
     endfor
 print,format='(a,f6.2,a,f6.2)','a) R using non-withheld subset is: ',mean(R1),' +/- ',stddev(R1)
 print,format='(a,f6.2,a,f6.2)','b) R using     withheld subset is: ',mean(R2),' +/- ',stddev(R2)
 print,'Note: in calculating b) regression coefficients from a) were used.'
 print,'----------------------------------------------------'
 return
 end
 
 PRO goplot,jd,ig
 plot,24.*(jd-jd(0)),ig,xtitle='JD',ytitle='IG',psym=-7,xrange=[0,4]
 ;a=get_kbrd()
 return
 end
 
 PRO testforequality,array
 oldJD=reform(array(0,*))
 newJD=oldJD(sort(oldJD))
 if (product(oldJD eq newJD) ne 1) then stop
 return
 end
 
 PRO getthedatawewanttoselect,array,vnams
 common poi,numpointers
 ; read  in original data
 files='kammas_folded_data.dat'
 array=get_data(files)
 ; check if JD is in time-order
 testforequality,array
 ; get rid of the JD
 l=size(array,/dimensions)
 nvars=l(0)
 vnams=['JD','ig','ca','boplusba','f1','f2','f3','ig0','frd','tsc','cat','tsb','igb','bob','istw']
 ;-------------------------------------------------------------
 ; select morning meals, healthy/ill days
 ; NOTE - columns refer to whole list of vnams above
 idx=where(array(8,*) gt 0.5)
 ;idx=where(array(8,*) gt 0.75 and array(8,*) lt 0.9583)
 ;idx=where(array(0,*) gt julday(11,1,2014) and array(8,*) gt 0.75 and array(8,*) lt 0.9583)
 array=array(*,idx)
 ;-------------------------------------------------------------
 ; Now indicate wanted columns - remember to use 'old column numbers' and include y=ig !!!
 pointers=[0,1,2,3,7,8,9,10,13]     ; these are regressors, still using original numbering
 numpointers=n_elements(pointers)
 ; now reformat into final set of vars and their names
 array=array([pointers],*)
 vnams=vnams(pointers)
 return
 end
 
 function poi,name,vnams
 idx=where(vnams eq name)
 if (idx(0) eq -1) then stop
 return,idx(0)
 end
 
 
 ;===========================================================================================
 ; Code to extract 4-hour bits of data, integrate it, and perform a simple regression of 
 ; ig against integrated carbs and insulin
 ;===========================================================================================
 getthedatawewanttoselect,data,vnams
 ifsmooig=1
 nsmoo=5
 if (ifsmooig eq 1) then data(poi('ig',vnams),*)=smooth(data(poi('ig',vnams),*),nsmoo,/edge_truncate)
 idx=where(data(7,*) eq 100 and abs(data(6,*)) lt 5.1/60./24.)
 ; fix the above JWP!!!
 ;
 nidx=n_elements(idx)
 openw,44,'extracted_4hourdeltas.dat'
 ic=0
 for i=0,nidx-1,1 do begin
     jdx=where(data(poi('JD',vnams),*) ge data(poi('JD',vnams),idx(i)) and data(poi('JD',vnams),*) lt data(poi('JD',vnams),idx(i))+4./24.)
     jd=data(poi('JD',vnams),jdx)
     ig=data(poi('ig',vnams),jdx)
     foldedcarbs=data(poi('ca',vnams),jdx)
     foldedboba=data(poi('boplusba',vnams),jdx)
     frd=data(poi('frd',vnams),jdx)
     bob=data(poi('bob',vnams),jdx)
     nrows=n_elements(jdx)
     ; now get delta-ig and integrals over folded columns
     ;deltaig=ig(nrows-1)
     deltaig=ig(nrows-1)-ig(0)
     ;deltaig=max(ig)-min(ig)
     delta_t=24.*(max(jd)-min(jd));jd(nrows-1)-jd(0))
     ;delta_t=24.*(jd(nrows-1)-jd(0))
     int_foldedcarbs=int_tabulated(jd,foldedcarbs,/double)
     int_foldedboba=int_tabulated(jd,foldedboba,/double)
     if (delta_t gt 3.6) then begin
         goplot,jd,ig
         printf,44,format='(i3,6(1x,f12.8))',nrows,ig(0),deltaig,int_foldedcarbs,int_foldedboba,frd(0),bob(0)
         print,format='(i3,6(1x,f12.8))',nrows,ig(0),deltaig,int_foldedcarbs,int_foldedboba,frd(0),bob(0)
	ic=ic+1
         endif
     endfor
 close,44
 print,ic,' segments were long enoughand therefore selected'
 plotdat=get_data('extracted_4hourdeltas.dat')
 idx=where(plotdat(1,*) gt 4.0 and plotdat(1,*) lt 8.0)
 plotdat=plotdat(*,idx)
 y=reform(plotdat(2,*))
 x=reform(plotdat([3,4,6],*))
 ;x=reform(plotdat(3:6,*))
 res=regress(x,y,/double,yfit=yhat,const=const,sigma=sigs)
 print,'Mean delta-ig: ',mean(y)
 print,'const=',const
 fmt='(f9.4,a,f9.4,a,f5.1)'
 print,format=fmt,res(0),' +/- ',sigs(0),' Z= ',abs(res(0)/sigs(0))
 print,format=fmt,res(1),' +/- ',sigs(1),' Z= ',abs(res(1)/sigs(1))
 print,format=fmt,res(2),' +/- ',sigs(2),' Z= ',abs(res(2)/sigs(2))
 ;print,format=fmt,res(3),' +/- ',sigs(3),' Z= ',abs(res(3)/sigs(3))
 R=correlate(y,yhat)
 print,format='(a,f6.2,a,f6.2)','R on COMPLETE dataset is = ',R
 validation,x,y
 !P.LINESTYLE=0
 plot,y,color=fsc_color('red'),xtitle='#',ytitle='delta_ig and model (dashed)',xstyle=3,charsize=1.7,title='4-hour segments and accumulated quantities'
 !P.LINESTYLE=2
 oplot,yhat
 end
