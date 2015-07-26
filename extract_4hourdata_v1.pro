PRO goplot,jd,ig
plot,24.*(jd-jd(0)),ig,xtitle='JD',ytitle='IG',psym=-7,xrange=[0,4]
a=get_kbrd()
return
end

PRO testforequality,array
 oldJD=reform(array(0,*))
 newJD=oldJD(sort(oldJD))
 if (product(oldJD eq newJD) ne 1) then stop
 return
 end
 
PRO getthedatawewanttoselect,array,varnames
 common poi,numpointers
 ; read  in original data
 files='kammas_folded_data.dat'
 array=get_data(files)
 ; check if JD is in time-order
 testforequality,array
 ; get rid of the JD
 l=size(array,/dimensions)
 nvars=l(0)
 varnames=['JD','ig','ca','boplusba','f1','f2','f3','ig0','frd','tsc','cat','tsb','igb','bob','istw']
 ;-------------------------------------------------------------
 ; select morning meals, healthy/ill days
 ; NOTE - columns refer to whole list of varnames above
 idx=where(array(8,*) gt 0.75)
 ;idx=where(array(8,*) gt 0.75 and array(8,*) lt 0.9583)
 ;idx=where(array(0,*) gt julday(11,1,2014) and array(8,*) gt 0.75 and array(8,*) lt 0.9583)
 array=array(*,idx)
 ;-------------------------------------------------------------
 ; Now indicate wanted columns - remember to use 'old column numbers' and include y=ig !!!
 pointers=[0,1,2,3,7,8,9,10,13]     ; these are regressors, still using original numbering
 numpointers=n_elements(pointers)
 ; now reformat into final set of vars and their names
 array=array([pointers],*)
 varnames=varnames(pointers)
 return
 end
 
 function poi,name,varnames
 idx=where(varnames eq name)
 if (idx(0) eq -1) then stop
 return,idx(0)
 end
  
 
 ;===========================================================================================
 ; Code to extract 4-hour bits of data, integrate it, and perform a simple regression of 
 ; ig against integrated carbs and insulin
 ;===========================================================================================
 getthedatawewanttoselect,data,varnames
 ifsmooig=0
 nsmoo=3
 if (ifsmooig eq 1) then data(poi('ig',varnames),*)=smooth(data(poi('ig',varnames),*),nsmoo,/edge_truncate)
 idx=where(data(7,*) eq 100 and abs(data(6,*)) lt 5.1/60./24.)
 ; fix the above JWP!!!
 ;
 nidx=n_elements(idx)
 openw,44,'extracted_4hourdeltas.dat'
 for i=0,nidx-1,1 do begin
     jdx=where(data(poi('JD',varnames),*) ge data(poi('JD',varnames),idx(i)) and data(poi('JD',varnames),*) lt data(poi('JD',varnames),idx(i))+4./24.)
     jd=data(poi('JD',varnames),jdx)
     ig=data(poi('ig',varnames),jdx)
	goplot,jd,ig
     foldedcarbs=data(poi('ca',varnames),jdx)
     foldedboba=data(poi('boplusba',varnames),jdx)
     frd=data(poi('frd',varnames),jdx)
     bob=data(poi('bob',varnames),jdx)
     nrows=n_elements(jdx)
     ; now get delta-ig and integrals over folded columns
     ;deltaig=ig(nrows-1)
     delta_t=24.*(jd(nrows-1)-jd(0))
     deltaig=ig(nrows-1)-ig(0)
     int_foldedcarbs=int_tabulated(jd,foldedcarbs,/double)
     int_foldedboba=int_tabulated(jd,foldedboba,/double)
	if (delta_t gt 2.5) then begin
     printf,44,format='(i3,6(1x,f12.8))',nrows,ig(0),deltaig,int_foldedcarbs,int_foldedboba,frd(0),bob(0)
     print,format='(i3,6(1x,f12.8))',nrows,ig(0),deltaig,int_foldedcarbs,int_foldedboba,frd(0),bob(0)
	endif
     endfor
 close,44
 help,nidx
 plotdat=get_data('extracted_4hourdeltas.dat')
 idx=where(plotdat(1,*) gt 4.0 and plotdat(1,*) lt 8.0)
 plotdat=plotdat(*,idx)
 y=reform(plotdat(2,*))
 x=reform(plotdat(3:6,*))
 res=regress(x,y,/double,yfit=yhat,const=const,sigma=sigs)
 print,'Mean delta-ig: ',mean(y)
 print,'const=',const
 print,res(0),' +/- ',sigs(0)
 print,res(1),' +/- ',sigs(1)
 print,res(2),' +/- ',sigs(2)
 print,res(3),' +/- ',sigs(3)
 R=correlate(y,yhat)
	print,'R= ',R
 	plot,y,color=fsc_color('red'),xtitle='#',ytitle='delta_ig and model'
	!P.LINESTYLE=2
	oplot,yhat
 end
