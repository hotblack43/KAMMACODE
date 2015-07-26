 PRO getthedatawewanttoselect,array,varnames
 common poi,numpointers
; read  in original data
 files='kammas_folded_data.dat'
 array=get_data(files)
; get rid of the JD
 l=size(array,/dimensions)
 nvars=l(0)
 varnames=['JD','ig','ca','boplusba','f1','f2','f3','ig0','frd','tsc','cat','tsb','igb','bob']
;-------------------------------------------------------------
 ; select morning meals, healthy/ill days
; NOTE - columns refer to whole list of varnames above
 idx=where(array(0,*) gt julday(11,1,2014) and array(7,*) ne -1 and array(6,*) eq 1 and array(8,*) gt 0.75 and array(8,*) lt 0.9583)
 array=array(*,idx)
;-------------------------------------------------------------
; Now indicate wanted columns - remember to use 'old column numbers' and include y=ig !!!
 pointers=[1,2,3,7,8,9,10,11,12,13]	; these are regressors, still using original numbering
 numpointers=n_elements(pointers)
; now reformat into final set of vars and their names
 array=array([pointers],*)	
 varnames=varnames(pointers)
 return
 end

PRO doablock,xx,yy_in,varnames,pointer1,pointer2,siglev
 common remember,Rmax
 yy=yy_in
 xx=[[xx],(xx(pointer1,*)*xx(pointer2,*))]
 varnames=[varnames,'('+varnames(pointer1+1)+')*('+varnames(pointer2+1)+')']
 dummy=backw_elim(xx,yy,siglev,varlist=varlist,const=const,yfit=yfit,sigma=sigs)
 residuals=yy-yfit
 sd=stddev(residuals)
 hdx=where(abs(residuals) gt 2.*sd)
 R=correlate(yy,yfit)
 print,' new Rmax = R ',correlate(yy,yfit),' old Rmax=',Rmax
 if (R gt Rmax) then begin
     Rmax=R
	outstr='ig(ca,boplusba,ig0,tsc,cat,frd):=('+string(const)+')+'
     print,'const=',const
     for k=0,n_elements(varlist)-1,1 do begin
         print,format='(a40,f10.5,a,g6.2,a,f5.1)',varnames(varlist(k)+1),dummy(k),' +/- ',abs(sigs(k)/dummy(k)*100.),'% or Z= ',abs(dummy(k)/sigs(k))
	 outstr=outstr+'+('+string(dummy(k),format='(f10.5)')+')*'+varnames(varlist(k)+1)
         endfor
	plot,yy,yfit,psym=1,xtitle='ig Observed',ytitle=' ig Model',charsize=1.7,title='Backwards elim.'
	oplot,yy(hdx),yfit(hdx),psym=1,color=fsc_color('red')
	openw,98,'formaxima.txt'
	printf,98,strcompress(outstr+';',/remove_all)
	printf,98,'diff(ig(ca,boplusba,ig0,tsc,cat,frd),ca);'
	printf,98,'diff(ig(ca,boplusba,ig0,tsc,cat,frd),boplusba);'
	printf,98,'diff(ig(ca,boplusba,ig0,tsc,cat,frd),ig0);'
	printf,98,'diff(ig(ca,boplusba,ig0,tsc,cat,frd),tsc);'
	printf,98,'diff(ig(ca,boplusba,ig0,tsc,cat,frd),cat);'
	printf,98,'diff(ig(ca,boplusba,ig0,tsc,cat,frd),frd);'

	printf,98,'quit();'
	close,98
     endif
 return
 end
 
 PRO dobackwardselimination,array,varnames_in
 common poi,numpointers
 varnames=varnames_in
 siglev=0.01
 yy=reform(array(0,*))
 l=size(array,/dimensions)
 nvars=l(0)-1
 xx=array(1:l(0)-1,*)
 dummy=backw_elim(xx,yy,siglev,varlist=varlist,const=const,yfit=yfit,sigmas=sigma)
 print,dummy
 print,'R: ',correlate(yy,yfit)
 print,'In pt_backw_REGRESS: varlist=',varnames(varlist)
 ;
 for i=0,numpointers-2,1 do begin
     for j=i+1,numpointers-2,1 do begin
         if (i ne j) then begin
             print,'---------------------------------------------'
             print,'VARNAMES inside i,j loop: ',varnames
             doablock,xx,yy,varnames,i,j,siglev
             print,'---------------------------------------------'
             endif
         endfor
     endfor
 return
 end
 
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
 
 PRO regress_with_validation,plotname,array,n_montecarlo,R_ensemblemean,R_ensemblemedian
 common varnames,varnames
 ;--------------------------------------------
 ; This routine performs multivariate regression with validation on bootstrapped withheld data
 ; ARRAY is INPUT and contains columns and rows - columns are the variables
 ; 	FIRST COLUMN is  y (dependent variable) and the rest are INdependent variables - i.e. x
 ; n_montecarlo is INPUT and is the number of Monte Carlo bootstraps to perform
 ; R_ensemblemean,R_ensemblemedian are OUTPUT and are the mean and median correlations obtained over the
 ; 	validation tsets
 ;-------------------------------------------
 ; First chop up array into Y and X
 l=size(array,/dimensions)
 nvarsplus1=l(0)
 nvars=nvarsplus1-1
 nrows=l(1)
 print,'nrows: ',nrows,' nvars: ',nvars
 r_list=[]
 r_list_2=[]
 ; Y
 y_orig=reform(array(0,*))
 ; X
 x_orig=reform(array(1:nvars,*))
 ; Loop over Monte Carlo trials
 z_list=[]
 for itrial=0L,n_montecarlo-1,1 do begin
     ; find pointer to the test sample and for the complementary withhold sample
     gofindpointers,nrows,idx,jdx
     coef=regress(x_orig(*,idx),y_orig(idx),/double,const=intercept,yfit=yhat,sigma=coef_si)
     z_list=[[z_list],[abs(coef_si/coef)]]
     r_list_2=[r_list_2,correlate(yhat,y_orig(idx))]
     ; build predicted y value for the withheld sample
     yhat_validate=intercept+reform(coef#x_orig(*,jdx))
     ; find and store correlation of the validation dataset
     r=correlate(y_orig(jdx),yhat_validate)
     r_list=[r_list,r]
     endfor
 R_ensemblemean=mean(r_list)
 R_ensemblemedian=median(r_list)
 !P.CHARSIZE=2
 histo,r_list,0,1,0.01,title=plotname
 pcolor=!P.color
 !P.color=fsc_color('red')
 histo,r_list_2,0,1,0.003,/overplot
 !P.color=pcolor
 print,'R mean on test sample  : ',mean(r_list_2)
 print,'R median on test sample: ',median(r_list_2)
 for j=0,n_elements(varnames)-2,1 do print,format='(a,a,a,f6.2)','Z for ',varnames(j+1),' : ',1./avg(z_list(j+1),1)
 return
 end
 
 ; Version 1 - uses the folded input file
 ; does backwards elimination of variables
 ; does original data vectors as well as all first-order cross terms
 common varnames,varnames
 common remember,Rmax
 Rmax=-1e22
 getthedatawewanttoselect,array,varnames
 n_montecarlo=500
 print,'-------------------------------------------------------------------------'
 regress_with_validation,files,array,n_montecarlo,R_ensemblemean,R_ensemblemedian
 print,'R mean on validation  : ',R_ensemblemean
 print,'R median on validation: ',R_ensemblemedian
 dobackwardselimination,array,varnames
 ;endfor
 print,'-------------------------------------------------------------------------'
 print,'Recall that results are means over MC trails - will change on next run!'
 print,'-------------------------------------------------------------------------'
 end
 
 
