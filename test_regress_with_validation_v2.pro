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
 
	; Version 2 - uses the folded input file
 common varnames,varnames
 varnames=['ig','ca','bo+ba','f1','f2','f3','ig0','frd','tsc','cat']
 pointers=[0,1,2,6,7,8,9]
 varnames=varnames(pointers)
 files='kammas_folded_data.dat'
 print,'-------------------------------------------------------------------------'
 array=get_data(files)
 l=size(array,/dimensions)
 nvars=l(0)
 array=[array(1:nvars-1,*)]	; skip the JD
; select healthy/ill days
 idx=where(array(5,*) ne 0)
 array=array(*,idx)
;--------------------------------
 l=size(array,/dimensions)
 nvars=l(0)
 nrows=l(1)
 array=array([pointers],*)
 n_montecarlo=100
 help,array
 regress_with_validation,files,array,n_montecarlo,R_ensemblemean,R_ensemblemedian
 print,'R mean on validation  : ',R_ensemblemean
 print,'R median on validation: ',R_ensemblemedian
;endfor
 print,'-------------------------------------------------------------------------'
 print,'Recall that results are means over MC trails - will change on next run!'
 print,'-------------------------------------------------------------------------'
 end
 
 
