PRO experiment,eidx,a,idx
 ; will return an index for a
 if (eidx eq 1) then begin
   print,'Morning experiment (6.00-11.00)'
   idx=where(a(5,*) eq 1 and a(7,*) gt 0.75 and a(7,*) lt 0.9583)
;   idx=where(a(7,*) gt 0.75 and a(7,*) lt 0.9583)
 endif
 if (eidx eq 2) then begin
   print,'Lunchtime experiment (11.00-15.00)'
   idx=where(a(5,*) eq 1 and (a(7,*) le 0.125 or a(7,*) gt 0.9583))
;   idx=where((a(7,*) le 0.125 or a(7,*) gt 0.9583))
 endif
 if (eidx eq 3) then begin
   print,'Afternoon/evening experiment (15.00-22.00)'
   idx=where(a(5,*) eq 1 and (a(7,*) gt 0.125 and a(7,*) lt 0.41667))
;   idx=where((a(7,*) gt 0.125 and a(7,*) lt 0.41667))
 endif
 if (eidx eq 4) then begin
   print,'Late evening/night experiment (22.00-06.00)'
   idx=where(a(5,*) eq 1 and (a(7,*) lt 0.75 and a(7,*) gt 0.41667))
;   idx=where((a(7,*) lt 0.75 and a(7,*) gt 0.41667))
 endif
 if (eidx eq 5) then begin
   print,'evening + late evening (20.00-24.00)'
   idx=where(a(5,*) eq 1 and (a(7,*) lt 0.5 and a(7,*) gt 0.333333))
;   idx=where((a(7,*) lt 0.75 and a(7,*) gt 0.41667))
 endif
 if (eidx eq 6) then begin
   print,'Breakfast experiment (cat=100)'
   idx=where(a(9,*) eq 100)
   help,idx
;   if (idx eq -1) then begin
;     stop
;   endif
 endif
 if (eidx eq 7) then begin
   print,'Lunch experiment (cat=200)'
   idx=where(a(9,*) eq 200)
   help,idx
 endif
 if (eidx eq 8) then begin
   print,'Dinner experiment (cat=300)'
   idx=where(a(9,*) eq 300)
   help,idx
 endif
 if (eidx eq 9) then begin
   print,'Evening/night experiment with NO CARBS (21.00-8.00)'
   delta=4./24.
   idx=where(a(7,*) gt 0.375 and a(7,*) le 0.83333 and a(8,*) gt delta)
   help,idx
stop
 endif
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
 !P.MULTI=[0,1,3]
 for itrial=0L,n_montecarlo-1,1 do begin
     ; find pointer to the test sample and for the complementary withhold sample
     gofindpointers,nrows,idx,jdx
     coef=regress(x_orig(*,idx),y_orig(idx),/double,const=intercept,yfit=yhat,sigma=coef_si)
     if (itrial eq 0) then begin
       delta=y_orig(idx)-yhat
       sd=stddev(delta)
idx_poor=where(abs(delta) gt sd*2)
idx_good=where(abs(delta) le sd*2)
       plot,xrange=[0,20],yrange=[0,20],xstyle=3,ystyle=3,y_orig(idx_good),yhat(idx_good),psym=7,xtitle='Observed IG',ytitle='Model'
       oplot,y_orig(idx_poor),yhat(idx_poor),color=fsc_color('red'),psym=7
       oplot,[0,20],[0,20],linestyle=2
       oplot,[0,20],[-SD*2,20-SD*2],linestyle=1
       oplot,[0,20],[SD*2,20+SD*2],linestyle=1

plot,yrange=[0,20],xstyle=3,ystyle=3,array(6,idx_good),yhat(idx_good),psym=7,xtitle='FIXME',ytitle='Model with first-order cross-products'
oplot,array(6,idx_poor),yhat(idx_poor),color=fsc_color('red'),psym=7
     endif
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
 
  ; Version 3 - uses the folded input file
  common varnames,varnames
  varnames=['ig','ca','bo+ba','f1','f2','f3','ig0','frd','tsc','cat','tsb','igb','bob']
  ;           0   1     2       3   4    5    6     7     8     9     10    11    12     
  pointers=[0,1,2,6,7,8,9,10,11,12]
  ;           0 1 2 3 4 5  6  7  8               
  varnames=varnames(pointers)
  files='kammas_folded_data.dat'
  print,'-------------------------------------------------------------------------'

  for expi=1,1,1 do begin
    array=get_data(files)
    l=size(array,/dimensions)
    nvars=l(0)
    idx=where (array(0,*) gt 2456963)
    arr=[array(1:nvars-1,idx)]	; skip the JD
    experiment,expi,arr,jdx
    arr=arr(*,jdx)
;--------------------------------
    l=size(arr,/dimensions)
    nvars=l(0)
    nrows=l(1)
    arr=arr([pointers],*)
    n_montecarlo=1000
    n_montecarlo=10
    help,arr
    regress_with_validation,files,arr,n_montecarlo,R_ensemblemean,R_ensemblemedian
    print,'R mean on validation  : ',R_ensemblemean
    print,'R median on validation: ',R_ensemblemedian
    print,'-------------------------------------------------------------------------'
    print,'Recall that results are means over MC trails - will change on next run!'
    print,'-------------------------------------------------------------------------'
  endfor
end
 
 
