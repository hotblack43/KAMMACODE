!P.CHARSIZE=1.0
str=''
openr,1,'nameoffile.txt'
readf,1,str
close,1
data_in=get_data('all_data')
names=['B','V','VE1','VE2','IRCUT']
colnames=['blue','green','orange','brown','red']
id=names(data_in(15,*))
data=data_in([1,3,10,14],*)
; JD,albedo,erralbedo,alfa1,rlimit,pedestal,xshift,yshift,corefactor,contrast,RMSE,totfl,zodi,SLcounts,flux
;  0   1       2        3     4      5       6      7        8         9       10   11   12    13       14
l=size(data,/dimensions)
nvars=l(0)
nobs=l(1)

if_want_center=1
if_want_stdize=1
; across all observations do centering or standardization of each variable
if (if_want_center eq 1) then for ivar=0,nvars-1,1 do data(ivar,*)=data(ivar,*)-mean(data(ivar,*),/double)
if (if_want_stdize eq 1) then for ivar=0,nvars-1,1 do data(ivar,*)=data(ivar,*)/stddev(data(ivar,*),/double)
;
; Compute the Singular Value Decomposition:
SVDC, data, W, pc_mat, eof_mat
; Print the singular values:
totsing=total(w)
PRINT,format='(a,'+string(nvars)+'(1x,f9.6),a)', 'Sing val: ',W/totsing*100.,' % of var.'

; Reconstruct the input data med havelåger
sv = FLTARR(nvars, nvars)
FOR K = 0, nvars-1 DO sv[K,K] = W[K]
result = pc_mat ## sv ## TRANSPOSE(eof_mat)
;print,'DATA - RECONSTRUCTION:'
;PRINT,format='('+string(nvars)+'(1x,f9.6))', data-result
; reconstrcut with do-loops
; it seems that columns of eof_mat are EOFs and rows of pc_mat are PCs
;for irow=0,nobs-1,1 do begin
;sum=(pc_mat(0,irow)*eof_mat(0,*)*w(0)+pc_mat(1,irow)*eof_mat(1,*)*w(1)+pc_mat(2,irow)*eof_mat(2,*)*w(2)+pc_mat(3,irow)*eof_mat(3,*)*w(3)+pc_mat(4,irow)*eof_mat(4,*)*w(4)+pc_mat(5,irow)*eof_mat(5,*)*w(5)+pc_mat(6,irow)*eof_mat(6,*)*w(6))
;print,format='(i4,a,7(1x,g12.6))',irow,'  error ',data(*,irow)-sum
;endfor
;
!P.MULTI=[0,3,2]
!P.charsize=0.9
plot,xstyle=3,ystyle=3,total(w/totsing,/cum)*100,title=str,ytitle='Cumulated var explained',xtitle='EOF #',psym=7

plot,/nodata,/isotropic,pc_mat(0,*),pc_mat(1,*),xtitle='PC0',ytitle='PC1',psym=7,xstyle=3,ystyle=3
for iobs=0,nobs-1,1 do xyouts,pc_mat(0,iobs),pc_mat(1,iobs),id(iobs),color=fsc_color(colnames(data_in(15,iobs)));strcompress(string(iobs,format='(a15)'),/remove_all)

plot,/nodata,/isotropic,pc_mat(0,*),pc_mat(2,*),xtitle='PC0',ytitle='PC2',psym=7,xstyle=3,ystyle=3
for iobs=0,nobs-1,1 do xyouts,pc_mat(0,iobs),pc_mat(2,iobs),id(iobs),color=fsc_color(colnames(data_in(15,iobs)));strcompress(string(iobs,format='(a15)'),/remove_all)

plot,/nodata,/isotropic,pc_mat(1,*),pc_mat(2,*),xtitle='PC1',ytitle='PC2',psym=7,xstyle=3,ystyle=3
for iobs=0,nobs-1,1 do xyouts,pc_mat(1,iobs),pc_mat(2,iobs),id(iobs),color=fsc_color(colnames(data_in(15,iobs)));strcompress(string(iobs,format='(a15)'),/remove_all)
end
