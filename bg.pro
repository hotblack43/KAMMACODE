FUNCTION BG,carbs,bolus
return,      4.3948784433864898E-4*carbs^3-1.3043349374444699E-1*bolus*carbs^2+2.0840515441355301E+0*bolus^2*carbs+1.76947946194235E+0*bolus*carbs-2.9022514387288401E-1*carbs-3.6079253706890697E+1*bolus^2+1.0474792459322E+1*bolus+6.4704219357098101E+0
end

FUNCTION BG2,carbs,bolus
return,  5.77485114474295*(carbs < 13.4410291155654) + 1.33597671950457*min([0.467644549638431*carbs, 5.69319650783746])*max([bolus, (3.39305956505394*bolus > tan(-8.12257249384799*bolus*!dtor))]) - 0.205757030232629*carbs 
end


!P.MULTI=[0,1,2]
carbs=findgen(30)
bolus=findgen(20)/10.
surf=fltarr(n_elements(carbs),n_elements(bolus))
for i=0,n_elements(carbs)-1,1 do begin
for j=0,n_elements(bolus)-1,1 do begin
surf(i,j)=bg(carbs(i),bolus(j))
endfor
endfor
surf=surf > 0.0
contour,xstyle=3,ystyle=3,/cell_fill,nlevels=11,surf,carbs,bolus,xtitle='Carbs',ytitle='Bolus',/downhill,title='BG(carbs,Bolus)'
levs=[-1,0,2,4,6,7,8,10,12,14,20]
nlevs=n_elements(levs)
contour,min=0.0,levels=levs,c_labels=findgen(nlevs)*0+1,surf,carbs,bolus,/overplot,/downhill
data=get_data('eureqa_kamma_train.noheader')
oplot,data(1,*),data(2,*),psym=7
;
surf=fltarr(n_elements(carbs),n_elements(bolus))
for i=0,n_elements(carbs)-1,1 do begin
for j=0,n_elements(bolus)-1,1 do begin
surf(i,j)=bg2(carbs(i),bolus(j))
endfor
endfor
surf=surf > 0.0
contour,xstyle=3,ystyle=3,/cell_fill,nlevels=11,surf,carbs,bolus,xtitle='Carbs',ytitle='Bolus',/downhill,title='BG(carbs,Bolus)'
levs=[-1,0,2,4,6,7,8,10,12,14,20]
nlevs=n_elements(levs)
contour,min=0.0,levels=levs,c_labels=findgen(nlevs)*0+1,surf,carbs,bolus,/overplot,/downhill
data=get_data('eureqa_kamma_train.noheader')
oplot,data(1,*),data(2,*),psym=7
end
