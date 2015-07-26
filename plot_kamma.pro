!P.MULTI=[0,2,2]
data=get_data('kamma.data')
dt=reform(data(0,*))
carbs=reform(data(1,*))
bolus=reform(data(2,*))
ig=reform(data(3,*))
plot,carbs,ig,xtitle='hat-weighted carbs',ytitle='ig',psym=7
plot,bolus,ig,xtitle='hat-weighted bolus',ytitle='ig',psym=7
;------------------
data2=get_data('delta_kamma.data')
dt=reform(data2(0,*))
carbs=reform(data2(1,*))
bolus=reform(data2(2,*))
delta_ig=reform(data2(3,*))
plot,carbs,delta_ig,xtitle='hat-weighted carbs',ytitle='!7D!3!d4!n ig',psym=7
plot,bolus,delta_ig,xtitle='hat-weighted bolus',ytitle='!7D!3!d4!n ig',psym=7
end
